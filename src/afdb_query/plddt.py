"""Statistics over per-residue pLDDT -- one definition of "mean pLDDT".

AFDB reports a summary field, ``confidence_avg_local_score``, that is already a mean
over the deposited model's residues. Other predictors (ColabFold, ESMFold) emit only
the per-residue array and leave the averaging to you. A codebase that reads the field
in one place and averages an array in another has two definitions of its headline
number and no way to tell them apart. :func:`mean_plddt` is the one definition.

Nothing here rounds. AFDB's own field arrives pre-rounded by the server, so
``mean_plddt(structure.plddt().scores)`` can differ from ``structure.global_plddt`` in
the last digit; that is a display difference, and rounding inside the computation
would bake it in. Round at the point you format, not here.

**Region means are the reason this takes bounds.** Comparing ``mean(whole protein)``
against ``mean(a sub-range of the same protein)`` measures the sub-range's composition
as much as anything else: drop a disordered N-terminal tail and the mean rises purely
because low scorers left the set. When two proteins share a region, compare over that
region on both sides -- ``mean_plddt(a, start=offset)`` against ``mean_plddt(b)`` --
so residue count cancels.
"""

from __future__ import annotations

from collections.abc import Sequence


def mean_plddt(
    scores: Sequence[float], start: int | None = None, stop: int | None = None
) -> float | None:
    """Mean of ``scores[start:stop]``, or ``None`` when that slice is empty.

    ``start`` / ``stop`` are ordinary Python slice bounds over the per-residue array,
    zero-based and ``stop``-exclusive. They index the ARRAY, not AFDB's residue
    numbering -- see :func:`residue_index`, and do not assume the two coincide.

    ``None`` rather than 0.0 or a raise: an empty selection has no mean, and a caller
    comparing two proteins needs to distinguish "no overlap" from "overlap that scored
    zero" without a try/except around every comparison.
    """
    window = scores[start:stop]
    if not window:
        return None
    return sum(window) / len(window)


def residue_index(residue_numbers: Sequence[int], residue: int) -> int | None:
    """Array position of AFDB residue number ``residue``, or ``None`` if absent.

    Per-residue arrays are parallel to ``residueNumber`` and are contiguous 1..N for
    every AFDB entry seen so far -- but that is an observation, not a guarantee the
    format makes, and a caller slicing by residue number on the assumption it holds
    would misalign silently rather than fail. Going through this function costs a
    lookup and removes the assumption.
    """
    try:
        return list(residue_numbers).index(residue)
    except ValueError:
        return None


def is_contiguous(residue_numbers: Sequence[int]) -> bool:
    """Whether ``residue_numbers`` is exactly ``1..N`` with no gaps.

    When true, array position ``i`` is residue ``i + 1`` and slicing by position is
    safe. Assert this before treating the two as interchangeable in bulk.
    """
    return list(residue_numbers) == list(range(1, len(residue_numbers) + 1))


def shared_suffix_means(
    long_scores: Sequence[float], short_scores: Sequence[float]
) -> dict[str, float | None]:
    """Length-controlled comparison of a protein against a suffix of itself.

    When one protein is the other with its N-terminus removed -- an alternative
    downstream start codon, say -- the naive comparison of two global means is
    confounded: it contrasts a mean over all residues with a mean over a subset of the
    same residues, and drops exactly the N-terminal region that scores lowest. The
    difference is then largely a readout of how much was removed.

    Returns three quantities that are not:

    * ``shared_long`` / ``shared_short`` -- means over the residues the two have in
      common, the same count on both sides, so length cancels. Their difference is
      whether truncation changed confidence in the part that survived.
    * ``displaced`` -- the mean over the residues only the longer protein has. A low
      value says the removed segment was the disordered tail, which is what makes the
      naive comparison move; a high one says real structure was removed.

    ``offset`` is ``len(long) - len(short)``. The caller is responsible for having
    established that the shorter sequence really is a suffix of the longer one -- this
    function aligns by length alone and cannot tell a suffix from an unrelated protein
    that happens to be shorter.
    """
    offset = len(long_scores) - len(short_scores)
    if offset < 0:
        raise ValueError("long_scores must be at least as long as short_scores")
    return {
        "offset": offset,
        "shared_long": mean_plddt(long_scores, start=offset),
        "shared_short": mean_plddt(short_scores),
        "displaced": mean_plddt(long_scores, stop=offset),
    }


def mean_per_residue(arrays, expected_length: int | None = None) -> list[float]:
    """Elementwise mean across a group's per-residue arrays. Raises on ragged input.

    The consensus array for a group of structures predicting the same sequence -- see
    ``afdb_query.selection``, which returns such a group rather than picking one member.
    Averaging across the group is unbiased and lower-variance than taking any single
    member, and it commutes with slicing: the mean over a region of this array equals
    the mean of each member's region mean, because both are unweighted means over the
    same rectangular index set.

    **That commutation is exactly why ragged input must raise rather than be handled.**
    A shorter member makes the index set non-rectangular, so position ``i`` stops
    denoting the same residue in every member and every downstream region mean silently
    describes a comparison that does not exist. ``zip`` would truncate to the shortest
    array and return a plausible number for it; that is the failure this refuses.

    ``expected_length``, when given, additionally requires the arrays to be that long --
    the query's residue count. ``coverage == 1.0`` from AFDB means the query is fully
    covered by the model, not that the model is the query's size, so a member can be
    uniformly longer and still align internally. Pair with
    ``afdb_query.selection.filter_by_length`` to drop such members before they get here.

    Raises ``ValueError`` on ragged input, on a length mismatch, or on no input at all.
    """
    arrays = [list(a) for a in arrays]
    if not arrays:
        raise ValueError("mean_per_residue needs at least one array")
    lengths = {len(a) for a in arrays}
    if len(lengths) > 1:
        raise ValueError(f"ragged per-residue arrays: lengths {sorted(lengths)}")
    length = lengths.pop()
    if expected_length is not None and length != expected_length:
        raise ValueError(f"arrays are {length} residues, expected {expected_length}")
    if length == 0:
        raise ValueError("per-residue arrays are empty")
    n = len(arrays)
    return [sum(a[i] for a in arrays) / n for i in range(length)]
