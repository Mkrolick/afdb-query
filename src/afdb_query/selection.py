"""Choosing which AFDB structures answer a sequence query.

A sequence lookup can return several structures, and they are not interchangeable.
Two things make a hit unsuitable even at ``sequence_identity == 1.0``:

* **Complexes.** AFDB's summary reports ``confidence_avg_local_score`` over the whole
  deposited model. For a HOMODIMER or HETERODIMER that average spans every chain, not
  just the one matching the query, so it is not comparable with a monomer's. Roughly
  one matched structure in nine is such a complex.
* **Numeric AB-INITIO models.** ``AF-0000000065764032`` style entries are not the
  canonical per-UniProt ``AF-<accession>-F1`` model.

Two rules follow, and both matter.

**Preference never reads the confidence scores.** Ranking candidates by pLDDT and
returning the winner makes the fetch layer pick the best-scoring structure. That is
not merely arbitrary, it is *biased*: the expected maximum of N draws rises with N, so
a comparison whose two arms match different numbers of candidates gets different
inflation on each side. The tiers below use identity and provenance only.

**What survives the tiers is a GROUP, not an element.** After preference there is
routinely more than one candidate left -- measured on a real cache, 25.5% of records
end in a tie. Breaking that tie by any rule at all discards N-1 predictions of the same
sequence for no reason, so :func:`select_group` returns all of them and
:func:`mean_global_plddt` averages across them. Averaging is unbiased for the same
reason an arbitrary pick is -- neither consults the value -- and it has strictly lower
variance, because it uses every prediction rather than throwing them away.

Tied candidates are near-always the same protein reached through different UniProt
entries: orthologs whose sequence is identical, isoform duplicates, TrEMBL redundancy.
On a real cache 4,288 of 4,289 tied sets consisted entirely of full-length ``-F1``
models, so their per-residue arrays are the same length and average elementwise.
"""

from __future__ import annotations

import re

# A canonical per-UniProt model id like "AF-P12345-F1" -- the accession starts with a
# letter. Numeric AB-INITIO ids ("AF-0000000065764032") do not match.
CANONICAL_F1 = re.compile(r"^AF-[A-Za-z]\w*-F\d+$")

MONOMER = "MONOMER"

# AFDB's sequence endpoint is an exact-sequence lookup: it returns structures whose sequence
# IS the query, reporting sequence_identity 1.0 and coverage 1.0 for every one, and 404s
# rather than offering a near neighbour. Verified two ways -- across 60,436 cached summaries
# from a real run, all 39,677 matched structures reported exactly these values; and against
# the live API across proteins of length 110-1273.
#
# Requiring them anyway, because that is a fact about a third-party API and not a guarantee
# this package can make. Every number downstream assumes the structure attached to a row IS
# that row's protein, and sequence_identity is the only field that ever asserts it. If AFDB
# ever adds fuzzy matching -- a new parameter default, a schema change -- an unchecked
# pipeline would attach a different protein's confidence to every affected row and look
# entirely normal doing it.
EXACT_IDENTITY = 1.0
EXACT_COVERAGE = 1.0


def is_exact(summary: dict) -> bool:
    """Whether this structure IS the queried sequence rather than a near neighbour.

    Requires ``sequence_identity == 1.0`` **and** ``coverage == 1.0``: identity alone says the
    aligned part matched, coverage says the alignment spanned the whole query. Both are needed
    for "this structure is my protein".

    A missing or non-numeric field is NOT exact. Unverifiable and verified are different
    states, and only one of them is safe to attach a pLDDT to.

    Note this still does not mean the model is the query's SIZE -- coverage is about the query
    being covered, not about the model carrying nothing else. See :func:`filter_by_length`.
    """
    return (
        summary.get("sequence_identity") == EXACT_IDENTITY
        and summary.get("coverage") == EXACT_COVERAGE
    )


def exact_matches(structures):
    """Split the endpoint's ``structures`` list into ``(exact, rejected)`` summaries.

    Exposed rather than kept internal so a caller can see what :func:`select_group` discarded.
    Today that list is always empty; if it ever is not, AFDB's search contract has changed and
    that is worth knowing loudly rather than inferring from a drop in coverage.
    """
    summaries = [item["summary"] for item in structures if item.get("summary")]
    exact = [s for s in summaries if is_exact(s)]
    rejected = [s for s in summaries if not is_exact(s)]
    return exact, rejected


def is_monomer(summary: dict) -> bool:
    """Whether this structure is a single-chain model of one entity.

    Checked from the summary alone -- no extra request. ``oligomeric_state`` carries
    AFDB's own verdict; the ``chain_ids`` count is the corroborating signal for entries
    that leave the state unset.
    """
    if summary.get("oligomeric_state") not in (MONOMER, None):
        return False
    entities = summary.get("entities") or []
    chains = sum(len(e.get("chain_ids") or []) for e in entities)
    return chains <= 1


def is_canonical_model(summary: dict) -> bool:
    """Whether the model id is a canonical ``AF-<accession>-F<n>`` rather than numeric."""
    return bool(CANONICAL_F1.match(summary.get("model_identifier") or ""))


def rank_tiers(summary: dict, accession: str | None = None) -> tuple[int, int, int]:
    """Preference tiers for one candidate (lower sorts first).

    1. the caller's ``accession`` (``AF-<accession>-F1``), then
    2. monomers over complexes, then
    3. canonical ``-F1`` ids over numeric AB-INITIO ids.

    There is no fourth tier. Confidence scores are absent from all three -- see the
    module docstring -- and so is any arbitrary tie-break, because candidates tying here
    are meant to be kept rather than whittled down to one.
    """
    mid = summary.get("model_identifier") or ""
    return (
        0 if (accession and mid == f"AF-{accession}-F1") else 1,
        0 if is_monomer(summary) else 1,
        0 if is_canonical_model(summary) else 1,
    )


def select_group(structures, accession: str | None = None) -> list[dict]:
    """Every candidate tied at the best :func:`rank_tiers`, ordered by model id.

    ``structures`` is the endpoint's list of ``{"summary": {...}}`` items; the return is
    the bare summary dicts. Empty when there are no candidates at all.

    The ordering exists for reproducible output only. Do NOT take ``[0]`` as "the"
    structure -- that reintroduces the arbitrary choice this function exists to avoid.
    Pass the whole group to :func:`mean_global_plddt`, or fetch each member's
    per-residue array and average those elementwise.

    **Non-exact matches are discarded before the tiers run** -- see :func:`is_exact`. That is
    the one thing this filters, and it filters rather than raises because a row with no exact
    match should end up with no pLDDT, which every consumer already handles, rather than
    halting a 60k-row batch on one anomaly. Use :func:`exact_matches` to see what was dropped.

    Nothing else is filtered. "No usable structure" and "a group whose average spans a
    complex" warrant different handling and only the caller knows which its analysis
    tolerates, so test :func:`is_monomer` on the members and decide.
    """
    summaries, _rejected = exact_matches(structures)
    if not summaries:
        return []
    best = min(rank_tiers(s, accession) for s in summaries)
    group = [s for s in summaries if rank_tiers(s, accession) == best]
    return sorted(group, key=lambda s: s.get("model_identifier") or "")


def filter_by_length(summaries, lengths, expected_length: int):
    """Split a group into members that model exactly ``expected_length`` residues, and the rest.

    Returns ``(kept, dropped)``. ``lengths`` maps ``model_identifier`` to that model's
    residue count -- which the summary does NOT carry, so the caller supplies it after
    fetching per-residue confidence. A member whose length is unknown is dropped, not
    assumed to conform.

    This exists because the tiers cannot see length and ``coverage == 1.0`` does not
    imply it. Coverage means the QUERY is fully covered by the model, not that the model
    is the query's size -- which is exactly how an 860-residue multi-chain entry reports
    coverage 1.0 against a 430-residue query. Filtering to monomers removes that case,
    but an ortholog carrying the query sequence plus a few extra terminal residues
    passes every visible test and still breaks two things:

    * **Positional slicing.** An offset computed from the query's amino-acid length
      indexes the wrong residues in a longer array -- silently, with a plausible result.
    * **The group average.** A longer member's ``confidence_avg_local_score`` spans
      residues the query does not contain, which is the complex-contamination problem
      in a subtler form.

    Dropped members are returned rather than discarded so a caller can report how many
    it lost. "3 of 5 matched entries were the wrong length" is a fact about the data,
    not a detail to absorb.
    """
    kept, dropped = [], []
    for summary in summaries:
        length = lengths.get(summary.get("model_identifier"))
        (kept if length == expected_length else dropped).append(summary)
    return kept, dropped


def mean_global_plddt(summaries) -> float | None:
    """Mean ``confidence_avg_local_score`` across a group, or ``None`` if none carry one.

    The unbiased replacement for picking one member. It costs no extra request -- the
    field is already in every summary -- and it is what a caller comparing one protein's
    confidence against another's should use in place of any single member's value.

    Members without the field are skipped rather than counted as zero: a missing score
    is an absent measurement, and averaging it in as 0 would drag the group toward a
    value no structure reported.

    **This cannot verify that the members model your query and nothing more.** Each
    member's field is a mean over that model's own residues, and the summary carries no
    residue count, so a member larger than the query contributes an average partly over
    residues you did not ask about. Run :func:`filter_by_length` first when the
    comparison depends on it.
    """
    values = [
        s.get("confidence_avg_local_score")
        for s in summaries
        if isinstance(s.get("confidence_avg_local_score"), (int, float))
    ]
    if not values:
        return None
    return sum(values) / len(values)
