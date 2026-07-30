"""Choosing which AFDB structure answers a sequence query.

A sequence lookup can return several structures, and they are not interchangeable.
Two things make a hit unsuitable even at ``sequence_identity == 1.0``:

* **Complexes.** AFDB's summary reports ``confidence_avg_local_score`` over the whole
  deposited model. For a HOMODIMER or HETERODIMER that average spans every chain, not
  just the one matching the query, so it is not comparable with a monomer's. Roughly
  one matched structure in nine is such a complex.
* **Numeric AB-INITIO models.** ``AF-0000000065764032`` style entries are not the
  canonical per-UniProt ``AF-<accession>-F1`` model.

Selection here is deliberately **independent of the confidence scores**. Ranking
candidates by pLDDT and returning the winner would make the fetch layer pick the
best-scoring structure, which silently biases any downstream comparison of one
protein's pLDDT against another's. The order below uses identity and provenance
only; where it cannot decide, it says so rather than reaching for the score.
"""

from __future__ import annotations

import re

# A canonical per-UniProt model id like "AF-P12345-F1" -- the accession starts with a
# letter. Numeric AB-INITIO ids ("AF-0000000065764032") do not match.
CANONICAL_F1 = re.compile(r"^AF-[A-Za-z]\w*-F\d+$")

MONOMER = "MONOMER"


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


def rank_key(summary: dict, accession: str | None = None):
    """Deterministic preference order for one candidate (lower sorts first).

    1. the caller's ``accession`` (``AF-<accession>-F1``), then
    2. monomers over complexes, then
    3. canonical ``-F1`` ids over numeric AB-INITIO ids, then
    4. ``model_identifier`` lexicographically.

    Confidence scores are deliberately absent -- see the module docstring.
    """
    mid = summary.get("model_identifier") or ""
    return (
        0 if (accession and mid == f"AF-{accession}-F1") else 1,
        0 if is_monomer(summary) else 1,
        0 if is_canonical_model(summary) else 1,
        mid,
    )


def select(structures, accession: str | None = None):
    """The best candidate summary by :func:`rank_key`, or ``None`` when there are none.

    ``structures`` is the endpoint's list of ``{"summary": {...}}`` items. Returns the
    bare summary dict. This never filters -- a caller that requires a monomer should
    test :func:`is_monomer` on the result and decide, because "no usable structure" and
    "a structure whose average spans a complex" warrant different handling and only the
    caller knows which its analysis can tolerate.
    """
    summaries = [item["summary"] for item in structures if item.get("summary")]
    if not summaries:
        return None
    return min(summaries, key=lambda s: rank_key(s, accession))
