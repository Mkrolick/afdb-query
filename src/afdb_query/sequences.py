"""Sequence validation for AFDB queries."""

from __future__ import annotations

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

# AFDB's sequence endpoint rejects queries shorter than this, so filtering here
# saves a request that could only 400. It is a property of the endpoint, not a
# judgement about which proteins are interesting -- a caller studying short ORFs
# should know these were never queried, which is why `filter_reason` names the
# reason rather than returning a bare bool, and why `search_many` reports
# `filtered["too_short"]` separately from `queried["misses"]`.
MIN_LENGTH = 20


def filter_reason(seq: str) -> str | None:
    """Why a sequence cannot be queried against AFDB, or None if it is queryable.

    Checked in priority order: internal stop, length, non-standard residues. The
    priority matters only for reporting -- a sequence failing two checks is counted
    under the first -- and the order runs cheapest and most-definitive first.
    """
    if "*" in seq:
        return "internal_stop"
    if len(seq) < MIN_LENGTH:
        return "too_short"
    if not set(seq) <= STANDARD_AA:
        return "nonstandard_aa"
    return None
