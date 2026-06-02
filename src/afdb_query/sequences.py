"""Sequence validation for AFDB queries (ported from the original pipeline)."""

from __future__ import annotations

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH = 20


def filter_reason(seq: str) -> str | None:
    """Why a sequence cannot be queried against AFDB, or None if it is queryable.

    Checked in priority order: internal stop, length, non-standard residues.
    """
    if "*" in seq:
        return "internal_stop"
    if len(seq) < MIN_LENGTH:
        return "too_short"
    if not set(seq) <= STANDARD_AA:
        return "nonstandard_aa"
    return None
