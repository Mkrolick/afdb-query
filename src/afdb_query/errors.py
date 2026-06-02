"""Exception types for afdb-query."""

from __future__ import annotations


class AFDBError(Exception):
    """Base class for all afdb-query errors."""


class InvalidSequenceError(AFDBError):
    """Raised when a sequence cannot be queried against AFDB.

    ``reason`` is one of ``"internal_stop"``, ``"too_short"``,
    ``"nonstandard_aa"`` (see ``afdb_query.sequences.filter_reason``).
    """

    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(f"sequence not queryable: {reason}")
