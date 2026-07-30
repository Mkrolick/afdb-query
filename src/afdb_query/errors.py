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


class AFDBHTTPError(AFDBError):
    """Raised when an AFDB request fails at the transport or HTTP level.

    The underlying ``httpx`` exception stays reachable as ``__cause__``, so a
    caller can handle everything this package raises by catching
    :class:`AFDBError` alone, without importing ``httpx``.

    A 404 from the sequence endpoint is NOT an error: it is AFDB's clean "no
    entry for this sequence" and surfaces as ``None`` / ``[]`` instead.
    """
