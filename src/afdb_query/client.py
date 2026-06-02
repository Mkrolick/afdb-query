"""The AlphaFold client: HTTP session and AFDB endpoint access."""

from __future__ import annotations

import httpx

from .batch import search_many as _search_many
from .errors import InvalidSequenceError
from .models import Structure, confidence_url
from .sequences import filter_reason

DEFAULT_BASE_URL = "https://alphafold.ebi.ac.uk"
SUMMARY_PATH = "/api/sequence/summary"


class AlphaFold:
    """Client for sequence-based access to the AlphaFold Protein Structure Database.

    Wraps a shared, thread-safe ``httpx.Client``. Use as a context manager, or
    call :meth:`close` when done.
    """

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        base_url: str = DEFAULT_BASE_URL,
        max_retries: int = 2,
    ) -> None:
        transport = httpx.HTTPTransport(retries=max_retries)
        self._client = httpx.Client(base_url=base_url, timeout=timeout, transport=transport)

    # -- lifecycle ---------------------------------------------------------
    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> "AlphaFold":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # -- low-level fetch ---------------------------------------------------
    def _get(self, url: str, params: dict | None = None) -> httpx.Response:
        return self._client.get(url, params=params, headers={"Accept": "application/json"})

    def _fetch_summary(self, sequence: str, rows: int = 10) -> dict | None:
        """Tier 1: query the sequence-summary endpoint.

        Returns the parsed ``{"entry": ..., "structures": [...]}`` document, or
        ``None`` when AFDB has no entry (HTTP 404 — a clean "not found"). Raises
        on any other HTTP error.
        """
        if rows < 2:
            raise ValueError("rows must be > 1 (AFDB rejects rows <= 1)")
        resp = self._get(
            SUMMARY_PATH, params={"id": sequence, "type": "sequence", "rows": rows}
        )
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        return resp.json()

    def _fetch_confidence(self, model_url: str) -> dict:
        """Tier 2: fetch the per-residue confidence JSON for a model URL."""
        resp = self._get(confidence_url(model_url))
        resp.raise_for_status()
        return resp.json()

    # -- public API --------------------------------------------------------
    def search(self, sequence: str, rows: int = 10) -> list[Structure]:
        """Tier 1: find AFDB structures matching ``sequence``, in AFDB's returned order.

        Results are ranked by sequence identity, but ``hits[0]`` is not guaranteed to
        be the canonical ``AF-<accession>-F1`` model — for some sequences a multi-chain
        or AB-INITIO model ranks first. Select by ``model_identifier`` if you need a
        specific entry.

        Raises :class:`InvalidSequenceError` if the sequence is not queryable.
        Returns ``[]`` when AFDB has no entry for it.
        """
        reason = filter_reason(sequence)
        if reason is not None:
            raise InvalidSequenceError(reason)
        data = self._fetch_summary(sequence, rows)
        if data is None:
            return []
        return [Structure(item["summary"], self) for item in data.get("structures", [])]

    def search_many(
        self,
        inputs,
        out_dir,
        *,
        concurrency: int = 6,
        rows: int = 10,
        plddt_first_n: int | None = None,
    ) -> dict:
        """Concurrent, resumable batch lookup. See ``afdb_query.batch.search_many``."""
        return _search_many(
            self,
            inputs,
            out_dir,
            concurrency=concurrency,
            rows=rows,
            plddt_first_n=plddt_first_n,
        )
