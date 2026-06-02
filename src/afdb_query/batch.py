"""Concurrent, resumable batch lookups over many sequences."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple

import httpx

from .sequences import filter_reason


def _normalize_inputs(inputs) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for item in inputs:
        if isinstance(item, dict):
            pairs.append((item["id"], item["sequence"]))
        else:
            id_, seq = item
            pairs.append((id_, seq))
    return pairs


def _chunked(items: list, size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


class _Result(NamedTuple):
    summary_path: Path
    plddt_path: Path | None
    outcome: str  # "found" | "notfound" | "error"
    summary_data: dict | None
    plddt_values: list | None


def search_many(
    client,
    inputs,
    out_dir,
    *,
    concurrency: int = 6,
    rows: int = 10,
    plddt_first_n: int | None = None,
) -> dict:
    """Query each queryable input's sequence concurrently, caching to disk.

    ``inputs`` is a list of ``(id, sequence)`` tuples or ``{"id":..., "sequence":...}``
    dicts. Results are cached under ``out_dir``:

    * ``out_dir/summaries/{id}.json`` — a hit stores the AFDB summary document; a
      404 miss stores ``{"structures": []}`` so re-runs skip it. An existing file
      is left untouched (resumability).
    * ``out_dir/plddt/{id}.json`` (only when ``plddt_first_n`` is set) — the raw
      first-n per-residue pLDDT array for the first/best structure (<= n values).

    A real per-query HTTP error is counted and NOT saved, so it retries next run.
    Returns a counts report (dict).

    Note: resumability keys on the summary file. If a record's summary file already
    exists, it is skipped entirely and ``plddt`` is not back-filled for it.
    """
    out_dir = Path(out_dir)
    summaries_dir = out_dir / "summaries"
    plddt_dir = out_dir / "plddt"

    pairs = _normalize_inputs(inputs)
    counts = {
        "internal_stop": 0,
        "too_short": 0,
        "nonstandard_aa": 0,
        "skipped": 0,
        "hits": 0,
        "misses": 0,
        "errors": 0,
    }

    pending: list[tuple[Path, Path | None, str]] = []
    for id_, seq in pairs:
        reason = filter_reason(seq)
        if reason is not None:
            counts[reason] += 1
            continue
        summary_path = summaries_dir / f"{id_}.json"
        if summary_path.exists():
            counts["skipped"] += 1
            continue
        plddt_path = (plddt_dir / f"{id_}.json") if plddt_first_n is not None else None
        pending.append((summary_path, plddt_path, seq))

    def _query(item: tuple[Path, Path | None, str]) -> _Result:
        summary_path, plddt_path, seq = item
        try:
            data = client._fetch_summary(seq, rows)
        except httpx.HTTPError:
            return _Result(summary_path, plddt_path, "error", None, None)
        if data is None:
            return _Result(summary_path, plddt_path, "notfound", None, None)
        plddt_values = None
        if plddt_first_n is not None:
            structures = data.get("structures") or []
            if structures:
                try:
                    conf = client._fetch_confidence(structures[0]["summary"]["model_url"])
                except httpx.HTTPError:
                    return _Result(summary_path, plddt_path, "error", None, None)
                plddt_values = conf.get("confidenceScore", [])[:plddt_first_n]
        return _Result(summary_path, plddt_path, "found", data, plddt_values)

    if pending:
        summaries_dir.mkdir(parents=True, exist_ok=True)
        if plddt_first_n is not None:
            plddt_dir.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=concurrency) as pool:
            for chunk in _chunked(pending, concurrency * 50):
                for res in pool.map(_query, chunk):
                    if res.outcome == "error":
                        counts["errors"] += 1
                    elif res.outcome == "notfound":
                        res.summary_path.write_text(json.dumps({"structures": []}))
                        counts["misses"] += 1
                    else:
                        res.summary_path.write_text(json.dumps(res.summary_data))
                        counts["hits"] += 1
                        if res.plddt_values is not None and res.plddt_path is not None:
                            res.plddt_path.write_text(json.dumps(res.plddt_values))

    return {
        "total": len(pairs),
        "filtered": counts["internal_stop"] + counts["too_short"] + counts["nonstandard_aa"],
        **counts,
        "queried": counts["hits"] + counts["misses"] + counts["errors"],
    }
