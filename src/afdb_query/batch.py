"""Concurrent, resumable batch lookups over many sequences."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple

from .errors import AFDBHTTPError
from .sequences import filter_reason


def _normalize_inputs(inputs) -> list[tuple[str, str]]:
    """``(id, sequence)`` pairs from either dicts or tuples, in the caller's order."""
    pairs: list[tuple[str, str]] = []
    for item in inputs:
        if isinstance(item, dict):
            pairs.append((item["id"], item["sequence"]))
        else:
            pairs.append((item[0], item[1]))
    return pairs


def _chunked(items: list, size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


class _Result(NamedTuple):
    summary_path: Path
    outcome: str  # "found" | "notfound" | "error"
    summary_data: dict | None


def search_many(client, inputs, out_dir, *, concurrency: int = 6, rows: int = 10) -> dict:
    """Query each queryable input's sequence concurrently, caching summaries to disk.

    ``inputs`` is a list of ``(id, sequence)`` tuples or ``{"id":..., "sequence":...}``
    dicts. Results are cached under ``out_dir/summaries/{id}.json``: a hit stores the
    AFDB summary document, a 404 miss stores ``{"structures": []}`` so re-runs skip it.
    An existing file is left untouched, which is what makes the run resumable.

    A per-query HTTP failure is counted and NOT written, so it is retried next run.

    This function fetches summaries only. It does not choose a structure and it does
    not fetch per-residue confidence: which of several exact-sequence matches is the
    right one is a question about the caller's analysis, not about AFDB, so making
    that choice here would hide it. Run :func:`afdb_query.selection.select_group` over
    the cached summaries and average across what it returns.

    Returns a report of disjoint counts::

        {
          "total":    int,                        # inputs seen
          "skipped":  int,                        # already cached, not re-queried
          "filtered": {"internal_stop", "too_short", "nonstandard_aa", "total"},
          "queried":  {"hits", "misses", "errors", "total"},
        }

    ``total == skipped + filtered["total"] + queried["total"]``. No count appears in
    more than one place.
    """
    out_dir = Path(out_dir)
    summaries_dir = out_dir / "summaries"

    pairs = _normalize_inputs(inputs)
    filtered = {"internal_stop": 0, "too_short": 0, "nonstandard_aa": 0}
    queried = {"hits": 0, "misses": 0, "errors": 0}
    skipped = 0

    pending: list[tuple[Path, str]] = []
    for id_, seq in pairs:
        reason = filter_reason(seq)
        if reason is not None:
            filtered[reason] += 1
            continue
        summary_path = summaries_dir / f"{id_}.json"
        if summary_path.exists():
            skipped += 1
            continue
        pending.append((summary_path, seq))

    def _query(item: tuple[Path, str]) -> _Result:
        summary_path, seq = item
        try:
            data = client.fetch_summary(seq, rows)
        except AFDBHTTPError:
            return _Result(summary_path, "error", None)
        if data is None:
            return _Result(summary_path, "notfound", None)
        return _Result(summary_path, "found", data)

    if pending:
        summaries_dir.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=concurrency) as pool:
            for chunk in _chunked(pending, concurrency * 50):
                for res in pool.map(_query, chunk):
                    if res.outcome == "error":
                        queried["errors"] += 1
                    elif res.outcome == "notfound":
                        res.summary_path.write_text(json.dumps({"structures": []}))
                        queried["misses"] += 1
                    else:
                        res.summary_path.write_text(json.dumps(res.summary_data))
                        queried["hits"] += 1

    return {
        "total": len(pairs),
        "skipped": skipped,
        "filtered": {**filtered, "total": sum(filtered.values())},
        "queried": {**queried, "total": sum(queried.values())},
    }
