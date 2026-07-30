"""Concurrent, resumable batch lookups over many sequences.

Two passes, cached side by side under one ``out_dir`` and resumable independently:

* :func:`search_many` -- sequence -> summary documents (``out_dir/summaries/{id}.json``)
* :func:`fetch_plddt_many` -- model URL -> per-residue pLDDT (``out_dir/plddt/{id}.json``)

They are separate because choosing WHICH of a summary's structures to fetch residues
for is the caller's decision (see :mod:`afdb_query.selection`), and because each keys
resumability on its OWN artifact. A summary-only run followed by a per-residue run
back-fills every already-cached record rather than skipping it.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple

from .errors import AFDBHTTPError
from .models import Plddt
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
    the cached summaries, then :func:`fetch_plddt_many`, when per-residue data is needed.

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


# -- per-residue pass ------------------------------------------------------


def _normalize_plddt_inputs(inputs) -> list[tuple[str, str]]:
    """``(id, model_url)`` pairs from either dicts or tuples, in the caller's order."""
    pairs: list[tuple[str, str]] = []
    for item in inputs:
        if isinstance(item, dict):
            pairs.append((item["id"], item["model_url"]))
        else:
            pairs.append((item[0], item[1]))
    return pairs


def plddt_path(out_dir, id_: str) -> Path:
    """Where :func:`fetch_plddt_many` caches record ``id_``'s per-residue array."""
    return Path(out_dir) / "plddt" / f"{id_}.json"


def load_plddt(out_dir, id_: str) -> Plddt | None:
    """The cached :class:`Plddt` for ``id_``, or ``None`` when it was never fetched.

    Reading goes through here rather than through ``json.load`` at each call site so
    the on-disk shape has exactly one reader, and so a record that was never fetched is
    distinguishable from one whose fetch returned nothing -- both are ``None`` to the
    caller, but only the first leaves no file, which is what makes a re-run resume it.
    """
    path = plddt_path(out_dir, id_)
    if not path.exists():
        return None
    return Plddt.from_dict(json.loads(path.read_text()))


def fetch_plddt_many(client, inputs, out_dir, *, concurrency: int = 6) -> dict:
    """Fetch and cache the FULL per-residue pLDDT array for each ``(id, model_url)``.

    ``inputs`` is a list of ``(id, model_url)`` tuples or ``{"id":..., "model_url":...}``
    dicts -- typically built by running :func:`afdb_query.selection.select` over
    summaries that :func:`search_many` already cached, which is what keeps the choice of
    structure visible in the caller rather than buried here.

    The whole array is stored, never a prefix: which residues a question needs is not
    knowable from this side, and a truncated array cannot be widened later without
    re-fetching.

    Resumability keys on ``out_dir/plddt/{id}.json``, NOT on the summary. Running
    summaries first and residues later back-fills every record.

    A failure is counted and not written, so it retries next run. Returns a report of
    disjoint counts::

        {"total": int, "skipped": int, "queried": {"hits", "errors", "total"}}
    """
    out_dir = Path(out_dir)
    pairs = _normalize_plddt_inputs(inputs)
    queried = {"hits": 0, "errors": 0}
    skipped = 0

    pending: list[tuple[Path, str]] = []
    for id_, model_url in pairs:
        path = plddt_path(out_dir, id_)
        if path.exists():
            skipped += 1
            continue
        pending.append((path, model_url))

    def _query(item: tuple[Path, str]):
        path, model_url = item
        try:
            doc = client.fetch_confidence(model_url)
        except AFDBHTTPError:
            return path, None
        return path, {
            "model_url": model_url,
            "confidenceScore": doc.get("confidenceScore") or [],
            "residueNumber": doc.get("residueNumber") or [],
        }

    if pending:
        (out_dir / "plddt").mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=concurrency) as pool:
            for chunk in _chunked(pending, concurrency * 50):
                for path, payload in pool.map(_query, chunk):
                    if payload is None:
                        queried["errors"] += 1
                    else:
                        path.write_text(json.dumps(payload))
                        queried["hits"] += 1

    return {
        "total": len(pairs),
        "skipped": skipped,
        "queried": {**queried, "total": sum(queried.values())},
    }
