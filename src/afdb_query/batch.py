"""Concurrent, resumable batch lookups over many sequences."""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple

import httpx

from .sequences import filter_reason

# A canonical per-UniProt model id like "AF-P12345-F1" — accession starts with a
# letter. Numeric AB-INITIO models ("AF-0000000065764032") do not match, so this
# lets us prefer the canonical model over numeric ones during tie-breaks.
_CANONICAL_F1 = re.compile(r"^AF-[A-Za-z]\w*-F\d+$")

# AFDB reports sequence_identity as a float; treat >= this as an exact match.
_EXACT_IDENTITY = 0.9995


def _normalize_inputs(inputs) -> list[tuple[str, str, str | None]]:
    pairs: list[tuple[str, str, str | None]] = []
    for item in inputs:
        if isinstance(item, dict):
            pairs.append((item["id"], item["sequence"], item.get("accession")))
        else:
            id_, seq = item[0], item[1]
            acc = item[2] if len(item) > 2 else None
            pairs.append((id_, seq, acc))
    return pairs


def _chunked(items: list, size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _is_exact(summary: dict) -> bool:
    si = summary.get("sequence_identity")
    return isinstance(si, (int, float)) and si >= _EXACT_IDENTITY


def _rank_key(summary: dict, accession: str | None):
    """Preference order for full-length selection (lower sorts first):

    1. caller's accession (``AF-<accession>-F1``), then
    2. canonical ``AF-<acc>-F1`` models over numeric AB-INITIO models, then
    3. highest ``global_plddt`` (confidence_avg_local_score), then
    4. ``model_identifier`` lexicographically (deterministic final tie-break).
    """
    mid = summary.get("model_identifier") or ""
    acc_match = 0 if (accession and mid == f"AF-{accession}-F1") else 1
    canonical = 0 if _CANONICAL_F1.match(mid) else 1
    gp = summary.get("confidence_avg_local_score")
    gp = gp if isinstance(gp, (int, float)) else -1.0
    return (acc_match, canonical, -gp, mid)


class _Result(NamedTuple):
    summary_path: Path
    plddt_path: Path | None
    outcome: str  # "found" | "notfound" | "no_full_length" | "error"
    summary_data: dict | None
    plddt_values: list | None
    ambiguous: bool


def search_many(
    client,
    inputs,
    out_dir,
    *,
    concurrency: int = 6,
    rows: int = 10,
    plddt_first_n: int | None = None,
    full_length: bool = False,
) -> dict:
    """Query each queryable input's sequence concurrently, caching to disk.

    ``inputs`` is a list of ``(id, sequence)`` tuples or ``{"id":..., "sequence":...}``
    dicts. A dict may also carry an optional ``"accession"`` (or a 3rd tuple element)
    used only by ``full_length`` selection. Results are cached under ``out_dir``:

    * ``out_dir/summaries/{id}.json`` — a hit stores the AFDB summary document; a
      404 miss stores ``{"structures": []}`` so re-runs skip it. An existing file
      is left untouched (resumability).
    * ``out_dir/plddt/{id}.json`` (only when ``plddt_first_n`` is set) — the raw
      first-n per-residue pLDDT array for the selected structure (<= n values).

    With ``full_length=False`` (default) the pLDDT comes from ``structures[0]`` —
    whatever AFDB ranks first, which is **not** guaranteed to be the canonical
    single-chain model (it may be a longer multi-chain / AB-INITIO model).

    With ``full_length=True`` the selected structure must have
    ``sequence_identity == 1.0`` **and** a per-residue length equal to the query
    length, guarding against multi-chain / fragment models. Among such structures
    the caller's ``accession`` wins (``AF-<accession>-F1``); otherwise selection
    falls back to canonical ``-F1`` over numeric models, then highest
    ``global_plddt``, deterministically. A record whose hits include no
    exact-length match is counted under ``no_full_length`` (its summary is still
    written so re-runs resume) and no pLDDT is cached. A hit selected by fallback
    while more than one exact-sequence candidate existed is counted under
    ``ambiguous`` (still cached). Because length is only knowable from the
    confidence JSON, ``full_length=True`` fetches confidence even when
    ``plddt_first_n`` is None, and may fetch more than one per record.

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
        "no_full_length": 0,
        "ambiguous": 0,
        "errors": 0,
    }

    pending: list[tuple[Path, Path | None, str, str | None]] = []
    for id_, seq, accession in pairs:
        reason = filter_reason(seq)
        if reason is not None:
            counts[reason] += 1
            continue
        summary_path = summaries_dir / f"{id_}.json"
        if summary_path.exists():
            counts["skipped"] += 1
            continue
        plddt_path = (plddt_dir / f"{id_}.json") if plddt_first_n is not None else None
        pending.append((summary_path, plddt_path, seq, accession))

    def _select_full_length(structures: list, query_len: int, accession: str | None):
        """Return (summary, scores, ambiguous) for the best exact full-length hit,
        or None if no exact-sequence, exact-length structure exists. May raise
        httpx.HTTPError while fetching confidence.

        ``ambiguous`` is True when the caller's accession did not pin the choice and
        more than one exact-sequence structure also matched the query length.
        """
        candidates = [s for s in structures if _is_exact(s.get("summary", {}))]
        candidates.sort(key=lambda s: _rank_key(s["summary"], accession))

        # Fast path: the caller's accession pins the choice — if its model matches
        # the query length there is no ambiguity and no need to fetch the rest.
        if accession:
            target = f"AF-{accession}-F1"
            for s in candidates:
                if (s["summary"].get("model_identifier") or "") == target:
                    scores = client._fetch_confidence(
                        s["summary"]["model_url"]
                    ).get("confidenceScore", [])
                    if len(scores) == query_len:
                        return s["summary"], scores, False
                    break  # accession model present but wrong length; fall through

        # General path: fetch every exact-sequence candidate, keep those whose
        # length matches, then pick the best by preference. More than one length
        # match means the choice was not uniquely determined.
        matches = []
        for s in candidates:
            scores = client._fetch_confidence(
                s["summary"]["model_url"]
            ).get("confidenceScore", [])
            if len(scores) == query_len:
                matches.append((s["summary"], scores))
        if not matches:
            return None
        summary, scores = matches[0]
        return summary, scores, len(matches) > 1

    def _query(item: tuple[Path, Path | None, str, str | None]) -> _Result:
        summary_path, plddt_path, seq, accession = item
        try:
            data = client._fetch_summary(seq, rows)
        except httpx.HTTPError:
            return _Result(summary_path, plddt_path, "error", None, None, False)
        if data is None:
            return _Result(summary_path, plddt_path, "notfound", None, None, False)

        if full_length:
            structures = data.get("structures") or []
            if not structures:
                return _Result(summary_path, plddt_path, "notfound", None, None, False)
            try:
                selected = _select_full_length(structures, len(seq), accession)
            except httpx.HTTPError:
                return _Result(summary_path, plddt_path, "error", None, None, False)
            if selected is None:
                return _Result(summary_path, plddt_path, "no_full_length", data, None, False)
            _summary, scores, ambiguous = selected
            plddt_values = scores[:plddt_first_n] if plddt_first_n is not None else None
            return _Result(summary_path, plddt_path, "found", data, plddt_values, ambiguous)

        plddt_values = None
        if plddt_first_n is not None:
            structures = data.get("structures") or []
            if structures:
                try:
                    conf = client._fetch_confidence(structures[0]["summary"]["model_url"])
                except httpx.HTTPError:
                    return _Result(summary_path, plddt_path, "error", None, None, False)
                plddt_values = conf.get("confidenceScore", [])[:plddt_first_n]
        return _Result(summary_path, plddt_path, "found", data, plddt_values, False)

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
                    elif res.outcome == "no_full_length":
                        res.summary_path.write_text(json.dumps(res.summary_data))
                        counts["no_full_length"] += 1
                    else:
                        res.summary_path.write_text(json.dumps(res.summary_data))
                        counts["hits"] += 1
                        if res.ambiguous:
                            counts["ambiguous"] += 1
                        if res.plddt_values is not None and res.plddt_path is not None:
                            res.plddt_path.write_text(json.dumps(res.plddt_values))

    return {
        "total": len(pairs),
        "filtered": counts["internal_stop"] + counts["too_short"] + counts["nonstandard_aa"],
        **counts,
        "queried": counts["hits"]
        + counts["misses"]
        + counts["no_full_length"]
        + counts["errors"],
    }
