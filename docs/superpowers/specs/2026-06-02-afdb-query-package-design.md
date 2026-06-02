# afdb-query — Design

**Date:** 2026-06-02
**Status:** Approved — endpoints validated against the live API 2026-06-02

## Purpose

A small, reusable Python package for programmatic access to the
[AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) (AFDB) **by raw
amino-acid sequence**. It exists to make two things easy that are currently painful:

1. **Query a protein by sequence** against AFDB and get back the matching structures
   with their metadata and global pLDDT.
2. **Retrieve per-residue pLDDT** for a matched structure — in particular, the
   *first n* per-residue pLDDT values — without hand-rolling URL derivation and JSON
   fetching each time.

It also ships a **first-class batch runner** (concurrency + resumable on-disk caching + a counts report) so the same operations run cleanly over many sequences in a data pipeline.

### Non-goals (v1)

- **UniProt-accession lookup** — designed for, not built. Sequence-only for now.
- **PAE (Predicted Aligned Error)** — reserved in the design (a future `Structure.pae()`
  slots in beside `plddt()`), but not implemented.
- **No statistics helpers.** The package returns raw values (e.g. the first-n pLDDT
  array). All downstream math (means, thresholds, aggregation) stays with the caller.
  There is explicitly **no `mean()`** function.
- Running AlphaFold predictions, parsing CIF/PDB geometry, structure visualization.

## Package identity

- **Distribution name (PyPI):** `afdb-query` (verified available 2026-06-02; `afquery`
  is taken).
- **Import package:** `afdb_query`.
- **Python:** 3.10+ (uses `X | None` syntax).
- **Runtime dependency:** `httpx` only.
- **Dev dependencies:** `pytest`, `respx`.
- **Layout:** `src/` layout, `pyproject.toml`.

## AFDB API facts (verified 2026-06-02)

These were confirmed against the live API and drive the design.

### Tier 1 — sequence summary

`GET https://alphafold.ebi.ac.uk/api/sequence/summary`
with params `id=<sequence>`, `type=sequence`, `rows=<n>` and header `Accept: application/json`.

- **`rows` must be > 1** (the API rejects `rows=1` with a 422 validation error).
- 404 → the sequence has no AFDB entry (a clean "not found", distinct from a real error).
- Response shape (3D-Beacons schema):
  ```json
  {
    "entry": {"sequence": "...", "checksum": "...", "checksum_type": "..."},
    "structures": [
      {"summary": {
        "model_identifier": "AF-0000000065764032",
        "model_url": "https://alphafold.ebi.ac.uk/files/AF-0000000065764032-model_v1.cif",
        "model_category": "AB-INITIO",
        "confidence_type": "pLDDT",
        "confidence_avg_local_score": 91.65,
        "sequence_identity": 1.0,
        "coverage": 1.0,
        "entities": [{
          "identifier": "P12345", "identifier_category": "UNIPROT",
          "description": "Aspartate aminotransferase, mitochondrial",
          "chain_ids": ["A"]
        }]
      }}
    ]
  }
  ```
- `confidence_avg_local_score` is the **global mean pLDDT** — available cheaply at Tier 1.
- The summary does **not** include a per-residue confidence URL.
- **Result order:** the API ranks matches by sequence identity but does **not** place the
  canonical `AF-<accession>-F1` model first (see the ordering caveat under *Public API*).

### Tier 2 — per-residue pLDDT

Derive the confidence-JSON URL directly from `model_url` by string substitution:
`-model_vN.cif` → `-confidence_vN.json` (verified: derived URL returns HTTP 200).

`GET <confidence_url>` returns:
```json
{
  "residueNumber":      [1, 2, 3, ...],
  "confidenceScore":    [18.17, 19.14, 21.67, ...],
  "confidenceCategory": ["D", "D", ...],
  "chains":             [...]
}
```
`confidenceScore` is the per-residue pLDDT array. **No CIF parser dependency is needed.**

- **`chains` is optional** (verified 2026-06-02): the multi-chain seq-search hit
  `AF-0000000065764032-confidence_v1.json` includes it; the canonical monomer
  `AF-P12345-F1-confidence_v6.json` omits it. The package reads only `confidenceScore` and
  `residueNumber` (stashing the full document in `Plddt.raw`), so this is non-breaking — but
  callers must not assume `chains` exists.
- Every `model_url` observed from the summary endpoint ends in `.cif` (v1 and v6). The
  URL-derivation also handles `.bcif` defensively, but that branch is **not exercised by live
  data** — it is speculative and may be dropped.

## Public API

```python
from afdb_query import AlphaFold

with AlphaFold() as af:                 # wraps a shared httpx.Client; context-managed
    hits = af.search(sequence)          # Tier 1 -> list[Structure], in AFDB's returned order
    s = hits[0]

    s.global_plddt        # float  — confidence_avg_local_score (cheap, from Tier 1)
    s.sequence_identity   # float  — 1.0 == exact match, <1.0 == near hit
    s.coverage            # float
    s.uniprot_accession   # str | None — from entities[] (identifier_category == UNIPROT)
    s.description         # str | None — from entities[]
    s.model_identifier    # str
    s.model_url           # str    — CIF url
    s.raw                 # dict   — full summary dict, escape hatch

    p = s.plddt()         # Tier 2: lazily fetches + caches the confidence JSON
    p.scores              # list[float] — full per-residue pLDDT array
    p.residue_numbers     # list[int]   — parallel to scores
    p.first(50)           # list[float] — first 50 values, or ALL if len < 50
```

> **Result-ordering caveat (verified 2026-06-02).** `search()` preserves the order AFDB
> returns. Matches are ranked by sequence identity, but `hits[0]` is **not** guaranteed to be
> the canonical `AF-<accession>-F1` monomer prediction. For the GOT2/P12345 sequence, `hits[0]`
> is the multi-chain model `AF-0000000065764032` (860 residues, mean pLDDT 91.65) while the
> canonical monomer `AF-P12345-F1` (430 residues, 94.12) is `hits[1]`. So `plddt()`/`first(n)`
> operate on whichever model the caller selects from `hits`; for a multi-chain model the
> per-residue array spans all chains and `residue_numbers` restart at each chain boundary.
> **v1 keeps `search()` order as-is** — callers who specifically want the canonical entry
> should pick the `hits[i]` whose `model_identifier` matches `AF-<accession>-F1`. Selecting the
> canonical entry automatically (e.g. a `canonical_only=` flag, or a `Structure.is_canonical`
> helper) is a possible future enhancement.

### Components

- **`AlphaFold`** (`client.py`) — holds the `httpx.Client`, timeout/retry config, base
  URL. Context-manager (`with AlphaFold() as af:`) and an explicit `close()`. Methods:
  - `search(sequence: str, rows: int = 10) -> list[Structure]`
  - `search_many(...)` (see batch section)
  - Construction options: `timeout`, `base_url`, `max_retries` (sensible defaults).
- **`Structure`** (`models.py`) — frozen dataclass wrapping one `summary` dict. Typed
  accessors as above; `raw` exposes the dict. Holds a reference to its `AlphaFold` client
  so `plddt()` can fetch. `plddt()` fetches lazily and **caches on the instance** (repeat
  calls do not re-hit the network). `pae()` is reserved (not implemented in v1).
- **`Plddt`** (`models.py`) — frozen dataclass holding `scores`, `residue_numbers`
  (and the raw dict). `first(n)` returns `scores[:n]` — i.e. `min(n, len(scores))`
  values, never padding, never erroring on short structures. **No `mean()`.**
- **`sequences.py`** — `filter_reason(seq) -> str | None` and constants
  (`STANDARD_AA`, min length). Returns one of `"internal_stop"`, `"too_short"`,
  `"nonstandard_aa"`, or `None` when queryable. Ported from existing pipeline code.
- **`errors.py`** — `AFDBError` (base) and `InvalidSequenceError`.

## Data flow

1. **`search(seq)`**: validate via `filter_reason`. If invalid → raise
   `InvalidSequenceError`. Otherwise `GET sequence/summary`. On 404 → return `[]`.
   On other HTTP error → raise (propagates). On success → wrap each
   `structures[i].summary` in a `Structure`.
2. **`structure.plddt()`**: derive confidence URL from `model_url`, `GET` it, build a
   `Plddt`, cache it on the `Structure`.

## Batch runner (first-class default)

The reusable core of the existing `run_lookup`, generalized so it is not tied to any one
pipeline's record shape.

```python
report = af.search_many(
    inputs,                 # list of (id, sequence) or {"id":..., "sequence":...}
    out_dir,                # disk cache root (Path)
    concurrency=6,
    rows=10,
    plddt_first_n=None,     # if set (e.g. 50), also fetch Tier-2 first-n pLDDT per hit
)
```

- **Caller-supplied generic `id`** keys every cached file and result, mapping back to the
  caller's own records (variants or anything else). The package does not impose a naming
  scheme beyond using the `id`.
- **Resumable on-disk cache:**
  - `out_dir/summaries/{id}.json` — a hit stores the AFDB summary; a 404 miss stores the
    marker `{"structures": []}` so re-runs skip it; an existing file is left untouched.
  - When `plddt_first_n` is set: `out_dir/plddt/{id}.json` stores the raw first-n pLDDT
    array (`<= n` values, all if the structure is shorter). Off by default → stays cheap.
    Uses the first/best structure of each hit.
- **Counts report** (returned `dict`): `total`, per-reason filter breakdown
  (`internal_stop`, `too_short`, `nonstandard_aa`), aggregate `filtered`, `skipped`,
  `hits`, `misses`, `errors`, and `queried` (= hits + misses + errors).
- **Error semantics:** a real per-query HTTP error is **counted and not saved**, so it is
  retried on the next run. Invalid sequences are counted by reason and skipped (never
  queried).
- **Concurrency:** `ThreadPoolExecutor(max_workers=concurrency)` over the shared,
  thread-safe `httpx.Client`, processed in chunks (mirrors the existing implementation).

## Error handling summary

| Situation                         | Single-protein (`search`)        | Batch (`search_many`)              |
|-----------------------------------|----------------------------------|------------------------------------|
| Sequence fails `filter_reason`    | raise `InvalidSequenceError`     | counted by reason, skipped         |
| AFDB has no entry (HTTP 404)      | return `[]`                      | write `{"structures": []}`, miss++ |
| Real HTTP / transport error       | propagate (raise)                | counted as error, **not saved**    |

## Testing strategy

- **Unit tests** with mocked HTTP (`respx`):
  - `filter_reason`: internal stop, too short (`< 20`), non-standard AA, valid.
  - Confidence-URL derivation from `model_url`.
  - `Plddt.first(n)`: normal, `n` larger than length (returns all), `n == 0`.
  - `search`: 404 → `[]`; success → typed `Structure` fields; invalid seq → raises.
  - Batch: resumability (existing file skipped), 404 marker written, counts correctness,
    real error counted and not saved, `plddt_first_n` writes first-n file.
- **Integration tests** against the live API, marked and skipped by default (network-gated).

## Module layout

```
src/afdb_query/
  __init__.py     # exports AlphaFold, Structure, Plddt, errors
  client.py       # AlphaFold (search, search_many, http/session mgmt)
  models.py       # Structure, Plddt
  sequences.py    # filter_reason, constants
  batch.py        # search_many implementation + report assembly
  errors.py       # AFDBError, InvalidSequenceError
tests/
  ...
pyproject.toml
README.md
```
