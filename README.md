# afdb-query

Sequence-based programmatic access to the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) (AFDB). Query a protein by its amino-acid sequence, then pull per-residue pLDDT without hand-rolling URL derivation and JSON fetching.

## Install

```bash
pip install afdb-query
```

## Quickstart

```python
from afdb_query import AlphaFold

with AlphaFold() as af:
    hits = af.search(sequence)        # Tier 1: list[Structure], in AFDB's returned order
    s = hits[0]

    s.global_plddt        # mean pLDDT AFDB reports for the model (from the summary)
    s.sequence_identity   # 1.0 == exact match
    s.oligomeric_state    # "MONOMER", "HOMODIMER", ...
    s.uniprot_accession   # e.g. "P12345", or None

    p = s.plddt()         # Tier 2: per-residue pLDDT (fetched once, then cached)
    p.scores              # full per-residue list[float]
    p.residue_numbers     # parallel residue numbers — do not assume 1..N
```

`search` raises `InvalidSequenceError` for sequences that cannot be queried
(internal stop `*`, shorter than 20 residues, or non-standard amino acids), and
returns `[]` when AFDB has no entry for a valid sequence. Every HTTP or transport
failure is raised as `AFDBHTTPError`, so catching `AFDBError` is enough — you never
need to import `httpx`.

## Choosing a structure

A sequence query can match several structures, and **they are not interchangeable**:

- `hits[0]` is *not* guaranteed to be the canonical `AF-<accession>-F1` model. For
  some sequences a multi-chain or numeric AB-INITIO model ranks first.
- `global_plddt` (`confidence_avg_local_score`) is averaged over the **whole
  deposited model**. For a HOMODIMER or HETERODIMER that spans every chain, not just
  the one matching your query, so it is not comparable with a monomer's.

`select` applies a deterministic preference order — caller's accession, then
monomers over complexes, then canonical `-F1` over numeric ids, then
`model_identifier` lexicographically:

```python
from afdb_query import AlphaFold, select, is_monomer

with AlphaFold() as af:
    doc = af.fetch_summary(sequence)
    chosen = select(doc["structures"])          # bare summary dict, or None
    if not is_monomer(chosen):
        ...                                     # your call: skip it, or use it knowingly
    scores = af.fetch_confidence(chosen["model_url"])["confidenceScore"]
```

**Selection never reads the confidence scores.** Ranking candidates by pLDDT and
returning the winner would make this library pick the best-scoring structure, which
silently biases any downstream comparison of one protein's pLDDT against another's.
Where identity and provenance cannot decide, `select` falls back to a lexicographic
tie-break rather than reaching for the score.

`select` does not filter, either. "No usable structure" and "a structure whose
average spans a complex" need different handling, and only the caller knows which
its analysis can tolerate — so test `is_monomer` on the result and decide.

## Batch lookups

`search_many` runs many sequences concurrently with resumable on-disk caching:

```python
report = af.search_many(
    [{"id": "rec1", "sequence": seq1}, {"id": "rec2", "sequence": seq2}],
    out_dir="afdb_cache",
    concurrency=6,
)
```

```python
{
  "total":    2,
  "skipped":  0,                                                   # already cached
  "filtered": {"internal_stop": 0, "too_short": 0, "nonstandard_aa": 0, "total": 0},
  "queried":  {"hits": 2, "misses": 0, "errors": 0, "total": 2},
}
```

`total == skipped + filtered["total"] + queried["total"]`; no count appears twice.

- You supply a generic `id` per sequence; it keys the cache file and maps back to
  your own records.
- `out_dir/summaries/{id}.json` stores each hit (a 404 miss stores
  `{"structures": []}`); existing files are left untouched, so re-runs resume.
- Real HTTP errors are counted but not saved, so they retry on the next run.

`search_many` fetches **summaries only**. It does not choose a structure and it does
not fetch per-residue confidence: which of several exact-sequence matches answers
your question is a property of your analysis, not of AFDB, and making that choice
inside a batch runner would hide it. Iterate the cached summaries with `select` and
`fetch_confidence` when you need per-residue data.

## Not (yet) supported

- UniProt-accession lookup (sequence-only for now)
- PAE (Predicted Aligned Error)
- No statistics helpers — the package returns raw values; downstream math is yours.

## Development

```bash
pip install -e ".[dev]"
ruff check . && ruff format --check .
pytest                    # unit suite; 100% branch coverage is enforced
pytest -m integration     # live AFDB tests (network required)
```

The integration suite includes a ground-truth check that the per-residue pLDDT in
AFDB's confidence JSON matches the B-factor column of the deposited mmCIF — two
independently generated files — so upstream API drift surfaces here rather than in
your results.
