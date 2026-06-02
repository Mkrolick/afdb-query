# afdb-query

Sequence-based programmatic access to the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) (AFDB). Query a protein by its amino-acid sequence, then pull per-residue pLDDT — including "the first n values" — without hand-rolling URL derivation and JSON fetching.

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

    s.global_plddt        # mean pLDDT for the model (cheap, from the summary)
    s.sequence_identity   # 1.0 == exact match, < 1.0 == near hit
    s.uniprot_accession   # e.g. "P12345", or None

    p = s.plddt()         # Tier 2: per-residue pLDDT (fetched once, then cached)
    p.scores              # full per-residue list[float]
    p.first(50)           # first 50 values — or all of them if the model is shorter
```

`search` raises `InvalidSequenceError` for sequences that cannot be queried
(internal stop `*`, shorter than 20 residues, or non-standard amino acids), and
returns `[]` when AFDB has no entry for a valid sequence.

Results come back in AFDB's returned order (ranked by sequence identity). Note that
`hits[0]` is **not** guaranteed to be the canonical `AF-<accession>-F1` model — for
some sequences a multi-chain or AB-INITIO model ranks first — so pick the hit whose
`model_identifier` you want if you need a specific entry.

## Batch lookups

`search_many` runs many sequences concurrently with resumable on-disk caching:

```python
report = af.search_many(
    [{"id": "rec1", "sequence": seq1}, {"id": "rec2", "sequence": seq2}],
    out_dir="afdb_cache",
    concurrency=6,
    plddt_first_n=50,   # optional: also save the first 50 per-residue pLDDT per hit
)
# report -> {"total":..., "hits":..., "misses":..., "errors":..., "skipped":..., ...}
```

- You supply a generic `id` per sequence; it keys the cache file and maps back to
  your own records.
- `out_dir/summaries/{id}.json` stores each hit (a 404 miss stores
  `{"structures": []}`); existing files are left untouched, so re-runs resume.
- With `plddt_first_n` set, `out_dir/plddt/{id}.json` stores the raw first-n
  per-residue pLDDT array for the best structure.
- Real HTTP errors are counted but not saved, so they retry on the next run.

  Note: resumability keys on the summary file. If you run once without
  `plddt_first_n` and again with it, already-cached records are skipped and their
  pLDDT is not back-filled.

## Not (yet) supported

- UniProt-accession lookup (sequence-only for now)
- PAE (Predicted Aligned Error)
- No statistics helpers — the package returns raw values; downstream math is yours.
