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

`select_group` applies three preference tiers — caller's accession, then monomers over
complexes, then canonical `-F1` over numeric ids — and returns **everything still tied**:

```python
from afdb_query import AlphaFold, select_group, mean_global_plddt, is_monomer

with AlphaFold() as af:
    doc = af.fetch_summary(sequence)           # None when AFDB has no entry
    group = select_group(doc["structures"]) if doc else []

    if not all(is_monomer(s) for s in group):
        ...                                   # your call: skip, or use it knowingly

    plddt = mean_global_plddt(group)           # average across the group, not one member
```

**No tier reads the confidence scores.** Ranking candidates by pLDDT and returning the
winner is not merely arbitrary, it is *biased*: the expected maximum of N draws rises
with N, so a comparison whose two arms match different numbers of candidates gets
different inflation on each side.

**There is no fourth tier, and that is deliberate.** After the three tiers a tie is
common — 25.5% of records on a real cache. Breaking it by any rule at all throws away
N−1 predictions of the same sequence for nothing. `select_group` keeps them and
`mean_global_plddt` averages across them: unbiased for the same reason an arbitrary
pick is (neither consults the value), with strictly lower variance because it uses
every prediction.

Tied candidates are near-always one protein reached through different UniProt entries —
identical-sequence orthologs, isoform duplicates, TrEMBL redundancy. On a real cache
4,288 of 4,289 tied sets were entirely full-length `-F1` models, so their per-residue
arrays are the same length and average elementwise too.

Do **not** take `group[0]`. The ordering is for reproducible output only; treating it
as "the" structure reintroduces exactly the arbitrary choice this API exists to avoid.

### Length is not visible to the tiers

`coverage == 1.0` means *the query is fully covered by the model*, **not** that the
model is the query's size — that is how an 860-residue multi-chain entry reports
coverage 1.0 against a 430-residue query. Filtering to monomers removes that case, but
an ortholog carrying the query sequence plus a few extra terminal residues passes every
visible tier, and then breaks two things: positional slicing (an offset from the
query's amino-acid length indexes the wrong residues) and the group average (its mean
spans residues your query does not contain).

The summary carries no residue count, so `filter_by_length` takes the lengths you
learned from fetching per-residue confidence:

```python
from afdb_query import filter_by_length

lengths = {s["model_identifier"]: len(fetched[s["model_identifier"]]) for s in group}
kept, dropped = filter_by_length(group, lengths, expected_length=len(sequence))
if dropped:
    log(f"{len(dropped)} of {len(group)} matched entries were the wrong length")
```

Unknown lengths are dropped, never assumed to conform, and dropped members are returned
rather than discarded so the loss is reportable.

`select_group` does not filter, either. "No usable structure" and "a group whose
average spans a complex" need different handling, and only the caller knows which its
analysis can tolerate — so test `is_monomer` on the members and decide.

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
inside a batch runner would hide it. Run `select_group` over the cached summaries and
average across what it returns.

## Not (yet) supported

- UniProt-accession lookup (sequence-only for now)
- PAE (Predicted Aligned Error)

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
