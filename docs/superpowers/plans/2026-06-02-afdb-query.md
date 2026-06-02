# afdb-query Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `afdb-query`, a small Python package for sequence-based programmatic access to the AlphaFold Protein Structure Database (AFDB), with a two-tier fetch (summary → per-residue pLDDT), ergonomic result objects, and a first-class resumable batch runner.

**Architecture:** An `AlphaFold` client wraps a shared `httpx.Client`. `search(sequence)` hits the AFDB `sequence/summary` endpoint (Tier 1) and returns `Structure` objects. `Structure.plddt()` lazily derives the confidence-JSON URL from the model URL and fetches per-residue pLDDT (Tier 2), exposing a `first(n)` accessor. `search_many(...)` runs the same lookups concurrently over many sequences with resumable on-disk caching and a counts report. UniProt lookup and PAE are designed-for-later, not built.

**Tech Stack:** Python 3.10+, `httpx` (runtime), `pytest` + `respx` (dev), `hatchling` build backend, `src/` layout.

**Spec:** `docs/superpowers/specs/2026-06-02-afdb-query-package-design.md`

**Commit convention:** Every commit message ends with the trailer:
```
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
```
(shown in each commit step below).

---

## File Structure

```
src/afdb_query/
  __init__.py     # public exports (Task 10)
  py.typed        # PEP 561 typed marker (Task 1)
  errors.py       # AFDBError, InvalidSequenceError (Task 2)
  sequences.py    # STANDARD_AA, MIN_LENGTH, filter_reason (Task 3)
  models.py       # confidence_url(), Plddt, Structure (Tasks 4, 5, 7)
  client.py       # AlphaFold: http session, _fetch_summary, _fetch_confidence, search, search_many (Tasks 6, 8)
  batch.py        # search_many implementation + report assembly (Task 9)
tests/
  test_smoke.py         # Tasks 1, 10
  test_sequences.py     # Task 3
  test_models.py        # Tasks 4, 5, 7
  test_client.py        # Tasks 6, 8
  test_batch.py         # Task 9
  test_integration.py   # Task 11 (network-gated)
pyproject.toml
README.md
```

---

### Task 1: Project scaffolding

**Files:**
- Create: `pyproject.toml`
- Create: `src/afdb_query/__init__.py` (empty for now)
- Create: `src/afdb_query/py.typed` (empty)
- Create: `tests/test_smoke.py`
- Create: `README.md` (stub; filled in Task 10)

- [ ] **Step 1: Create `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "afdb-query"
version = "0.1.0"
description = "Sequence-based programmatic access to the AlphaFold Protein Structure Database"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "MIT" }
dependencies = ["httpx>=0.27"]

[project.optional-dependencies]
dev = ["pytest>=8", "respx>=0.21"]

[tool.hatch.build.targets.wheel]
packages = ["src/afdb_query"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-m 'not integration'"
markers = [
  "integration: live network tests against AFDB (run with: pytest -m integration)",
]
```

- [ ] **Step 2: Create empty package files**

Create `src/afdb_query/__init__.py` with a single line:

```python
"""afdb-query: sequence-based programmatic access to the AlphaFold DB."""
```

Create `src/afdb_query/py.typed` as an empty file (PEP 561 marker so type checkers treat the package as typed).

- [ ] **Step 3: Create stub `README.md`**

```markdown
# afdb-query

Sequence-based programmatic access to the AlphaFold Protein Structure Database (AFDB).

Status: in development. See `docs/superpowers/specs/2026-06-02-afdb-query-package-design.md`.
```

- [ ] **Step 4: Write the smoke test**

Create `tests/test_smoke.py`:

```python
def test_package_imports():
    import afdb_query

    assert afdb_query.__doc__
```

- [ ] **Step 5: Create and activate a virtualenv, install in editable mode**

Run:
```bash
python3 -m venv .venv && . .venv/bin/activate && pip install -e '.[dev]'
```
Expected: ends with `Successfully installed ... afdb-query-0.1.0 httpx-... respx-... pytest-...`.

- [ ] **Step 6: Run the smoke test**

Run: `. .venv/bin/activate && pytest -q`
Expected: `1 passed`.

- [ ] **Step 7: Create `.gitignore`**

Create `.gitignore`:
```
.venv/
__pycache__/
*.egg-info/
.pytest_cache/
dist/
build/
```

- [ ] **Step 8: Commit**

```bash
git add pyproject.toml src/afdb_query/__init__.py src/afdb_query/py.typed tests/test_smoke.py README.md .gitignore
git commit -m "chore: scaffold afdb-query package" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Errors module

**Files:**
- Create: `src/afdb_query/errors.py`
- Test: `tests/test_models.py` (start the file here)

- [ ] **Step 1: Write the failing test**

Create `tests/test_models.py`:

```python
import pytest

from afdb_query.errors import AFDBError, InvalidSequenceError


def test_invalid_sequence_error_is_afdb_error():
    err = InvalidSequenceError("too_short")
    assert isinstance(err, AFDBError)
    assert err.reason == "too_short"
    assert "too_short" in str(err)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_models.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'afdb_query.errors'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/afdb_query/errors.py`:

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_models.py -q`
Expected: `1 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/errors.py tests/test_models.py
git commit -m "feat: add AFDBError and InvalidSequenceError" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Sequence validation (`filter_reason`)

**Files:**
- Create: `src/afdb_query/sequences.py`
- Test: `tests/test_sequences.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_sequences.py`:

```python
from afdb_query.sequences import MIN_LENGTH, STANDARD_AA, filter_reason

VALID = "ACDEFGHIKLMNPQRSTVWY"  # 20 residues, all standard


def test_valid_sequence_returns_none():
    assert filter_reason(VALID) is None


def test_internal_stop():
    assert filter_reason("ACD*EFGHIKLMNPQRSTVWY") == "internal_stop"


def test_too_short():
    assert filter_reason("ACDEF") == "too_short"


def test_nonstandard_aa():
    # 'B' and 'Z' are not standard amino acids
    assert filter_reason("ACDEFGHIKLMNPQRSTVWB") == "nonstandard_aa"


def test_internal_stop_takes_priority_over_length():
    # A short sequence containing '*' reports the stop first
    assert filter_reason("AC*") == "internal_stop"


def test_constants():
    assert MIN_LENGTH == 20
    assert "A" in STANDARD_AA and "B" not in STANDARD_AA
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_sequences.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'afdb_query.sequences'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/afdb_query/sequences.py`:

```python
"""Sequence validation for AFDB queries (ported from the original pipeline)."""

from __future__ import annotations

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH = 20


def filter_reason(seq: str) -> str | None:
    """Why a sequence cannot be queried against AFDB, or None if it is queryable.

    Checked in priority order: internal stop, length, non-standard residues.
    """
    if "*" in seq:
        return "internal_stop"
    if len(seq) < MIN_LENGTH:
        return "too_short"
    if not set(seq) <= STANDARD_AA:
        return "nonstandard_aa"
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_sequences.py -q`
Expected: `6 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/sequences.py tests/test_sequences.py
git commit -m "feat: add filter_reason sequence validation" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Confidence-URL derivation

**Files:**
- Create: `src/afdb_query/models.py`
- Test: `tests/test_models.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_models.py`:

```python
from afdb_query.models import confidence_url


def test_confidence_url_v1():
    assert (
        confidence_url("https://x/files/AF-1-model_v1.cif")
        == "https://x/files/AF-1-confidence_v1.json"
    )


def test_confidence_url_v6():
    assert (
        confidence_url("https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif")
        == "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-confidence_v6.json"
    )


def test_confidence_url_bcif():
    assert (
        confidence_url("https://x/files/AF-1-model_v4.bcif")
        == "https://x/files/AF-1-confidence_v4.json"
    )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_models.py -q`
Expected: FAIL with `ImportError: cannot import name 'confidence_url'` (or `ModuleNotFoundError`).

- [ ] **Step 3: Write minimal implementation**

Create `src/afdb_query/models.py`:

```python
"""Result objects and helpers for afdb-query."""

from __future__ import annotations


def confidence_url(model_url: str) -> str:
    """Derive the per-residue confidence-JSON URL from a model (CIF) URL.

    AFDB names files ``...-model_vN.cif`` and ``...-confidence_vN.json`` in the
    same directory, so the per-residue pLDDT URL is a pure string transform of
    the model URL (verified against the live API).
    """
    url = model_url.replace("-model_", "-confidence_")
    if url.endswith(".bcif"):
        return url[: -len(".bcif")] + ".json"
    if url.endswith(".cif"):
        return url[: -len(".cif")] + ".json"
    return url
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_models.py -q`
Expected: all passing (includes the Task 2 error test).

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/models.py tests/test_models.py
git commit -m "feat: derive confidence-JSON URL from model URL" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: `Plddt` result object

**Files:**
- Modify: `src/afdb_query/models.py`
- Test: `tests/test_models.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_models.py`:

```python
from afdb_query.models import Plddt


def _plddt(scores):
    return Plddt(scores=scores, residue_numbers=list(range(1, len(scores) + 1)), raw={})


def test_first_normal():
    assert _plddt([1.0, 2.0, 3.0, 4.0]).first(2) == [1.0, 2.0]


def test_first_more_than_len_returns_all():
    assert _plddt([1.0, 2.0]).first(10) == [1.0, 2.0]


def test_first_zero():
    assert _plddt([1.0, 2.0]).first(0) == []


def test_from_dict():
    p = Plddt.from_dict(
        {"confidenceScore": [5.0, 6.0], "residueNumber": [1, 2], "confidenceCategory": ["D", "D"]}
    )
    assert p.scores == [5.0, 6.0]
    assert p.residue_numbers == [1, 2]
    assert p.raw["confidenceCategory"] == ["D", "D"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_models.py -q`
Expected: FAIL with `ImportError: cannot import name 'Plddt'`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/afdb_query/models.py` (after `confidence_url`, and add the `dataclasses` import at the top of the file):

```python
from dataclasses import dataclass, field
```

```python
@dataclass(frozen=True)
class Plddt:
    """Per-residue pLDDT for one structure.

    ``scores`` and ``residue_numbers`` are parallel lists. ``raw`` is the full
    confidence-JSON document (escape hatch).
    """

    scores: list[float]
    residue_numbers: list[int]
    raw: dict = field(repr=False)

    @classmethod
    def from_dict(cls, data: dict) -> "Plddt":
        return cls(
            scores=data["confidenceScore"],
            residue_numbers=data["residueNumber"],
            raw=data,
        )

    def first(self, n: int) -> list[float]:
        """First ``n`` per-residue pLDDT values, or all of them if fewer than ``n``.

        Never pads and never raises on short structures: returns ``scores[:n]``.
        """
        return self.scores[:n]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_models.py -q`
Expected: all passing.

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/models.py tests/test_models.py
git commit -m "feat: add Plddt result object with first(n)" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: `AlphaFold` client core (session + low-level fetch)

**Files:**
- Create: `src/afdb_query/client.py`
- Test: `tests/test_client.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_client.py`:

```python
import httpx
import pytest
import respx

from afdb_query.client import AlphaFold

SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/sequence/summary"
VALID = "ACDEFGHIKLMNPQRSTVWY"


@respx.mock
def test_fetch_summary_success():
    respx.get(SUMMARY_URL).mock(
        return_value=httpx.Response(200, json={"entry": {}, "structures": []})
    )
    with AlphaFold() as af:
        data = af._fetch_summary(VALID)
    assert data == {"entry": {}, "structures": []}


@respx.mock
def test_fetch_summary_404_returns_none():
    respx.get(SUMMARY_URL).mock(return_value=httpx.Response(404))
    with AlphaFold() as af:
        assert af._fetch_summary(VALID) is None


@respx.mock
def test_fetch_summary_500_raises():
    respx.get(SUMMARY_URL).mock(return_value=httpx.Response(500))
    with AlphaFold() as af:
        with pytest.raises(httpx.HTTPStatusError):
            af._fetch_summary(VALID)


def test_fetch_summary_rows_guard():
    # API rejects rows <= 1; we fail fast with a clear ValueError before requesting.
    with AlphaFold() as af:
        with pytest.raises(ValueError):
            af._fetch_summary(VALID, rows=1)


@respx.mock
def test_fetch_confidence():
    model_url = "https://alphafold.ebi.ac.uk/files/AF-X-model_v4.cif"
    conf_url = "https://alphafold.ebi.ac.uk/files/AF-X-confidence_v4.json"
    respx.get(conf_url).mock(
        return_value=httpx.Response(
            200, json={"residueNumber": [1, 2], "confidenceScore": [10.0, 20.0]}
        )
    )
    with AlphaFold() as af:
        data = af._fetch_confidence(model_url)
    assert data["confidenceScore"] == [10.0, 20.0]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_client.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'afdb_query.client'`.

- [ ] **Step 3: Write minimal implementation**

Create `src/afdb_query/client.py`:

```python
"""The AlphaFold client: HTTP session and AFDB endpoint access."""

from __future__ import annotations

import httpx

from .models import confidence_url

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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_client.py -q`
Expected: `5 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/client.py tests/test_client.py
git commit -m "feat: add AlphaFold client core (summary + confidence fetch)" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: `Structure` result object (accessors + lazy `plddt()`)

**Files:**
- Modify: `src/afdb_query/models.py`
- Test: `tests/test_models.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_models.py` (add `httpx` and `respx` imports at the top of the file):

```python
import httpx
import respx

from afdb_query.client import AlphaFold
from afdb_query.models import Structure

SUMMARY = {
    "model_identifier": "AF-X",
    "model_url": "https://alphafold.ebi.ac.uk/files/AF-X-model_v1.cif",
    "confidence_avg_local_score": 91.65,
    "sequence_identity": 1.0,
    "coverage": 1.0,
    "entities": [
        {
            "identifier": "P12345",
            "identifier_category": "UNIPROT",
            "description": "Aspartate aminotransferase, mitochondrial",
        }
    ],
}


def test_structure_accessors():
    s = Structure(SUMMARY, None)
    assert s.model_identifier == "AF-X"
    assert s.model_url.endswith("AF-X-model_v1.cif")
    assert s.global_plddt == 91.65
    assert s.sequence_identity == 1.0
    assert s.coverage == 1.0
    assert s.uniprot_accession == "P12345"
    assert s.description == "Aspartate aminotransferase, mitochondrial"
    assert s.raw is SUMMARY


def test_structure_uniprot_missing_returns_none():
    s = Structure({"entities": [{"identifier": "x", "identifier_category": "PDB"}]}, None)
    assert s.uniprot_accession is None


@respx.mock
def test_structure_plddt_lazy_and_cached():
    conf_url = "https://alphafold.ebi.ac.uk/files/AF-X-confidence_v1.json"
    route = respx.get(conf_url).mock(
        return_value=httpx.Response(
            200, json={"residueNumber": [1, 2, 3], "confidenceScore": [10.0, 20.0, 30.0]}
        )
    )
    with AlphaFold() as af:
        s = Structure(SUMMARY, af)
        p1 = s.plddt()
        p2 = s.plddt()
    assert p1.scores == [10.0, 20.0, 30.0]
    assert p1.first(2) == [10.0, 20.0]
    assert p1 is p2  # cached on the instance
    assert route.call_count == 1  # fetched once
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_models.py -q`
Expected: FAIL with `ImportError: cannot import name 'Structure'`.

- [ ] **Step 3: Write minimal implementation**

Append to `src/afdb_query/models.py`:

```python
@dataclass(frozen=True)
class Structure:
    """One AFDB structure match for a queried sequence.

    Thin typed wrapper over the endpoint's ``summary`` dict. ``raw`` is the full
    summary (escape hatch). ``plddt()`` lazily fetches per-residue pLDDT.
    """

    raw: dict
    _client: "AlphaFold" = field(repr=False, compare=False)  # noqa: F821
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def model_identifier(self) -> str | None:
        return self.raw.get("model_identifier")

    @property
    def model_url(self) -> str | None:
        return self.raw.get("model_url")

    @property
    def global_plddt(self) -> float | None:
        return self.raw.get("confidence_avg_local_score")

    @property
    def sequence_identity(self) -> float | None:
        return self.raw.get("sequence_identity")

    @property
    def coverage(self) -> float | None:
        return self.raw.get("coverage")

    @property
    def uniprot_accession(self) -> str | None:
        for entity in self.raw.get("entities") or []:
            if entity.get("identifier_category") == "UNIPROT":
                return entity.get("identifier")
        return None

    @property
    def description(self) -> str | None:
        for entity in self.raw.get("entities") or []:
            if entity.get("description"):
                return entity["description"]
        return None

    def plddt(self) -> Plddt:
        """Tier 2: per-residue pLDDT for this structure (fetched once, then cached)."""
        if "plddt" not in self._cache:
            data = self._client._fetch_confidence(self.model_url)
            self._cache["plddt"] = Plddt.from_dict(data)
        return self._cache["plddt"]
```

Note: `_client` is typed as a forward reference (`"AlphaFold"`) and never imported here — `from __future__ import annotations` (already at the top) keeps this a string, avoiding a circular import with `client.py`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_models.py -q`
Expected: all passing.

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/models.py tests/test_models.py
git commit -m "feat: add Structure result object with lazy plddt()" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: `AlphaFold.search()` (Tier 1 → `list[Structure]`)

**Files:**
- Modify: `src/afdb_query/client.py`
- Test: `tests/test_client.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_client.py` (add the import at the top):

```python
from afdb_query.errors import InvalidSequenceError

STRUCT_SUMMARY = {
    "model_identifier": "AF-X",
    "model_url": "https://alphafold.ebi.ac.uk/files/AF-X-model_v1.cif",
    "confidence_avg_local_score": 91.65,
    "entities": [{"identifier": "P12345", "identifier_category": "UNIPROT"}],
}


@respx.mock
def test_search_success_returns_structures():
    respx.get(SUMMARY_URL).mock(
        return_value=httpx.Response(
            200, json={"entry": {}, "structures": [{"summary": STRUCT_SUMMARY}]}
        )
    )
    with AlphaFold() as af:
        hits = af.search(VALID)
    assert len(hits) == 1
    assert hits[0].global_plddt == 91.65
    assert hits[0].uniprot_accession == "P12345"


@respx.mock
def test_search_404_returns_empty_list():
    respx.get(SUMMARY_URL).mock(return_value=httpx.Response(404))
    with AlphaFold() as af:
        assert af.search(VALID) == []


def test_search_invalid_sequence_raises():
    with AlphaFold() as af:
        with pytest.raises(InvalidSequenceError) as excinfo:
            af.search("ACD*EFGHIKLMNPQRSTVWY")
    assert excinfo.value.reason == "internal_stop"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_client.py -q`
Expected: FAIL with `AttributeError: 'AlphaFold' object has no attribute 'search'`.

- [ ] **Step 3: Write minimal implementation**

In `src/afdb_query/client.py`, add these imports near the top (below the existing `import httpx`):

```python
from .errors import InvalidSequenceError
from .models import Structure, confidence_url
from .sequences import filter_reason
```

(Replace the existing `from .models import confidence_url` line with the combined import above.)

Then add the `search` method to the `AlphaFold` class (after `_fetch_confidence`):

```python
    # -- public API --------------------------------------------------------
    def search(self, sequence: str, rows: int = 10) -> list[Structure]:
        """Tier 1: find AFDB structures matching ``sequence``, best matches first.

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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_client.py -q`
Expected: all passing (8 tests).

- [ ] **Step 5: Commit**

```bash
git add src/afdb_query/client.py tests/test_client.py
git commit -m "feat: add AlphaFold.search (sequence -> Structures)" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Batch runner (`search_many`)

**Files:**
- Create: `src/afdb_query/batch.py`
- Modify: `src/afdb_query/client.py` (add `search_many` method delegating to `batch`)
- Test: `tests/test_batch.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_batch.py`:

```python
import json

import httpx
import pytest
import respx

from afdb_query.client import AlphaFold

SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/sequence/summary"
CONF_URL = "https://alphafold.ebi.ac.uk/files/AF-X-confidence_v1.json"

SEQ_HIT = "ACDEFGHIKLMNPQRSTVWY"   # valid -> hit
SEQ_MISS = "MNPQRSTVWYACDEFGHIKL"  # valid -> 404
SEQ_ERR = "GHIKLMNPQRSTVWYACDEFG"  # valid -> 500
SEQ_SHORT = "ACDEF"                # invalid -> too_short

HIT_SUMMARY = {
    "model_identifier": "AF-X",
    "model_url": "https://alphafold.ebi.ac.uk/files/AF-X-model_v1.cif",
    "confidence_avg_local_score": 91.65,
}


def _summary_side_effect(request):
    seq = request.url.params["id"]
    if seq == SEQ_HIT:
        return httpx.Response(200, json={"entry": {}, "structures": [{"summary": HIT_SUMMARY}]})
    if seq == SEQ_MISS:
        return httpx.Response(404)
    return httpx.Response(500)


@respx.mock
def test_search_many_counts_and_files(tmp_path):
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    inputs = [
        {"id": "hit", "sequence": SEQ_HIT},
        {"id": "miss", "sequence": SEQ_MISS},
        {"id": "err", "sequence": SEQ_ERR},
        {"id": "short", "sequence": SEQ_SHORT},
    ]
    with AlphaFold() as af:
        report = af.search_many(inputs, tmp_path, concurrency=2)

    assert report["total"] == 4
    assert report["hits"] == 1
    assert report["misses"] == 1
    assert report["errors"] == 1
    assert report["too_short"] == 1
    assert report["filtered"] == 1
    assert report["skipped"] == 0
    assert report["queried"] == 3

    assert (tmp_path / "summaries" / "hit.json").exists()
    assert json.loads((tmp_path / "summaries" / "miss.json").read_text()) == {"structures": []}
    assert not (tmp_path / "summaries" / "err.json").exists()  # errors are not saved


@respx.mock
def test_search_many_accepts_tuples(tmp_path):
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    with AlphaFold() as af:
        report = af.search_many([("hit", SEQ_HIT)], tmp_path)
    assert report["hits"] == 1


@respx.mock(assert_all_called=False)
def test_search_many_resumable_skips_existing(tmp_path):
    (tmp_path / "summaries").mkdir(parents=True)
    (tmp_path / "summaries" / "hit.json").write_text(
        json.dumps({"structures": [{"summary": HIT_SUMMARY}]})
    )
    route = respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    with AlphaFold() as af:
        report = af.search_many([{"id": "hit", "sequence": SEQ_HIT}], tmp_path)
    assert report["skipped"] == 1
    assert report["hits"] == 0
    assert route.call_count == 0  # existing file left untouched, no request made


@respx.mock
def test_search_many_plddt_first_n(tmp_path):
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    respx.get(CONF_URL).mock(
        return_value=httpx.Response(
            200, json={"residueNumber": [1, 2, 3, 4], "confidenceScore": [1.0, 2.0, 3.0, 4.0]}
        )
    )
    with AlphaFold() as af:
        report = af.search_many(
            [{"id": "hit", "sequence": SEQ_HIT}], tmp_path, plddt_first_n=2
        )
    assert report["hits"] == 1
    assert json.loads((tmp_path / "plddt" / "hit.json").read_text()) == [1.0, 2.0]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_batch.py -q`
Expected: FAIL with `AttributeError: 'AlphaFold' object has no attribute 'search_many'`.

- [ ] **Step 3: Write the batch implementation**

Create `src/afdb_query/batch.py`:

```python
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
```

- [ ] **Step 4: Wire the method onto the client**

In `src/afdb_query/client.py`, add this import near the other local imports:

```python
from .batch import search_many as _search_many
```

Then add the method to the `AlphaFold` class (after `search`):

```python
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
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_batch.py -q`
Expected: `4 passed`.

- [ ] **Step 6: Run the full suite**

Run: `pytest -q`
Expected: all passing, no failures.

- [ ] **Step 7: Commit**

```bash
git add src/afdb_query/batch.py src/afdb_query/client.py tests/test_batch.py
git commit -m "feat: add resumable concurrent search_many batch runner" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: Public exports + README

**Files:**
- Modify: `src/afdb_query/__init__.py`
- Modify: `tests/test_smoke.py`
- Modify: `README.md`

- [ ] **Step 1: Write the failing test**

Replace `tests/test_smoke.py` with:

```python
def test_public_exports():
    import afdb_query as m

    for name in [
        "AlphaFold",
        "Structure",
        "Plddt",
        "filter_reason",
        "confidence_url",
        "AFDBError",
        "InvalidSequenceError",
    ]:
        assert hasattr(m, name), f"missing export: {name}"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_smoke.py -q`
Expected: FAIL with `AssertionError: missing export: AlphaFold`.

- [ ] **Step 3: Write the exports**

Replace `src/afdb_query/__init__.py` with:

```python
"""afdb-query: sequence-based programmatic access to the AlphaFold DB."""

from .client import AlphaFold
from .errors import AFDBError, InvalidSequenceError
from .models import Plddt, Structure, confidence_url
from .sequences import filter_reason

__all__ = [
    "AlphaFold",
    "Structure",
    "Plddt",
    "filter_reason",
    "confidence_url",
    "AFDBError",
    "InvalidSequenceError",
]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_smoke.py -q`
Expected: `1 passed`.

- [ ] **Step 5: Write the README**

Replace `README.md` with:

````markdown
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
    hits = af.search(sequence)        # Tier 1: list[Structure], best matches first
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
````

- [ ] **Step 6: Run the full suite**

Run: `pytest -q`
Expected: all passing.

- [ ] **Step 7: Commit**

```bash
git add src/afdb_query/__init__.py tests/test_smoke.py README.md
git commit -m "feat: public exports and README" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 11: Network-gated integration tests

**Files:**
- Create: `tests/test_integration.py`

- [ ] **Step 1: Write the integration tests**

Create `tests/test_integration.py`. These hit the live AFDB API and are skipped by default (`addopts = -m 'not integration'`); run them with `pytest -m integration`.

```python
import pytest

from afdb_query import AlphaFold

# Rabbit GOT2 (UniProt P12345) — a stable AFDB entry.
GOT2 = (
    "MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYKRDTNSKKMNLGVGAYRDDNGKPYVLPSVRKAEAQ"
    "IAAKGLDKEYLPIGGLAEFCRASAELALGENSEVVKSGRFVTVQTISGTGALRIGASFLQRFFKFSRDVFLPKPSWGNHTPIFRDA"
    "GMQLQSYRYYDPKTCGFDFTGALEDISKIPEQSVLLLHACAHNPTGVDPRPEQWKEIATVVKKRNLFAFFDMAYQGFASGDGDKDA"
    "WAVRHFIEQGINVCLCQSYAKNMGLYGERVGAFTVICKDADEAKRVESQLKILIRPMYSNPPIHGARIASTILTSPDLRKQWLQEV"
    "KGMADRIIGMRTQLVSNLKKEGSTHSWQHITDQIGMFCFTGLKPEQVERLTKEFSIYMTKDGRISVAGVTSGNVGYLAHAIHQVTK"
)


@pytest.mark.integration
def test_live_search_returns_hits():
    with AlphaFold() as af:
        hits = af.search(GOT2)
    assert hits
    assert isinstance(hits[0].global_plddt, float)


@pytest.mark.integration
def test_live_plddt_first_n():
    with AlphaFold() as af:
        hits = af.search(GOT2)
        plddt = hits[0].plddt()
    assert len(plddt.scores) > 0
    first5 = plddt.first(5)
    assert len(first5) == min(5, len(plddt.scores))
    assert all(isinstance(x, float) for x in first5)


@pytest.mark.integration
def test_live_unknown_sequence_returns_empty():
    # A plausible-looking but synthetic sequence AFDB will not have.
    bogus = "ACDEFGHIKLMNPQRSTVWY" * 3
    with AlphaFold() as af:
        assert af.search(bogus) == []
```

- [ ] **Step 2: Verify they are deselected by default**

Run: `pytest -q`
Expected: integration tests are NOT run (deselected); all other tests pass. Output includes a `deselected` count.

- [ ] **Step 3: Run them against the live API**

Run: `pytest -m integration -q`
Expected: `3 passed` (requires network access to `alphafold.ebi.ac.uk`).

If `test_live_unknown_sequence_returns_empty` fails because the synthetic sequence
unexpectedly matches, replace `bogus` with a different synthetic standard-AA string
of length ≥ 20 and re-run.

- [ ] **Step 4: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: add network-gated integration tests" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] Run the full unit suite: `pytest -q` → all pass, integration deselected.
- [ ] Run integration once: `pytest -m integration -q` → all pass (network required).
- [ ] Build the package: `python -m build` (after `pip install build`) → produces
      `dist/afdb_query-0.1.0-py3-none-any.whl` and the sdist with no errors.
