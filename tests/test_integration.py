import httpx
import pytest

from afdb_query import AlphaFold, confidence_url, is_monomer, mean_global_plddt, select_group

BASE_URL = "https://alphafold.ebi.ac.uk"
SUMMARY_URL = f"{BASE_URL}/api/sequence/summary"

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
def test_live_per_residue_plddt():
    with AlphaFold() as af:
        hits = af.search(GOT2)
        plddt = hits[0].plddt()
    assert len(plddt.scores) > 0
    assert len(plddt.scores) == len(plddt.residue_numbers)
    assert all(isinstance(x, float) for x in plddt.scores)


@pytest.mark.integration
def test_live_unknown_sequence_returns_empty():
    # A plausible-looking but synthetic sequence AFDB will not have.
    bogus = "ACDEFGHIKLMNPQRSTVWY" * 3
    with AlphaFold() as af:
        assert af.search(bogus) == []


# -- raw AFDB URL contract tests -------------------------------------------
#
# These hit the live AFDB endpoints directly (not through the client) to catch
# upstream API drift: the URL scheme, status codes, and JSON field names the
# package depends on. If AFDB changes any of these, the package breaks — and
# these tests are what tells us before users do.


@pytest.mark.integration
def test_live_select_avoids_multichain_trap():
    # For P12345, AFDB ranks an 860-residue multi-chain model first. `select` must
    # return the 430-residue canonical monomer instead -- and its per-residue array
    # must be the query's length, which is the property the global average depends on.
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": GOT2, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    group = select_group(resp.json()["structures"])
    assert group, "expected at least one candidate"
    assert all(is_monomer(s) for s in group), "selection must not return a complex"

    # every member must be a full-length model of the query, or averaging across the
    # group would mix different residue counts
    with AlphaFold() as af:
        for s in group:
            scores = af.fetch_confidence(s["model_url"])["confidenceScore"]
            assert len(scores) == len(GOT2), s["model_identifier"]


@pytest.mark.integration
def test_live_selection_ignores_pldd_ranking():
    """Selection must not depend on the confidence scores AFDB reports."""
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": GOT2, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    structures = resp.json()["structures"]
    baseline = [s["model_identifier"] for s in select_group(structures)]
    # Perturbing every score must not move the group.
    bumped = [
        {"summary": {**s["summary"], "confidence_avg_local_score": 99.99}} for s in structures
    ]
    assert [s["model_identifier"] for s in select_group(bumped)] == baseline


@pytest.mark.integration
def test_live_group_mean_is_not_the_max():
    """The group average must not coincide with the old max-over-candidates value."""
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": GOT2, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    group = select_group(resp.json()["structures"])
    scores = [s["confidence_avg_local_score"] for s in group]
    mean = mean_global_plddt(group)
    assert mean is not None
    assert min(scores) <= mean <= max(scores)
    if len(set(scores)) > 1:
        assert mean < max(scores)


@pytest.mark.integration
def test_raw_summary_endpoint_returns_expected_shape():
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": GOT2, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    assert resp.status_code == 200
    doc = resp.json()
    assert doc.get("structures"), "expected at least one structure"
    summary = doc["structures"][0]["summary"]
    # the exact fields the Structure wrapper reads
    for key in (
        "model_identifier",
        "model_url",
        "confidence_avg_local_score",
        "sequence_identity",
    ):
        assert key in summary, f"missing summary field: {key}"


@pytest.mark.integration
def test_raw_summary_endpoint_404s_for_unknown_sequence():
    bogus = "ACDEFGHIKLMNPQRSTVWY" * 3
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": bogus, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    # the client relies on a clean 404 to distinguish "no entry" from a real error
    assert resp.status_code == 404


@pytest.mark.integration
def test_raw_model_and_derived_confidence_urls_are_live():
    # The package derives the confidence-JSON URL from the model URL by a pure
    # string transform. Verify both the model file and that derived URL resolve
    # to live files on AFDB.
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": GOT2, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    model_url = resp.json()["structures"][0]["summary"]["model_url"]
    assert model_url.startswith(BASE_URL) and model_url.endswith(".cif")

    model = httpx.get(model_url, timeout=60, follow_redirects=True)
    assert model.status_code == 200
    assert model.text.lstrip().startswith("data_")  # mmCIF marker

    conf_url = confidence_url(model_url)
    assert conf_url.endswith(".json") and conf_url != model_url
    conf = httpx.get(conf_url, timeout=60, follow_redirects=True)
    assert conf.status_code == 200
    doc = conf.json()
    assert {"confidenceScore", "residueNumber"} <= doc.keys()
    assert len(doc["confidenceScore"]) == len(doc["residueNumber"]) > 0


@pytest.mark.integration
def test_raw_confidence_json_matches_cif_bfactor():
    # Ground-truth correctness: the per-residue pLDDT the package would return
    # (from the confidence JSON) equals the B-factor column of the deposited
    # mmCIF — two independently generated files served by AFDB.
    resp = httpx.get(
        SUMMARY_URL,
        params={"id": GOT2, "type": "sequence", "rows": 10},
        headers={"Accept": "application/json"},
        timeout=30,
    )
    summary = next(
        s["summary"]
        for s in resp.json()["structures"]
        if s["summary"]["model_identifier"] == "AF-P12345-F1"
    )
    model_url = summary["model_url"]

    json_scores = httpx.get(confidence_url(model_url), timeout=60, follow_redirects=True).json()[
        "confidenceScore"
    ]

    cif_text = httpx.get(model_url, timeout=60, follow_redirects=True).text
    cif_bfactors = _cif_ca_bfactors(cif_text)

    assert len(json_scores) == len(cif_bfactors) > 0
    assert all(abs(a - b) <= 0.01 for a, b in zip(json_scores, cif_bfactors, strict=True))


def _cif_ca_bfactors(text):
    """Per-residue pLDDT from an mmCIF: the B-factor of each CA atom, parsing the
    _atom_site loop generically (no external CIF dependency)."""
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "loop_":
            cols, j = [], i + 1
            while j < n and lines[j].lstrip().startswith("_"):
                cols.append(lines[j].strip())
                j += 1
            if any(c.startswith("_atom_site.") for c in cols):
                idx = {c: k for k, c in enumerate(cols)}
                ca = idx["_atom_site.label_atom_id"]
                cb = idx["_atom_site.B_iso_or_equiv"]
                out = []
                while j < n:
                    s = lines[j].strip()
                    if s in ("#", "loop_", "") or s.startswith("_") or s.startswith("data_"):
                        break
                    f = lines[j].split()
                    if f and f[0] in ("ATOM", "HETATM") and f[ca] == "CA":
                        out.append(float(f[cb]))
                    j += 1
                return out
            i = j
        else:
            i += 1
    return []
