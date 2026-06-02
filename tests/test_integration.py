import json

import httpx
import pytest

from afdb_query import AlphaFold, confidence_url

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


# -- raw AFDB URL contract tests -------------------------------------------
#
# These hit the live AFDB endpoints directly (not through the client) to catch
# upstream API drift: the URL scheme, status codes, and JSON field names the
# package depends on. If AFDB changes any of these, the package breaks — and
# these tests are what tells us before users do.


@pytest.mark.integration
def test_live_full_length_avoids_multichain_trap(tmp_path):
    # For P12345, AFDB ranks an 860-residue multi-chain model first; full_length
    # must instead cache the 430-residue canonical AF-P12345-F1.
    with AlphaFold() as af:
        report = af.search_many(
            [{"id": "got2", "sequence": GOT2, "accession": "P12345"}],
            tmp_path,
            plddt_first_n=9999999,
            full_length=True,
        )
    assert report["hits"] == 1
    assert report["no_full_length"] == 0
    cached = json.loads((tmp_path / "plddt" / "got2.json").read_text())
    assert len(cached) == len(GOT2)  # exact-length canonical model, not the 2x trap


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

    json_scores = httpx.get(
        confidence_url(model_url), timeout=60, follow_redirects=True
    ).json()["confidenceScore"]

    cif_text = httpx.get(model_url, timeout=60, follow_redirects=True).text
    cif_bfactors = _cif_ca_bfactors(cif_text)

    assert len(json_scores) == len(cif_bfactors) > 0
    assert all(abs(a - b) <= 0.01 for a, b in zip(json_scores, cif_bfactors))


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
