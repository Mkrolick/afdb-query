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


# -- full_length selection -------------------------------------------------
#
# SEQ_HIT is 20 residues. A "good" full-length model has 20 confidence scores;
# a multi-chain trap has 40. Each structure gets its own model/confidence URL.

def _struct(model_id, n_scores, *, identity=1.0, gp=90.0):
    cif = f"https://alphafold.ebi.ac.uk/files/{model_id}-model_v4.cif"
    return {
        "summary": {
            "model_identifier": model_id,
            "model_url": cif,
            "sequence_identity": identity,
            "confidence_avg_local_score": gp,
        }
    }, [float(i) for i in range(n_scores)]


def _route_conf(model_id, scores):
    url = f"https://alphafold.ebi.ac.uk/files/{model_id}-confidence_v4.json"
    respx.get(url).mock(
        return_value=httpx.Response(
            200, json={"residueNumber": list(range(1, len(scores) + 1)), "confidenceScore": scores}
        )
    )


@respx.mock
def test_full_length_picks_exact_length_over_multichain(tmp_path):
    # structures[0] is the 40-residue multi-chain trap; the exact-length model is second.
    trap, trap_scores = _struct("AF-0000000000000001", 40)        # numeric, multi-chain
    good, good_scores = _struct("AF-P00001-F1", len(SEQ_HIT))      # canonical, exact length
    respx.get(SUMMARY_URL).mock(
        return_value=httpx.Response(200, json={"structures": [trap, good]})
    )
    _route_conf("AF-0000000000000001", trap_scores)
    _route_conf("AF-P00001-F1", good_scores)

    with AlphaFold() as af:
        report = af.search_many(
            [{"id": "rec", "sequence": SEQ_HIT}], tmp_path, plddt_first_n=999, full_length=True
        )
    assert report["hits"] == 1
    assert report["no_full_length"] == 0
    assert report["ambiguous"] == 0
    # cached pLDDT is the exact-length model's, not the 40-residue trap's
    assert json.loads((tmp_path / "plddt" / "rec.json").read_text()) == good_scores


@respx.mock
def test_full_length_no_match_is_reported_and_resumes(tmp_path):
    trap, trap_scores = _struct("AF-0000000000000001", 40)
    respx.get(SUMMARY_URL).mock(
        return_value=httpx.Response(200, json={"structures": [trap]})
    )
    _route_conf("AF-0000000000000001", trap_scores)

    with AlphaFold() as af:
        report = af.search_many(
            [{"id": "rec", "sequence": SEQ_HIT}], tmp_path, plddt_first_n=999, full_length=True
        )
    assert report["hits"] == 0
    assert report["no_full_length"] == 1
    assert report["queried"] == 1
    # summary written (so re-runs resume), but no pLDDT cached
    assert (tmp_path / "summaries" / "rec.json").exists()
    assert not (tmp_path / "plddt" / "rec.json").exists()


@respx.mock
def test_full_length_accession_breaks_tie(tmp_path):
    # Two exact-length, identity-1.0 canonical models with different pLDDT.
    a, a_scores = _struct("AF-P00001-F1", len(SEQ_HIT), gp=90.0)
    b, b_scores = _struct("AF-P00002-F1", len(SEQ_HIT), gp=95.0)
    b_scores = [s + 50 for s in b_scores]  # make them distinguishable
    respx.get(SUMMARY_URL).mock(
        return_value=httpx.Response(200, json={"structures": [a, b]})
    )
    _route_conf("AF-P00001-F1", a_scores)
    _route_conf("AF-P00002-F1", b_scores)

    with AlphaFold() as af:
        report = af.search_many(
            [{"id": "rec", "sequence": SEQ_HIT, "accession": "P00001"}],
            tmp_path,
            plddt_first_n=999,
            full_length=True,
        )
    assert report["hits"] == 1
    assert report["ambiguous"] == 0  # accession disambiguated
    assert json.loads((tmp_path / "plddt" / "rec.json").read_text()) == a_scores


@respx.mock
def test_full_length_fallback_is_ambiguous_and_deterministic(tmp_path):
    # No accession: among two exact canonical models, pick higher global_plddt,
    # and flag the record ambiguous.
    a, a_scores = _struct("AF-P00001-F1", len(SEQ_HIT), gp=90.0)
    b, b_scores = _struct("AF-P00002-F1", len(SEQ_HIT), gp=95.0)
    b_scores = [s + 50 for s in b_scores]
    respx.get(SUMMARY_URL).mock(
        return_value=httpx.Response(200, json={"structures": [a, b]})
    )
    _route_conf("AF-P00001-F1", a_scores)
    _route_conf("AF-P00002-F1", b_scores)

    with AlphaFold() as af:
        report = af.search_many(
            [{"id": "rec", "sequence": SEQ_HIT}], tmp_path, plddt_first_n=999, full_length=True
        )
    assert report["hits"] == 1
    assert report["ambiguous"] == 1
    # higher global_plddt (b) wins the fallback
    assert json.loads((tmp_path / "plddt" / "rec.json").read_text()) == b_scores
