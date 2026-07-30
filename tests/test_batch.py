import json

import httpx
import respx

from afdb_query.batch import fetch_plddt_many, load_plddt, plddt_path
from afdb_query.client import AlphaFold

SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/sequence/summary"

SEQ_HIT = "ACDEFGHIKLMNPQRSTVWY"  # valid -> hit
SEQ_MISS = "MNPQRSTVWYACDEFGHIKL"  # valid -> 404
SEQ_ERR = "GHIKLMNPQRSTVWYACDEFG"  # valid -> 500
SEQ_SHORT = "ACDEF"  # invalid -> too_short

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
    assert report["skipped"] == 0
    assert report["filtered"] == {
        "internal_stop": 0,
        "too_short": 1,
        "nonstandard_aa": 0,
        "total": 1,
    }
    assert report["queried"] == {"hits": 1, "misses": 1, "errors": 1, "total": 3}

    assert (tmp_path / "summaries" / "hit.json").exists()
    assert json.loads((tmp_path / "summaries" / "miss.json").read_text()) == {"structures": []}
    assert not (tmp_path / "summaries" / "err.json").exists()  # errors are not saved


@respx.mock
def test_search_many_report_counts_are_disjoint(tmp_path):
    """total is exactly skipped + filtered + queried -- no count is reported twice."""
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    (tmp_path / "summaries").mkdir(parents=True)
    (tmp_path / "summaries" / "cached.json").write_text(json.dumps({"structures": []}))
    inputs = [
        {"id": "cached", "sequence": SEQ_HIT},
        {"id": "hit", "sequence": SEQ_HIT},
        {"id": "miss", "sequence": SEQ_MISS},
        {"id": "short", "sequence": SEQ_SHORT},
        {"id": "stop", "sequence": "ACD*EFGHIKLMNPQRSTVWY"},
    ]
    with AlphaFold() as af:
        report = af.search_many(inputs, tmp_path)

    assert report["total"] == 5
    assert (
        report["skipped"] + report["filtered"]["total"] + report["queried"]["total"]
        == report["total"]
    )
    assert report["filtered"]["internal_stop"] == 1
    assert report["filtered"]["too_short"] == 1


@respx.mock
def test_search_many_accepts_tuples(tmp_path):
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    with AlphaFold() as af:
        report = af.search_many([("hit", SEQ_HIT)], tmp_path)
    assert report["queried"]["hits"] == 1


@respx.mock
def test_search_many_nonstandard_aa_is_filtered(tmp_path):
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    with AlphaFold() as af:
        report = af.search_many([("bad", "ACDEFGHIKLMNPQRSTVWB")], tmp_path)
    assert report["filtered"]["nonstandard_aa"] == 1
    assert report["queried"]["total"] == 0


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
    assert report["queried"]["hits"] == 0
    assert route.call_count == 0  # existing file left untouched, no request made


@respx.mock
def test_search_many_chunks_larger_than_concurrency(tmp_path):
    """The chunking loop runs more than once when pending exceeds concurrency * 50."""
    respx.get(SUMMARY_URL).mock(side_effect=_summary_side_effect)
    inputs = [{"id": f"r{i}", "sequence": SEQ_HIT} for i in range(3)]
    with AlphaFold() as af:
        report = af.search_many(inputs, tmp_path, concurrency=1)
    assert report["queried"]["hits"] == 3


def test_search_many_empty_input_makes_no_directory(tmp_path):
    with AlphaFold() as af:
        report = af.search_many([], tmp_path)
    assert report["total"] == 0
    assert report["queried"]["total"] == 0
    assert not (tmp_path / "summaries").exists()


# -- per-residue pass ------------------------------------------------------

CONF_URL = "https://alphafold.ebi.ac.uk/files/AF-X-confidence_v1.json"
MODEL_URL = "https://alphafold.ebi.ac.uk/files/AF-X-model_v1.cif"
CONF_DOC = {"residueNumber": [1, 2, 3], "confidenceScore": [10.0, 20.0, 30.0]}


@respx.mock
def test_fetch_plddt_many_caches_full_array(tmp_path):
    respx.get(CONF_URL).mock(return_value=httpx.Response(200, json=CONF_DOC))
    with AlphaFold() as af:
        report = fetch_plddt_many(af, [{"id": "rec", "model_url": MODEL_URL}], tmp_path)
    assert report == {"total": 1, "skipped": 0, "queried": {"hits": 1, "errors": 0, "total": 1}}
    got = load_plddt(tmp_path, "rec")
    assert got.scores == [10.0, 20.0, 30.0]  # whole array, never a prefix
    assert got.residue_numbers == [1, 2, 3]
    assert got.raw["model_url"] == MODEL_URL  # provenance kept


@respx.mock
def test_fetch_plddt_many_accepts_tuples(tmp_path):
    respx.get(CONF_URL).mock(return_value=httpx.Response(200, json=CONF_DOC))
    with AlphaFold() as af:
        report = fetch_plddt_many(af, [("rec", MODEL_URL)], tmp_path)
    assert report["queried"]["hits"] == 1


@respx.mock
def test_fetch_plddt_many_resumes_on_its_own_artifact(tmp_path):
    """Resumability keys on the plddt file, NOT the summary -- so a summary-only
    run followed by a per-residue run back-fills instead of skipping."""
    (tmp_path / "summaries").mkdir(parents=True)
    (tmp_path / "summaries" / "rec.json").write_text(json.dumps({"structures": []}))
    route = respx.get(CONF_URL).mock(return_value=httpx.Response(200, json=CONF_DOC))
    with AlphaFold() as af:
        report = fetch_plddt_many(af, [("rec", MODEL_URL)], tmp_path)
    assert report["skipped"] == 0  # the cached SUMMARY does not block it
    assert report["queried"]["hits"] == 1
    assert route.call_count == 1


@respx.mock(assert_all_called=False)
def test_fetch_plddt_many_skips_already_fetched(tmp_path):
    (tmp_path / "plddt").mkdir(parents=True)
    (tmp_path / "plddt" / "rec.json").write_text(json.dumps(CONF_DOC))
    route = respx.get(CONF_URL).mock(return_value=httpx.Response(200, json=CONF_DOC))
    with AlphaFold() as af:
        report = fetch_plddt_many(af, [("rec", MODEL_URL)], tmp_path)
    assert report["skipped"] == 1
    assert route.call_count == 0


@respx.mock
def test_fetch_plddt_many_error_is_counted_not_written(tmp_path):
    respx.get(CONF_URL).mock(return_value=httpx.Response(500))
    with AlphaFold() as af:
        report = fetch_plddt_many(af, [("rec", MODEL_URL)], tmp_path)
    assert report["queried"]["errors"] == 1
    assert load_plddt(tmp_path, "rec") is None  # retried next run


@respx.mock
def test_fetch_plddt_many_missing_fields_become_empty(tmp_path):
    respx.get(CONF_URL).mock(return_value=httpx.Response(200, json={}))
    with AlphaFold() as af:
        fetch_plddt_many(af, [("rec", MODEL_URL)], tmp_path)
    got = load_plddt(tmp_path, "rec")
    assert got.scores == [] and got.residue_numbers == []
    assert got.mean() is None


def test_fetch_plddt_many_empty_input_makes_no_directory(tmp_path):
    with AlphaFold() as af:
        report = fetch_plddt_many(af, [], tmp_path)
    assert report["total"] == 0
    assert not (tmp_path / "plddt").exists()


def test_load_plddt_absent_is_none(tmp_path):
    assert load_plddt(tmp_path, "nope") is None


def test_plddt_path_shape(tmp_path):
    assert plddt_path(tmp_path, "rec") == tmp_path / "plddt" / "rec.json"
