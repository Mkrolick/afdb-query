import json

import httpx
import respx

from afdb_query.client import AlphaFold

SUMMARY_URL = "https://alphafold.ebi.ac.uk/api/sequence/summary"

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
