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
