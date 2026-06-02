import pytest

from afdb_query.errors import AFDBError, InvalidSequenceError
from afdb_query.models import confidence_url


def test_invalid_sequence_error_is_afdb_error():
    err = InvalidSequenceError("too_short")
    assert isinstance(err, AFDBError)
    assert err.reason == "too_short"
    assert "too_short" in str(err)


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
