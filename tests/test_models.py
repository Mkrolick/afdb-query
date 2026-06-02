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
