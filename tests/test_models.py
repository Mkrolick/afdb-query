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


def test_from_dict():
    p = Plddt.from_dict(
        {"confidenceScore": [5.0, 6.0], "residueNumber": [1, 2], "confidenceCategory": ["D", "D"]}
    )
    assert p.scores == [5.0, 6.0]
    assert p.residue_numbers == [1, 2]
    assert p.raw["confidenceCategory"] == ["D", "D"]


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


def test_structure_oligomeric_state_and_monomer():
    assert Structure(SUMMARY, None).oligomeric_state is None  # SUMMARY leaves it unset
    assert Structure(SUMMARY, None).is_monomer
    dimer = Structure(
        {"oligomeric_state": "HOMODIMER", "entities": [{"chain_ids": ["A", "B"]}]}, None
    )
    assert dimer.oligomeric_state == "HOMODIMER"
    assert not dimer.is_monomer


def test_structure_no_entities_returns_none():
    s = Structure({}, None)
    assert s.uniprot_accession is None
    assert s.description is None


def test_structure_description_missing_returns_none():
    s = Structure({"entities": [{"identifier": "x", "identifier_category": "PDB"}]}, None)
    assert s.description is None


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
    assert p1.residue_numbers == [1, 2, 3]
    assert p1 is p2  # cached on the instance
    assert route.call_count == 1  # fetched once
