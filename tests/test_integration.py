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
