from afdb_query.selection import is_canonical_model, is_monomer, rank_key, select


def _s(model_id, *, oligomeric_state="MONOMER", chains=1, gp=90.0):
    return {
        "model_identifier": model_id,
        "model_url": f"https://alphafold.ebi.ac.uk/files/{model_id}-model_v6.cif",
        "sequence_identity": 1.0,
        "coverage": 1.0,
        "confidence_avg_local_score": gp,
        "oligomeric_state": oligomeric_state,
        "entities": [{"chain_ids": [chr(65 + i) for i in range(chains)]}],
    }


def _wrap(*summaries):
    return [{"summary": s} for s in summaries]


# -- is_monomer ------------------------------------------------------------


def test_monomer_is_monomer():
    assert is_monomer(_s("AF-P00001-F1"))


def test_homodimer_is_not_monomer():
    assert not is_monomer(_s("AF-P00001-F1", oligomeric_state="HOMODIMER", chains=2))


def test_heterodimer_is_not_monomer():
    assert not is_monomer(_s("AF-P00001-F1", oligomeric_state="HETERODIMER"))


def test_missing_state_falls_back_to_chain_count():
    assert is_monomer(_s("AF-P00001-F1", oligomeric_state=None, chains=1))
    assert not is_monomer(_s("AF-P00001-F1", oligomeric_state=None, chains=2))


def test_no_entities_is_monomer():
    assert is_monomer({"oligomeric_state": "MONOMER"})


# -- is_canonical_model ----------------------------------------------------


def test_canonical_f1_id():
    assert is_canonical_model(_s("AF-P12345-F1"))


def test_numeric_ab_initio_id_is_not_canonical():
    assert not is_canonical_model(_s("AF-0000000065764032"))


def test_missing_id_is_not_canonical():
    assert not is_canonical_model({})


# -- select ----------------------------------------------------------------


def test_select_empty_returns_none():
    assert select([]) is None


def test_select_skips_items_without_summary():
    assert select([{"summary": None}]) is None


def test_select_prefers_monomer_over_complex():
    dimer = _s("AF-P00001-F1", oligomeric_state="HOMODIMER", chains=2)
    mono = _s("AF-P00002-F1")
    assert select(_wrap(dimer, mono))["model_identifier"] == "AF-P00002-F1"


def test_select_prefers_canonical_over_numeric():
    numeric = _s("AF-0000000065764032")
    canonical = _s("AF-P99999-F1")
    assert select(_wrap(numeric, canonical))["model_identifier"] == "AF-P99999-F1"


def test_select_accession_wins_over_everything():
    other = _s("AF-P00001-F1")
    mine = _s("AF-P77777-F1")
    got = select(_wrap(other, mine), accession="P77777")
    assert got["model_identifier"] == "AF-P77777-F1"


def test_select_ignores_confidence_score():
    """The whole point: a higher pLDDT must not pull a candidate to the front."""
    low = _s("AF-P00001-F1", gp=10.0)
    high = _s("AF-P00002-F1", gp=99.0)
    # Same monomer/canonical tier, so only the lexicographic id decides -- not the score.
    assert select(_wrap(high, low))["model_identifier"] == "AF-P00001-F1"
    assert select(_wrap(low, high))["model_identifier"] == "AF-P00001-F1"


def test_select_is_order_independent():
    a, b, c = _s("AF-P00003-F1"), _s("AF-P00001-F1"), _s("AF-P00002-F1")
    assert select(_wrap(a, b, c)) == select(_wrap(c, a, b))


def test_rank_key_is_a_total_order_without_scores():
    a = _s("AF-P00001-F1", gp=1.0)
    b = _s("AF-P00001-F1", gp=99.0)
    assert rank_key(a) == rank_key(b)
