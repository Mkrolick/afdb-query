from afdb_query.selection import (
    exact_matches,
    filter_by_length,
    is_canonical_model,
    is_exact,
    is_monomer,
    mean_global_plddt,
    rank_tiers,
    select_group,
)


def _s(model_id, *, oligomeric_state="MONOMER", chains=1, gp=90.0):
    s = {
        "model_identifier": model_id,
        "model_url": f"https://alphafold.ebi.ac.uk/files/{model_id}-model_v6.cif",
        "sequence_identity": 1.0,
        "coverage": 1.0,
        "oligomeric_state": oligomeric_state,
        "entities": [{"chain_ids": [chr(65 + i) for i in range(chains)]}],
    }
    if gp is not None:
        s["confidence_avg_local_score"] = gp
    return s


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


# -- rank_tiers ------------------------------------------------------------


def test_rank_tiers_ignores_confidence_score():
    """The tiers must be a function of identity and provenance only."""
    assert rank_tiers(_s("AF-P00001-F1", gp=1.0)) == rank_tiers(_s("AF-P00001-F1", gp=99.0))


def test_rank_tiers_has_exactly_three_tiers():
    """A fourth component would be an arbitrary tie-break -- what group-return avoids."""
    assert len(rank_tiers(_s("AF-P00001-F1"))) == 3


# -- select_group ----------------------------------------------------------


def test_select_group_empty_returns_empty_list():
    assert select_group([]) == []


def test_select_group_skips_items_without_summary():
    assert select_group([{"summary": None}]) == []


def test_select_group_drops_complexes_when_a_monomer_exists():
    dimer = _s("AF-P00001-F1", oligomeric_state="HOMODIMER", chains=2)
    mono = _s("AF-P00002-F1")
    got = select_group(_wrap(dimer, mono))
    assert [s["model_identifier"] for s in got] == ["AF-P00002-F1"]


def test_select_group_drops_numeric_when_a_canonical_exists():
    numeric = _s("AF-0000000065764032")
    canonical = _s("AF-P99999-F1")
    got = select_group(_wrap(numeric, canonical))
    assert [s["model_identifier"] for s in got] == ["AF-P99999-F1"]


def test_select_group_keeps_every_tied_candidate():
    """The whole point: equal-tier candidates are all returned, not whittled to one."""
    a, b, c = _s("AF-P00001-F1"), _s("AF-P00002-F1"), _s("AF-P00003-F1")
    got = select_group(_wrap(a, b, c))
    assert [s["model_identifier"] for s in got] == [
        "AF-P00001-F1",
        "AF-P00002-F1",
        "AF-P00003-F1",
    ]


def test_select_group_membership_ignores_confidence_score():
    """A higher pLDDT must neither promote nor demote a candidate."""
    low = _s("AF-P00001-F1", gp=10.0)
    high = _s("AF-P00002-F1", gp=99.0)
    assert len(select_group(_wrap(low, high))) == 2
    assert select_group(_wrap(low, high)) == select_group(_wrap(high, low))


def test_select_group_falls_back_to_complexes_when_nothing_better():
    dimer = _s("AF-P00001-F1", oligomeric_state="HOMODIMER", chains=2)
    got = select_group(_wrap(dimer))
    assert len(got) == 1  # never filters -- the caller decides what to do with it
    assert not is_monomer(got[0])


def test_select_group_accession_pins_to_one():
    other = _s("AF-P00001-F1")
    mine = _s("AF-P77777-F1")
    got = select_group(_wrap(other, mine), accession="P77777")
    assert [s["model_identifier"] for s in got] == ["AF-P77777-F1"]


def test_select_group_order_is_stable_and_input_independent():
    a, b, c = _s("AF-P00003-F1"), _s("AF-P00001-F1"), _s("AF-P00002-F1")
    assert select_group(_wrap(a, b, c)) == select_group(_wrap(c, a, b))


# -- mean_global_plddt -----------------------------------------------------


def test_mean_global_plddt_averages_the_group():
    got = select_group(_wrap(_s("AF-P00001-F1", gp=80.0), _s("AF-P00002-F1", gp=90.0)))
    assert mean_global_plddt(got) == 85.0


def test_mean_global_plddt_is_not_the_max():
    """Regression guard: the old max-over-structures behaviour returned 99.0 here."""
    group = select_group(_wrap(_s("AF-P00001-F1", gp=10.0), _s("AF-P00002-F1", gp=99.0)))
    assert mean_global_plddt(group) == 54.5


def test_mean_global_plddt_single_member():
    assert mean_global_plddt(select_group(_wrap(_s("AF-P00001-F1", gp=77.0)))) == 77.0


def test_mean_global_plddt_empty_is_none():
    assert mean_global_plddt([]) is None


def test_mean_global_plddt_skips_missing_scores():
    """A missing score is an absent measurement, not a zero."""
    group = [_s("AF-P00001-F1", gp=80.0), _s("AF-P00002-F1", gp=None)]
    assert mean_global_plddt(group) == 80.0


def test_mean_global_plddt_all_missing_is_none():
    assert mean_global_plddt([_s("AF-P00001-F1", gp=None)]) is None


def test_mean_global_plddt_ignores_non_numeric():
    assert mean_global_plddt([{"confidence_avg_local_score": "high"}]) is None


# -- filter_by_length ------------------------------------------------------


def test_filter_by_length_keeps_exact_matches():
    a, b = _s("AF-P00001-F1"), _s("AF-P00002-F1")
    kept, dropped = filter_by_length([a, b], {"AF-P00001-F1": 430, "AF-P00002-F1": 430}, 430)
    assert [s["model_identifier"] for s in kept] == ["AF-P00001-F1", "AF-P00002-F1"]
    assert dropped == []


def test_filter_by_length_drops_longer_models():
    """An ortholog with extra terminal residues passes every visible tier but is not the query."""
    good, long_ = _s("AF-P00001-F1"), _s("AF-P00002-F1")
    kept, dropped = filter_by_length([good, long_], {"AF-P00001-F1": 430, "AF-P00002-F1": 437}, 430)
    assert [s["model_identifier"] for s in kept] == ["AF-P00001-F1"]
    assert [s["model_identifier"] for s in dropped] == ["AF-P00002-F1"]


def test_filter_by_length_drops_shorter_models():
    kept, dropped = filter_by_length([_s("AF-P00001-F1")], {"AF-P00001-F1": 12}, 430)
    assert kept == []
    assert len(dropped) == 1


def test_filter_by_length_drops_unknown_length():
    """Unknown is dropped, never assumed to conform."""
    kept, dropped = filter_by_length([_s("AF-P00001-F1")], {}, 430)
    assert kept == []
    assert len(dropped) == 1


def test_filter_by_length_empty_group():
    assert filter_by_length([], {}, 430) == ([], [])


def test_filter_by_length_then_mean_excludes_the_outlier():
    """End to end: the wrong-length member must not reach the average."""
    good = _s("AF-P00001-F1", gp=70.0)
    long_ = _s("AF-P00002-F1", gp=95.0)  # higher score, wrong protein extent
    kept, dropped = filter_by_length([good, long_], {"AF-P00001-F1": 430, "AF-P00002-F1": 860}, 430)
    assert mean_global_plddt(kept) == 70.0
    assert mean_global_plddt(dropped) == 95.0  # visible, not silently absorbed


# -- is_exact / exact_matches ----------------------------------------------


def _inexact(model_id, *, identity=1.0, coverage=1.0):
    s = _s(model_id)
    s["sequence_identity"] = identity
    s["coverage"] = coverage
    return s


def test_exact_requires_identity_and_coverage():
    assert is_exact(_inexact("AF-P1-F1"))


def test_near_match_is_not_exact():
    """The failure this guards: a different protein's confidence attached to a row."""
    assert not is_exact(_inexact("AF-P1-F1", identity=0.42, coverage=0.31))


def test_high_but_imperfect_identity_is_not_exact():
    assert not is_exact(_inexact("AF-P1-F1", identity=0.9995))


def test_full_identity_partial_coverage_is_not_exact():
    """Identity says the aligned part matched; coverage says the alignment spanned the query."""
    assert not is_exact(_inexact("AF-P1-F1", coverage=0.5))


def test_missing_fields_are_not_exact():
    """Unverifiable and verified are different states."""
    assert not is_exact({})
    assert not is_exact({"sequence_identity": 1.0})
    assert not is_exact({"coverage": 1.0})


def test_non_numeric_identity_is_not_exact():
    assert not is_exact(_inexact("AF-P1-F1", identity="1.0"))


def test_integer_one_counts_as_exact():
    assert is_exact({"sequence_identity": 1, "coverage": 1})


def test_exact_matches_splits_and_reports():
    good, bad = _inexact("AF-P1-F1"), _inexact("AF-P2-F1", identity=0.42)
    exact, rejected = exact_matches(_wrap(good, bad))
    assert [s["model_identifier"] for s in exact] == ["AF-P1-F1"]
    assert [s["model_identifier"] for s in rejected] == ["AF-P2-F1"]


def test_exact_matches_empty():
    assert exact_matches([]) == ([], [])


def test_exact_matches_skips_items_without_summary():
    assert exact_matches([{"summary": None}]) == ([], [])


def test_select_group_discards_near_matches():
    """Regression: before this, a 0.42-identity hit sailed through and supplied the pLDDT."""
    junk = _inexact("AF-WRONG-F1", identity=0.42, coverage=0.31)
    junk["confidence_avg_local_score"] = 95.0
    assert select_group(_wrap(junk)) == []
    assert mean_global_plddt(select_group(_wrap(junk))) is None


def test_select_group_keeps_exact_drops_inexact():
    good = _inexact("AF-P1-F1")
    good["confidence_avg_local_score"] = 70.0
    junk = _inexact("AF-P2-F1", identity=0.5)
    junk["confidence_avg_local_score"] = 99.0
    group = select_group(_wrap(good, junk))
    assert [s["model_identifier"] for s in group] == ["AF-P1-F1"]
    assert mean_global_plddt(group) == 70.0  # the near match's 99.0 never enters
