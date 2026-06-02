from afdb_query.sequences import MIN_LENGTH, STANDARD_AA, filter_reason

VALID = "ACDEFGHIKLMNPQRSTVWY"  # 20 residues, all standard


def test_valid_sequence_returns_none():
    assert filter_reason(VALID) is None


def test_internal_stop():
    assert filter_reason("ACD*EFGHIKLMNPQRSTVWY") == "internal_stop"


def test_too_short():
    assert filter_reason("ACDEF") == "too_short"


def test_nonstandard_aa():
    # 'B' and 'Z' are not standard amino acids
    assert filter_reason("ACDEFGHIKLMNPQRSTVWB") == "nonstandard_aa"


def test_internal_stop_takes_priority_over_length():
    # A short sequence containing '*' reports the stop first
    assert filter_reason("AC*") == "internal_stop"


def test_constants():
    assert MIN_LENGTH == 20
    assert "A" in STANDARD_AA and "B" not in STANDARD_AA
