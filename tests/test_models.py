import pytest

from afdb_query.errors import AFDBError, InvalidSequenceError


def test_invalid_sequence_error_is_afdb_error():
    err = InvalidSequenceError("too_short")
    assert isinstance(err, AFDBError)
    assert err.reason == "too_short"
    assert "too_short" in str(err)
