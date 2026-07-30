import pytest

from afdb_query.plddt import (
    is_contiguous,
    mean_per_residue,
    mean_plddt,
    residue_index,
    shared_suffix_means,
)

# -- mean_plddt ------------------------------------------------------------


def test_mean_of_whole_array():
    assert mean_plddt([10.0, 20.0, 30.0]) == 20.0


def test_mean_with_start():
    assert mean_plddt([0.0, 30.0, 50.0], start=1) == 40.0


def test_mean_with_stop():
    assert mean_plddt([30.0, 50.0, 0.0], stop=2) == 40.0


def test_mean_with_both_bounds():
    assert mean_plddt([0.0, 10.0, 20.0, 99.0], start=1, stop=3) == 15.0


def test_mean_of_empty_is_none():
    assert mean_plddt([]) is None


def test_mean_of_empty_slice_is_none():
    """None, not 0.0 -- 'no overlap' must be distinguishable from 'overlap scoring 0'."""
    assert mean_plddt([1.0, 2.0], start=5) is None
    assert mean_plddt([1.0, 2.0], stop=0) is None


def test_mean_does_not_round():
    assert mean_plddt([1.0, 2.0]) == 1.5
    assert mean_plddt([1.0, 1.0, 2.0]) == pytest.approx(4 / 3)


# -- residue_index ---------------------------------------------------------


def test_residue_index_found():
    assert residue_index([1, 2, 3], 2) == 1


def test_residue_index_non_contiguous_numbering():
    """The point of the helper: position != residue number when there are gaps."""
    assert residue_index([5, 6, 9], 9) == 2


def test_residue_index_missing_is_none():
    assert residue_index([1, 2, 3], 99) is None


# -- is_contiguous ---------------------------------------------------------


def test_contiguous_true():
    assert is_contiguous([1, 2, 3, 4])


def test_contiguous_false_on_gap():
    assert not is_contiguous([1, 2, 4])


def test_contiguous_false_when_not_one_based():
    assert not is_contiguous([0, 1, 2])


def test_contiguous_empty_is_true():
    assert is_contiguous([])


# -- shared_suffix_means ---------------------------------------------------


def test_shared_suffix_means_length_cancels():
    """A disordered N-terminus drags the global mean but not the shared-region mean."""
    long_scores = [10.0, 10.0, 90.0, 90.0]  # global mean 50
    short_scores = [90.0, 90.0]  # global mean 90 -- the naive "wins by 40"
    got = shared_suffix_means(long_scores, short_scores)
    assert got["offset"] == 2
    assert got["shared_long"] == 90.0
    assert got["shared_short"] == 90.0
    assert got["shared_long"] == got["shared_short"]  # no difference once length cancels
    assert got["displaced"] == 10.0  # the removed segment is what moved it


def test_shared_suffix_means_detects_real_difference():
    """When the shared region genuinely folds better, the comparison still says so."""
    got = shared_suffix_means([50.0, 20.0, 20.0], [80.0, 80.0])
    assert got["shared_long"] == 20.0
    assert got["shared_short"] == 80.0
    assert got["displaced"] == 50.0


def test_shared_suffix_means_equal_lengths():
    got = shared_suffix_means([1.0, 3.0], [10.0, 20.0])
    assert got["offset"] == 0
    assert got["shared_long"] == 2.0
    assert got["shared_short"] == 15.0
    assert got["displaced"] is None  # nothing was displaced


def test_shared_suffix_means_rejects_inverted_arguments():
    with pytest.raises(ValueError):
        shared_suffix_means([1.0], [1.0, 2.0])


# -- mean_per_residue ------------------------------------------------------


def test_mean_per_residue_averages_elementwise():
    assert mean_per_residue([[10.0, 20.0], [20.0, 40.0]]) == [15.0, 30.0]


def test_mean_per_residue_single_array_is_itself():
    assert mean_per_residue([[1.0, 2.0, 3.0]]) == [1.0, 2.0, 3.0]


def test_mean_per_residue_refuses_ragged():
    """zip would truncate and return a plausible number for a comparison that does not exist."""
    with pytest.raises(ValueError, match="ragged"):
        mean_per_residue([[1.0, 2.0, 3.0], [1.0, 2.0]])


def test_mean_per_residue_enforces_expected_length():
    with pytest.raises(ValueError, match="expected"):
        mean_per_residue([[1.0, 2.0]], expected_length=430)


def test_mean_per_residue_accepts_matching_expected_length():
    assert mean_per_residue([[1.0, 2.0]], expected_length=2) == [1.0, 2.0]


def test_mean_per_residue_refuses_empty_input():
    with pytest.raises(ValueError):
        mean_per_residue([])


def test_mean_per_residue_refuses_empty_arrays():
    with pytest.raises(ValueError):
        mean_per_residue([[], []])


def test_mean_per_residue_commutes_with_region_mean():
    """The property the whole design rests on: slice-then-mean == mean-then-slice."""
    a = [10.0, 10.0, 90.0, 80.0]
    b = [20.0, 30.0, 70.0, 60.0]
    consensus = mean_per_residue([a, b])
    via_consensus = mean_plddt(consensus, start=2)
    via_members = (mean_plddt(a, start=2) + mean_plddt(b, start=2)) / 2
    assert via_consensus == pytest.approx(via_members)
