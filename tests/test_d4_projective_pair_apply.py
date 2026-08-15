import pytest

from d4_projective_pair_apply import validate_pair_scan


def test_validate_pair_scan_maps_sparse_candidate_indices() -> None:
    payload = {
        "missing_color": 0,
        "pair_edges": [[7, 11], [11, 19]],
        "used_candidates": [
            {"candidate_index": 19},
            {"candidate_index": 7},
            {"candidate_index": 11},
        ],
    }

    records, edges = validate_pair_scan(payload)

    assert sorted(records) == [7, 11, 19]
    assert edges == ((7, 11), (11, 19))


def test_validate_pair_scan_rejects_missing_pair_endpoint() -> None:
    payload = {
        "missing_color": 0,
        "pair_edges": [[7, 11]],
        "used_candidates": [{"candidate_index": 7}],
    }

    with pytest.raises(ValueError, match="missing candidate"):
        validate_pair_scan(payload)
