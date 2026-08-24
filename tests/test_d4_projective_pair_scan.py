from d4_exact import D4ProjectivePoint
from d4_pair_search import ForcedPoint, ForcedPointStats
from d4_projective_pair_scan import build_pair_scan_payload
from exact_geometry import Point


def _projective(marker: int) -> D4ProjectivePoint:
    coefficients = (marker,) + (0,) * 31
    denominator = (1,) + (0,) * 31
    return D4ProjectivePoint(coefficients, coefficients, denominator, -1)


def test_build_pair_scan_payload_keeps_only_incident_candidates() -> None:
    candidates = [
        ForcedPoint(
            point=Point(0, 0),
            projective=_projective(index + 1),
            neighbors=(index, index + 1, index + 2, index + 3),
            forced_color=4,
            generator_pairs=((index, index + 1), (index + 2, index + 3)),
        )
        for index in range(3)
    ]
    stats = ForcedPointStats(4, 10, 4, 6, 3, 3, 3)

    payload = build_pair_scan_payload(
        source_report="source.json",
        missing_color=4,
        candidates=candidates,
        pair_edges={(0, 2)},
        stats=stats,
        timings_seconds={"forced_recovery": 1.25, "pair_filter": 0.5},
    )

    assert payload["forced_candidates"] == 3
    assert payload["pair_edges"] == [[0, 2]]
    assert [record["candidate_index"] for record in payload["used_candidates"]] == [0, 2]
    assert payload["search_stats"]["exact_candidates"] == 3
    assert payload["timings_seconds"]["pair_filter"] == 0.5
    assert payload["used_candidates"][0]["neighbors_zero_based"] == [0, 1, 2, 3]
    assert payload["used_candidates"][0]["projective_certificate"]["fifth_neighbor_zero_based"] == -1
