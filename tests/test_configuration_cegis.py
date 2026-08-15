import sympy as sp
import pytest

from cegis import CertifiedCandidate, candidate_extension
from configuration_cegis import (
    augmented_edges_from_pool,
    certified_candidate_unit_edges,
    configuration_cegis_step,
    deletion_minimal_blocker,
    map_selected_candidates,
    resume_phase1_state,
)
from exact_geometry import Point


def _candidate(x: sp.Expr, y: sp.Expr) -> CertifiedCandidate:
    return CertifiedCandidate(
        point=Point(x, y),
        neighbors=(),
        generator_pair=(0, 1),
    )


def test_candidate_interaction_edges_are_exactly_certified() -> None:
    candidates = [
        _candidate(sp.Integer(0), sp.Integer(0)),
        _candidate(sp.Integer(1), sp.Integer(0)),
        _candidate(sp.Integer(0), sp.Integer(1)),
        _candidate(sp.Rational(2001, 1000), sp.Integer(0)),
    ]

    edges, stats = certified_candidate_unit_edges(candidates)

    assert edges == {(0, 1), (0, 2)}
    assert stats.exact_unit_edges == 2
    assert stats.numeric_pairs >= 2


def test_v_configuration_is_a_deletion_minimal_blocker() -> None:
    base_coloring = (0, 1, 2, 3, 4)
    neighbors = {
        "forced_four": (0, 1, 2, 3),
        "forced_three": (0, 1, 2, 4),
        "center": (0, 1, 2),
    }
    edges = {
        ("forced_four", "center"),
        ("forced_three", "center"),
    }

    core = deletion_minimal_blocker(neighbors, edges, base_coloring, 5)

    assert core is not None
    assert set(core) == set(neighbors)
    assert candidate_extension(
        neighbors, edges, base_coloring, 5, subset=core
    ) is None
    for removed in core:
        subset = [candidate for candidate in core if candidate != removed]
        assert candidate_extension(
            neighbors, edges, base_coloring, 5, subset=subset
        ) is not None


def test_selected_candidates_are_mapped_by_exact_point() -> None:
    old = [_candidate(sp.Integer(0), sp.Integer(0)), _candidate(1, 0)]
    expanded = [_candidate(1, 0), _candidate(0, 0), _candidate(0, 1)]

    assert map_selected_candidates([0, 1], old, expanded) == {0, 1}


def test_selected_candidate_mapping_rejects_missing_exact_point() -> None:
    old = [_candidate(0, 0), _candidate(2, 0)]
    expanded = [_candidate(0, 0), _candidate(0, 1)]

    with pytest.raises(ValueError, match="absent from expanded candidate pool"):
        map_selected_candidates([1], old, expanded)


def test_phase1_state_uses_exact_report_points_and_last_stagnation_coloring() -> None:
    report = {
        "candidate_pool": [
            {
                "x_exact": "0",
                "y_exact": "0",
                "neighbors_zero_based": [0, 1],
                "generator_pair_zero_based": [0, 1],
            },
            {
                "x_exact": "1",
                "y_exact": "0",
                "neighbors_zero_based": [1, 2],
                "generator_pair_zero_based": [1, 2],
            },
        ],
        "cegis": {
            "selected_candidate_indices": [1],
            "stagnation_analysis": [
                {"full_pool_coloring": [0, 1, 2, 3, 4, 0, 1]},
                {"full_pool_coloring": [0, 1, 2, 4, 3, 1, 0]},
            ],
        },
    }
    expanded = [_candidate(1, 0), _candidate(0, 0), _candidate(0, 1)]

    selected, coloring = resume_phase1_state(
        report, expanded, base_vertex_count=5
    )

    assert selected == {0}
    assert coloring == (0, 1, 2, 4, 3)


def test_augmented_edges_use_a_deterministic_selected_vertex_order() -> None:
    candidates = [
        CertifiedCandidate(Point(0, 1), (0,), (0, 1)),
        CertifiedCandidate(Point(1, 1), (1,), (0, 1)),
    ]

    order, edges = augmented_edges_from_pool(
        2,
        {(0, 1)},
        candidates,
        {(0, 1)},
        {1, 0},
    )

    assert order == (0, 1)
    assert edges == {(0, 1), (0, 2), (1, 3), (2, 3)}


def test_configuration_step_blocks_fixed_coloring_then_finds_new_model() -> None:
    candidates = [
        CertifiedCandidate(Point(0, 0), (0, 1, 2, 3), (0, 1)),
        CertifiedCandidate(Point(1, 0), (0, 1, 2, 4), (0, 1)),
        CertifiedCandidate(Point(sp.Rational(1, 2), 0), (0, 1, 2), (0, 1)),
    ]

    step = configuration_cegis_step(
        base_vertex_count=5,
        base_edges=set(),
        candidates=candidates,
        candidate_edges={(0, 2), (1, 2)},
        selected_indices=set(),
        base_coloring=(0, 1, 2, 3, 4),
        color_count=5,
    )

    assert step.result == "SAT_COUNTEREXAMPLE"
    assert set(step.blocking_core) == {0, 1, 2}
    assert step.selected_indices == (0, 1, 2)
    assert step.new_base_coloring is not None
    assert step.new_base_coloring != (0, 1, 2, 3, 4)
