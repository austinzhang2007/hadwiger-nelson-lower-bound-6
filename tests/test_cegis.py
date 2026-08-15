import sympy as sp

from cegis import (
    CertifiedCandidate,
    blocking_candidates,
    candidate_pool_extension_conclusion,
    candidate_extension,
    generate_certified_unit_circle_candidates,
    greedy_cover,
    interacting_pair_coverage,
    recertify_candidate_neighbors,
)
from exact_geometry import Point


def test_candidate_blocks_coloring_only_when_all_colors_are_neighbors() -> None:
    candidate_neighbors = {
        Point(sp.Integer(0), sp.Integer(0)): (0, 1, 2, 3, 4),
        Point(sp.Integer(2), sp.Integer(0)): (0, 1, 2, 3),
    }
    colorings = [
        (0, 1, 2, 3, 4),
        (0, 1, 2, 3, 3),
    ]
    coverage = blocking_candidates(candidate_neighbors, colorings, 5)
    assert coverage[Point(sp.Integer(0), sp.Integer(0))] == frozenset({0})
    assert coverage[Point(sp.Integer(2), sp.Integer(0))] == frozenset()


def test_greedy_cover_reports_uncovered_colorings() -> None:
    coverage = {"a": frozenset({0, 2}), "b": frozenset({1})}
    chosen, uncovered = greedy_cover(coverage, universe={0, 1, 2, 3})
    assert chosen == ["a", "b"]
    assert uncovered == {3}


def test_unit_circle_candidate_generation_recovers_exact_intersections() -> None:
    points = [
        Point(sp.Integer(0), sp.Integer(0)),
        Point(sp.Integer(1), sp.Integer(0)),
    ]
    candidates, stats = generate_certified_unit_circle_candidates(
        points, min_neighbors=2
    )
    assert stats.center_pairs == 1
    assert len(candidates) == 2
    assert {candidate.neighbors for candidate in candidates} == {(0, 1)}
    assert {candidate.point.y for candidate in candidates} == {
        -sp.sqrt(3) / 2,
        sp.sqrt(3) / 2,
    }


def test_adjacent_candidate_pair_can_jointly_block_a_coloring() -> None:
    candidate_neighbors = {
        "left": (0, 1, 2, 3),
        "right": (0, 1, 2, 3),
    }
    colorings = [(0, 1, 2, 3)]
    singles = blocking_candidates(candidate_neighbors, colorings, 5)
    assert singles["left"] == singles["right"] == frozenset()
    pairs = interacting_pair_coverage(
        candidate_neighbors,
        [("left", "right")],
        colorings,
        5,
    )
    assert pairs[("left", "right")] == frozenset({0})


def test_three_candidate_configuration_can_block_without_pair_blocking() -> None:
    # Base colors 2,3,4 forbid those colors at every candidate, leaving {0,1}.
    candidate_neighbors = {0: (0, 1, 2), 1: (0, 1, 2), 2: (0, 1, 2)}
    base_coloring = (2, 3, 4)
    assert candidate_extension(
        candidate_neighbors,
        [(0, 1), (1, 2)],
        base_coloring,
        5,
    ) is not None
    assert candidate_extension(
        candidate_neighbors,
        [(0, 1), (1, 2), (0, 2)],
        base_coloring,
        5,
    ) is None


def test_candidate_pool_extension_conclusion_uses_dynamic_pool_size() -> None:
    assert candidate_pool_extension_conclusion(248) == (
        "No subset of this 248-candidate pool can block this base coloring."
    )


def test_candidate_neighbors_are_rebuilt_from_all_exact_base_points() -> None:
    base = [Point(0, 0), Point(2, 0), Point(0, 1)]
    candidate = CertifiedCandidate(
        point=Point(1, 0),
        neighbors=(0,),
        generator_pair=(0, 1),
    )

    rebuilt = recertify_candidate_neighbors(base, [candidate])

    assert rebuilt[0].neighbors == (0, 1)
