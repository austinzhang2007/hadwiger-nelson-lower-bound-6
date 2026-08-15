import json
from pathlib import Path

import sympy as sp

from cegis import candidate_extension
from coloring_sat import parse_dimacs_edge
from exact_geometry import Point, parse_mathematica_vertices
from graph_verifier import load_graph_json
from translated_copy_cegis import (
    certified_unit_shift_copy_edges,
    circle_translation,
    enumerate_structured_cross_edges,
    primary_d4_candidate_pool,
    translated_points,
)
from translated_copy_verifier import verify_translated_copy_artifact


ROOT = Path(__file__).resolve().parents[1]
BASE_GRAPH = (
    ROOT
    / "artifacts/cegis/heule1061-second-closure-recertified-complete.graph.json"
)
BASE_REPORT = (
    ROOT / "artifacts/cegis/heule1061-second-closure-recertified.json"
)


def test_circle_translation_recovers_exact_nonabelian_extension() -> None:
    points, _ = load_graph_json(BASE_GRAPH)
    r = points[97]

    translation = circle_translation(r, branch=1)

    assert sp.expand(translation.h_squared - (25 + 3 * sp.sqrt(33)) / 8) == 0
    x = sp.symbols("x")
    assert translation.minimal_polynomial == sp.Poly(
        8 * x**4 - 50 * x**2 + 41, x, domain=sp.QQ
    )
    assert translation.galois_group_is_nonabelian
    assert translation.point.squared_distance(Point(0, 0)) == 1
    assert translation.point.squared_distance(r) == 1


def test_structured_cross_edges_are_exact_on_two_anchor_example() -> None:
    points, _ = load_graph_json(BASE_GRAPH)
    r = points[97]
    translation = circle_translation(r, branch=1)
    source = [Point(0, 0)]
    base = [Point(0, 0), r, Point(3, 0)]
    shifted = translated_points(source, translation.point)

    cross_edges = enumerate_structured_cross_edges(
        base,
        source,
        r,
        shifted,
        translation,
    )

    assert cross_edges == {(0, 0), (1, 0)}
    assert all(
        base[base_index].squared_distance(shifted[source_index]) == 1
        for base_index, source_index in cross_edges
    )


def test_heule_copy_has_772_cross_edges_and_blocks_saved_coloring() -> None:
    base, _ = load_graph_json(BASE_GRAPH)
    with BASE_REPORT.open(encoding="utf-8") as handle:
        report = json.load(handle)
    base_coloring = tuple(
        report["whole_pool_failure_certificate"]["coloring"]
    )
    source = parse_mathematica_vertices(
        ROOT / "third_party/CNP-SAT/vtx/529.vtx"
    )
    _, source_edges = parse_dimacs_edge(
        ROOT / "third_party/CNP-SAT/edge/529.edge"
    )
    r = base[97]
    translation = circle_translation(r, branch=1)
    shifted = translated_points(source, translation.point)

    cross_edges = enumerate_structured_cross_edges(
        base,
        source,
        r,
        shifted,
        translation,
    )
    neighbors = {index: [] for index in range(len(source))}
    for base_index, source_index in cross_edges:
        neighbors[source_index].append(base_index)

    assert len(cross_edges) == 772
    assert sum(len(indices) == 2 for indices in neighbors.values()) == 243
    assert candidate_extension(
        neighbors,
        source_edges,
        base_coloring,
        5,
    ) is None


def test_primary_d4_candidate_pool_is_exact_and_reproducible() -> None:
    base, _ = load_graph_json(BASE_GRAPH)
    source = parse_mathematica_vertices(
        ROOT / "third_party/CNP-SAT/vtx/529.vtx"
    )

    pool = primary_d4_candidate_pool(base, source)

    assert len(pool) == 24
    assert [(candidate.second_neighbor_count, candidate.source_r_index) for candidate in pool[:3]] == [
        (243, 201),
        (243, 97),
        (241, 206),
    ]
    assert all(
        sum(len(indices) == 2 for indices in candidate.neighbors.values())
        == candidate.second_neighbor_count
        for candidate in pool
    )


def test_independent_verifier_accepts_primary_d4_artifact() -> None:
    result = verify_translated_copy_artifact(
        ROOT / "artifacts/cegis/heule2510-primary-d4-cegis.graph.json",
        ROOT / "artifacts/cegis/heule2510-primary-d4-cegis.json",
        BASE_GRAPH,
        ROOT / "third_party/CNP-SAT/vtx/529.vtx",
        ROOT / "third_party/CNP-SAT/edge/529.edge",
    )

    assert result["valid"]
    assert result["vertices"] == 8329
    assert result["listed_edges"] == 51378
    assert result["duplicate_vertices"] == 0
    assert result["coloring_conflicts"] == 0


def test_three_selected_copy_pairs_have_1954_certified_interaction_edges() -> None:
    source = parse_mathematica_vertices(
        ROOT / "third_party/CNP-SAT/vtx/529.vtx"
    )
    selected = [201, 101, 216, 97, 104, 135, 184, 118, 178, 206, 124]

    interactions, pair_counts = certified_unit_shift_copy_edges(
        source,
        selected,
    )

    assert pair_counts == {(0, 6): 648, (1, 10): 655, (3, 7): 651}
    assert len(interactions) == 1954
