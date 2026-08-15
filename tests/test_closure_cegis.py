import json
from pathlib import Path

import pytest

from cegis import CertifiedCandidate
from closure_cegis import (
    finalize_closure_checkpoint,
    load_exact_colored_graph,
    run_direct_closure_audit,
    solve_full_candidate_pool,
    write_checkpoint,
)
from exact_geometry import Point


def _write_graph_and_report(tmp_path, coloring):
    graph_path = tmp_path / "graph.json"
    report_path = tmp_path / "report.json"
    graph_path.write_text(
        json.dumps(
            {
                "schema": "udg-exact-v1",
                "coordinate_format": "sympy-radical-v1",
                "name": "unit edge",
                "vertices": [
                    {"id": 1, "x": "0", "y": "0"},
                    {"id": 2, "x": "1", "y": "0"},
                ],
                "edges": [[1, 2]],
            }
        ),
        encoding="utf-8",
    )
    report_path.write_text(
        json.dumps(
            {"whole_pool_failure_certificate": {"coloring": coloring}}
        ),
        encoding="utf-8",
    )
    return graph_path, report_path


def test_exact_colored_graph_loader_accepts_valid_coloring(tmp_path) -> None:
    graph_path, report_path = _write_graph_and_report(tmp_path, [0, 1])

    points, edges, coloring = load_exact_colored_graph(
        graph_path, report_path, color_count=5
    )

    assert len(points) == 2
    assert edges == {(0, 1)}
    assert coloring == (0, 1)


def test_exact_colored_graph_loader_rejects_edge_conflict(tmp_path) -> None:
    graph_path, report_path = _write_graph_and_report(tmp_path, [2, 2])

    with pytest.raises(ValueError, match="invalid on the exact base graph"):
        load_exact_colored_graph(graph_path, report_path, color_count=5)


def test_checkpoint_round_trips_resume_state(tmp_path) -> None:
    path = tmp_path / "checkpoint.json"
    report = {
        "status": "RUNNING",
        "resume_state": {
            "next_iteration": 7,
            "selected_candidate_indices": [2, 5],
            "current_base_coloring": [0, 1, 2],
        },
    }

    write_checkpoint(path, report)

    assert json.loads(path.read_text(encoding="utf-8")) == report
    assert not path.with_suffix(".json.tmp").exists()


def test_full_candidate_pool_is_solved_without_fixing_base_colors() -> None:
    candidates = [
        CertifiedCandidate(Point(0, 1), (0,), (0, 1)),
        CertifiedCandidate(Point(1, 1), (1,), (0, 1)),
    ]

    edges, coloring = solve_full_candidate_pool(
        base_vertex_count=2,
        base_edges={(0, 1)},
        candidates=candidates,
        candidate_edges={(0, 1)},
        color_count=3,
    )

    assert edges == {(0, 1), (0, 2), (1, 3), (2, 3)}
    assert coloring is not None
    assert all(coloring[u] != coloring[v] for u, v in edges)


def test_checkpoint_finalizer_saves_full_pool_coloring(tmp_path) -> None:
    graph_path, _ = _write_graph_and_report(tmp_path, [0, 1])
    checkpoint_path = tmp_path / "closure.json"
    checkpoint_path.write_text(
        json.dumps(
            {
                "status": "ITERATION_LIMIT",
                "candidate_pool": [
                    {
                        "x_exact": "0",
                        "y_exact": "1",
                        "neighbors_zero_based": [0],
                        "generator_pair_zero_based": [0, 1],
                    },
                    {
                        "x_exact": "1",
                        "y_exact": "1",
                        "neighbors_zero_based": [1],
                        "generator_pair_zero_based": [0, 1],
                    },
                ],
                "candidate_unit_edges_zero_based": [[0, 1]],
            }
        ),
        encoding="utf-8",
    )

    report = finalize_closure_checkpoint(
        graph_path, checkpoint_path, color_count=3
    )

    assert report["status"] == "EXACT_LISTED_EDGE_WHOLE_POOL_5_COLORING"
    certificate = report["whole_pool_failure_certificate"]
    assert certificate["vertices"] == 4
    assert certificate["certified_listed_edges"] == 4
    assert Path(certificate["listed_graph_path"]).exists()


def test_direct_closure_audit_builds_and_solves_complete_exact_graph(
    tmp_path,
) -> None:
    graph_path, _ = _write_graph_and_report(tmp_path, [0, 1])
    output_path = tmp_path / "direct.json"

    report = run_direct_closure_audit(
        graph_path,
        output_path,
        minimum_candidate_degree=2,
        color_count=3,
    )

    assert report["candidate_generation"]["certified_candidates"] == 2
    assert report["complete_graph"]["vertices"] == 4
    assert report["complete_graph"]["edges"] == 5
    assert report["status"] == "EXACT_COMPLETE_WHOLE_POOL_COLORING"
    assert Path(report["complete_graph"]["path"]).exists()
