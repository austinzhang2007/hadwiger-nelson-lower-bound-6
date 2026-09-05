import json
from pathlib import Path

import pytest

from attach_validated_coloring import attach_validated_coloring
from coloring_sat import write_dimacs_edge


def write_inputs(tmp_path: Path, coloring):
    edge = tmp_path / "graph.edge"
    report = tmp_path / "report.json"
    result = tmp_path / "coloring.json"
    write_dimacs_edge(edge, 3, {(0, 1), (1, 2)})
    report.write_text(
        json.dumps(
            {
                "status": "PENDING",
                "vertices": 3,
                "edges": 2,
                "edge_path": str(edge),
                "lower_bound_6_proved": False,
                "whole_graph_5_sat": None,
                "whole_graph_coloring": None,
            }
        ),
        encoding="utf-8",
    )
    result.write_text(
        json.dumps(
            {
                "status": "VALIDATED_SAT_COLORING",
                "solver": "test-solver",
                "source_log": "test.log",
                "coloring": coloring,
            }
        ),
        encoding="utf-8",
    )
    return report, result


def test_attach_revalidates_and_updates_report_atomically(tmp_path: Path):
    report, result = write_inputs(tmp_path, [0, 1, 0])

    payload = attach_validated_coloring(
        target_report_path=report,
        coloring_result_path=result,
    )

    assert payload["status"] == "EXACT_INTERACTIONS_AUGMENTED_WITH_5_COLORING"
    assert payload["whole_graph_5_sat"] is True
    assert payload["whole_graph_coloring"] == [0, 1, 0]
    assert payload["validated_5_coloring"]["solver"] == "test-solver"
    assert json.loads(report.read_text()) == payload


def test_attach_rejects_a_monochromatic_target_edge(tmp_path: Path):
    report, result = write_inputs(tmp_path, [0, 0, 1])

    with pytest.raises(ValueError, match="violates"):
        attach_validated_coloring(
            target_report_path=report,
            coloring_result_path=result,
        )
