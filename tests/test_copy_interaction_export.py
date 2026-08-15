import json
from pathlib import Path

from coloring_sat import parse_dimacs_edge
from copy_interactions import export_interaction_stages
from graph_verifier import load_graph_json


ROOT = Path(__file__).resolve().parents[1]


def test_export_interaction_stages_writes_exact_graph_edge_and_cnf(
    tmp_path: Path,
) -> None:
    artifacts = export_interaction_stages(
        graph_path=(
            ROOT
            / "artifacts/cegis/heule2510-primary-d4-cegis.graph.json"
        ),
        report_path=(
            ROOT / "artifacts/cegis/heule2510-primary-d4-cegis.json"
        ),
        source_vertex_path=ROOT / "third_party/CNP-SAT/vtx/529.vtx",
        output_directory=tmp_path,
    )

    assert [artifact.edge_count for artifact in artifacts] == [
        52026,
        52681,
        53332,
    ]
    assert [artifact.added_interaction_edges for artifact in artifacts] == [
        648,
        1303,
        1954,
    ]
    for stage_number, artifact in enumerate(artifacts, 1):
        assert artifact.graph_path.is_file()
        assert artifact.edge_path.is_file()
        assert artifact.cnf_path.is_file()
        assert artifact.report_path.is_file()
        points, graph_edges = load_graph_json(artifact.graph_path)
        edge_vertices, dimacs_edges = parse_dimacs_edge(artifact.edge_path)
        assert len(points) == edge_vertices == 8329
        assert graph_edges == dimacs_edges
        with artifact.cnf_path.open(encoding="ascii") as handle:
            assert handle.readline().startswith("p cnf 41645 ")
        with artifact.report_path.open(encoding="utf-8") as handle:
            report = json.load(handle)
        assert report["status"] == "UNSOLVED"
        assert report["complete_unit_distance_graph"] is (stage_number == 3)
