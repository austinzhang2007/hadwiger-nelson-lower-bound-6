from pathlib import Path

from copy_interaction_verifier import verify_complete_interaction_artifact


ROOT = Path(__file__).resolve().parents[1]


def test_independent_verifier_accepts_complete_stage3_graph() -> None:
    result = verify_complete_interaction_artifact(
        graph_path=(
            ROOT
            / "artifacts/interactions/"
            "heule8329-interactions-stage3.graph.json"
        ),
        stage_report_path=(
            ROOT
            / "artifacts/interactions/heule8329-interactions-stage3.json"
        ),
        primary_report_path=(
            ROOT / "artifacts/cegis/heule2510-primary-d4-cegis.json"
        ),
        base_graph_path=(
            ROOT
            / "artifacts/cegis/"
            "heule1061-second-closure-recertified-complete.graph.json"
        ),
        source_vertex_path=ROOT / "third_party/CNP-SAT/vtx/529.vtx",
        source_edge_path=ROOT / "third_party/CNP-SAT/edge/529.edge",
    )

    assert result["valid"]
    assert result["vertices"] == 8329
    assert result["edges"] == 53332
    assert result["interaction_edges"] == 1954
    assert result["pair_certificates"] == 55
    assert result["duplicate_vertices"] == 0
    assert result["missing_edges"] == 0
    assert result["unexpected_edges"] == 0
