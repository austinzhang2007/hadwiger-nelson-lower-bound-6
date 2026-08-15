from pathlib import Path

from graph_verifier import verify_author_graph, verify_drat


ROOT = Path(__file__).resolve().parents[1]


def test_heule_529_edge_file_matches_all_exact_unit_distances() -> None:
    report = verify_author_graph(
        ROOT / "third_party/CNP-SAT/vtx/529.vtx",
        ROOT / "third_party/CNP-SAT/edge/529.edge",
    )
    assert report.vertex_count == 529
    assert report.edge_count == 2670
    assert report.missing_edges == ()
    assert report.extra_edges == ()


def test_published_heule_drat_verifies_with_reconstructed_sbp_input() -> None:
    report = verify_drat(
        ROOT / "third_party/drat-trim/drat-trim",
        ROOT / "data/heule529/529-4-sbp.cnf",
        ROOT / "third_party/CNP-SAT/proof/529-4-sbp.drat",
    )
    assert report.valid, report.summary
