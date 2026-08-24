import json

from tabu_coloring import load_initial_coloring, summary_without_colorings, tabu_repair


def test_load_initial_coloring_accepts_edge_ladder_checkpoint(tmp_path) -> None:
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(json.dumps({"coloring": [0, 1, 2]}), encoding="utf-8")

    assert load_initial_coloring(checkpoint) == [0, 1, 2]


def test_tabu_iteration_limit_preserves_best_coloring(tmp_path) -> None:
    coloring, stats = tabu_repair(
        2,
        {(0, 1)},
        [0, 0],
        colors=2,
        iterations=0,
    )
    assert coloring is None
    assert stats["best_conflicts"] == 1
    assert stats["best_coloring"] == [0, 0]

    report = tmp_path / "tabu.json"
    report.write_text(json.dumps(stats), encoding="utf-8")
    assert load_initial_coloring(report) == [0, 0]


def test_tabu_summary_omits_full_color_arrays() -> None:
    assert summary_without_colorings(
        {
            "status": "UNKNOWN_ITERATION_LIMIT",
            "coloring": None,
            "best_coloring": [0, 1],
            "best_conflicts": 1,
        }
    ) == {"status": "UNKNOWN_ITERATION_LIMIT", "best_conflicts": 1}


def test_tabu_stagnation_kick_is_recorded() -> None:
    coloring, stats = tabu_repair(
        3,
        {(0, 1), (1, 2), (0, 2)},
        [0, 0, 0],
        colors=2,
        seed=67,
        iterations=20,
        sample_size=3,
        stagnation_limit=2,
        kick_size=1,
    )

    assert coloring is None
    assert stats["perturbations"] >= 1
    assert len(stats["best_coloring"]) == 3
