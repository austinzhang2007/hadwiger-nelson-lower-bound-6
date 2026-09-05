from core_guided_coloring_repair import repair_coloring_core_guided


def test_core_guided_repair_unfreezes_an_unsat_boundary_core():
    edges = {(0, 1), (0, 2), (1, 3)}
    initial = [0, 0, 1, 1]

    coloring, stats = repair_coloring_core_guided(
        4,
        edges,
        initial,
        colors=2,
        solver_name="glucose4",
    )

    assert coloring is not None
    assert all(coloring[left] != coloring[right] for left, right in edges)
    assert stats["status"] == "SAT"
    assert stats["initial_conflicts"] == 1
    assert stats["final_conflicts"] == 0
    assert stats["rounds"] >= 1


def test_core_guided_full_unsat_is_not_reported_as_a_proved_certificate():
    coloring, stats = repair_coloring_core_guided(
        3,
        {(0, 1), (1, 2), (0, 2)},
        [0, 0, 1],
        colors=2,
        solver_name="glucose4",
    )

    assert coloring is None
    assert stats["status"] == "UNKNOWN_FULL_UNSAT_UNCERTIFIED"
