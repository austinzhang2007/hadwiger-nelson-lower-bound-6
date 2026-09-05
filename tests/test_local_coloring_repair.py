from local_coloring_repair import conflict_edges, repair_coloring_locally


def test_conflict_edges_normalizes_and_finds_monochromatic_edges():
    edges = {(1, 0), (1, 2), (2, 3)}
    coloring = [0, 0, 1, 2]

    assert conflict_edges(edges, coloring) == ((0, 1),)


def test_radius_zero_repairs_a_single_conflicting_edge():
    edges = {(0, 1), (1, 2), (2, 3)}
    initial = [0, 0, 1, 0]

    coloring, stats = repair_coloring_locally(
        4,
        edges,
        initial,
        colors=2,
        radius=0,
        solver_name="glucose4",
    )

    assert coloring is not None
    assert all(coloring[left] != coloring[right] for left, right in edges)
    assert stats["status"] == "SAT"
    assert stats["initial_conflicts"] == 1
    assert stats["final_conflicts"] == 0
    assert stats["active_vertices"] == 2


def test_frozen_boundary_can_make_local_problem_unsatisfiable():
    # The two conflict endpoints see opposite frozen colors, so neither can
    # change at radius zero.  Expanding to radius one admits a repair.
    edges = {(0, 1), (0, 2), (1, 3)}
    initial = [0, 0, 1, 1]

    coloring0, stats0 = repair_coloring_locally(
        4,
        edges,
        initial,
        colors=2,
        radius=0,
        solver_name="glucose4",
    )
    coloring1, stats1 = repair_coloring_locally(
        4,
        edges,
        initial,
        colors=2,
        radius=1,
        solver_name="glucose4",
    )

    assert coloring0 is None
    assert stats0["status"] == "UNKNOWN_LOCAL_UNSAT"
    assert coloring1 is not None
    assert stats1["status"] == "SAT"


def test_already_valid_coloring_is_returned_without_sat_call():
    coloring, stats = repair_coloring_locally(
        3,
        {(0, 1), (1, 2)},
        [0, 1, 0],
        colors=2,
        radius=0,
    )

    assert coloring == (0, 1, 0)
    assert stats["status"] == "SAT"
    assert stats["active_vertices"] == 0
