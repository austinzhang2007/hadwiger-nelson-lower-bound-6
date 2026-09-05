import json

import pytest


def test_audit_checks_every_edge_and_persists_partial_progress(tmp_path):
    from primary_checkpoint_audit import audit_coordinates
    from primary_field_audit import PrimaryField

    f = PrimaryField()
    def point(x):
        return (f.element([x] + [0] * 31), f.zero, f.one)
    output = tmp_path / "audit.json"
    snapshots = []
    result = audit_coordinates(
        [point(0), point(1), point(3)], [(0, 1), (1, 2)], [0, 1, 0],
        output=output, batch_size=1,
        progress=lambda r: snapshots.append(json.loads(output.read_text())),
    )
    assert result["invalid_edges"] == [[1, 2]]
    assert result["checked_edges"] == 2
    assert result["listed_edges_exactly_unit"] is False
    assert snapshots[0]["status"] == "RUNNING"
    assert snapshots[0]["checked_edges"] == 1
    assert snapshots[-1]["status"] == "FINISHED"


def test_audit_rejects_color_coercions_and_zero_denominators(tmp_path):
    from primary_checkpoint_audit import audit_coordinates
    from primary_field_audit import PrimaryField

    f = PrimaryField()
    p = (f.zero, f.zero, f.one)
    for bad in (0.5, True, "0", -1, 5):
        with pytest.raises(ValueError, match="color"):
            audit_coordinates([p], [], [bad], output=tmp_path / "bad.json")
    with pytest.raises(ValueError, match="denominator"):
        audit_coordinates([(f.zero, f.zero, f.zero)], [], [0],
                          output=tmp_path / "bad.json")


def test_sat_audit_does_not_claim_complete_geometry_or_lower_bound(tmp_path):
    from primary_checkpoint_audit import audit_coordinates
    from primary_field_audit import PrimaryField

    f = PrimaryField()
    result = audit_coordinates([(f.zero, f.zero, f.one)], [], [0],
                               output=tmp_path / "audit.json")
    assert result["coloring_valid"] is True
    assert result["listed_edges_exactly_unit"] is True
    assert result["complete_unit_graph_certified"] is False
    assert result["duplicate_free_certified"] is False
    assert result["lower_bound_6_proved"] is False
