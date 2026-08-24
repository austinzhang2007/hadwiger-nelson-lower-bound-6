import json
from pathlib import Path

from phase_normalized_kissat_model import (
    validate_phase_normalized_kissat_coloring,
)


def test_unflips_then_decodes_and_checks_all_edges(tmp_path: Path):
    log = tmp_path / "kissat.log"
    mapping = tmp_path / "mapping.json"
    log.write_text(
        "s SATISFIABLE\nv -1 2 3 -4 -5 6 0\n",
        encoding="ascii",
    )
    mapping.write_text(
        json.dumps({"true_variables": [1, 3]}),
        encoding="utf-8",
    )

    coloring, stats = validate_phase_normalized_kissat_coloring(
        log_path=log,
        mapping_path=mapping,
        vertex_count=2,
        edges={(0, 1)},
        colors=3,
    )

    assert coloring == (0, 2)
    assert stats["transformed_model_literals"] == 6
    assert stats["flipped_variables"] == 2
