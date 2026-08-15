from pathlib import Path

from copy_interactions import complete_copy_interactions
from exact_geometry import parse_mathematica_vertices


ROOT = Path(__file__).resolve().parents[1]
SELECTED = [201, 101, 216, 97, 104, 135, 184, 118, 178, 206, 124]


def test_complete_copy_interactions_certifies_all_55_pairs() -> None:
    source = parse_mathematica_vertices(
        ROOT / "third_party/CNP-SAT/vtx/529.vtx"
    )

    result = complete_copy_interactions(source, SELECTED)

    assert len(result.pair_certificates) == 55
    assert len(result.edges) == 1954
    nonempty = {
        (certificate.left_copy, certificate.right_copy): certificate.edge_count
        for certificate in result.pair_certificates
        if certificate.edge_count
    }
    assert nonempty == {(0, 6): 648, (1, 10): 655, (3, 7): 651}
    assert sum(
        certificate.method == "same_extension_basis_separation"
        for certificate in result.pair_certificates
    ) == 25
    assert sum(
        certificate.method == "mixed_extension_product_coefficient"
        for certificate in result.pair_certificates
    ) == 30
    assert all(certificate.complete for certificate in result.pair_certificates)
