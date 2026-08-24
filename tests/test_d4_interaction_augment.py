import numpy as np
from scipy.spatial import cKDTree

from d4_interaction_augment import numeric_annulus_pairs_batched


def test_batched_annulus_matches_full_query_pairs() -> None:
    coordinates = np.asarray(
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (0.6, 0.8),
            (2.0, 0.0),
            (4.0, 4.0),
        ]
    )
    tolerance = 1e-8
    possible = cKDTree(coordinates).query_pairs(
        r=1.0 + tolerance,
        output_type="ndarray",
    )
    delta = coordinates[possible[:, 0]] - coordinates[possible[:, 1]]
    squared = np.einsum("ij,ij->i", delta, delta)
    expected = {
        tuple(map(int, edge))
        for edge in possible[np.abs(squared - 1.0) <= 4.0 * tolerance]
    }

    actual, within = numeric_annulus_pairs_batched(
        coordinates,
        tolerance=tolerance,
        chunk_size=2,
    )

    assert set(actual) == expected
    assert within == len(possible)
    assert len(actual) == len(set(actual))
