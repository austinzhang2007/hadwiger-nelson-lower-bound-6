import sympy as sp

from d4_blocker_cegis import _projective_from_payload, _projective_payload
from d4_exact import D4ProjectivePoint


def test_projective_checkpoint_round_trip() -> None:
    zero = tuple([sp.S.Zero] * 32)
    one = (sp.S.One, *([sp.S.Zero] * 31))
    point = D4ProjectivePoint(one, zero, one, 17)

    recovered = _projective_from_payload(
        {
            **_projective_payload(point),
            "fifth_neighbor_zero_based": 17,
        }
    )

    assert recovered == point
