import pytest

from primary_field_audit import PrimaryField


@pytest.mark.parametrize("value", [0.5, "0.5", "1e-6", True])
def test_field_coefficients_require_explicit_exact_rationals(value):
    with pytest.raises(ValueError, match="rational"):
        PrimaryField.element([value] + [0] * 31)
