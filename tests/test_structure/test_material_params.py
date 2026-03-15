"""Tests for material parameter sanity checks."""

from __future__ import annotations

from conftest import DENSITY_LSHT, N_DENSITY_LSHT


class TestMaterialParams:
    """Test that material parameters are self-consistent."""

    def test_density_positive(self):
        assert DENSITY_LSHT > 0

    def test_n_density_positive(self):
        assert N_DENSITY_LSHT > 0

    def test_density_reasonable(self):
        assert 3000 < DENSITY_LSHT < 10000, f"Density = {DENSITY_LSHT:.1f} kg/m³"
