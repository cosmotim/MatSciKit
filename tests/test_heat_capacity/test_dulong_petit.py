"""Tests for Dulong-Petit limit calculation."""

from __future__ import annotations

from MatSciKit_COSMOTIM.heat_capacity import dulong_petit


class TestDulongPetit:
    """Test Dulong-Petit limit calculation."""

    def test_lsht_limit(self, lsht_params):
        dp = dulong_petit.calculate(lsht_params["n_density"], lsht_params["density"])
        assert dp > 0
        assert 0.1 < dp < 2.0, f"Dulong-Petit limit = {dp:.4f} seems unreasonable"
