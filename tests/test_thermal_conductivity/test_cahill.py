"""Tests for Cahill minimum thermal conductivity."""

from __future__ import annotations

import numpy as np
from conftest import N_DENSITY_LSHT

from MatSciKit_COSMOTIM.thermal_conductivity import cahill


class TestCahillMinimum:
    """Test Cahill minimum thermal conductivity."""

    def test_positive_values(self, lsht_debye_params):
        theta_D, v_s = lsht_debye_params
        T = np.linspace(10, 300, 50)
        kappa_min = cahill.minimum_tc(T, N_DENSITY_LSHT, theta_D, v_s)
        assert len(kappa_min) == 50
        assert np.all(kappa_min > 0)

    def test_increases_with_T(self, lsht_debye_params):
        theta_D, v_s = lsht_debye_params
        T = np.array([10, 50, 100, 200, 300])
        kappa_min = cahill.minimum_tc(T, N_DENSITY_LSHT, theta_D, v_s)
        assert np.all(np.diff(kappa_min) > 0)

    def test_scalar_input(self, lsht_debye_params):
        theta_D, v_s = lsht_debye_params
        result = cahill.minimum_tc(300.0, N_DENSITY_LSHT, theta_D, v_s)
        assert isinstance(result, float)
        assert result > 0

    def test_reasonable_magnitude(self, lsht_debye_params):
        theta_D, v_s = lsht_debye_params
        kappa_300 = cahill.minimum_tc(300.0, N_DENSITY_LSHT, theta_D, v_s)
        assert 0.1 < kappa_300 < 5.0, f"κ_min(300K) = {kappa_300:.3f} W/(m·K)"
