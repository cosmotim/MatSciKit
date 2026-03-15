"""Tests for LFA-DSC thermal conductivity calculation."""

from __future__ import annotations

import numpy as np

from MatSciKit_COSMOTIM.thermal_conductivity import lfa_dsc


class TestLFADSC:
    """Test LFA-DSC thermal conductivity calculation."""

    def test_simple_case(self):
        cp_T = np.array([300, 350, 400, 450, 500], dtype=float)
        cp = np.array([0.5, 0.55, 0.6, 0.62, 0.63], dtype=float)
        cp_err = cp * 0.02

        diff_T = np.array([325, 400, 475], dtype=float)
        diff = np.array([0.8, 0.9, 1.0], dtype=float)
        diff_err = diff * 0.05

        tc_data, cp_avg = lfa_dsc.calculate(cp_T, cp, cp_err, diff_T, diff, diff_err, 6.0, 0.1)

        assert tc_data.shape == (3, 3)
        assert cp_avg.shape == (3, 3)
        assert np.all(tc_data[:, 1] > 0), "κ should be positive"
        assert np.all(tc_data[:, 2] > 0), "κ error should be positive"
