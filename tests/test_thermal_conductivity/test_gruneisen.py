"""Tests for Grüneisen parameter calculation."""

from __future__ import annotations

from MatSciKit_COSMOTIM.thermal_conductivity import gruneisen


class TestGruneisen:
    """Test Grüneisen parameter calculation."""

    def test_nzp_example(self):
        gamma, gamma_err = gruneisen.calculate(
            thermal_expansion=2 * 4.2e-6 + 9.8e-6,
            bulk_modulus=87.7e9,
            cp=549.2,
            density=3.254e3,
            thermal_expansion_r_error=0.1,
            bulk_modulus_r_error=5.9 / 87.7,
            cp_r_error=11.3 / 549.2,
        )
        assert gamma > 0
        assert gamma_err > 0
        assert gamma_err < gamma

    def test_no_errors(self):
        gamma, gamma_err = gruneisen.calculate(
            thermal_expansion=1e-5,
            bulk_modulus=100e9,
            cp=500,
            density=5000,
        )
        assert gamma > 0
        assert gamma_err == 0.0
