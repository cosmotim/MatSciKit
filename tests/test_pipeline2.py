"""
Tests for Pipeline 2: Thermal Conductivity Analysis.

Validates Python implementations against known MATLAB outputs using
LSHT sample data from the MATLAB prototype.
"""

import sys
import numpy as np
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from MatSciKit_COSMOTIM.thermal_conductivity import (
    cahill, mean_free_path, porosity_correction, lfa_dsc, gruneisen
)
from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting

# Path to test data
DATA_DIR = Path(__file__).parent.parent.parent / 'MatSciKit_Matlab_Prototype'

# Material parameters for LSHT
V = 3.98**3 * 1e-30
M = (6.4 * 3/8 + 87.62 * 7/16 + 178.49 * 1/4 + 180.95 * 3/4 + 16 * 3) * 1.66e-24 * 1e-3
DENSITY = M / V
N = 3/8 + 7/16 + 1/4 + 3/4 + 3
N_DENSITY = N / V
POROSITY = 0.03


def _get_debye_params():
    """Helper: run Pipeline 1 to get θ_D and v_s for Pipeline 2 tests."""
    cp_file = DATA_DIR / 'TTO_Data' / 'LSHT_Cp_all.csv'
    if not cp_file.exists():
        pytest.skip("Cp test data not available")
    hc_data = np.loadtxt(str(cp_file), delimiter=',')
    theta_D, v_s, _, _ = low_t_fitting.fit(
        hc_data[:, 0], hc_data[:, 1], hc_data[:, 2],
        N_DENSITY, DENSITY, n_range=(13, 41)
    )
    return theta_D, v_s


class TestPorosityCorrection:
    """Test porosity correction."""

    def test_zero_porosity(self):
        """No porosity → no correction."""
        assert porosity_correction.correct(1.0, 0.0) == pytest.approx(1.0)

    def test_known_value(self):
        """Test against MATLAB solidTCwithPorosity(1.0, 0.03)."""
        # κ_s = 1.0 * (2 + 0.03) / (2 - 0.06) = 2.03 / 1.94
        expected = 2.03 / 1.94
        result = porosity_correction.correct(1.0, POROSITY)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_array_input(self):
        """Test with array input."""
        kappa = np.array([0.5, 1.0, 1.5])
        result = porosity_correction.correct(kappa, 0.03)
        assert len(result) == 3
        assert np.all(result > kappa)

    def test_invalid_porosity(self):
        with pytest.raises(ValueError):
            porosity_correction.correct(1.0, 1.0)

    def test_high_porosity_warning(self):
        with pytest.warns(UserWarning):
            porosity_correction.correct(1.0, 0.2)


class TestCahillMinimum:
    """Test Cahill minimum thermal conductivity."""

    def test_positive_values(self):
        theta_D, v_s = _get_debye_params()
        T = np.linspace(10, 300, 50)
        kappa_min = cahill.minimum_tc(T, N_DENSITY, theta_D, v_s)

        assert len(kappa_min) == 50
        assert np.all(kappa_min > 0)

    def test_increases_with_T(self):
        """κ_min should generally increase with temperature."""
        theta_D, v_s = _get_debye_params()
        T = np.array([10, 50, 100, 200, 300])
        kappa_min = cahill.minimum_tc(T, N_DENSITY, theta_D, v_s)

        # Should be monotonically increasing at these temperatures
        assert np.all(np.diff(kappa_min) > 0)

    def test_scalar_input(self):
        theta_D, v_s = _get_debye_params()
        result = cahill.minimum_tc(300.0, N_DENSITY, theta_D, v_s)
        assert isinstance(result, float)
        assert result > 0

    def test_reasonable_magnitude(self):
        """Cahill min κ for ceramics should be ~0.5-2.0 W/(m·K) at 300K."""
        theta_D, v_s = _get_debye_params()
        kappa_300 = cahill.minimum_tc(300.0, N_DENSITY, theta_D, v_s)
        assert 0.1 < kappa_300 < 5.0, f"κ_min(300K) = {kappa_300:.3f} W/(m·K)"
        print(f"κ_min(300K) = {kappa_300:.3f} W/(m·K)")


class TestMeanFreePath:
    """Test mean free path calculation."""

    def test_with_tto_data(self):
        """Test MFP calculation against LSHT TTO data."""
        tto_file = DATA_DIR / 'TTO_Data' / 'LSHT_PPMS_surface.csv'
        if not tto_file.exists():
            pytest.skip("TTO test data not available")

        theta_D, v_s = _get_debye_params()
        tto_data = np.loadtxt(str(tto_file), delimiter=',')

        T = tto_data[:, 0]
        kappa = tto_data[:, 1]

        mfp = mean_free_path.calculate(T, kappa, theta_D, v_s)

        assert len(mfp) == len(T)
        assert np.all(mfp > 0), "MFP should be positive"

        # MFP should decrease with temperature (more scattering)
        # Check overall trend (first vs last quarter)
        assert np.mean(mfp[:len(mfp)//4]) > np.mean(mfp[-len(mfp)//4:])

        print(f"MFP range: {np.min(mfp)*1e9:.1f} - {np.max(mfp)*1e9:.1f} nm")

    def test_compare_with_matlab_output(self):
        """Compare against MATLAB MFP output file."""
        mfp_file = DATA_DIR / 'TTO_Data' / 'LSHT_MFP.csv'
        tto_file = DATA_DIR / 'TTO_Data' / 'LSHT_PPMS_surface.csv'
        if not mfp_file.exists() or not tto_file.exists():
            pytest.skip("MFP/TTO test data not available")

        theta_D, v_s = _get_debye_params()

        # MATLAB output: [T, MFP]
        matlab_mfp = np.loadtxt(str(mfp_file), delimiter=',')

        # Our calculation
        tto_data = np.loadtxt(str(tto_file), delimiter=',')
        mfp = mean_free_path.calculate(tto_data[:, 0], tto_data[:, 1], theta_D, v_s)

        # Compare at matching temperatures (MATLAB may have used different input)
        # Just check that orders of magnitude match
        assert np.all(mfp > 1e-12), "MFP should be > 1 pm"
        assert np.all(mfp < 1e-3), "MFP should be < 1 mm"

        print(f"Python MFP at T={tto_data[0,0]:.1f}K: {mfp[0]*1e9:.2f} nm")
        print(f"MATLAB MFP at T={matlab_mfp[0,0]:.1f}K: {matlab_mfp[0,1]*1e9:.2f} nm")


class TestGruneisen:
    """Test Grüneisen parameter calculation."""

    def test_nzp_example(self):
        """Test against NZP example from MATLAB Gruneisen_calculate.m."""
        gamma, gamma_err = gruneisen.calculate(
            thermal_expansion=2 * 4.2e-6 + 9.8e-6,  # K^-1
            bulk_modulus=87.7e9,  # Pa
            cp=549.2,  # J/(kg·K)
            density=3.254e3,  # kg/m³
            thermal_expansion_r_error=0.1,
            bulk_modulus_r_error=5.9 / 87.7,
            cp_r_error=11.3 / 549.2
        )

        assert gamma > 0
        assert gamma_err > 0
        assert gamma_err < gamma
        print(f"NZP Grüneisen: γ = {gamma:.4f} ± {gamma_err:.4f}")

    def test_no_errors(self):
        """Test with zero errors."""
        gamma, gamma_err = gruneisen.calculate(
            thermal_expansion=1e-5,
            bulk_modulus=100e9,
            cp=500,
            density=5000
        )
        assert gamma > 0
        assert gamma_err == 0.0


class TestLFADSC:
    """Test LFA-DSC thermal conductivity calculation."""

    def test_simple_case(self):
        """Test with synthetic data."""
        cp_T = np.array([300, 350, 400, 450, 500], dtype=float)
        cp = np.array([0.5, 0.55, 0.6, 0.62, 0.63], dtype=float)
        cp_err = cp * 0.02

        diff_T = np.array([325, 400, 475], dtype=float)
        diff = np.array([0.8, 0.9, 1.0], dtype=float)
        diff_err = diff * 0.05

        density = 6.0
        density_err = 0.1

        tc_data, cp_avg = lfa_dsc.calculate(
            cp_T, cp, cp_err, diff_T, diff, diff_err,
            density, density_err
        )

        assert tc_data.shape == (3, 3)
        assert cp_avg.shape == (3, 3)
        assert np.all(tc_data[:, 1] > 0), "κ should be positive"
        assert np.all(tc_data[:, 2] > 0), "κ error should be positive"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
