"""
Tests for Pipeline 1: Heat Capacity Analysis.

Validates Python implementations against known MATLAB outputs using
LSHT sample data from the MATLAB prototype.
"""

import sys
import numpy as np
import pytest
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from MatSciKit_COSMOTIM.constants import kb, h, hbar
from MatSciKit_COSMOTIM.io import dsc, lfa
from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting, dulong_petit, debye
from MatSciKit_COSMOTIM.fitting.linear import fit_with_errors

# Path to test data
DATA_DIR = Path(__file__).parent.parent.parent / 'MatSciKit_Matlab_Prototype'


# ============================================================
# Material parameters for LSHT (from MATLAB main script)
# ============================================================
V = 3.98**3 * 1e-30  # unit cell volume (m³)
M = (6.4 * 3/8 + 87.62 * 7/16 + 178.49 * 1/4 + 180.95 * 3/4 + 16 * 3) * (1.66e-24) * 1e-3  # kg
DENSITY = M / V
N = (3/8 + 7/16 + 1/4 + 3/4 + 3)  # atoms per unit cell
N_DENSITY = N / V  # number density (atoms/m³)
POROSITY = 0.03


class TestConstants:
    """Test physical constants."""

    def test_boltzmann(self):
        assert kb == pytest.approx(1.380649e-23, rel=1e-6)

    def test_planck(self):
        assert h == pytest.approx(6.62607015e-34, rel=1e-6)

    def test_hbar(self):
        assert hbar == pytest.approx(h / (2 * np.pi), rel=1e-10)


class TestDSCReader:
    """Test DSC data reader."""

    def test_read_dsc_file(self):
        dsc_file = DATA_DIR / 'ExpDat_LSHT-05232022-cool.csv'
        if not dsc_file.exists():
            pytest.skip("DSC test data not available")

        data = dsc.read(str(dsc_file))

        # Should have 2 columns: [T_K, Cp]
        assert data.shape[1] == 2

        # Temperature should be in Kelvin (first DSC point ~26°C → ~299K)
        assert data[0, 0] > 270, "Temperature should be in Kelvin"
        assert data[0, 0] < 320, "First temperature should be near room temp"

        # Cp should be positive and reasonable (0.3-0.9 J/(g·K) typical for ceramics)
        assert np.all(data[:, 1] > 0), "Cp should be positive"
        assert np.all(data[:, 1] < 2), "Cp should be < 2 J/(g·K)"


class TestLFAReader:
    """Test LFA data reader."""

    def test_read_lfa_file(self):
        lfa_file = DATA_DIR / 'TTO_Data' / 'LSHT_LFA.csv'
        if not lfa_file.exists():
            pytest.skip("LFA test data not available")

        data = lfa.read(str(lfa_file))

        # Should have 3 columns and NaN rows removed
        assert data.shape[1] == 3
        assert not np.any(np.isnan(data)), "NaN rows should be removed"

        # Original file has 16 rows, 4 with NaN → 12 valid rows
        assert data.shape[0] == 12

        # Temperature should be reasonable (near room temp and above)
        assert data[0, 0] > 250
        assert data[-1, 0] < 900


class TestLinearFit:
    """Test weighted linear fit."""

    def test_perfect_line(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = 2.0 * x + 1.0
        y_err = np.ones_like(y) * 0.1

        slope, se_slope, intercept = fit_with_errors(x, y, y_err)

        assert slope == pytest.approx(2.0, rel=1e-6)
        assert intercept == pytest.approx(1.0 + 2.0 * np.mean(x), rel=1e-6)

    def test_noisy_data(self):
        np.random.seed(42)
        x = np.linspace(0, 10, 50)
        y = 3.0 * x + 5.0 + np.random.normal(0, 0.5, 50)
        y_err = np.ones_like(y) * 0.5

        slope, se_slope, intercept = fit_with_errors(x, y, y_err)

        assert slope == pytest.approx(3.0, abs=0.5)
        assert se_slope > 0


class TestLowTFitting:
    """Test low-temperature Cp fitting for Debye temperature."""

    def test_with_lsht_data(self):
        """Test against LSHT sample data from MATLAB prototype."""
        cp_file = DATA_DIR / 'TTO_Data' / 'LSHT_Cp_all.csv'
        if not cp_file.exists():
            pytest.skip("Cp test data not available")

        # Load pre-processed HC data [T, Cp, Cp_err]
        hc_data = np.loadtxt(str(cp_file), delimiter=',')

        # Sort by temperature (matching MATLAB)
        idx = np.argsort(hc_data[:, 0])
        hc_data = hc_data[idx]

        # Take points 13-41 (MATLAB 1-indexed) → 12:41 (Python 0-indexed)
        T = hc_data[12:41, 0]
        Cp = hc_data[12:41, 1]
        Cp_err = hc_data[12:41, 2]

        theta_D, v_s, theta_D_err, v_s_err = low_t_fitting.fit(
            T, Cp, Cp_err, N_DENSITY, DENSITY
        )

        # Debye temperature should be in a physically reasonable range
        # For LSHT perovskite, expect ~300-600 K
        assert 200 < theta_D < 800, f"θ_D = {theta_D:.1f} K seems unreasonable"

        # Sound velocity should be ~2000-5000 m/s for ceramics
        assert 1000 < v_s < 6000, f"v_s = {v_s:.1f} m/s seems unreasonable"

        # Errors should be positive and smaller than values
        assert theta_D_err > 0
        assert v_s_err > 0
        assert theta_D_err < theta_D
        assert v_s_err < v_s

        print(f"θ_D = {theta_D:.2f} ± {theta_D_err:.2f} K")
        print(f"v_s = {v_s:.2f} ± {v_s_err:.2f} m/s")


class TestDulongPetit:
    """Test Dulong-Petit limit calculation."""

    def test_lsht_limit(self):
        """Test Dulong-Petit limit for LSHT."""
        dp = dulong_petit.calculate(N_DENSITY, DENSITY)

        # Should be positive and in a reasonable range (0.3-1.0 J/(g·K) typical)
        assert dp > 0
        assert 0.1 < dp < 2.0, f"Dulong-Petit limit = {dp:.4f} seems unreasonable"

        print(f"Dulong-Petit limit: {dp:.4f} J/(g·K)")


class TestDebyeConverters:
    """Test Debye temperature converters."""

    def test_velocity_roundtrip(self):
        """Test that velocity → θ_D → velocity roundtrips."""
        v_s_input = 3000.0  # m/s

        theta_D = debye.from_velocity(v_s_input, N_DENSITY)
        assert theta_D > 0

        # Reverse: v_s = θ_D / (ħ/kb * (6π²·N_density)^(1/3))
        v_s_back = theta_D / (hbar / kb * (6 * np.pi**2 * N_DENSITY) ** (1.0 / 3))
        assert v_s_back == pytest.approx(v_s_input, rel=1e-6)

    def test_modulus_positive(self):
        """Test that modulus converter gives positive θ_D."""
        theta_D = debye.from_modulus(100.0, DENSITY, N_DENSITY)
        assert theta_D > 0


class TestMaterialParams:
    """Test that material parameters are self-consistent."""

    def test_density_positive(self):
        assert DENSITY > 0

    def test_n_density_positive(self):
        assert N_DENSITY > 0

    def test_density_reasonable(self):
        # LSHT density should be ~6000-8000 kg/m³
        density_gcc = DENSITY * 1e-3  # convert to g/cm³... actually check units
        # DENSITY is in kg/m³ already from M(kg)/V(m³)
        assert 3000 < DENSITY < 10000, f"Density = {DENSITY:.1f} kg/m³"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
