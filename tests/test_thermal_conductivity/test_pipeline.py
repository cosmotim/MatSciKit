"""Tests for Pipeline 3: thermal conductivity pipeline functions."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.thermal_conductivity import pipeline


class TestTTOFromPPMS:
    """Test PPMS TTO pipeline."""

    def test_read_tto_csv(self, data_dir):
        """Test reading pre-processed TTO CSV data."""
        tto_file = data_dir / "TTO_Data" / "LSHT_PPMS_surface.csv"
        if not tto_file.exists():
            pytest.skip("TTO test data not available")

        # Pre-processed CSV has 3 columns [T, κ, κ_err] — read directly
        data = np.loadtxt(str(tto_file), delimiter=",")
        assert data.shape[1] == 3
        assert np.all(data[:, 0] > 0)
        assert np.all(data[:, 1] > 0)


class TestKappaFromLFA:
    """Test LFA → κ conversion pipeline."""

    def test_synthetic_data(self):
        """Test with synthetic Cp and diffusivity data."""
        from MatSciKit_COSMOTIM.thermal_conductivity import lfa_dsc

        # Synthetic DSC data
        cp_T = np.linspace(300, 800, 50)
        cp = 0.5 + 0.0003 * (cp_T - 300)  # J/(g·K), linearly increasing
        cp_err = cp * 0.02

        # Synthetic LFA data
        diff_T = np.array([325, 425, 525, 625, 725])
        diff = np.array([0.8, 0.85, 0.9, 0.88, 0.86])  # mm²/s
        diff_err = diff * 0.05

        density = 6.87
        density_err = density * 0.041

        tc_data, cp_avg = lfa_dsc.calculate(
            cp_T, cp, cp_err, diff_T, diff, diff_err, density, density_err
        )

        assert tc_data.shape == (5, 3)
        assert np.all(tc_data[:, 1] > 0), "κ should be positive"
        assert np.all(tc_data[:, 2] > 0), "κ error should be positive"
        assert cp_avg.shape == (5, 3)

    def test_with_real_data(self, data_dir):
        """Test LFA→κ with real LSHT data files."""
        dsc_file = data_dir / "ExpDat_LSHT-05232022-cool.csv"
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"

        if not dsc_file.exists() or not lfa_file.exists():
            pytest.skip("DSC/LFA test data not available")

        result = pipeline.kappa_from_lfa(
            cp_file=str(dsc_file),
            lfa_file=str(lfa_file),
            density=6.87,
            density_error=6.87 * 0.041,
            cp_format="dsc",
            lfa_format="csv",
        )

        assert len(result["temperature"]) > 0
        assert np.all(result["kappa"] > 0)
        assert np.all(result["kappa_error"] > 0)
        assert result["kappa_solid"] is None  # no porosity correction
        assert result["cp_averaged"].shape[1] == 3

    def test_with_porosity(self, data_dir):
        """Test LFA→κ with porosity correction."""
        dsc_file = data_dir / "ExpDat_LSHT-05232022-cool.csv"
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"

        if not dsc_file.exists() or not lfa_file.exists():
            pytest.skip("DSC/LFA test data not available")

        result = pipeline.kappa_from_lfa(
            cp_file=str(dsc_file),
            lfa_file=str(lfa_file),
            density=6.87,
            density_error=6.87 * 0.041,
            porosity=0.03,
        )

        assert result["kappa_solid"] is not None
        assert np.all(result["kappa_solid"] > result["kappa"])


class TestCompareTTOLFA:
    """Test TTO vs LFA comparison."""

    def test_comparison(self):
        """Test comparison with synthetic overlapping data."""
        tto_result = {
            "temperature": np.linspace(5, 500, 50),
            "kappa": np.linspace(0.2, 1.8, 50),
            "kappa_error": np.ones(50) * 0.05,
        }
        lfa_result = {
            "temperature": np.array([300, 400, 500, 600, 700]),
            "kappa": np.array([1.6, 1.65, 1.7, 1.68, 1.65]),
            "kappa_error": np.array([0.1, 0.1, 0.1, 0.1, 0.1]),
        }

        comp = pipeline.compare_tto_lfa(tto_result, lfa_result)

        assert len(comp["tto_temperature"]) == 50
        assert len(comp["lfa_temperature"]) == 5
        # T=300, 400, 500 overlap (TTO goes up to 500)
        assert len(comp["overlap_temps"]) == 3
        assert np.all(comp["ratio"] > 0)

    def test_no_overlap(self):
        """Test when temperature ranges don't overlap."""
        tto_result = {
            "temperature": np.linspace(5, 50, 10),
            "kappa": np.ones(10),
        }
        lfa_result = {
            "temperature": np.linspace(300, 800, 5),
            "kappa": np.ones(5),
        }

        comp = pipeline.compare_tto_lfa(tto_result, lfa_result)
        assert len(comp["overlap_temps"]) == 0
        assert len(comp["ratio"]) == 0
