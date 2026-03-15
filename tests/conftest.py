"""Shared test fixtures for MatSciKit."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# Path to MATLAB prototype test data
DATA_DIR = Path(__file__).parent.parent.parent / "MatSciKit_Matlab_Prototype"

# LSHT material parameters (from MATLAB main script)
V_LSHT = 3.98**3 * 1e-30  # unit cell volume (m³)
M_LSHT = (
    (6.4 * 3 / 8 + 87.62 * 7 / 16 + 178.49 * 1 / 4 + 180.95 * 3 / 4 + 16 * 3) * (1.66e-24) * 1e-3
)  # kg
DENSITY_LSHT = M_LSHT / V_LSHT  # kg/m³
N_LSHT = 3 / 8 + 7 / 16 + 1 / 4 + 3 / 4 + 3  # atoms per unit cell
N_DENSITY_LSHT = N_LSHT / V_LSHT  # atoms/m³
POROSITY_LSHT = 0.03


@pytest.fixture
def lsht_params():
    """LSHT material parameters dict."""
    return {
        "volume": V_LSHT,
        "mass": M_LSHT,
        "density": DENSITY_LSHT,
        "n_atoms": N_LSHT,
        "n_density": N_DENSITY_LSHT,
        "porosity": POROSITY_LSHT,
    }


@pytest.fixture
def data_dir():
    """Path to MATLAB prototype data directory."""
    return DATA_DIR


@pytest.fixture
def lsht_cp_data(data_dir):
    """Load LSHT Cp data, skip if unavailable."""
    cp_file = data_dir / "TTO_Data" / "LSHT_Cp_all.csv"
    if not cp_file.exists():
        pytest.skip("Cp test data not available")
    data = np.loadtxt(str(cp_file), delimiter=",")
    idx = np.argsort(data[:, 0])
    return data[idx]


@pytest.fixture
def lsht_tto_data(data_dir):
    """Load LSHT TTO data, skip if unavailable."""
    tto_file = data_dir / "TTO_Data" / "LSHT_PPMS_surface.csv"
    if not tto_file.exists():
        pytest.skip("TTO test data not available")
    return np.loadtxt(str(tto_file), delimiter=",")


@pytest.fixture
def lsht_debye_params(lsht_cp_data, lsht_params):
    """Run Pipeline 1 to get θ_D and v_s."""
    from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting

    theta_D, v_s, _, _ = low_t_fitting.fit(
        lsht_cp_data[:, 0],
        lsht_cp_data[:, 1],
        lsht_cp_data[:, 2],
        lsht_params["n_density"],
        lsht_params["density"],
        n_range=(13, 41),
    )
    return theta_D, v_s
