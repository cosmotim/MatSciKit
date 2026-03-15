"""Tests for Debye temperature converters."""

from __future__ import annotations

import numpy as np
import pytest
from conftest import DENSITY_LSHT, N_DENSITY_LSHT

from MatSciKit_COSMOTIM.constants import hbar, kb
from MatSciKit_COSMOTIM.heat_capacity import debye


class TestDebyeConverters:
    """Test Debye temperature converters."""

    def test_velocity_roundtrip(self):
        """Test that velocity → θ_D → velocity roundtrips."""
        v_s_input = 3000.0

        theta_D = debye.from_velocity(v_s_input, N_DENSITY_LSHT)
        assert theta_D > 0

        v_s_back = theta_D / (hbar / kb * (6 * np.pi**2 * N_DENSITY_LSHT) ** (1.0 / 3))
        assert v_s_back == pytest.approx(v_s_input, rel=1e-6)

    def test_modulus_positive(self):
        theta_D = debye.from_modulus(100.0, DENSITY_LSHT, N_DENSITY_LSHT)
        assert theta_D > 0
