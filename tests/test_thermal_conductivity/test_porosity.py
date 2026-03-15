"""Tests for porosity correction."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.thermal_conductivity import porosity_correction


class TestPorosityCorrection:
    """Test porosity correction."""

    def test_zero_porosity(self):
        assert porosity_correction.correct(1.0, 0.0) == pytest.approx(1.0)

    def test_known_value(self):
        expected = 2.03 / 1.94
        result = porosity_correction.correct(1.0, 0.03)
        assert result == pytest.approx(expected, rel=1e-6)

    def test_array_input(self):
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
