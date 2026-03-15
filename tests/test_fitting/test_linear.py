"""Tests for linear fitting utilities."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.fitting.linear import fit_with_errors


class TestLinearFit:
    """Test weighted linear fit."""

    def test_perfect_line(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        y = 2.0 * x + 1.0
        y_err = np.ones_like(y) * 0.1

        slope, _se_slope, intercept = fit_with_errors(x, y, y_err)

        assert slope == pytest.approx(2.0, rel=1e-6)
        assert intercept == pytest.approx(1.0 + 2.0 * np.mean(x), rel=1e-6)

    def test_noisy_data(self):
        np.random.seed(42)
        x = np.linspace(0, 10, 50)
        y = 3.0 * x + 5.0 + np.random.normal(0, 0.5, 50)
        y_err = np.ones_like(y) * 0.5

        slope, se_slope, _intercept = fit_with_errors(x, y, y_err)

        assert slope == pytest.approx(3.0, abs=0.5)
        assert se_slope > 0
