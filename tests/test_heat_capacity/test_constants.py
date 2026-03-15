"""Tests for physical constants."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.constants import h, hbar, kb


class TestConstants:
    """Test physical constants."""

    def test_boltzmann(self):
        assert kb == pytest.approx(1.380649e-23, rel=1e-6)

    def test_planck(self):
        assert h == pytest.approx(6.62607015e-34, rel=1e-6)

    def test_hbar(self):
        assert hbar == pytest.approx(h / (2 * np.pi), rel=1e-10)
