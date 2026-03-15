"""Tests for LFA data reader."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.io import lfa


class TestLFAReader:
    """Test LFA data reader."""

    def test_read_lfa_file(self, data_dir):
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"
        if not lfa_file.exists():
            pytest.skip("LFA test data not available")

        data = lfa.read(str(lfa_file))

        assert data.shape[1] == 3
        assert not np.any(np.isnan(data)), "NaN rows should be removed"
        assert data.shape[0] == 12
        assert data[0, 0] > 250
        assert data[-1, 0] < 900
