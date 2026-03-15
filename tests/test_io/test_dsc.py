"""Tests for DSC data reader."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.io import dsc


class TestDSCReader:
    """Test DSC data reader."""

    def test_read_dsc_file(self, data_dir):
        dsc_file = data_dir / "ExpDat_LSHT-05232022-cool.csv"
        if not dsc_file.exists():
            pytest.skip("DSC test data not available")

        data = dsc.read(str(dsc_file))

        assert data.shape[1] == 2
        assert data[0, 0] > 270, "Temperature should be in Kelvin"
        assert data[0, 0] < 320, "First temperature should be near room temp"
        assert np.all(data[:, 1] > 0), "Cp should be positive"
        assert np.all(data[:, 1] < 2), "Cp should be < 2 J/(g·K)"
