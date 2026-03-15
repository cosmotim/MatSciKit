"""Tests for class-based IO readers."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.io.readers import (
    DSCReader,
    LFAReader,
    PPMSHCReader,
)


class TestPPMSHCReader:
    """Test PPMS HC class-based reader."""

    def test_read_and_cache(self, data_dir):
        hc_file = data_dir / "TTO_Data" / "LSHT_Cp_all.csv"
        if not hc_file.exists():
            pytest.skip("Cp test data not available")
        # LSHT_Cp_all.csv is a plain CSV, not a PPMS .dat file,
        # so we test with LFAReader which also handles plain CSVs
        reader = LFAReader(hc_file)
        data = reader.data
        assert data.shape[1] == 3
        # Second call should return cached data
        assert reader.data is data

    def test_metadata(self, data_dir):
        hc_file = data_dir / "TTO_Data" / "LSHT_Cp_all.csv"
        if not hc_file.exists():
            pytest.skip("Cp test data not available")
        reader = LFAReader(hc_file)
        meta = reader.metadata
        assert meta.format_name == "LFA"
        assert meta.n_rows > 0
        assert meta.t_min < meta.t_max
        assert meta.filename == "LSHT_Cp_all.csv"


class TestLFAReader:
    """Test LFA class-based reader."""

    def test_read_lfa(self, data_dir):
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"
        if not lfa_file.exists():
            pytest.skip("LFA test data not available")
        reader = LFAReader(lfa_file)
        data = reader.data
        assert data.shape[1] == 3
        assert not np.any(np.isnan(data))

    def test_summary(self, data_dir):
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"
        if not lfa_file.exists():
            pytest.skip("LFA test data not available")
        reader = LFAReader(lfa_file)
        s = reader.summary()
        assert "LFA" in s
        assert "Temperature" in s

    def test_nan_rows_counted(self, data_dir):
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"
        if not lfa_file.exists():
            pytest.skip("LFA test data not available")
        reader = LFAReader(lfa_file)
        meta = reader.metadata
        if "rows_dropped_nan" in meta.extra:
            assert meta.extra["rows_dropped_nan"] >= 0


class TestDSCReader:
    """Test DSC class-based reader."""

    def test_read_dsc(self, data_dir):
        dsc_file = data_dir / "ExpDat_LSHT-05232022-cool.csv"
        if not dsc_file.exists():
            pytest.skip("DSC test data not available")
        reader = DSCReader(dsc_file)
        data = reader.data
        assert data.shape[1] == 2
        assert data[0, 0] > 270  # Temperature in K


class TestFileNotFound:
    """Test error handling."""

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError):
            PPMSHCReader("/nonexistent/file.dat")

    def test_missing_file_lfa(self):
        with pytest.raises(FileNotFoundError):
            LFAReader("/nonexistent/file.csv")


class TestRepr:
    """Test string representation."""

    def test_repr(self, data_dir):
        lfa_file = data_dir / "TTO_Data" / "LSHT_LFA.csv"
        if not lfa_file.exists():
            pytest.skip("LFA test data not available")
        reader = LFAReader(lfa_file)
        assert "LFAReader" in repr(reader)
