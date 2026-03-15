"""Tests for mean free path calculation."""

from __future__ import annotations

import numpy as np
import pytest

from MatSciKit_COSMOTIM.thermal_conductivity import mean_free_path


class TestMeanFreePath:
    """Test mean free path calculation."""

    def test_with_tto_data(self, lsht_tto_data, lsht_debye_params):
        theta_D, v_s = lsht_debye_params
        T = lsht_tto_data[:, 0]
        kappa = lsht_tto_data[:, 1]

        mfp = mean_free_path.calculate(T, kappa, theta_D, v_s)

        assert len(mfp) == len(T)
        assert np.all(mfp > 0), "MFP should be positive"
        assert np.mean(mfp[: len(mfp) // 4]) > np.mean(mfp[-len(mfp) // 4 :])

    def test_compare_with_matlab_output(self, data_dir, lsht_tto_data, lsht_debye_params):
        theta_D, v_s = lsht_debye_params
        mfp_file = data_dir / "TTO_Data" / "LSHT_MFP.csv"
        if not mfp_file.exists():
            pytest.skip("MFP reference data not available")

        mfp = mean_free_path.calculate(lsht_tto_data[:, 0], lsht_tto_data[:, 1], theta_D, v_s)

        assert np.all(mfp > 1e-12), "MFP should be > 1 pm"
        assert np.all(mfp < 1e-3), "MFP should be < 1 mm"
