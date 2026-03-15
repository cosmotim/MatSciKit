"""Tests for low-temperature Cp fitting."""

from __future__ import annotations


class TestLowTFitting:
    """Test low-temperature Cp fitting for Debye temperature."""

    def test_with_lsht_data(self, lsht_cp_data, lsht_params):
        """Test against LSHT sample data from MATLAB prototype."""
        from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting

        T = lsht_cp_data[:, 0]
        Cp = lsht_cp_data[:, 1]
        Cp_err = lsht_cp_data[:, 2]

        theta_D, v_s, theta_D_err, v_s_err = low_t_fitting.fit(
            T,
            Cp,
            Cp_err,
            lsht_params["n_density"],
            lsht_params["density"],
            t_range=(3.0, 10.0),
        )

        assert 200 < theta_D < 800, f"θ_D = {theta_D:.1f} K seems unreasonable"
        assert 1000 < v_s < 6000, f"v_s = {v_s:.1f} m/s seems unreasonable"
        assert theta_D_err > 0
        assert v_s_err > 0
        assert theta_D_err < theta_D
        assert v_s_err < v_s
