"""
Low-temperature heat capacity fitting for Debye temperature and sound velocity.

Fits Cp/T vs T² in the low-temperature regime to extract the Debye temperature
(θ_D) and average sound velocity (v_s).

Translated from LowT_Cp_fitting.m
"""

import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple

from MatSciKit_COSMOTIM.constants import kb, h, hbar


def _linear(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Linear model: y = a*x + b."""
    return a * x + b


def fit(T: np.ndarray, Cp: np.ndarray, Cp_err: np.ndarray,
        n_density: float, density: float
        ) -> Tuple[float, float, float, float]:
    """
    Extract Debye temperature and sound velocity from low-T Cp data.

    Performs a weighted linear fit of Cp/T vs T² in the low-temperature
    regime. The slope relates to the Debye temperature through:

        Cp/T = β·T² + γ

    where β = (12π⁴/5) · N_density · kb / (Density · θ_D³)

    Parameters
    ----------
    T : np.ndarray
        Temperature values (K).
    Cp : np.ndarray
        Heat capacity values (J/(g·K)).
    Cp_err : np.ndarray
        Heat capacity errors (J/(g·K)).
    n_density : float
        Number density N/V (atoms/m³).
    density : float
        Mass density (kg/m³).

    Returns
    -------
    theta_D : float
        Debye temperature (K).
    v_s : float
        Average sound velocity (m/s).
    theta_D_error : float
        Uncertainty in Debye temperature (K).
    v_s_error : float
        Uncertainty in sound velocity (m/s).

    Notes
    -----
    The fit is weighted by (Cp_err/Cp)^(-2), matching the MATLAB implementation
    which uses relative errors as weights for the ``poly1`` fit.
    """
    T = np.asarray(T, dtype=float)
    Cp = np.asarray(Cp, dtype=float)
    Cp_err = np.asarray(Cp_err, dtype=float)

    # Prepare data: Cp/T vs T²
    x = T**2
    y = Cp / T

    # Weights: inverse of relative error squared, matching MATLAB W = (Cp_err./Cp).^(-2)
    sigma = Cp_err / Cp  # relative errors used as sigma for curve_fit

    # Weighted linear fit
    popt, pcov = curve_fit(_linear, x, y, sigma=sigma, absolute_sigma=False)
    slope = popt[0]

    # 95% confidence interval for slope error
    slope_std = np.sqrt(pcov[0, 0])
    slope_95 = slope_std * 1.96  # approximate 95% CI
    r_error = slope_95 / abs(slope)

    # Debye temperature: θ_D = (slope * Density * 1e3 / (12π⁴/5 * N_density * kb))^(-1/3)
    theta_D = (1.0 / (12 * np.pi**4 / 5 * n_density * kb) * slope * density * 1e3) ** (-1.0 / 3)
    theta_D_error = theta_D * r_error / 3

    # Sound velocity: v_s = θ_D / (ħ/kb * (6·N_density·π²)^(1/3))
    v_s = theta_D / (hbar / kb * (6 * n_density * np.pi**2) ** (1.0 / 3))
    v_s_error = v_s * r_error / 3

    return theta_D, v_s, theta_D_error, v_s_error
