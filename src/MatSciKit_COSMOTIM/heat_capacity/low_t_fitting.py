"""
Low-temperature heat capacity fitting for Debye temperature and sound velocity.

Fits Cp/T vs T² in the low-temperature regime to extract the Debye temperature
(θ_D) and average sound velocity (v_s).

Translated from LowT_Cp_fitting.m
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

from MatSciKit_COSMOTIM.constants import hbar, kb


def _linear(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Linear model: y = a*x + b."""
    return a * x + b


def fit(
    T: np.ndarray,
    Cp: np.ndarray,
    Cp_err: np.ndarray,
    n_density: float,
    density: float,
    t_range: tuple[float, float] | None = None,
    n_range: tuple[int, int] | None = None,
) -> tuple[float, float, float, float]:
    """
    Extract Debye temperature and sound velocity from low-T Cp data.

    Performs a weighted linear fit of Cp/T vs T² in the low-temperature
    regime. The slope relates to the Debye temperature through:

        Cp/T = β·T² + γ

    where β = (12π⁴/5) · N_density · kb / (Density · θ_D³)

    The fitting region can be specified in two ways:

    - **Temperature range:** ``t_range=(T_min, T_max)`` selects all points
      within that range.
    - **Index range:** ``n_range=(n_start, n_end)`` uses data points at
      indices n_start through n_end (1-based, inclusive) after sorting by T.
      This matches the MATLAB convention ``data(n_start:n_end, :)``.

    If both are given, ``t_range`` takes priority. If neither is given,
    all data is used.

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
    t_range : tuple of (float, float), optional
        Temperature range (T_min, T_max) in Kelvin for the fitting region.
    n_range : tuple of (int, int), optional
        Index range (n_start, n_end) using 1-based inclusive indexing,
        applied after sorting by temperature. Matches MATLAB convention
        ``data(n_start:n_end, :)``. For example, ``n_range=(13, 41)``
        selects the 13th through 41st data points.

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

    Examples
    --------
    >>> # Method 1: Fit using temperature range
    >>> theta_D, v_s, err_D, err_v = fit(T, Cp, Cp_err, n_density, density,
    ...                                   t_range=(3.0, 10.0))
    >>> # Method 2: Use data points 13 through 41 (1-based, like MATLAB)
    >>> theta_D, v_s, err_D, err_v = fit(T, Cp, Cp_err, n_density, density,
    ...                                   n_range=(13, 41))
    """
    T = np.asarray(T, dtype=float)
    Cp = np.asarray(Cp, dtype=float)
    Cp_err = np.asarray(Cp_err, dtype=float)

    # Sort by temperature
    sort_idx = np.argsort(T)
    T = T[sort_idx]
    Cp = Cp[sort_idx]
    Cp_err = Cp_err[sort_idx]

    # Select fitting region
    if t_range is not None:
        # Method 1: Temperature range
        t_min, t_max = t_range
        mask = (t_min <= T) & (t_max >= T)
        T = T[mask]
        Cp = Cp[mask]
        Cp_err = Cp_err[mask]
    elif n_range is not None:
        # Method 2: Index range (1-based inclusive, like MATLAB)
        n_start, n_end = n_range
        # Convert 1-based inclusive to 0-based Python slice
        T = T[n_start - 1 : n_end]
        Cp = Cp[n_start - 1 : n_end]
        Cp_err = Cp_err[n_start - 1 : n_end]

    if len(T) < 3:
        raise ValueError(
            f"Fewer than 3 data points in the selected range. "
            f"Need at least 3 for a meaningful linear fit. Got {len(T)}."
        )

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
