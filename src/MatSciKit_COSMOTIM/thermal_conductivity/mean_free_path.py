"""
Phonon mean free path calculation from thermal conductivity.

Uses kinetic theory to extract the average phonon mean free path
from measured thermal conductivity data.

Translated from MFP_from_TC.m
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from MatSciKit_COSMOTIM.constants import hbar, kb


def _integrand(x: float) -> float:
    """Debye integrand for MFP: x⁴·eˣ / (eˣ - 1)²."""
    if x <= 0:
        return 0.0
    if x > 500:
        return 0.0
    ex = np.exp(x)
    return x**4 * ex / (ex - 1) ** 2


def _debye_integral(T: float, theta_D: float) -> float:
    """
    Evaluate the Debye integral ∫₀^(θ_D/T) x⁴·eˣ/(eˣ-1)² dx.

    Parameters
    ----------
    T : float
        Temperature (K).
    theta_D : float
        Debye temperature (K).

    Returns
    -------
    value : float
        Integral value.
    """
    if T <= 0:
        return 0.0
    upper = theta_D / T
    result, _ = quad(_integrand, 0, upper)
    return result


def calculate(
    T: float | np.ndarray, kappa: float | np.ndarray, theta_D: float, v_s: float
) -> float | np.ndarray:
    """
    Calculate phonon mean free path from thermal conductivity.

    Uses the kinetic theory relation κ = (1/3)·C·v·ℓ, inverted
    with the Debye model for the heat capacity.

    Parameters
    ----------
    T : float or np.ndarray
        Temperature(s) (K).
    kappa : float or np.ndarray
        Thermal conductivity values (W/(m·K)).
    theta_D : float
        Debye temperature (K).
    v_s : float
        Average sound velocity (m/s).

    Returns
    -------
    mfp : float or np.ndarray
        Mean free path (m).

    Notes
    -----
    Formula:
        ℓ = κ / (kb⁴·T³ / (2π²·v_s²·ℏ³) · I(θ_D/T))

    where I = ∫₀^(θ_D/T) x⁴·eˣ/(eˣ-1)² dx.

    Note the integrand here uses x⁴ (not x³ as in Cahill), because
    this comes from the volumetric heat capacity expression in the
    Debye model.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    kappa = np.atleast_1d(np.asarray(kappa, dtype=float))

    mfp = np.zeros(len(T))
    for i in range(len(T)):
        integral_val = _debye_integral(T[i], theta_D)
        # Prefactor: kb⁴·T³ / (2π²·v_s²·ℏ³)
        prefactor = (kb**4 * T[i] ** 3) / (2 * np.pi**2 * v_s**2 * hbar**3)
        mfp[i] = kappa[i] / (prefactor * integral_val)

    if mfp.size == 1:
        return float(mfp[0])
    return mfp
