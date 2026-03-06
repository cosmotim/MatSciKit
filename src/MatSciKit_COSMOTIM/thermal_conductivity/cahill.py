"""
Cahill minimum thermal conductivity model.

Computes the theoretical minimum thermal conductivity based on the
Cahill-Watson-Pohl model for amorphous solids.

Translated from calculate_cahill_TC.m
"""

import numpy as np
from scipy.integrate import quad
from typing import Union

from MatSciKit_COSMOTIM.constants import kb


def _integrand(x: float) -> float:
    """Debye integrand: x³·eˣ / (eˣ - 1)²."""
    if x <= 0:
        return 0.0
    if x > 500:
        # For very large x, eˣ dominates: x³·eˣ/(eˣ)² = x³·e⁻ˣ → 0
        return 0.0
    ex = np.exp(x)
    return x**3 * ex / (ex - 1)**2


def _debye_integral(T: float, theta_D: float) -> float:
    """
    Evaluate the Debye integral from 0 to θ_D/T.

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


def minimum_tc(T: Union[float, np.ndarray],
               n_density: float,
               theta_D: float,
               v_s: float) -> Union[float, np.ndarray]:
    """
    Calculate the Cahill minimum thermal conductivity.

    Parameters
    ----------
    T : float or np.ndarray
        Temperature(s) (K).
    n_density : float
        Number density N/V (atoms/m³).
    theta_D : float
        Debye temperature (K).
    v_s : float
        Average sound velocity (m/s).

    Returns
    -------
    kappa_min : float or np.ndarray
        Minimum thermal conductivity (W/(m·K)).

    Notes
    -----
    Formula:
        κ_min = (π/6)^(1/3) · kb · n^(2/3) · 3·v_s · (T/θ_D)² · I(θ_D/T)

    where I is the Debye integral ∫₀^(θ_D/T) x³·eˣ/(eˣ-1)² dx.

    Reference: Cahill, Watson, and Pohl, Phys. Rev. B 46, 6131 (1992).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))

    prefactor = (np.pi / 6) ** (1.0 / 3) * kb * n_density ** (2.0 / 3) * 3 * v_s

    # Vectorized Debye integral evaluation
    integrals = np.array([_debye_integral(t, theta_D) for t in T])

    kappa_min = prefactor * (T / theta_D) ** 2 * integrals

    # Return scalar if input was scalar
    if kappa_min.size == 1:
        return float(kappa_min[0])
    return kappa_min
