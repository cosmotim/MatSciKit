"""
Cahill minimum thermal conductivity model.

Computes the theoretical minimum thermal conductivity based on the
Cahill-Watson-Pohl model for amorphous solids.

Translated from calculate_cahill_TC.m
"""

from __future__ import annotations

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


def minimum_tc_high_T(n_density: float,
                      v_L: float,
                      v_T: float) -> float:
    """
    Cahill minimum thermal conductivity in the high-temperature limit.

    Uses the simplified Cheng et al. formula valid when T >> θ_D,
    where the Debye integrals saturate to their classical limit.

    Parameters
    ----------
    n_density : float
        Number density N/V (atoms/m³).
    v_L : float
        Longitudinal sound velocity (m/s).
    v_T : float
        Transverse sound velocity (m/s).

    Returns
    -------
    kappa_min : float
        High-temperature minimum thermal conductivity (W/(m·K)).

    Notes
    -----
    Formula:

    .. math::

        \\kappa_{\\min} = \\left(\\frac{\\pi}{48}\\right)^{1/3}
        k_B \\, n^{2/3} (v_L + 2v_T)

    The coefficient (π/48)^(1/3) relates to the full formula as:

    .. math::

        \\left(\\frac{\\pi}{48}\\right)^{1/3} =
        \\left(\\frac{\\pi}{6}\\right)^{1/3} \\times \\frac{1}{2}

    The factor of 1/2 arises from the high-T limit of the Debye integral:
    as T → ∞, ∫₀^(θ/T) x³eˣ/(eˣ-1)² dx → (θ/T)²/2, which cancels the
    (T/θ)² prefactor and leaves a factor of 1/2.

    This is the formula used in Cheng et al., Small 17, 2101693 (2021).
    It treats the three acoustic branches separately (1 longitudinal +
    2 transverse) rather than using a single average velocity.

    References
    ----------
    Cahill, D. G., Watson, S. K., & Pohl, R. O., Phys. Rev. B 46, 6131 (1992).
    Cheng, Z. et al., Small 17, 2101693 (2021).
    """
    return (np.pi / 48) ** (1.0 / 3) * kb * n_density ** (2.0 / 3) * (v_L + 2 * v_T)
