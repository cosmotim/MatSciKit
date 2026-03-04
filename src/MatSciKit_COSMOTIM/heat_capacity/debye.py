"""
Debye temperature converters.

Convert between Debye temperature, sound velocity, and elastic modulus.

Translated from modulus2debyeT.m and velocity2debyeT.m
"""

import numpy as np

from MatSciKit_COSMOTIM.constants import h, kb, hbar


def from_velocity(v_s: float, n_density: float) -> float:
    """
    Calculate Debye temperature from sound velocity.

    Parameters
    ----------
    v_s : float
        Sound velocity (m/s).
    n_density : float
        Number density N/V (atoms/m³).

    Returns
    -------
    theta_D : float
        Debye temperature (K).

    Notes
    -----
    Formula: θ_D = v_s · (ħ/kb) · (6π²·N_density)^(1/3)
    """
    return v_s * (hbar / kb) * (6 * np.pi**2 * n_density) ** (1.0 / 3)


def from_modulus(modulus_gpa: float, density: float, n_density: float) -> float:
    """
    Calculate Debye temperature from bulk modulus.

    First computes the sound velocity as v = sqrt(B/ρ), then converts
    to Debye temperature.

    Parameters
    ----------
    modulus_gpa : float
        Bulk modulus (GPa).
    density : float
        Mass density (kg/m³).
    n_density : float
        Number density N/V (atoms/m³).

    Returns
    -------
    theta_D : float
        Debye temperature (K).
    """
    v_s = np.sqrt(modulus_gpa * 1e9 / density)
    return from_velocity(v_s, n_density)
