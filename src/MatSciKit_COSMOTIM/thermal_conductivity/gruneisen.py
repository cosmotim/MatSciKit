"""
Grüneisen parameter calculation.

The Grüneisen parameter relates thermal expansion to the elastic
and thermal properties of a material.

Translated from Gruneisen_calculate.m
"""

from __future__ import annotations

import numpy as np


def calculate(
    thermal_expansion: float,
    bulk_modulus: float,
    cp: float,
    density: float,
    thermal_expansion_r_error: float = 0.0,
    bulk_modulus_r_error: float = 0.0,
    cp_r_error: float = 0.0,
) -> tuple[float, float]:
    """
    Calculate the Grüneisen parameter.

    Parameters
    ----------
    thermal_expansion : float
        Volumetric thermal expansion coefficient (K⁻¹).
    bulk_modulus : float
        Bulk modulus (Pa).
    cp : float
        Specific heat capacity (J/(kg·K)).
    density : float
        Mass density (kg/m³).
    thermal_expansion_r_error : float, optional
        Relative error in thermal expansion. Default is 0.
    bulk_modulus_r_error : float, optional
        Relative error in bulk modulus. Default is 0.
    cp_r_error : float, optional
        Relative error in heat capacity. Default is 0.

    Returns
    -------
    gamma : float
        Grüneisen parameter (dimensionless).
    gamma_error : float
        Uncertainty in Grüneisen parameter.

    Notes
    -----
    Formula: γ = α_V · B / (Cp · ρ)

    where α_V is the volumetric thermal expansion coefficient,
    B is the bulk modulus, Cp is the specific heat, and ρ is the density.
    """
    gamma = thermal_expansion * bulk_modulus / (cp * density)

    gamma_error = (
        np.sqrt(thermal_expansion_r_error**2 + bulk_modulus_r_error**2 + cp_r_error**2) * gamma
    )

    return gamma, gamma_error
