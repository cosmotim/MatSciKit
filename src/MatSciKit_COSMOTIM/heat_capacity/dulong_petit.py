"""
Dulong-Petit limit calculation.

The classical upper limit of heat capacity per unit mass.
"""

from MatSciKit_COSMOTIM.constants import kb


def calculate(n_density: float, density: float) -> float:
    """
    Calculate the Dulong-Petit heat capacity limit.

    The Dulong-Petit law states that the molar heat capacity of a solid
    approaches 3R per mole of atoms at high temperatures.

    Parameters
    ----------
    n_density : float
        Number density N/V (atoms/m³).
    density : float
        Mass density (kg/m³).

    Returns
    -------
    limit : float
        Dulong-Petit limit in J/(g·K).

    Notes
    -----
    Formula: C_DP = 3 · N_density · kb / Density × 1e-3

    The factor 1e-3 converts from J/(kg·K) to J/(g·K) to match
    typical specific heat capacity units in materials science.
    """
    return 3 * n_density * kb / density * 1e-3
