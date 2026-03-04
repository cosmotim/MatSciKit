"""
Porosity correction for thermal conductivity.

Corrects measured thermal conductivity of a porous sample to obtain
the intrinsic solid thermal conductivity.

Translated from solidTCwithPorosity.m
"""

import numpy as np
from typing import Union


def correct(kappa: Union[float, np.ndarray],
            porosity: float) -> Union[float, np.ndarray]:
    """
    Calculate solid thermal conductivity from measured porous value.

    Assumes negligible pore thermal conductivity (κ_pore ≈ 0) and uses
    the Maxwell-Eucken effective medium approximation.

    Parameters
    ----------
    kappa : float or np.ndarray
        Measured thermal conductivity (W/(m·K)).
    porosity : float
        Volume fraction of pores (0 < φ < 1).

    Returns
    -------
    kappa_s : float or np.ndarray
        Intrinsic solid thermal conductivity (W/(m·K)).

    Raises
    ------
    ValueError
        If porosity is not in (0, 1).

    Notes
    -----
    Formula: κ_s = κ × (2 + φ) / (2 − 2φ)

    This is valid for dilute porosity (φ < ~0.15). A warning is issued
    for higher porosity values.
    """
    if porosity >= 1 or porosity < 0:
        raise ValueError(f"Porosity must be in [0, 1), got {porosity}")

    if porosity > 0.15:
        import warnings
        warnings.warn(
            f"Porosity {porosity:.2f} > 0.15; Maxwell-Eucken correction "
            "may not be accurate."
        )

    return kappa * (2 + porosity) / (2 - 2 * porosity)
