"""
Thermal conductivity from LFA diffusivity and DSC heat capacity.

Combines laser flash analysis (thermal diffusivity) with differential
scanning calorimetry (heat capacity) to compute thermal conductivity:

    κ = Cp × α × ρ

Translated from LFA_DSC_TC.m
"""

import numpy as np
from typing import Tuple


def calculate(cp_T: np.ndarray, cp: np.ndarray, cp_error: np.ndarray,
              diff_T: np.ndarray, diffusivity: np.ndarray,
              diff_error: np.ndarray,
              density: float, density_error: float,
              t_window: float = 5.0
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate thermal conductivity from LFA diffusivity and DSC Cp.

    For each LFA measurement temperature, averages nearby DSC Cp values
    (within ±t_window K), then computes κ = Cp × α × ρ with full
    error propagation.

    Parameters
    ----------
    cp_T : np.ndarray
        DSC/PPMS temperature values (K).
    cp : np.ndarray
        Heat capacity values (J/(g·K)).
    cp_error : np.ndarray
        Heat capacity errors (J/(g·K)).
    diff_T : np.ndarray
        LFA measurement temperatures (K).
    diffusivity : np.ndarray
        Thermal diffusivity values (mm²/s or appropriate units).
    diff_error : np.ndarray
        Diffusivity errors.
    density : float
        Sample density (g/cm³ or consistent with Cp and diffusivity units).
    density_error : float
        Density error.
    t_window : float, optional
        Temperature window (K) for averaging DSC data at each LFA point.
        Default is 5.0 K.

    Returns
    -------
    tc_data : np.ndarray
        Array with columns [Temperature, κ, κ_error].
    cp_averaged : np.ndarray
        Array with columns [Temperature, Cp_avg, Cp_avg_error] at LFA temperatures.

    Notes
    -----
    Error propagation: (δκ/κ)² = (δCp/Cp)² + (δα/α)² + (δρ/ρ)²
    """
    n_lfa = len(diff_T)

    cp_avg = np.zeros(n_lfa)
    cp_avg_err = np.zeros(n_lfa)

    for i in range(n_lfa):
        # Find DSC points within ±t_window of this LFA temperature
        idx = np.where(np.abs(cp_T - diff_T[i]) <= t_window)[0]

        if len(idx) == 0:
            # No nearby Cp data — interpolate instead
            cp_avg[i] = np.interp(diff_T[i], cp_T, cp)
            # Estimate error from nearest point
            nearest = np.argmin(np.abs(cp_T - diff_T[i]))
            cp_avg_err[i] = cp_error[nearest]
        else:
            cp_avg[i] = np.mean(cp[idx])
            cp_avg_err[i] = np.sqrt(np.sum(cp_error[idx]**2))

    # Calculate thermal conductivity: κ = Cp × α × ρ
    tc = cp_avg * diffusivity * density

    # Error propagation: (δκ/κ)² = (δCp/Cp)² + (δα/α)² + (δρ/ρ)²
    relative_error = np.sqrt(
        (cp_avg_err / cp_avg)**2 +
        (diff_error / diffusivity)**2 +
        (density_error / density)**2
    )
    tc_error = tc * relative_error

    tc_data = np.column_stack([diff_T, tc, tc_error])
    cp_averaged = np.column_stack([diff_T, cp_avg, cp_avg_err])

    return tc_data, cp_averaged
