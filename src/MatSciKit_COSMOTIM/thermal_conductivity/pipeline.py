"""
Pipeline 3: Thermal conductivity analysis.

High-level functions that chain IO reading, LFA→κ conversion,
porosity correction, and PPMS TTO comparison.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from MatSciKit_COSMOTIM.thermal_conductivity import lfa_dsc, porosity_correction


def tto_from_ppms(
    filepath: str,
    drop_temps: list[float] | None = None,
    porosity: float = 0.0,
    skip_rows: int = 27,
) -> dict[str, Any]:
    """Read PPMS TTO data and optionally apply porosity correction.

    Parameters
    ----------
    filepath : str
        Path to the PPMS TTO .dat file.
    drop_temps : list of float, optional
        Temperatures to exclude.
    porosity : float, optional
        Sample porosity for solid-state correction (default 0, no correction).
    skip_rows : int, optional
        Header rows to skip (default 27).

    Returns
    -------
    result : dict
        Dictionary with keys:

        - ``"temperature"`` : np.ndarray, temperature in K
        - ``"kappa"`` : np.ndarray, thermal conductivity in W/(m·K)
        - ``"kappa_error"`` : np.ndarray, conductivity error
        - ``"kappa_solid"`` : np.ndarray or None, porosity-corrected κ
        - ``"raw"`` : np.ndarray, raw [T, κ, κ_err] array
        - ``"porosity"`` : float
    """
    from MatSciKit_COSMOTIM.io import ppms_tto

    data = ppms_tto.read(filepath, drop_temps=drop_temps, skip_rows=skip_rows)
    T = data[:, 0]
    kappa = data[:, 1]
    kappa_err = data[:, 2]

    kappa_solid = None
    if porosity > 0:
        kappa_solid = porosity_correction.correct(kappa, porosity)

    return {
        "temperature": T,
        "kappa": kappa,
        "kappa_error": kappa_err,
        "kappa_solid": kappa_solid,
        "raw": data,
        "porosity": porosity,
    }


def kappa_from_lfa(
    cp_file: str,
    lfa_file: str,
    density: float,
    density_error: float,
    cp_format: str = "dsc",
    lfa_format: str = "csv",
    lfa_sheet: str | int = 0,
    porosity: float = 0.0,
    t_window: float = 5.0,
    cp_skip_rows: int | None = None,
) -> dict[str, Any]:
    """Convert LFA diffusivity to thermal conductivity using Cp data.

    Loads heat capacity (from DSC or PPMS HC) and LFA diffusivity data,
    computes κ = Cp × α × ρ with error propagation.

    Parameters
    ----------
    cp_file : str
        Path to heat capacity data file.
    lfa_file : str
        Path to LFA diffusivity data file.
    density : float
        Sample density in g/cm³.
    density_error : float
        Density error in g/cm³.
    cp_format : str, optional
        Format of Cp file: ``"dsc"`` or ``"ppms_hc"`` (default ``"dsc"``).
    lfa_format : str, optional
        Format of LFA file: ``"csv"`` or ``"excel"`` (default ``"csv"``).
    lfa_sheet : str or int, optional
        Sheet name/index for Excel LFA files (default 0).
    porosity : float, optional
        Sample porosity for correction (default 0).
    t_window : float, optional
        Temperature window (K) for Cp averaging (default 5.0).
    cp_skip_rows : int or None, optional
        Override skip_rows for Cp reader.

    Returns
    -------
    result : dict
        Dictionary with keys:

        - ``"temperature"`` : np.ndarray, LFA temperatures (K)
        - ``"kappa"`` : np.ndarray, thermal conductivity
        - ``"kappa_error"`` : np.ndarray, conductivity error
        - ``"kappa_solid"`` : np.ndarray or None, porosity-corrected
        - ``"cp_averaged"`` : np.ndarray, [T, Cp_avg, Cp_avg_err]
        - ``"diffusivity"`` : np.ndarray, raw LFA [T, α, α_err]
        - ``"density"`` : float
    """
    # Load heat capacity data
    if cp_format == "dsc":
        from MatSciKit_COSMOTIM.io import dsc

        kwargs = {}
        if cp_skip_rows is not None:
            kwargs["skip_rows"] = cp_skip_rows
        cp_data = dsc.read(cp_file, **kwargs)
        cp_T = cp_data[:, 0]
        cp_vals = cp_data[:, 1]
        # DSC doesn't provide error; estimate as 2%
        cp_err = cp_vals * 0.02
    elif cp_format == "ppms_hc":
        from MatSciKit_COSMOTIM.io import ppms_hc

        kwargs = {}
        if cp_skip_rows is not None:
            kwargs["skip_rows"] = cp_skip_rows
        cp_data = ppms_hc.read(cp_file, **kwargs)
        cp_T = cp_data[:, 0]
        cp_vals = cp_data[:, 1]
        cp_err = cp_data[:, 2]
    else:
        raise ValueError(f"Unknown cp_format: {cp_format!r}. Use 'dsc' or 'ppms_hc'.")

    # Load LFA diffusivity data
    if lfa_format == "csv":
        from MatSciKit_COSMOTIM.io import lfa

        lfa_data = lfa.read(lfa_file)
    elif lfa_format == "excel":
        from MatSciKit_COSMOTIM.io import lfa_excel

        lfa_data = lfa_excel.read(lfa_file, sheet_name=lfa_sheet)
    else:
        raise ValueError(f"Unknown lfa_format: {lfa_format!r}. Use 'csv' or 'excel'.")

    diff_T = lfa_data[:, 0]
    diffusivity = lfa_data[:, 1]
    diff_err = lfa_data[:, 2]

    # Compute κ = Cp × α × ρ
    tc_data, cp_averaged = lfa_dsc.calculate(
        cp_T,
        cp_vals,
        cp_err,
        diff_T,
        diffusivity,
        diff_err,
        density,
        density_error,
        t_window=t_window,
    )

    kappa = tc_data[:, 1]
    kappa_err = tc_data[:, 2]

    kappa_solid = None
    if porosity > 0:
        kappa_solid = porosity_correction.correct(kappa, porosity)

    return {
        "temperature": tc_data[:, 0],
        "kappa": kappa,
        "kappa_error": kappa_err,
        "kappa_solid": kappa_solid,
        "cp_averaged": cp_averaged,
        "diffusivity": lfa_data,
        "density": density,
    }


def compare_tto_lfa(
    tto_result: dict[str, Any],
    lfa_result: dict[str, Any],
) -> dict[str, Any]:
    """Compare PPMS TTO direct measurement with LFA-derived conductivity.

    Parameters
    ----------
    tto_result : dict
        Output from :func:`tto_from_ppms`.
    lfa_result : dict
        Output from :func:`kappa_from_lfa`.

    Returns
    -------
    comparison : dict
        Dictionary with keys:

        - ``"tto_temperature"`` : np.ndarray
        - ``"tto_kappa"`` : np.ndarray
        - ``"lfa_temperature"`` : np.ndarray
        - ``"lfa_kappa"`` : np.ndarray
        - ``"overlap_temps"`` : np.ndarray, temperatures where both have data
        - ``"tto_at_overlap"`` : np.ndarray, TTO κ interpolated at overlap temps
        - ``"lfa_at_overlap"`` : np.ndarray, LFA κ interpolated at overlap temps
        - ``"ratio"`` : np.ndarray, TTO/LFA ratio at overlap points
    """
    tto_T = tto_result["temperature"]
    tto_k = tto_result["kappa"]
    lfa_T = lfa_result["temperature"]
    lfa_k = lfa_result["kappa"]

    # Find overlapping temperature range
    t_min = max(np.min(tto_T), np.min(lfa_T))
    t_max = min(np.max(tto_T), np.max(lfa_T))

    if t_min >= t_max:
        # No overlap
        return {
            "tto_temperature": tto_T,
            "tto_kappa": tto_k,
            "lfa_temperature": lfa_T,
            "lfa_kappa": lfa_k,
            "overlap_temps": np.array([]),
            "tto_at_overlap": np.array([]),
            "lfa_at_overlap": np.array([]),
            "ratio": np.array([]),
        }

    # Use LFA temperatures that fall within TTO range as overlap points
    overlap_mask = (lfa_T >= t_min) & (lfa_T <= t_max)
    overlap_temps = lfa_T[overlap_mask]

    # Interpolate TTO at overlap temperatures
    tto_interp = np.interp(overlap_temps, tto_T, tto_k)
    lfa_at_overlap = lfa_k[overlap_mask]

    ratio = tto_interp / lfa_at_overlap

    return {
        "tto_temperature": tto_T,
        "tto_kappa": tto_k,
        "lfa_temperature": lfa_T,
        "lfa_kappa": lfa_k,
        "overlap_temps": overlap_temps,
        "tto_at_overlap": tto_interp,
        "lfa_at_overlap": lfa_at_overlap,
        "ratio": ratio,
    }
