"""
Thermal conductivity plotting with optional κ_min overlay.

Generates publication-ready thermal conductivity plots from
TTO, LFA, or combined data, with the Cahill minimum model
as a reference curve.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure


def plot_kappa(
    datasets: list[dict[str, Any]],
    kappa_min_params: dict[str, Any] | None = None,
    title: str = "",
    xlabel: str = "Temperature (K)",
    ylabel: str = "Thermal Conductivity (W m⁻¹ K⁻¹)",
    figsize: tuple[float, float] = (6, 4.5),
    loglog: bool = False,
    t_range: tuple[float, float] | None = None,
    show: bool = True,
    save: str | Path | None = None,
    dpi: int = 300,
) -> Figure:
    """Plot thermal conductivity data with optional κ_min curve.

    Parameters
    ----------
    datasets : list of dict
        Each dict defines a dataset to plot with keys:

        - ``"temperature"`` : np.ndarray (required)
        - ``"kappa"`` : np.ndarray (required)
        - ``"kappa_error"`` : np.ndarray (optional)
        - ``"label"`` : str (optional, for legend)
        - ``"marker"`` : str (optional, default ``"o"``)
        - ``"color"`` : str (optional)
        - ``"linestyle"`` : str (optional, default ``"none"`` for data)
        - ``"linewidth"`` : float (optional)
    kappa_min_params : dict, optional
        Parameters for computing the Cahill κ_min curve:

        - ``"n_density"`` : float (required) — atoms/m³
        - ``"theta_D"`` : float (required) — Debye temperature (K)
        - ``"v_s"`` : float (required) — sound velocity (m/s)
        - ``"label"`` : str (optional, default ``"κ_min (Cahill)"``)
        - ``"color"`` : str (optional, default ``"black"``)
        - ``"linestyle"`` : str (optional, default ``"--"``)
        - ``"t_range"`` : tuple (optional) — override T range for curve
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    figsize : tuple, optional
        Figure size (width, height) in inches.
    loglog : bool, optional
        Use log-log axes (default False).
    t_range : tuple of (float, float), optional
        Temperature range for the κ_min curve. If None, derived from data.
    show : bool, optional
        Call ``plt.show()`` (default True).
    save : str or Path, optional
        Save figure to this path.
    dpi : int, optional
        Resolution for saved figure (default 300).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure.

    Examples
    --------
    >>> plot_kappa(
    ...     datasets=[
    ...         {"temperature": tto_T, "kappa": tto_k, "kappa_error": tto_err,
    ...          "label": "PPMS TTO", "marker": "o"},
    ...         {"temperature": lfa_T, "kappa": lfa_k, "kappa_error": lfa_err,
    ...          "label": "LFA", "marker": "s"},
    ...     ],
    ...     kappa_min_params={
    ...         "n_density": 7.63e28,
    ...         "theta_D": 312.4,
    ...         "v_s": 2841,
    ...     },
    ... )
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Default colors cycle
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # Track T range for κ_min curve
    all_t_min = []
    all_t_max = []

    # Plot each dataset
    for i, ds in enumerate(datasets):
        T = ds["temperature"]
        kappa = ds["kappa"]
        kappa_err = ds.get("kappa_error")
        label = ds.get("label", f"Dataset {i + 1}")
        marker = ds.get("marker", "o")
        color = ds.get("color", colors[i % len(colors)])
        ls = ds.get("linestyle", "none")
        lw = ds.get("linewidth", 1.5)

        all_t_min.append(np.min(T))
        all_t_max.append(np.max(T))

        if kappa_err is not None:
            ax.errorbar(
                T,
                kappa,
                yerr=kappa_err,
                marker=marker,
                linestyle=ls,
                color=color,
                label=label,
                capsize=3,
                markersize=5,
                linewidth=lw,
            )
        else:
            ax.plot(
                T,
                kappa,
                marker=marker,
                linestyle=ls,
                color=color,
                label=label,
                markersize=5,
                linewidth=lw,
            )

    # Plot κ_min curve
    if kappa_min_params is not None:
        from MatSciKit_COSMOTIM.thermal_conductivity import cahill

        n_density = kappa_min_params["n_density"]
        theta_D = kappa_min_params["theta_D"]
        v_s = kappa_min_params["v_s"]
        km_label = kappa_min_params.get("label", "κ_min (Cahill)")
        km_color = kappa_min_params.get("color", "black")
        km_ls = kappa_min_params.get("linestyle", "--")
        km_lw = kappa_min_params.get("linewidth", 1.5)

        # Determine T range
        if "t_range" in kappa_min_params:
            t_lo, t_hi = kappa_min_params["t_range"]
        elif t_range is not None:
            t_lo, t_hi = t_range
        elif all_t_min:
            t_lo = max(1.0, min(all_t_min) * 0.8)
            t_hi = max(all_t_max) * 1.1
        else:
            t_lo, t_hi = 5.0, 300.0

        T_curve = np.linspace(t_lo, t_hi, 300)
        kappa_min = np.array([cahill.minimum_tc(t, n_density, theta_D, v_s) for t in T_curve])

        ax.plot(
            T_curve,
            kappa_min,
            linestyle=km_ls,
            color=km_color,
            label=km_label,
            linewidth=km_lw,
        )

    # Formatting
    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.legend(frameon=True, fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()

    return fig


def plot_tto_vs_lfa(
    tto_result: dict[str, Any],
    lfa_result: dict[str, Any],
    kappa_min_params: dict[str, Any] | None = None,
    tto_label: str = "PPMS TTO",
    lfa_label: str = "LFA",
    **kwargs: Any,
) -> Figure:
    """Plot PPMS TTO data alongside LFA-derived κ, with optional κ_min.

    Convenience wrapper around :func:`plot_kappa` for comparing
    TTO and LFA results.

    Parameters
    ----------
    tto_result : dict
        Output from ``pipeline.tto_from_ppms()``.
    lfa_result : dict
        Output from ``pipeline.kappa_from_lfa()``.
    kappa_min_params : dict, optional
        Cahill κ_min parameters (see :func:`plot_kappa`).
    tto_label : str, optional
        Legend label for TTO data.
    lfa_label : str, optional
        Legend label for LFA data.
    **kwargs
        Additional keyword arguments passed to :func:`plot_kappa`.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    datasets = [
        {
            "temperature": tto_result["temperature"],
            "kappa": tto_result.get("kappa_solid", tto_result["kappa"]),
            "kappa_error": tto_result.get("kappa_error"),
            "label": tto_label,
            "marker": "o",
        },
        {
            "temperature": lfa_result["temperature"],
            "kappa": lfa_result.get("kappa_solid", lfa_result["kappa"]),
            "kappa_error": lfa_result.get("kappa_error"),
            "label": lfa_label,
            "marker": "s",
        },
    ]

    return plot_kappa(datasets, kappa_min_params=kappa_min_params, **kwargs)
