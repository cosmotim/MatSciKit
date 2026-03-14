"""
PPMS HC (Heat Capacity) data reader.

Reads heat capacity data from Quantum Design PPMS HC .dat files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def read(
    filepath: str, drop_temps: list[float] | np.ndarray | None = None, skip_rows: int = 15
) -> np.ndarray:
    """
    Read and process heat capacity data from a PPMS HC .dat file.

    Parameters
    ----------
    filepath : str
        Full path to the HC .dat file.
    drop_temps : list of float, optional
        Temperatures at which to drop data points. Default is None.
    skip_rows : int, optional
        Number of header rows to skip. Default is 15.

    Returns
    -------
    data : np.ndarray
        Array with columns [Temperature (K), Cp (J/(g·K)), Cp_error].

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    if drop_temps is None:
        drop_temps = []
    elif isinstance(drop_temps, np.ndarray):
        drop_temps = drop_temps.tolist()

    # Read data, skipping header rows
    # PPMS HC files can be comma or tab separated
    try:
        data = np.genfromtxt(filepath, delimiter=",", skip_header=skip_rows)
    except ValueError:
        data = np.genfromtxt(filepath, skip_header=skip_rows)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Extract columns: Temperature (col 4), Cp (col 9), Cp_error (col 10) — 0-indexed
    # MATLAB: Data(:, 5), Data(:, 10), Data(:, 11) — 1-indexed
    if data.shape[1] < 11:
        raise ValueError(
            f"Data file has insufficient columns. Expected at least 11, got {data.shape[1]}"
        )

    result = data[:, [4, 9, 10]]

    # Drop specified temperature points
    rows_to_drop = []
    for temp in drop_temps:
        idx = np.argmin(np.abs(result[:, 0] - temp))
        rows_to_drop.append(idx)

    if rows_to_drop:
        result = np.delete(result, rows_to_drop, axis=0)

    return result
