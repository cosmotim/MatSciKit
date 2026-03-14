"""
PPMS TTO (Thermal Transport Option) data reader.

Reads thermal conductivity data from Quantum Design PPMS TTO .dat files.
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Optional, Union, List


def read(filepath: str,
         drop_temps: Optional[Union[List[float], np.ndarray]] = None,
         skip_rows: int = 27) -> np.ndarray:
    """
    Read and process thermal conductivity data from a PPMS TTO .dat file.

    Parameters
    ----------
    filepath : str
        Full path to the TTO .dat file.
    drop_temps : list of float, optional
        Temperatures at which to drop data points. Default is None.
    skip_rows : int, optional
        Number of header rows to skip. Default is 27.

    Returns
    -------
    data : np.ndarray
        Array with columns [Temperature (K), Conductivity (W/m/K), Error].

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If drop_temps contains non-numeric values or data has insufficient columns.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    if drop_temps is None:
        drop_temps = []
    elif isinstance(drop_temps, np.ndarray):
        drop_temps = drop_temps.tolist()

    # Read data, skipping header rows
    data = np.genfromtxt(filepath, delimiter=',', skip_header=skip_rows)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 8:
        raise ValueError(
            f"Data file has insufficient columns. Expected at least 8, got {data.shape[1]}"
        )

    # Check for non-positive temperatures (column index 5, 0-based)
    if np.any(data[:, 5] <= 0):
        import warnings
        warnings.warn(
            "This function uses K as a unit of temperature. "
            "Some temperature data is below or equal to 0."
        )

    # Drop specified temperature points
    rows_to_drop = []
    for temp in drop_temps:
        idx = np.argmin(np.abs(data[:, 5] - temp))
        rows_to_drop.append(idx)

    if rows_to_drop:
        data = np.delete(data, rows_to_drop, axis=0)

    # Extract columns: Temperature, Conductivity, Error (0-indexed: 5, 6, 7)
    return data[:, 5:8]
