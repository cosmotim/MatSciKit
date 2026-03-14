"""
DSC (Differential Scanning Calorimetry) data reader.

Reads heat capacity data from DSC CSV export files (e.g., Netzsch DSC 214).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def read(filepath: str, skip_rows: int = 34) -> np.ndarray:
    """
    Read DSC heat capacity data from a CSV export file.

    The file is expected to have header rows followed by data columns:
    [Temperature (°C), Time (min), Cp (J/(g·K))].
    Temperature is automatically converted from °C to K.

    Parameters
    ----------
    filepath : str
        Full path to the DSC CSV file.
    skip_rows : int, optional
        Number of header rows to skip. Default is 34 (Netzsch format).

    Returns
    -------
    data : np.ndarray
        Array with columns [Temperature (K), Cp (J/(g·K))].

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    # Read CSV data
    # DSC files may contain non-UTF-8 characters (e.g., µ in µV)
    # Read with latin-1 encoding, then parse with numpy
    with open(filepath, encoding="latin-1") as f:
        lines = f.readlines()[skip_rows:]

    # Parse numeric data from remaining lines
    raw = np.genfromtxt(lines, delimiter=",")

    if raw.ndim == 1:
        raw = raw.reshape(1, -1)

    # Extract Temperature (°C) and Cp columns, convert T to Kelvin
    temp_k = raw[:, 0] + 273.15
    cp = raw[:, 2]  # Third column is Cp

    return np.column_stack([temp_k, cp])
