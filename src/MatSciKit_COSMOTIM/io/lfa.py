"""
LFA (Laser Flash Analysis) data reader.

Reads thermal diffusivity or conductivity data from LFA CSV files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def read(filepath: str) -> np.ndarray:
    """
    Read LFA data from a CSV file.

    Expected format: columns [Temperature (K), Value, Error].
    Rows containing NaN values are automatically dropped.

    Parameters
    ----------
    filepath : str
        Full path to the LFA CSV file.

    Returns
    -------
    data : np.ndarray
        Array with columns [Temperature (K), Value, Error], NaN rows removed.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    data = np.genfromtxt(filepath, delimiter=",")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Drop rows with any NaN
    mask = ~np.any(np.isnan(data), axis=1)
    data = data[mask]

    return data
