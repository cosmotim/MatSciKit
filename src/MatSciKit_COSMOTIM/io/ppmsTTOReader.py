"""
PPMS TTO (Thermal Transport Option) Data Reader

This module provides functionality to read and process thermal conductivity data
from Quantum Design PPMS (Physical Property Measurement System) TTO.dat files.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Union, List
import warnings


def ttodataplot(filename: str, datadrop: Optional[Union[List[float], np.ndarray]] = None,
                base_directory: str = "/Users/ywang/Documents/QdDynacool/TTO/Yitian Wang") -> np.ndarray:
    """
    Read and process thermal conductivity data from a PPMS TTO.dat file.
    
    This function reads thermal transport data from Quantum Design PPMS files,
    with options to drop specific temperature data points.
    
    Args:
        filename (str): TTO.dat filename (without .dat extension)
        datadrop (list or ndarray, optional): Array of temperatures at which to drop data points.
            Default is None (no data dropped).
        base_directory (str, optional): Base directory path where TTO data files are stored.
            Default is "/Users/ywang/Documents/QdDynacool/TTO/Yitian Wang"
    
    Returns:
        np.ndarray: Processed data array with columns [Temperature, Conductivity, Error]
        
    Raises:
        FileNotFoundError: If the specified file doesn't exist
        ValueError: If datadrop is not numeric
        
    Example:
        >>> data = ttodataplot("sample_data", datadrop=[50.0, 100.0])
        >>> temperature = data[:, 0]
        >>> conductivity = data[:, 1]
        >>> error = data[:, 2]
    """
    # Validate inputs
    if datadrop is None:
        datadrop = []
    else:
        # Convert to list if it's a numpy array
        if isinstance(datadrop, np.ndarray):
            datadrop = datadrop.tolist()
        # Validate that datadrop is numeric
        if not all(isinstance(x, (int, float)) for x in datadrop):
            raise ValueError('datadrop must be a numeric array.')
    
    # Construct file path
    filepath = Path(base_directory) / f"{filename}.dat"
    
    # Check if file exists
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Read data from file (skip first 27 rows, start from row 28)
    # Using pandas to read the data
    try:
        data = pd.read_csv(filepath, sep=r'\s+', skiprows=27, header=None)
        Data = data.values
    except Exception as e:
        raise IOError(f"Error reading file {filepath}: {e}")
    
    # Check if we have enough columns (need at least columns 0-7 for indices 6-7 in MATLAB, 5-7 in Python)
    if Data.shape[1] < 8:
        raise ValueError(f"Data file has insufficient columns. Expected at least 8, got {Data.shape[1]}")
    
    # Check if all temperature data is above 0 (column 6 in MATLAB = column 5 in Python, 0-indexed)
    # Note: MATLAB uses 1-based indexing, Python uses 0-based
    if np.any(Data[:, 5] <= 0):
        warnings.warn('This function uses K as a unit of temperature. '
                     'Some temperature data is below or equal to 0.')
    
    # Drop specified data points
    rows_to_drop = []
    for temp in datadrop:
        # Find the row index with temperature closest to the specified value
        dx = np.argmin(np.abs(Data[:, 5] - temp))
        rows_to_drop.append(dx)
        print(f"Data of T={temp:.6f} is dropped")
    
    # Remove rows (in reverse order to avoid index shifting)
    Data = np.delete(Data, rows_to_drop, axis=0)
    
    # Output processed data (columns 6:8 in MATLAB = columns 5:8 in Python, 0-indexed)
    Data_out = Data[:, 5:8]
    
    return Data_out

