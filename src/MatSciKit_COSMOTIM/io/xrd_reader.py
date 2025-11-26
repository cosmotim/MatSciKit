import numpy as np
import pandas as pd
from scipy.signal import medfilt
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple

class XRDDataReader:
    """
    A class for reading and processing XRD data files.
    """
    
    def __init__(self, data_directory: str):
        """
        Initialize the XRD data reader.
        
        Args:
            data_directory (str): Path to directory containing XRD data files
        """
        self.data_directory = Path(data_directory)

    def load_thermal_series(self, file_pattern: str, temperatures: List[int], 
                           smooth: bool = True, smooth_window: int = 5) -> Dict[int, Dict[str, np.ndarray]]:
        """
        Load thermal XRD data series based on a filename pattern.
        
        Args:
            file_pattern (str): Pattern for filenames containing {temp}, e.g., "{temp}C.txt"
            temperatures (List[int]): List of temperatures to load
            smooth (bool): Whether to apply smoothing
            smooth_window (int): Window size for median filtering
            
        Returns:
            Dict: Dictionary mapping temperature to data dictionary {'two_theta': ..., 'intensity': ...}
        """
        thermal_data = {}
        
        for temp in temperatures:
            filename = file_pattern.format(temp=temp)
            file_path = self.data_directory / filename
            
            if file_path.exists():
                try:
                    # Try reading with tab separator first (common for these files)
                    try:
                        df = pd.read_csv(file_path, sep='\t', header=None)
                        if df.shape[1] < 2:
                             df = pd.read_csv(file_path, sep=None, engine='python', header=None)
                    except:
                        df = pd.read_csv(file_path, sep=None, engine='python', header=None)
                        
                    data = df.values
                    if data.shape[1] >= 2:
                        two_theta = data[:, 0]
                        intensity = data[:, 1]
                        
                        # Normalize intensity
                        if np.max(intensity) > 0:
                            intensity_norm = intensity / np.max(intensity)
                        else:
                            intensity_norm = intensity
                        
                        # Apply smoothing if requested
                        if smooth:
                            intensity_norm = medfilt(intensity_norm, smooth_window)
                        
                        thermal_data[temp] = {
                            'two_theta': two_theta,
                            'intensity': intensity_norm
                        }
                        
                except Exception as e:
                    print(f"Error loading {file_path}: {e}")
            else:
                print(f"Warning: File {filename} not found in {self.data_directory}")
                
        return thermal_data

    def load_xrd_dataset(self, filename: str, columns: List[int] = [0, 1], 
                    normalize: bool = True, scale_factor: float = 1.0, 
                    header: Optional[int] = None) -> Optional[np.ndarray]:
        """
        Load a generic XRD dataset (sample, reference, etc.).
        
        Args:
            filename (str): Filename in the data directory
            columns (List[int]): List of column indices to extract [x, y1, y2, ...]
            normalize (bool): Whether to normalize intensities by max of first y column
            scale_factor (float): Scaling factor if not normalizing
            header (int): Row number to use as header (None for no header)
            
        Returns:
            np.ndarray: Processed data array or None if loading failed
        """
        file_path = self.data_directory / filename
        if not file_path.exists():
            print(f"Error: File {filename} not found.")
            return None

        try:
            # Auto-detect separator or try common ones
            try:
                df = pd.read_csv(file_path, sep=None, engine='python', header=header)
            except:
                df = pd.read_csv(file_path, sep='\t', header=header)
                
            data = df.values
            
            # Check if we have enough columns
            if data.shape[1] <= max(columns):
                print(f"Error: File {filename} has {data.shape[1]} columns, but requested index {max(columns)}")
                return None

            selected_data = data[:, columns]
            
            # Normalize or scale
            # Assuming column 0 is x (2theta), and rest are y (intensities)
            if normalize:
                # Normalize all y columns by the max of the first y column (index 1 in selected_data)
                max_val = np.max(selected_data[:, 1])
                if max_val > 0:
                    selected_data[:, 1:] = selected_data[:, 1:] / max_val
            elif scale_factor != 1.0:
                selected_data[:, 1:] = selected_data[:, 1:] * scale_factor
                
            print(f"Loaded dataset from {filename}")
            return selected_data
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None
