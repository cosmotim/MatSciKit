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


    def load_xrd_dataset(self, filename: str, columns: List[int] = [0, 1], 
                    normalize: bool = False, scale_factor: float = 1.0
                    ) -> Optional[np.ndarray]:
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
            # Try reading with whitespace delimiter (common for XRD data)
            try:
                df = pd.read_csv(file_path, sep=r'\s+', engine='python')
            except:
                # Fall back to auto-detect separator
                try:
                    df = pd.read_csv(file_path, sep=None, engine='python')
                except:
                    # Last resort: try tab-separated
                    df = pd.read_csv(file_path, sep='\t')
                
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

