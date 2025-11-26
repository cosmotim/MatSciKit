"""
XRD Multiline Plot Module

This module provides functionality for plotting XRD (X-ray Diffraction) data
with multiple patterns, particularly for thermal XRD analysis of LLZTO materials.

Features:
- Plot multiple XRD patterns with vertical offsets
- Support for reference patterns (stick plots)
- Peak and impurity markers
- Text tags or legend labeling
- Refinement data visualization (observed, calculated, background, residual)
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional, Union, Any





class XRDMultilinePlotter:
    """
    A class for plotting multiple XRD patterns with various visualization options.

    This class handles visualizing XRD data from thermal measurements,
    reference patterns, and aged samples.
    """

    def __init__(self):
        """
        Initialize the XRD plotter.
        """
        # Default color cycle similar to MATLAB
        self.colors = plt.cm.tab10.colors

    def plot_patterns(self, data_list: Union[np.ndarray, List[np.ndarray]],
                     labels: Optional[List[str]] = None,
                     y_shift: float = 1.0,
                     markers: Optional[List[float]] = None
                     ) -> plt.Figure:
        """
        General method to plot multiple XRD patterns.

        Args:
            data_list (Union[np.ndarray, List[np.ndarray]]): Single numpy array or list of numpy arrays from XRDDataReader.
            labels (List[str]): List of labels for each pattern. If None, generic labels are used.
            y_shift (float): Absolute vertical offset between patterns (default: 1.0).
            markers (list): List of 2theta positions to mark with triangular markers.
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        if not isinstance(data_list, list):
            data_list = [data_list]

        # Set default labels if not provided
        if labels is None:
            labels = [f"Pattern {i+1}" for i in range(len(data_list))]
        elif len(labels) != len(data_list):
            print(f"Warning: Number of labels ({len(labels)}) doesn't match number of datasets ({len(data_list)}). Using default labels.")
            labels = [f"Pattern {i+1}" for i in range(len(data_list))]

        # Use y_shift as absolute shift amount
        shift_amount = y_shift

        # Plot patterns
        for i, data_item in enumerate(data_list):
            if data_item is None or data_item.shape[1] < 2:
                print(f"Warning: Could not extract data from item {i}.")
                continue

            color_idx = i % len(self.colors)
            
            # Extract two_theta and intensity from numpy array
            two_theta = data_item[:, 0]
            intensity = data_item[:, 1]
            
            # Apply vertical offset for multiple patterns
            y_offset = i * shift_amount

            # Plot Main Intensity
            ax.plot(two_theta, intensity + y_offset,
                   '-',
                   color=self.colors[color_idx],
                   label=labels[i])

        # Plot Peak Markers (Stars)
        if markers and len(data_list) >= 1:
            # Use the first dataset for height reference
            data_item = data_list[0]
            
            if data_item is not None and data_item.shape[1] >= 2:
                two_theta = data_item[:, 0]
                intensity = data_item[:, 1]
                
                for pos in markers:
                    idx = np.argmin(np.abs(two_theta - pos))
                    peak_intensity = intensity[idx] * 1.5  # Markers on first pattern only
                    ax.plot(two_theta[idx], peak_intensity, 'v', color='black', markersize=8)

        # Finalize Plot
        ax.set_xlabel('2θ (degree)')
        ax.set_ylabel('Intensity (a.u.)')
        ax.set_yticks([])
        ax.grid(True, alpha=0.3)

        # Legend handling
        handles, labels_legend = ax.get_legend_handles_labels()
        by_label = dict(zip(labels_legend, handles))
        ax.legend(by_label.values(), by_label.keys())

        return fig


def main():
    """Example usage of the XRDMultilinePlotter class."""
    import sys
    import os
    from pathlib import Path
    
    # Add src to path to allow imports
    current_dir = Path(__file__).resolve().parent
    src_root = current_dir.parent.parent
    if str(src_root) not in sys.path:
        sys.path.append(str(src_root))
    
    from MatSciKit_COSMOTIM.io.xrd_reader import XRDDataReader
    from MatSciKit_COSMOTIM.visualization.export_journal_style import export_journal_figure

    # Try find this data file and other files in that folder ExampleData/XRD_data/Li_08.txt
    project_root = src_root.parent  # Go up from src/ to project root
    data_dir = project_root / "ExampleData" / "XRD_data"
    if not data_dir.exists():
        print(f"Error: Data directory {data_dir} not found.")
        return

    

    # Initialize reader and plotter
    reader = XRDDataReader(str(data_dir))
    plotter = XRDMultilinePlotter()
    
    # Create a name list based on filenames in the ExampleData/XRD_data directory which can be used for plotting later
    data_filenames = [fname for fname in os.listdir(str(data_dir)) if fname.endswith(".txt")]
    print("Available data files:")
    for fname in data_filenames:
        print(f" - {fname}")

    # Create output directory for figures
    output_dir = project_root / "ExampleFigures"
    output_dir.mkdir(exist_ok=True)

    # Example 1: Load and plot single or multiple XRD patterns using available files
    print("\nExample 1: Single or Multiple XRD patterns")
    data_list = []
    name_list = []
    # Load first 3 files from available data files (or all if less than 3)
    for fname in data_filenames[:3]:
        pattern = reader.load_xrd_dataset(fname, columns=[0, 1], normalize=False)
        if pattern is not None:
            data_list.append(pattern)
            # Use filename without extension as label
            name_list.append(fname.replace('.txt', ''))
    
    if data_list:
        fig1 = plotter.plot_patterns(data_list, labels=name_list)
        export_journal_figure(fig1, output_dir / "example1_patterns", format='png', dpi=600)
        # plt.savefig(output_dir / "example1_patterns.png", dpi=300)
        print(f"Saved {output_dir / 'example1_patterns.tiff'} and .png")
    
    # Example 2: Plot with peak and impurity markers (using first available file)
    print("\nExample 2: XRD pattern with markers")
    if data_filenames:
        sample_data = reader.load_xrd_dataset(data_filenames[0], columns=[0, 1], normalize=True)
        
        if sample_data is not None:
            peak_positions = [25.5, 32.1, 40.8, 46.2]  # Example 2theta positions
            fig2 = plotter.plot_patterns(
                sample_data, 
                markers = peak_positions
            )
            export_journal_figure(fig2, output_dir / "example2_with_markers", format='png', dpi=600)
            # plt.savefig(output_dir / "example2_with_markers.png", dpi=300)
            print(f"Saved {output_dir / 'example2_with_markers.tiff'} and .png")
    


if __name__ == "__main__":
    main()
