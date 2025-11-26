"""
XRD Multiline Plot Module

This module provides functionality for plotting XRD (X-ray Diffraction) data
with multiple patterns, particularly for thermal XRD analysis of LLZTO materials.

Features:
- Load thermal XRD data from multiple temperature files
- Normalize and smooth XRD patterns
- Plot reference patterns (ICSD database entries)
- Create temperature-dependent XRD plots with vertical offsets
- Support for Rietveld refinement data visualization
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path
import os
from typing import Dict, List, Tuple, Optional, Union
import sys

# Add src to path to allow imports if running as script
current_dir = Path(__file__).resolve().parent
src_root = current_dir.parent.parent
if str(src_root) not in sys.path:
    sys.path.append(str(src_root))



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
        self.thermal_data = {}
        self.datasets = {}
        
        # Default color cycle similar to MATLAB
        self.colors = plt.cm.tab10.colors
        
    def set_thermal_data(self, thermal_data: Dict[int, Dict[str, np.ndarray]]) -> None:
        """
        Set the thermal XRD data series.
        
        Args:
            thermal_data (Dict): Dictionary mapping temperature to data dictionary
        """
        self.thermal_data = thermal_data

    def add_dataset(self, name: str, data: np.ndarray) -> None:
        """
        Add a generic XRD dataset (sample, reference, etc.).
        
        Args:
            name (str): Unique identifier for the dataset
            data (np.ndarray): The XRD data array
        """
        self.datasets[name] = data
    
    def plot_pattern_comparison(self, sample_name: str, reference_name: Optional[str] = None,
                              figsize: Tuple[int, int] = (10, 6),
                              xlim: Tuple[float, float] = (20, 70),
                              ylim: Tuple[float, float] = (0, 2.0),
                              impurity_markers: Optional[List[Tuple[float, float]]] = None,
                              title: str = "") -> plt.Figure:
        """
        Plot a sample pattern, optionally with a reference stick pattern and impurity markers.
        
        Args:
            sample_name (str): Name of the sample dataset to plot
            reference_name (str): Name of the reference dataset to plot (optional)
            figsize (tuple): Figure size
            xlim (tuple): X-axis limits
            ylim (tuple): Y-axis limits
            impurity_markers (list): List of (x, y_scale) tuples for markers
            title (str): Plot title
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot Sample
        if sample_name in self.datasets:
            data = self.datasets[sample_name]
            # Determine if we have refinement data (obs, calc, bkg, diff) or just simple data
            has_refinement = data.shape[1] >= 5
            
            spacer = 0.8 if has_refinement else 0.6
            
            # Plot main intensity (Observed)
            ax.plot(data[:, 0], data[:, 1] + spacer, 
                   'x' if has_refinement else '-', 
                   color=self.colors[0], label=f'{sample_name} Obs', markersize=4)
            
            if has_refinement:
                ax.plot(data[:, 0], data[:, 2] + spacer, color=self.colors[1], label='Refinement')
                ax.plot(data[:, 0], data[:, 3] + spacer, color=self.colors[2], label='Background')
                ax.plot(data[:, 0], data[:, 4] + spacer - 0.1, color=self.colors[3], label='Residual')
        else:
            ax.text(0.5, 0.5, f"Sample '{sample_name}' not found", ha='center', transform=ax.transAxes)
        
        # Plot Reference (Stick Pattern)
        if reference_name and reference_name in self.datasets:
            ref_data = self.datasets[reference_name]
            y_base = 0.1
            n_peaks = min(52, len(ref_data))
            
            for i in range(n_peaks):
                ax.plot([ref_data[i, 0], ref_data[i, 0]], 
                       [y_base, ref_data[i, 1] / 2 + y_base], 'k-', linewidth=1)
            
            ax.axhline(y=y_base, color='k', linewidth=1)
            ax.text(xlim[0] + (xlim[1]-xlim[0])*0.4, 0.3 + y_base, reference_name, fontsize=8, color='black')

        # Plot Impurity Markers
        if impurity_markers:
            for x, y_scale in impurity_markers:
                # Adjust y position based on scale and fixed offset
                ax.plot(x + 0.2, y_scale * 0.1 + 0.9, 'v', color='black', label='Impurity')
            
            # Deduplicate legend labels
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())
        else:
            ax.legend()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('2θ (degree)')
        ax.set_ylabel('Intensity (a.u.)')
        ax.set_yticks([])
        ax.grid(True, alpha=0.3)
        if title:
            ax.set_title(title)
        
        return fig
    
    def plot_thermal_evolution(self, figsize: Tuple[int, int] = (10, 8),
                             reference_name: Optional[str] = None,
                             peak_markers: Optional[List[float]] = None) -> plt.Figure:
        """
        Plot thermal evolution of XRD patterns with vertical offset.
        
        Args:
            figsize (tuple): Figure size
            reference_name (str): Name of reference dataset to plot at bottom
            peak_markers (list): List of 2theta positions to mark with stars
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        if not self.thermal_data:
            print("No thermal XRD data loaded.")
            return fig
        
        spacer = 1.1
        temperatures = sorted(self.thermal_data.keys())
        
        # Plot thermal evolution patterns
        for i, temp in enumerate(temperatures):
            data = self.thermal_data[temp]
            color_idx = i % len(self.colors)
            
            ax.plot(data['two_theta'], data['intensity'] + spacer * i, 
                   color=self.colors[color_idx], 
                   label=f'{temp} K' if temp > 200 else f'{temp + 273} K') # Simple heuristic for K vs C
            
            # Add temperature label
            label_txt = f'{temp} K' if temp > 200 else f'{temp + 273} K'
            ax.text(18, 1.1 * (i + 1) - 1.3, label_txt, 
                   color=self.colors[color_idx], fontsize=8)
        
        # Mark specific peaks if requested
        if peak_markers and len(temperatures) >= 5:
            # Use the 5th temperature dataset for finding peak heights (arbitrary choice from original code)
            ref_temp = temperatures[min(4, len(temperatures)-1)]
            temp_data = self.thermal_data[ref_temp]
            two_theta = temp_data['two_theta']
            intensity = temp_data['intensity']
            
            for pos in peak_markers:
                idx = np.argmin(np.abs(two_theta - pos))
                # Plot marker at the height of the 5th plot (index 4)
                plot_idx = min(4, len(temperatures)-1)
                peak_intensity = intensity[idx] + spacer * plot_idx + 0.2
                ax.plot(two_theta[idx], peak_intensity, '*', color='black', markersize=8)
        
        # Plot Reference at bottom
        if reference_name and reference_name in self.datasets:
            ref_data = self.datasets[reference_name]
            y_base = -1.1
            scale_factor = 1.0
            n_peaks = min(52, len(ref_data))
            
            for i in range(n_peaks):
                ax.plot([ref_data[i, 0], ref_data[i, 0]], 
                       [y_base, ref_data[i, 1] * scale_factor + y_base], 
                       'k-', linewidth=1)
            
            ax.axhline(y=y_base, color='k', linewidth=1)
            ax.text(40, 0.6 + y_base, reference_name, fontsize=8, color='black')
        
        ax.set_xlim([16, 60])
        ax.set_ylim([-1.5, 5.5])
        ax.set_xlabel('2θ (degree)')
        ax.set_ylabel('Intensity (a.u.)')
        ax.set_yticks([])
        ax.grid(True, alpha=0.3)
        
        return fig


# Example usage
def main():
    """Example usage of the XRDMultilinePlotter class."""
    # Import here to avoid dependency in the module level if desired, 
    # but since we are in the same package structure, top-level import is also fine.
    # However, to demonstrate separation, we use the reader explicitly here.
    from MatSciKit_COSMOTIM.io.xrd_reader import XRDDataReader
    
    reader = XRDDataReader("XRD_plot_data")
    plotter = XRDMultilinePlotter()
    
    # 1. Load Reference Data (ICSD-894)
    icsd_data = reader.load_xrd_dataset("LLZTO_CollCode894.txt", 
                                      columns=[7, 8], normalize=False, scale_factor=0.01)
    if icsd_data is not None:
        plotter.add_dataset("ICSD-894", icsd_data)
    
    # 2. Load Pristine Sample Data
    pristine_data = reader.load_xrd_dataset("LLZTO_ploy1_powder_bed_sintered_GSAS.csv", 
                                          columns=[1, 2, 3, 4, 5])
    if pristine_data is not None:
        plotter.add_dataset("Pristine LLZTO", pristine_data)
    
    # 3. Load Aged Sample Data
    aged_data = reader.load_xrd_dataset("LLZTO_UCR_poly1_Nonsinter_GSAS.csv", 
                                      columns=[1, 2, 3, 4, 5])
    if aged_data is not None:
        plotter.add_dataset("Aged LLZTO", aged_data)
    
    # 4. Load Thermal Data Series
    temps = list(range(50, 251, 50)) # 50, 100, 150, 200, 250
    thermal_data = reader.load_thermal_series("{temp}C.txt", temps)
    plotter.set_thermal_data(thermal_data)
    
    # Create Plots
    
    # Figure 1: Pristine with Reference
    fig1 = plotter.plot_pattern_comparison("Pristine LLZTO", "ICSD-894", 
                                         title="Pristine LLZTO Analysis")
    
    # Figure 2: Aged with Reference and Impurities
    impurities = [(21.348, 0.853), (31.794, 1.0)]
    fig2 = plotter.plot_pattern_comparison("Aged LLZTO", "ICSD-894", 
                                         title="Aged LLZTO Analysis",
                                         xlim=(16, 70), ylim=(0, 1.7),
                                         impurity_markers=impurities)
    
    # Figure 3: Thermal Evolution
    fig3 = plotter.plot_thermal_evolution(reference_name="ICSD-894", 
                                        peak_markers=[40.28, 46.86])
    
    plt.show()


if __name__ == "__main__":
    main()

