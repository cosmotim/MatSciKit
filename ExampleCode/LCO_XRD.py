"""
Example script for plotting XRD patterns using MatSciKit_COSMOTIM.

This script demonstrates how to:
1. Load multiple XRD data files
2. Plot them with custom labels
3. Add peak markers
4. Export in journal-style format
"""

import os
from pathlib import Path

from MatSciKit_COSMOTIM.io.xrd_reader import XRDDataReader
from MatSciKit_COSMOTIM.core.XRD_multiline_plot import XRDMultilinePlotter
from MatSciKit_COSMOTIM.visualization.export_journal_style import export_journal_figure


def main():
    """Main function to plot XRD patterns."""
    
    # Set up paths
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "ExampleData" / "XRD_data"
    if not data_dir.exists():
        print(f"Error: Data directory {data_dir} not found.")
        return

    # Initialize reader and plotter
    reader = XRDDataReader(str(data_dir))
    plotter = XRDMultilinePlotter()
    
    # Get available data files
    data_filenames = [fname for fname in os.listdir(str(data_dir)) if fname.endswith(".txt")]
    print("Available data files:")
    for fname in data_filenames:
        print(f" - {fname}")

    # Create output directory for figures
    output_dir = project_root / "ExampleFigures"

    # Example 1: Load and plot multiple XRD patterns
    print("\nExample 1: Multiple XRD patterns")
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
        print(f"Saved {output_dir / 'example1_patterns.png'} and .tiff")
    
    # Example 2: Plot with peak markers
    print("\nExample 2: XRD pattern with markers")
    if data_filenames:
        sample_data = reader.load_xrd_dataset(data_filenames[0], columns=[0, 1], normalize=True)
        
        if sample_data is not None:
            peak_positions = [25.5, 32.1, 40.8, 46.2]  # Example 2theta positions
            fig2 = plotter.plot_patterns(
                sample_data, 
                markers=peak_positions
            )
            export_journal_figure(fig2, output_dir / "example2_with_markers", format='png', dpi=600)
            print(f"Saved {output_dir / 'example2_with_markers.png'} and .tiff")


if __name__ == "__main__":
    main()
