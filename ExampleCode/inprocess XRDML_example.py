"""
Example script demonstrating how to use the XRDML reader.

This script shows how to read Panalytical .xrdml files and plot the data.
"""

import numpy as np
import matplotlib.pyplot as plt
from MatSciKit_COSMOTIM.io import XRDMLReader, read_xrdml

# Example 1: Using the convenience function
def example_simple_read():
    """Simple example using the convenience function."""
    print("Example 1: Simple read using convenience function")
    print("-" * 50)
    
    # Read a .xrdml file (replace with your actual file path)
    filepath = "path/to/your/data.xrdml"
    
    try:
        two_theta, intensity, metadata = read_xrdml(filepath, normalize=True)
        
        print(f"Data points: {len(two_theta)}")
        print(f"2θ range: {two_theta[0]:.2f}° - {two_theta[-1]:.2f}°")
        print("\nMetadata:")
        for key, value in metadata.items():
            print(f"  {key}: {value}")
        
        # Plot the data
        plt.figure(figsize=(10, 6))
        plt.plot(two_theta, intensity, 'b-', linewidth=1)
        plt.xlabel('2θ (degrees)', fontsize=12)
        plt.ylabel('Normalized Intensity', fontsize=12)
        plt.title(f"XRD Pattern: {metadata.get('sample_id', 'Unknown')}", fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        print("Please update the filepath in the script.")


# Example 2: Using the XRDMLReader class for multiple scans
def example_multi_scan():
    """Example reading multiple scans from a file."""
    print("\nExample 2: Reading multiple scans")
    print("-" * 50)
    
    filepath = "path/to/your/data.xrdml"
    
    try:
        # Create reader instance
        reader = XRDMLReader(filepath)
        
        print(f"File: {reader.filepath.name}")
        print(f"Number of scans: {reader.get_num_scans()}")
        
        # Read all scans
        all_scans = reader.read_all_scans()
        
        # Plot all scans
        fig, axes = plt.subplots(reader.get_num_scans(), 1, 
                                figsize=(10, 4*reader.get_num_scans()))
        
        if reader.get_num_scans() == 1:
            axes = [axes]
        
        for i, (two_theta, intensity, metadata) in enumerate(all_scans):
            axes[i].plot(two_theta, intensity, 'b-', linewidth=1)
            axes[i].set_xlabel('2θ (degrees)', fontsize=11)
            axes[i].set_ylabel('Intensity (counts)', fontsize=11)
            axes[i].set_title(f"Scan {i+1}", fontsize=12)
            axes[i].grid(True, alpha=0.3)
            
            # Add temperature info if available
            if 'temperature' in metadata:
                axes[i].text(0.02, 0.98, f"T = {metadata['temperature']:.1f} K",
                           transform=axes[i].transAxes,
                           verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.show()
        
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        print("Please update the filepath in the script.")


# Example 3: Extracting and saving data
def example_export_data():
    """Example showing how to export data to numpy array or text file."""
    print("\nExample 3: Exporting data")
    print("-" * 50)
    
    filepath = "path/to/your/data.xrdml"
    
    try:
        reader = XRDMLReader(filepath)
        
        # Get data as numpy array
        data_array = reader.to_array(scan_index=0, normalize=False)
        
        print(f"Data shape: {data_array.shape}")
        print(f"First few rows:\n{data_array[:5]}")
        
        # Save to text file
        output_file = "xrd_data_exported.txt"
        np.savetxt(output_file, data_array, 
                  header="2theta(deg) Intensity(counts)",
                  fmt='%.4f')
        print(f"\nData saved to: {output_file}")
        
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        print("Please update the filepath in the script.")


# Example 4: Comparing wavelengths
def example_metadata():
    """Example focusing on metadata extraction."""
    print("\nExample 4: Detailed metadata extraction")
    print("-" * 50)
    
    filepath = "path/to/your/data.xrdml"
    
    try:
        reader = XRDMLReader(filepath)
        metadata = reader.get_metadata(0)
        
        print("\nInstrument Configuration:")
        print(f"  X-ray tube: {metadata.get('tube_type', 'N/A')}")
        print(f"  Voltage: {metadata.get('voltage', 'N/A')}")
        print(f"  Current: {metadata.get('current', 'N/A')}")
        print(f"  Goniometer radius: {metadata.get('gonio_radius', 'N/A')} mm")
        
        print("\nWavelength Information:")
        print(f"  Kα1: {metadata.get('wavelength_ka1', 'N/A')} Å")
        print(f"  Kα2: {metadata.get('wavelength_ka2', 'N/A')} Å")
        print(f"  Kβ: {metadata.get('wavelength_kb', 'N/A')} Å")
        print(f"  Kα2/Kα1 ratio: {metadata.get('ka2_ka1_ratio', 'N/A')}")
        
        print("\nMeasurement Details:")
        print(f"  Sample ID: {metadata.get('sample_id', 'N/A')}")
        print(f"  Start time: {metadata.get('start_time', 'N/A')}")
        print(f"  End time: {metadata.get('end_time', 'N/A')}")
        
        if 'temperature' in metadata:
            print(f"  Temperature: {metadata['temperature']} K")
        
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        print("Please update the filepath in the script.")


if __name__ == "__main__":
    print("XRDML Reader Examples")
    print("=" * 50)
    
    # Run examples (uncomment the ones you want to try)
    example_simple_read()
    # example_multi_scan()
    # example_export_data()
    # example_metadata()
