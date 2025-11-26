"""
XRDML Reader for Panalytical X-ray Diffraction Data

This module provides functionality to read and parse Panalytical .xrdml (XML) files,
which are commonly used for storing X-ray diffraction powder data.

Reference:
    Based on GSAS-II Panalytical reader implementation
    https://github.com/AdvancedPhotonSource/GSAS-II
"""

import xml.etree.ElementTree as ET
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union


class XRDMLReader:
    """
    A class for reading Panalytical .xrdml (XML) format XRD data files.
    
    The .xrdml format is an XML-based file format used by Panalytical (now Malvern Panalytical)
    X-ray diffraction instruments to store measurement data and metadata.
    
    Attributes:
        filepath (Path): Path to the .xrdml file
        root (ET.Element): Root element of the parsed XML tree
        metadata (dict): Dictionary containing measurement metadata
        scans (list): List of available scans in the file
    """
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize the XRDML reader.
        
        Args:
            filepath (str or Path): Path to the .xrdml file
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file is not a valid .xrdml file
        """
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        self.root = None
        self.metadata = {}
        self.scans = []
        self._namespace = ""
        
        # Validate and parse the file
        self._parse_file()
    
    def _parse_file(self):
        """
        Parse the XML file and extract the namespace.
        
        Raises:
            ValueError: If the file is not a valid .xrdml file
        """
        try:
            tree = ET.parse(self.filepath)
            self.root = tree.getroot()
            
            # Extract namespace from root tag
            tag = self.root.tag
            if '}' in tag:
                self._namespace = tag.split('}')[0] + '}'
            
            # Validate it's an XRDML file by checking for comment element
            if self.root.find(f'{self._namespace}comment') is None:
                # Try without namespace
                if self.root.find('comment') is None:
                    raise ValueError("Not a valid .xrdml file")
                self._namespace = ""
            
            # Find all scans
            xrd_measurements = self.root.findall(f'{self._namespace}xrdMeasurement')
            for measurement in xrd_measurements:
                scans = measurement.findall(f'{self._namespace}scan')
                self.scans.extend(scans)
                
        except ET.ParseError as e:
            raise ValueError(f"Invalid XML file: {e}")
    
    def get_metadata(self, scan_index: int = 0) -> Dict:
        """
        Extract metadata from the specified scan.
        
        Args:
            scan_index (int): Index of the scan to extract metadata from (0-based)
            
        Returns:
            dict: Dictionary containing metadata including:
                - sample_id: Sample identifier
                - start_time: Measurement start timestamp
                - end_time: Measurement end timestamp
                - wavelength_ka1: Kα1 wavelength
                - wavelength_ka2: Kα2 wavelength
                - wavelength_kb: Kβ wavelength
                - ka2_ka1_ratio: Ratio of Kα2 to Kα1
                - tube_type: X-ray tube type
                - voltage: Tube voltage
                - current: Tube current
                - gonio_radius: Goniometer radius
                - temperature: Sample temperature (if available)
        """
        if scan_index >= len(self.scans):
            raise IndexError(f"Scan index {scan_index} out of range. Available scans: {len(self.scans)}")
        
        metadata = {}
        scan = self.scans[scan_index]
        ns = self._namespace
        
        # Sample information
        sample = self.root.find(f'{ns}sample')
        if sample is not None:
            sample_id = sample.find(f'{ns}id')
            if sample_id is not None and sample_id.text:
                metadata['sample_id'] = sample_id.text
        
        # Measurement information
        xrd_measurement = self.root.find(f'{ns}xrdMeasurement')
        if xrd_measurement is not None:
            # Wavelength information
            wavelength = xrd_measurement.find(f'{ns}usedWavelength')
            if wavelength is not None:
                ka1 = wavelength.find(f'{ns}kAlpha1')
                ka2 = wavelength.find(f'{ns}kAlpha2')
                kb = wavelength.find(f'{ns}kBeta')
                ratio = wavelength.find(f'{ns}ratioKAlpha2KAlpha1')
                
                if ka1 is not None and ka1.text:
                    metadata['wavelength_ka1'] = float(ka1.text)
                if ka2 is not None and ka2.text:
                    metadata['wavelength_ka2'] = float(ka2.text)
                if kb is not None and kb.text:
                    metadata['wavelength_kb'] = float(kb.text)
                if ratio is not None and ratio.text:
                    metadata['ka2_ka1_ratio'] = float(ratio.text)
            
            # Incident beam path (goniometer radius and tube info)
            try:
                incident_beam = xrd_measurement.find(f'{ns}incidentBeamPath')
                if incident_beam is not None:
                    radius = incident_beam.find(f'{ns}radius')
                    if radius is not None and radius.text:
                        metadata['gonio_radius'] = float(radius.text)
                    
                    tube = incident_beam.find(f'{ns}xRayTube')
                    if tube is not None:
                        if 'name' in tube.attrib:
                            metadata['tube_type'] = tube.attrib['name']
                        
                        voltage = tube.find(f'{ns}tension')
                        current = tube.find(f'{ns}current')
                        if voltage is not None and voltage.text:
                            metadata['voltage'] = voltage.text
                        if current is not None and current.text:
                            metadata['current'] = current.text
            except:
                metadata['gonio_radius'] = 300.0  # Default value
        
        # Scan header information
        header = scan.find(f'{ns}header')
        if header is not None:
            start_time = header.find(f'{ns}startTimeStamp')
            end_time = header.find(f'{ns}endTimeStamp')
            
            if start_time is not None and start_time.text:
                metadata['start_time'] = start_time.text
            if end_time is not None and end_time.text:
                metadata['end_time'] = end_time.text
        
        # Non-ambient conditions (temperature, etc.)
        conditions = scan.find(f'{ns}nonAmbientPoints')
        if conditions is not None:
            kind = conditions.attrib.get('type')
            if kind == 'Temperature':
                values = conditions.find(f'{ns}nonAmbientValues')
                if values is not None and values.text:
                    # Take the last temperature value
                    temp_values = values.text.split()
                    if temp_values:
                        metadata['temperature'] = float(temp_values[-1])
        
        return metadata
    
    def read_scan(self, scan_index: int = 0) -> Tuple[np.ndarray, np.ndarray, Dict]:
        """
        Read XRD data from the specified scan.
        
        Args:
            scan_index (int): Index of the scan to read (0-based)
            
        Returns:
            tuple: (two_theta, intensities, metadata)
                - two_theta (np.ndarray): Array of 2θ values (degrees)
                - intensities (np.ndarray): Array of intensity counts
                - metadata (dict): Dictionary of measurement metadata
                
        Raises:
            IndexError: If scan_index is out of range
            ValueError: If intensity data cannot be found
        """
        if scan_index >= len(self.scans):
            raise IndexError(f"Scan index {scan_index} out of range. Available scans: {len(self.scans)}")
        
        scan = self.scans[scan_index]
        ns = self._namespace
        
        # Get data points
        data_points = scan.find(f'{ns}dataPoints')
        if data_points is None:
            raise ValueError("No dataPoints element found in scan")
        
        # Get positions (2-theta values)
        positions = data_points.find(f'{ns}positions')
        if positions is None:
            raise ValueError("No positions element found in dataPoints")
        
        start_pos_elem = positions.find(f'{ns}startPosition')
        end_pos_elem = positions.find(f'{ns}endPosition')
        
        if start_pos_elem is None or end_pos_elem is None:
            raise ValueError("Start or end position not found")
        
        start_pos = float(start_pos_elem.text)
        end_pos = float(end_pos_elem.text)
        
        # Get intensities (try both 'intensities' and 'counts' tags)
        intensities = None
        for tag in ['intensities', 'counts']:
            intensity_elem = data_points.find(f'{ns}{tag}')
            if intensity_elem is not None and intensity_elem.text:
                intensities = np.fromstring(intensity_elem.text, sep=' ')
                break
        
        if intensities is None:
            raise ValueError("Intensity data could not be located")
        
        # Generate 2-theta values
        n_points = len(intensities)
        two_theta = np.linspace(start_pos, end_pos, n_points)
        
        # Get metadata
        metadata = self.get_metadata(scan_index)
        
        return two_theta, intensities, metadata
    
    def read_all_scans(self) -> List[Tuple[np.ndarray, np.ndarray, Dict]]:
        """
        Read all scans from the file.
        
        Returns:
            list: List of tuples, each containing (two_theta, intensities, metadata)
        """
        results = []
        for i in range(len(self.scans)):
            results.append(self.read_scan(i))
        return results
    
    def get_num_scans(self) -> int:
        """
        Get the number of scans in the file.
        
        Returns:
            int: Number of scans available
        """
        return len(self.scans)
    
    def to_array(self, scan_index: int = 0, normalize: bool = False) -> np.ndarray:
        """
        Read scan data and return as a 2D numpy array.
        
        Args:
            scan_index (int): Index of the scan to read (0-based)
            normalize (bool): If True, normalize intensities to [0, 1]
            
        Returns:
            np.ndarray: Array of shape (N, 2) where columns are [2θ, intensity]
        """
        two_theta, intensities, _ = self.read_scan(scan_index)
        
        if normalize:
            max_intensity = np.max(intensities)
            if max_intensity > 0:
                intensities = intensities / max_intensity
        
        return np.column_stack((two_theta, intensities))
    
    def __repr__(self) -> str:
        """String representation of the reader."""
        return f"XRDMLReader('{self.filepath.name}', scans={len(self.scans)})"


def read_xrdml(filepath: Union[str, Path], scan_index: int = 0, 
               normalize: bool = False) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Convenience function to read a .xrdml file.
    
    Args:
        filepath (str or Path): Path to the .xrdml file
        scan_index (int): Index of the scan to read (0-based)
        normalize (bool): If True, normalize intensities to [0, 1]
        
    Returns:
        tuple: (two_theta, intensities, metadata)
            - two_theta (np.ndarray): Array of 2θ values (degrees)
            - intensities (np.ndarray): Array of intensity counts
            - metadata (dict): Dictionary of measurement metadata
            
    Example:
        >>> two_theta, intensity, metadata = read_xrdml('sample.xrdml')
        >>> print(f"Wavelength: {metadata['wavelength_ka1']} Å")
        >>> print(f"Sample: {metadata['sample_id']}")
    """
    reader = XRDMLReader(filepath)
    two_theta, intensities, metadata = reader.read_scan(scan_index)
    
    if normalize:
        max_intensity = np.max(intensities)
        if max_intensity > 0:
            intensities = intensities / max_intensity
    
    return two_theta, intensities, metadata
