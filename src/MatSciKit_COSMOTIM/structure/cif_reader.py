"""
CIF file reader for extracting crystallographic material properties.

Reads Crystallographic Information Files (.cif) and computes material
properties needed for thermal transport analysis (Pipeline 0).

Based on cif_analyzer.py by the researcher agent, adapted to MatSciKit
conventions with bug fixes to the basic parser.

Notes
-----
Uses a built-in parser by default. If ``pymatgen`` is installed, it will
be used for more robust parsing of complex CIF files.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from ..constants import AMU_TO_KG


# Atomic masses (amu) — common elements in solid electrolytes and oxides
ATOMIC_MASSES: Dict[str, float] = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81,
    'C': 12.01, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906,
    'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87,
    'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60,
    'I': 126.90, 'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12,
    'Pr': 140.91, 'Nd': 144.24, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25,
    'Tb': 158.93, 'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93,
    'Yb': 173.04, 'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84,
    'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97,
    'Pb': 207.2, 'Bi': 208.98,
}

# Try importing pymatgen for advanced parsing
try:
    from pymatgen.io.cif import CifParser as _PmgCifParser
    from pymatgen.core.periodic_table import Element as _PmgElement
    _PYMATGEN_AVAILABLE = True
except ImportError:
    _PYMATGEN_AVAILABLE = False


def _parse_float(s: str) -> float:
    """Parse a CIF numeric value, stripping parenthesized uncertainties."""
    return float(s.split('(')[0])


def _clean_symbol(s: str) -> str:
    """Extract element symbol from a CIF atom label (e.g. 'Sr1' → 'Sr')."""
    alpha = ''.join(c for c in s if c.isalpha())
    # Try 2-char first, then 1-char
    if len(alpha) >= 2 and alpha[:2] in ATOMIC_MASSES:
        return alpha[:2]
    if len(alpha) >= 1 and alpha[0] in ATOMIC_MASSES:
        return alpha[0]
    return alpha


def _parse_cif_basic(path: str) -> Dict:
    """
    Basic CIF parser (no pymatgen dependency).

    Handles standard CIF loop_ blocks for atom sites.
    """
    data: Dict = {}

    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # --- Extract scalar cell parameters ---
    for line in lines:
        line = line.strip()
        parts = line.split()
        if len(parts) < 2:
            continue
        tag, val = parts[0], parts[1]
        if tag == '_cell_length_a':
            data['a'] = _parse_float(val)
        elif tag == '_cell_length_b':
            data['b'] = _parse_float(val)
        elif tag == '_cell_length_c':
            data['c'] = _parse_float(val)
        elif tag == '_cell_angle_alpha':
            data['alpha'] = _parse_float(val)
        elif tag == '_cell_angle_beta':
            data['beta'] = _parse_float(val)
        elif tag == '_cell_angle_gamma':
            data['gamma'] = _parse_float(val)
        elif tag == '_cell_volume':
            data['volume'] = _parse_float(val)

    # --- Parse atom_site loop_ blocks ---
    atoms: List[Dict] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Find a loop_ block
        if line == 'loop_':
            i += 1
            # Collect header tags
            headers: List[str] = []
            while i < len(lines):
                hline = lines[i].strip()
                if hline.startswith('_'):
                    headers.append(hline.split()[0])
                    i += 1
                else:
                    break

            # Check if this is an atom_site loop
            if not any('_atom_site' in h for h in headers):
                continue

            # Determine column indices
            col_label = None
            col_symbol = None
            col_occ = None
            for ci, h in enumerate(headers):
                if h == '_atom_site_type_symbol':
                    col_symbol = ci
                elif h == '_atom_site_label':
                    col_label = ci
                elif h == '_atom_site_occupancy':
                    col_occ = ci

            # Read data rows
            while i < len(lines):
                dline = lines[i].strip()
                if not dline or dline.startswith('_') or dline.startswith('loop_') or dline.startswith('#'):
                    break
                parts = dline.split()
                if len(parts) < len(headers):
                    break

                # Get element symbol
                symbol = None
                if col_symbol is not None:
                    symbol = _clean_symbol(parts[col_symbol])
                elif col_label is not None:
                    symbol = _clean_symbol(parts[col_label])

                # Get occupancy
                occ = 1.0
                if col_occ is not None:
                    try:
                        occ = _parse_float(parts[col_occ])
                    except (ValueError, IndexError):
                        occ = 1.0

                if symbol and symbol in ATOMIC_MASSES:
                    atoms.append({'symbol': symbol, 'occupancy': occ})

                i += 1
        else:
            i += 1

    data['atoms'] = atoms
    return data


def _calculate_volume(a: float, b: float, c: float,
                      alpha: float, beta: float, gamma: float) -> float:
    """Calculate unit cell volume from lattice parameters."""
    ar = math.radians(alpha)
    br = math.radians(beta)
    gr = math.radians(gamma)
    return a * b * c * math.sqrt(
        1 + 2 * math.cos(ar) * math.cos(br) * math.cos(gr)
        - math.cos(ar)**2 - math.cos(br)**2 - math.cos(gr)**2
    )


def read(path: Union[str, Path], use_pymatgen: Optional[bool] = None) -> Dict:
    """
    Read a CIF file and extract material properties.

    Parameters
    ----------
    path : str or Path
        Path to the CIF file.
    use_pymatgen : bool, optional
        Force pymatgen parser (True) or basic parser (False).
        If None (default), uses pymatgen when available.

    Returns
    -------
    result : dict
        Dictionary with keys:

        - ``'file'`` : str — filename
        - ``'unit_cell'`` : dict — {a, b, c, alpha, beta, gamma}
        - ``'volume'`` : float — unit cell volume (Å³)
        - ``'atoms'`` : list of dict — [{symbol, occupancy}, ...]
        - ``'total_atoms'`` : float — total atoms per unit cell
          (sum of occupancies)
        - ``'total_mass'`` : float — total mass per unit cell (kg)
        - ``'density'`` : float — mass density (kg/m³)
        - ``'n_density'`` : float — number density (atoms/m³)
        - ``'average_mass'`` : float — average atomic mass (amu)
        - ``'parser'`` : str — 'pymatgen' or 'basic'

    Raises
    ------
    FileNotFoundError
        If the CIF file does not exist.
    ValueError
        If essential crystallographic data cannot be extracted.

    Examples
    --------
    >>> from MatSciKit_COSMOTIM.structure import cif_reader
    >>> mat = cif_reader.read("SrTiO3.cif")
    >>> print(f"Density: {mat['density']:.0f} kg/m³")
    Density: 5117 kg/m³
    >>> print(f"N/V: {mat['n_density']:.2e} atoms/m³")
    N/V: 8.39e+28 atoms/m³
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"CIF file not found: {path}")

    should_use_pmg = use_pymatgen if use_pymatgen is not None else _PYMATGEN_AVAILABLE

    result: Dict = {
        'file': path.name,
        'parser': 'pymatgen' if (should_use_pmg and _PYMATGEN_AVAILABLE) else 'basic',
    }

    if should_use_pmg and _PYMATGEN_AVAILABLE:
        parser = _PmgCifParser(str(path))
        # Use parse_structures (new API) if available, else get_structures
        if hasattr(parser, 'parse_structures'):
            structure = parser.parse_structures(primitive=True)[0]
        else:
            structure = parser.get_structures()[0]
        lattice = structure.lattice

        result['unit_cell'] = {
            'a': lattice.a, 'b': lattice.b, 'c': lattice.c,
            'alpha': lattice.alpha, 'beta': lattice.beta, 'gamma': lattice.gamma,
        }
        result['volume'] = lattice.volume

        atoms = []
        total_atoms = 0.0
        total_mass_amu = 0.0
        for element, amount in structure.composition.items():
            sym = str(element)
            atoms.append({'symbol': sym, 'occupancy': amount})
            total_atoms += amount
            mass = ATOMIC_MASSES.get(sym, float(_PmgElement(sym).atomic_mass))
            total_mass_amu += amount * mass

        result['atoms'] = atoms

    elif should_use_pmg and not _PYMATGEN_AVAILABLE:
        raise ImportError(
            "pymatgen requested but not installed. "
            "Install with: pip install pymatgen"
        )

    else:
        # Basic parser
        data = _parse_cif_basic(str(path))

        a = data.get('a', 0)
        b = data.get('b', 0)
        c = data.get('c', 0)
        alpha = data.get('alpha', 90)
        beta = data.get('beta', 90)
        gamma = data.get('gamma', 90)

        result['unit_cell'] = {
            'a': a, 'b': b, 'c': c,
            'alpha': alpha, 'beta': beta, 'gamma': gamma,
        }
        result['volume'] = data.get('volume', _calculate_volume(a, b, c, alpha, beta, gamma))

        atoms = data.get('atoms', [])
        result['atoms'] = atoms

        total_atoms = sum(a['occupancy'] for a in atoms)
        total_mass_amu = sum(
            a['occupancy'] * ATOMIC_MASSES.get(a['symbol'], 0)
            for a in atoms
        )

    result['total_atoms'] = total_atoms

    if total_atoms == 0:
        raise ValueError(
            f"No atoms found in {path.name}. "
            "Check CIF format or install pymatgen for better parsing."
        )

    result['average_mass'] = total_mass_amu / total_atoms  # amu
    result['total_mass'] = total_mass_amu * AMU_TO_KG  # kg

    volume_m3 = result['volume'] * 1e-30  # Å³ → m³
    result['density'] = result['total_mass'] / volume_m3  # kg/m³
    result['n_density'] = total_atoms / volume_m3  # atoms/m³

    return result


def get_material_params(path: Union[str, Path],
                        use_pymatgen: Optional[bool] = None) -> Tuple[float, float, float]:
    """
    Convenience function: extract (volume_m3, density, n_density) from a CIF.

    These are the three material parameters needed by Pipelines 1 and 2.

    Parameters
    ----------
    path : str or Path
        Path to the CIF file.
    use_pymatgen : bool, optional
        Force parser choice (see :func:`read`).

    Returns
    -------
    volume : float
        Unit cell volume (m³).
    density : float
        Mass density (kg/m³).
    n_density : float
        Number density (atoms/m³).

    Examples
    --------
    >>> V, rho, n = cif_reader.get_material_params("SrTiO3.cif")
    """
    result = read(path, use_pymatgen=use_pymatgen)
    volume_m3 = result['volume'] * 1e-30
    return volume_m3, result['density'], result['n_density']
