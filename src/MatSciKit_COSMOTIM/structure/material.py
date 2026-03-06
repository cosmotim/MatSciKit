"""
Material properties container and database loader.

Provides the :class:`Material` dataclass that holds crystallographic and
elastic properties needed by Pipelines 1 and 2, plus utilities for
loading materials from JSON databases or CIF files.

Part of Pipeline 0 (Structure).
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np

from ..constants import kb, hbar, AMU_TO_KG


# Average atomic masses (amu) for common solid electrolytes
# Calculated from stoichiometric composition
_DEFAULT_AVG_MASSES: Dict[str, float] = {
    'LAGP': (6.941*1.5 + 26.982*0.5 + 72.630*1.5 + 30.974*3 + 15.999*12) / (1.5+0.5+1.5+3+12),
    'NZP': (22.990 + 91.224*2 + 30.974*3 + 15.999*12) / (1+2+3+12),
    'LLZTO': (6.941*6.5 + 138.91*3 + 91.224*1.5 + 180.95*0.5 + 15.999*12) / (6.5+3+1.5+0.5+12),
    'LSHT': (6.941*3/8 + 87.62*7/16 + 178.49*1/4 + 180.95*3/4 + 15.999*3) / (3/8+7/16+1/4+3/4+3),
    'LGPS': (6.941*10 + 72.630 + 30.974*2 + 32.06*12) / (10+1+2+12),
    'LPSCl': (6.941*6 + 30.974 + 32.06*5 + 35.45) / (6+1+5+1),
    'Na3PS4': (22.990*3 + 30.974 + 32.06*4) / (3+1+4),
    'Li3InCl6': (6.941*3 + 114.82 + 35.45*6) / (3+1+6),
}


@dataclass
class Material:
    """
    Container for material properties.

    Holds crystallographic, elastic, and thermal properties needed
    for Cahill minimum κ, mean free path, and other calculations.

    Parameters
    ----------
    name : str
        Short name (e.g. 'LLZTO', 'LGPS').
    formula : str
        Chemical formula.
    density : float
        Mass density (kg/m³).
    n_density : float
        Number density (atoms/m³).
    v_avg : float, optional
        Average sound velocity (m/s).
    v_L : float, optional
        Longitudinal sound velocity (m/s).
    v_T : float, optional
        Transverse sound velocity (m/s).
    theta_D : float, optional
        Debye temperature (K). Estimated from v_avg if not provided.
    volume : float, optional
        Unit cell volume (m³).
    category : str, optional
        Material class (e.g. 'oxide', 'sulfide', 'halide').
    sources : list of str, optional
        Literature references.

    Examples
    --------
    >>> mat = Material(
    ...     name='LSHT',
    ...     formula='Li3/8Sr7/16Hf1/4Ta3/4O3',
    ...     density=6870,
    ...     n_density=4.8125 / (3.98e-10)**3,
    ...     v_avg=3461.3,
    ... )
    >>> print(f"θ_D = {mat.theta_D:.0f} K")
    """

    name: str
    formula: str
    density: float                     # kg/m³
    n_density: float                   # atoms/m³
    v_avg: Optional[float] = None      # m/s
    v_L: Optional[float] = None        # m/s
    v_T: Optional[float] = None        # m/s
    theta_D: Optional[float] = None    # K
    volume: Optional[float] = None     # m³
    category: Optional[str] = None
    sources: List[str] = field(default_factory=list)

    def __post_init__(self):
        # Estimate θ_D from v_avg and n_density if not provided
        if self.theta_D is None and self.v_avg is not None:
            self.theta_D = estimate_debye_temperature(self.v_avg, self.n_density)

    @classmethod
    def from_cif(cls, path: Union[str, Path], name: str = '',
                 v_avg: Optional[float] = None, **kwargs) -> 'Material':
        """
        Create a Material from a CIF file.

        Parameters
        ----------
        path : str or Path
            Path to the CIF file.
        name : str, optional
            Short name for the material.
        v_avg : float, optional
            Average sound velocity (m/s). Must be provided externally
            (not in CIF files).
        **kwargs
            Additional keyword arguments passed to Material.

        Returns
        -------
        Material
        """
        from . import cif_reader

        result = cif_reader.read(path)

        return cls(
            name=name or Path(path).stem,
            formula=name,
            density=result['density'],
            n_density=result['n_density'],
            volume=result['volume'] * 1e-30,
            v_avg=v_avg,
            **kwargs,
        )

    @classmethod
    def from_params(cls, name: str, formula: str,
                    volume_m3: float, total_atoms: float,
                    total_mass_kg: float,
                    v_avg: Optional[float] = None, **kwargs) -> 'Material':
        """
        Create a Material from manual parameters.

        Parameters
        ----------
        name : str
            Short name.
        formula : str
            Chemical formula.
        volume_m3 : float
            Unit cell volume (m³).
        total_atoms : float
            Atoms per unit cell (may be fractional for partial occupancy).
        total_mass_kg : float
            Total mass per unit cell (kg).
        v_avg : float, optional
            Average sound velocity (m/s).
        """
        density = total_mass_kg / volume_m3
        n_density = total_atoms / volume_m3
        return cls(
            name=name,
            formula=formula,
            density=density,
            n_density=n_density,
            volume=volume_m3,
            v_avg=v_avg,
            **kwargs,
        )


def estimate_debye_temperature(v_s: float, n_density: float) -> float:
    """
    Estimate Debye temperature from sound velocity and number density.

    .. math::

        \\theta_D = \\frac{\\hbar}{k_B} v_s (6\\pi^2 n)^{1/3}

    Parameters
    ----------
    v_s : float
        Average sound velocity (m/s).
    n_density : float
        Number density (atoms/m³).

    Returns
    -------
    theta_D : float
        Estimated Debye temperature (K).
    """
    return (hbar / kb) * v_s * (6 * np.pi**2 * n_density) ** (1.0 / 3)


def estimate_n_density(density: float, avg_mass_amu: float) -> float:
    """
    Estimate number density from mass density and average atomic mass.

    Parameters
    ----------
    density : float
        Mass density (kg/m³).
    avg_mass_amu : float
        Average atomic mass (amu).

    Returns
    -------
    n_density : float
        Number density (atoms/m³).
    """
    return density / (avg_mass_amu * AMU_TO_KG)


def load_database(path: Union[str, Path]) -> List[Material]:
    """
    Load materials from a sound velocity database JSON file.

    Reads the database format used by the researcher agent and
    returns a list of :class:`Material` objects for all materials
    with complete data (density + sound velocity).

    Parameters
    ----------
    path : str or Path
        Path to the JSON database file.

    Returns
    -------
    materials : list of Material
        Materials with complete data.

    Examples
    --------
    >>> materials = load_database("sound_velocity_database.json")
    >>> for m in materials:
    ...     print(f"{m.name}: v_avg={m.v_avg} m/s, ρ={m.density} kg/m³")
    """
    path = Path(path)
    with open(path) as f:
        db = json.load(f)

    materials = []

    for category, mats in db.get('materials', {}).items():
        for name, data in mats.items():
            sv = data.get('sound_velocities', {})
            v_avg = sv.get('average_m_per_s')
            density = data.get('density_kg_per_m3')

            if v_avg is None or density is None:
                continue

            # Get number density from avg mass if available
            avg_mass = _DEFAULT_AVG_MASSES.get(name)
            if avg_mass is not None:
                n_density = estimate_n_density(density, avg_mass)
            else:
                # Fallback: rough estimate assuming avg mass ~30 amu
                n_density = estimate_n_density(density, 30.0)

            # Collect sources
            sources = []
            for src in data.get('sources', []):
                if isinstance(src, dict):
                    sources.append(src.get('citation', ''))
                elif isinstance(src, str):
                    sources.append(src)

            mat = Material(
                name=name,
                formula=data.get('full_name', name),
                density=density,
                n_density=n_density,
                v_avg=v_avg,
                v_L=sv.get('longitudinal_m_per_s'),
                v_T=sv.get('transverse_m_per_s'),
                category=category,
                sources=sources,
            )
            materials.append(mat)

    return materials


def batch_cahill(materials: List[Material],
                 T: Union[float, np.ndarray] = 300.0) -> Dict[str, Union[float, np.ndarray]]:
    """
    Compute Cahill minimum κ for a list of materials.

    Parameters
    ----------
    materials : list of Material
        Materials with v_avg and n_density set.
    T : float or np.ndarray
        Temperature(s) (K).

    Returns
    -------
    results : dict
        Mapping of material name → κ_min value(s).

    Examples
    --------
    >>> mats = load_database("sound_velocity_database.json")
    >>> results = batch_cahill(mats, T=300)
    >>> for name, kappa in sorted(results.items(), key=lambda x: -x[1]):
    ...     print(f"{name}: κ_min = {kappa:.3f} W/(m·K)")
    """
    from ..thermal_conductivity import cahill

    results = {}
    for mat in materials:
        if mat.v_avg is None or mat.theta_D is None:
            continue
        kappa = cahill.minimum_tc(T, mat.n_density, mat.theta_D, mat.v_avg)
        results[mat.name] = kappa

    return results
