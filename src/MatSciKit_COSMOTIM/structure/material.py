"""
Material properties container and database loader.

Provides the :class:`Material` dataclass that holds crystallographic and
elastic properties needed by Pipelines 1 and 2, plus utilities for
loading materials from JSON databases, CIF files, or the Materials
Project API.

Field naming follows the Materials Project convention where applicable
(e.g. ``bulk_modulus``, ``shear_modulus``, ``sound_velocity``).

Part of Pipeline 0 (Structure).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..constants import AMU_TO_KG, hbar, kb

# Average atomic masses (amu) for common solid electrolytes
_DEFAULT_AVG_MASSES: dict[str, float] = {
    "LAGP": (6.941 * 1.5 + 26.982 * 0.5 + 72.630 * 1.5 + 30.974 * 3 + 15.999 * 12)
    / (1.5 + 0.5 + 1.5 + 3 + 12),
    "NZP": (22.990 + 91.224 * 2 + 30.974 * 3 + 15.999 * 12) / (1 + 2 + 3 + 12),
    "LLZTO": (6.941 * 6.5 + 138.91 * 3 + 91.224 * 1.5 + 180.95 * 0.5 + 15.999 * 12)
    / (6.5 + 3 + 1.5 + 0.5 + 12),
    "LSHT": (6.941 * 3 / 8 + 87.62 * 7 / 16 + 178.49 * 1 / 4 + 180.95 * 3 / 4 + 15.999 * 3)
    / (3 / 8 + 7 / 16 + 1 / 4 + 3 / 4 + 3),
    "LGPS": (6.941 * 10 + 72.630 + 30.974 * 2 + 32.06 * 12) / (10 + 1 + 2 + 12),
    "LPSCl": (6.941 * 6 + 30.974 + 32.06 * 5 + 35.45) / (6 + 1 + 5 + 1),
    "Na3PS4": (22.990 * 3 + 30.974 + 32.06 * 4) / (3 + 1 + 4),
    "Li3InCl6": (6.941 * 3 + 114.82 + 35.45 * 6) / (3 + 1 + 6),
}


@dataclass
class Material:
    """
    Container for material properties.

    Field names follow the Materials Project convention where applicable.
    Holds crystallographic, elastic, and thermal properties needed
    for Cahill minimum κ, mean free path, and other calculations.

    Parameters
    ----------
    name : str
        Short name (e.g. 'LLZTO', 'LGPS').
    formula : str
        Chemical formula (pretty-printed).
    material_id : str, optional
        Materials Project ID (e.g. 'mp-149').
    density : float
        Mass density (kg/m³). Same as MP ``density``.
    density_atomic : float, optional
        Atomic density (atoms/Å³). Same as MP ``density_atomic``.
    n_density : float
        Number density (atoms/m³). Derived from ``density_atomic``
        or calculated from structure.
    volume : float, optional
        Unit cell volume (ų). Same as MP ``volume``.
    nsites : int, optional
        Number of sites in the unit cell. Same as MP ``nsites``.

    sound_velocity : dict, optional
        Sound velocities (m/s). Keys follow MP ElasticityDoc:

        - ``'transverse'``: transverse (shear) wave velocity
        - ``'longitudinal'``: longitudinal (pressure) wave velocity
        - ``'snyder_acoustic'``: Snyder acoustic velocity
        - ``'snyder_optical'``: Snyder optical velocity
        - ``'snyder_total'``: Snyder total velocity

    bulk_modulus : dict, optional
        Bulk modulus (GPa). Keys: ``'voigt'``, ``'reuss'``, ``'vrh'``.
    shear_modulus : dict, optional
        Shear modulus (GPa). Keys: ``'voigt'``, ``'reuss'``, ``'vrh'``.
    young_modulus : float, optional
        Young's modulus (GPa).
    homogeneous_poisson : float, optional
        Poisson's ratio.
    universal_anisotropy : float, optional
        Universal anisotropy index.

    debye_temperature : float, optional
        Debye temperature (K). Estimated from v_s if not provided.
    thermal_conductivity : dict, optional
        Thermal conductivity estimates (W/(m·K)).
        Keys: ``'clarke'``, ``'cahill'``.

    category : str, optional
        Material class (e.g. 'oxide', 'sulfide', 'halide').
    symmetry : dict, optional
        Symmetry info (space group, crystal system, etc.).
    elements : list of str, optional
        Constituent elements.
    sources : list of str, optional
        Literature references.

    Examples
    --------
    >>> mat = Material(
    ...     name='LSHT',
    ...     formula='Li₃/₈Sr₇/₁₆Hf₁/₄Ta₃/₄O₃',
    ...     density=6870,
    ...     n_density=4.8125 / (3.98e-10)**3,
    ...     sound_velocity={'transverse': 1900, 'longitudinal': 3400},
    ... )
    >>> print(f"θ_D = {mat.debye_temperature:.0f} K")
    """

    # --- Identity ---
    name: str
    formula: str
    material_id: str | None = None

    # --- Structure ---
    density: float = 0.0  # kg/m³
    density_atomic: float | None = None  # atoms/ų
    n_density: float = 0.0  # atoms/m³
    volume: float | None = None  # ų
    nsites: int | None = None

    # --- Sound velocities (m/s) ---
    sound_velocity: dict[str, float] | None = None

    # --- Elastic moduli ---
    bulk_modulus: dict[str, float] | None = None  # GPa
    shear_modulus: dict[str, float] | None = None  # GPa
    young_modulus: float | None = None  # GPa
    homogeneous_poisson: float | None = None
    universal_anisotropy: float | None = None

    # --- Thermal ---
    debye_temperature: float | None = None  # K
    thermal_conductivity: dict[str, float] | None = None  # W/(m·K)

    # --- Metadata ---
    category: str | None = None
    symmetry: dict[str, Any] | None = None
    elements: list[str] | None = None
    sources: list[str] = field(default_factory=list)

    def __post_init__(self):
        # Derive n_density from density_atomic if available
        if self.n_density == 0 and self.density_atomic is not None:
            self.n_density = self.density_atomic * 1e30  # atoms/ų → atoms/m³

        # Estimate Debye temperature from sound velocity if not provided
        if self.debye_temperature is None and self.v_avg is not None and self.n_density > 0:
            self.debye_temperature = estimate_debye_temperature(self.v_avg, self.n_density)

    # --- Convenience properties ---

    @property
    def v_avg(self) -> float | None:
        """Average (Debye) sound velocity (m/s).

        Priority:
        1. Explicit ``'average'`` key (from local database)
        2. Debye average computed from v_L and v_T:
           v_avg = [1/3 (1/v_L³ + 2/v_T³)]^(-1/3)
        3. ``'snyder_total'`` as fallback

        Note: MP's ``snyder_acoustic`` is the Snyder model velocity
        for thermoelectric κ_L, NOT the Debye average. We compute
        the Debye average from v_L/v_T instead.
        """
        if self.sound_velocity is None:
            return None

        # Explicit average from our database format
        if self.sound_velocity.get("average"):
            return self.sound_velocity["average"]

        # Compute Debye average from v_L and v_T (standard for Cahill model)
        v_L = self.sound_velocity.get("longitudinal")
        v_T = self.sound_velocity.get("transverse")
        if v_L and v_T:
            return (1.0 / 3 * (1 / v_L**3 + 2 / v_T**3)) ** (-1.0 / 3)

        return None

    @property
    def v_L(self) -> float | None:
        """Longitudinal sound velocity (m/s)."""
        if self.sound_velocity is None:
            return None
        return self.sound_velocity.get("longitudinal")

    @property
    def v_T(self) -> float | None:
        """Transverse sound velocity (m/s)."""
        if self.sound_velocity is None:
            return None
        return self.sound_velocity.get("transverse")

    @property
    def K_vrh(self) -> float | None:
        """Voigt-Reuss-Hill bulk modulus (GPa)."""
        if self.bulk_modulus is None:
            return None
        return self.bulk_modulus.get("vrh")

    @property
    def G_vrh(self) -> float | None:
        """Voigt-Reuss-Hill shear modulus (GPa)."""
        if self.shear_modulus is None:
            return None
        return self.shear_modulus.get("vrh")

    # --- Factory methods ---

    @classmethod
    def from_mp(cls, mp_doc: Any, name: str | None = None) -> Material:
        """
        Create a Material from a Materials Project API document.

        Works with both ``SummaryDoc`` (from ``mpr.materials.summary``)
        and ``ElasticityDoc`` (from ``mpr.materials.elasticity``).

        Parameters
        ----------
        mp_doc : SummaryDoc or ElasticityDoc
            A document returned by the mp-api client.
        name : str, optional
            Short name. Defaults to the formula.

        Returns
        -------
        Material

        Examples
        --------
        >>> from mp_api.client import MPRester
        >>> with MPRester("your_api_key") as mpr:
        ...     docs = mpr.materials.elasticity.search(
        ...         material_ids=["mp-5827"],  # SrTiO3
        ...         fields=["material_id", "formula_pretty", "density",
        ...                 "volume", "nsites", "sound_velocity",
        ...                 "bulk_modulus", "shear_modulus",
        ...                 "debye_temperature", "thermal_conductivity",
        ...                 "homogeneous_poisson", "universal_anisotropy"]
        ...     )
        ...     mat = Material.from_mp(docs[0], name="STO")
        """
        formula = getattr(mp_doc, "formula_pretty", "") or ""
        mpid = str(getattr(mp_doc, "material_id", "") or "")
        density = getattr(mp_doc, "density", None) or 0.0
        density_atomic = getattr(mp_doc, "density_atomic", None)
        volume = getattr(mp_doc, "volume", None)
        nsites = getattr(mp_doc, "nsites", None)

        # Sound velocity
        sv_obj = getattr(mp_doc, "sound_velocity", None)
        sound_velocity = None
        if sv_obj is not None:
            sound_velocity = {}
            for key in [
                "transverse",
                "longitudinal",
                "snyder_acoustic",
                "snyder_optical",
                "snyder_total",
            ]:
                val = getattr(sv_obj, key, None)
                if val is not None:
                    sound_velocity[key] = val

        # Elastic moduli
        bulk = getattr(mp_doc, "bulk_modulus", None)
        bulk_modulus = None
        if bulk is not None:
            bulk_modulus = {}
            for key in ["voigt", "reuss", "vrh"]:
                val = getattr(bulk, key, None)
                if val is not None:
                    bulk_modulus[key] = val

        shear = getattr(mp_doc, "shear_modulus", None)
        shear_modulus = None
        if shear is not None:
            shear_modulus = {}
            for key in ["voigt", "reuss", "vrh"]:
                val = getattr(shear, key, None)
                if val is not None:
                    shear_modulus[key] = val

        # Thermal conductivity
        tc_obj = getattr(mp_doc, "thermal_conductivity", None)
        thermal_conductivity = None
        if tc_obj is not None:
            thermal_conductivity = {}
            for key in ["clarke", "cahill"]:
                val = getattr(tc_obj, key, None)
                if val is not None:
                    thermal_conductivity[key] = val

        # Number density
        n_density = 0.0
        if density_atomic is not None:
            n_density = density_atomic * 1e30
        elif nsites and volume:
            n_density = nsites / (volume * 1e-30)

        # Elements
        elements = None
        elem_list = getattr(mp_doc, "elements", None)
        if elem_list:
            elements = [str(e) for e in elem_list]

        # Symmetry
        sym_obj = getattr(mp_doc, "symmetry", None)
        symmetry = None
        if sym_obj is not None:
            symmetry = {
                "crystal_system": getattr(sym_obj, "crystal_system", None),
                "symbol": getattr(sym_obj, "symbol", None),
                "number": getattr(sym_obj, "number", None),
            }

        return cls(
            name=name or formula,
            formula=formula,
            material_id=mpid,
            density=density * 1000 if density and density < 100 else density,
            # MP density is in g/cm³, convert to kg/m³
            density_atomic=density_atomic,
            n_density=n_density,
            volume=volume,
            nsites=nsites,
            sound_velocity=sound_velocity,
            bulk_modulus=bulk_modulus,
            shear_modulus=shear_modulus,
            young_modulus=(
                getattr(mp_doc, "youngs_modulus", None) or getattr(mp_doc, "young_modulus", None)
            ),
            homogeneous_poisson=getattr(mp_doc, "homogeneous_poisson", None),
            universal_anisotropy=getattr(mp_doc, "universal_anisotropy", None),
            debye_temperature=getattr(mp_doc, "debye_temperature", None),
            thermal_conductivity=thermal_conductivity,
            elements=elements,
            symmetry=symmetry,
            sources=[f"Materials Project: {mpid}"],
        )

    @classmethod
    def from_cif(
        cls,
        path: str | Path,
        name: str = "",
        sound_velocity: dict[str, float] | None = None,
        **kwargs,
    ) -> Material:
        """
        Create a Material from a CIF file.

        Parameters
        ----------
        path : str or Path
            Path to the CIF file.
        name : str, optional
            Short name for the material.
        sound_velocity : dict, optional
            Sound velocities (m/s). CIF files don't contain these;
            must be provided from experiment or literature.
        **kwargs
            Additional keyword arguments passed to Material.

        Returns
        -------
        Material
        """
        from . import cif_reader

        result = cif_reader.read(path)
        volume_ang3 = result["volume"]

        return cls(
            name=name or Path(path).stem,
            formula=name or Path(path).stem,
            density=result["density"],
            n_density=result["n_density"],
            volume=volume_ang3,
            nsites=round(result["total_atoms"]),
            sound_velocity=sound_velocity,
            **kwargs,
        )

    @classmethod
    def from_params(
        cls,
        name: str,
        formula: str,
        volume_m3: float,
        total_atoms: float,
        total_mass_kg: float,
        sound_velocity: dict[str, float] | None = None,
        **kwargs,
    ) -> Material:
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
            Atoms per unit cell (may be fractional).
        total_mass_kg : float
            Total mass per unit cell (kg).
        sound_velocity : dict, optional
            Sound velocities (m/s).
        """
        density = total_mass_kg / volume_m3
        n_density = total_atoms / volume_m3
        return cls(
            name=name,
            formula=formula,
            density=density,
            n_density=n_density,
            volume=volume_m3 * 1e30,  # m³ → ų
            sound_velocity=sound_velocity,
            **kwargs,
        )


def estimate_debye_temperature(v_s: float, n_density: float) -> float:
    r"""
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


def load_database(path: str | Path) -> list[Material]:
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

    for category, mats in db.get("materials", {}).items():
        for name, data in mats.items():
            sv = data.get("sound_velocities", {})
            v_avg = sv.get("average_m_per_s")
            density = data.get("density_kg_per_m3")

            if v_avg is None or density is None:
                continue

            # Build sound_velocity dict in MP format
            sound_velocity = {}
            if sv.get("longitudinal_m_per_s"):
                sound_velocity["longitudinal"] = sv["longitudinal_m_per_s"]
            if sv.get("transverse_m_per_s"):
                sound_velocity["transverse"] = sv["transverse_m_per_s"]
            if v_avg:
                sound_velocity["average"] = v_avg

            # Get number density from avg mass if available
            avg_mass = _DEFAULT_AVG_MASSES.get(name)
            if avg_mass is not None:
                n_density = estimate_n_density(density, avg_mass)
            else:
                n_density = estimate_n_density(density, 30.0)

            # Elastic moduli
            data.get("elastic_constants", {})
            em = data.get("elastic_moduli", {})
            bulk_modulus = None
            shear_modulus = None
            if em.get("bulk_modulus_GPa"):
                bulk_modulus = {"vrh": em["bulk_modulus_GPa"]}
            if em.get("shear_modulus_GPa"):
                shear_modulus = {"vrh": em["shear_modulus_GPa"]}

            # Sources
            sources = []
            for src in data.get("sources", []):
                if isinstance(src, dict):
                    sources.append(src.get("citation", ""))
                elif isinstance(src, str):
                    sources.append(src)

            mat = Material(
                name=name,
                formula=data.get("full_name", name),
                density=density,
                n_density=n_density,
                sound_velocity=sound_velocity,
                bulk_modulus=bulk_modulus,
                shear_modulus=shear_modulus,
                young_modulus=em.get("youngs_modulus_GPa"),
                category=category,
                sources=sources,
            )
            materials.append(mat)

    return materials


def fetch_mp(material_ids: list[str], api_key: str | None = None) -> list[Material]:
    """
    Fetch materials from the Materials Project API.

    Requires the ``mp-api`` package and a valid API key.

    Parameters
    ----------
    material_ids : list of str
        Materials Project IDs (e.g. ``["mp-5827", "mp-149"]``).
    api_key : str, optional
        MP API key. If None, uses the ``MP_API_KEY`` environment variable.

    Returns
    -------
    materials : list of Material
        Materials with elasticity data from MP.

    Examples
    --------
    >>> materials = fetch_mp(["mp-5827"])  # SrTiO3
    >>> print(materials[0].v_avg)
    """
    try:
        from mp_api.client import MPRester
    except ImportError as exc:
        raise ImportError(
            "mp-api is required for Materials Project queries. Install with: pip install mp-api"
        ) from exc

    kwargs = {}
    if api_key:
        kwargs["api_key"] = api_key

    materials = []
    with MPRester(**kwargs) as mpr:
        # Fetch elasticity data (includes sound velocity, moduli, θ_D)
        docs = mpr.materials.elasticity.search(
            material_ids=material_ids,
            fields=[
                "material_id",
                "formula_pretty",
                "density",
                "density_atomic",
                "volume",
                "nsites",
                "elements",
                "sound_velocity",
                "bulk_modulus",
                "shear_modulus",
                "youngs_modulus",
                "homogeneous_poisson",
                "universal_anisotropy",
                "debye_temperature",
                "thermal_conductivity",
                "symmetry",
            ],
        )
        for doc in docs:
            materials.append(Material.from_mp(doc))

    return materials


def batch_cahill(
    materials: list[Material], T: float | np.ndarray = 300.0
) -> dict[str, float | np.ndarray]:
    """
    Compute Cahill minimum κ for a list of materials.

    Parameters
    ----------
    materials : list of Material
        Materials with sound velocity and n_density set.
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
        if mat.v_avg is None or mat.debye_temperature is None:
            continue
        kappa = cahill.minimum_tc(T, mat.n_density, mat.debye_temperature, mat.v_avg)
        results[mat.name] = kappa

    return results
