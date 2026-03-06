#!/usr/bin/env python3
"""
Cahill Minimum Thermal Conductivity Comparison for Solid Electrolytes.

Computes the Cahill minimum κ for multiple solid electrolyte materials
using the sound velocity database compiled by the researcher agent.

This script demonstrates the full MatSciKit workflow:
  CIF / manual params → Pipeline 1 (v_s) → Pipeline 2 (κ_min)

Usage:
    python cahill_comparison.py [--database path/to/database.json]
    python cahill_comparison.py --cif path/to/file.cif --v_avg 3000

References:
    Cahill, D. G., Watson, S. K., & Pohl, R. O.
    Lower limit to the thermal conductivity of disordered crystals.
    Phys. Rev. B 46, 6131 (1992).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Add src to path if running from examples/
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from MatSciKit_COSMOTIM.thermal_conductivity import cahill
from MatSciKit_COSMOTIM.constants import kb


# Default database path (researcher agent workspace)
DEFAULT_DB = Path(__file__).parent.parent.parent / (
    'MatSciKit_Matlab_Prototype'  # fallback
)
RESEARCHER_DB = Path.home() / '.openclaw/agents/researcher/sound_velocity_database.json'


def load_database(db_path: str) -> Dict:
    """Load the sound velocity database JSON."""
    with open(db_path) as f:
        return json.load(f)


def extract_materials(db: Dict) -> List[Dict]:
    """
    Extract materials with complete data from the database.

    Returns a list of dicts with keys:
        name, full_name, category, v_avg, density, n_density_approx
    """
    materials = []

    for category, mats in db.get('materials', {}).items():
        for name, data in mats.items():
            sv = data.get('sound_velocities', {})
            v_avg = sv.get('average_m_per_s')
            density = data.get('density_kg_per_m3')

            if v_avg is None or density is None:
                continue

            materials.append({
                'name': name,
                'full_name': data.get('full_name', name),
                'category': category,
                'v_avg': v_avg,
                'v_L': sv.get('longitudinal_m_per_s'),
                'v_T': sv.get('transverse_m_per_s'),
                'density': density,
                'search_status': data.get('search_status', 'unknown'),
            })

    return materials


def compute_cahill_for_material(
    T: np.ndarray,
    n_density: float,
    theta_D: float,
    v_s: float,
) -> np.ndarray:
    """
    Compute Cahill minimum κ(T) for a material.

    Uses MatSciKit's cahill.minimum_tc which implements the full
    Debye integral model.

    Parameters
    ----------
    T : np.ndarray
        Temperature array (K).
    n_density : float
        Number density (atoms/m³).
    theta_D : float
        Debye temperature (K).
    v_s : float
        Average sound velocity (m/s).

    Returns
    -------
    kappa_min : np.ndarray
        Minimum thermal conductivity (W/(m·K)).
    """
    return cahill.minimum_tc(T, n_density, theta_D, v_s)


def estimate_theta_D(v_s: float, n_density: float) -> float:
    """
    Estimate Debye temperature from sound velocity and number density.

    θ_D = (ħ/k_B) * v_s * (6π² n)^(1/3)

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
    from MatSciKit_COSMOTIM.constants import hbar

    theta_D = (hbar / kb) * v_s * (6 * np.pi**2 * n_density) ** (1/3)
    return theta_D


def estimate_n_density(density: float, avg_mass_amu: float) -> float:
    """
    Estimate number density from mass density and average atomic mass.

    n = ρ / (M_avg * AMU_TO_KG)
    """
    from MatSciKit_COSMOTIM.constants import AMU_TO_KG
    return density / (avg_mass_amu * AMU_TO_KG)


# Average atomic masses (amu) for each material
# Calculated from composition
AVG_MASSES = {
    'LAGP': (6.941*1.5 + 26.982*0.5 + 72.630*1.5 + 30.974*3 + 15.999*12) / (1.5+0.5+1.5+3+12),
    'NZP': (22.990 + 91.224*2 + 30.974*3 + 15.999*12) / (1+2+3+12),
    'LLZTO': (6.941*6.5 + 138.91*3 + 91.224*1.5 + 180.95*0.5 + 15.999*12) / (6.5+3+1.5+0.5+12),
    'LGPS': (6.941*10 + 72.630 + 30.974*2 + 32.06*12) / (10+1+2+12),
    'LPSCl': (6.941*6 + 30.974 + 32.06*5 + 35.45) / (6+1+5+1),
    'Na3PS4': (22.990*3 + 30.974 + 32.06*4) / (3+1+4),
    'Li3InCl6': (6.941*3 + 114.82 + 35.45*6) / (3+1+6),
}


def main():
    parser = argparse.ArgumentParser(
        description='Compute Cahill minimum κ for solid electrolytes'
    )
    parser.add_argument(
        '--database', '-d', type=str,
        default=str(RESEARCHER_DB),
        help='Path to sound_velocity_database.json'
    )
    parser.add_argument(
        '--T_min', type=float, default=1.0,
        help='Minimum temperature (K)'
    )
    parser.add_argument(
        '--T_max', type=float, default=1000.0,
        help='Maximum temperature (K)'
    )
    parser.add_argument(
        '--T_eval', type=float, default=300.0,
        help='Temperature for comparison table (K)'
    )
    parser.add_argument(
        '--output', '-o', type=str,
        help='Output CSV file'
    )
    parser.add_argument(
        '--plot', action='store_true',
        help='Generate comparison plot'
    )

    args = parser.parse_args()

    # Load database
    db_path = Path(args.database)
    if not db_path.exists():
        print(f"Database not found: {db_path}", file=sys.stderr)
        print("Using built-in material data.", file=sys.stderr)
        db = None
    else:
        db = load_database(str(db_path))

    # Extract materials with data
    if db:
        materials = extract_materials(db)
    else:
        materials = []

    if not materials:
        print("No materials with complete data found in database.")
        return

    # Temperature range
    T = np.linspace(args.T_min, args.T_max, 500)

    # Compute Cahill minimum κ for each material
    print(f"\n{'='*75}")
    print(f"  Cahill Minimum Thermal Conductivity at T = {args.T_eval:.0f} K")
    print(f"{'='*75}")
    print(f"{'Material':<12} {'Formula':<25} {'v_avg':>8} {'ρ':>8} {'θ_D':>8} {'κ_min':>10}")
    print(f"{'':12} {'':25} {'(m/s)':>8} {'(kg/m³)':>8} {'(K)':>8} {'(W/m·K)':>10}")
    print(f"{'-'*75}")

    results = []

    for mat in materials:
        name = mat['name']
        v_avg = mat['v_avg']
        density = mat['density']

        # Get average atomic mass
        avg_mass = AVG_MASSES.get(name)
        if avg_mass is None:
            print(f"  {name:<12} — skipped (no avg mass data)")
            continue

        # Estimate number density and Debye temperature
        n_density = estimate_n_density(density, avg_mass)
        theta_D = estimate_theta_D(v_avg, n_density)

        # Compute κ_min at evaluation temperature
        kappa_min = cahill.minimum_tc(args.T_eval, n_density, theta_D, v_avg)

        # Compute full curve
        kappa_curve = cahill.minimum_tc(T, n_density, theta_D, v_avg)

        print(f"  {name:<12} {mat['full_name']:<25} {v_avg:>8.0f} {density:>8.0f} "
              f"{theta_D:>8.0f} {kappa_min:>10.3f}")

        results.append({
            'name': name,
            'full_name': mat['full_name'],
            'category': mat['category'],
            'v_avg': v_avg,
            'density': density,
            'n_density': n_density,
            'theta_D': theta_D,
            'kappa_min': kappa_min,
            'kappa_curve': kappa_curve,
        })

    print(f"{'-'*75}")
    print(f"  {len(results)} materials computed")
    print()

    # Save CSV output
    if args.output:
        header = 'Material,Formula,Category,v_avg(m/s),density(kg/m3),theta_D(K),kappa_min(W/mK)\n'
        with open(args.output, 'w') as f:
            f.write(header)
            for r in results:
                f.write(f"{r['name']},{r['full_name']},{r['category']},"
                        f"{r['v_avg']:.0f},{r['density']:.0f},"
                        f"{r['theta_D']:.1f},{r['kappa_min']:.4f}\n")
        print(f"Results saved to {args.output}")

    # Generate plot
    if args.plot:
        try:
            import matplotlib.pyplot as plt
            from MatSciKit_COSMOTIM.visualization.journal_style import (
                export_journal_figure,
            )

            fig, ax = plt.subplots(figsize=(8, 5))

            # Color scheme by category
            colors = {
                'oxides': '#1f77b4',
                'sulfides': '#ff7f0e',
                'halides': '#2ca02c',
                'others': '#d62728',
            }
            markers = {
                'oxides': 'o',
                'sulfides': 's',
                'halides': '^',
                'others': 'D',
            }

            for r in results:
                color = colors.get(r['category'], 'gray')
                marker = markers.get(r['category'], 'o')
                ax.plot(T, r['kappa_curve'], '-', color=color, alpha=0.7)
                # Mark at T_eval
                ax.plot(args.T_eval, r['kappa_min'], marker,
                        color=color, markersize=8,
                        label=f"{r['name']} ({r['kappa_min']:.2f})")

            ax.set_xlabel('Temperature (K)')
            ax.set_ylabel(r'$\kappa_{\min}$ (W m$^{-1}$ K$^{-1}$)')
            ax.set_title('Cahill Minimum Thermal Conductivity — Solid Electrolytes')
            ax.legend(fontsize=8, ncol=2)
            ax.set_xlim(0, args.T_max)
            ax.set_ylim(bottom=0)

            out_path = args.output.replace('.csv', '') if args.output else 'cahill_comparison'
            export_journal_figure(fig, out_path, format='png')
            print(f"Plot saved to {out_path}.png")

        except ImportError as e:
            print(f"Plotting skipped: {e}", file=sys.stderr)


if __name__ == '__main__':
    main()
