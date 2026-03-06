#!/usr/bin/env python3
"""
Plot measured thermal conductivity vs Cahill minimum for solid electrolytes.

Generates a multi-panel or overlay comparison of experimental κ(T) data
(PPMS TTO and LFA) against the Cahill minimum thermal conductivity model
for oxide solid electrolytes.

Data sources:
  - PPMS TTO measurements from Yitian's lab
  - LFA measurements from Sally Jia
  - Cahill model computed using MatSciKit

For the SE thermal review paper.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from MatSciKit_COSMOTIM.thermal_conductivity import cahill, porosity_correction
from MatSciKit_COSMOTIM.structure.material import Material, estimate_debye_temperature, estimate_n_density

# ─── Data paths ─────────────────────────────────────────────────────────
DATA_DIR = Path("/Users/cosmotim/Documents/SE thermal review/TC of SE/TTO Plot/Export_Data")

# ─── Material definitions ───────────────────────────────────────────────
# Each material: name, formula, density (kg/m³), avg_mass (amu),
#                v_avg (m/s), porosity, PPMS file, LFA file, color, marker

MATERIALS = {
    'LSHT': {
        'formula': r'Li$_{3/8}$Sr$_{7/16}$Hf$_{1/4}$Ta$_{3/4}$O$_3$',
        'density': None,  # computed from structure
        'V': 3.98e-10 ** 3,  # m³
        'N': 3/8 + 7/16 + 1/4 + 3/4 + 3,
        'M_amu': 6.941*3/8 + 87.62*7/16 + 178.49*1/4 + 180.95*3/4 + 15.999*3,
        'v_avg': 3461.3,  # from Pipeline 1 fit
        'theta_D': 437.14,  # from Pipeline 1 fit
        'porosity': 0.03,
        'ppms_file': 'LSHT_PPMS.csv',
        'lfa_file': 'LSHT_LFA.csv',
        'color': '#d62728',  # red
        'marker': 'o',
    },
    'LAGP': {
        'formula': r'Li$_{1.5}$Al$_{0.5}$Ge$_{1.5}$(PO$_4$)$_3$',
        'density': 3090,  # kg/m³
        'avg_mass_amu': (6.941*1.5 + 26.982*0.5 + 72.630*1.5 + 30.974*3 + 15.999*12) / (1.5+0.5+1.5+3+12),
        'v_avg': 2430,
        'porosity': 0.0,
        'ppms_file': 'LAGP_PPMS.csv',
        'lfa_file': 'LAGP_LFA.csv',
        'color': '#1f77b4',  # blue
        'marker': 's',
    },
    'NZP': {
        'formula': r'NaZr$_2$(PO$_4$)$_3$',
        'density': 3800,  # kg/m³
        'avg_mass_amu': (22.990 + 91.224*2 + 30.974*3 + 15.999*12) / (1+2+3+12),
        'v_avg': 2900,
        'porosity': 0.0,
        'ppms_file': 'NZP_PPMS.csv',
        'lfa_file': 'NZP_LFA.csv',
        'color': '#2ca02c',  # green
        'marker': '^',
    },
    'LLZTO_poly': {
        'formula': r'Li$_{6.5}$La$_3$Zr$_{1.5}$Ta$_{0.5}$O$_{12}$ (poly)',
        'density': 5100,  # kg/m³
        'avg_mass_amu': (6.941*6.5 + 138.91*3 + 91.224*1.5 + 180.95*0.5 + 15.999*12) / (6.5+3+1.5+0.5+12),
        'v_avg': 2550,
        'porosity': 0.0,
        'ppms_file': 'LLZTO poly_PPMS.csv',
        'lfa_file': 'LLZTO poly_LFA.csv',
        'color': '#ff7f0e',  # orange
        'marker': 'D',
    },
    'LLZTO_SC': {
        'formula': r'Li$_{6.5}$La$_3$Zr$_{1.5}$Ta$_{0.5}$O$_{12}$ (SC)',
        'density': 5100,  # kg/m³
        'avg_mass_amu': (6.941*6.5 + 138.91*3 + 91.224*1.5 + 180.95*0.5 + 15.999*12) / (6.5+3+1.5+0.5+12),
        'v_avg': 2550,
        'porosity': 0.0,
        'ppms_file': 'LLZTO SC_PPMS.csv',
        'lfa_file': 'LLZTO SC_LFA.csv',
        'color': '#9467bd',  # purple
        'marker': 'v',
    },
}


def load_ppms(path: Path) -> tuple:
    """Load PPMS TTO data: [T, κ, κ_err]."""
    data = np.loadtxt(str(path), delimiter=',')
    T = data[:, 0]
    kappa = data[:, 1]
    kappa_err = data[:, 2] if data.shape[1] > 2 else None
    # Sort by temperature
    idx = np.argsort(T)
    T, kappa = T[idx], kappa[idx]
    if kappa_err is not None:
        kappa_err = kappa_err[idx]
    return T, kappa, kappa_err


def load_lfa(path: Path) -> tuple:
    """Load LFA data: [T, κ]."""
    data = np.loadtxt(str(path), delimiter=',')
    T = data[:, 0]
    kappa = data[:, 1]
    idx = np.argsort(T)
    return T[idx], kappa[idx]


def get_n_density_theta(mat_info: dict) -> tuple:
    """Compute n_density and theta_D for a material."""
    if 'theta_D' in mat_info and mat_info['theta_D']:
        # Use fitted values (e.g. LSHT)
        n_density = mat_info['N'] / mat_info['V']
        return n_density, mat_info['theta_D']
    else:
        # Estimate from v_avg and density
        n_density = estimate_n_density(mat_info['density'], mat_info['avg_mass_amu'])
        theta_D = estimate_debye_temperature(mat_info['v_avg'], n_density)
        return n_density, theta_D


def main():
    # Temperature range for Cahill curves
    T_theory = np.linspace(1, 1000, 500)

    # ─── Figure: κ vs T with Cahill minimum ─────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5.5))

    for name, info in MATERIALS.items():
        # Compute Cahill minimum
        n_density, theta_D = get_n_density_theta(info)
        kappa_min = cahill.minimum_tc(T_theory, n_density, theta_D, info['v_avg'])

        # Plot Cahill minimum curve
        ax.plot(T_theory, kappa_min, '--', color=info['color'], alpha=0.6, linewidth=1.2)

        # Load and plot PPMS data
        ppms_path = DATA_DIR / info['ppms_file']
        if ppms_path.exists():
            T_ppms, kappa_ppms, kappa_err = load_ppms(ppms_path)

            # Apply porosity correction if needed
            if info.get('porosity', 0) > 0:
                kappa_ppms = porosity_correction.correct(kappa_ppms, info['porosity'])

            ax.scatter(T_ppms, kappa_ppms, c=info['color'], marker=info['marker'],
                       s=15, alpha=0.7, edgecolors='none',
                       label=f"{name} (PPMS)")
        else:
            print(f"  Warning: PPMS file not found: {ppms_path}")

        # Load and plot LFA data
        lfa_path = DATA_DIR / info['lfa_file']
        if lfa_path.exists():
            T_lfa, kappa_lfa = load_lfa(lfa_path)

            if info.get('porosity', 0) > 0:
                kappa_lfa = porosity_correction.correct(kappa_lfa, info['porosity'])

            ax.scatter(T_lfa, kappa_lfa, c=info['color'], marker=info['marker'],
                       s=40, alpha=0.9, edgecolors='black', linewidths=0.5,
                       label=f"{name} (LFA)")
        else:
            print(f"  Warning: LFA file not found: {lfa_path}")

        print(f"{name}: θ_D={theta_D:.0f} K, v_avg={info['v_avg']:.0f} m/s, "
              f"κ_min(300K)={cahill.minimum_tc(300.0, n_density, theta_D, info['v_avg']):.3f} W/(m·K)")

    # Formatting
    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel(r'$\kappa$ (W m$^{-1}$ K$^{-1}$)', fontsize=12)
    ax.set_xlim(0, 700)
    ax.set_ylim(bottom=0)
    ax.tick_params(direction='in', which='both', top=True, right=True)
    ax.legend(fontsize=8, ncol=2, loc='upper left', framealpha=0.8)
    ax.set_title('Thermal Conductivity vs Cahill Minimum — Oxide Solid Electrolytes', fontsize=11)

    plt.tight_layout()

    out_path = Path(__file__).parent / 'kappa_min_vs_measured'
    fig.savefig(str(out_path) + '.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved to {out_path}.png")

    # Also save a log-scale version for low-T behavior
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1, 1000)
    ax.set_ylim(0.01, 10)
    ax.set_title('Thermal Conductivity vs Cahill Minimum (log-log)', fontsize=11)
    fig.savefig(str(out_path) + '_loglog.png', dpi=300, bbox_inches='tight')
    print(f"Saved to {out_path}_loglog.png")

    plt.close(fig)


if __name__ == '__main__':
    main()
