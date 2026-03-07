#!/usr/bin/env python3
"""
κ_min vs κ_measured at 300 K — scatter plot with linear fit.

Recreates Figure 6 from Cheng et al., Small 17, 2101693 (2021),
extended with additional oxide/sulfide/halide solid electrolyte data.

κ_measured values are from the authoritative paper table (ver7_AIP_APR.tex)
to ensure figure-table consistency for peer review.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

from MatSciKit_COSMOTIM.thermal_conductivity import cahill
from MatSciKit_COSMOTIM.structure.material import estimate_debye_temperature, estimate_n_density


def kmin_300K(density, avg_mass, v_avg, theta_D=None, V=None, N=None):
    """Compute Cahill κ_min at 300 K."""
    if V is not None:
        n = N / V
        theta = theta_D if theta_D else estimate_debye_temperature(v_avg, n)
    else:
        n = estimate_n_density(density, avg_mass)
        theta = theta_D if theta_D else estimate_debye_temperature(v_avg, n)
    return cahill.minimum_tc(300.0, n, theta, v_avg)


# ═══════════════════════════════════════════════════════════════════════
# Material database — authoritative κ_measured values from paper table
# (ver7_AIP_APR.tex, Table I)
# ═══════════════════════════════════════════════════════════════════════

materials = []

def add(name, cat, kappa_meas, ref, density, avg_mass, v_avg,
        theta_D=None, **kw):
    km = kmin_300K(density, avg_mass, v_avg, theta_D, **kw)
    materials.append({
        'name': name, 'category': cat, 'kappa_meas': kappa_meas,
        'kappa_min': km, 'ref': ref,
    })


# ─── Oxides ─────────────────────────────────────────────────────────────

# LAGP — Rohde 2020
add('LAGP', 'oxide', 1.4, 'Rohde 2020',
    density=3090, avg_mass=28.95, v_avg=2430)

# NZP — Böger 2023
add('NZP', 'oxide', 1.13, 'Böger 2023',
    density=3800, avg_mass=24.46, v_avg=2900)

# NZSP — Wang 2025
add('NZSP', 'oxide', 1.0, 'Wang 2025',
    density=3400, avg_mass=28.0, v_avg=2700)

# NZS — Böger 2023
add('NZS', 'oxide', 1.0, 'Böger 2023',
    density=3300, avg_mass=26.0, v_avg=2600)

# LLZTO — Wang 2025 (SC, intrinsic)
add('LLZTO', 'oxide', 1.6, 'Wang 2025',
    density=5100, avg_mass=37.87, v_avg=2550)

# LSHT — Wang 2025 (our PRX Energy paper)
V_lsht = 3.98e-10 ** 3
N_lsht = 3/8 + 7/16 + 1/4 + 3/4 + 3
add('LSHT', 'oxide', 1.7, 'Wang 2025',
    density=None, avg_mass=None, v_avg=3461.3, theta_D=437.14,
    V=V_lsht, N=N_lsht)

# ─── Sulfides ───────────────────────────────────────────────────────────

# LGPS — Böger 2023
add('LGPS', 'sulfide', 0.8, 'Böger 2023',
    density=2900, avg_mass=24.32, v_avg=2000)

# Li₁₀SnP₂S₁₂ — Bron 2013
add('LSPS', 'sulfide', 0.7, 'Bron 2013',
    density=2700, avg_mass=25.0, v_avg=1900)

# LPSCl — Böger 2023
add('LPSCl', 'sulfide', 0.6, 'Böger 2023',
    density=2400, avg_mass=26.05, v_avg=2200)

# 3(Li₂S)-P₂S₅ — Cheng 2021
add('LPS', 'sulfide', 0.6, 'Cheng 2021',
    density=1900, avg_mass=20.0, v_avg=2000)

# Na₃PS₄ — Cheng/Bernges
add('Na₃PS₄', 'sulfide', 0.5, 'Cheng/Bernges',
    density=2500, avg_mass=28.78, v_avg=1900)

# NaSbS₂ — Gunatilleke 2023
add('NaSbS₂', 'sulfide', 0.6, 'Gunatilleke 2023',
    density=4200, avg_mass=38.43, v_avg=1500)

# ─── Halides ────────────────────────────────────────────────────────────

# Li₃InCl₆ — Cheng 2021
add('Li₃InCl₆', 'halide', 0.63, 'Cheng 2021',
    density=2700, avg_mass=33.62, v_avg=1800)

# Li₃YCl₆ — Cheng 2021
add('Li₃YCl₆', 'halide', 0.55, 'Cheng 2021',
    density=2600, avg_mass=32.26, v_avg=1800)


# ═══════════════════════════════════════════════════════════════════════
# Print results table
# ═══════════════════════════════════════════════════════════════════════

print(f"\n{'Material':15s} {'Category':10s} {'Reference':20s} {'κ_min':>8s} {'κ_meas':>8s} {'Ratio':>6s}")
print('-' * 72)
for m in sorted(materials, key=lambda x: x['kappa_meas']/x['kappa_min']):
    ratio = m['kappa_meas'] / m['kappa_min']
    print(f"{m['name']:15s} {m['category']:10s} {m['ref']:20s} "
          f"{m['kappa_min']:8.3f} {m['kappa_meas']:8.2f} {ratio:6.2f}")


# ═══════════════════════════════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(6, 6))

# Diagonal: κ_meas = κ_min
diag = np.linspace(0.1, 10, 100)
ax.plot(diag, diag, 'k-', linewidth=1, alpha=0.4, zorder=1)

# Style by category
cat_style = {
    'oxide':   {'color': '#d62728', 'marker': 'o', 'label': 'Oxide'},
    'sulfide': {'color': '#1f77b4', 'marker': 's', 'label': 'Sulfide'},
    'halide':  {'color': '#2ca02c', 'marker': '^', 'label': 'Halide'},
}

plotted_cats = set()
for m in materials:
    cs = cat_style[m['category']]

    label = cs['label'] if m['category'] not in plotted_cats else None
    plotted_cats.add(m['category'])

    ax.scatter(m['kappa_min'], m['kappa_meas'], zorder=3,
               c=cs['color'], marker=cs['marker'], edgecolors='black',
               s=100, linewidths=0.8, label=label, alpha=0.85)

    # Label positioning
    ha, va, ox, oy = 'left', 'bottom', 6, 4
    # Avoid overlaps
    if m['name'] in ('NZS', 'LPS'):
        va, oy = 'top', -8
    if m['name'] in ('Li₃YCl₆',):
        ha, ox = 'right', -6
    if m['name'] in ('LSPS',):
        va, oy = 'top', -8

    ax.annotate(m['name'], (m['kappa_min'], m['kappa_meas']),
                textcoords="offset points", xytext=(ox, oy),
                fontsize=7, ha=ha, va=va, color='#333')

# Linear fit in log space
log_kmin = np.log10([m['kappa_min'] for m in materials])
log_kmeas = np.log10([m['kappa_meas'] for m in materials])
slope, intercept, r, p, se = stats.linregress(log_kmin, log_kmeas)

x_fit = np.linspace(-0.6, 0.8, 100)
y_fit = slope * x_fit + intercept
ax.plot(10**x_fit, 10**y_fit, '--', color='#888', linewidth=1.5, alpha=0.7, zorder=2)

ax.text(0.55, 0.08,
        f'log-log fit: slope = {slope:.2f}, R² = {r**2:.2f}',
        transform=ax.transAxes, fontsize=9, color='#555',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ccc', alpha=0.8))

# Formatting
ax.set_xlabel(r'$\kappa_{\min}$ at 300 K (W m$^{-1}$ K$^{-1}$)', fontsize=12)
ax.set_ylabel(r'$\kappa_{\mathrm{measured}}$ at 300 K (W m$^{-1}$ K$^{-1}$)', fontsize=12)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.2, 3)
ax.set_ylim(0.2, 3)
ax.set_aspect('equal')
ax.tick_params(direction='in', which='both', top=True, right=True)
ax.legend(fontsize=10, loc='upper left', framealpha=0.9, edgecolor='#ccc')

# κ = κ_min annotation
ax.text(1.5, 1.3, r'$\kappa = \kappa_{\min}$',
        fontsize=9, color='gray', rotation=45, ha='center', va='center')

plt.tight_layout()

out = Path(__file__).parent / 'kappa_min_vs_kappa_300K'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
print(f"\nSaved to {out}.png")
print(f"Linear fit (log-log): slope={slope:.3f}, R²={r**2:.3f}, p={p:.1e}")
plt.close(fig)
