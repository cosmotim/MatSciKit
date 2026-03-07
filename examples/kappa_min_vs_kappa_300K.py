#!/usr/bin/env python3
"""
κ_min vs κ_measured at 300 K — scatter plot with linear fit.

Recreates Figure 6 from Cheng et al., Small 17, 2101693 (2021),
extended with additional oxide/sulfide/halide solid electrolyte data.
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

DATA = Path("/Users/cosmotim/Documents/SE thermal review")
EXPORT = DATA / "TC of SE/TTO Plot/Export_Data"
CHENG = DATA / "Data"


def interp_300K(T, kappa):
    """Interpolate κ at 300 K."""
    idx = np.argsort(T)
    T, kappa = T[idx], kappa[idx]
    return float(np.interp(300, T, kappa))


def kmin_300K(density, avg_mass, v_avg, theta_D=None, V=None, N=None, M_amu=None):
    """Compute Cahill κ_min at 300 K."""
    if V is not None:
        n = N / V
        theta = theta_D if theta_D else estimate_debye_temperature(v_avg, n)
    else:
        n = estimate_n_density(density, avg_mass)
        theta = theta_D if theta_D else estimate_debye_temperature(v_avg, n)
    return cahill.minimum_tc(300.0, n, theta, v_avg)


def load_kappa_300(filepath, col_T=0, col_k=1):
    """Load CSV and interpolate at 300 K."""
    data = np.loadtxt(str(filepath), delimiter=',')
    return interp_300K(data[:, col_T], data[:, col_k])


# ═══════════════════════════════════════════════════════════════════════
# Material database: name, category, source, κ_measured(300K), κ_min params
# ═══════════════════════════════════════════════════════════════════════

materials = []

def add(name, cat, source, kappa_meas, density, avg_mass, v_avg, theta_D=None, **kw):
    km = kmin_300K(density, avg_mass, v_avg, theta_D, **kw)
    materials.append({'name': name, 'category': cat, 'source': source,
                      'kappa_meas': kappa_meas, 'kappa_min': km})


# ─── Our data (PPMS/LFA) ───────────────────────────────────────────────

# LSHT (from Pipeline 1 fit)
V_lsht = 3.98e-10 ** 3
N_lsht = 3/8 + 7/16 + 1/4 + 3/4 + 3
M_lsht = 6.941*3/8 + 87.62*7/16 + 178.49*1/4 + 180.95*3/4 + 15.999*3
add('LSHT', 'oxide', 'This work',
    load_kappa_300(EXPORT / 'LSHT_LFA.csv'),
    None, None, 3461.3, 437.14, V=V_lsht, N=N_lsht)

# LAGP (our LFA)
add('LAGP', 'oxide', 'This work',
    load_kappa_300(EXPORT / 'LAGP_LFA.csv'),
    3090, 28.95, 2430)

# NZP (our LFA)
add('NZP', 'oxide', 'This work',
    load_kappa_300(EXPORT / 'NZP_LFA.csv'),
    3800, 24.46, 2900)

# LLZTO poly (our LFA)
add('LLZTO (poly)', 'oxide', 'This work',
    load_kappa_300(EXPORT / 'LLZTO poly_LFA.csv'),
    5100, 37.87, 2550)

# LLZTO SC (our LFA)
add('LLZTO (SC)', 'oxide', 'This work',
    load_kappa_300(EXPORT / 'LLZTO SC_LFA.csv'),
    5100, 37.87, 2550)

# ─── Cheng et al. (TDTR) ───────────────────────────────────────────────

add('LAGP*', 'oxide', 'Cheng 2021',
    load_kappa_300(CHENG / 'LAGP_Cheng.csv'),
    3090, 28.95, 2430)

add('LLZTO*', 'oxide', 'Cheng 2021',
    load_kappa_300(CHENG / 'LLZTO_ZheCheng.csv'),
    5100, 37.87, 2550)

add('LPSC', 'sulfide', 'Cheng 2021',
    load_kappa_300(CHENG / 'LPSC_ZheCheng.csv'),
    2400, 26.05, 2200)

add('LPS', 'sulfide', 'Cheng 2021',
    load_kappa_300(CHENG / 'LPS_ZheCheng.csv'),
    1900, 20.0, 2000)

add('NPS*', 'sulfide', 'Cheng 2021',
    load_kappa_300(CHENG / 'NPS_ZheCheng.csv'),
    2500, 28.78, 1900)

add('LIC', 'halide', 'Cheng 2021',
    load_kappa_300(CHENG / 'LIC.csv'),
    2700, 33.62, 1800)

add('LYC', 'halide', 'Cheng 2021',
    load_kappa_300(CHENG / 'LYC.csv'),
    2600, 32.26, 1800)

# ─── Other literature data ─────────────────────────────────────────────

# LGPS
add('LGPS', 'sulfide', 'Literature',
    load_kappa_300(CHENG / 'LGPS.csv'),
    2900, 24.32, 2000)

# LPSC (our PPMS)
add('LPSC (PPMS)', 'sulfide', 'Literature',
    load_kappa_300(CHENG / 'LPSC.csv'),
    2400, 26.05, 2200)

# NPS (Bernges)
add('NPS (Bernges)', 'sulfide', 'Literature',
    load_kappa_300(CHENG / 'NPS_Bernges.csv'),
    2500, 28.78, 1900)

# NZSP
add('NZSP', 'oxide', 'Literature',
    load_kappa_300(CHENG / 'NZSP.csv'),
    3400, 28.0, 2700)

# NATP (NaTi₂(PO₄)₃)
add('NATP', 'oxide', 'Literature',
    load_kappa_300(CHENG / 'NATP.csv'),
    3200, 24.0, 2800)

# Li₂TiO₃
add('Li₂TiO₃', 'oxide', 'Literature',
    load_kappa_300(CHENG / 'Li2TiO3.csv'),
    3430, 16.35, 4000)

# NaSbS₂
add('NaSbS₂', 'sulfide', 'Literature',
    load_kappa_300(CHENG / 'NaSbS2.csv'),
    4200, 38.43, 1500)

# NZP (Boger)
add('NZP (Böger)', 'oxide', 'Literature',
    load_kappa_300(CHENG / 'NZP_Boger.csv'),
    3800, 24.46, 2900)

# NZS (Boger) — Na₄Zr₂(SiO₄)₃
add('NZS', 'oxide', 'Literature',
    load_kappa_300(CHENG / 'NZS_Boger.csv'),
    3300, 26.0, 2600)


# ═══════════════════════════════════════════════════════════════════════
# Print results table
# ═══════════════════════════════════════════════════════════════════════

print(f"\n{'Material':20s} {'Category':10s} {'Source':15s} {'κ_min':>8s} {'κ_meas':>8s} {'Ratio':>6s}")
print('-' * 72)
for m in sorted(materials, key=lambda x: x['kappa_meas']/x['kappa_min']):
    ratio = m['kappa_meas'] / m['kappa_min']
    print(f"{m['name']:20s} {m['category']:10s} {m['source']:15s} "
          f"{m['kappa_min']:8.3f} {m['kappa_meas']:8.3f} {ratio:6.2f}")


# ═══════════════════════════════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 7))

# Diagonal: κ_meas = κ_min
diag = np.linspace(0.1, 10, 100)
ax.plot(diag, diag, 'k-', linewidth=1, alpha=0.4, zorder=1)

# Style by category
cat_style = {
    'oxide':   {'color': '#d62728', 'marker': 'o', 'label': 'Oxide'},
    'sulfide': {'color': '#1f77b4', 'marker': 's', 'label': 'Sulfide'},
    'halide':  {'color': '#2ca02c', 'marker': '^', 'label': 'Halide'},
}

# Marker size by source
size_map = {'This work': 120, 'Cheng 2021': 90, 'Literature': 70}
edge_map = {'This work': 'black', 'Cheng 2021': '#333', 'Literature': '#666'}

# Plot each point
plotted_cats = set()
for m in materials:
    cs = cat_style[m['category']]
    sz = size_map.get(m['source'], 70)
    ec = edge_map.get(m['source'], 'gray')
    lw = 1.2 if m['source'] == 'This work' else 0.8

    label = cs['label'] if m['category'] not in plotted_cats else None
    plotted_cats.add(m['category'])

    ax.scatter(m['kappa_min'], m['kappa_meas'], zorder=3,
               c=cs['color'], marker=cs['marker'], edgecolors=ec,
               s=sz, linewidths=lw, label=label, alpha=0.85)

    # Label
    ha, va, ox, oy = 'left', 'bottom', 6, 4
    if m['name'] in ('LLZTO (SC)', 'LYC', 'NPS (Bernges)', 'LPS'):
        va, oy = 'top', -6
    if m['name'] in ('LAGP*', 'LLZTO*'):
        ha, ox = 'right', -6

    ax.annotate(m['name'], (m['kappa_min'], m['kappa_meas']),
                textcoords="offset points", xytext=(ox, oy),
                fontsize=6.5, ha=ha, va=va, color='#333')

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
ax.set_xlim(0.15, 5)
ax.set_ylim(0.15, 5)
ax.set_aspect('equal')
ax.tick_params(direction='in', which='both', top=True, right=True)
ax.legend(fontsize=10, loc='upper left', framealpha=0.9,
          edgecolor='#ccc', markerscale=1.2)

# κ = κ_min annotation
ax.text(2.5, 2.1, r'$\kappa = \kappa_{\min}$',
        fontsize=9, color='gray', rotation=45, ha='center', va='center')

plt.tight_layout()

out = Path(__file__).parent / 'kappa_min_vs_kappa_300K'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
print(f"\nSaved to {out}.png")
print(f"Linear fit (log-log): slope={slope:.3f}, R²={r**2:.3f}, p={p:.1e}")
plt.close(fig)
