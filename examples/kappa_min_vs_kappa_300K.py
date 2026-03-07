#!/usr/bin/env python3
"""
κ_min vs κ_measured at 300 K — Cheng et al. (2021) style scatter plot.

Recreates Figure 6 from:
  Cheng et al., "Good Solid-State Electrolytes Have Low, Glass-Like
  Thermal Conductivity," Small 17, 2101693 (2021).

Combines our experimental data (LSHT, LAGP, NZP, LLZTO) with
Cheng et al.'s TDTR data (LPSC, LPS, NPS, LIC, LYC).
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from MatSciKit_COSMOTIM.thermal_conductivity import cahill
from MatSciKit_COSMOTIM.structure.material import estimate_debye_temperature, estimate_n_density

# ─── Data directories ───────────────────────────────────────────────────
OUR_DATA = Path("/Users/cosmotim/Documents/SE thermal review/TC of SE/TTO Plot/Export_Data")
CHENG_DATA = Path("/Users/cosmotim/Documents/SE thermal review/Data")


def interp_at_300K(T, kappa):
    """Interpolate κ at 300 K from measured data."""
    idx = np.argsort(T)
    T, kappa = T[idx], kappa[idx]
    # Find points bracketing 300 K
    mask = (T > 250) & (T < 350)
    if np.sum(mask) > 0:
        return float(np.interp(300, T[mask], kappa[mask]))
    return float(np.interp(300, T, kappa))


def compute_kappa_min_300K(density, avg_mass_amu, v_avg, theta_D=None):
    """Compute Cahill κ_min at 300 K."""
    n = estimate_n_density(density, avg_mass_amu)
    if theta_D is None:
        theta_D = estimate_debye_temperature(v_avg, n)
    return cahill.minimum_tc(300.0, n, theta_D, v_avg)


# ─── Our materials (measured PPMS + LFA) ────────────────────────────────
# Get κ_measured at 300 K from experimental data
our_materials = {}

# LSHT — use PPMS + LFA
lsht_ppms = np.loadtxt(str(OUR_DATA / 'LSHT_PPMS.csv'), delimiter=',')
lsht_lfa = np.loadtxt(str(OUR_DATA / 'LSHT_LFA.csv'), delimiter=',')
# LFA is closer to 300K
our_materials['LSHT'] = {
    'kappa_300': interp_at_300K(lsht_lfa[:, 0], lsht_lfa[:, 1]),
    'density': None,  # compute from structure
    'V': 3.98e-10 ** 3,
    'N': 3/8 + 7/16 + 1/4 + 3/4 + 3,
    'M_amu': 6.941*3/8 + 87.62*7/16 + 178.49*1/4 + 180.95*3/4 + 15.999*3,
    'v_avg': 3461.3,
    'theta_D': 437.14,
    'category': 'oxide',
    'source': 'This work',
}

# LAGP
lagp = np.loadtxt(str(OUR_DATA / 'LAGP_LFA.csv'), delimiter=',')
our_materials['LAGP'] = {
    'kappa_300': interp_at_300K(lagp[:, 0], lagp[:, 1]),
    'density': 3090, 'avg_mass_amu': 28.95, 'v_avg': 2430,
    'category': 'oxide', 'source': 'This work',
}

# NZP
nzp = np.loadtxt(str(OUR_DATA / 'NZP_LFA.csv'), delimiter=',')
our_materials['NZP'] = {
    'kappa_300': interp_at_300K(nzp[:, 0], nzp[:, 1]),
    'density': 3800, 'avg_mass_amu': 24.46, 'v_avg': 2900,
    'category': 'oxide', 'source': 'This work',
}

# LLZTO poly
llzto_p = np.loadtxt(str(OUR_DATA / 'LLZTO poly_LFA.csv'), delimiter=',')
our_materials['LLZTO\n(poly)'] = {
    'kappa_300': interp_at_300K(llzto_p[:, 0], llzto_p[:, 1]),
    'density': 5100, 'avg_mass_amu': 37.87, 'v_avg': 2550,
    'category': 'oxide', 'source': 'This work',
}

# LLZTO SC
llzto_sc = np.loadtxt(str(OUR_DATA / 'LLZTO SC_LFA.csv'), delimiter=',')
our_materials['LLZTO\n(SC)'] = {
    'kappa_300': interp_at_300K(llzto_sc[:, 0], llzto_sc[:, 1]),
    'density': 5100, 'avg_mass_amu': 37.87, 'v_avg': 2550,
    'category': 'oxide', 'source': 'This work',
}

# ─── Cheng et al. data (TDTR measurements) ──────────────────────────────
cheng_materials = {}

# LAGP (Cheng)
lagp_c = np.loadtxt(str(CHENG_DATA / 'LAGP_Cheng.csv'), delimiter=',')
cheng_materials['LAGP\n(Cheng)'] = {
    'kappa_300': interp_at_300K(lagp_c[:, 0], lagp_c[:, 1]),
    'density': 3090, 'avg_mass_amu': 28.95, 'v_avg': 2430,
    'category': 'oxide', 'source': 'Cheng 2021',
}

# LLZTO (Cheng)
llzto_c = np.loadtxt(str(CHENG_DATA / 'LLZTO_ZheCheng.csv'), delimiter=',')
cheng_materials['LLZTO\n(Cheng)'] = {
    'kappa_300': interp_at_300K(llzto_c[:, 0], llzto_c[:, 1]),
    'density': 5100, 'avg_mass_amu': 37.87, 'v_avg': 2550,
    'category': 'oxide', 'source': 'Cheng 2021',
}

# LPSC (Cheng)
lpsc = np.loadtxt(str(CHENG_DATA / 'LPSC_ZheCheng.csv'), delimiter=',')
cheng_materials['LPSC'] = {
    'kappa_300': interp_at_300K(lpsc[:, 0], lpsc[:, 1]),
    'density': 2400, 'avg_mass_amu': 26.05, 'v_avg': 2200,
    'category': 'sulfide', 'source': 'Cheng 2021',
}

# LPS (Cheng) — amorphous Li₂S-P₂S₅
lps = np.loadtxt(str(CHENG_DATA / 'LPS_ZheCheng.csv'), delimiter=',')
cheng_materials['LPS'] = {
    'kappa_300': interp_at_300K(lps[:, 0], lps[:, 1]),
    'density': 1900, 'avg_mass_amu': 20.0, 'v_avg': 2000,
    'category': 'sulfide', 'source': 'Cheng 2021',
}

# NPS (Cheng) — Na₃PS₄
nps = np.loadtxt(str(CHENG_DATA / 'NPS_ZheCheng.csv'), delimiter=',')
cheng_materials['NPS'] = {
    'kappa_300': interp_at_300K(nps[:, 0], nps[:, 1]),
    'density': 2500, 'avg_mass_amu': 28.78, 'v_avg': 1900,
    'category': 'sulfide', 'source': 'Cheng 2021',
}


# ─── Compute κ_min for all ──────────────────────────────────────────────
all_data = []

for name, info in {**our_materials, **cheng_materials}.items():
    if 'V' in info:
        # LSHT — compute from unit cell
        n = info['N'] / info['V']
        density = info['M_amu'] * 1.66053906660e-27 / info['V']
        kmin = cahill.minimum_tc(300.0, n, info['theta_D'], info['v_avg'])
    else:
        n = estimate_n_density(info['density'], info['avg_mass_amu'])
        theta_D = estimate_debye_temperature(info['v_avg'], n)
        kmin = cahill.minimum_tc(300.0, n, theta_D, info['v_avg'])

    all_data.append({
        'name': name,
        'kappa_min': kmin,
        'kappa_meas': info['kappa_300'],
        'category': info['category'],
        'source': info['source'],
    })
    print(f"{name.replace(chr(10),' '):20s}  κ_min={kmin:.3f}  κ_meas={info['kappa_300']:.3f}  "
          f"ratio={info['kappa_300']/kmin:.2f}")


# ─── Plot ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))

# Diagonal: κ_meas = κ_min
diag = np.linspace(0.1, 10, 100)
ax.plot(diag, diag, 'k-', linewidth=1, alpha=0.5, zorder=1)
ax.fill_between(diag, diag, 10, alpha=0.04, color='gray')

# Color/marker by category and source
style = {
    ('oxide', 'This work'):    {'color': '#d62728', 'marker': 'o', 'edgecolor': 'black', 's': 100},
    ('oxide', 'Cheng 2021'):   {'color': '#ff9896', 'marker': 's', 'edgecolor': 'black', 's': 80},
    ('sulfide', 'Cheng 2021'): {'color': '#aec7e8', 'marker': '^', 'edgecolor': 'black', 's': 80},
    ('halide', 'Cheng 2021'):  {'color': '#98df8a', 'marker': 'D', 'edgecolor': 'black', 's': 80},
}

plotted_labels = set()
for d in all_data:
    key = (d['category'], d['source'])
    s = style.get(key, {'color': 'gray', 'marker': 'x', 'edgecolor': 'black', 's': 60})

    label = None
    legend_key = f"{d['category'].capitalize()} ({d['source']})"
    if legend_key not in plotted_labels:
        label = legend_key
        plotted_labels.add(legend_key)

    ax.scatter(d['kappa_min'], d['kappa_meas'], zorder=3,
               c=s['color'], marker=s['marker'], edgecolors=s['edgecolor'],
               s=s['s'], linewidths=0.8, label=label)

    # Label each point
    offset = (5, 5)
    if 'LLZTO' in d['name'] and 'SC' in d['name']:
        offset = (5, -12)
    elif 'Cheng' in d['name']:
        offset = (-10, -14)
    ax.annotate(d['name'], (d['kappa_min'], d['kappa_meas']),
                textcoords="offset points", xytext=offset,
                fontsize=7, ha='left', va='bottom')

# Formatting
ax.set_xlabel(r'$\kappa_{\min}$ at 300 K (W m$^{-1}$ K$^{-1}$)', fontsize=12)
ax.set_ylabel(r'$\kappa_{\mathrm{measured}}$ at 300 K (W m$^{-1}$ K$^{-1}$)', fontsize=12)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.2, 5)
ax.set_ylim(0.2, 5)
ax.set_aspect('equal')
ax.tick_params(direction='in', which='both', top=True, right=True)
ax.legend(fontsize=9, loc='upper left', framealpha=0.8)

# Add annotation
ax.text(0.95, 0.05, r'$\kappa_{\mathrm{meas}} = \kappa_{\min}$',
        transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
        color='gray', style='italic')

plt.tight_layout()

out_path = Path(__file__).parent / 'kappa_min_vs_kappa_300K'
fig.savefig(str(out_path) + '.png', dpi=300, bbox_inches='tight')
print(f"\nSaved to {out_path}.png")
plt.close(fig)
