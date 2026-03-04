# MatSciKit Design Document

**Date:** 2026-03-03  
**Authors:** Yitian Wang, Coder Agent  
**Status:** Agreed — ready to implement

---

## Overview

MatSciKit is a Python package for materials science research data processing, with a focus on thermal transport analysis. It translates a MATLAB prototype (originally developed for LSHT perovskite solid electrolyte analysis) into a clean, reusable, modular Python package.

## MATLAB Prototype

The original MATLAB codebase lives in Google Drive ("MatSciKit Matlab Prototype") and centers on `LSHT_Cp_functionized_main.m`, which orchestrates the full thermal analysis workflow. The MATLAB code combines heat capacity analysis and thermal conductivity analysis into a single monolithic function (`Cp_functionized`), with hardcoded material parameters, paths, and tightly coupled plotting.

### MATLAB Function Inventory

| Function | Purpose | Python Target |
|---|---|---|
| `ttodataplot.m` | Read PPMS TTO .dat files | `io/ppms_tto.py` |
| `hcdataplot.m` | Read PPMS HC .dat files | `io/ppms_hc.py` |
| `LowT_Cp_fitting.m` | Debye T & v_s from low-T Cp fit | `heat_capacity/low_t_fitting.py` |
| `calculate_cahill_TC.m` | Cahill minimum thermal conductivity | `thermal_conductivity/cahill.py` |
| `solidTCwithPorosity.m` | Porosity correction for κ | `thermal_conductivity/porosity_correction.py` |
| `MFP_from_TC.m` | Mean free path from κ | `thermal_conductivity/mean_free_path.py` |
| `LFA_DSC_TC.m` | κ = Cp × α × ρ | `thermal_conductivity/lfa_dsc.py` |
| `Gruneisen_calculate.m` | Grüneisen parameter | `thermal_conductivity/gruneisen.py` |
| `modulus2debyeT.m` | Bulk modulus → Debye T | `heat_capacity/debye.py` |
| `velocity2debyeT.m` | Sound velocity → Debye T | `heat_capacity/debye.py` |
| `linear_fit_with_errors.m` | Weighted linear regression | `fitting/linear.py` |
| `getAtomicMass.m` | Periodic table lookup | `structure/atomic_masses.py` |
| `calculateDensityFromCIF.m` | Density from CIF file | `structure/cif_reader.py` (future) |
| `plotExportJournalStyle.m` | Journal figure formatting | `visualization/journal_style.py` |
| `plotExportInsetStyle.m` | Inset figure formatting | `visualization/inset_style.py` |
| `Porosity_calculator.m` | Archimedes density/porosity | `structure/porosity.py` |
| `Cp_functionized.m` | Combined Cp+κ analysis (to be split) | `heat_capacity/` + `thermal_conductivity/` |
| `XRDDataReader` (existing Python) | XRD data file reader | `structure/xrd_reader.py` |
| `XRDMultilinePlotter` (existing Python) | Multi-pattern XRD plots | `structure/xrd_plot.py` |

## Architecture: Three Pipelines

```
Pipeline 0: Structure (XRD + CIF)  →  Material properties (V, ρ, N/V)
     ↓
Pipeline 1: Heat Capacity (Cp)     →  θ_D, v_s, Dulong-Petit limit
     ↓
Pipeline 2: Thermal Conductivity   →  κ_solid, κ_min (Cahill), MFP
```

### Pipeline 0 — Structure Analysis (future)
- **Input:** CIF files, XRD data
- **Output:** Material object with V, ρ, N/V, N
- **Status:** XRD reader/plotter exists. CIF reader is future work. For now, material parameters are provided manually.

### Pipeline 1 — Heat Capacity Analysis
- **Input:** PPMS HC data [T, Cp, Cp_err] + optional DSC data [T, Cp]
- **Processing:**
  1. Low-T region: Cp/T vs T² weighted linear fit
  2. Extract Debye temperature (θ_D) and average sound velocity (v_s)
  3. Compute Dulong-Petit limit: `3 × N_density × kb / Density`
  4. Combine PPMS + DSC for full temperature range
- **Output:** θ_D, v_s (with errors), Cp(T)

### Pipeline 2 — Thermal Conductivity Analysis
- **Input:** PPMS TTO data [T, κ, κ_err] + optional LFA data
- **Depends on:** θ_D, v_s from Pipeline 1; Cp(T) for LFA path
- **Processing:**
  1. Porosity correction: κ_s = κ × (2+φ)/(2−2φ)
  2. Cahill minimum κ: Debye integral model
  3. Mean free path: MFP from kinetic theory
  4. LFA→κ conversion: κ = Cp × α × ρ (needs P1 Cp data)
- **Output:** κ_solid(T), κ_min(T), MFP(T)

## Package Structure

```
src/MatSciKit_COSMOTIM/
├── __init__.py
├── constants.py                     # kb, h, ħ — physical constants
│
├── io/                              # Instrument-specific data readers
│   ├── __init__.py
│   ├── ppms_tto.py                  # PPMS TTO .dat reader
│   ├── ppms_hc.py                   # PPMS HC .dat reader
│   ├── dsc.py                       # DSC CSV reader (°C→K)
│   └── lfa.py                       # LFA Excel/CSV reader
│
├── structure/                       # Crystallography & XRD (Pipeline 0)
│   ├── __init__.py
│   ├── material.py                  # Material class (manual or from CIF)
│   ├── atomic_masses.py             # Periodic table dict
│   ├── cif_reader.py                # CIF file parser (future)
│   ├── porosity.py                  # Archimedes method density/porosity
│   ├── xrd_reader.py               # XRDDataReader class
│   └── xrd_plot.py                  # XRDMultilinePlotter class
│
├── heat_capacity/                   # Pipeline 1
│   ├── __init__.py
│   ├── low_t_fitting.py             # Debye T, v_s from Cp/T vs T² fit
│   ├── dulong_petit.py              # Dulong-Petit limit calculation
│   └── debye.py                     # modulus→θ_D, velocity→θ_D converters
│
├── thermal_conductivity/            # Pipeline 2
│   ├── __init__.py
│   ├── cahill.py                    # Cahill minimum κ
│   ├── mean_free_path.py            # MFP from κ
│   ├── porosity_correction.py       # κ_solid from porous κ
│   ├── lfa_dsc.py                   # κ = Cp × α × ρ (P1→P2 bridge)
│   └── gruneisen.py                 # Grüneisen parameter
│
├── fitting/                         # Curve fitting utilities
│   ├── __init__.py
│   └── linear.py                    # Weighted linear fit with errors
│
└── visualization/                   # Plotting (matplotlib)
    ├── __init__.py
    ├── journal_style.py             # Journal figure export
    └── inset_style.py               # Inset figure export
```

## Dependencies

```toml
dependencies = ["numpy", "scipy", "matplotlib"]
```

- **NumPy:** Core arrays, linear algebra, vectorized operations
- **SciPy:** Integration (Debye integrals), curve fitting (weighted least squares), signal processing
- **Matplotlib:** Figure generation and export
- No Pandas — CSV files are clean enough for `numpy.loadtxt`/`numpy.genfromtxt`

## Design Principles

1. **Separate I/O from analysis from plotting** — analysis functions never import matplotlib
2. **Analysis functions return data** — they compute and return, never plot
3. **Physical constants in one place** — `constants.py`, not redefined per function
4. **Material parameters from crystallography** — CIF reader (future) or manual input, not hardcoded
5. **No hardcoded paths** — all file paths are user-provided arguments
6. **Porosity is sample-specific** — not a property of the Material class
7. **Pipelines are independently runnable** — P2 accepts θ_D/v_s as arguments, doesn't require running P1 in the same session

## Future Work

- [ ] CIF file reader for automated Material construction (Pipeline 0)
- [ ] Two-channel thermal conductivity fitting model
- [ ] Vectorized variables for inner product and matrix operations (for two-channel fitting)
- [ ] XRDML file format support
- [ ] Additional instrument readers as needed

## Example Usage (Target)

```python
from MatSciKit_COSMOTIM.io import ppms_hc, ppms_tto
from MatSciKit_COSMOTIM.structure import Material
from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting, dulong_petit
from MatSciKit_COSMOTIM.thermal_conductivity import cahill, mean_free_path, porosity_correction

# Material definition (manual until CIF reader exists)
lsht = Material(a=3.98, atoms={'Li': 3/8, 'Sr': 7/16, 'Hf': 1/4, 'Ta': 3/4, 'O': 3})

# Pipeline 1: Heat Capacity
hc_data = ppms_hc.read("LSHT_HC.dat")
theta_D, v_s, errors = low_t_fitting.fit(hc_data, lsht.n_density, lsht.density)
dp_limit = dulong_petit.calculate(lsht.n_density, lsht.density)

# Pipeline 2: Thermal Conductivity
tto_data = ppms_tto.read("LSHT_TTO.dat", drop_temps=[302, 303])
kappa_solid = porosity_correction.correct(tto_data[:, 1], porosity=0.03)
kappa_min = cahill.minimum_tc(T, lsht.n_density, theta_D, v_s)
mfp = mean_free_path.calculate(tto_data, theta_D, v_s)
```
