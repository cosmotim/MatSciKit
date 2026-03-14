# MatSciKit

A Python toolkit for materials science research data processing, with a focus on thermal transport analysis.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

MatSciKit provides a modular pipeline for analyzing thermal properties of materials, from raw instrument data to publication-ready results. It was designed to replace a collection of MATLAB scripts used in thermal transport research on solid-state ionic conductors and functional oxides.

### Three-Pipeline Architecture

```
Pipeline 0: Structure (XRD + CIF)  →  Material properties (V, ρ, N/V)
     ↓
Pipeline 1: Heat Capacity (Cp)     →  Debye temperature (θ_D), sound velocity (v_s)
     ↓
Pipeline 2: Thermal Conductivity   →  κ_solid, κ_min (Cahill), mean free path
```

## Installation

### From PyPI (recommended)

```bash
pip install matscikit-cosmotim
```

With interactive terminal interface:
```bash
pip install "matscikit-cosmotim[tui]"
```

### From GitHub (latest development version)

```bash
pip install git+https://github.com/cosmotim/MatSciKit.git
```

### For Development

```bash
git clone https://github.com/cosmotim/MatSciKit.git
cd MatSciKit
pip install -e ".[dev,tui]"
```

### Dependencies

- **NumPy** — array operations and linear algebra
- **SciPy** — numerical integration and curve fitting
- **Matplotlib** — publication-quality figure export

## Quick Start

### Pipeline 1: Heat Capacity Analysis

Extract Debye temperature and sound velocity from low-temperature heat capacity data:

```python
import numpy as np
from MatSciKit_COSMOTIM.io import ppms_hc, dsc
from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting, dulong_petit

# Load PPMS heat capacity data
hc_data = np.loadtxt("LSHT_Cp_all.csv", delimiter=",")
T, Cp, Cp_err = hc_data[:, 0], hc_data[:, 1], hc_data[:, 2]

# Material parameters
V = 3.98**3 * 1e-30          # unit cell volume (m³)
N = 3/8 + 7/16 + 1/4 + 3/4 + 3  # atoms per unit cell
n_density = N / V             # number density (atoms/m³)
density = 6870                # mass density (kg/m³)

# Fit low-T region: Cp/T vs T²
# Method 1: by temperature range
theta_D, v_s, err_D, err_v = low_t_fitting.fit(
    T, Cp, Cp_err, n_density, density,
    t_range=(3.0, 10.0)
)

# Method 2: by index range (1-based, like MATLAB)
theta_D, v_s, err_D, err_v = low_t_fitting.fit(
    T, Cp, Cp_err, n_density, density,
    n_range=(13, 41)
)

print(f"Debye temperature: {theta_D:.1f} ± {err_D:.1f} K")
print(f"Sound velocity: {v_s:.1f} ± {err_v:.1f} m/s")

# Dulong-Petit limit
dp_limit = dulong_petit.calculate(n_density, density)
print(f"Dulong-Petit limit: {dp_limit:.4f} J/(g·K)")
```

### Pipeline 2: Thermal Conductivity Analysis

Compute Cahill minimum thermal conductivity, porosity correction, and mean free path:

```python
from MatSciKit_COSMOTIM.io import ppms_tto, lfa
from MatSciKit_COSMOTIM.thermal_conductivity import (
    cahill, mean_free_path, porosity_correction
)

# Load TTO data, dropping bad temperature points
tto_data = ppms_tto.read("sample_TTO.dat", drop_temps=[302, 303])
T_tto, kappa, kappa_err = tto_data[:, 0], tto_data[:, 1], tto_data[:, 2]

# Porosity correction: κ_solid from measured κ
kappa_solid = porosity_correction.correct(kappa, porosity=0.03)

# Cahill minimum thermal conductivity
T_plot = np.linspace(1, 1000, 500)
kappa_min = cahill.minimum_tc(T_plot, n_density, theta_D, v_s)

# Mean free path
mfp = mean_free_path.calculate(T_tto, kappa, theta_D, v_s)
```

### LFA + DSC → Thermal Conductivity

When LFA measures diffusivity and DSC measures heat capacity:

```python
from MatSciKit_COSMOTIM.thermal_conductivity import lfa_dsc

# κ = Cp × α × ρ (with error propagation)
tc_data, cp_avg = lfa_dsc.calculate(
    cp_T, cp, cp_error,           # DSC heat capacity
    diff_T, diffusivity, diff_error,  # LFA diffusivity
    density=6.87, density_error=0.28
)
# tc_data columns: [Temperature, κ, κ_error]
```

### Visualization

Export publication-ready figures:

```python
from MatSciKit_COSMOTIM.visualization.journal_style import export_journal_figure
from MatSciKit_COSMOTIM.visualization.inset_style import export_inset_figure

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(T, kappa, 'o')
ax.set_xlabel("T (K)")
ax.set_ylabel("κ (W m⁻¹ K⁻¹)")

# Export with journal formatting (Arial 11pt, 600 DPI, inward ticks)
export_journal_figure(fig, "thermal_conductivity", format="png")
```

## Package Structure

```
src/MatSciKit_COSMOTIM/
├── constants.py                 # Physical constants (kb, h, ħ)
├── io/                          # Instrument data readers
│   ├── ppms_tto.py              #   PPMS Thermal Transport Option
│   ├── ppms_hc.py               #   PPMS Heat Capacity
│   ├── dsc.py                   #   Differential Scanning Calorimetry
│   └── lfa.py                   #   Laser Flash Analysis
├── structure/                   # Crystallography & XRD (Pipeline 0)
│   ├── xrd_reader.py            #   XRD data file reader
│   └── xrd_plot.py              #   Multi-pattern XRD plotter
├── heat_capacity/               # Pipeline 1
│   ├── low_t_fitting.py         #   Debye T & v_s from Cp/T vs T²
│   ├── dulong_petit.py          #   Classical heat capacity limit
│   └── debye.py                 #   θ_D from velocity or modulus
├── thermal_conductivity/        # Pipeline 2
│   ├── cahill.py                #   Cahill minimum κ model
│   ├── mean_free_path.py        #   Phonon MFP from κ
│   ├── porosity_correction.py   #   Solid κ from porous samples
│   ├── lfa_dsc.py               #   κ = Cp × α × ρ
│   └── gruneisen.py             #   Grüneisen parameter
├── fitting/                     # Curve fitting utilities
│   └── linear.py                #   Weighted linear regression
└── visualization/               # Publication figures
    ├── journal_style.py         #   Journal figure formatting
    └── inset_style.py           #   Inset panel formatting
```

## IO Readers

| Module | Instrument | Input Format | Output |
|--------|-----------|-------------|--------|
| `ppms_tto` | Quantum Design PPMS TTO | `.dat` file | [T, κ, κ_err] |
| `ppms_hc` | Quantum Design PPMS HC | `.dat` file | [T, Cp, Cp_err] |
| `dsc` | Netzsch DSC 214 | CSV export | [T(K), Cp] |
| `lfa` | Laser Flash Analysis | CSV | [T, value, error] |

## Interactive Terminal Interface

MatSciKit includes an interactive TUI for loading and analyzing data without writing code:

```bash
# Launch the TUI (requires tui extra)
pip install "matscikit-cosmotim[tui]"
matscikit
```

Features:
- File browser to select data files
- Auto-detection of PPMS HC, PPMS TTO, DSC, and LFA formats
- Loading progress and data preview
- Pipeline analysis selection

## Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific pipeline
python -m pytest tests/test_pipeline1.py -v
python -m pytest tests/test_pipeline2.py -v
```

## Roadmap

- [x] Pipeline 1: Heat capacity analysis
- [x] Pipeline 2: Thermal conductivity analysis
- [ ] Pipeline 0: CIF file reader for automated material properties
- [ ] Two-channel thermal conductivity fitting model
- [ ] XRDML file format support
- [ ] Example Jupyter notebooks

## Citation

If you use MatSciKit in your research, please cite:

```bibtex
@software{matscikit,
  author = {Wang, Yitian},
  title = {MatSciKit: Python Toolkit for Materials Science Research},
  url = {https://github.com/cosmotim/MatSciKit},
  year = {2026}
}
```

The methods implemented in this package are based on the following work:

```bibtex
@article{wang2025thermal,
  title = {Origin of Intrinsically Low Thermal Conductivity in a Garnet-Type Solid Electrolyte},
  author = {Wang, Yitian and Su, Yaokun and Carrete, Jes\'us and Zhang, Huanyu and Wu, Nan and Li, Yutao and Li, Hongze and He, Jiaming and Xu, Youming and Guo, Shucheng and Cai, Qingan and Abernathy, Douglas L. and Williams, Travis and Kravchyk, Kostiantyn V. and Kovalenko, Maksym V. and Madsen, Georg K.H. and Li, Chen and Chen, Xi},
  journal = {PRX Energy},
  volume = {4},
  pages = {033004},
  year = {2025},
  doi = {10.1103/6wj2-kzhh}
}
```

## License

MIT License — see [LICENSE](LICENSE) for details.
