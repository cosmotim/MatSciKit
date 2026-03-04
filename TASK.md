# Build Task: MatSciKit Pipeline 1 + Scaffold

## Context
You are building a Python package called MatSciKit_COSMOTIM. Read `docs/DESIGN.md` first for the full architecture.

The MATLAB prototype source is at `../MatSciKit_Matlab_Prototype/` — use it as reference for translating functions.
Test data CSVs are at `../MatSciKit_Matlab_Prototype/TTO_Data/` and `../MatSciKit_Matlab_Prototype/ExpDat_LSHT-05232022-cool.csv`.

## Step 1: Scaffold the Package Structure

Restructure `src/MatSciKit_COSMOTIM/` to match the design in `docs/DESIGN.md`. Create all directories, `__init__.py` files, and move existing code:
- `xrdDataReader.py` → `structure/xrd_reader.py`
- `core/XRD_multiline_plot.py` → `structure/xrd_plot.py`
- `visualization/export_journal_style.py` → `visualization/journal_style.py`
- Delete empty files: `core/TTO`, `io/ppmsHCReader`
- Delete old empty directories: `core/`
- Update imports in moved files

## Step 2: Build `constants.py`

Physical constants (SI units):
```python
kb = 1.380649e-23   # Boltzmann constant (J/K)
h = 6.62607015e-34  # Planck constant (J·s)
hbar = h / (2 * np.pi)  # reduced Planck constant
```

## Step 3: Build `io/` Readers

### `io/ppms_hc.py`
Translate from `../MatSciKit_Matlab_Prototype/hcdataplot.m`:
- Read PPMS HC .dat files (skip first 15 header rows)
- Extract columns: [col5=Temperature, col10=Cp, col11=Cp_error] (MATLAB 1-indexed → Python 0-indexed: cols 4, 9, 10)
- Support dropping specific temperature points
- Return numpy array [T, Cp, Cp_err]
- NO plotting, NO hardcoded paths

### `io/ppms_tto.py`
Translate from `../MatSciKit_Matlab_Prototype/ttodataplot.m`:
- Read PPMS TTO .dat files (skip first 27 header rows)
- Extract columns 5,6,7 (0-indexed) = [Temperature, Conductivity, Error]
- Support dropping specific temperature points
- Return numpy array [T, κ, κ_err]
- NO plotting, NO hardcoded paths

### `io/dsc.py`
Based on how DSC data is read in `LSHT_Cp_functionized_main.m`:
- Read DSC CSV files (skip 34 header rows)
- Columns: [Temp_C, Time, Cp]
- Convert temperature from °C to K (add 273.15)
- Return numpy array [T_K, Cp]
- NO plotting

### `io/lfa.py`
- Read LFA CSV files (simple CSV, no header)
- Return numpy array [T, diffusivity/conductivity, error]
- Handle NaN rows (drop them)

## Step 4: Build `fitting/linear.py`

Translate from `../MatSciKit_Matlab_Prototype/Cp_and_TTO_functions/linear_fit_with_errors.m`:
- Weighted linear regression with error propagation
- Inputs: x, y, y_err
- Returns: slope, slope_error, intercept

## Step 5: Build `heat_capacity/` (Pipeline 1)

### `heat_capacity/low_t_fitting.py`
Translate from `../MatSciKit_Matlab_Prototype/Cp_and_TTO_functions/LowT_Cp_fitting.m`:
- Input: T, Cp, Cp_err, N_density, Density
- Perform weighted linear fit of Cp/T vs T²
- Calculate Debye temperature: θ_D = (slope * Density * 1e3 / (12π⁴/5 * N_density * kb))^(-1/3)
- Calculate sound velocity: v_s = θ_D / (ħ/kb * (6*N_density*π²)^(1/3))
- Return: theta_D, v_s, theta_D_error, v_s_error

### `heat_capacity/dulong_petit.py`
- Calculate Dulong-Petit limit: 3 * N_density * kb / Density * 1e-3
- Input: N_density, Density
- Return: limit value (J g⁻¹ K⁻¹)

### `heat_capacity/debye.py`
Translate from `modulus2debyeT.m` and `velocity2debyeT.m`:
- `from_modulus(modulus_GPa, density, N_density)` → θ_D
- `from_velocity(v_s, N_density)` → θ_D

## Step 6: Build `visualization/inset_style.py`

Translate from `../MatSciKit_Matlab_Prototype/Cp_and_TTO_functions/plotExportInsetStyle.m`.

## Step 7: Write Tests

Create `tests/test_pipeline1.py`:
- Load `../MatSciKit_Matlab_Prototype/TTO_Data/LSHT_Cp_all.csv` as HC data
- Load `../MatSciKit_Matlab_Prototype/ExpDat_LSHT-05232022-cool.csv` as DSC data
- Material params: V = 3.98^3 * 1e-30, N = (3/8 + 7/16 + 1/4 + 3/4 + 3), N_density = N/V, same M/Density calc as MATLAB
- Run low_t_fitting with points 13-41 (0-indexed: 12:41)
- Verify θ_D and v_s are reasonable (θ_D should be a few hundred K, v_s a few thousand m/s)
- Test Dulong-Petit limit calculation
- Test DSC reader (°C→K conversion)
- Test io readers with sample data

## Step 8: Update `pyproject.toml`

```toml
dependencies = ["numpy", "scipy", "matplotlib"]
```
Remove pandas from dependencies.

## Rules
- Use numpy.loadtxt / numpy.genfromtxt for CSV reading, NOT pandas
- Analysis functions must NEVER import matplotlib
- All functions need NumPy-style docstrings
- Add type hints
- MATLAB 1-based indexing → Python 0-based indexing
- Preserve original algorithm logic exactly
- Use `from MatSciKit_COSMOTIM.constants import kb, h, hbar` — don't redefine constants

## IMPORTANT
If you hit any issues with:
- Understanding the MATLAB code
- Data file format problems
- Test failures you can't resolve
STOP and describe the issue clearly. Do not guess or skip.

When completely finished, run this command to notify:
openclaw system event --text "Done: MatSciKit Pipeline 1 scaffold + io readers + heat capacity analysis + tests" --mode now
