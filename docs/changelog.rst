Changelog
=========

v0.2.0 (2026-03-03)
--------------------

Major restructuring and new features:

**New Package Structure:**

- Reorganized into ``io/``, ``structure/``, ``heat_capacity/``,
  ``thermal_conductivity/``, ``fitting/``, ``visualization/``
- Separated I/O from analysis from plotting
- Dropped Pandas dependency (NumPy-only for data reading)

**Pipeline 1 — Heat Capacity:**

- ``heat_capacity.low_t_fitting`` — Debye temperature and sound velocity
  from low-T Cp/T vs T² fitting (supports ``t_range`` and ``n_range``)
- ``heat_capacity.dulong_petit`` — Classical heat capacity limit
- ``heat_capacity.debye`` — θ_D from velocity or bulk modulus

**Pipeline 2 — Thermal Conductivity:**

- ``thermal_conductivity.cahill`` — Cahill minimum κ model
- ``thermal_conductivity.mean_free_path`` — Phonon MFP from kinetic theory
- ``thermal_conductivity.porosity_correction`` — Maxwell-Eucken correction
- ``thermal_conductivity.lfa_dsc`` — κ = Cp × α × ρ with error propagation
- ``thermal_conductivity.gruneisen`` — Grüneisen parameter

**IO Readers:**

- ``io.ppms_tto`` — PPMS TTO .dat reader
- ``io.ppms_hc`` — PPMS HC .dat reader
- ``io.dsc`` — DSC CSV reader (auto °C→K conversion)
- ``io.lfa`` — LFA CSV reader (auto NaN removal)

**Other:**

- ``constants.py`` — Physical constants (kb, h, ħ)
- ``fitting.linear`` — Weighted linear regression with errors
- ``visualization.inset_style`` — Inset panel figure export

v0.1.0 (2025-12)
-----------------

Initial release:

- XRD data reader and multi-line plotter
- PPMS TTO reader (basic)
- Journal-style figure export
