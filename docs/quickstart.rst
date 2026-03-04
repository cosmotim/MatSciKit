Quick Start
===========

This guide walks through a typical thermal analysis workflow using MatSciKit.

Define Material Parameters
--------------------------

Before running any analysis, you need to define the material's structural
parameters. These can be calculated from crystallographic data (CIF file
support is planned for Pipeline 0).

.. code-block:: python

   import numpy as np

   # Unit cell volume (m³)
   V = 3.98**3 * 1e-30

   # Atoms per unit cell and number density
   N = 3/8 + 7/16 + 1/4 + 3/4 + 3
   n_density = N / V  # atoms/m³

   # Mass density (kg/m³)
   M = (6.4*3/8 + 87.62*7/16 + 178.49*1/4 + 180.95*3/4 + 16*3) * 1.66e-27
   density = M / V

Pipeline 1: Heat Capacity
--------------------------

Extract Debye temperature and sound velocity from low-temperature Cp data.

.. code-block:: python

   from MatSciKit_COSMOTIM.heat_capacity import low_t_fitting, dulong_petit

   # Load heat capacity data [T, Cp, Cp_err]
   hc_data = np.loadtxt("LSHT_Cp_all.csv", delimiter=",")
   T, Cp, Cp_err = hc_data[:, 0], hc_data[:, 1], hc_data[:, 2]

   # Fit low-T region using temperature range
   theta_D, v_s, err_D, err_v = low_t_fitting.fit(
       T, Cp, Cp_err, n_density, density,
       t_range=(3.0, 10.0)
   )

   # Or using index range (1-based inclusive, like MATLAB)
   theta_D, v_s, err_D, err_v = low_t_fitting.fit(
       T, Cp, Cp_err, n_density, density,
       n_range=(13, 41)
   )

   print(f"θ_D = {theta_D:.1f} ± {err_D:.1f} K")
   print(f"v_s = {v_s:.1f} ± {err_v:.1f} m/s")

   # Dulong-Petit limit
   dp = dulong_petit.calculate(n_density, density)
   print(f"Dulong-Petit limit: {dp:.4f} J/(g·K)")

Pipeline 2: Thermal Conductivity
----------------------------------

Use θ_D and v_s from Pipeline 1 to analyze thermal conductivity data.

.. code-block:: python

   from MatSciKit_COSMOTIM.io import ppms_tto
   from MatSciKit_COSMOTIM.thermal_conductivity import (
       cahill, mean_free_path, porosity_correction
   )

   # Load TTO data
   tto_data = ppms_tto.read("sample_TTO.dat", drop_temps=[302, 303])
   T_tto = tto_data[:, 0]
   kappa = tto_data[:, 1]

   # Correct for porosity
   kappa_solid = porosity_correction.correct(kappa, porosity=0.03)

   # Cahill minimum thermal conductivity
   T_theory = np.linspace(1, 1000, 500)
   kappa_min = cahill.minimum_tc(T_theory, n_density, theta_D, v_s)

   # Mean free path
   mfp = mean_free_path.calculate(T_tto, kappa, theta_D, v_s)

LFA + DSC Bridge
-----------------

When thermal conductivity must be calculated from LFA diffusivity and DSC
heat capacity (κ = Cp × α × ρ):

.. code-block:: python

   from MatSciKit_COSMOTIM.thermal_conductivity import lfa_dsc

   tc_data, cp_avg = lfa_dsc.calculate(
       cp_T, cp, cp_error,
       diff_T, diffusivity, diff_error,
       density=6.87, density_error=0.28
   )
   # tc_data: [Temperature, κ, κ_error]

Visualization
-------------

Export figures with journal-style formatting:

.. code-block:: python

   import matplotlib.pyplot as plt
   from MatSciKit_COSMOTIM.visualization.journal_style import export_journal_figure

   fig, ax = plt.subplots()
   ax.plot(T_tto, kappa_solid, 'o', label='Measurement')
   ax.plot(T_theory, kappa_min, '--', label='Cahill minimum')
   ax.set_xlabel("T (K)")
   ax.set_ylabel("κ (W m⁻¹ K⁻¹)")
   ax.legend()

   # Exports with Arial 11pt, inward ticks, 600 DPI
   export_journal_figure(fig, "thermal_conductivity", format="png")
