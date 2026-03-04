Pipeline Architecture
=====================

MatSciKit organizes thermal analysis into three sequential pipelines,
each building on the outputs of the previous one.

.. code-block:: text

   Pipeline 0: Structure (XRD + CIF)  →  Material properties (V, ρ, N/V)
        ↓
   Pipeline 1: Heat Capacity (Cp)     →  θ_D, v_s, Dulong-Petit limit
        ↓
   Pipeline 2: Thermal Conductivity   →  κ_solid, κ_min, MFP

Pipeline 0: Structure Analysis
-------------------------------

**Status:** Partially implemented (XRD reader/plotter available; CIF reader planned)

Provides material structural properties needed by Pipelines 1 and 2:

- Unit cell volume (*V*)
- Mass density (*ρ*)
- Number density (*N/V*)
- Atoms per unit cell (*N*)

Currently these must be specified manually. A future CIF file reader will
automate this step.

Pipeline 1: Heat Capacity
--------------------------

**Input:** PPMS HC data and/or DSC data

**Processing:**

1. Load heat capacity data [T, Cp, Cp_err]
2. Select low-temperature fitting region (by temperature range or index range)
3. Perform weighted linear fit of Cp/T vs T²
4. Extract Debye temperature (θ_D) and average sound velocity (v_s)
5. Compute Dulong-Petit classical limit

**Output:** θ_D, v_s (with uncertainties), Cp(T)

**Key equation:**

.. math::

   C_p/T = \beta T^2 + \gamma

where :math:`\beta = \frac{12\pi^4}{5} \frac{N_{density} \cdot k_B}{\rho \cdot \theta_D^3}`

Pipeline 2: Thermal Conductivity
----------------------------------

**Input:** PPMS TTO data and/or LFA data

**Depends on:** θ_D and v_s from Pipeline 1

**Processing:**

1. Load thermal conductivity data [T, κ, κ_err]
2. Apply porosity correction: :math:`\kappa_s = \kappa \cdot \frac{2+\phi}{2-2\phi}`
3. Compute Cahill minimum κ using Debye integral model
4. Calculate phonon mean free path from kinetic theory
5. (Optional) Compute κ from LFA diffusivity + DSC Cp: :math:`\kappa = C_p \cdot \alpha \cdot \rho`

**Output:** κ_solid(T), κ_min(T), MFP(T)

**Key equations:**

Cahill minimum:

.. math::

   \kappa_{min} = \left(\frac{\pi}{6}\right)^{1/3} k_B n^{2/3} \cdot 3 v_s \left(\frac{T}{\theta_D}\right)^2 \int_0^{\theta_D/T} \frac{x^3 e^x}{(e^x - 1)^2} dx

Mean free path:

.. math::

   \ell = \kappa \left[ \frac{k_B^4 T^3}{2\pi^2 v_s^2 \hbar^3} \int_0^{\theta_D/T} \frac{x^4 e^x}{(e^x - 1)^2} dx \right]^{-1}
