IO Readers
==========

MatSciKit provides readers for common materials science instruments.
All readers return NumPy arrays and use ``numpy.loadtxt``/``numpy.genfromtxt``
internally (no Pandas dependency).

PPMS TTO Reader
----------------

Reads thermal conductivity data from Quantum Design PPMS TTO ``.dat`` files.

.. code-block:: python

   from MatSciKit_COSMOTIM.io import ppms_tto

   # Basic read
   data = ppms_tto.read("sample_TTO.dat")
   # Returns: [Temperature (K), Conductivity (W/m/K), Error]

   # Drop specific temperature points
   data = ppms_tto.read("sample_TTO.dat", drop_temps=[302, 303])

   # Custom header skip
   data = ppms_tto.read("sample_TTO.dat", skip_rows=27)

PPMS HC Reader
---------------

Reads heat capacity data from Quantum Design PPMS HC ``.dat`` files.

.. code-block:: python

   from MatSciKit_COSMOTIM.io import ppms_hc

   data = ppms_hc.read("sample_HC.dat")
   # Returns: [Temperature (K), Cp (J/(g·K)), Cp_error]

   data = ppms_hc.read("sample_HC.dat", drop_temps=[50.0])

DSC Reader
-----------

Reads heat capacity data from Netzsch DSC CSV export files.
Automatically converts temperature from °C to K.

.. code-block:: python

   from MatSciKit_COSMOTIM.io import dsc

   data = dsc.read("ExpDat_sample.csv")
   # Returns: [Temperature (K), Cp (J/(g·K))]

   # Custom header rows (default is 34 for Netzsch format)
   data = dsc.read("other_dsc.csv", skip_rows=20)

LFA Reader
-----------

Reads data from Laser Flash Analysis CSV files.
Automatically drops rows containing NaN values.

.. code-block:: python

   from MatSciKit_COSMOTIM.io import lfa

   data = lfa.read("sample_LFA.csv")
   # Returns: [Temperature (K), Value, Error] with NaN rows removed

Data Format Summary
--------------------

.. list-table::
   :header-rows: 1
   :widths: 20 25 20 35

   * - Reader
     - Instrument
     - Input
     - Output columns
   * - ``ppms_tto``
     - QD PPMS TTO
     - ``.dat``
     - [T(K), κ(W/m/K), κ_err]
   * - ``ppms_hc``
     - QD PPMS HC
     - ``.dat``
     - [T(K), Cp(J/g/K), Cp_err]
   * - ``dsc``
     - Netzsch DSC
     - CSV
     - [T(K), Cp(J/g/K)]
   * - ``lfa``
     - LFA
     - CSV
     - [T(K), value, error]
