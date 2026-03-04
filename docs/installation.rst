Installation
============

Requirements
------------

- Python 3.8 or later
- NumPy
- SciPy
- Matplotlib

Install from Source
-------------------

.. code-block:: bash

   git clone https://github.com/cosmotim/MatSciKit.git
   cd MatSciKit
   pip install -e .

For development (includes pytest):

.. code-block:: bash

   pip install -e ".[dev]"

Dependencies
------------

MatSciKit uses only three core dependencies:

- **NumPy** — array operations, linear algebra, and data I/O
- **SciPy** — numerical integration (Debye integrals) and curve fitting
- **Matplotlib** — publication-quality figure export

Pandas is intentionally not required. All data readers use
``numpy.loadtxt`` / ``numpy.genfromtxt`` since instrument data files
are well-structured.
