MatSciKit Documentation
=======================

**MatSciKit** is a Python toolkit for materials science research data processing,
with a focus on thermal transport analysis.

It provides a modular pipeline for analyzing thermal properties of materials,
from raw instrument data to publication-ready results.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   pipelines
   io_readers

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/constants
   api/io
   api/structure
   api/heat_capacity
   api/thermal_conductivity
   api/fitting
   api/visualization

.. toctree::
   :maxdepth: 1
   :caption: Development

   design
   changelog

Citation
========

If you use MatSciKit in your research, please cite:

.. code-block:: bibtex

   @software{matscikit,
     author = {Wang, Yitian},
     title = {MatSciKit: Python Toolkit for Materials Science Research},
     url = {https://github.com/cosmotim/MatSciKit},
     year = {2026}
   }

The methods implemented in this package are based on:

   Y. Wang *et al.*, "Origin of Intrinsically Low Thermal Conductivity in a
   Garnet-Type Solid Electrolyte: Linking Lattice and Ionic Dynamics with
   Thermal Transport," *PRX Energy* **4**, 033004 (2025).
   `DOI: 10.1103/6wj2-kzhh <https://link.aps.org/doi/10.1103/6wj2-kzhh>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
