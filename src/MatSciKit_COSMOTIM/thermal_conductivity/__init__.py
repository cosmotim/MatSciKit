"""
Thermal conductivity analysis (Pipeline 2).

Modules
-------
cahill : Cahill minimum thermal conductivity model
mean_free_path : Phonon mean free path from thermal conductivity
porosity_correction : Solid κ from porous κ (Maxwell-Eucken)
lfa_dsc : κ = Cp × α × ρ (bridge between Pipeline 1 and 2)
gruneisen : Grüneisen parameter calculation
"""

from __future__ import annotations

from . import cahill, gruneisen, lfa_dsc, mean_free_path, porosity_correction

__all__ = ["cahill", "gruneisen", "lfa_dsc", "mean_free_path", "porosity_correction"]
