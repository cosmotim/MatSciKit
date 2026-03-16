"""
Thermal conductivity analysis (Pipeline 3).

Modules
-------
cahill : Cahill minimum thermal conductivity model
mean_free_path : Phonon mean free path from thermal conductivity
porosity_correction : Solid κ from porous κ (Maxwell-Eucken)
lfa_dsc : κ = Cp × α × ρ (LFA diffusivity → conductivity)
gruneisen : Grüneisen parameter calculation
pipeline : High-level pipeline (TTO reading, LFA→κ, comparison)
"""

from __future__ import annotations

from . import cahill, gruneisen, lfa_dsc, mean_free_path, pipeline, porosity_correction

__all__ = [
    "cahill",
    "gruneisen",
    "lfa_dsc",
    "mean_free_path",
    "pipeline",
    "porosity_correction",
]
