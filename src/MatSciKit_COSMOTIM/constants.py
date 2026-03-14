"""
Physical constants in SI units.

Values are sourced from :mod:`scipy.constants` (CODATA 2018).
Aliases are provided for backward compatibility.
"""

from __future__ import annotations

from scipy import constants as _c

# Boltzmann constant (J/K)
kb: float = _c.Boltzmann

# Planck constant (J·s)
h: float = _c.Planck

# Reduced Planck constant ħ (J·s)
hbar: float = _c.hbar

# Avogadro constant (mol⁻¹)
NA: float = _c.Avogadro

# Atomic mass unit (kg)
AMU_TO_KG: float = _c.atomic_mass
