"""
Instrument-specific data readers.

Function-based readers (quick use):
    ``ppms_hc.read()``, ``ppms_tto.read()``, ``dsc.read()``, ``lfa.read()``

Class-based readers (metadata + caching):
    ``PPMSHCReader``, ``PPMSTTOReader``, ``DSCReader``, ``LFAReader``,
    ``auto_reader()``
"""

from __future__ import annotations

from . import dsc, lfa, lfa_excel, ppms_hc, ppms_tto
from .readers import (
    BaseReader,
    DSCReader,
    LFAReader,
    PPMSHCReader,
    PPMSTTOReader,
    auto_reader,
)

__all__ = [
    "BaseReader",
    "DSCReader",
    "LFAReader",
    "PPMSHCReader",
    "PPMSTTOReader",
    "auto_reader",
    "dsc",
    "lfa",
    "lfa_excel",
    "ppms_hc",
    "ppms_tto",
]
