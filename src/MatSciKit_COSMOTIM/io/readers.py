"""
Class-based IO readers for MatSciKit data files.

These readers wrap the function-based readers, adding metadata extraction
and a unified interface. The simple ``read()`` functions remain available
for quick one-liner usage.

Examples
--------
>>> reader = PPMSHCReader("data/sample_hc.dat")
>>> reader.metadata
{'filepath': 'data/sample_hc.dat', 'format': 'PPMS HC', ...}
>>> data = reader.get_data()
>>> data.shape
(200, 3)
"""

from __future__ import annotations

import abc
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class ReaderMetadata:
    """Metadata extracted from a data file.

    Attributes
    ----------
    filepath : str
        Path to the source file.
    format_name : str
        Human-readable format name.
    filename : str
        Base filename.
    n_rows : int
        Number of data rows.
    n_cols : int
        Number of data columns.
    columns : list of str
        Column names/descriptions.
    t_min : float
        Minimum temperature (K).
    t_max : float
        Maximum temperature (K).
    extra : dict
        Format-specific metadata.
    """

    filepath: str = ""
    format_name: str = ""
    filename: str = ""
    n_rows: int = 0
    n_cols: int = 0
    columns: list[str] = field(default_factory=list)
    t_min: float = 0.0
    t_max: float = 0.0
    extra: dict[str, Any] = field(default_factory=dict)


class BaseReader(abc.ABC):
    """Abstract base class for data file readers.

    Parameters
    ----------
    filepath : str or Path
        Path to the data file.
    """

    FORMAT_NAME: str = "Unknown"
    COLUMNS: list[str] = []

    def __init__(self, filepath: str | Path) -> None:
        self._filepath = Path(filepath)
        if not self._filepath.exists():
            raise FileNotFoundError(f"File not found: {self._filepath}")
        self._data: np.ndarray | None = None
        self._raw_header: list[str] = []

    @property
    def filepath(self) -> Path:
        """Path to the data file."""
        return self._filepath

    @abc.abstractmethod
    def _parse(self) -> np.ndarray:
        """Parse the file and return data array."""

    def get_data(self, **kwargs: Any) -> np.ndarray:
        """Load and return the data array.

        Results are cached after first call.

        Returns
        -------
        data : np.ndarray
            Parsed data array.
        """
        if self._data is None:
            self._data = self._parse(**kwargs)
        return self._data

    @property
    def data(self) -> np.ndarray:
        """Shorthand for ``get_data()``."""
        return self.get_data()

    @property
    def metadata(self) -> ReaderMetadata:
        """Extract metadata from the loaded data.

        Returns
        -------
        meta : ReaderMetadata
            File and data metadata.
        """
        d = self.get_data()
        return ReaderMetadata(
            filepath=str(self._filepath),
            format_name=self.FORMAT_NAME,
            filename=self._filepath.name,
            n_rows=d.shape[0],
            n_cols=d.shape[1],
            columns=self.COLUMNS[: d.shape[1]],
            t_min=float(np.nanmin(d[:, 0])),
            t_max=float(np.nanmax(d[:, 0])),
            extra=self._extra_metadata(),
        )

    def _extra_metadata(self) -> dict[str, Any]:
        """Override to provide format-specific metadata."""
        return {}

    def _read_header(self, n_lines: int) -> list[str]:
        """Read the first n lines of the file as header."""
        if not self._raw_header:
            with open(self._filepath, encoding="latin-1") as f:
                self._raw_header = [f.readline() for _ in range(n_lines)]
        return self._raw_header

    def summary(self) -> str:
        """Return a human-readable summary string.

        Returns
        -------
        text : str
            Formatted summary of the loaded data.
        """
        m = self.metadata
        lines = [
            f"File:        {m.filename}",
            f"Format:      {m.format_name}",
            f"Rows:        {m.n_rows}",
            f"Columns:     {', '.join(m.columns)}",
            f"Temperature: {m.t_min:.2f} – {m.t_max:.2f} K",
        ]
        for k, v in m.extra.items():
            lines.append(f"{k}: {v}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({str(self._filepath)!r})"


class PPMSHCReader(BaseReader):
    """Reader for Quantum Design PPMS Heat Capacity (.dat) files.

    Parameters
    ----------
    filepath : str or Path
        Path to the PPMS HC .dat file.
    skip_rows : int, optional
        Number of header rows to skip. Default is 15.
    drop_temps : list of float, optional
        Temperatures at which to drop data points.

    Examples
    --------
    >>> reader = PPMSHCReader("LSHT_HC.dat")
    >>> print(reader.summary())
    >>> data = reader.data  # [T, Cp, Cp_err]
    """

    FORMAT_NAME = "PPMS HC"
    COLUMNS = ["Temperature (K)", "Cp (J/(g·K))", "Cp error"]

    def __init__(
        self,
        filepath: str | Path,
        skip_rows: int = 15,
        drop_temps: list[float] | None = None,
    ) -> None:
        super().__init__(filepath)
        self._skip_rows = skip_rows
        self._drop_temps = drop_temps

    def _parse(self, **kwargs: Any) -> np.ndarray:
        from MatSciKit_COSMOTIM.io import ppms_hc

        return ppms_hc.read(
            str(self._filepath),
            drop_temps=self._drop_temps,
            skip_rows=self._skip_rows,
        )

    def _extra_metadata(self) -> dict[str, Any]:
        header = self._read_header(self._skip_rows)
        info: dict[str, Any] = {"header_rows": self._skip_rows}
        for line in header:
            if "BYAPP" in line:
                info["application"] = line.strip().split(",")[-1]
            elif "SAMPLE" in line or "Title" in line.lower():
                info["sample"] = line.strip()
        return info


class PPMSTTOReader(BaseReader):
    """Reader for Quantum Design PPMS Thermal Transport (.dat) files.

    Parameters
    ----------
    filepath : str or Path
        Path to the PPMS TTO .dat file.
    skip_rows : int, optional
        Number of header rows to skip. Default is 27.
    drop_temps : list of float, optional
        Temperatures at which to drop data points.

    Examples
    --------
    >>> reader = PPMSTTOReader("LSHT_TTO.dat")
    >>> print(reader.summary())
    >>> data = reader.data  # [T, κ, κ_err]
    """

    FORMAT_NAME = "PPMS TTO"
    COLUMNS = ["Temperature (K)", "κ (W/m/K)", "κ error"]

    def __init__(
        self,
        filepath: str | Path,
        skip_rows: int = 27,
        drop_temps: list[float] | None = None,
    ) -> None:
        super().__init__(filepath)
        self._skip_rows = skip_rows
        self._drop_temps = drop_temps

    def _parse(self, **kwargs: Any) -> np.ndarray:
        from MatSciKit_COSMOTIM.io import ppms_tto

        return ppms_tto.read(
            str(self._filepath),
            drop_temps=self._drop_temps,
            skip_rows=self._skip_rows,
        )

    def _extra_metadata(self) -> dict[str, Any]:
        header = self._read_header(self._skip_rows)
        info: dict[str, Any] = {"header_rows": self._skip_rows}
        for line in header:
            if "BYAPP" in line:
                info["application"] = line.strip().split(",")[-1]
        return info


class DSCReader(BaseReader):
    """Reader for DSC (Differential Scanning Calorimetry) CSV files.

    Parameters
    ----------
    filepath : str or Path
        Path to the DSC CSV file.
    skip_rows : int, optional
        Number of header rows to skip. Default is 34 (Netzsch format).

    Examples
    --------
    >>> reader = DSCReader("LSHT_DSC.csv")
    >>> print(reader.summary())
    >>> data = reader.data  # [T (K), Cp (J/(g·K))]
    """

    FORMAT_NAME = "DSC"
    COLUMNS = ["Temperature (K)", "Cp (J/(g·K))"]

    def __init__(self, filepath: str | Path, skip_rows: int = 34) -> None:
        super().__init__(filepath)
        self._skip_rows = skip_rows

    def _parse(self, **kwargs: Any) -> np.ndarray:
        from MatSciKit_COSMOTIM.io import dsc

        return dsc.read(str(self._filepath), skip_rows=self._skip_rows)

    def _extra_metadata(self) -> dict[str, Any]:
        header = self._read_header(min(self._skip_rows, 10))
        info: dict[str, Any] = {"header_rows": self._skip_rows}
        for line in header:
            stripped = line.strip()
            if "instrument" in stripped.lower() or "netzsch" in stripped.lower():
                info["instrument"] = stripped
            elif "sample" in stripped.lower():
                info["sample"] = stripped
        return info


class LFAReader(BaseReader):
    """Reader for LFA (Laser Flash Analysis) CSV files.

    Parameters
    ----------
    filepath : str or Path
        Path to the LFA CSV file.

    Examples
    --------
    >>> reader = LFAReader("LSHT_LFA.csv")
    >>> print(reader.summary())
    >>> data = reader.data  # [T (K), Value, Error]
    """

    FORMAT_NAME = "LFA"
    COLUMNS = ["Temperature (K)", "Value", "Error"]

    def _parse(self, **kwargs: Any) -> np.ndarray:
        from MatSciKit_COSMOTIM.io import lfa

        return lfa.read(str(self._filepath))

    def _extra_metadata(self) -> dict[str, Any]:
        d = self.get_data()
        original = np.genfromtxt(self._filepath, delimiter=",")
        n_nan = original.shape[0] - d.shape[0]
        return {"rows_dropped_nan": n_nan} if n_nan > 0 else {}


def auto_reader(filepath: str | Path) -> BaseReader:
    """Auto-detect format and return the appropriate reader.

    Tries to determine the file format from the extension and header content.

    Parameters
    ----------
    filepath : str or Path
        Path to the data file.

    Returns
    -------
    reader : BaseReader
        An instance of the appropriate reader class.

    Raises
    ------
    ValueError
        If the format cannot be determined.
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()

    if suffix == ".csv":
        # Peek at header to distinguish DSC vs LFA
        with open(filepath, encoding="latin-1") as f:
            first_lines = [f.readline() for _ in range(5)]
        header_text = "\n".join(first_lines).lower()

        if "netzsch" in header_text or "dsc" in header_text or "cp" in header_text:
            return DSCReader(filepath)
        return LFAReader(filepath)

    if suffix == ".dat":
        # Peek at header to distinguish HC vs TTO
        with open(filepath, encoding="latin-1") as f:
            first_lines = [f.readline() for _ in range(30)]
        header_text = "\n".join(first_lines).lower()

        if "thermal transport" in header_text or "tto" in header_text:
            return PPMSTTOReader(filepath)
        return PPMSHCReader(filepath)

    raise ValueError(
        f"Cannot auto-detect format for '{filepath.name}'. Use a specific reader class instead."
    )
