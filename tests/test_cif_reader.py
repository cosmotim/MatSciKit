"""
Tests for Pipeline 0: CIF reader (structure/cif_reader.py).

Tests the basic CIF parser against known crystal structures with
independently verified properties.
"""

import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from MatSciKit_COSMOTIM.structure import cif_reader


# ---- Test CIF content ----

SRTIO3_CIF = textwrap.dedent("""\
    data_SrTiO3
    _cell_length_a                    3.905
    _cell_length_b                    3.905
    _cell_length_c                    3.905
    _cell_angle_alpha                 90.000
    _cell_angle_beta                  90.000
    _cell_angle_gamma                 90.000
    _cell_volume                      59.55

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    Sr1 Sr 0.5 0.5 0.5 1.0
    Ti1 Ti 0.0 0.0 0.0 1.0
    O1  O  0.5 0.0 0.0 1.0
    O2  O  0.0 0.5 0.0 1.0
    O3  O  0.0 0.0 0.5 1.0
""")

NACL_CIF = textwrap.dedent("""\
    data_NaCl
    _cell_length_a                    5.6402
    _cell_length_b                    5.6402
    _cell_length_c                    5.6402
    _cell_angle_alpha                 90.000
    _cell_angle_beta                  90.000
    _cell_angle_gamma                 90.000
    _cell_volume                      179.43

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    Na1 Na 0.0 0.0 0.0 1.0
    Na2 Na 0.5 0.5 0.0 1.0
    Na3 Na 0.5 0.0 0.5 1.0
    Na4 Na 0.0 0.5 0.5 1.0
    Cl1 Cl 0.5 0.5 0.5 1.0
    Cl2 Cl 0.0 0.0 0.5 1.0
    Cl3 Cl 0.0 0.5 0.0 1.0
    Cl4 Cl 0.5 0.0 0.0 1.0
""")

# Partial occupancy test (like LSHT perovskite)
PARTIAL_OCC_CIF = textwrap.dedent("""\
    data_LSHT_like
    _cell_length_a                    3.98
    _cell_length_b                    3.98
    _cell_length_c                    3.98
    _cell_angle_alpha                 90.000
    _cell_angle_beta                  90.000
    _cell_angle_gamma                 90.000

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_occupancy
    Li1 Li 0.0 0.0 0.0 0.375
    Sr1 Sr 0.0 0.0 0.0 0.4375
    Hf1 Hf 0.5 0.5 0.5 0.25
    Ta1 Ta 0.5 0.5 0.5 0.75
    O1  O  0.5 0.5 0.0 1.0
    O2  O  0.5 0.0 0.5 1.0
    O3  O  0.0 0.5 0.5 1.0
""")

# No-occupancy column test
NO_OCC_CIF = textwrap.dedent("""\
    data_simple
    _cell_length_a 4.0
    _cell_length_b 4.0
    _cell_length_c 4.0
    _cell_angle_alpha 90.0
    _cell_angle_beta 90.0
    _cell_angle_gamma 90.0
    _cell_volume 64.0

    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    Fe1 Fe 0.0 0.0 0.0
    O1  O  0.5 0.5 0.5
""")


def _write_temp_cif(content: str) -> str:
    """Write CIF content to a temp file and return the path."""
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False)
    f.write(content)
    f.close()
    return f.name


class TestCIFReaderBasicParser:
    """Test the basic (non-pymatgen) CIF parser."""

    def test_srtio3_unit_cell(self):
        """SrTiO3: cubic perovskite, a = 3.905 Å."""
        path = _write_temp_cif(SRTIO3_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        assert result['unit_cell']['a'] == pytest.approx(3.905)
        assert result['unit_cell']['b'] == pytest.approx(3.905)
        assert result['unit_cell']['c'] == pytest.approx(3.905)
        assert result['volume'] == pytest.approx(59.55)

    def test_srtio3_atoms(self):
        """SrTiO3: 5 atoms (1 Sr + 1 Ti + 3 O)."""
        path = _write_temp_cif(SRTIO3_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        assert result['total_atoms'] == pytest.approx(5.0)

        symbols = [a['symbol'] for a in result['atoms']]
        assert symbols.count('Sr') == 1
        assert symbols.count('Ti') == 1
        assert symbols.count('O') == 3

    def test_srtio3_density(self):
        """SrTiO3 density should be ~5100 kg/m³."""
        path = _write_temp_cif(SRTIO3_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        # Known: ρ(SrTiO3) ≈ 5110-5130 kg/m³
        assert 4900 < result['density'] < 5300, \
            f"SrTiO3 density = {result['density']:.0f} kg/m³"
        print(f"SrTiO3 density: {result['density']:.0f} kg/m³")

    def test_srtio3_number_density(self):
        """SrTiO3 number density: 5 atoms / 59.55 Å³."""
        path = _write_temp_cif(SRTIO3_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        expected = 5.0 / (59.55e-30)
        assert result['n_density'] == pytest.approx(expected, rel=0.01)
        print(f"SrTiO3 n_density: {result['n_density']:.2e} atoms/m³")

    def test_nacl_8_atoms(self):
        """NaCl: 8 atoms in conventional cell (4 Na + 4 Cl)."""
        path = _write_temp_cif(NACL_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        assert result['total_atoms'] == pytest.approx(8.0)
        assert result['volume'] == pytest.approx(179.43)

        # NaCl density ~ 2165 kg/m³
        assert 2000 < result['density'] < 2300, \
            f"NaCl density = {result['density']:.0f} kg/m³"
        print(f"NaCl density: {result['density']:.0f} kg/m³")

    def test_partial_occupancy(self):
        """Partial occupancy site: total atoms = sum of occupancies."""
        path = _write_temp_cif(PARTIAL_OCC_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        # 0.375 + 0.4375 + 0.25 + 0.75 + 1.0 + 1.0 + 1.0 = 4.8125
        expected_atoms = 0.375 + 0.4375 + 0.25 + 0.75 + 3.0
        assert result['total_atoms'] == pytest.approx(expected_atoms)

        print(f"LSHT-like total atoms: {result['total_atoms']:.4f}")
        print(f"LSHT-like density: {result['density']:.0f} kg/m³")

    def test_no_occupancy_column(self):
        """CIF without occupancy column defaults to 1.0."""
        path = _write_temp_cif(NO_OCC_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        assert result['total_atoms'] == pytest.approx(2.0)
        symbols = [a['symbol'] for a in result['atoms']]
        assert 'Fe' in symbols
        assert 'O' in symbols

    def test_volume_calculation(self):
        """Volume calculated from lattice params matches _cell_volume."""
        path = _write_temp_cif(SRTIO3_CIF)
        result = cif_reader.read(path, use_pymatgen=False)

        # Manually calculate cubic volume
        calc_vol = 3.905**3
        assert result['volume'] == pytest.approx(59.55, abs=0.1)
        assert calc_vol == pytest.approx(59.55, abs=0.1)

    def test_get_material_params(self):
        """get_material_params returns (V, ρ, n) tuple."""
        path = _write_temp_cif(SRTIO3_CIF)
        V, rho, n = cif_reader.get_material_params(path, use_pymatgen=False)

        assert V == pytest.approx(59.55e-30)
        assert 4900 < rho < 5300
        assert n > 0

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            cif_reader.read("/nonexistent/file.cif")

    def test_parser_label(self):
        """Verify parser field is set correctly."""
        path = _write_temp_cif(SRTIO3_CIF)
        result = cif_reader.read(path, use_pymatgen=False)
        assert result['parser'] == 'basic'


class TestCIFReaderIntegration:
    """Integration tests: CIF reader → Pipeline 1/2 inputs."""

    def test_lsht_like_feeds_pipeline(self):
        """Material params from CIF can feed into Pipeline 1 functions."""
        from MatSciKit_COSMOTIM.heat_capacity import dulong_petit

        path = _write_temp_cif(PARTIAL_OCC_CIF)
        V, rho, n = cif_reader.get_material_params(path, use_pymatgen=False)

        # Dulong-Petit should work with CIF-derived params
        dp = dulong_petit.calculate(n, rho)
        assert dp > 0
        print(f"LSHT-like Dulong-Petit: {dp:.4f} J/(g·K)")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
