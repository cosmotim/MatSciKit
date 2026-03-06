"""
Tests for Pipeline 0: Material class, database loader, MP integration.
"""

import json
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from MatSciKit_COSMOTIM.structure.material import (
    Material, estimate_debye_temperature, estimate_n_density,
    load_database, batch_cahill,
)


# LSHT known values for validation
LSHT_V = 3.98e-10 ** 3  # m³
LSHT_N = 3/8 + 7/16 + 1/4 + 3/4 + 3
LSHT_M = (6.941*3/8 + 87.62*7/16 + 178.49*1/4 + 180.95*3/4 + 15.999*3) * 1.66053906660e-27
LSHT_DENSITY = LSHT_M / LSHT_V
LSHT_N_DENSITY = LSHT_N / LSHT_V


class TestMaterial:
    """Test the Material dataclass."""

    def test_basic_creation(self):
        mat = Material(
            name='LSHT',
            formula='Li3/8Sr7/16Hf1/4Ta3/4O3',
            density=LSHT_DENSITY,
            n_density=LSHT_N_DENSITY,
            sound_velocity={'transverse': 1900, 'longitudinal': 3400, 'average': 3461.3},
        )
        assert mat.name == 'LSHT'
        assert mat.density > 0
        assert mat.debye_temperature is not None
        assert mat.debye_temperature > 0
        assert mat.v_avg == pytest.approx(3461.3)
        assert mat.v_L == 3400
        assert mat.v_T == 1900
        print(f"LSHT θ_D (estimated) = {mat.debye_temperature:.1f} K")

    def test_v_avg_from_vl_vt(self):
        """v_avg computed from v_L and v_T when no average provided."""
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28,
            sound_velocity={'longitudinal': 5000, 'transverse': 3000},
        )
        assert mat.v_avg is not None
        assert mat.v_avg > 0
        # Debye average: (1/3 * (1/5000³ + 2/3000³))^(-1/3)
        expected = (1.0/3 * (1/5000**3 + 2/3000**3)) ** (-1.0/3)
        assert mat.v_avg == pytest.approx(expected, rel=1e-6)

    def test_v_avg_debye_over_snyder(self):
        """Debye average from v_L/v_T takes priority over snyder."""
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28,
            sound_velocity={'longitudinal': 5000, 'transverse': 3000,
                            'snyder_acoustic': 178},
        )
        # Should compute Debye average, not use snyder
        expected = (1.0/3 * (1/5000**3 + 2/3000**3)) ** (-1.0/3)
        assert mat.v_avg == pytest.approx(expected, rel=1e-6)
        assert mat.v_avg > 3000  # Not 178 (snyder)

    def test_theta_D_auto_estimate(self):
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28,
            sound_velocity={'average': 3000},
        )
        assert mat.debye_temperature is not None
        assert 200 < mat.debye_temperature < 800

    def test_theta_D_explicit(self):
        """Explicit θ_D should not be overwritten."""
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28,
            sound_velocity={'average': 3000},
            debye_temperature=500.0,
        )
        assert mat.debye_temperature == 500.0

    def test_no_sound_velocity(self):
        """Without sound_velocity, v_avg and θ_D are None."""
        mat = Material(name='test', formula='test', density=5000, n_density=8e28)
        assert mat.v_avg is None
        assert mat.v_T is None
        assert mat.debye_temperature is None

    def test_moduli_vrh(self):
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28,
            bulk_modulus={'voigt': 120, 'reuss': 110, 'vrh': 115},
            shear_modulus={'voigt': 50, 'reuss': 45, 'vrh': 47.5},
        )
        assert mat.K_vrh == 115
        assert mat.G_vrh == 47.5

    def test_n_density_from_density_atomic(self):
        """n_density derived from density_atomic."""
        mat = Material(
            name='test', formula='test',
            density=5000, density_atomic=0.084,  # atoms/ų
        )
        assert mat.n_density == pytest.approx(0.084 * 1e30)

    def test_from_params(self):
        mat = Material.from_params(
            name='LSHT',
            formula='Li3/8Sr7/16Hf1/4Ta3/4O3',
            volume_m3=LSHT_V,
            total_atoms=LSHT_N,
            total_mass_kg=LSHT_M,
            sound_velocity={'average': 3461.3},
        )
        assert mat.density == pytest.approx(LSHT_DENSITY, rel=1e-6)
        assert mat.n_density == pytest.approx(LSHT_N_DENSITY, rel=1e-6)
        assert mat.volume == pytest.approx(LSHT_V * 1e30)  # stored in ų

    def test_from_cif(self):
        cif_content = textwrap.dedent("""\
            data_SrTiO3
            _cell_length_a 3.905
            _cell_length_b 3.905
            _cell_length_c 3.905
            _cell_angle_alpha 90.0
            _cell_angle_beta 90.0
            _cell_angle_gamma 90.0
            _cell_volume 59.55

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
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False)
        f.write(cif_content)
        f.close()

        mat = Material.from_cif(
            f.name, name='STO',
            sound_velocity={'longitudinal': 7800, 'transverse': 4500},
        )
        assert mat.name == 'STO'
        assert 4900 < mat.density < 5300
        assert mat.v_avg is not None
        assert mat.debye_temperature is not None
        assert mat.nsites == 5
        assert mat.volume == pytest.approx(59.55, rel=0.01)


class TestMPIntegration:
    """Test Materials Project integration (mocked)."""

    def test_from_mp_mock(self):
        """Test from_mp with a mock document object."""
        class MockSV:
            transverse = 4500.0
            longitudinal = 7800.0
            snyder_acoustic = 5100.0
            snyder_optical = None
            snyder_total = None

        class MockBulk:
            voigt = 180.0
            reuss = 170.0
            vrh = 175.0

        class MockShear:
            voigt = 110.0
            reuss = 105.0
            vrh = 107.5

        class MockTC:
            clarke = 3.5
            cahill = 2.8

        class MockSym:
            crystal_system = 'Cubic'
            symbol = 'Pm-3m'
            number = 221

        class MockDoc:
            material_id = 'mp-5827'
            formula_pretty = 'SrTiO3'
            density = 5.117  # g/cm³
            density_atomic = 0.084
            volume = 59.55
            nsites = 5
            elements = ['Sr', 'Ti', 'O']
            sound_velocity = MockSV()
            bulk_modulus = MockBulk()
            shear_modulus = MockShear()
            young_modulus = 270.0
            homogeneous_poisson = 0.25
            universal_anisotropy = 0.0
            debye_temperature = 520.0
            thermal_conductivity = MockTC()
            symmetry = MockSym()

        mat = Material.from_mp(MockDoc(), name='STO')

        assert mat.name == 'STO'
        assert mat.material_id == 'mp-5827'
        assert mat.formula == 'SrTiO3'
        assert mat.density == pytest.approx(5117, rel=0.01)  # converted to kg/m³
        # v_avg = Debye average from v_L/v_T, not snyder_acoustic
        expected_v_avg = (1.0/3 * (1/7800**3 + 2/4500**3)) ** (-1.0/3)
        assert mat.v_avg == pytest.approx(expected_v_avg, rel=1e-6)
        assert mat.v_L == 7800.0
        assert mat.v_T == 4500.0
        assert mat.K_vrh == 175.0
        assert mat.G_vrh == 107.5
        assert mat.debye_temperature == 520.0
        assert mat.thermal_conductivity['cahill'] == 2.8
        assert mat.symmetry['crystal_system'] == 'Cubic'
        assert 'Sr' in mat.elements

        print(f"MP STO: θ_D={mat.debye_temperature} K, "
              f"v_avg={mat.v_avg} m/s, K={mat.K_vrh} GPa")


class TestEstimateFunctions:

    def test_debye_temperature(self):
        theta = estimate_debye_temperature(3461.3, LSHT_N_DENSITY)
        assert 350 < theta < 550, f"θ_D = {theta:.0f} K"
        print(f"Estimated θ_D = {theta:.0f} K (Pipeline 1 fit: 437 K)")

    def test_n_density(self):
        avg_mass = LSHT_M / LSHT_N / 1.66053906660e-27
        n = estimate_n_density(LSHT_DENSITY, avg_mass)
        assert n == pytest.approx(LSHT_N_DENSITY, rel=0.01)


class TestDatabaseLoader:

    def _make_test_db(self) -> str:
        db = {
            "materials": {
                "oxides": {
                    "TestOxide": {
                        "full_name": "Li3La3Zr2O12",
                        "sound_velocities": {
                            "longitudinal_m_per_s": 4100,
                            "transverse_m_per_s": 2300,
                            "average_m_per_s": 2550,
                        },
                        "elastic_moduli": {
                            "bulk_modulus_GPa": 110,
                            "shear_modulus_GPa": 55,
                        },
                        "density_kg_per_m3": 5100,
                        "sources": [{"citation": "Test et al. 2026"}],
                    },
                    "Incomplete": {
                        "full_name": "MissingData",
                        "sound_velocities": {"average_m_per_s": None},
                        "density_kg_per_m3": None,
                    },
                },
            }
        }
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False)
        json.dump(db, f)
        f.close()
        return f.name

    def test_load_database(self):
        path = self._make_test_db()
        materials = load_database(path)

        assert len(materials) == 1
        mat = materials[0]
        assert mat.name == 'TestOxide'
        assert mat.v_avg == 2550
        assert mat.v_L == 4100
        assert mat.v_T == 2300
        assert mat.density == 5100
        assert mat.category == 'oxides'
        assert mat.K_vrh == 110
        assert mat.G_vrh == 55

    def test_load_researcher_db(self):
        db_path = Path.home() / '.openclaw/agents/researcher/sound_velocity_database.json'
        if not db_path.exists():
            pytest.skip("Researcher database not available")

        materials = load_database(db_path)
        assert len(materials) >= 7
        names = [m.name for m in materials]
        assert 'LLZTO' in names

        for mat in materials:
            assert mat.density > 0
            assert mat.v_avg > 0
            assert mat.debye_temperature is not None
            print(f"  {mat.name}: v_avg={mat.v_avg} m/s, θ_D={mat.debye_temperature:.0f} K")


class TestBatchCahill:

    def test_batch_at_300K(self):
        materials = [
            Material('A', 'A', density=5000, n_density=8e28,
                     sound_velocity={'average': 3000}),
            Material('B', 'B', density=3000, n_density=5e28,
                     sound_velocity={'average': 2000}),
        ]
        results = batch_cahill(materials, T=300.0)
        assert len(results) == 2
        assert results['A'] > results['B']
        assert all(v > 0 for v in results.values())

    def test_batch_array_T(self):
        materials = [
            Material('X', 'X', density=5000, n_density=8e28,
                     sound_velocity={'average': 3000}),
        ]
        T = np.linspace(10, 300, 50)
        results = batch_cahill(materials, T=T)
        assert len(results['X']) == 50
        assert np.all(np.diff(results['X']) > 0)

    def test_batch_with_researcher_db(self):
        db_path = Path.home() / '.openclaw/agents/researcher/sound_velocity_database.json'
        if not db_path.exists():
            pytest.skip("Researcher database not available")

        materials = load_database(db_path)
        results = batch_cahill(materials, T=300.0)
        assert len(results) >= 7
        for name, kappa in sorted(results.items(), key=lambda x: -x[1]):
            print(f"  {name}: κ_min(300K) = {kappa:.3f} W/(m·K)")
            assert 0.1 < kappa < 3.0


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
