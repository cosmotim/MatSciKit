"""
Tests for Pipeline 0: Material class and database loader.
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
            v_avg=3461.3,
        )
        assert mat.name == 'LSHT'
        assert mat.density > 0
        assert mat.theta_D is not None  # auto-estimated
        assert mat.theta_D > 0
        print(f"LSHT θ_D (estimated) = {mat.theta_D:.1f} K")

    def test_theta_D_auto_estimate(self):
        """θ_D should be auto-estimated when v_avg is provided."""
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28, v_avg=3000,
        )
        assert mat.theta_D is not None
        assert 200 < mat.theta_D < 800

    def test_theta_D_explicit(self):
        """Explicit θ_D should not be overwritten."""
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28, v_avg=3000,
            theta_D=500.0,
        )
        assert mat.theta_D == 500.0

    def test_no_v_avg_no_theta(self):
        """Without v_avg, θ_D stays None."""
        mat = Material(
            name='test', formula='test',
            density=5000, n_density=8e28,
        )
        assert mat.theta_D is None

    def test_from_params(self):
        mat = Material.from_params(
            name='LSHT',
            formula='Li3/8Sr7/16Hf1/4Ta3/4O3',
            volume_m3=LSHT_V,
            total_atoms=LSHT_N,
            total_mass_kg=LSHT_M,
            v_avg=3461.3,
        )
        assert mat.density == pytest.approx(LSHT_DENSITY, rel=1e-6)
        assert mat.n_density == pytest.approx(LSHT_N_DENSITY, rel=1e-6)

    def test_from_cif(self):
        """Create Material from a CIF file."""
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

        mat = Material.from_cif(f.name, name='STO', v_avg=5000)
        assert mat.name == 'STO'
        assert 4900 < mat.density < 5300
        assert mat.v_avg == 5000
        assert mat.theta_D is not None


class TestEstimateFunctions:

    def test_debye_temperature(self):
        theta = estimate_debye_temperature(3461.3, LSHT_N_DENSITY)
        # LSHT θ_D from Pipeline 1 fit = 437 K; estimate should be similar
        assert 350 < theta < 550, f"θ_D = {theta:.0f} K"
        print(f"Estimated θ_D = {theta:.0f} K (Pipeline 1 fit: 437 K)")

    def test_n_density(self):
        avg_mass = LSHT_M / LSHT_N / 1.66053906660e-27  # back to amu
        n = estimate_n_density(LSHT_DENSITY, avg_mass)
        assert n == pytest.approx(LSHT_N_DENSITY, rel=0.01)


class TestDatabaseLoader:

    def _make_test_db(self) -> str:
        """Create a minimal test database JSON."""
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

        # Should load 1 material (skip incomplete)
        assert len(materials) == 1
        mat = materials[0]
        assert mat.name == 'TestOxide'
        assert mat.v_avg == 2550
        assert mat.density == 5100
        assert mat.category == 'oxides'
        assert len(mat.sources) == 1

    def test_load_researcher_db(self):
        """Load the actual researcher database if available."""
        db_path = Path.home() / '.openclaw/agents/researcher/sound_velocity_database.json'
        if not db_path.exists():
            pytest.skip("Researcher database not available")

        materials = load_database(db_path)
        assert len(materials) >= 7
        names = [m.name for m in materials]
        assert 'LLZTO' in names
        assert 'LGPS' in names

        for mat in materials:
            assert mat.density > 0
            assert mat.v_avg > 0
            assert mat.theta_D is not None
            print(f"  {mat.name}: v_avg={mat.v_avg} m/s, θ_D={mat.theta_D:.0f} K")


class TestBatchCahill:

    def test_batch_at_300K(self):
        """Batch Cahill calculation at 300 K."""
        materials = [
            Material('A', 'A', density=5000, n_density=8e28, v_avg=3000),
            Material('B', 'B', density=3000, n_density=5e28, v_avg=2000),
        ]

        results = batch_cahill(materials, T=300.0)
        assert len(results) == 2
        assert results['A'] > results['B']  # higher v → higher κ_min
        assert all(v > 0 for v in results.values())

    def test_batch_array_T(self):
        """Batch with array temperature."""
        materials = [
            Material('X', 'X', density=5000, n_density=8e28, v_avg=3000),
        ]
        T = np.linspace(10, 300, 50)
        results = batch_cahill(materials, T=T)

        assert 'X' in results
        assert len(results['X']) == 50
        assert np.all(np.diff(results['X']) > 0)  # monotonic increase

    def test_batch_with_researcher_db(self):
        """Full integration test with researcher database."""
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
