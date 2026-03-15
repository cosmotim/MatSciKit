"""Tests for Material serialization (as_dict / from_dict / JSON)."""

from __future__ import annotations

import json

import pytest

from MatSciKit_COSMOTIM.structure.material import Material


class TestSerialization:
    """Test Material serialization round-trips."""

    def _make_material(self):
        return Material(
            name="LSHT",
            formula="Li₃/₈Sr₇/₁₆Hf₁/₄Ta₃/₄O₃",
            density=6870.0,
            n_density=7.63e28,
            sound_velocity={"transverse": 1900.0, "longitudinal": 3400.0},
            bulk_modulus={"vrh": 120.0},
        )

    def test_as_dict(self):
        mat = self._make_material()
        d = mat.as_dict()
        assert d["__class__"] == "Material"
        assert d["name"] == "LSHT"
        assert d["density"] == 6870.0
        assert d["sound_velocity"]["transverse"] == 1900.0

    def test_roundtrip_dict(self):
        mat = self._make_material()
        d = mat.as_dict()
        mat2 = Material.from_dict(d)
        assert mat2.name == mat.name
        assert mat2.density == mat.density
        assert mat2.sound_velocity == mat.sound_velocity
        assert mat2.bulk_modulus == mat.bulk_modulus

    def test_roundtrip_json(self, tmp_path):
        mat = self._make_material()
        json_path = tmp_path / "test_material.json"
        mat.to_json(json_path)

        # File should exist and be valid JSON
        assert json_path.exists()
        with open(json_path) as f:
            data = json.load(f)
        assert data["name"] == "LSHT"

        # Round-trip
        mat2 = Material.from_json(json_path)
        assert mat2.name == mat.name
        assert mat2.density == mat.density
        assert mat2.n_density == mat.n_density

    def test_dict_does_not_mutate(self):
        """Ensure from_dict doesn't fail with extra keys."""
        mat = self._make_material()
        d = mat.as_dict()
        d["extra_field"] = "should be ignored"
        # This should raise TypeError for unexpected kwarg
        with pytest.raises(TypeError):
            Material.from_dict(d)
