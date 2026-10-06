"""Tests for .mat generator."""

import tempfile
from pathlib import Path

from hstar_yaml.schema import (
    HstarProject, MeshConfig, MaterialConfig, MaterialType,
    MaterialProperties, ElementGroupConfig, ElementType, Phase,
)
from hstar_yaml.generators.mat import MatGenerator


def _make_project(materials):
    return HstarProject(
        title="MAT Test",
        problem_name="test",
        mesh=MeshConfig(file="test.msh", dimension=2),
        materials=materials,
        element_groups=[
            ElementGroupConfig(group=1, material=1, element_type=ElementType.Q4),
        ],
    )


class TestMatGenerator:
    def test_elastic_material(self):
        project = _make_project([
            MaterialConfig(
                id=1, name="Steel", type=MaterialType.ELASTIC_ISOTROPIC,
                properties=MaterialProperties(E=210e9, nu=0.3, density=7800),
            )
        ])
        gen = MatGenerator(project)
        content = gen.build()
        assert "ELASTIC_ISOTROPIC" in content
        assert "MECHANICAL" in content
        assert "SOLID" in content

    def test_multiple_materials(self):
        project = _make_project([
            MaterialConfig(
                id=1, name="Soil", type=MaterialType.ELASTIC_ISOTROPIC,
                properties=MaterialProperties(E=1e9, nu=0.3),
            ),
            MaterialConfig(
                id=2, name="Water", type=MaterialType.FLUID, phase=Phase.FLUID,
                properties=MaterialProperties(density=1000, bulk_modulus=2.2e9),
            ),
        ])
        gen = MatGenerator(project)
        content = gen.build()
        assert "material_serial         1" in content
        assert "material_serial         2" in content
        assert "FLUID" in content

    def test_goodman_material(self):
        project = _make_project([
            MaterialConfig(
                id=1, name="Joint", type=MaterialType.GOODMAN, phase=Phase.CONTACT,
                properties=MaterialProperties(kn=1e8, ks=1e6, friction=0.3),
            )
        ])
        gen = MatGenerator(project)
        content = gen.build()
        assert "GOODMAN" in content
        assert "CONTACT" in content

    def test_writes_file(self):
        project = _make_project([
            MaterialConfig(
                id=1, type=MaterialType.ELASTIC_ISOTROPIC,
                properties=MaterialProperties(E=1e9, nu=0.3),
            )
        ])
        gen = MatGenerator(project)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = gen.generate(Path(tmpdir))
            assert path.name == "test.mat"
            content = path.read_text(encoding="utf-8")
            assert "material_serial" in content
