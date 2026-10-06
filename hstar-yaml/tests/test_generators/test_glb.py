"""Tests for .glb generator."""

import tempfile
from pathlib import Path

from hstar_yaml.schema import (
    HstarProject, MeshConfig, MaterialConfig, MaterialType,
    MaterialProperties, ElementGroupConfig, ElementType, DofType,
    AnalysisConfig, AnalysisType, SolverConfig, SolverType,
    OutputConfig,
)
from hstar_yaml.generators.glb import GlbGenerator
from hstar_yaml.mesh_reader import MeshData, MeshNode, MeshElement, MeshGroup


def _make_mesh():
    mesh = MeshData(npoin=4, nelem=1, ndimn=2)
    for i in range(1, 5):
        mesh.nodes[i] = MeshNode(id=i, x=float(i), y=0.0)
    mesh.elements[1] = MeshElement(id=1, nodes=[1, 2, 3, 4], group=1)
    mesh.groups[1] = MeshGroup(id=1, element_ids=[1], nelem=1)
    return mesh


def _make_project():
    return HstarProject(
        title="GLB Test",
        problem_name="test",
        mesh=MeshConfig(file="test.msh", dimension=2),
        materials=[
            MaterialConfig(id=1, type=MaterialType.ELASTIC_ISOTROPIC,
                           properties=MaterialProperties(E=1e9, nu=0.3)),
        ],
        element_groups=[
            ElementGroupConfig(group=1, material=1, element_type=ElementType.Q4),
        ],
        output=OutputConfig(results=["displacement", "stress"]),
    )


class TestGlbGenerator:
    def test_generates_content(self):
        project = _make_project()
        mesh = _make_mesh()
        gen = GlbGenerator(project, mesh)
        content = gen.build()
        assert content is not None
        assert "NPOIN" in content
        assert "PARDISO" in content

    def test_npoin_nelem(self):
        project = _make_project()
        mesh = _make_mesh()
        gen = GlbGenerator(project, mesh)
        content = gen.build()
        lines = content.split("\n")
        # Line 2 should contain npoin=4 and nelem=1
        assert "4" in lines[1]
        assert "1" in lines[1]

    def test_static_type(self):
        project = _make_project()
        mesh = _make_mesh()
        gen = GlbGenerator(project, mesh)
        content = gen.build()
        assert "F         PARDISO" in content

    def test_group_definition(self):
        project = _make_project()
        mesh = _make_mesh()
        gen = GlbGenerator(project, mesh)
        content = gen.build()
        assert "GROUP1" in content
        assert "Q4" in content

    def test_gid_output_flags(self):
        project = _make_project()
        mesh = _make_mesh()
        gen = GlbGenerator(project, mesh)
        content = gen.build()
        # displacement and stress should be enabled
        assert "gid_u" in content

    def test_write_to_file(self):
        project = _make_project()
        mesh = _make_mesh()
        gen = GlbGenerator(project, mesh)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = gen.generate(Path(tmpdir))
            assert path is not None
            assert path.name == "test.glb"
            assert path.read_text(encoding="utf-8").startswith("NPOIN")
