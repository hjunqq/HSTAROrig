"""Tests for input validation."""

import pytest
from hstar_yaml.schema import (
    HstarProject, MeshConfig, MaterialConfig, MaterialType,
    MaterialProperties, ElementGroupConfig, ElementType, DofType,
    AnalysisConfig, AnalysisType, BoundaryCondition, DofName,
    Phase,
)
from hstar_yaml.validator import validate_project
from hstar_yaml.mesh_reader import MeshData, MeshNode


def _project(**overrides):
    defaults = dict(
        title="Test",
        problem_name="1",
        mesh=MeshConfig(file="test.msh", dimension=2),
        materials=[
            MaterialConfig(id=1, name="Soil", type=MaterialType.ELASTIC_ISOTROPIC,
                           properties=MaterialProperties(E=1e9, nu=0.3)),
        ],
        element_groups=[
            ElementGroupConfig(group=1, material=1, element_type=ElementType.Q4),
        ],
    )
    defaults.update(overrides)
    return HstarProject(**defaults)


def _mesh_with_nodes(*node_ids):
    mesh = MeshData(npoin=len(node_ids), ndimn=2)
    for nid in node_ids:
        mesh.nodes[nid] = MeshNode(id=nid, x=0.0, y=0.0)
    return mesh


class TestMaterialValidation:
    def test_valid_elastic(self):
        p = _project()
        errors, warnings = validate_project(p)
        assert len(errors) == 0

    def test_invalid_material_ref(self):
        p = _project(
            element_groups=[
                ElementGroupConfig(group=1, material=99, element_type=ElementType.Q4),
            ]
        )
        errors, _ = validate_project(p)
        assert any("material 99" in e for e in errors)

    def test_zero_E_elastic(self):
        p = _project(
            materials=[
                MaterialConfig(id=1, type=MaterialType.ELASTIC_ISOTROPIC,
                               properties=MaterialProperties(E=0, nu=0.3)),
            ]
        )
        errors, _ = validate_project(p)
        assert any("E must be > 0" in e for e in errors)

    def test_high_poisson_warning(self):
        p = _project(
            materials=[
                MaterialConfig(id=1, type=MaterialType.ELASTIC_ISOTROPIC,
                               properties=MaterialProperties(E=1e9, nu=0.499)),
            ]
        )
        _, warnings = validate_project(p)
        assert any("nu=" in w for w in warnings)

    def test_unused_material_warning(self):
        p = _project(
            materials=[
                MaterialConfig(id=1, type=MaterialType.ELASTIC_ISOTROPIC,
                               properties=MaterialProperties(E=1e9, nu=0.3)),
                MaterialConfig(id=2, name="Unused", type=MaterialType.ELASTIC_ISOTROPIC,
                               properties=MaterialProperties(E=1e9, nu=0.3)),
            ]
        )
        _, warnings = validate_project(p)
        assert any("Unused" in w or "material 2" in w.lower() for w in warnings)


class TestNodeValidation:
    def test_bc_node_not_in_mesh(self):
        mesh = _mesh_with_nodes(1, 2, 3)
        p = _project(
            boundary_conditions=[
                BoundaryCondition(name="Test", nodes=[999], dofs=[DofName.UX]),
            ]
        )
        errors, _ = validate_project(p, mesh)
        assert any("999" in e for e in errors)

    def test_bc_node_in_mesh(self):
        mesh = _mesh_with_nodes(1, 2, 3)
        p = _project(
            boundary_conditions=[
                BoundaryCondition(name="Test", nodes=[1, 2], dofs=[DofName.UX]),
            ]
        )
        errors, _ = validate_project(p, mesh)
        assert len(errors) == 0


class TestDimensionValidation:
    def test_3d_element_in_2d(self):
        p = _project(
            element_groups=[
                ElementGroupConfig(group=1, material=1, element_type=ElementType.B8),
            ]
        )
        errors, _ = validate_project(p)
        assert any("3D element" in e for e in errors)

    def test_2d_element_in_2d_ok(self):
        p = _project()
        errors, _ = validate_project(p)
        assert not any("element" in e.lower() and "dimension" in e.lower() for e in errors)
