"""Tests for YAML schema validation."""

import pytest
from hstar_yaml.schema import (
    HstarProject, MeshConfig, MaterialConfig, MaterialType, Phase,
    MaterialProperties, ElementGroupConfig, ElementType, DofType,
    AnalysisConfig, AnalysisType, NonlinearMethod, BoundaryCondition,
    BcType, DofName, SolverConfig, SolverType, TimeCurve, TimeCurveType,
    LoadConfig, GravityLoad, OutputConfig,
)


def _minimal_project(**overrides):
    """Create a minimal valid project config."""
    defaults = dict(
        title="Test",
        problem_name="1",
        mesh=MeshConfig(file="test.msh", dimension=2),
        materials=[
            MaterialConfig(
                id=1,
                name="Soil",
                type=MaterialType.ELASTIC_ISOTROPIC,
                properties=MaterialProperties(E=1e9, nu=0.3),
            )
        ],
        element_groups=[
            ElementGroupConfig(group=1, material=1, element_type=ElementType.Q4),
        ],
    )
    defaults.update(overrides)
    return HstarProject(**defaults)


class TestMinimalProject:
    def test_create_minimal(self):
        p = _minimal_project()
        assert p.title == "Test"
        assert p.problem_name == "1"
        assert len(p.materials) == 1
        assert len(p.element_groups) == 1

    def test_default_analysis(self):
        p = _minimal_project()
        assert p.analysis.type == AnalysisType.STATIC
        assert p.analysis.nonlinear == NonlinearMethod.NEWTON_RAPHSON

    def test_default_solver(self):
        p = _minimal_project()
        assert p.solver.type == SolverType.PARDISO
        assert p.solver.num_threads == 4


class TestMaterials:
    def test_elastic(self):
        m = MaterialConfig(
            id=1, name="Steel", type=MaterialType.ELASTIC_ISOTROPIC,
            properties=MaterialProperties(E=210e9, nu=0.3, density=7800),
        )
        assert m.properties.E == 210e9
        assert m.phase == Phase.SOLID

    def test_mohr_coulomb(self):
        m = MaterialConfig(
            id=2, name="Clay", type=MaterialType.MOHR_COULOMB,
            properties=MaterialProperties(
                E=30e6, nu=0.35, cohesion=25e3,
                friction_angle=20.0, dilation_angle=0.0,
            ),
        )
        assert m.properties.cohesion == 25e3

    def test_fluid(self):
        m = MaterialConfig(
            id=3, name="Water", type=MaterialType.FLUID, phase=Phase.FLUID,
            properties=MaterialProperties(
                density=1000, bulk_modulus=2.2e9,
            ),
        )
        assert m.phase == Phase.FLUID

    def test_goodman(self):
        m = MaterialConfig(
            id=4, name="Joint", type=MaterialType.GOODMAN, phase=Phase.CONTACT,
            properties=MaterialProperties(kn=1e8, ks=1e6, friction=0.3),
        )
        assert m.properties.kn == 1e8


class TestBoundaryConditions:
    def test_fixed_bc(self):
        bc = BoundaryCondition(
            name="Bottom", type=BcType.FIXED,
            nodes=[1, 2, 3], dofs=[DofName.UX, DofName.UY], value=0.0,
        )
        assert len(bc.nodes) == 3
        assert len(bc.dofs) == 2

    def test_prescribed_bc(self):
        bc = BoundaryCondition(
            name="Top", type=BcType.PRESCRIBED,
            nodes=[10], dofs=[DofName.UY], value=-0.01, time_curve=1,
        )
        assert bc.time_curve == 1

    def test_must_have_nodes(self):
        with pytest.raises(ValueError, match="nodes.*node_range"):
            BoundaryCondition(
                name="Empty", dofs=[DofName.UX],
            )

    def test_node_range(self):
        bc = BoundaryCondition(
            name="Range", dofs=[DofName.UX],
            node_range={"start": 1, "end": 10, "step": 1},
        )
        assert bc.node_range["start"] == 1


class TestAnalysis:
    def test_dynamic_requires_dynamic_section(self):
        with pytest.raises(ValueError, match="Dynamic.*dynamic"):
            AnalysisConfig(type=AnalysisType.DYNAMIC)

    def test_static_ok_without_dynamic(self):
        a = AnalysisConfig(type=AnalysisType.STATIC)
        assert a.dynamic is None


class TestTimeCurves:
    def test_linear_curve(self):
        tc = TimeCurve(id=1, type=TimeCurveType.LINEAR, points=[[0, 0], [1, 1]])
        assert len(tc.points) == 2

    def test_seismic_curve(self):
        tc = TimeCurve(
            id=2, type=TimeCurveType.SEISMIC,
            file="eq.dat", direction="X", dt=0.01, total_time=30.0,
        )
        assert tc.file == "eq.dat"


class TestYamlParsing:
    def test_from_dict(self):
        """Test creating project from a dict (simulating YAML load)."""
        data = {
            "title": "Test Project",
            "problem_name": "test1",
            "mesh": {"file": "test.msh", "dimension": 2},
            "materials": [
                {
                    "id": 1, "name": "Soil",
                    "type": "ELASTIC_ISOTROPIC",
                    "properties": {"E": 1e9, "nu": 0.3},
                }
            ],
            "element_groups": [
                {"group": 1, "material": 1, "element_type": "Q4"},
            ],
            "analysis": {"type": "STATIC", "nonlinear": "NEWTON_RAPHSON"},
            "solver": {"type": "PARDISO", "matrix_type": -2},
        }
        p = HstarProject.model_validate(data)
        assert p.title == "Test Project"
        assert p.materials[0].type == MaterialType.ELASTIC_ISOTROPIC
        assert p.solver.matrix_type == -2
