"""Pydantic models for HSTAR YAML input schema."""

from __future__ import annotations

from enum import Enum
from typing import Literal

from pydantic import BaseModel, Field, model_validator


# ── Enums ──────────────────────────────────────────────────────────────────

class MaterialType(str, Enum):
    ELASTIC_ISOTROPIC = "ELASTIC_ISOTROPIC"
    MOHR_COULOMB = "MOHR_COULOMB"
    DRUCKER_PRAGER = "DRUCKER_PRAGER"
    DUNCAN_CHANG = "DUNCAN_CHANG"
    CAMCLAY = "CAMCLAY"
    GOODMAN = "GOODMAN"
    CONCRETE = "CONCRETE"
    VON_MISES = "VON_MISES"
    TRESCA = "TRESCA"
    FLUID = "FLUID"
    CREEP = "CREEP"


class Phase(str, Enum):
    SOLID = "SOLID"
    FLUID = "FLUID"
    AIR = "AIR"
    CONTACT = "CONTACT"


class ElementType(str, Enum):
    L2 = "L2"
    L3 = "L3"
    T3 = "T3"
    T6 = "T6"
    Q4 = "Q4"
    Q8 = "Q8"
    H4 = "H4"
    H10 = "H10"
    B8 = "B8"
    B20 = "B20"
    BEAM = "BEAM"
    CONTACT = "CONTACT"


class DofType(str, Enum):
    U = "U"         # displacement only
    UW = "UW"       # displacement + pore pressure (Biot)
    UPC = "UPC"     # u-p coupled
    T = "T"         # temperature
    W = "W"         # pore pressure only


class AnalysisType(str, Enum):
    STATIC = "STATIC"
    DYNAMIC = "DYNAMIC"
    FREQUENCY = "FREQUENCY"
    THERMAL = "THERMAL"
    SEEPAGE = "SEEPAGE"


class NonlinearMethod(str, Enum):
    NEWTON_RAPHSON = "NEWTON_RAPHSON"
    MODIFIED_NR = "MODIFIED_NR"
    LINEAR = "LINEAR"


class SolverType(str, Enum):
    PARDISO = "PARDISO"
    PROFILE = "PROFILE"
    JPCG = "JPCG"
    PBCG = "PBCG"


class DampingType(str, Enum):
    RAYLEIGH = "RAYLEIGH"
    NONE = "NONE"


class TimeCurveType(str, Enum):
    LINEAR = "LINEAR"
    SEISMIC = "SEISMIC"
    TABULAR = "TABULAR"


class BcType(str, Enum):
    FIXED = "FIXED"
    PRESCRIBED = "PRESCRIBED"
    SPRING = "SPRING"


class DofName(str, Enum):
    UX = "UX"
    UY = "UY"
    UZ = "UZ"
    RX = "RX"
    RY = "RY"
    RZ = "RZ"
    PW = "PW"     # pore water pressure
    T = "T"       # temperature


class OutputFormat(str, Enum):
    GID = "GID"
    COSMOS = "COSMOS"


class OutputMode(str, Enum):
    REPLACE = "REPLACE"
    APPEND = "APPEND"


class BackAnalysisMethod(str, Enum):
    TRUST_REGION = "TRUST_REGION"


# ── Mesh ───────────────────────────────────────────────────────────────────

class MeshConfig(BaseModel):
    file: str = Field(default="", description="Path to GiD .msh file")
    cor_file: str = Field(default="", description="Path to .cor coordinate file (alternative to .msh)")
    ele_file: str = Field(default="", description="Path to .ele element file (used with cor_file)")
    dimension: Literal[2, 3] = Field(description="2D or 3D problem")


# ── Materials ──────────────────────────────────────────────────────────────

class MaterialProperties(BaseModel):
    """Material properties — different fields used depending on material type."""
    # Common
    density: float = 0.0
    E: float = 0.0
    nu: float = 0.0
    thickness: float = 1.0
    thermal_expansion: float = 0.0

    # Plasticity (MC/DP)
    cohesion: float = 0.0
    friction_angle: float = 0.0
    dilation_angle: float = 0.0
    tensile_strength: float = 0.0

    # Hardening
    hardening: float = 0.0

    # Goodman joint
    kn: float = 0.0               # normal stiffness
    ks: float = 0.0               # shear stiffness
    friction: float = 0.0
    janbu_modulus: float = 0.0
    janbu_exponent: float = 0.0
    janbu_k0: float = 0.0
    janbu_phi: float = 0.0

    # Duncan-Chang
    Rf: float = 0.0
    K_dc: float = 0.0
    n_dc: float = 0.0
    pa: float = 101.325e3         # atmospheric pressure

    # Cam-Clay
    Pc: float = 0.0
    lamda: float = 0.0
    kappa: float = 0.0
    Mg: float = 0.0
    Mf: float = 0.0
    D0: float = 0.0
    D1: float = 0.0

    # Fluid
    bulk_modulus: float = 0.0
    permeability: list[float] = Field(default_factory=list)

    # Creep
    creep_A: float = 0.0
    creep_n: float = 0.0
    creep_m: float = 0.0


class MaterialConfig(BaseModel):
    id: int = Field(ge=1, description="Material ID (1-based)")
    name: str = ""
    type: MaterialType
    phase: Phase = Phase.SOLID
    properties: MaterialProperties = Field(default_factory=MaterialProperties)


# ── Element Groups ─────────────────────────────────────────────────────────

class ElementGroupConfig(BaseModel):
    group: int = Field(ge=1, description="Group number (1-based)")
    material: int = Field(ge=1, description="Material ID reference")
    element_type: ElementType = ElementType.Q4
    gauss_order: int = Field(default=2, ge=1, le=5)
    dof_type: DofType = DofType.U
    # Additional group params from .glb
    element_class: str = "CO"       # CO = continuum, ST = structural
    algo_type: str = "ST"           # ST = standard
    stiffness_type: str = "PE"      # PE = standard
    mass_type: int = 0              # 0 = lumped, 1 = consistent
    rayleigh_alpha: float = 0.0
    rayleigh_beta: float = 0.0
    elcod_local: float = 0.0       # local coordinate system rotation


# ── Boundary Conditions ───────────────────────────────────────────────────

class BoundaryCondition(BaseModel):
    name: str = ""
    type: BcType = BcType.FIXED
    nodes: list[int] = Field(default_factory=list)
    node_range: dict | None = None  # {start, end, step}
    dofs: list[DofName]
    value: float = 0.0
    time_curve: int | None = None

    @model_validator(mode="after")
    def check_nodes_or_range(self):
        if not self.nodes and self.node_range is None:
            raise ValueError("Must specify either 'nodes' or 'node_range'")
        return self


# ── Loads ──────────────────────────────────────────────────────────────────

class GravityLoad(BaseModel):
    enabled: bool = True
    direction: list[float] = Field(default_factory=lambda: [0, -1])
    magnitude: float = 9.81


class PointLoad(BaseModel):
    nodes: list[int]
    values: list[float]       # Fx, Fy [, Fz]
    time_curve: int | None = None


class DistributedLoad(BaseModel):
    elements: list[int]
    face: int = 1
    pressure: float = 0.0
    time_curve: int | None = None


class LoadConfig(BaseModel):
    gravity: GravityLoad | None = None
    point_loads: list[PointLoad] = Field(default_factory=list)
    distributed_loads: list[DistributedLoad] = Field(default_factory=list)


# ── Time Curves ────────────────────────────────────────────────────────────

class TimeCurve(BaseModel):
    id: int = Field(ge=1)
    type: TimeCurveType = TimeCurveType.LINEAR
    points: list[list[float]] = Field(default_factory=list)  # [[t, val], ...]
    # For SEISMIC type
    file: str | None = None
    scale: float = 1.0
    direction: str = "X"          # X, Y, Z
    dt: float = 0.01
    total_time: float = 0.0
    nlines: int = 0               # number of data lines in seismic file


# ── Analysis ───────────────────────────────────────────────────────────────

class ToleranceConfig(BaseModel):
    force: float = 1.0e-3
    displacement: float = 1.0e-3


class IncrementBlock(BaseModel):
    """One increment block (nincs group)."""
    steps: int = 10
    dt: float = 1.0
    stiffness_update: int = 1     # istif
    load_curve: int = 1           # icurve
    max_iterations: int = 50
    print_interval: int = 1       # iprin
    output_interval: int = 1      # outinc
    tolerance: ToleranceConfig = Field(default_factory=ToleranceConfig)


class IncrementConfig(BaseModel):
    total_steps: int = 1
    max_iterations: int = 50
    tolerance: ToleranceConfig = Field(default_factory=ToleranceConfig)
    dt: float = 1.0               # pseudo time step for static
    blocks: list[IncrementBlock] = Field(default_factory=list,
        description="Multiple increment blocks (for multi-stage analysis)")


class DampingConfig(BaseModel):
    type: DampingType = DampingType.NONE
    alpha: float = 0.0
    beta: float = 0.0


class DynamicConfig(BaseModel):
    method: str = "NEWMARK"
    beta1: float = 0.25
    beta2: float = 0.50
    theta1: float = 1.0
    total_time: float = 1.0
    dt: float = 0.01
    damping: DampingConfig = Field(default_factory=DampingConfig)


class AnalysisConfig(BaseModel):
    type: AnalysisType = AnalysisType.STATIC
    nonlinear: NonlinearMethod = NonlinearMethod.NEWTON_RAPHSON
    large_deformation: bool = False

    increments: IncrementConfig = Field(default_factory=IncrementConfig)
    dynamic: DynamicConfig | None = None

    @model_validator(mode="after")
    def check_dynamic_config(self):
        if self.type == AnalysisType.DYNAMIC and self.dynamic is None:
            raise ValueError("Dynamic analysis requires 'dynamic' section")
        return self


# ── Solver ─────────────────────────────────────────────────────────────────

class SolverConfig(BaseModel):
    type: SolverType = SolverType.PARDISO
    matrix_type: int = Field(default=-2, description="1=SPD, 2=symmetric indefinite, -2=symmetric indefinite (PARDISO)")
    num_threads: int = 4
    message_level: int = 0
    # PARDISO advanced params (when isdefault=1)
    isdefault: int = Field(default=1, description="0=use PARDISO defaults, 1=use custom params below")
    reducing_order: int = 2
    PreCGS: int = 0
    permutation: int = 0
    maxiter: int = 20
    out_of_core: int = 0
    eps_pivot: int = 20
    iparm11: int = 0
    iparm13: int = 0


# ── Construction Stages ────────────────────────────────────────────────────

class MaterialChange(BaseModel):
    group: int
    new_material: int


class ConstructionStage(BaseModel):
    stage: int = Field(ge=1)
    name: str = ""
    active_groups: list[int] = Field(default_factory=list)
    deactivate_groups: list[int] = Field(default_factory=list)
    activate_groups: list[int] = Field(default_factory=list)
    steps: list[int] = Field(min_length=2, max_length=2)  # [start, end]
    material_change: list[MaterialChange] = Field(default_factory=list)


# ── Contact ────────────────────────────────────────────────────────────────

class ContactGapGroup(BaseModel):
    id: int = Field(ge=1)
    name: str = ""
    master_nodes: list[int]
    slave_nodes: list[int]
    initial_gap: float = 0.0
    friction: float = 0.0
    cohesion: float = 0.0
    normal_stiffness: float = 1.0e8
    tangential_stiffness: float = 1.0e6


class ContactConfig(BaseModel):
    gap_groups: list[ContactGapGroup] = Field(default_factory=list)
    penalty_stiffness: float = 1.0
    max_iterations: int = 500
    tolerance: float = 1.0e-5
    solver_type: SolverType = SolverType.PROFILE


# ── Output ─────────────────────────────────────────────────────────────────

class OutputConfig(BaseModel):
    format: OutputFormat = OutputFormat.GID
    mode: OutputMode = OutputMode.REPLACE
    results: list[str] = Field(
        default_factory=lambda: ["displacement", "stress"]
    )
    monitor_nodes: list[int] = Field(default_factory=list)
    monitor_interval: int = 1


# ── Back Analysis ──────────────────────────────────────────────────────────

class BackAnalysisParameter(BaseModel):
    material: int
    property: str
    initial: float
    bounds: list[float] = Field(min_length=2, max_length=2)


class ObservationPoint(BaseModel):
    node: int
    dof: DofName
    weight: float = 1.0


class BackAnalysisConfig(BaseModel):
    enabled: bool = False
    method: BackAnalysisMethod = BackAnalysisMethod.TRUST_REGION
    parameters: list[BackAnalysisParameter] = Field(default_factory=list)
    observations: dict | None = None


# ── Surface Force (FTR) ───────────────────────────────────────────────────

class SurfaceForce(BaseModel):
    group: int
    face_nodes: int = 2           # number of nodes per face
    elements: list[int]
    appear_type: int = 3          # 1=safety, 2=force, 3=both


class SurfaceForceConfig(BaseModel):
    forces: list[SurfaceForce] = Field(default_factory=list)
    safety_gaps: int = 0


# ── Water / IFS ────────────────────────────────────────────────────────────

class WaterConfig(BaseModel):
    """Water level and interface settings (IFS)."""
    enabled: bool = False
    add_mass: int = 1             # Icaddmass
    surface_id: int = 0           # swlifs2006
    total_height: float = 0.0     # toth
    ifswater: int = 0
    gravity: float = 9.8
    absorb: float = 0.0
    alfa_p4: float = 0.0
    stiff_p4: float = 0.0
    dam_height: float = 0.0
    water_level: float = -99.0


# ── Global Control Params ─────────────────────────────────────────────────

class GlbControl(BaseModel):
    """Advanced .glb parameters that most users won't change."""
    ninit: int = 0                # initial stress computation steps
    kinit: int = 0                # initial stress type (0=none, 1=K0)
    winit: int = -1               # initial water flag (-1=no, 1=yes)
    kstab: int = 0                # stabilization method
    kstat: int = 0                # static condensation
    nmass: int = 999              # mass matrix control (999=auto)
    nsmat: int = 999              # stiffness matrix control
    nhmat: int = 999              # h matrix control
    nqmat: int = 999              # q matrix control
    nswkw: int = 999              # water stiffness
    ngrav: int = 999              # gravity (999=auto, 0=off, 1=on)
    ntsmat: int = 999
    nthmat: int = 999
    # IFS/water params
    icaddmass: int = 0
    swlifs2006: int = 0
    toth: float = 0.0
    ifswater: int = 0
    ifsgravity: float = 9.8
    absorb: float = 0.0
    alfa_p4: float = 0.0
    stiff_p4: float = 1.0e20
    # Crack/steel params
    ftcrack: float = 1.5e6
    coefMpa: float = 1.0e6
    ikindks: int = 0
    doubsig: int = 2
    ktan1: float = 1.0e8
    ktan2: float = 1.0e8
    # Newmark params (override analysis.dynamic defaults)
    beta1: float | None = None
    beta2: float | None = None
    theta1: float | None = None


# ── Root Model ─────────────────────────────────────────────────────────────

class HstarProject(BaseModel):
    title: str = "HSTAR Analysis"
    problem_name: str = "1"

    mesh: MeshConfig
    materials: list[MaterialConfig]
    element_groups: list[ElementGroupConfig]

    boundary_conditions: list[BoundaryCondition] = Field(default_factory=list)
    loads: LoadConfig = Field(default_factory=LoadConfig)
    time_curves: list[TimeCurve] = Field(default_factory=list)

    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    solver: SolverConfig = Field(default_factory=SolverConfig)

    construction_stages: list[ConstructionStage] = Field(default_factory=list)
    contact: ContactConfig | None = None
    output: OutputConfig = Field(default_factory=OutputConfig)
    back_analysis: BackAnalysisConfig | None = None
    surface_forces: SurfaceForceConfig | None = None
    water: WaterConfig | None = None
    glb_control: GlbControl = Field(default_factory=GlbControl)
