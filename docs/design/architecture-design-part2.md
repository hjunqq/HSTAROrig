# HSTAR-Next 架构设计文档 (续)
# Architecture Design Document (Continued)

---

## 7. 输入输出规范

### 7.1 输入文件格式

#### 7.1.1 主控制文件 (*.hstar)

```yaml
# HSTAR-Next 控制文件
# 采用类YAML格式，易于解析和阅读

*HEADER
  VERSION: 2.0
  TITLE: "Dam Stress Analysis"
  DATE: 2024-01-15

*PROBLEM
  TYPE: STATIC                  # STATIC | DYNAMIC | THERMAL | COUPLED_HM
  DIMENSION: 3                  # 2 | 3
  LARGE_DEFORMATION: FALSE      # TRUE | FALSE

*FILES
  MESH: mesh/dam.msh            # 网格文件
  MATERIALS: materials.mat      # 材料文件
  BOUNDARY: boundary.bnd        # 边界条件
  LOADS: loads.lod              # 荷载文件
  OUTPUT_DIR: results/          # 输出目录

*SOLVER
  TYPE: PARDISO                 # PARDISO | CG | GMRES
  TOLERANCE: 1.0E-10
  MAX_ITERATIONS: 10000

*NONLINEAR
  METHOD: NEWTON_RAPHSON        # NEWTON_RAPHSON | MODIFIED_NR | BFGS
  MAX_ITERATIONS: 25
  FORCE_TOLERANCE: 1.0E-6
  DISPLACEMENT_TOLERANCE: 1.0E-8
  LINE_SEARCH: TRUE

*TIME  # 仅动力分析
  INTEGRATION: NEWMARK          # NEWMARK | HHT | EXPLICIT
  TOTAL_TIME: 10.0
  TIME_STEP: 0.01
  BETA: 0.25
  GAMMA: 0.5

*OUTPUT
  FORMAT: GID                   # GID | VTK | BOTH
  FREQUENCY: 1                  # 每N步输出
  VARIABLES:
    - DISPLACEMENT
    - STRESS
    - STRAIN
    - PLASTIC_STRAIN

*LOAD_BLOCKS
  INCLUDE: loadblocks.ldb       # 或直接定义

*END
```

#### 7.1.2 网格文件格式 (*.msh)

```
# HSTAR-Next Mesh File
# Version 2.0

*NODES, COUNT=1000
# ID, X, Y, Z
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
...

*ELEMENTS, COUNT=500, TYPE=HEX8
# ID, N1, N2, N3, N4, N5, N6, N7, N8, MATERIAL_ID, GROUP_ID
1, 1, 2, 3, 4, 5, 6, 7, 8, 1, 1
2, 5, 6, 7, 8, 9, 10, 11, 12, 1, 1
...

*ELEMENTS, COUNT=100, TYPE=TET4
# ID, N1, N2, N3, N4, MATERIAL_ID, GROUP_ID
501, 100, 101, 102, 103, 2, 2
...

*NODE_SETS
SET, NAME=bottom
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10
END_SET

SET, NAME=top
  991, 992, 993, 994, 995, 996, 997, 998, 999, 1000
END_SET

*ELEMENT_SETS
SET, NAME=concrete
  1, 2, 3, 4, 5, ..., 400
END_SET

SET, NAME=rock
  401, 402, ..., 600
END_SET

*END
```

#### 7.1.3 材料文件格式 (*.mat)

```
# HSTAR-Next Materials File

*MATERIAL, ID=1, NAME="Concrete"
  TYPE: ELASTIC
  E: 30.0E9
  NU: 0.2
  DENSITY: 2400.0
  THERMAL_EXPANSION: 1.0E-5

*MATERIAL, ID=2, NAME="Rock_Foundation"
  TYPE: MOHR_COULOMB
  E: 20.0E9
  NU: 0.25
  DENSITY: 2600.0
  COHESION: 1.0E6
  FRICTION_ANGLE: 35.0      # degrees
  DILATION_ANGLE: 10.0      # degrees
  TENSILE_STRENGTH: 0.5E6

*MATERIAL, ID=3, NAME="Joint"
  TYPE: GOODMAN
  KN: 1.0E10               # 法向刚度
  KS: 1.0E9                # 切向刚度
  COHESION: 0.1E6
  FRICTION_ANGLE: 30.0
  TENSILE_STRENGTH: 0.0

*MATERIAL, ID=4, NAME="Concrete_Creep"
  TYPE: CREEP_POWER
  E: 30.0E9
  NU: 0.2
  DENSITY: 2400.0
  # 蠕变参数
  CREEP_COEFFICIENT: 2.5
  CREEP_EXPONENT: 0.3
  REFERENCE_TIME: 28.0     # days

*END
```

#### 7.1.4 边界条件文件 (*.bnd)

```
# HSTAR-Next Boundary Conditions File

*DISPLACEMENT_BC
  # 固定底面
  SET: bottom
  DOF: ALL                  # ALL | X | Y | Z | XY | XZ | YZ
  VALUE: 0.0

  # 对称面
  SET: symmetry_plane
  DOF: X
  VALUE: 0.0

*PRESCRIBED_DISPLACEMENT
  # 指定位移
  SET: top_center
  DOF: Z
  VALUE: -0.01              # 或 CURVE=1

*SPRING_BC
  SET: elastic_support
  DOF: Z
  STIFFNESS: 1.0E8

*CONTACT
  MASTER_SET: dam_base
  SLAVE_SET: foundation_top
  TYPE: PENALTY             # PENALTY | LAGRANGE | AUGMENTED
  NORMAL_STIFFNESS: 1.0E12
  FRICTION: 0.6

*END
```

#### 7.1.5 荷载文件 (*.lod)

```
# HSTAR-Next Load File

*CURVES
  # 时间-荷载曲线
  CURVE, ID=1, NAME="ramp"
    0.0, 0.0
    1.0, 1.0
    10.0, 1.0

  CURVE, ID=2, NAME="earthquake"
    INCLUDE: earthquake_acc.dat

*LOAD_BLOCK, ID=1, NAME="Gravity"
  INCREMENTS: 10
  *BODY_FORCE
    SET: ALL
    TYPE: GRAVITY
    DIRECTION: 0.0, 0.0, -9.81
    CURVE: 1

*LOAD_BLOCK, ID=2, NAME="Water_Pressure"
  INCREMENTS: 5
  AFTER_BLOCK: 1
  *PRESSURE
    SET: upstream_face
    VALUE: HYDROSTATIC
    WATER_LEVEL: 150.0
    WATER_DENSITY: 1000.0

*LOAD_BLOCK, ID=3, NAME="Earthquake"
  INCREMENTS: 1000
  AFTER_BLOCK: 2
  TIME_STEP: 0.01
  *ACCELERATION
    DIRECTION: X
    CURVE: 2
    SCALE: 1.0

*END
```

### 7.2 输出格式

#### 7.2.1 GiD结果格式 (*.res)

```
GiD Post Results File 1.0

Result "Displacement" "Analysis" 1 Vector OnNodes
ComponentNames "Ux" "Uy" "Uz"
Values
1  0.000000E+00  0.000000E+00  -1.234567E-03
2  1.234567E-04  0.000000E+00  -1.345678E-03
...
End Values

Result "Stress" "Analysis" 1 Matrix OnGaussPoints "HEX8_GP"
ComponentNames "Sxx" "Syy" "Szz" "Sxy" "Syz" "Sxz"
Values
1  1.234E+06  2.345E+06  3.456E+06  1.234E+05  2.345E+05  3.456E+05
   1.234E+06  2.345E+06  3.456E+06  1.234E+05  2.345E+05  3.456E+05
   ...
End Values
```

#### 7.2.2 VTK格式 (*.vtu)

```xml
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="1000" NumberOfCells="500">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0.0 0.0 0.0
          1.0 0.0 0.0
          ...
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 1 2 3 4 5 6 7
          ...
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          8 16 24 ...
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          12 12 12 ...
        </DataArray>
      </Cells>
      <PointData Vectors="Displacement">
        <DataArray type="Float64" Name="Displacement" NumberOfComponents="3">
          0.0 0.0 -0.001234
          ...
        </DataArray>
      </PointData>
      <CellData Tensors="Stress">
        <DataArray type="Float64" Name="Stress" NumberOfComponents="6">
          1.234e6 2.345e6 3.456e6 1.234e5 2.345e5 3.456e5
          ...
        </DataArray>
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
```

### 7.3 输入解析器设计

```fortran
module input_parser_mod
    use precision_mod
    implicit none

    ! 解析器状态
    type :: parser_state
        character(len=512) :: filename = ''
        integer(IP) :: unit = -1
        integer(IP) :: line_number = 0
        character(len=256) :: current_line = ''
        logical :: eof = .false.
        integer(IP) :: error_count = 0
    end type

    ! 主解析器
    type :: input_parser
        type(parser_state) :: state
    contains
        procedure :: open => parser_open
        procedure :: close => parser_close
        procedure :: read_line => parser_read_line
        procedure :: parse_keyword => parser_parse_keyword
        procedure :: parse_value => parser_parse_value
        procedure :: skip_to => parser_skip_to
        procedure :: expect => parser_expect
        procedure :: error => parser_error
    end type

    ! 关键字处理接口
    abstract interface
        subroutine keyword_handler(parser, context, ierr)
            import :: input_parser, IP
            class(input_parser), intent(inout) :: parser
            class(*), intent(inout) :: context
            integer(IP), intent(out) :: ierr
        end subroutine
    end interface

    ! 关键字注册表
    type :: keyword_entry
        character(len=32) :: keyword
        procedure(keyword_handler), pointer, nopass :: handler
    end type

end module input_parser_mod
```

---

## 8. 错误处理策略

### 8.1 错误码定义

```fortran
module error_codes_mod
    use precision_mod
    implicit none

    ! 成功
    integer(IP), parameter :: ERR_OK = 0

    ! 通用错误 (1-99)
    integer(IP), parameter :: ERR_UNKNOWN = 1
    integer(IP), parameter :: ERR_NOT_IMPLEMENTED = 2
    integer(IP), parameter :: ERR_INVALID_ARGUMENT = 3

    ! 内存错误 (100-199)
    integer(IP), parameter :: ERR_ALLOC_FAILED = 100
    integer(IP), parameter :: ERR_DEALLOC_FAILED = 101
    integer(IP), parameter :: ERR_OUT_OF_MEMORY = 102

    ! 文件错误 (200-299)
    integer(IP), parameter :: ERR_FILE_NOT_FOUND = 200
    integer(IP), parameter :: ERR_FILE_OPEN = 201
    integer(IP), parameter :: ERR_FILE_READ = 202
    integer(IP), parameter :: ERR_FILE_WRITE = 203
    integer(IP), parameter :: ERR_FILE_FORMAT = 204

    ! 输入错误 (300-399)
    integer(IP), parameter :: ERR_PARSE_ERROR = 300
    integer(IP), parameter :: ERR_MISSING_KEYWORD = 301
    integer(IP), parameter :: ERR_INVALID_VALUE = 302
    integer(IP), parameter :: ERR_INCONSISTENT_DATA = 303

    ! 网格错误 (400-499)
    integer(IP), parameter :: ERR_INVALID_MESH = 400
    integer(IP), parameter :: ERR_INVALID_NODE = 401
    integer(IP), parameter :: ERR_INVALID_ELEMENT = 402
    integer(IP), parameter :: ERR_NEGATIVE_JACOBIAN = 403

    ! 材料错误 (500-599)
    integer(IP), parameter :: ERR_INVALID_MATERIAL = 500
    integer(IP), parameter :: ERR_MATERIAL_NOT_FOUND = 501
    integer(IP), parameter :: ERR_YIELD_SURFACE = 502

    ! 求解器错误 (600-699)
    integer(IP), parameter :: ERR_SOLVER_FAILED = 600
    integer(IP), parameter :: ERR_SINGULAR_MATRIX = 601
    integer(IP), parameter :: ERR_NO_CONVERGENCE = 602
    integer(IP), parameter :: ERR_DIVERGENCE = 603

    ! 数值错误 (700-799)
    integer(IP), parameter :: ERR_NAN_DETECTED = 700
    integer(IP), parameter :: ERR_INF_DETECTED = 701
    integer(IP), parameter :: ERR_NEGATIVE_STIFFNESS = 702

end module error_codes_mod
```

### 8.2 错误处理类型

```fortran
module error_handling_mod
    use precision_mod
    use error_codes_mod
    implicit none

    ! 错误信息类型
    type :: error_info
        integer(IP) :: code = ERR_OK
        character(len=256) :: message = ''
        character(len=64) :: source_file = ''
        character(len=64) :: source_proc = ''
        integer(IP) :: source_line = 0
        type(error_info), pointer :: cause => null()  ! 链式错误
    contains
        procedure :: is_error => error_is_error
        procedure :: set => error_set
        procedure :: chain => error_chain
        procedure :: print => error_print
        procedure :: clear => error_clear
    end type

    ! 全局错误处理配置
    type :: error_config
        logical :: stop_on_error = .true.
        logical :: print_stack_trace = .true.
        integer(IP) :: log_unit = 6
    end type

    type(error_config), save :: global_error_config

contains

    function is_error(self) result(flag)
        class(error_info), intent(in) :: self
        logical :: flag
        flag = (self%code /= ERR_OK)
    end function

    subroutine error_set(self, code, message, file, proc, line)
        class(error_info), intent(inout) :: self
        integer(IP), intent(in) :: code
        character(len=*), intent(in) :: message
        character(len=*), intent(in), optional :: file, proc
        integer(IP), intent(in), optional :: line

        self%code = code
        self%message = message
        if (present(file)) self%source_file = file
        if (present(proc)) self%source_proc = proc
        if (present(line)) self%source_line = line

        if (global_error_config%stop_on_error .and. code /= ERR_OK) then
            call self%print()
            error stop code
        end if
    end subroutine

    subroutine error_print(self)
        class(error_info), intent(in) :: self
        type(error_info), pointer :: current

        write(global_error_config%log_unit, '(A)') repeat('=', 70)
        write(global_error_config%log_unit, '(A)') 'ERROR REPORT'
        write(global_error_config%log_unit, '(A)') repeat('=', 70)

        write(global_error_config%log_unit, '(A,I0)') 'Error Code: ', self%code
        write(global_error_config%log_unit, '(A,A)') 'Message: ', trim(self%message)

        if (len_trim(self%source_file) > 0) then
            write(global_error_config%log_unit, '(A,A,A,I0)') &
                'Location: ', trim(self%source_file), ':', self%source_line
        end if

        if (len_trim(self%source_proc) > 0) then
            write(global_error_config%log_unit, '(A,A)') &
                'Procedure: ', trim(self%source_proc)
        end if

        ! 打印错误链
        current => self%cause
        do while (associated(current))
            write(global_error_config%log_unit, '(A)') 'Caused by:'
            write(global_error_config%log_unit, '(A,I0)') '  Code: ', current%code
            write(global_error_config%log_unit, '(A,A)') '  Message: ', trim(current%message)
            current => current%cause
        end do

        write(global_error_config%log_unit, '(A)') repeat('=', 70)
    end subroutine

end module error_handling_mod
```

### 8.3 使用模式

```fortran
! 宏定义 (通过预处理器)
#define CHECK_ERROR(ierr) if (ierr /= ERR_OK) return
#define SET_ERROR(err, code, msg) call err%set(code, msg, __FILE__, __PROCEDURE__, __LINE__)

! 使用示例
subroutine calculate_stiffness(elem, coords, D, Ke, err)
    class(element_base), intent(in) :: elem
    real(WP), intent(in) :: coords(:,:), D(:,:)
    real(WP), intent(out) :: Ke(:,:)
    type(error_info), intent(out) :: err

    real(WP) :: detJ
    integer(IP) :: igauss, ierr

    err%code = ERR_OK

    ! 检查输入
    if (size(coords, 2) /= elem%n_nodes) then
        call err%set(ERR_INVALID_ARGUMENT, &
            'Coordinate array size mismatch', &
            __FILE__, 'calculate_stiffness', __LINE__)
        return
    end if

    Ke = ZERO

    do igauss = 1, elem%n_gauss
        call elem%b_matrix(coords, igauss, B, detJ, ierr)

        ! 检查雅可比
        if (detJ <= ZERO) then
            call err%set(ERR_NEGATIVE_JACOBIAN, &
                'Negative Jacobian determinant detected', &
                __FILE__, 'calculate_stiffness', __LINE__)
            return
        end if

        ! ... 继续计算
    end do

end subroutine
```

---

## 9. 扩展机制

### 9.1 单元扩展

添加新单元类型的步骤：

```fortran
! 1. 创建新单元模块 (例如 wedge6_element_mod.f90)
module wedge6_element_mod
    use precision_mod
    use element_interface_mod
    implicit none

    type, extends(element_base) :: wedge6_element
    contains
        procedure :: shape_functions => wedge6_shape
        procedure :: b_matrix => wedge6_bmatrix
        procedure :: stiffness => wedge6_stiffness
        procedure :: internal_force => wedge6_internal_force
    end type

contains
    ! 实现所有必需的方法
    ! ...
end module

! 2. 在工厂中注册
module element_factory_mod
    use wedge6_element_mod  ! 添加新模块
    ! ...

    function create_element(type_id) result(elem)
        integer(IP), intent(in) :: type_id
        class(element_base), allocatable :: elem

        select case (type_id)
        ! ... 现有单元
        case (ELEM_WEDGE6)
            allocate(wedge6_element :: elem)
        end select
    end function
end module

! 3. 在配置中添加单元类型
! element_types_mod 中添加 ELEM_WEDGE6 = 13
```

### 9.2 材料扩展

```fortran
! 添加新材料类型的步骤

! 1. 创建新材料模块
module hyperbolic_material_mod
    use precision_mod
    use material_interface_mod
    implicit none

    type, extends(material_base) :: hyperbolic_material
        ! 材料参数
        real(WP) :: K, n, Rf     ! 双曲线参数
        real(WP) :: c, phi       ! 强度参数
        real(WP) :: Pa           ! 大气压力
    contains
        procedure :: tangent_modulus => hyper_tangent
        procedure :: stress_update => hyper_update
        procedure :: init_from_params => hyper_init
    end type

contains
    ! 实现方法
end module

! 2. 在工厂中注册
! 3. 在输入解析器中添加解析逻辑
```

### 9.3 求解器扩展

```fortran
! 添加新求解器

! 1. 实现求解器接口
module gmres_solver_mod
    use solver_interface_mod
    implicit none

    type, extends(solver_base) :: gmres_solver
        ! GMRES特定参数
        integer(IP) :: restart = 30
        real(WP), allocatable :: basis(:,:)
    contains
        procedure :: setup => gmres_setup
        procedure :: factorize => gmres_factorize  ! 可能为空操作
        procedure :: solve => gmres_solve
        procedure :: cleanup => gmres_cleanup
    end type

end module

! 2. 在求解器管理器中注册
```

### 9.4 输出格式扩展

```fortran
! 添加新输出格式

! 1. 实现输出接口
module exodus_writer_mod
    use output_interface_mod
    implicit none

    type, extends(output_writer_base) :: exodus_writer
        integer(IP) :: exoid = -1  ! Exodus文件ID
    contains
        procedure :: init => exodus_init
        procedure :: write_mesh => exodus_write_mesh
        ! ... 其他方法
    end type

end module

! 2. 在输出管理器中注册
```

---

## 10. 文件组织

### 10.1 目录结构

```
HSTAR-Next/
├── CMakeLists.txt                 # 主CMake文件
├── README.md
├── LICENSE
├── VERSION                        # 版本号文件
│
├── src/                           # 源代码
│   ├── CMakeLists.txt
│   │
│   ├── core/                      # 核心层
│   │   ├── CMakeLists.txt
│   │   ├── precision_mod.f90
│   │   ├── tensor_mod.f90
│   │   ├── sparse_matrix_mod.f90
│   │   ├── gauss_quadrature_mod.f90
│   │   ├── error_handling_mod.f90
│   │   └── error_codes_mod.f90
│   │
│   ├── data/                      # 数据层
│   │   ├── CMakeLists.txt
│   │   ├── mesh_mod.f90
│   │   ├── node_mod.f90
│   │   ├── state_manager_mod.f90
│   │   └── result_store_mod.f90
│   │
│   ├── domain/                    # 领域层
│   │   ├── CMakeLists.txt
│   │   │
│   │   ├── elements/              # 单元库
│   │   │   ├── element_interface_mod.f90
│   │   │   ├── element_types_mod.f90
│   │   │   ├── element_factory_mod.f90
│   │   │   ├── line2_element_mod.f90
│   │   │   ├── tri3_element_mod.f90
│   │   │   ├── tri6_element_mod.f90
│   │   │   ├── quad4_element_mod.f90
│   │   │   ├── quad8_element_mod.f90
│   │   │   ├── tet4_element_mod.f90
│   │   │   ├── tet10_element_mod.f90
│   │   │   ├── hex8_element_mod.f90
│   │   │   └── hex20_element_mod.f90
│   │   │
│   │   ├── materials/             # 材料库
│   │   │   ├── material_interface_mod.f90
│   │   │   ├── material_types_mod.f90
│   │   │   ├── material_factory_mod.f90
│   │   │   ├── elastic_mod.f90
│   │   │   ├── mohr_coulomb_mod.f90
│   │   │   ├── drucker_prager_mod.f90
│   │   │   ├── von_mises_mod.f90
│   │   │   ├── duncan_chang_mod.f90
│   │   │   ├── cam_clay_mod.f90
│   │   │   ├── goodman_joint_mod.f90
│   │   │   ├── creep_mod.f90
│   │   │   └── contact_mod.f90
│   │   │
│   │   ├── boundary/              # 边界条件
│   │   │   ├── bc_interface_mod.f90
│   │   │   ├── dirichlet_bc_mod.f90
│   │   │   ├── neumann_bc_mod.f90
│   │   │   └── contact_bc_mod.f90
│   │   │
│   │   └── loads/                 # 荷载
│   │       ├── load_interface_mod.f90
│   │       ├── body_force_mod.f90
│   │       ├── surface_load_mod.f90
│   │       └── concentrated_load_mod.f90
│   │
│   ├── service/                   # 服务层
│   │   ├── CMakeLists.txt
│   │   ├── assembler_mod.f90
│   │   ├── boundary_manager_mod.f90
│   │   ├── load_manager_mod.f90
│   │   │
│   │   ├── solver/                # 求解器
│   │   │   ├── solver_interface_mod.f90
│   │   │   ├── solver_factory_mod.f90
│   │   │   ├── pardiso_solver_mod.f90
│   │   │   ├── cg_solver_mod.f90
│   │   │   └── gmres_solver_mod.f90
│   │   │
│   │   └── output/                # 输出
│   │       ├── output_interface_mod.f90
│   │       ├── output_manager_mod.f90
│   │       ├── gid_writer_mod.f90
│   │       └── vtk_writer_mod.f90
│   │
│   ├── analysis/                  # 分析层
│   │   ├── CMakeLists.txt
│   │   ├── fem_engine_mod.f90
│   │   ├── newton_raphson_mod.f90
│   │   ├── static_analysis_mod.f90
│   │   ├── dynamic_analysis_mod.f90
│   │   ├── thermal_analysis_mod.f90
│   │   └── coupled_analysis_mod.f90
│   │
│   ├── io/                        # 输入/输出
│   │   ├── CMakeLists.txt
│   │   ├── input_parser_mod.f90
│   │   ├── mesh_reader_mod.f90
│   │   ├── material_reader_mod.f90
│   │   └── load_reader_mod.f90
│   │
│   └── app/                       # 应用层
│       ├── CMakeLists.txt
│       ├── main.f90
│       ├── fem_app_mod.f90
│       └── cli_mod.f90
│
├── tests/                         # 测试
│   ├── CMakeLists.txt
│   │
│   ├── unit/                      # 单元测试
│   │   ├── test_tensor.f90
│   │   ├── test_sparse_matrix.f90
│   │   ├── test_elements.f90
│   │   ├── test_materials.f90
│   │   └── test_solver.f90
│   │
│   ├── integration/               # 集成测试
│   │   ├── test_assembly.f90
│   │   └── test_analysis.f90
│   │
│   └── benchmarks/                # 基准测试
│       ├── patch_test.f90
│       ├── cantilever_beam.f90
│       └── cook_membrane.f90
│
├── examples/                      # 示例
│   ├── elastic_block/
│   ├── dam_analysis/
│   └── slope_stability/
│
├── docs/                          # 文档
│   ├── design/                    # 设计文档
│   │   └── architecture-design.md
│   ├── user/                      # 用户手册
│   └── api/                       # API文档
│
├── tools/                         # 工具脚本
│   ├── mesh_converter.py
│   └── result_extractor.py
│
└── third_party/                   # 第三方库
    └── CMakeLists.txt
```

### 10.2 模块依赖规则

```
┌─────────────────────────────────────────────────────────────────┐
│                     DEPENDENCY RULES                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│ 1. 层间只允许向下依赖                                            │
│    app -> analysis -> service -> domain -> data -> core        │
│                                                                 │
│ 2. 同层模块通过接口通信                                          │
│    materials <--> elements (通过 core 类型)                     │
│                                                                 │
│ 3. 接口模块可被任何层使用                                        │
│    *_interface_mod 放在对应层的顶部                              │
│                                                                 │
│ 4. 工厂模块知道所有具体实现                                       │
│    element_factory 依赖所有 *_element_mod                       │
│                                                                 │
│ 5. 禁止循环依赖                                                 │
│    使用接口打破循环                                              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 10.3 命名规范

| 类型 | 规范 | 示例 |
|------|------|------|
| 模块 | 小写_下划线_mod | `elastic_material_mod` |
| 类型 | 小写_下划线 | `elastic_material`, `mesh_type` |
| 子程序 | 小写_下划线 | `calculate_stiffness` |
| 函数 | 小写_下划线 | `get_jacobian` |
| 常量 | 大写_下划线 | `ERR_OK`, `ELEM_HEX8` |
| 变量 | 小写_下划线 | `n_nodes`, `element_id` |
| 类型绑定过程 | 小写_下划线 | `init`, `cleanup`, `calc_stress` |
| 文件 | 小写_下划线.f90 | `elastic_material_mod.f90` |

---

## 附录 A: 快速参考

### A.1 核心类型表

| 类型 | 模块 | 用途 |
|------|------|------|
| `WP`, `IP`, `LP` | precision_mod | 精度定义 |
| `vec3` | tensor_mod | 3D向量 |
| `stress_tensor` | tensor_mod | 应力张量 |
| `strain_tensor` | tensor_mod | 应变张量 |
| `csr_matrix` | sparse_matrix_mod | 稀疏矩阵 |
| `mesh_type` | mesh_mod | 网格容器 |
| `element_base` | element_interface_mod | 单元基类 |
| `material_base` | material_interface_mod | 材料基类 |
| `material_state` | material_interface_mod | 材料状态 |
| `solver_base` | solver_interface_mod | 求解器基类 |
| `error_info` | error_handling_mod | 错误信息 |

### A.2 常用操作

```fortran
! 创建网格
type(mesh_type) :: mesh
call mesh%init(dimension=3)
call mesh%add_node(1, [0.0_WP, 0.0_WP, 0.0_WP])
call mesh%build_dof_numbering()

! 创建单元
class(element_base), allocatable :: elem
elem = create_element(ELEM_HEX8)

! 创建材料
class(material_base), allocatable :: mat
mat = create_material(MAT_MOHR_COULOMB, [E, nu, c, phi, psi])

! 创建求解器
class(solver_base), allocatable :: solver
solver = create_solver(SOLVER_PARDISO)
call solver%setup(n_equations, ierr)
call solver%solve(K, b, x, ierr)

! 错误处理
type(error_info) :: err
call some_operation(args, err)
if (err%is_error()) then
    call err%print()
    return
end if
```

---

## 11. 迁移实施路线

### 11.1 迁移原则

- 先建护栏、后改内核。
- 先迁横切能力、后迁高风险物理能力。
- 全程可灰度、可回退，确保随时恢复到生产状态，防止一次性整体替换。

### 11.2 推荐周期表 (12周)

- M1（第1-2周）：基线与质量护栏——算例集与基准建立
- M2（第3-4周）：`legacy_facade` + `adapter` + `comparator` 影子运行，将影响降到最小
- M3（第5-6周）：迁移横切能力（输入解析、错误处理、输出管理）
- M4（第7-8周）：统一非线性迭代驱动 `newton_driver`，覆盖静力+线弹性路径
- M5（第9-10周）：分批迁移元素/材料：HEX8/TET4, Elastic/Mohr-Coulomb
- M6（第11-12周）：灰度切换、双线治理固化、迁移看板

### 11.3 双线维护规则

1. 新功能默认落在新接口层；紧急修老代码需同步补适配器契约测试。
2. 每个PR必须附"老跑 + 影子跑"差异报告。
3. 禁止在 `Fem.f90` 继续扩散重复迭代分支，统一收敛到 `newton_driver`。
4. 每周固定审查：数值偏差、性能回退阈值、高风险未迁移点。

### 11.4 详细执行方案

- 细分阶段、交付物、验收标准及具体执行清单，见 `docs/design/migration-roadmap.md`

---
**文档版本历史:**

| 版本 | 日期 | 修改内容 |
|------|------|---------|
| 1.0 | 2024-01 | 初始版本 |

---

*本文档为HSTAR-Next有限元程序的架构设计规范，所有新代码应遵循此规范。*
