# HSTAR-Next 架构设计文档
# Architecture Design Document

**Version:** 1.0
**Date:** 2024-01
**Status:** Draft

---

## 目录

1. [设计目标](#1-设计目标)
2. [系统架构](#2-系统架构)
3. [核心数据类型](#3-核心数据类型)
4. [模块规范](#4-模块规范)
5. [接口定义](#5-接口定义)
6. [数据流](#6-数据流)
7. [输入输出规范](#7-输入输出规范)
8. [错误处理策略](#8-错误处理策略)
9. [扩展机制](#9-扩展机制)
10. [文件组织](#10-文件组织)
11. [Migration Guide](#11-migration-guide)

---

## 1. 设计目标

### 1.1 核心目标

| 目标 | 描述 | 衡量标准 |
|------|------|---------|
| **可测试** | 每个模块可独立进行单元测试 | 100%模块有对应测试 |
| **可扩展** | 添加新单元、材料无需修改现有代码 | 符合开放封闭原则 |
| **可维护** | 代码清晰，易于理解和修改 | 遵循命名规范，有文档 |
| **高性能** | 保持与原程序相当或更优的性能 | 基准测试对比 |
| **兼容性** | 支持原有输入格式 | 读取旧格式文件 |

### 1.2 设计原则

1. **单一职责原则 (SRP)**: 每个模块只负责一件事
2. **依赖倒置原则 (DIP)**: 高层模块依赖抽象，不依赖具体实现
3. **接口隔离原则 (ISP)**: 接口小而专一
4. **显式依赖**: 所有依赖通过参数传递，禁止全局变量
5. **不可变优先**: 优先使用 `intent(in)`，状态变更明确

### 1.3 约束条件

- 语言: Fortran 2008 (部分2018特性)
- 编译器: Intel Fortran 2019+, GFortran 9+
- 外部依赖: Intel MKL (PARDISO, BLAS, LAPACK)
- 平台: Windows, Linux

---

## 2. 系统架构

### 2.1 分层架构

```
┌────────────────────────────────────────────────────────────────────┐
│                       APPLICATION LAYER                            │
│ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐               │
│ │  CLI App │ │Batch Mode│ │ Library  │ │ Python   │               │
│ │          │ │          │ │  API     │ │ Binding  │               │
│ └──────────┘ └──────────┘ └──────────┘ └──────────┘               │
└────────────────────────────────────────────────────────────────────┘
          │              │              │              │
          └──────────────┴──────────────┴──────────────┘
                                  │
┌─────────────────────────────────▼──────────────────────────────────┐
│                       ANALYSIS LAYER                               │
│ ┌──────────────┐ ┌──────────────┐ ┌──────────────┐                │
│ │StaticAnalysis│ │DynamicAnalysis│ │ThermalAnalysis│               │
│ └──────────────┘ └──────────────┘ └──────────────┘                │
│          └───────────────┴───────────────┘                         │
│                              │                                     │
│                   ┌──────────▼──────────┐                          │
│                   │   FEM_Engine        │                          │
│                   │ (核心计算引擎)       │                          │
│                   └─────────────────────┘                          │
└────────────────────────────────│───────────────────────────────────┘
                                 │
┌────────────────────────────────▼───────────────────────────────────┐
│                       SERVICE LAYER                                │
│ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐               │
│ │ Assembler│ │  Solver  │ │ Boundary │ │  Output  │               │
│ │          │ │  Manager │ │  Manager │ │  Manager │               │
│ └──────────┘ └──────────┘ └──────────┘ └──────────┘               │
└────────────────────────────────────────────────────────────────────┘
          │              │              │              │
┌─────────▼──────────────▼──────────────▼──────────────▼─────────────┐
│                       DOMAIN LAYER                                 │
│ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐               │
│ │ Elements │ │ Materials│ │ Boundary │ │  Loads   │               │
│ │ (Factory)│ │ (Factory)│ │Conditions│ │          │               │
│ └──────────┘ └──────────┘ └──────────┘ └──────────┘               │
└────────────────────────────────────────────────────────────────────┘
          │              │              │              │
          └──────────────┴──────────────┴──────────────┘
                                  │
┌─────────────────────────────────▼──────────────────────────────────┐
│                       DATA LAYER                                   │
│ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐          │
│ │     Mesh       │ │  StateManager  │ │  ResultStore   │          │
│ │  (节点/单元)    │ │ (高斯点状态)    │ │  (结果存储)    │          │
│ └────────────────┘ └────────────────┘ └────────────────┘          │
└────────────────────────────────────────────────────────────────────┘
            │                  │                  │
            └──────────────────┴──────────────────┘
                                          │
┌─────────────────────────────────────────▼──────────────────────────┐
│                       CORE LAYER                                   │
│ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐            │
│ │Precision│ │Tensor  │ │Sparse  │ │ Gauss  │ │ Error  │            │
│ │        │ │ Math   │ │Matrix  │ │ Quad   │ │Handle  │            │
│ └────────┘ └────────┘ └────────┘ └────────┘ └────────┘            │
└────────────────────────────────────────────────────────────────────┘
```

### 2.2 层职责

| 层 | 职责 | 依赖 |
|---|------|-----|
| **Application** | 用户交互、命令解析、批处理 | Analysis |
| **Analysis** | 分析流程控制、时间步进、收敛判断 | Service |
| **Service** | 组装、求解、边界处理、输出管理 | Domain, Data |
| **Domain** | 单元、材料、边界条件的业务逻辑 | Data, Core |
| **Data** | 网格存储、状态管理、结果存储 | Core |
| **Core** | 数学工具、精度定义、基础设施 | 无 |

### 2.3 模块依赖图

```
                                 ┌───────────┐
                                 │  fem_app  │
                                 └─────┬─────┘
                                       │
                    ┌──────────────────┼──────────────────┐
                    │                  │                  │
             ┌──────▼──────┐    ┌──────▼──────┐    ┌──────▼──────┐
             │  static     │    │  dynamic    │    │  thermal    │
             │ _analysis   │    │ _analysis   │    │ _analysis   │
             └──────┬──────┘    └──────┬──────┘    └──────┬──────┘
                    │                  │                  │
                    └──────────────────┼──────────────────┘
                                       │
                                 ┌─────▼─────┐
                                 │fem_engine │
                                 └─────┬─────┘
                                       │
          ┌────────────┬───────────────┼───────────────┬────────────┐
          │            │               │               │            │
   ┌──────▼──────┐┌────▼────┐┌─────────▼─────────┐┌────▼────┐┌──────▼──────┐
   │ assembler   ││ solver  ││  boundary_mgr     ││output_mgr││ input_mgr   │
   └──────┬──────┘└────┬────┘└─────────┬─────────┘└────┬────┘└──────┬──────┘
          │            │               │               │            │
          │     ┌──────┴──────┐        │               │            │
          │     │             │        │               │            │
   ┌──────▼─────▼────┐  ┌─────▼────┐   │        ┌──────▼──────┐     │
   │element_factory  │  │pardiso_  │   │        │gid_writer   │     │
   │                 │  │solver    │   │        │vtk_writer   │     │
   └────────┬────────┘  └──────────┘   │        └─────────────┘     │
            │                          │                            │
   ┌────────▼────────┐          ┌──────▼──────┐              ┌──────▼──────┐
   │material_        │          │   mesh      │              │  input_     │
   │factory          │          │             │              │  parser     │
   └────────┬────────┘          └──────┬──────┘              └─────────────┘
            │                          │
            │     ┌────────────────────┤
            │     │                    │
   ┌────────▼─────▼────┐         ┌─────▼─────┐
   │state_manager      │         │   node    │
   └────────┬──────────┘         └───────────┘
            │
   ┌────────▼────────┐
   │   core          │
   │(precision,      │
   │ tensor,         │
   │ sparse,         │
   │ gauss)          │
   └─────────────────┘
```

---

## 3. 核心数据类型

### 3.1 精度定义 (precision_mod)

```fortran
module precision_mod
    use, intrinsic :: iso_fortran_env
    implicit none

    ! 工作精度
    integer, parameter :: WP = real64      ! 双精度 (8字节)
    integer, parameter :: IP = int32       ! 整数精度 (4字节)
    integer, parameter :: LP = int64       ! 长整数 (索引用)

    ! 数学常量
    real(WP), parameter :: PI = 3.141592653589793238_WP
    real(WP), parameter :: ZERO = 0.0_WP
    real(WP), parameter :: ONE = 1.0_WP
    real(WP), parameter :: HALF = 0.5_WP

    ! 容差
    real(WP), parameter :: EPS_MACHINE = epsilon(1.0_WP)
    real(WP), parameter :: EPS_GEOMETRY = 1.0e-10_WP
    real(WP), parameter :: EPS_STRESS = 1.0e-6_WP
    real(WP), parameter :: EPS_CONVERGENCE = 1.0e-8_WP

end module precision_mod
```

### 3.2 向量与张量 (tensor_mod)

```fortran
module tensor_mod
    use precision_mod
    implicit none

    ! 向量类型 (用于小向量，避免allocatable)
    type :: vec3
        real(WP) :: v(3) = ZERO
    contains
        procedure :: norm => vec3_norm
        procedure :: dot => vec3_dot
        procedure :: cross => vec3_cross
    end type

    ! 应力/应变向量 (Voigt表示)
    type :: stress_tensor
        real(WP) :: s(6) = ZERO  ! [s11, s22, s33, s12, s23, s13]
    contains
        procedure :: von_mises => stress_von_mises
        procedure :: mean => stress_mean
        procedure :: deviatoric => stress_deviatoric
        procedure :: invariants => stress_invariants
        procedure :: principal => stress_principal
    end type

    type :: strain_tensor
        real(WP) :: e(6) = ZERO  ! [e11, e22, e33, 2*e12, 2*e23, 2*e13]
    contains
        procedure :: volumetric => strain_volumetric
        procedure :: deviatoric => strain_deviatoric
    end type

    ! 本构矩阵 (6x6 for 3D)
    type :: constitutive_matrix
        real(WP) :: D(6,6) = ZERO
    contains
        procedure :: apply => const_apply  ! stress = D * strain
    end type

end module tensor_mod
```

### 3.3 稀疏矩阵 (sparse_matrix_mod)

```fortran
module sparse_matrix_mod
    use precision_mod
    implicit none

    ! CSR格式稀疏矩阵
    type :: csr_matrix
        integer(IP) :: n_rows = 0
        integer(IP) :: n_cols = 0
        integer(LP) :: nnz = 0          ! 非零元素数

        integer(LP), allocatable :: row_ptr(:)   ! 行指针 (n_rows + 1)
        integer(IP), allocatable :: col_idx(:)   ! 列索引 (nnz)
        real(WP), allocatable :: values(:)       ! 值 (nnz)

        logical :: is_symmetric = .false.
        logical :: is_assembled = .false.
    contains
        procedure :: init => csr_init
        procedure :: add_value => csr_add_value
        procedure :: add_matrix => csr_add_matrix
        procedure :: finalize => csr_finalize
        procedure :: clear_values => csr_clear_values
        procedure :: matvec => csr_matvec
        procedure :: cleanup => csr_cleanup
    end type

    ! 用于组装的COO格式
    type :: coo_entry
        integer(IP) :: row, col
        real(WP) :: val
    end type

    type :: coo_matrix
        type(coo_entry), allocatable :: entries(:)
        integer(LP) :: count = 0
        integer(LP) :: capacity = 0
    contains
        procedure :: add => coo_add
        procedure :: to_csr => coo_to_csr
    end type

end module sparse_matrix_mod
```

### 3.4 网格数据结构 (mesh_mod)

```fortran
module mesh_mod
    use precision_mod
    implicit none

    !=========================================
    ! 节点
    !=========================================
    type :: node_data
        integer(IP) :: id
        real(WP) :: coord(3)           ! 坐标 (始终3D)
        integer(IP) :: dof_start = 0   ! 第一个自由度的全局编号
        integer(IP) :: n_dof = 0       ! 该节点的自由度数
    end type

    !=========================================
    ! 单元连接
    !=========================================
    type :: element_connectivity
        integer(IP) :: id
        integer(IP) :: element_type    ! 单元类型枚举
        integer(IP) :: material_id
        integer(IP) :: group_id
        integer(IP), allocatable :: node_ids(:)
    end type

    !=========================================
    ! 节点集合 (用于边界条件)
    !=========================================
    type :: node_set
        character(len=32) :: name
        integer(IP), allocatable :: node_ids(:)
    end type

    !=========================================
    ! 单元集合 (用于材料分配)
    !=========================================
    type :: element_set
        character(len=32) :: name
        integer(IP), allocatable :: element_ids(:)
    end type

    !=========================================
    ! 网格容器
    !=========================================
    type :: mesh_type
        ! 基本信息
        integer(IP) :: dimension = 3
        integer(IP) :: n_nodes = 0
        integer(IP) :: n_elements = 0
        integer(IP) :: n_total_dof = 0

        ! 数据存储
        type(node_data), allocatable :: nodes(:)
        type(element_connectivity), allocatable :: elements(:)

        ! 集合
        type(node_set), allocatable :: node_sets(:)
        type(element_set), allocatable :: element_sets(:)

        ! 自由度映射
        integer(IP), allocatable :: dof_map(:,:)  ! (n_dof_per_node, n_nodes)

    contains
        procedure :: init => mesh_init
        procedure :: add_node => mesh_add_node
        procedure :: add_element => mesh_add_element
        procedure :: add_node_set => mesh_add_node_set
        procedure :: build_dof_numbering => mesh_build_dof
        procedure :: get_element_coords => mesh_get_elem_coords
        procedure :: get_element_dofs => mesh_get_elem_dofs
        procedure :: cleanup => mesh_cleanup
    end type

end module mesh_mod
```

---

## 4. 模块规范

### 4.1 单元类型枚举

```fortran
module element_types_mod
    use precision_mod
    implicit none

    ! 单元类型枚举
    integer(IP), parameter :: ELEM_LINE2 = 1      ! 2节点线单元
    integer(IP), parameter :: ELEM_LINE3 = 2      ! 3节点线单元
    integer(IP), parameter :: ELEM_TRI3 = 3       ! 3节点三角形
    integer(IP), parameter :: ELEM_TRI6 = 4       ! 6节点三角形
    integer(IP), parameter :: ELEM_QUAD4 = 5      ! 4节点四边形
    integer(IP), parameter :: ELEM_QUAD8 = 6      ! 8节点四边形
    integer(IP), parameter :: ELEM_QUAD9 = 7      ! 9节点四边形
    integer(IP), parameter :: ELEM_TET4 = 8       ! 4节点四面体
    integer(IP), parameter :: ELEM_TET10 = 9      ! 10节点四面体
    integer(IP), parameter :: ELEM_HEX8 = 10      ! 8节点六面体
    integer(IP), parameter :: ELEM_HEX20 = 11     ! 20节点六面体
    integer(IP), parameter :: ELEM_HEX27 = 12     ! 27节点六面体
    integer(IP), parameter :: ELEM_WEDGE6 = 13    ! 6节点楔形
    integer(IP), parameter :: ELEM_PYRAMID5 = 14  ! 5节点金字塔

    ! 特殊单元
    integer(IP), parameter :: ELEM_BEAM2 = 20     ! 2节点梁
    integer(IP), parameter :: ELEM_CONTACT_P2P = 30  ! 点对点接触
    integer(IP), parameter :: ELEM_INTERFACE = 40    ! 界面单元

    ! 单元信息查询
    type :: element_info
        integer(IP) :: type_id
        character(len=16) :: name
        integer(IP) :: n_nodes
        integer(IP) :: dimension        ! 1D, 2D, 3D
        integer(IP) :: n_gauss_default  ! 默认高斯点数
        integer(IP) :: n_edges
        integer(IP) :: n_faces
    end type

    type(element_info), parameter :: ELEM_INFO(14) = [ &
        element_info(ELEM_LINE2, 'LINE2', 2, 1, 2, 0, 0), &
        element_info(ELEM_TRI3, 'TRI3', 3, 2, 1, 3, 0), &
        element_info(ELEM_TRI6, 'TRI6', 6, 2, 3, 3, 0), &
        element_info(ELEM_QUAD4, 'QUAD4', 4, 2, 4, 4, 0), &
        element_info(ELEM_QUAD8, 'QUAD8', 8, 2, 9, 4, 0), &
        element_info(ELEM_TET4, 'TET4', 4, 3, 1, 6, 4), &
        element_info(ELEM_TET10, 'TET10', 10, 3, 4, 6, 4), &
        element_info(ELEM_HEX8, 'HEX8', 8, 3, 8, 12, 6), &
        element_info(ELEM_HEX20, 'HEX20', 20, 3, 27, 12, 6), &
        ! ... 其他单元
    ]

end module element_types_mod
```

### 4.2 材料类型枚举

```fortran
module material_types_mod
    use precision_mod
    implicit none

    ! 材料类型枚举
    integer(IP), parameter :: MAT_ELASTIC = 1
    integer(IP), parameter :: MAT_ELASTIC_ORTHOTROPIC = 2
    integer(IP), parameter :: MAT_MOHR_COULOMB = 10
    integer(IP), parameter :: MAT_DRUCKER_PRAGER = 11
    integer(IP), parameter :: MAT_VON_MISES = 12
    integer(IP), parameter :: MAT_TRESCA = 13
    integer(IP), parameter :: MAT_CAM_CLAY = 20
    integer(IP), parameter :: MAT_DUNCAN_CHANG = 30
    integer(IP), parameter :: MAT_HYPERBOLIC = 31
    integer(IP), parameter :: MAT_GOODMAN = 40      ! 节理
    integer(IP), parameter :: MAT_CONTACT = 50      ! 接触
    integer(IP), parameter :: MAT_CREEP_POWER = 60  ! 幂律蠕变
    integer(IP), parameter :: MAT_CREEP_BURGERS = 61 ! Burgers蠕变
    integer(IP), parameter :: MAT_CONCRETE_DAMAGE = 70 ! 混凝土损伤

    ! 材料行为标志
    integer(IP), parameter :: BEHAVIOR_LINEAR = 0
    integer(IP), parameter :: BEHAVIOR_NONLINEAR = 1
    integer(IP), parameter :: BEHAVIOR_RATE_DEPENDENT = 2

end module material_types_mod
```

### 4.3 分析类型枚举

```fortran
module analysis_types_mod
    use precision_mod
    implicit none

    ! 分析类型
    integer(IP), parameter :: ANALYSIS_STATIC = 1
    integer(IP), parameter :: ANALYSIS_DYNAMIC = 2
    integer(IP), parameter :: ANALYSIS_MODAL = 3
    integer(IP), parameter :: ANALYSIS_THERMAL = 10
    integer(IP), parameter :: ANALYSIS_SEEPAGE = 11
    integer(IP), parameter :: ANALYSIS_COUPLED_TM = 20   ! 热力耦合
    integer(IP), parameter :: ANALYSIS_COUPLED_HM = 21   ! 水力耦合
    integer(IP), parameter :: ANALYSIS_COUPLED_THM = 22  ! 热水力耦合

    ! 时间积分方法
    integer(IP), parameter :: TIME_STATIC = 0
    integer(IP), parameter :: TIME_NEWMARK = 1
    integer(IP), parameter :: TIME_HHT = 2
    integer(IP), parameter :: TIME_GENERALIZED_ALPHA = 3
    integer(IP), parameter :: TIME_EXPLICIT_CD = 10  ! 中心差分

    ! 非线性求解方法
    integer(IP), parameter :: NONLIN_NEWTON_RAPHSON = 1
    integer(IP), parameter :: NONLIN_MODIFIED_NR = 2
    integer(IP), parameter :: NONLIN_BFGS = 3
    integer(IP), parameter :: NONLIN_ARC_LENGTH = 4

end module analysis_types_mod
```

---

## 5. 接口定义

### 5.1 单元抽象接口

```fortran
module element_interface_mod
    use precision_mod
    use tensor_mod
    implicit none

    !=========================================
    ! 单元抽象基类
    !=========================================
    type, abstract :: element_base
        integer(IP) :: id = 0
        integer(IP) :: type_id = 0
        integer(IP) :: n_nodes = 0
        integer(IP) :: n_dof = 0
        integer(IP) :: n_gauss = 0
        integer(IP) :: n_stress = 0        ! 应力分量数
    contains
        ! 必须实现的方法
        procedure(shape_if), deferred :: shape_functions
        procedure(bmatrix_if), deferred :: b_matrix
        procedure(stiffness_if), deferred :: stiffness
        procedure(internal_force_if), deferred :: internal_force

        ! 可覆盖的方法
        procedure :: mass_matrix => element_mass_default
        procedure :: get_gauss_points => element_gauss_default
    end type

    !=========================================
    ! 抽象接口定义
    !=========================================
    abstract interface

        ! 形函数计算
        pure subroutine shape_if(self, xi, N, dN_dxi)
            import :: element_base, WP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: xi(:)           ! 自然坐标
            real(WP), intent(out) :: N(:)           ! 形函数值
            real(WP), intent(out) :: dN_dxi(:,:)    ! 形函数导数
        end subroutine

        ! B矩阵计算
        pure subroutine bmatrix_if(self, coords, igauss, B, detJ, ierr)
            import :: element_base, WP, IP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: coords(:,:)     ! 节点坐标
            integer(IP), intent(in) :: igauss      ! 高斯点编号
            real(WP), intent(out) :: B(:,:)         ! B矩阵
            real(WP), intent(out) :: detJ           ! 雅可比行列式
            integer(IP), intent(out) :: ierr       ! 错误码
        end subroutine

        ! 刚度矩阵计算
        pure subroutine stiffness_if(self, coords, D, Ke, ierr)
            import :: element_base, WP, IP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: coords(:,:)     ! 节点坐标
            real(WP), intent(in) :: D(:,:)          ! 本构矩阵
            real(WP), intent(out) :: Ke(:,:)        ! 单元刚度
            integer(IP), intent(out) :: ierr
        end subroutine

        ! 内力向量计算
        pure subroutine internal_force_if(self, coords, stress, Fe, ierr)
            import :: element_base, WP, IP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: coords(:,:)     ! 节点坐标
            real(WP), intent(in) :: stress(:,:)     ! 高斯点应力
            real(WP), intent(out) :: Fe(:)          ! 单元内力
            integer(IP), intent(out) :: ierr
        end subroutine

    end interface

end module element_interface_mod
```

### 5.2 材料抽象接口

```fortran
module material_interface_mod
    use precision_mod
    use tensor_mod
    implicit none

    !=========================================
    ! 材料状态 (高斯点历史变量)
    !=========================================
    type :: material_state
        real(WP), allocatable :: stress(:)           ! 应力
        real(WP), allocatable :: strain(:)           ! 总应变
        real(WP), allocatable :: plastic_strain(:)   ! 塑性应变
        real(WP), allocatable :: internal_vars(:)    ! 内变量
        integer(IP) :: status = 0                    ! 0=弹性, 1=塑性
    contains
        procedure :: init => state_init
        procedure :: copy_from => state_copy
        procedure :: reset => state_reset
    end type

    !=========================================
    ! 材料抽象基类
    !=========================================
    type, abstract :: material_base
        integer(IP) :: id = 0
        integer(IP) :: type_id = 0
        character(len=32) :: name = ''
        integer(IP) :: n_stress = 6        ! 应力分量数
        integer(IP) :: n_internal_vars = 0 ! 内变量数
        real(WP) :: density = 0.0_WP
    contains
        ! 必须实现
        procedure(tangent_if), deferred :: tangent_modulus
        procedure(stress_update_if), deferred :: stress_update

        ! 可覆盖实现
        procedure :: is_linear => material_is_linear_default
        procedure :: get_elastic_props => material_elastic_default
    end type

    !=========================================
    ! 抽象接口
    !=========================================
    abstract interface

        ! 切线刚度矩阵
        pure subroutine tangent_if(self, state, D)
            import :: material_base, material_state, WP
            class(material_base), intent(in) :: self
            type(material_state), intent(in) :: state
            real(WP), intent(out) :: D(:,:)
        end subroutine

        ! 应力更新 (本构积分)
        subroutine stress_update_if(self, d_strain, state, stress, D, ierr)
            import :: material_base, material_state, WP, IP
            class(material_base), intent(in) :: self
            real(WP), intent(in) :: d_strain(:)
            type(material_state), intent(inout) :: state
            real(WP), intent(out) :: stress(:)
            real(WP), intent(out) :: D(:,:)     ! 一致切线刚度
            integer(IP), intent(out) :: ierr
        end subroutine

    end interface

end module material_interface_mod
```

### 5.3 求解器抽象接口

```fortran
module solver_interface_mod
    use precision_mod
    use sparse_matrix_mod
    implicit none

    !=========================================
    ! 求解器抽象基类
    !=========================================
    type, abstract :: solver_base
        integer(IP) :: n_equations = 0
        integer(IP) :: max_iterations = 1000
        real(WP) :: tolerance = 1.0e-10_WP
        logical :: is_initialized = .false.
        logical :: is_factorized = .false.
    contains
        procedure(setup_if), deferred :: setup
        procedure(factorize_if), deferred :: factorize
        procedure(solve_if), deferred :: solve
        procedure(cleanup_if), deferred :: cleanup
    end type

    abstract interface

        subroutine setup_if(self, n_eq, ierr)
            import :: solver_base, IP
            class(solver_base), intent(inout) :: self
            integer(IP), intent(in) :: n_eq
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine factorize_if(self, A, ierr)
            import :: solver_base, csr_matrix, IP
            class(solver_base), intent(inout) :: self
            type(csr_matrix), intent(in) :: A
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine solve_if(self, A, b, x, ierr)
            import :: solver_base, csr_matrix, WP, IP
            class(solver_base), intent(inout) :: self
            type(csr_matrix), intent(in) :: A
            real(WP), intent(in) :: b(:)
            real(WP), intent(inout) :: x(:)
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine cleanup_if(self)
            import :: solver_base
            class(solver_base), intent(inout) :: self
        end subroutine

    end interface

end module solver_interface_mod
```

### 5.4 输出器抽象接口

```fortran
module output_interface_mod
    use precision_mod
    use mesh_mod
    implicit none

    !=========================================
    ! 输出器抽象基类
    !=========================================
    type, abstract :: output_writer_base
        character(len=256) :: base_path = ''
        integer(IP) :: current_step = 0
        real(WP) :: current_time = 0.0_WP
    contains
        procedure(init_output_if), deferred :: init
        procedure(write_mesh_if), deferred :: write_mesh
        procedure(begin_step_if), deferred :: begin_step
        procedure(write_nodal_if), deferred :: write_nodal_vector
        procedure(write_nodal_if), deferred :: write_nodal_scalar
        procedure(write_elemental_if), deferred :: write_elemental_tensor
        procedure(end_step_if), deferred :: end_step
        procedure(finalize_if), deferred :: finalize
    end type

    abstract interface

        subroutine init_output_if(self, path, mesh, ierr)
            import :: output_writer_base, mesh_type, IP
            class(output_writer_base), intent(inout) :: self
            character(len=*), intent(in) :: path
            type(mesh_type), intent(in) :: mesh
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine write_mesh_if(self, mesh, ierr)
            import :: output_writer_base, mesh_type, IP
            class(output_writer_base), intent(inout) :: self
            type(mesh_type), intent(in) :: mesh
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine begin_step_if(self, step, time)
            import :: output_writer_base, IP, WP
            class(output_writer_base), intent(inout) :: self
            integer(IP), intent(in) :: step
            real(WP), intent(in) :: time
        end subroutine

        subroutine write_nodal_if(self, name, data, ierr)
            import :: output_writer_base, WP, IP
            class(output_writer_base), intent(inout) :: self
            character(len=*), intent(in) :: name
            real(WP), intent(in) :: data(:,:)
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine write_elemental_if(self, name, data, ierr)
            import :: output_writer_base, WP, IP
            class(output_writer_base), intent(inout) :: self
            character(len=*), intent(in) :: name
            real(WP), intent(in) :: data(:,:,:)  ! (n_components, n_gauss, n_elem)
            integer(IP), intent(out) :: ierr
        end subroutine

        subroutine end_step_if(self)
            import :: output_writer_base
            class(output_writer_base), intent(inout) :: self
        end subroutine

        subroutine finalize_if(self)
            import :: output_writer_base
            class(output_writer_base), intent(inout) :: self
        end subroutine

    end interface

end module output_interface_mod
```

---

## 6. 数据流

### 6.1 分析流程数据流

```
┌────────────────────────────────────────────────────────────────────┐
│                        INPUT PHASE                                 │
│                                                                    │
│ ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐        │
│ │Control   │-> │ Mesh     │-> │Materials │-> │ BCs &    │        │
│ │ File     │   │ File     │   │ File     │   │ Loads    │        │
│ └──────────┘   └──────────┘   └──────────┘   └──────────┘        │
│      │             │             │             │                  │
│      │             │             │             │                  │
│ ┌────▼─────────────▼─────────────▼─────────────▼────────────────┐ │
│ │                   INPUT PARSER                                │ │
│ └───────────────────────────────────────────────────────────────┘ │
└────────────────────────────────────────│───────────────────────────┘
                              │
                              │
┌─────────────────────────────▼──────────────────────────────────────┐
│                     INITIALIZATION PHASE                           │
│                                                                    │
│ ┌──────────────┐ ┌──────────────┐ ┌──────────────┐                │
│ │   Mesh       │ │  Element     │ │  Material    │                │
│ │  Builder     │ │  Factory     │ │  Factory     │                │
│ └──────────────┘ └──────────────┘ └──────────────┘                │
│        │                │                │                        │
│        │                │                │                        │
│ ┌──────▼────────────────▼────────────────▼──────────────────────┐ │
│ │                    FEM ENGINE                                 │ │
│ │ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐                  │ │
│ │ │ Mesh   │ │Elements│ │Materials│ │State   │                  │ │
│ │ │        │ │Array   │ │ Array  │ │Manager │                  │ │
│ │ └────────┘ └────────┘ └────────┘ └────────┘                  │ │
│ └───────────────────────────────────────────────────────────────┘ │
└────────────────────────────────│───────────────────────────────────┘
                            │
                            │
┌───────────────────────────▼────────────────────────────────────────┐
│                      ANALYSIS PHASE                                │
│                                                                    │
│ FOR each load_block:                                              │
│   FOR each increment:                                             │
│     ┌─────────────────────────────────────────────────────────┐   │
│     │             NEWTON-RAPHSON LOOP                         │   │
│     │                                                         │   │
│     │ ┌──────────┐    ┌──────────┐    ┌──────────┐           │   │
│     │ │Assemble  │--> │ Apply    │--> │ Solve    │           │   │
│     │ │   K      │    │  BCs     │    │ K*du=R   │           │   │
│     │ └──────────┘    └──────────┘    └──────────┘           │   │
│     │      │                                 │                │   │
│     │      │          ┌──────────┐           │                │   │
│     │      │--------->│ Check    │-----------│                │   │
│     │                 │Convergence│                           │   │
│     │                 └──────────┘                           │   │
│     │                       │                                │   │
│     │           ┌───────────┴───────────┐                    │   │
│     │           │                       │                    │   │
│     │    [Not Converged]           [Converged]               │   │
│     │     Update state              Commit state             │   │
│     │     Continue loop             Exit loop                │   │
│     └─────────────────────────────────────────────────────────┘   │
│                                                                    │
│     Write output if needed                                        │
└────────────────────────────────────────────────────────────────────┘
                            │
                            │
┌───────────────────────────▼────────────────────────────────────────┐
│                      OUTPUT PHASE                                  │
│                                                                    │
│ ┌────────────┐ ┌────────────┐ ┌────────────┐                      │
│ │ GiD Writer │ │ VTK Writer │ │ Text Writer│                      │
│ └────────────┘ └────────────┘ └────────────┘                      │
│                                                                    │
│ Output files: .res, .msh, .vtk, .dat                              │
└────────────────────────────────────────────────────────────────────┘
```

### 6.2 Newton-Raphson迭代数据流

```
┌─────────────────────────────────────────────────────────────────┐
│                   ITERATION START                               │
│                                                                 │
│ Input:                                                          │
│   u_n        : 上一收敛步位移                                   │
│   state_n    : 上一收敛步状态                                   │
│   F_ext      : 外力向量                                         │
│   load_factor: 荷载因子                                         │
└────────────────────────────────────│────────────────────────────┘
                             │
                             │
┌────────────────────────────▼────────────────────────────────────┐
│ STEP 1: Initialize iteration                                    │
│                                                                 │
│   u = u_n                                                       │
│   state = copy(state_n)                                         │
│   iter = 0                                                      │
└────────────────────────────────────│────────────────────────────┘
                             │
              ┌──────────────▼──────────────┐
              │    ITERATION LOOP           │
              │                             │
              │ ┌────────────────────────┐  │
              │ │STEP 2: Assemble K      │  │
              │ │                        │  │
              │ │For each element:       │  │
              │ │  D = mat.tangent()     │  │
              │ │  Ke = elem.stiff()     │  │
              │ │  K += assemble(Ke)     │  │
              │ └────────────────────────┘  │
              │            │                │
              │ ┌──────────▼─────────────┐  │
              │ │STEP 3: Internal F      │  │
              │ │                        │  │
              │ │For each element:       │  │
              │ │  σ = state.stress      │  │
              │ │  Fe = elem.Fint()      │  │
              │ │  F_int += Fe           │  │
              │ └────────────────────────┘  │
              │            │                │
              │ ┌──────────▼─────────────┐  │
              │ │STEP 4: Residual        │  │
              │ │                        │  │
              │ │R = λ*F_ext - F_int     │  │
              │ └────────────────────────┘  │
              │            │                │
              │ ┌──────────▼─────────────┐  │
              │ │STEP 5: Apply BC        │  │
              │ │                        │  │
              │ │K, R = apply_bc()       │  │
              │ └────────────────────────┘  │
              │            │                │
              │ ┌──────────▼─────────────┐  │
              │ │STEP 6: Solve           │  │
              │ │                        │  │
              │ │K * Δu = R              │  │
              │ └────────────────────────┘  │
              │            │                │
              │ ┌──────────▼─────────────┐  │
              │ │STEP 7: Update          │  │
              │ │                        │  │
              │ │u = u + Δu              │  │
              │ │ε = B * u               │  │
              │ │σ, D = mat.update()     │  │
              │ │state.update()          │  │
              │ └────────────────────────┘  │
              │            │                │
              │ ┌──────────▼─────────────┐  │
              │ │STEP 8: Check           │  │
              │ │                        │  │
              │ │||R|| < tol ?           │  │
              │ │iter < max_iter ?       │  │
              │ └────────────────────────┘  │
              │            │                │
              │     ┌──────┴──────┐         │
              │     │             │         │
              │ [Continue]    [Exit]        │
              │  iter += 1    Converged     │
              └─────────────────────────────┘
                     │            │
                     └────────────┘
                            │
                            │
┌───────────────────────────▼─────────────────────────────────────┐
│                   ITERATION END                                 │
│                                                                 │
│ If converged:                                                   │
│   u_n+1 = u                                                     │
│   state_n+1 = commit(state)                                     │
│                                                                 │
│ If not converged:                                               │
│   Report error                                                  │
│   Optionally: reduce step, try again                            │
└─────────────────────────────────────────────────────────────────┘
```

---

(继续 architecture-design-part2.md ...)


## 11. Migration Guide

- Progressive migration plan (12 weeks): `docs/design/migration-roadmap.md`
- Integrated migration chapter: `docs/design/architecture-design-part2.md`

---

