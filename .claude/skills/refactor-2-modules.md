# HSTAR Refactoring Guide - Part 2: Core Module Patterns
# HSTAR重构指南 - 第二部分：核心模块设计模式

## 1. Precision Module (精度模块)

### Current (现有)
```fortran
module variable_types
    integer,parameter::ink=k_int,irk=k_double
end module
```

### Refactored (重构后)
```fortran
module precision_mod
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
    implicit none
    private

    ! 导出精度类型
    integer, parameter, public :: WP = real64     ! 工作精度
    integer, parameter, public :: IP = int32      ! 整数精度
    integer, parameter, public :: LP = int64      ! 长整数(索引)

    ! 常用数学常量
    real(WP), parameter, public :: PI = 3.141592653589793238_WP
    real(WP), parameter, public :: ZERO = 0.0_WP
    real(WP), parameter, public :: ONE = 1.0_WP
    real(WP), parameter, public :: TOLERANCE = 1.0e-10_WP

end module precision_mod
```

---

## 2. Element Interface (单元抽象接口)

### Current Problem
```fortran
! 26种单元类型硬编码在一个大文件里
call l2_define(elkn(1))
call t3_define(elkn(3))
...
! 每种单元的代码结构不统一
```

### Refactored Design
```fortran
module element_interface_mod
    use precision_mod
    implicit none
    private

    !===============================================
    ! 单元抽象类型 - 所有单元的基类
    !===============================================
    type, abstract, public :: element_base
        integer(IP) :: id              ! 单元编号
        integer(IP) :: n_nodes         ! 节点数
        integer(IP) :: n_dof           ! 自由度数
        integer(IP) :: n_gauss         ! 高斯点数
        integer(IP) :: dimension       ! 问题维度 (2D/3D)
        integer(IP), allocatable :: node_ids(:)     ! 节点列表
        integer(IP), allocatable :: dof_indices(:)  ! 自由度索引
    contains
        ! 必须实现的方法
        procedure(shape_func_interface), deferred :: calc_shape_functions
        procedure(b_matrix_interface), deferred :: calc_b_matrix
        procedure(stiffness_interface), deferred :: calc_stiffness
        procedure(internal_force_interface), deferred :: calc_internal_force

        ! 通用方法
        procedure :: get_gauss_points => element_get_gauss_points
        procedure :: get_jacobian => element_get_jacobian
    end type element_base

    !===============================================
    ! 抽象接口定义
    !===============================================
    abstract interface
        subroutine shape_func_interface(self, xi, eta, zeta, N, dN)
            import :: element_base, WP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: xi, eta, zeta  ! 自然坐标
            real(WP), intent(out) :: N(:)          ! 形函数值
            real(WP), intent(out) :: dN(:,:)       ! 形函数导数
        end subroutine

        subroutine b_matrix_interface(self, coords, igauss, B, detJ)
            import :: element_base, WP, IP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: coords(:,:)    ! 节点坐标
            integer(IP), intent(in) :: igauss     ! 高斯点索引
            real(WP), intent(out) :: B(:,:)        ! B矩阵
            real(WP), intent(out) :: detJ          ! 雅可比行列式
        end subroutine

        subroutine stiffness_interface(self, coords, D, Ke)
            import :: element_base, WP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: coords(:,:)    ! 节点坐标
            real(WP), intent(in) :: D(:,:)         ! 本构矩阵
            real(WP), intent(out) :: Ke(:,:)       ! 单元刚度矩阵
        end subroutine

        subroutine internal_force_interface(self, coords, stress, Fe)
            import :: element_base, WP
            class(element_base), intent(in) :: self
            real(WP), intent(in) :: coords(:,:)    ! 节点坐标
            real(WP), intent(in) :: stress(:,:)    ! 高斯点应力
            real(WP), intent(out) :: Fe(:)         ! 单元内力向量
        end subroutine
    end interface

end module element_interface_mod
```

### Concrete Element Example (具体单元实现)
```fortran
module quad4_element_mod
    use precision_mod
    use element_interface_mod
    implicit none
    private

    type, extends(element_base), public :: quad4_element
    contains
        procedure :: calc_shape_functions => quad4_shape
        procedure :: calc_b_matrix => quad4_b_matrix
        procedure :: calc_stiffness => quad4_stiffness
        procedure :: calc_internal_force => quad4_internal_force
    end type

contains

    subroutine quad4_shape(self, xi, eta, zeta, N, dN)
        class(quad4_element), intent(in) :: self
        real(WP), intent(in) :: xi, eta, zeta
        real(WP), intent(out) :: N(4), dN(4,2)

        ! 4节点四边形形函数
        N(1) = 0.25_WP * (1.0_WP - xi) * (1.0_WP - eta)
        N(2) = 0.25_WP * (1.0_WP + xi) * (1.0_WP - eta)
        N(3) = 0.25_WP * (1.0_WP + xi) * (1.0_WP + eta)
        N(4) = 0.25_WP * (1.0_WP - xi) * (1.0_WP + eta)

        ! 形函数对自然坐标的导数
        dN(1,1) = -0.25_WP * (1.0_WP - eta)  ! dN1/dxi
        dN(1,2) = -0.25_WP * (1.0_WP - xi)   ! dN1/deta
        ! ... 其他节点
    end subroutine

    ! ... 其他方法实现

end module quad4_element_mod
```

---

## 3. Material Interface (材料抽象接口)

### Current Problem
```fortran
! 材料类型通过字符串判断，本构计算分散在各处
if (material=='ELASTIC_ISOTROPIC') then
    ! ...
elseif (material=='MOHRCOLUMB') then
    ! ...
```

### Refactored Design
```fortran
module material_interface_mod
    use precision_mod
    implicit none
    private

    !===============================================
    ! 材料状态类型 - 存储高斯点的历史变量
    !===============================================
    type, public :: material_state
        real(WP), allocatable :: stress(:)      ! 应力 [nstre]
        real(WP), allocatable :: strain(:)      ! 应变 [nstre]
        real(WP), allocatable :: plastic_strain(:) ! 塑性应变
        real(WP), allocatable :: internal_vars(:)  ! 内变量
        logical :: is_plastic = .false.         ! 塑性状态标志
    contains
        procedure :: copy => state_copy
        procedure :: reset => state_reset
    end type

    !===============================================
    ! 材料抽象类型
    !===============================================
    type, abstract, public :: material_base
        integer(IP) :: id                       ! 材料编号
        character(len=32) :: name               ! 材料名称
        integer(IP) :: n_stress_components      ! 应力分量数
        real(WP) :: density = 0.0_WP           ! 密度
    contains
        ! 必须实现的方法
        procedure(tangent_interface), deferred :: calc_tangent_modulus
        procedure(stress_interface), deferred :: calc_stress
        procedure(update_interface), deferred :: update_state

        ! 通用方法
        procedure :: get_elastic_modulus => material_get_E
        procedure :: get_poisson_ratio => material_get_nu
    end type

    !===============================================
    ! 抽象接口
    !===============================================
    abstract interface
        subroutine tangent_interface(self, state, D)
            import :: material_base, material_state, WP
            class(material_base), intent(in) :: self
            type(material_state), intent(in) :: state
            real(WP), intent(out) :: D(:,:)     ! 切线刚度矩阵
        end subroutine

        subroutine stress_interface(self, strain, state, stress)
            import :: material_base, material_state, WP
            class(material_base), intent(in) :: self
            real(WP), intent(in) :: strain(:)
            type(material_state), intent(inout) :: state
            real(WP), intent(out) :: stress(:)
        end subroutine

        subroutine update_interface(self, d_strain, state)
            import :: material_base, material_state, WP
            class(material_base), intent(in) :: self
            real(WP), intent(in) :: d_strain(:)
            type(material_state), intent(inout) :: state
        end subroutine
    end interface

end module material_interface_mod
```

### Elastic Material Example
```fortran
module elastic_material_mod
    use precision_mod
    use material_interface_mod
    implicit none
    private

    type, extends(material_base), public :: elastic_isotropic
        real(WP) :: E       ! 弹性模量
        real(WP) :: nu      ! 泊松比
        real(WP) :: G       ! 剪切模量 (derived)
        real(WP) :: K       ! 体积模量 (derived)
    contains
        procedure :: calc_tangent_modulus => elastic_tangent
        procedure :: calc_stress => elastic_stress
        procedure :: update_state => elastic_update
        procedure :: initialize => elastic_init
    end type

contains

    subroutine elastic_init(self, E, nu)
        class(elastic_isotropic), intent(inout) :: self
        real(WP), intent(in) :: E, nu

        self%E = E
        self%nu = nu
        self%G = E / (2.0_WP * (1.0_WP + nu))
        self%K = E / (3.0_WP * (1.0_WP - 2.0_WP * nu))
        self%n_stress_components = 6  ! 3D case
    end subroutine

    subroutine elastic_tangent(self, state, D)
        class(elastic_isotropic), intent(in) :: self
        type(material_state), intent(in) :: state
        real(WP), intent(out) :: D(:,:)

        real(WP) :: lambda, mu

        lambda = self%E * self%nu / ((1.0_WP + self%nu) * (1.0_WP - 2.0_WP * self%nu))
        mu = self%G

        D = 0.0_WP
        D(1,1) = lambda + 2.0_WP * mu
        D(2,2) = lambda + 2.0_WP * mu
        D(3,3) = lambda + 2.0_WP * mu
        D(1,2) = lambda; D(2,1) = lambda
        D(1,3) = lambda; D(3,1) = lambda
        D(2,3) = lambda; D(3,2) = lambda
        D(4,4) = mu
        D(5,5) = mu
        D(6,6) = mu
    end subroutine

    ! ... 其他方法

end module elastic_material_mod
```

### Mohr-Coulomb Example
```fortran
module mohr_coulomb_mod
    use precision_mod
    use material_interface_mod
    use elastic_material_mod
    implicit none
    private

    type, extends(material_base), public :: mohr_coulomb
        ! 弹性参数
        real(WP) :: E, nu
        ! 强度参数
        real(WP) :: cohesion      ! 粘聚力 c
        real(WP) :: friction_angle ! 摩擦角 phi
        real(WP) :: dilation_angle ! 剪胀角 psi
        ! 派生参数
        real(WP) :: sin_phi, cos_phi
        real(WP) :: sin_psi, cos_psi
    contains
        procedure :: calc_tangent_modulus => mc_tangent
        procedure :: calc_stress => mc_stress
        procedure :: update_state => mc_update
        procedure :: yield_function => mc_yield
        procedure :: plastic_potential => mc_potential
    end type

contains

    pure function mc_yield(self, stress) result(f)
        class(mohr_coulomb), intent(in) :: self
        real(WP), intent(in) :: stress(:)
        real(WP) :: f

        real(WP) :: sigma1, sigma3, mean_stress, deviatoric

        ! 计算主应力
        ! ... 主应力计算

        ! Mohr-Coulomb屈服函数
        f = (sigma1 - sigma3) + (sigma1 + sigma3) * self%sin_phi &
            - 2.0_WP * self%cohesion * self%cos_phi
    end function

    ! ... 其他方法

end module mohr_coulomb_mod
```

---

## 4. Solver Interface (求解器接口)

```fortran
module solver_interface_mod
    use precision_mod
    use sparse_matrix_mod
    implicit none
    private

    !===============================================
    ! 求解器抽象类型
    !===============================================
    type, abstract, public :: solver_base
        integer(IP) :: n_equations           ! 方程数
        integer(IP) :: max_iterations = 1000 ! 最大迭代次数
        real(WP) :: tolerance = 1.0e-10_WP   ! 收敛容差
        logical :: is_initialized = .false.
    contains
        procedure(setup_interface), deferred :: setup
        procedure(solve_interface), deferred :: solve
        procedure(cleanup_interface), deferred :: cleanup
        procedure :: check_convergence => solver_check_convergence
    end type

    abstract interface
        subroutine setup_interface(self, K)
            import :: solver_base, csr_matrix
            class(solver_base), intent(inout) :: self
            type(csr_matrix), intent(in) :: K
        end subroutine

        subroutine solve_interface(self, K, b, x, info)
            import :: solver_base, csr_matrix, WP, IP
            class(solver_base), intent(inout) :: self
            type(csr_matrix), intent(in) :: K
            real(WP), intent(in) :: b(:)
            real(WP), intent(inout) :: x(:)
            integer(IP), intent(out) :: info
        end subroutine

        subroutine cleanup_interface(self)
            import :: solver_base
            class(solver_base), intent(inout) :: self
        end subroutine
    end interface

end module solver_interface_mod
```

---

## 5. Mesh Container (网格容器)

```fortran
module mesh_mod
    use precision_mod
    use node_mod
    use element_interface_mod
    implicit none
    private

    type, public :: mesh
        ! 基本信息
        integer(IP) :: n_nodes = 0
        integer(IP) :: n_elements = 0
        integer(IP) :: n_dimension = 3

        ! 节点数据
        real(WP), allocatable :: coordinates(:,:)  ! (ndim, nnodes)

        ! 单元数据 - 多态单元列表
        class(element_base), allocatable :: elements(:)

        ! 材料分配
        integer(IP), allocatable :: element_material(:)

        ! 节点自由度编号
        integer(IP), allocatable :: dof_map(:,:)   ! (ndof, nnodes)
        integer(IP) :: total_dof = 0

    contains
        procedure :: initialize => mesh_init
        procedure :: add_node => mesh_add_node
        procedure :: add_element => mesh_add_element
        procedure :: build_dof_map => mesh_build_dof
        procedure :: get_element_coords => mesh_get_elem_coords
        procedure :: cleanup => mesh_cleanup
    end type

end module mesh_mod
```
