# HSTAR Refactoring Guide - Part 4: Migration Strategy
# HSTAR重构指南 - 第四部分：迁移策略与步骤

## Phase 1: Foundation (基础层) - 2-3 weeks

### Step 1.1: Create New Project Structure
```bash
mkdir -p HSTAR-Modern/{src/{core,mesh,elements,materials,boundary,solver,analysis,io,app},tests/{unit,integration,benchmarks},docs,examples}
```

### Step 1.2: Implement Core Modules
```
优先级: ★★★★★ (最高)
依赖: 无
可测试性: 独立测试

实现内容:
├── precision_mod.f90      # 1天
├── tensor_mod.f90         # 2天
├── sparse_matrix_mod.f90  # 3天
└── vector_mod.f90         # 1天

验收标准:
- 所有单元测试通过
- 与原代码数值结果一致
- 文档完整
```

### Step 1.3: Extract and Test Gauss Integration
```fortran
! 从 Elements.f90 提取高斯积分
module gauss_quadrature_mod
    use precision_mod
    implicit none
    private

    type, public :: gauss_rule
        integer(IP) :: n_points
        real(WP), allocatable :: points(:,:)   ! 积分点坐标
        real(WP), allocatable :: weights(:)    ! 权重
    contains
        procedure :: initialize => gauss_init
    end type

    ! 预定义规则
    type(gauss_rule), public, save :: GAUSS_LINE_2, GAUSS_LINE_3
    type(gauss_rule), public, save :: GAUSS_QUAD_2x2, GAUSS_QUAD_3x3
    type(gauss_rule), public, save :: GAUSS_TRI_1, GAUSS_TRI_3
    type(gauss_rule), public, save :: GAUSS_TET_1, GAUSS_TET_4
    type(gauss_rule), public, save :: GAUSS_HEX_2x2x2, GAUSS_HEX_3x3x3

contains
    subroutine initialize_all_rules()
        ! 一次性初始化所有预定义规则
        call GAUSS_LINE_2%initialize('LINE', 2)
        call GAUSS_QUAD_2x2%initialize('QUAD', 4)
        ! ...
    end subroutine
end module
```

---

## Phase 2: Mesh Layer (网格层) - 2 weeks

### Step 2.1: Node and Coordinate Module
```fortran
module node_mod
    use precision_mod
    implicit none
    private

    type, public :: node
        integer(IP) :: id
        real(WP) :: coords(3)   ! 始终使用3D坐标
        integer(IP), allocatable :: dof_ids(:)
    end type

    type, public :: node_set
        character(len=32) :: name
        integer(IP), allocatable :: node_ids(:)
    end type
end module
```

### Step 2.2: Element Factory Pattern
```fortran
module element_factory_mod
    use precision_mod
    use element_interface_mod
    ! 导入所有具体单元类型
    use line2_element_mod
    use tri3_element_mod
    use quad4_element_mod
    use tet4_element_mod
    use hex8_element_mod
    implicit none
    private

    public :: create_element

contains

    function create_element(element_type, n_nodes) result(elem)
        character(len=*), intent(in) :: element_type
        integer(IP), intent(in) :: n_nodes
        class(element_base), allocatable :: elem

        select case (trim(element_type))
        case ('LINE2', 'L2')
            allocate(line2_element :: elem)

        case ('TRI3', 'T3')
            allocate(tri3_element :: elem)

        case ('QUAD4', 'Q4')
            allocate(quad4_element :: elem)

        case ('TET4', 'H4')
            allocate(tet4_element :: elem)

        case ('HEX8', 'B8')
            allocate(hex8_element :: elem)

        case default
            error stop 'Unknown element type: ' // trim(element_type)
        end select
    end function

end module element_factory_mod
```

### Step 2.3: Mesh I/O (支持多种格式)
```fortran
module mesh_io_mod
    use precision_mod
    use mesh_mod
    implicit none
    private

    public :: read_mesh, write_mesh

    interface read_mesh
        module procedure read_mesh_native
        module procedure read_mesh_gmsh
        module procedure read_mesh_vtk
    end interface

contains

    subroutine read_mesh_gmsh(filename, msh, ierr)
        character(len=*), intent(in) :: filename
        type(mesh), intent(out) :: msh
        integer(IP), intent(out) :: ierr

        integer(IP) :: unit_num, n_nodes, n_elements
        character(len=256) :: line

        ! GMSH格式解析
        open(newunit=unit_num, file=filename, status='old', iostat=ierr)
        if (ierr /= 0) return

        ! 读取$MeshFormat
        ! 读取$Nodes
        ! 读取$Elements
        ! ...

        close(unit_num)
    end subroutine

end module mesh_io_mod
```

---

## Phase 3: Material Layer (材料层) - 3 weeks

### Step 3.1: Material Factory
```fortran
module material_factory_mod
    use precision_mod
    use material_interface_mod
    use elastic_material_mod
    use mohr_coulomb_mod
    use drucker_prager_mod
    use duncan_chang_mod
    implicit none
    private

    public :: create_material, material_registry

    type :: material_registry_type
        class(material_base), allocatable :: materials(:)
        integer(IP) :: count = 0
    contains
        procedure :: add => registry_add
        procedure :: get => registry_get
        procedure :: get_by_name => registry_get_by_name
    end type

    type(material_registry_type), save :: material_registry

contains

    function create_material(material_type, params) result(mat)
        character(len=*), intent(in) :: material_type
        real(WP), intent(in) :: params(:)
        class(material_base), allocatable :: mat

        select case (trim(material_type))
        case ('ELASTIC', 'ELASTIC_ISOTROPIC')
            allocate(elastic_isotropic :: mat)
            select type (mat)
            type is (elastic_isotropic)
                call mat%initialize(E=params(1), nu=params(2))
            end select

        case ('MOHR_COULOMB', 'MC')
            allocate(mohr_coulomb :: mat)
            select type (mat)
            type is (mohr_coulomb)
                call mat%initialize(E=params(1), nu=params(2), &
                    cohesion=params(3), phi=params(4), psi=params(5))
            end select

        ! ... 其他材料类型
        end select
    end function

end module material_factory_mod
```

### Step 3.2: State Variables Manager
```fortran
module state_manager_mod
    use precision_mod
    use material_interface_mod
    implicit none
    private

    public :: state_manager

    type :: state_manager_type
        type(material_state), allocatable :: current(:,:)    ! (n_gauss, n_elem)
        type(material_state), allocatable :: previous(:,:)
        type(material_state), allocatable :: converged(:,:)
    contains
        procedure :: initialize => sm_init
        procedure :: save_state => sm_save
        procedure :: restore_state => sm_restore
        procedure :: commit => sm_commit
        procedure :: get_state => sm_get
        procedure :: set_state => sm_set
    end type

    type(state_manager_type), save :: state_manager

contains

    subroutine sm_commit(self)
        ! 收敛后将当前状态提交为已收敛状态
        class(state_manager_type), intent(inout) :: self
        self%converged = self%current
    end subroutine

    subroutine sm_restore(self)
        ! 迭代不收敛时恢复到上一收敛状态
        class(state_manager_type), intent(inout) :: self
        self%current = self%converged
    end subroutine

end module state_manager_mod
```

---

## Phase 4: Solver Layer (求解器层) - 2 weeks

### Step 4.1: Abstract Solver with PARDISO Implementation
```fortran
module pardiso_solver_mod
    use precision_mod
    use solver_interface_mod
    use sparse_matrix_mod
    implicit none
    private

    type, extends(solver_base), public :: pardiso_solver
        integer(IP) :: pt(64)           ! PARDISO内部指针
        integer(IP) :: iparm(64)        ! PARDISO参数
        integer(IP) :: mtype = 11       ! 矩阵类型: 11=非对称
        integer(IP) :: phase = 0        ! 当前阶段
        logical :: factorized = .false.
    contains
        procedure :: setup => pardiso_setup
        procedure :: solve => pardiso_solve
        procedure :: cleanup => pardiso_cleanup
        procedure :: factorize => pardiso_factorize
    end type

contains

    subroutine pardiso_setup(self, K)
        class(pardiso_solver), intent(inout) :: self
        type(csr_matrix), intent(in) :: K

        ! 初始化PARDISO
        self%pt = 0
        call pardisoinit(self%pt, self%mtype, self%iparm)

        ! 设置参数
        self%iparm(1) = 1   ! 不使用默认值
        self%iparm(2) = 3   ! 并行嵌套剖分
        self%iparm(10) = 13 ! 扰动
        self%iparm(11) = 1  ! 缩放
        self%iparm(13) = 1  ! 改进精度

        self%n_equations = K%n_rows
        self%is_initialized = .true.
    end subroutine

    subroutine pardiso_solve(self, K, b, x, info)
        class(pardiso_solver), intent(inout) :: self
        type(csr_matrix), intent(in) :: K
        real(WP), intent(in) :: b(:)
        real(WP), intent(inout) :: x(:)
        integer(IP), intent(out) :: info

        integer(IP) :: nrhs, maxfct, mnum, msglvl, error
        real(WP) :: ddum

        nrhs = 1
        maxfct = 1
        mnum = 1
        msglvl = 0

        ! 分析 + 数值分解 + 求解
        self%phase = 13

        call pardiso(self%pt, maxfct, mnum, self%mtype, self%phase, &
            self%n_equations, K%values, K%row_ptr, K%col_idx, &
            idum, nrhs, self%iparm, msglvl, b, x, error)

        info = error
    end subroutine

end module pardiso_solver_mod
```

### Step 4.2: Newton-Raphson Controller
```fortran
module newton_raphson_mod
    use precision_mod
    use solver_interface_mod
    use sparse_matrix_mod
    implicit none
    private

    public :: newton_raphson_solver

    type :: newton_raphson_solver
        class(solver_base), allocatable :: linear_solver
        integer(IP) :: max_iterations = 25
        real(WP) :: force_tolerance = 1.0e-6_WP
        real(WP) :: displacement_tolerance = 1.0e-8_WP
        integer(IP) :: current_iteration = 0
        real(WP) :: current_residual = 0.0_WP
        logical :: converged = .false.
    contains
        procedure :: initialize => nr_init
        procedure :: iterate => nr_iterate
        procedure :: check_convergence => nr_check
        procedure :: report => nr_report
    end type

contains

    subroutine nr_iterate(self, K, R, du, info)
        class(newton_raphson_solver), intent(inout) :: self
        type(csr_matrix), intent(in) :: K    ! 切线刚度
        real(WP), intent(in) :: R(:)         ! 残差向量
        real(WP), intent(out) :: du(:)       ! 位移增量
        integer(IP), intent(out) :: info

        self%current_iteration = self%current_iteration + 1

        ! 求解线性系统 K * du = R
        call self%linear_solver%solve(K, R, du, info)

        ! 计算残差范数
        self%current_residual = norm2(R)

        ! 检查收敛
        call self%check_convergence()
    end subroutine

    subroutine nr_check(self)
        class(newton_raphson_solver), intent(inout) :: self

        if (self%current_residual < self%force_tolerance) then
            self%converged = .true.
        else if (self%current_iteration >= self%max_iterations) then
            self%converged = .false.
            print '(A,I0,A)', 'Warning: Newton-Raphson did not converge after ', &
                self%max_iterations, ' iterations'
        end if
    end subroutine

end module newton_raphson_mod
```

---

## Phase 5: Integration (集成) - 2 weeks

### Step 5.1: FEM Engine
```fortran
module fem_engine_mod
    use precision_mod
    use mesh_mod
    use material_factory_mod
    use element_factory_mod
    use solver_interface_mod
    use newton_raphson_mod
    use sparse_matrix_mod
    use state_manager_mod
    implicit none
    private

    public :: fem_engine

    type :: fem_engine
        type(mesh) :: msh
        type(csr_matrix) :: K_global      ! 全局刚度矩阵
        real(WP), allocatable :: F_ext(:) ! 外力向量
        real(WP), allocatable :: F_int(:) ! 内力向量
        real(WP), allocatable :: u(:)     ! 位移向量
        real(WP), allocatable :: du(:)    ! 位移增量
        type(newton_raphson_solver) :: nr_solver
    contains
        procedure :: initialize => engine_init
        procedure :: assemble_stiffness => engine_assemble_K
        procedure :: assemble_internal_force => engine_assemble_Fint
        procedure :: apply_bc => engine_apply_bc
        procedure :: solve_step => engine_solve_step
        procedure :: update_state => engine_update
    end type

contains

    subroutine engine_solve_step(self, load_factor)
        class(fem_engine), intent(inout) :: self
        real(WP), intent(in) :: load_factor

        real(WP), allocatable :: R(:)  ! 残差
        integer(IP) :: info

        allocate(R(self%msh%total_dof))

        ! Newton-Raphson迭代
        self%nr_solver%current_iteration = 0
        self%nr_solver%converged = .false.

        do while (.not. self%nr_solver%converged)
            ! 1. 组装刚度矩阵
            call self%assemble_stiffness()

            ! 2. 计算内力
            call self%assemble_internal_force()

            ! 3. 计算残差 R = F_ext - F_int
            R = load_factor * self%F_ext - self%F_int

            ! 4. 应用边界条件
            call self%apply_bc(self%K_global, R)

            ! 5. 求解位移增量
            call self%nr_solver%iterate(self%K_global, R, self%du, info)

            ! 6. 更新位移和状态
            self%u = self%u + self%du
            call self%update_state()
        end do

        ! 提交收敛状态
        call state_manager%commit()

        deallocate(R)
    end subroutine

end module fem_engine_mod
```

### Step 5.2: Analysis Controller
```fortran
module static_analysis_mod
    use precision_mod
    use fem_engine_mod
    use input_config_mod
    use output_writer_mod
    implicit none
    private

    public :: static_analysis

    type :: static_analysis
        type(fem_engine) :: engine
        type(analysis_config) :: config
        class(output_writer_base), allocatable :: output
        integer(IP) :: current_step = 0
        real(WP) :: current_time = 0.0_WP
    contains
        procedure :: initialize => analysis_init
        procedure :: run => analysis_run
        procedure :: step => analysis_step
        procedure :: finalize => analysis_finalize
    end type

contains

    subroutine analysis_run(self)
        class(static_analysis), intent(inout) :: self
        integer(IP) :: istep
        real(WP) :: load_factor

        do istep = 1, self%config%n_load_blocks
            load_factor = real(istep, WP) / real(self%config%n_load_blocks, WP)

            call self%step(load_factor)

            ! 输出结果
            if (mod(istep, self%config%output_frequency) == 0) then
                call self%output%write_results(istep, self%current_time, &
                    self%engine%u)
            end if
        end do
    end subroutine

end module static_analysis_mod
```

---

## Phase 6: Migration of Legacy Code (遗留代码迁移) - 4+ weeks

### Gradual Migration Strategy
```
┌────────────────────────────────────────────────────┐
│                  Migration Steps                    │
├────────────────────────────────────────────────────┤
│                                                    │
│  Week 1-2: 创建适配层(Adapter)                      │
│  ├── 旧代码调用新模块的桥接接口                      │
│  └── 保持旧代码可运行                               │
│                                                    │
│  Week 3-4: 迁移简单单元                             │
│  ├── L2, T3, Q4 (线性单元)                         │
│  └── 单元测试验证                                   │
│                                                    │
│  Week 5-6: 迁移弹性材料                             │
│  ├── ELASTIC_ISOTROPIC                            │
│  └── 基准测试验证                                   │
│                                                    │
│  Week 7-8: 迁移塑性材料                             │
│  ├── MOHR_COULOMB, DRUCKER                        │
│  └── 复杂案例验证                                   │
│                                                    │
│  Week 9+: 迁移特殊功能                              │
│  ├── 接触、蠕变、温度场                             │
│  └── 耦合分析                                      │
│                                                    │
└────────────────────────────────────────────────────┘
```

### Adapter Pattern for Gradual Migration
```fortran
module legacy_adapter_mod
    ! 允许新代码调用旧代码，或旧代码调用新代码
    use precision_mod
    use global_var           ! 旧的全局变量
    use material_interface_mod  ! 新的材料接口
    implicit none

contains

    subroutine call_legacy_dmatrix(matno, strain, stress, dmatx)
        ! 从新代码调用旧的dmatrix计算
        integer(IP), intent(in) :: matno
        real(WP), intent(in) :: strain(:)
        real(WP), intent(out) :: stress(:), dmatx(:,:)

        ! 设置旧代码需要的全局变量
        ! ... 调用旧的 dmatrix 子程序
    end subroutine

    subroutine wrap_new_material(new_mat, old_matno, old_strain, old_dmatx)
        ! 从旧代码调用新的材料模块
        class(material_base), intent(in) :: new_mat
        integer, intent(in) :: old_matno
        real(irk), intent(in) :: old_strain(:)
        real(irk), intent(out) :: old_dmatx(:,:)

        type(material_state) :: state
        real(WP) :: D(size(old_dmatx,1), size(old_dmatx,2))

        call new_mat%calc_tangent_modulus(state, D)
        old_dmatx = D
    end subroutine

end module legacy_adapter_mod
```

---

## Naming Conventions (命名规范)

### Variables (变量)
```fortran
! ✓ Good - Descriptive names
integer(IP) :: n_nodes, n_elements, n_gauss_points
real(WP) :: young_modulus, poisson_ratio, cohesion
real(WP) :: displacement(:), stress(:,:), strain(:,:)
type(mesh) :: current_mesh
type(material_state) :: gauss_point_state

! ✗ Bad - Cryptic abbreviations (current code)
integer :: npoin, nelem, ngaus
real :: e, nu, c
real :: deltafi(:), gpvar(:,:)
```

### Procedures (子程序)
```fortran
! ✓ Good - verb_noun pattern
subroutine calculate_stiffness_matrix(...)
subroutine assemble_global_matrix(...)
subroutine apply_boundary_conditions(...)
function compute_strain_from_displacement(...) result(strain)

! ✗ Bad - Unclear names
subroutine STIFF_U(...)
subroutine RESIDU_F(...)
```

### Types (类型)
```fortran
! ✓ Good - Descriptive with suffix
type :: element_base         ! 基类用 _base
type :: elastic_isotropic    ! 具体类型用描述性名称
type :: analysis_config      ! 配置类型用 _config
type :: material_state       ! 状态类型用 _state

! ✗ Bad - Unclear
type :: element_lib
type :: material_1, material_2, ...
type :: property_solid
```

### Modules (模块)
```fortran
! ✓ Good - Ends with _mod
module precision_mod
module elastic_material_mod
module sparse_matrix_mod

! ✗ Bad
module variable_types
module global_var
module elements
```
