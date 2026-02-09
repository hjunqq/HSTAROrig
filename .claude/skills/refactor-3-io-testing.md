# HSTAR Refactoring Guide - Part 3: I/O and Testing
# HSTAR重构指南 - 第三部分：输入输出与测试

## 1. Input System Refactoring (输入系统重构)

### Current Problem
```fortran
! 输入代码分散在各处，格式硬编码
open(inpunit,file='inp')
read (inpunit,*) text
read (inpunit,*) restart,relis,sysrelis,ADINA,Uopt_R,gamamax
read (inpunit,*) text
read (inpunit,*) probn

! 魔法数字作为文件单元号
integer (ink) gunit,cunit,eunit,punit,loadunit,munit,...
```

### Refactored Design

#### Input Configuration Module
```fortran
module input_config_mod
    use precision_mod
    implicit none
    private

    !===============================================
    ! 分析配置类型
    !===============================================
    type, public :: analysis_config
        ! 问题识别
        character(len=256) :: problem_name = ''
        character(len=32) :: analysis_type = 'STATIC'  ! STATIC/DYNAMIC/THERMAL

        ! 求解控制
        integer(IP) :: n_load_blocks = 1
        integer(IP) :: max_iterations = 20
        real(WP) :: convergence_tol = 1.0e-6_WP

        ! 时间控制
        real(WP) :: total_time = 0.0_WP
        real(WP) :: time_step = 0.0_WP
        character(len=16) :: time_integration = 'NEWMARK'

        ! 求解器选项
        character(len=16) :: solver_type = 'PARDISO'

        ! 输出控制
        integer(IP) :: output_frequency = 1
        logical :: output_displacement = .true.
        logical :: output_stress = .true.
        logical :: output_strain = .false.

        ! 重启控制
        logical :: is_restart = .false.
        integer(IP) :: restart_step = 0
    contains
        procedure :: validate => config_validate
        procedure :: print_summary => config_print
    end type

end module input_config_mod
```

#### Input Parser Module
```fortran
module input_parser_mod
    use precision_mod
    use input_config_mod
    use mesh_mod
    use material_interface_mod
    implicit none
    private

    public :: input_parser

    !===============================================
    ! 输入解析器类型
    !===============================================
    type :: input_parser
        character(len=512) :: input_file = ''
        integer(IP) :: current_line = 0
        integer(IP) :: unit_number = -1
        logical :: is_open = .false.
    contains
        procedure :: open_file => parser_open
        procedure :: close_file => parser_close
        procedure :: read_config => parser_read_config
        procedure :: read_mesh => parser_read_mesh
        procedure :: read_materials => parser_read_materials
        procedure :: read_boundary_conditions => parser_read_bc
        procedure :: read_loads => parser_read_loads

        ! 辅助方法
        procedure, private :: skip_comments => parser_skip_comments
        procedure, private :: read_line => parser_read_line
        procedure, private :: parse_keyword => parser_parse_keyword
        procedure, private :: report_error => parser_error
    end type

contains

    subroutine parser_read_config(self, config, ierr)
        class(input_parser), intent(inout) :: self
        type(analysis_config), intent(out) :: config
        integer(IP), intent(out) :: ierr

        character(len=256) :: line, keyword, value
        logical :: found_end

        ierr = 0
        found_end = .false.

        do while (.not. found_end)
            call self%read_line(line, ierr)
            if (ierr /= 0) exit

            call self%parse_keyword(line, keyword, value)

            select case (trim(keyword))
            case ('PROBLEM_NAME')
                config%problem_name = trim(value)

            case ('ANALYSIS_TYPE')
                config%analysis_type = trim(value)

            case ('MAX_ITERATIONS')
                read(value, *, iostat=ierr) config%max_iterations

            case ('CONVERGENCE_TOL')
                read(value, *, iostat=ierr) config%convergence_tol

            case ('SOLVER')
                config%solver_type = trim(value)

            case ('END_CONFIG')
                found_end = .true.

            case default
                call self%report_error('Unknown keyword: ' // trim(keyword))
            end select
        end do
    end subroutine

    ! ... 其他解析方法

end module input_parser_mod
```

#### Modern Input File Format
```
# HSTAR Modern Input File Format
# Lines starting with # are comments

*PROBLEM
  NAME = dam_analysis
  TYPE = STATIC
  DIMENSION = 3
*END

*CONTROL
  MAX_ITERATIONS = 25
  TOLERANCE = 1.0E-8
  SOLVER = PARDISO
*END

*MESH_FILE
  FORMAT = MSH
  PATH = mesh/dam.msh
*END

*MATERIAL, ID=1, NAME=Concrete
  TYPE = ELASTIC
  E = 30.0E9
  NU = 0.2
  DENSITY = 2400.0
*END

*MATERIAL, ID=2, NAME=Rock
  TYPE = MOHR_COULOMB
  E = 20.0E9
  NU = 0.25
  COHESION = 1.0E6
  FRICTION_ANGLE = 35.0
  DILATION_ANGLE = 10.0
*END

*BOUNDARY
  TYPE = FIXED
  NODE_SET = bottom
  DOF = ALL
*END

*LOAD_BLOCK, ID=1
  *INCREMENT, N=10, DT=1.0
    *GRAVITY
      DIRECTION = 0.0, 0.0, -9.81
      GROUPS = ALL
    *END
  *END
*END

*OUTPUT
  FREQUENCY = 1
  VARIABLES = DISPLACEMENT, STRESS, STRAIN
  FORMAT = GID
*END
```

---

## 2. Output System Refactoring (输出系统重构)

### Output Writer Module
```fortran
module output_writer_mod
    use precision_mod
    use mesh_mod
    implicit none
    private

    !===============================================
    ! 输出格式枚举
    !===============================================
    integer(IP), parameter, public :: OUTPUT_FORMAT_TEXT = 1
    integer(IP), parameter, public :: OUTPUT_FORMAT_BINARY = 2
    integer(IP), parameter, public :: OUTPUT_FORMAT_GID = 3
    integer(IP), parameter, public :: OUTPUT_FORMAT_VTK = 4

    !===============================================
    ! 输出器抽象类型
    !===============================================
    type, abstract, public :: output_writer_base
        character(len=256) :: output_path = ''
        integer(IP) :: current_step = 0
    contains
        procedure(init_interface), deferred :: initialize
        procedure(write_mesh_interface), deferred :: write_mesh
        procedure(write_results_interface), deferred :: write_results
        procedure(finalize_interface), deferred :: finalize
    end type

    !===============================================
    ! GiD输出器
    !===============================================
    type, extends(output_writer_base), public :: gid_writer
        integer(IP) :: mesh_unit = -1
        integer(IP) :: result_unit = -1
    contains
        procedure :: initialize => gid_init
        procedure :: write_mesh => gid_write_mesh
        procedure :: write_results => gid_write_results
        procedure :: finalize => gid_finalize
    end type

    !===============================================
    ! VTK输出器 (现代可视化)
    !===============================================
    type, extends(output_writer_base), public :: vtk_writer
        logical :: binary_output = .true.
    contains
        procedure :: initialize => vtk_init
        procedure :: write_mesh => vtk_write_mesh
        procedure :: write_results => vtk_write_results
        procedure :: finalize => vtk_finalize
    end type

contains

    subroutine gid_write_results(self, step, time, displacement, stress)
        class(gid_writer), intent(inout) :: self
        integer(IP), intent(in) :: step
        real(WP), intent(in) :: time
        real(WP), intent(in) :: displacement(:,:)
        real(WP), intent(in), optional :: stress(:,:)

        ! 写入位移结果
        write(self%result_unit, '(A)') 'Result "Displacement" "Analysis" ' // &
            trim(adjustl(int_to_str(step))) // ' Vector OnNodes'
        write(self%result_unit, '(A)') 'Values'
        ! ... 写入数据
        write(self%result_unit, '(A)') 'End Values'

        ! 写入应力结果
        if (present(stress)) then
            write(self%result_unit, '(A)') 'Result "Stress" "Analysis" ' // &
                trim(adjustl(int_to_str(step))) // ' Matrix OnGaussPoints'
            ! ...
        end if
    end subroutine

end module output_writer_mod
```

---

## 3. Testing Framework (测试框架)

### Unit Test Module
```fortran
module test_framework_mod
    use precision_mod
    implicit none
    private

    public :: test_suite, test_case, run_tests
    public :: assert_equal, assert_near, assert_true, assert_false

    !===============================================
    ! 测试用例类型
    !===============================================
    type :: test_case
        character(len=128) :: name = ''
        logical :: passed = .false.
        character(len=256) :: message = ''
        real(WP) :: elapsed_time = 0.0_WP
    end type

    !===============================================
    ! 测试套件类型
    !===============================================
    type :: test_suite
        character(len=128) :: name = ''
        type(test_case), allocatable :: tests(:)
        integer(IP) :: n_tests = 0
        integer(IP) :: n_passed = 0
        integer(IP) :: n_failed = 0
    contains
        procedure :: add_test => suite_add_test
        procedure :: run => suite_run
        procedure :: report => suite_report
    end type

contains

    !===============================================
    ! 断言函数
    !===============================================
    subroutine assert_equal(actual, expected, test_name, passed, message)
        class(*), intent(in) :: actual, expected
        character(len=*), intent(in) :: test_name
        logical, intent(out) :: passed
        character(len=*), intent(out) :: message

        passed = .false.
        message = ''

        select type (actual)
        type is (integer(IP))
            select type (expected)
            type is (integer(IP))
                passed = (actual == expected)
                if (.not. passed) then
                    write(message, '(A,I0,A,I0)') 'Expected ', expected, ', got ', actual
                end if
            end select

        type is (real(WP))
            select type (expected)
            type is (real(WP))
                passed = (abs(actual - expected) < TOLERANCE)
                if (.not. passed) then
                    write(message, '(A,ES12.5,A,ES12.5)') 'Expected ', expected, ', got ', actual
                end if
            end select
        end select
    end subroutine

    subroutine assert_near(actual, expected, tolerance, test_name, passed, message)
        real(WP), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: test_name
        logical, intent(out) :: passed
        character(len=*), intent(out) :: message

        passed = abs(actual - expected) <= tolerance
        if (.not. passed) then
            write(message, '(A,ES12.5,A,ES12.5,A,ES12.5)') &
                'Expected ', expected, ' +/- ', tolerance, ', got ', actual
        end if
    end subroutine

    subroutine suite_report(self)
        class(test_suite), intent(in) :: self
        integer(IP) :: i

        print '(A)', repeat('=', 60)
        print '(A,A)', 'Test Suite: ', trim(self%name)
        print '(A)', repeat('=', 60)

        do i = 1, self%n_tests
            if (self%tests(i)%passed) then
                print '(A,A,A)', '[PASS] ', trim(self%tests(i)%name)
            else
                print '(A,A,A)', '[FAIL] ', trim(self%tests(i)%name)
                print '(A,A)', '       ', trim(self%tests(i)%message)
            end if
        end do

        print '(A)', repeat('-', 60)
        print '(A,I0,A,I0,A,I0)', 'Total: ', self%n_tests, &
            ' | Passed: ', self%n_passed, ' | Failed: ', self%n_failed
        print '(A)', repeat('=', 60)
    end subroutine

end module test_framework_mod
```

### Example Test File
```fortran
program test_elastic_material
    use precision_mod
    use test_framework_mod
    use elastic_material_mod
    use material_interface_mod
    implicit none

    type(test_suite) :: suite
    type(elastic_isotropic) :: mat
    type(material_state) :: state
    real(WP) :: D(6,6), expected_D(6,6)
    real(WP) :: E, nu, lambda, mu
    logical :: passed
    character(len=256) :: msg

    suite%name = 'Elastic Material Tests'

    ! 初始化材料
    E = 200.0e9_WP  ! 200 GPa
    nu = 0.3_WP
    call mat%initialize(E, nu)

    !-------------------------------------------------
    ! Test 1: Check elastic constants
    !-------------------------------------------------
    call assert_near(mat%E, E, 1.0_WP, 'Elastic modulus', passed, msg)
    call suite%add_test('Elastic modulus storage', passed, msg)

    call assert_near(mat%nu, nu, 1.0e-10_WP, 'Poisson ratio', passed, msg)
    call suite%add_test('Poisson ratio storage', passed, msg)

    !-------------------------------------------------
    ! Test 2: Check D matrix symmetry
    !-------------------------------------------------
    call state%reset(6)
    call mat%calc_tangent_modulus(state, D)

    passed = all(abs(D - transpose(D)) < TOLERANCE)
    msg = 'D matrix is not symmetric'
    call suite%add_test('D matrix symmetry', passed, msg)

    !-------------------------------------------------
    ! Test 3: Check D matrix values
    !-------------------------------------------------
    lambda = E * nu / ((1.0_WP + nu) * (1.0_WP - 2.0_WP * nu))
    mu = E / (2.0_WP * (1.0_WP + nu))

    call assert_near(D(1,1), lambda + 2.0_WP * mu, 1.0_WP, 'D(1,1)', passed, msg)
    call suite%add_test('D matrix diagonal value', passed, msg)

    !-------------------------------------------------
    ! Report results
    !-------------------------------------------------
    call suite%report()

    if (suite%n_failed > 0) then
        stop 1
    end if

end program test_elastic_material
```

---

## 4. Benchmark Tests (基准测试)

### Patch Test
```fortran
module patch_test_mod
    use precision_mod
    use mesh_mod
    use element_interface_mod
    use elastic_material_mod
    implicit none

contains

    subroutine run_patch_test_quad4(passed, message)
        ! 四边形单元片检验 - 线性位移场
        logical, intent(out) :: passed
        character(len=*), intent(out) :: message

        type(mesh) :: test_mesh
        real(WP) :: coords(2,4), displacement(8)
        real(WP) :: u_exact(8), error

        ! 创建单个四边形单元
        coords(:,1) = [0.0_WP, 0.0_WP]
        coords(:,2) = [1.0_WP, 0.0_WP]
        coords(:,3) = [1.0_WP, 1.0_WP]
        coords(:,4) = [0.0_WP, 1.0_WP]

        ! 施加线性位移场 u = a*x + b*y
        ! 对于线性位移场，所有单元应该精确通过
        u_exact(1) = 0.0_WP; u_exact(2) = 0.0_WP  ! node 1
        u_exact(3) = 0.001_WP; u_exact(4) = 0.0_WP  ! node 2
        u_exact(5) = 0.001_WP; u_exact(6) = 0.001_WP  ! node 3
        u_exact(7) = 0.0_WP; u_exact(8) = 0.001_WP  ! node 4

        ! 运行分析并比较结果
        ! ... (详细实现)

        error = maxval(abs(displacement - u_exact))
        passed = error < 1.0e-12_WP

        if (.not. passed) then
            write(message, '(A,ES12.5)') 'Patch test failed with error: ', error
        else
            message = 'Patch test passed'
        end if
    end subroutine

end module patch_test_mod
```

### Benchmark Against Analytical Solution
```fortran
module benchmark_tests_mod
    use precision_mod
    implicit none

contains

    subroutine benchmark_cantilever_beam(passed, rel_error)
        ! 悬臂梁弯曲 - 与解析解比较
        logical, intent(out) :: passed
        real(WP), intent(out) :: rel_error

        real(WP) :: L, h, E, nu, P  ! 梁参数
        real(WP) :: I, delta_analytical, delta_fem

        L = 1.0_WP      ! 长度
        h = 0.1_WP      ! 高度
        E = 200.0e9_WP  ! 弹性模量
        nu = 0.3_WP     ! 泊松比
        P = 1000.0_WP   ! 端部载荷

        I = h**3 / 12.0_WP  ! 惯性矩

        ! 解析解
        delta_analytical = P * L**3 / (3.0_WP * E * I)

        ! 有限元解 (调用FEM求解器)
        ! delta_fem = solve_cantilever_fem(L, h, E, nu, P)

        rel_error = abs(delta_fem - delta_analytical) / delta_analytical
        passed = rel_error < 0.01_WP  ! 1% tolerance
    end subroutine

end module benchmark_tests_mod
```

---

## 5. CMake Build System

```cmake
cmake_minimum_required(VERSION 3.16)
project(HSTAR VERSION 2.0 LANGUAGES Fortran)

# 编译器选项
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -check all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
endif()

# 查找依赖
find_package(MKL REQUIRED)
find_package(OpenMP REQUIRED)

# 核心库
add_library(hstar_core
    src/core/precision_mod.f90
    src/core/tensor_mod.f90
    src/core/sparse_matrix_mod.f90
)

# 网格库
add_library(hstar_mesh
    src/mesh/node_mod.f90
    src/mesh/element_interface.f90
    src/mesh/mesh_mod.f90
)
target_link_libraries(hstar_mesh PUBLIC hstar_core)

# 单元库
add_library(hstar_elements
    src/elements/line2_mod.f90
    src/elements/tri3_mod.f90
    src/elements/quad4_mod.f90
    src/elements/tet4_mod.f90
    src/elements/hex8_mod.f90
)
target_link_libraries(hstar_elements PUBLIC hstar_mesh)

# 材料库
add_library(hstar_materials
    src/materials/material_interface.f90
    src/materials/elastic_mod.f90
    src/materials/mohr_coulomb_mod.f90
    src/materials/drucker_prager_mod.f90
)
target_link_libraries(hstar_materials PUBLIC hstar_core)

# 求解器库
add_library(hstar_solver
    src/solver/solver_interface.f90
    src/solver/direct_solver_mod.f90
    src/solver/iterative_solver_mod.f90
)
target_link_libraries(hstar_solver PUBLIC hstar_core MKL::MKL)

# 主程序
add_executable(hstar src/app/main.f90)
target_link_libraries(hstar PRIVATE
    hstar_elements
    hstar_materials
    hstar_solver
    OpenMP::OpenMP_Fortran
)

# 测试
enable_testing()

add_executable(test_materials tests/unit/test_materials.f90)
target_link_libraries(test_materials PRIVATE hstar_materials)
add_test(NAME MaterialTests COMMAND test_materials)

add_executable(test_elements tests/unit/test_elements.f90)
target_link_libraries(test_elements PRIVATE hstar_elements)
add_test(NAME ElementTests COMMAND test_elements)

add_executable(test_benchmark tests/benchmarks/test_benchmark.f90)
target_link_libraries(test_benchmark PRIVATE hstar)
add_test(NAME BenchmarkTests COMMAND test_benchmark)
```
