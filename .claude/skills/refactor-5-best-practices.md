# HSTAR Refactoring Guide - Part 5: Modern Fortran Best Practices
# HSTAR重构指南 - 第五部分：现代Fortran最佳实践

## 1. Modern Fortran Features to Use (应使用的现代Fortran特性)

### 1.1 Modules and Submodules (模块与子模块)
```fortran
! 主模块定义接口
module material_mod
    implicit none
    private

    type, abstract, public :: material_base
    contains
        procedure(calc_stress_if), deferred :: calc_stress
    end type

    abstract interface
        subroutine calc_stress_if(self, strain, stress)
            import :: material_base
            class(material_base), intent(in) :: self
            real(8), intent(in) :: strain(:)
            real(8), intent(out) :: stress(:)
        end subroutine
    end interface

end module

! 子模块提供实现 (分离接口与实现)
submodule (material_mod) material_impl
contains
    ! 实现细节
end submodule
```

### 1.2 Object-Oriented Features (面向对象特性)
```fortran
! 继承与多态
type, extends(material_base) :: elastic_material
    real(8) :: E, nu
contains
    procedure :: calc_stress => elastic_calc_stress
    procedure :: calc_tangent => elastic_calc_tangent
    final :: elastic_cleanup  ! 析构函数
end type

! 类型绑定过程
type :: mesh
contains
    procedure :: initialize
    procedure :: add_node
    procedure :: add_element
    generic :: add => add_node, add_element  ! 泛型接口
end type
```

### 1.3 Pure and Elemental Functions (纯函数与逐元函数)
```fortran
! 纯函数 - 无副作用，可并行化
pure function von_mises_stress(stress) result(vm)
    real(8), intent(in) :: stress(6)
    real(8) :: vm
    real(8) :: s1, s2, s3, s12, s23, s13

    s1 = stress(1); s2 = stress(2); s3 = stress(3)
    s12 = stress(4); s23 = stress(5); s13 = stress(6)

    vm = sqrt(0.5d0 * ((s1-s2)**2 + (s2-s3)**2 + (s3-s1)**2 &
              + 6.0d0 * (s12**2 + s23**2 + s13**2)))
end function

! 逐元函数 - 可作用于数组
elemental function deg_to_rad(degrees) result(radians)
    real(8), intent(in) :: degrees
    real(8) :: radians
    real(8), parameter :: PI = 3.141592653589793d0
    radians = degrees * PI / 180.0d0
end function
```

### 1.4 Allocatable and Pointer (可分配与指针)
```fortran
! 优先使用allocatable而非pointer
type :: element
    real(8), allocatable :: stiffness(:,:)  ! ✓ 自动管理内存
    real(8), pointer :: legacy_ptr(:)        ! ✗ 需手动管理
end type

! 自动重分配
real(8), allocatable :: array(:)
array = [1.0d0, 2.0d0, 3.0d0]  ! 自动分配
array = [array, 4.0d0, 5.0d0]  ! 自动扩展
```

### 1.5 Assumed-Shape and Assumed-Rank (假定形状与假定秩)
```fortran
! 假定形状 - 灵活的数组参数
subroutine process_matrix(A)
    real(8), intent(inout) :: A(:,:)  ! 任意大小的2D数组
    integer :: n, m
    n = size(A, 1)
    m = size(A, 2)
end subroutine

! 假定秩 (Fortran 2018)
subroutine generic_print(array)
    real(8), intent(in) :: array(..)  ! 任意维度
    select rank(array)
    rank(1)
        print *, 'Vector:', array
    rank(2)
        print *, 'Matrix:', array
    end select
end subroutine
```

### 1.6 Coarrays for Parallel (协数组并行)
```fortran
! 并行分布式计算
program parallel_fem
    real(8), allocatable :: local_stiffness(:,:)[:]  ! 协数组
    integer :: me, num_images

    me = this_image()
    num_images = num_images()

    allocate(local_stiffness(100,100)[*])

    ! 每个镜像计算自己的部分
    call compute_local_stiffness(me, local_stiffness)

    ! 同步
    sync all

    ! 收集结果
    if (me == 1) then
        do i = 2, num_images
            ! 从其他镜像获取数据
            global_stiffness = global_stiffness + local_stiffness(:,:)[i]
        end do
    end if
end program
```

---

## 2. Error Handling Patterns (错误处理模式)

### 2.1 Error Codes with Message
```fortran
module error_handling_mod
    implicit none
    private

    public :: error_info, check_error, report_error

    ! 错误代码
    integer, parameter, public :: ERR_NONE = 0
    integer, parameter, public :: ERR_MEMORY = 1
    integer, parameter, public :: ERR_FILE_NOT_FOUND = 2
    integer, parameter, public :: ERR_INVALID_INPUT = 3
    integer, parameter, public :: ERR_CONVERGENCE = 4
    integer, parameter, public :: ERR_SINGULAR_MATRIX = 5

    type :: error_info
        integer :: code = ERR_NONE
        character(len=256) :: message = ''
        character(len=64) :: source = ''
        integer :: line = 0
    contains
        procedure :: is_error => error_is_error
        procedure :: print => error_print
    end type

contains

    subroutine check_error(ierr, message, source)
        integer, intent(in) :: ierr
        character(len=*), intent(in) :: message, source
        if (ierr /= 0) then
            print '(A,A,A,A)', 'ERROR in ', trim(source), ': ', trim(message)
            error stop ierr
        end if
    end subroutine

    subroutine report_error(err)
        type(error_info), intent(in) :: err
        if (err%code /= ERR_NONE) then
            print '(A)', '================================'
            print '(A,I0)', 'Error code: ', err%code
            print '(A,A)', 'Source: ', trim(err%source)
            print '(A,A)', 'Message: ', trim(err%message)
            print '(A)', '================================'
        end if
    end subroutine

end module
```

### 2.2 Result Types (结果类型)
```fortran
module result_types_mod
    implicit none

    ! 类似Rust的Result类型
    type :: result_real
        logical :: success = .false.
        real(8) :: value = 0.0d0
        character(len=256) :: error_message = ''
    end type

    type :: result_array
        logical :: success = .false.
        real(8), allocatable :: value(:)
        character(len=256) :: error_message = ''
    end type

contains

    function ok_real(value) result(res)
        real(8), intent(in) :: value
        type(result_real) :: res
        res%success = .true.
        res%value = value
    end function

    function err_real(message) result(res)
        character(len=*), intent(in) :: message
        type(result_real) :: res
        res%success = .false.
        res%error_message = message
    end function

end module

! 使用示例
function safe_divide(a, b) result(res)
    real(8), intent(in) :: a, b
    type(result_real) :: res

    if (abs(b) < 1.0d-15) then
        res = err_real('Division by zero')
    else
        res = ok_real(a / b)
    end if
end function

! 调用
res = safe_divide(x, y)
if (res%success) then
    z = res%value
else
    print *, 'Error: ', trim(res%error_message)
end if
```

---

## 3. Memory Management (内存管理)

### 3.1 RAII Pattern (资源获取即初始化)
```fortran
module raii_example_mod
    implicit none

    type :: auto_array
        real(8), allocatable :: data(:)
    contains
        procedure :: init => array_init
        procedure :: resize => array_resize
        final :: array_cleanup  ! 自动清理
    end type

contains

    subroutine array_init(self, n)
        class(auto_array), intent(inout) :: self
        integer, intent(in) :: n
        if (allocated(self%data)) deallocate(self%data)
        allocate(self%data(n))
        self%data = 0.0d0
    end subroutine

    subroutine array_cleanup(self)
        type(auto_array), intent(inout) :: self
        if (allocated(self%data)) then
            deallocate(self%data)
            print *, 'Array deallocated'
        end if
    end subroutine

end module
```

### 3.2 Memory Pool (内存池)
```fortran
module memory_pool_mod
    implicit none
    private

    public :: memory_pool, get_work_array, release_work_array

    type :: memory_pool
        real(8), allocatable :: pool(:)
        integer :: pool_size = 0
        integer :: next_free = 1
        logical :: initialized = .false.
    contains
        procedure :: init => pool_init
        procedure :: allocate => pool_allocate
        procedure :: reset => pool_reset
    end type

    type(memory_pool), save :: global_pool

contains

    subroutine pool_init(self, size)
        class(memory_pool), intent(inout) :: self
        integer, intent(in) :: size
        allocate(self%pool(size))
        self%pool_size = size
        self%next_free = 1
        self%initialized = .true.
    end subroutine

    function pool_allocate(self, n) result(ptr)
        class(memory_pool), intent(inout) :: self
        integer, intent(in) :: n
        real(8), pointer :: ptr(:)

        if (self%next_free + n - 1 > self%pool_size) then
            error stop 'Memory pool exhausted'
        end if

        ptr => self%pool(self%next_free : self%next_free + n - 1)
        self%next_free = self%next_free + n
    end function

    subroutine pool_reset(self)
        class(memory_pool), intent(inout) :: self
        self%next_free = 1
    end subroutine

end module
```

---

## 4. Debugging and Profiling (调试与性能分析)

### 4.1 Debug Build Configuration
```fortran
! 编译时启用调试
! Intel: ifort -g -check all -warn all -traceback
! GNU:   gfortran -g -fcheck=all -Wall -fbacktrace

module debug_mod
    implicit none

#ifdef DEBUG
    logical, parameter :: DEBUG_MODE = .true.
#else
    logical, parameter :: DEBUG_MODE = .false.
#endif

contains

    subroutine debug_print(message, file, line)
        character(len=*), intent(in) :: message
        character(len=*), intent(in), optional :: file
        integer, intent(in), optional :: line

        if (DEBUG_MODE) then
            if (present(file) .and. present(line)) then
                print '(A,A,A,I0,A,A)', '[DEBUG] ', trim(file), ':', line, ' - ', trim(message)
            else
                print '(A,A)', '[DEBUG] ', trim(message)
            end if
        end if
    end subroutine

    subroutine assert(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (DEBUG_MODE .and. .not. condition) then
            print '(A,A)', 'Assertion failed: ', trim(message)
            error stop 1
        end if
    end subroutine

end module
```

### 4.2 Timing Utilities
```fortran
module timing_mod
    use iso_fortran_env, only: int64
    implicit none
    private

    public :: timer, start_timer, stop_timer, report_timers

    type :: timer
        character(len=64) :: name = ''
        integer(int64) :: start_count = 0
        integer(int64) :: total_count = 0
        integer :: call_count = 0
        logical :: running = .false.
    end type

    type(timer), save :: timers(100)
    integer, save :: n_timers = 0

contains

    subroutine start_timer(name)
        character(len=*), intent(in) :: name
        integer :: idx
        integer(int64) :: count

        idx = find_or_create_timer(name)
        call system_clock(count)
        timers(idx)%start_count = count
        timers(idx)%running = .true.
    end subroutine

    subroutine stop_timer(name)
        character(len=*), intent(in) :: name
        integer :: idx
        integer(int64) :: count

        idx = find_timer(name)
        if (idx > 0 .and. timers(idx)%running) then
            call system_clock(count)
            timers(idx)%total_count = timers(idx)%total_count + &
                (count - timers(idx)%start_count)
            timers(idx)%call_count = timers(idx)%call_count + 1
            timers(idx)%running = .false.
        end if
    end subroutine

    subroutine report_timers()
        integer :: i
        integer(int64) :: rate
        real(8) :: elapsed

        call system_clock(count_rate=rate)

        print '(A)', repeat('=', 60)
        print '(A)', 'Timer Report'
        print '(A)', repeat('=', 60)
        print '(A20,A15,A15,A10)', 'Name', 'Total (s)', 'Average (s)', 'Calls'
        print '(A)', repeat('-', 60)

        do i = 1, n_timers
            elapsed = real(timers(i)%total_count, 8) / real(rate, 8)
            print '(A20,F15.6,F15.6,I10)', trim(timers(i)%name), &
                elapsed, elapsed/max(1,timers(i)%call_count), timers(i)%call_count
        end do
        print '(A)', repeat('=', 60)
    end subroutine

end module

! 使用示例
call start_timer('stiffness_assembly')
call assemble_stiffness(...)
call stop_timer('stiffness_assembly')

call start_timer('solve')
call solve(...)
call stop_timer('solve')

call report_timers()
```

### 4.3 Logging System
```fortran
module logging_mod
    implicit none
    private

    public :: logger, log_debug, log_info, log_warning, log_error

    integer, parameter, public :: LOG_DEBUG = 1
    integer, parameter, public :: LOG_INFO = 2
    integer, parameter, public :: LOG_WARNING = 3
    integer, parameter, public :: LOG_ERROR = 4

    type :: logger_type
        integer :: level = LOG_INFO
        integer :: unit = 6  ! stdout
        logical :: to_file = .false.
        character(len=256) :: filename = ''
    contains
        procedure :: set_level => logger_set_level
        procedure :: set_file => logger_set_file
        procedure :: log => logger_log
    end type

    type(logger_type), save :: logger

contains

    subroutine logger_log(self, level, message, module_name)
        class(logger_type), intent(in) :: self
        integer, intent(in) :: level
        character(len=*), intent(in) :: message
        character(len=*), intent(in), optional :: module_name

        character(len=8) :: level_str
        character(len=20) :: timestamp

        if (level < self%level) return

        select case(level)
        case(LOG_DEBUG);   level_str = 'DEBUG'
        case(LOG_INFO);    level_str = 'INFO'
        case(LOG_WARNING); level_str = 'WARNING'
        case(LOG_ERROR);   level_str = 'ERROR'
        end select

        call date_and_time(time=timestamp)

        if (present(module_name)) then
            write(self%unit, '(A,A,A,A,A,A,A)') &
                '[', timestamp(1:6), '] [', trim(level_str), '] [', &
                trim(module_name), '] ', trim(message)
        else
            write(self%unit, '(A,A,A,A,A,A)') &
                '[', timestamp(1:6), '] [', trim(level_str), '] ', trim(message)
        end if
    end subroutine

end module

! 便捷函数
subroutine log_info(message)
    character(len=*), intent(in) :: message
    call logger%log(LOG_INFO, message)
end subroutine
```

---

## 5. Documentation Standards (文档标准)

### 5.1 Doxygen-Style Comments
```fortran
!> @brief Calculate the element stiffness matrix
!>
!> This subroutine computes the element stiffness matrix using
!> numerical integration (Gauss quadrature).
!>
!> @param[in]  elem    Element object containing geometry and DOF info
!> @param[in]  coords  Node coordinates (ndim x nnode)
!> @param[in]  D       Constitutive matrix (nstre x nstre)
!> @param[out] Ke      Element stiffness matrix (ndof x ndof)
!>
!> @note The stiffness matrix is symmetric. Only the upper triangle
!>       is computed and then copied to the lower triangle.
!>
!> @warning This routine assumes the element has been properly initialized.
!>
!> @see assemble_global_stiffness
!> @see calc_b_matrix
!>
subroutine calc_element_stiffness(elem, coords, D, Ke)
    class(element_base), intent(in) :: elem
    real(8), intent(in) :: coords(:,:)
    real(8), intent(in) :: D(:,:)
    real(8), intent(out) :: Ke(:,:)
    ! Implementation...
end subroutine
```

### 5.2 Module Header Template
```fortran
!===============================================================================
!> @file mohr_coulomb_mod.f90
!> @brief Mohr-Coulomb plasticity material model
!>
!> @details
!> This module implements the Mohr-Coulomb yield criterion for geomaterials.
!> The yield function is:
!>
!> \f$ F = \sigma_1 - \sigma_3 + (\sigma_1 + \sigma_3) \sin\phi - 2c\cos\phi \f$
!>
!> @author Your Name
!> @date 2024-01-15
!> @version 1.0
!>
!> @copyright (c) 2024 Your Organization
!>
!> @par References:
!> - Chen, W.F. (1975) Limit Analysis and Soil Plasticity
!> - de Souza Neto et al. (2008) Computational Methods for Plasticity
!===============================================================================
module mohr_coulomb_mod
    ! ...
end module
```

---

## 6. Quick Reference Card (快速参考卡)

| 旧代码模式 | 现代替代方案 |
|-----------|-------------|
| `integer i0, ic, f0` | `integer :: index, count, factor` |
| `real*8 x` | `real(real64) :: x` 或自定义 `WP` |
| `common /block/` | `module` 变量 |
| `include 'file.inc'` | `use module_name` |
| `goto 100` | `exit`, `cycle`, `return` |
| `character*20` | `character(len=20)` |
| `dimension(100)` | `allocatable` 或 `dimension(:)` |
| 手动内存管理 | `allocatable` + `final` |
| 字符串比较条件 | 多态 + 类型绑定过程 |
| 全局变量 | 参数传递 + 类型封装 |
