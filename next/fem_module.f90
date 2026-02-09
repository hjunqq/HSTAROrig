    include 'mkl_rci.f90'

    module fem_module
    !
    ! Module wrapper for HSTAR FEM subroutines.
    ! Phase 1 decomposition: pure file-level split, no logic changes.
    !
    ! Original: all 133 subroutines were in PROGRAM FEM90 CONTAINS block.
    ! Now: subroutines are in 14 include files, accessed via this module.
    !

    use variable_types
    use elements
    use global_var
    use materials
    use prescribed
    use applied_load
    use stiffness_matrix
    use internal_force
    use solver
    use output
    use temperature
    use meshfine
    !use portlib ! for what?
    use levelset !levelset
    use MKL_RCI  !20190810
    use MKL_RCI_TYPE !20190810
    use vsl_gauss_module

    implicit none

    ! === Shared variables (moved from PROGRAM scope) ===
    ! These were program-local variables accessed by CONTAINS subroutines
    ! via host association. Now they are module variables with the same
    ! implicit SAVE semantics.

    ! Reliability analysis (accessed by STATIC_U_reli, system_reliability, etc.)
    type beta_resultm
        real(irk),pointer:: ga(:)
        real(irk) beta
    end type beta_resultm
    type(beta_resultm),allocatable::betas(:)

    ! MKL DTRNLSP variables (accessed by parameter_back_analysis, trust_region)
    integer             RCI_REQUEST, res
    double precision    EPS(6), JAC_EPS
    double precision,allocatable:: FVEC(:), FJAC(:,:), Fvec1(:), Fvec2(:), &
        value_vc(:,:,:), fjac22(:,:,:)
    integer             ITER, ST_CR, INFO(6), SUCCESSFUL, ITER1, ITER2
    double precision    RS, R1, R2
    TYPE(HANDLE_TR) :: HANDLE

    ! Back-analysis verification variables
    integer             nback_point, nstoch, istoch, tbstep

    contains

    include 'fem_process.f90'
    include 'fem_static.f90'
    include 'fem_dynamic.f90'
    include 'fem_frequency.f90'
    include 'fem_rigid.f90'
    include 'fem_reliability.f90'
    include 'fem_backanalysis.f90'
    include 'fem_compliance.f90'
    include 'fem_force.f90'
    include 'fem_nr_update.f90'
    include 'fem_state.f90'
    include 'fem_init.f90'
    include 'fem_contact.f90'
    include 'fem_misc.f90'

    end module fem_module
