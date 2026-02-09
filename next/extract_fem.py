"""
Phase 1: Extract subroutines from Fem.f90 into 14 include files.

This script reads the original HSTAR/Fem.f90, extracts all subroutines from
the CONTAINS block, and writes them to categorized include files.
It also creates fem_module.f90 (module wrapper) and Fem.f90 (slimmed main program).

NO logic changes - pure file-level decomposition.
"""

import os

HSTAR_DIR = os.path.join(os.path.dirname(__file__), '..', 'HSTAR')
NEXT_DIR = os.path.dirname(__file__)
FEM_PATH = os.path.join(HSTAR_DIR, 'Fem.f90')

# Subroutine line ranges: (start_line, end_line) - 1-indexed, inclusive
# Grouped by target include file

FILE_ASSIGNMENTS = {
    'fem_init.f90': [
        (545, 585),    # modify_coord
        (588, 598),    # update_coord_blarge
        (601, 663),    # modify_element_information
        (11665, 12218), # modf_element_lib
        (12220, 12269), # normal_local
        (12271, 12508), # modf_var_prescribed
        (12510, 12560), # dfact_temp_pre
        (12563, 12571), # modf_var_prescribed_w
        (12573, 12709), # dfact_time_curve
        (14797, 14833), # read_permanent_strain
        (14836, 15252), # read_initial
        (15255, 15376), # placement_temperature
        (15378, 15413), # placement_temperature0
        (16292, 16346), # modf_time_order
    ],
    'fem_backanalysis.f90': [
        (668, 1024),    # parameter_back_analysis
        (1027, 1276),   # trust_region_back_analysis
        (1278, 1333),   # solve_dx
        (1335, 1360),   # update_bk
        (1362, 1463),   # rigid_dis_back_analysis
        (1465, 1559),   # nodal_value_back_analysis
        (1561, 1630),   # parameter_back_analysis_verify
        (2001, 2029),   # parameter_back_analysis_read
        (2031, 2064),   # trust_region_back_analysis_read
        (2066, 2198),   # parameter_back_analysis_verify_read
        (2201, 2329),   # observe_back_analysis_read
        (2368, 2955),   # back_analysis
        (2957, 3555),   # back_d_analysis
    ],
    'fem_process.f90': [
        (1633, 1999),   # process_analysis
        (2332, 2366),   # stab_initialize
    ],
    'fem_static.f90': [
        (3560, 4246),   # STATIC_U
        (7400, 7571),   # STATIC_U_P
        (7574, 8258),   # STATIC_U_PW
        (8260, 8467),   # STATIC_U_PWm
    ],
    'fem_compliance.f90': [
        (4249, 4289),   # modf_element_local_direction
        (4291, 4390),   # dudx
        (4392, 4488),   # cmatrix_c_formation
        (4491, 4584),   # Tcmatrix_c_formation
        (4587, 4703),   # cmatrix_dtv_formation
        (6427, 6721),   # forAdirect_back_analysis
        (6724, 6984),   # forAdirect
        (6986, 7113),   # matrix_rigid_dis
        (7115, 7178),   # matrix_nodal_value
    ],
    'fem_reliability.f90': [
        (4708, 5364),   # STATIC_U_reli
        (5367, 5386),   # betaindex
        (5388, 5427),   # RI3
        (5429, 5480),   # DANGLI
        (5483, 5651),   # system_reliability
        (5653, 5718),   # aft_gauss
        (5720, 5739),   # af_normalx
        (5742, 5762),   # af_normal
        (5764, 5789),   # af_beta
        (16707, 16801), # stab_rcandgy_reli
    ],
    'fem_contact.f90': [
        (5793, 5843),   # unit_force_trans
        (5846, 5913),   # unit_dis_force_trans
        (5917, 5953),   # result_node_to_center
        (12711, 12997), # contact_state
        (13000, 13065), # crack_state
        (13068, 13130), # local_stress
    ],
    'fem_rigid.f90': [
        (5957, 6214),   # STATIC_rigid_reli
        (6217, 6424),   # STATIC_rigid_1
    ],
    'fem_misc.f90': [
        (7182, 7206),   # write_stiff_u
        (7208, 7298),   # neuman_expan
        (7301, 7352),   # kdelt
        (7356, 7396),   # EIGV
        (16348, 16523), # RESTA_READ_WRITE
        (16526, 16705), # safety_factor
        (16802, 16915), # FORCE_interface
        (16919, 16946), # write_force_interface
        (16948, 16983), # judge_fine_mesh
        (16985, 17012), # find_remesh_element1
        (17014, 17041), # find_remesh_element2
        (17044, 17190), # find_remesh_stran
        (17193, 17231), # modf_inpwav
        (17234, 17415), # liquifaction_judge
        (17417, 17536), # permdeform_judge
        (17538, 17589), # gamamaxupdate
        (17591, 17621), # readgamamax
        (17623, 17655), # writegamamax
        (17660, 17738), # inivdval
        (17741, 17776), # strain_for_steel_bar
        (17780, 17853), # change_list
        (17855, 17873), # change4
        (17876, 18256), # GHM2ADINA
    ],
    'fem_dynamic.f90': [
        (8470, 9408),   # time_dependent
        (9410, 9524),   # value_submodel_boundary
        (9527, 9598),   # force_unit_rigid_accs
        (9601, 9701),   # dfat_rigid
        (9703, 9803),   # dfat_rigid_ctfor
        (9941, 10112),  # explicit
    ],
    'fem_frequency.f90': [
        (9806, 9938),   # frequency_analysis
        (10118, 10393), # response_spectrum
        (10396, 10476), # base_frequency_analysis
        (10479, 10524), # load_of_mass_base
        (10528, 10583), # load_of_mass
        (10585, 10613), # load_of_stiff
        (10615, 10707), # load_of_addtional_mass
    ],
    'fem_nr_update.f90': [
        (10711, 10863), # PREDICT
        (10865, 11010), # varupdate
        (11013, 11082), # relative_dis_watertight
        (11086, 11146), # construction_dis_modify
        (11152, 11170), # acc_modify
        (11172, 11189), # acc_check
        (11192, 11220), # acc_rigid
        (11223, 11242), # vel_modify
        (11245, 11275), # dis_modify
        (11278, 11299), # varupdate_w
        (11302, 11347), # updalfa
        (15417, 15563), # ALGORT
    ],
    'fem_state.f90': [
        (11350, 11431), # gpvarupdate
        (11434, 11486), # gpvarupdate1
        (11488, 11536), # gpvarupdate2
        (11538, 11577), # gpvar_initial
        (11580, 11619), # gpvar1_initial
        (11621, 11660), # gpvar2_initial
        (15566, 16007), # conver_load
        (16009, 16073), # conver_load_w
        (16075, 16206), # conver_nodal_value
    ],
    'fem_force.f90': [
        (13132, 13935), # FORCE_EXTERNAL
        (13940, 13987), # FORCE_EXTERNAL_w
        (13991, 14021), # assemble_boundt_estif
        (14023, 14057), # assemble_boundt_eload
        (14063, 14158), # eload_ifs2006
        (14160, 14216), # eload_ifs2006_w
        (14222, 14270), # eload_interface_fluid_solid
        (14272, 14305), # eload_interface_fs_w
        (14307, 14334), # eload_absorb_fluid
        (14336, 14362), # eload_absorb_fluid_w
        (14364, 14393), # eload_absorb_solid
        (14396, 14408), # eload_back_spring
        (14410, 14443), # eload_absorb_solid_w
        (14447, 14685), # FORCE_INTERNAL
        (14689, 14742), # force_release
        (14744, 14794), # ELOAD_INITIALIZE
    ],
}

# Include file order in the module (matters for readability, not compilation)
INCLUDE_ORDER = [
    'fem_process.f90',
    'fem_static.f90',
    'fem_dynamic.f90',
    'fem_frequency.f90',
    'fem_rigid.f90',
    'fem_reliability.f90',
    'fem_backanalysis.f90',
    'fem_compliance.f90',
    'fem_force.f90',
    'fem_nr_update.f90',
    'fem_state.f90',
    'fem_init.f90',
    'fem_contact.f90',
    'fem_misc.f90',
]


def read_fem():
    """Read the entire Fem.f90 into a list of lines (0-indexed)."""
    with open(FEM_PATH, 'r', encoding='utf-8', errors='replace') as f:
        return f.readlines()


def extract_include_files(lines):
    """Extract subroutines into include files."""
    total_extracted = 0
    for fname in INCLUDE_ORDER:
        ranges = FILE_ASSIGNMENTS[fname]
        outpath = os.path.join(NEXT_DIR, fname)
        with open(outpath, 'w', encoding='utf-8') as f:
            for i, (start, end) in enumerate(ranges):
                if i > 0:
                    f.write('\n')  # blank line between subroutines
                # Extract lines (1-indexed to 0-indexed)
                for lineno in range(start - 1, end):
                    f.write(lines[lineno])
                total_extracted += (end - start + 1)

        # Count subroutines
        n_subs = len(ranges)
        actual_lines = sum(e - s + 1 for s, e in ranges)
        print(f'  {fname}: {n_subs} subroutines, {actual_lines} lines')

    print(f'  Total extracted: {total_extracted} lines')
    return total_extracted


def create_fem_module(lines):
    """Create fem_module.f90 - the module wrapper."""
    outpath = os.path.join(NEXT_DIR, 'fem_module.f90')

    content = """    include 'mkl_rci.f90'

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
"""

    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f'  fem_module.f90: created')


def create_slim_fem(lines):
    """Create slimmed Fem.f90 - main program only, no CONTAINS."""
    outpath = os.path.join(NEXT_DIR, 'Fem.f90')

    with open(outpath, 'w', encoding='utf-8') as f:
        # Line 1 was "include 'mkl_rci.f90'" - now handled by fem_module
        # Lines 4-8: PROGRAM header
        f.write('\n')
        f.write('    PROGRAM FEM90\n')
        f.write('\n')
        f.write('    !                 M. PASTOR, TONCHUN LI AND P. MIRA\n')
        f.write('    !                          May, 1997\n')
        f.write('    !                      All rights reserved.\n')
        f.write('\n')
        f.write('    use fem_module\n')
        f.write('\n')
        f.write('    implicit none\n')
        f.write('\n')
        # Local variables only used in main body (lines 49-52 from original)
        f.write('    ! Local variables (only used in main program body)\n')
        f.write('    integer             ie, ig, nel_sub, npoin_sub, ngroup_sub\n')
        f.write('    integer,allocatable :: list_nel_sub(:), icpoin_sub(:), listpoin_sub(:), &\n')
        f.write('                           listgroup_sub(:), list_group_sub(:)\n')
        f.write('    integer,allocatable :: icgroup(:)\n')
        f.write('    real(irk),allocatable :: coord_sub(:,:)\n')
        f.write('\n')

        # Main body: lines 92-541 (before "contains" at line 543)
        # Note: line 542 is blank, line 543 is "    contains"
        for lineno in range(91, 541):  # 0-indexed: lines 92-541
            f.write(lines[lineno])

        f.write('\n')
        f.write('    END PROGRAM FEM90\n')

    print(f'  Fem.f90 (slim): created')


def verify_coverage(lines):
    """Verify all subroutine lines are accounted for."""
    # Collect all extracted line numbers
    extracted = set()
    for fname, ranges in FILE_ASSIGNMENTS.items():
        for start, end in ranges:
            for i in range(start, end + 1):
                extracted.add(i)

    # Check for subroutine lines in CONTAINS block not extracted
    # CONTAINS is at line 543, END PROGRAM at line 18258
    missing_subs = []
    in_sub = False
    sub_name = ''
    sub_start = 0

    for i in range(543, 18258):  # 1-indexed
        line = lines[i - 1].strip().lower()
        if line.startswith('subroutine ') or line.startswith('function '):
            if i not in extracted:
                sub_name = line.split('(')[0].split()[-1] if '(' in line else line.split()[-1]
                missing_subs.append((i, sub_name))

    if missing_subs:
        print(f'\n  WARNING: {len(missing_subs)} subroutines NOT extracted:')
        for lineno, name in missing_subs:
            print(f'    Line {lineno}: {name}')
    else:
        print(f'\n  Verification: All subroutine start lines covered.')

    # Count total subroutines
    total = sum(len(v) for v in FILE_ASSIGNMENTS.values())
    print(f'  Total subroutine blocks: {total}')

    return len(missing_subs) == 0


def main():
    print('Phase 1: Fem.f90 decomposition')
    print(f'Reading {FEM_PATH}...')
    lines = read_fem()
    print(f'  Total lines: {len(lines)}')

    print('\nCreating include files...')
    extract_include_files(lines)

    print('\nCreating fem_module.f90...')
    create_fem_module(lines)

    print('\nCreating slim Fem.f90...')
    create_slim_fem(lines)

    print('\nVerifying coverage...')
    verify_coverage(lines)

    print('\nDone!')


if __name__ == '__main__':
    main()
