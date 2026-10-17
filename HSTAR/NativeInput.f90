! TOML-03/04: two bounded profiles. No solver state mirror or legacy-card emitter.
module native_input
  use iso_fortran_env, only: real64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use yl_authoring_toml
  use yl_sets
  implicit none
  character(1024) :: native_mesh_prefix='1'
  logical :: native_mode = .false.
  logical :: native_dam = .false.
  ! GLB pilot (glb-pilot/DESIGN.md, G1-G12): bounded sizes instead of the Cook/dam fixed (2) arrays.
  integer, parameter :: NATIVE_MAXG=64, NATIVE_MAXB=64, NATIVE_MAXN=200000
  integer :: native_points=289, native_elements=256, native_materials=1, native_blocks=1
  integer :: native_elcounts(NATIVE_MAXG)=0, native_appearance(NATIVE_MAXG,NATIVE_MAXB)=1, native_reset(NATIVE_MAXB)=0
  integer :: native_bc(50,2)=0, native_nbc(2)=[17,17], native_edges(3,40)
  real(real64) :: native_e(NATIVE_MAXG), native_nu(NATIVE_MAXG), native_density(NATIVE_MAXG), native_gravity(NATIVE_MAXB)
  real(real64) :: native_pressure=0.
  integer :: native_substeps(NATIVE_MAXB), native_iterations(NATIVE_MAXB)
  ! HSTAR_NATIVE_GLB=contract (default, Cook/dam profile constants) | toml (global control from TOML;
  ! counts derived from 1.cor/1.ele and the TOML tables; .sol implied by solver.linear="profile")
  logical :: native_glb_toml=.false., native_glb_receipt_on=.false.
  integer :: native_ndim=2, native_gid_f=0
  integer :: native_dof_list(3)=[1,2,3]
  character(10) :: native_gname(NATIVE_MAXG)='Q4'
  character(20) :: native_kname(NATIVE_MAXG)='Default'
  integer :: native_gindex(NATIVE_MAXG)=5, native_gnnode(NATIVE_MAXG)=4
  ! beam pilot C2: B2 = legacy index 20 (b2, 2-node beam). point_direct from section.orientation.node_pair (0 0 = default).
  integer :: native_pdir(2,NATIVE_MAXG)=0
  logical, save :: native_gps(NATIVE_MAXG)=.false.   ! plane stress (legacy group SPTYPE 'PS'), section formulation "plane_stress"
  logical :: native_beam=.false.                  ! any B2 section: legacy mdofn 6 (3D) / 4 (2D, rotation in slot 4)
  ! steel25 pilot: [[section]] element "bond_slip" = legacy index 25 "steel" (2-node bond-slip spring), and the
  ! .glb global bond line ftcrack coefMpa ikindks doubsig ktan1 ktan2 nlocalbeam ndimnrt listglocbeam lelenrt
  ! from top-level [bond_slip] / [crack]. Defaults = the values native mode used before this pilot.
  logical :: native_bond=.false.
  integer :: native_ikindks=0, native_doubsig=2, native_nlocalbeam=0, native_ndimnrt=0
  real(real64) :: native_coefmpa=1.e6_real64, native_ktan1=1.e8_real64, native_ktan2=1.e8_real64
  real(real64) :: native_ftcrack=1.5e6_real64
  integer :: native_listglocbeam(NATIVE_MAXG)=0
  integer, allocatable :: native_lelenrt(:)
  integer, private :: force_unit, force_frame=0
  integer, parameter :: native_nodes(17) = [1,18,35,52,69,86,103,120,137,154,171,188,205,222,239,256,273]
  type(toml_doc_t), private :: doc
  logical, allocatable, private :: seen(:)
  ! PRE pilot: HSTAR_NATIVE_PRE=contract (default, profile constants) | toml (boundary data from TOML)
  logical :: native_pre_toml=.false., native_pre_receipt_on=.false.
  type, private :: pre_set_t
    integer :: ifixvar=0, itcurve=0, outfix=1
    real(real64) :: value=0._real64
    integer, allocatable :: nodes(:)
  end type
  type, private :: pre_step_t
    type(pre_set_t), allocatable :: sets(:)
  end type
  type, private :: pre_nset_t
    character(len=TOML_LEN_TEXT) :: name=""
    integer, allocatable :: nodes(:)
  end type
  type, private :: pre_cand_t
    integer :: dof, nset, bnd, itcurve, outfix
    real(real64) :: val
  end type
  type(pre_step_t), allocatable, private :: native_pre(:)
  ! LOA pilot: HSTAR_NATIVE_LOA=contract (default, profile constants) | toml (load data from TOML)
  logical :: native_loa_toml=.false., native_loa_receipt_on=.false., native_pre_curve1=.false.
  logical, private :: pre_used_curve1=.false.
  type :: loa_amp_t
    character(len=TOML_LEN_TEXT) :: name=""
    character(len=10) :: typ="LINEAR"
    real(real64), allocatable :: t(:), f(:)
    logical :: used=.false.
  end type
  type(loa_amp_t), allocatable :: native_amp(:)
  type, private :: loa_surf_t
    character(len=TOML_LEN_TEXT) :: name=""
    integer :: first=0, last=0, kind=2           ! kind 2 = edge2 (LOA edges), 4 = face4 (TEM faces)
    integer :: nnode=2, index=1                  ! legacy edge group header: Q4 edge (2, index 1) / B8 face (4, index 5)
    logical :: used=.false.
  end type
  type(loa_surf_t), allocatable, private :: native_surf(:)
  integer :: native_loa_nedge=0
  integer, allocatable :: native_loa_edges(:,:)     ! (5,n): nodes 1..nnode (in legacy order), row 5 = element
  ! B8 face node order of the 1.set "face" form for face loads: outward normal for the HSTAR B8 numbering
  ! (face 1 as train04 lists it, 4-3-2-1; faces 2..6 = yl_sets B8_FACE, outward on train04's 160 elements)
  integer, parameter :: PR_B8_FACE(4,6)=reshape([4,3,2,1, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8],[4,6])
  type :: loa_hyd_t
    integer :: first, last, itcurve
    real(real64) :: cor0, cor1, p0, p1, fact
    integer :: water=2                             ! legacy water: |water| = coordinate axis (hydrostatic); 0 = pressure
    real(real64) :: press=0                        ! water=0: value on the local normal at every face node
  end type
  type :: loa_step_t
    real(real64) :: g=0, dir(3)=0, fdir(3)=0
    integer :: gamp=0
    type(loa_hyd_t), allocatable :: hyd(:)
  end type
  type(loa_step_t), allocatable :: native_loa(:)
  type :: loa_pl_t
    integer :: itcurve=0
    real(real64) :: comp(6)=0
    integer :: nudofn=2
    integer, allocatable :: nodes(:)
  end type
  type(loa_pl_t), allocatable :: native_pl(:)
  integer, private :: loa_edge_group=0
  character(len=TOML_LEN_TEXT), allocatable, private :: loa_nset_refs(:)
  ! MAT pilot: HSTAR_NATIVE_MAT=contract (default, profile constants) | toml (material/section data from TOML)
  ! Contract defaults equal the constants Material.f90 used before the seam (ratio=1, thickness=1, alfa=1e-5).
  logical :: native_mat_toml=.false., native_mat_receipt_on=.false.
  real(real64) :: native_alfa(NATIVE_MAXG)=1.e-5_real64, native_thick(NATIVE_MAXG)=1._real64, native_ratio(NATIVE_MAXG)=1._real64
  integer :: native_group_mat(NATIVE_MAXG)=0   ! set to 1..NATIVE_MAXG at native_start (was [1,2])
  ! matline pilot (matline-pilot/DESIGN.md M1-M6, L1-L4): [[stress_strain_curve]] (legacy .mat curve table)
  ! and L2 line-element section data (legacy GEOMETRY ROD_or_BEAM: aera J Iy Iz), per material index.
  integer, parameter :: NATIVE_MAXC=1000
  type :: mat_curve_t
    character(len=TOML_LEN_TEXT) :: name=""
    real(real64), allocatable :: strain(:), stress(:)
  end type
  type(mat_curve_t), allocatable :: native_scurve(:)
  integer :: native_nscurve=0
  logical :: native_has_geom(NATIVE_MAXG)=.false.
  real(real64) :: native_area(NATIVE_MAXG)=0._real64, native_tj(NATIVE_MAXG)=0._real64, &
    native_iy(NATIVE_MAXG)=0._real64, native_iz(NATIVE_MAXG)=0._real64
  integer, private :: mat_receipt_n=0
  ! TEM minimal step (T1/M8): HSTAR_NATIVE_TEM=contract (default) | toml (no temperature loading allowed)
  logical :: native_tem_toml=.false., native_tem_receipt_on=.false.
  integer, private :: tem_receipt_n=0
  ! CREEP pilot phase 1 (creep-tem-pilot/DESIGN.md, C1-C8): HSTAR_NATIVE_MAN=contract (default,
  ! time_increment 1 / tolerances 1e-5) | toml (step[i].controls.time_increment, tolerance_force,
  ! tolerance_dof from TOML). Creep (Kelvin/Burgers = legacy icreep=4) lives under [[material]].creep
  ! with HSTAR_NATIVE_MAT=toml. Contract defaults equal the constants Fem.f90 used before the seam.
  logical :: native_man_toml=.false., native_man_receipt_on=.false.
  real(real64) :: native_ditime(NATIVE_MAXB)=1._real64, native_tol_force(NATIVE_MAXB)=1.e-5_real64, native_tol_dof(NATIVE_MAXB)=1.e-5_real64
  integer :: native_icreep(NATIVE_MAXG)=0
  real(real64) :: native_ek(NATIVE_MAXG)=0._real64, native_etak(NATIVE_MAXG)=0._real64, native_etam(NATIVE_MAXG)=0._real64
  real(real64), parameter :: native_kelvin_etam=1.e30_real64
  integer, private :: man_receipt_n=0
  ! SETS pilot (sets-pilot/DESIGN.md, S1-S10): [mesh] sets_file = "1.set" (needs HSTAR_NATIVE_GLB=toml).
  ! Without it the inline [[nset]]/[[surface]] edges/[[elset]] element_count forms are still read for
  ! one version, with a deprecation warning on stderr.
  logical :: native_sets_on=.false.
  ! TEMPERATURE pilot phase 2 (creep-tem-pilot/DESIGN.md section 9, T1-T8 confirmed). analysis.type =
  ! "transient_heat" under HSTAR_NATIVE_GLB=toml: B8 temperature field (dof 10, legacy type_problem S),
  ! time_dependent .man increments ([[step.increment]]), PARDISO (.sol implied), HEAT material, DABT curve,
  ! convection on face4 surfaces, cooling pipes (algo 3) and [[tie]] (HSTAR_NATIVE_NRT=toml, .nrt TRAL).
  logical :: native_heat=.false.
  logical :: native_nrt_toml=.false., native_nrt_receipt_on=.false.
  ! SOL seam (sol-pilot V1-V8): HSTAR_NATIVE_SOL=contract (default: implied .sol records, PROFILE 0 0 1 1 /
  ! heat PARDISO mtype 11 ncpu 2) | toml ([solver] linear = "profile" | "pardiso", [solver.pardiso]
  ! matrix_type -2|11, threads (default 1, V5), iparm = "default"). Receipt native-sol.txt.
  logical :: native_sol_toml=.false., native_sol_receipt_on=.false., native_pardiso=.false.
  integer :: native_sol_mtype=0, native_sol_threads=1
  integer, private :: sol_receipt_n=0
  integer, parameter :: NATIVE_MAXI=100
  integer :: native_ninc(NATIVE_MAXB)=0
  real(real64) :: native_inc_dt(NATIVE_MAXI,NATIVE_MAXB)=0._real64, native_inc_tol(NATIVE_MAXI,NATIVE_MAXB)=0._real64
  integer :: native_inc_nstep(NATIVE_MAXI,NATIVE_MAXB)=0, native_inc_miter(NATIVE_MAXI,NATIVE_MAXB)=0
  logical :: native_mat_heat(NATIVE_MAXG)=.false.
  ! joint pilot D3a: GOODMAN node-pair joints (section element="joint", formulation="goodman"), FCM sub-model only
  logical :: native_mat_goodman(NATIVE_MAXG)=.false., native_mat_janbu(NATIVE_MAXG)=.false., native_gjoint(NATIVE_MAXG)=.false.
  integer :: native_fcm_law(NATIVE_MAXG)=0
  ! joint pilot D3b: mc_joint thin layer (legacy CLASSICALEP / MCJOINT). native_mcj(1:5,i) = sigma0, hardening,
  ! frict_angle, dilan_angle, ft; native_mcj_curve(1:4,i) = csigma0, cfrict, cdilan, cft (scurves index, 0 = none)
  logical :: native_mat_mcj(NATIVE_MAXG)=.false.
  ! joint pilot D3c: goodman_janbu CONTACT joints (legacy MECHANICAL CONTACT group + GOODMAN / JANBU card).
  ! native_janbu(1:11,i) = K1 n Kzz Rf phi Kzx pa gamaw cohes Ft Kzy; native_gap(1:2,i) = gap0, ftcontact
  logical :: native_mat_contact(NATIVE_MAXG)=.false.
  ! output.stress_averaging_exclude: element sets left out of nodal stress averaging (legacy average_appear = 0)
  logical :: native_avg_off(NATIVE_MAXG)=.false.
  real(real64) :: native_janbu(11,NATIVE_MAXG)=0._real64, native_gap(2,NATIVE_MAXG)=0._real64
  integer :: native_igap0(NATIVE_MAXG)=0, native_icft(NATIVE_MAXG)=0
  real(real64) :: native_mcj(5,NATIVE_MAXG)=0._real64, native_elcod(NATIVE_MAXG)=0._real64, native_kstab=0._real64
  integer :: native_mcj_curve(4,NATIVE_MAXG)=0, native_mat_curve=0, native_type_nl=5
  integer :: native_gid_ms=0, native_gid_ep=0, native_gid_y=0, native_stab_matde=99999
  logical, private :: sr_amp_on=.false., any_mcj=.false.
  logical, allocatable, private :: scurve_ref(:)
  real(real64) :: native_fcm_gf(NATIVE_MAXG)=0._real64, native_fcm_ft(NATIVE_MAXG)=0._real64
  real(real64) :: native_diff(3,NATIVE_MAXG)=0._real64
  integer :: native_src_curve(NATIVE_MAXG)=0, native_place_curve(NATIVE_MAXG)=0
  character(len=TOML_LEN_TEXT), private :: mat_src_name(NATIVE_MAXG)='', mat_place_name(NATIVE_MAXG)=''
  integer :: native_nface=0
  integer, allocatable :: native_faces(:,:)          ! (5,nface): 4 face nodes + element
  type :: tem_conv_t
    integer :: first=0, last=0, itcurve=0, surf=0
    real(real64) :: beta=0._real64
  end type
  type :: tem_step_t
    type(tem_conv_t), allocatable :: conv(:)
  end type
  type(tem_step_t), allocatable :: native_tem_step(:)
  type :: tem_pipe_t
    integer :: igroup=0, twcurve=0
    real(real64) :: coef=1._real64, qw=0._real64, lamda=0._real64, dens=0._real64, cw=0._real64, t0=0._real64, t1=0._real64
    integer, allocatable :: seg(:,:)
  end type
  type(tem_pipe_t), allocatable :: native_pipe(:)
  integer :: native_pipe_algo=3
  type :: nrt_tie_t
    integer :: slave=0
    integer :: axis=0                             ! 0 = TRAL; 1..3 = ROTAX/ROTAY/ROTAZ (C2, rotation_tie set)
    integer, allocatable :: comps(:)              ! ROTA: displacement component of each term node
    integer, allocatable :: masters(:)
    real(real64), allocatable :: weights(:)
  end type
  type(nrt_tie_t), allocatable :: native_tie(:)
  logical :: native_nrt_dof_u=.false.             ! D1: [[tie]] dof = "U" -> legacy translg=88 (TRAL rows only)
  integer, private :: tem_face_blk=0, tem_conv_grp=0, tem_receipt2_n=0
  integer, allocatable, private :: mesh_ele(:,:)       ! (8,nelem) kept by validate_mesh_general for face checks
contains
  subroutine native_start()
    character(1024) :: flag, filename
    character(256) :: pre_env, loa_env, mat_env, tem_env, man_env, glb_env, nrt_env, sol_env
    integer :: n, i, j, file_size, ios
    n = command_argument_count()
    if(n == 0) return
    native_group_mat=[(i,i=1,NATIVE_MAXG)]
    call get_command_argument(1,flag)
    if(n /= 2 .or. trim(flag) /= '--input') call native_error('usage: hstar [--input case.toml]')
    call get_command_argument(2,filename,length=n)
    if(n > len(filename)) call native_error('input path too long')
    file_size=-1
    inquire(file=trim(filename),size=file_size,iostat=ios)
    if(ios /= 0 .or. file_size < 0) call native_error('cannot open TOML input')
    if(file_size > 65536) call native_error('TOML exceeds profile limit (64 KiB)')
    j=max(index(trim(filename),'/',back=.true.),index(trim(filename),achar(92),back=.true.))
    native_mesh_prefix=filename(:j)//'1'
    call get_environment_variable("HSTAR_NATIVE_PRE",pre_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_PRE too long")
    if(ios==0 .and. len_trim(pre_env)>0) then
      select case(trim(pre_env))
      case("contract"); native_pre_receipt_on=.true.
      case("toml"); native_pre_toml=.true.; native_pre_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_PRE must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_LOA",loa_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_LOA too long")
    if(ios==0 .and. len_trim(loa_env)>0) then
      select case(trim(loa_env))
      case("contract"); native_loa_receipt_on=.true.
      case("toml"); native_loa_toml=.true.; native_loa_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_LOA must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_MAT",mat_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_MAT too long")
    if(ios==0 .and. len_trim(mat_env)>0) then
      select case(trim(mat_env))
      case("contract"); native_mat_receipt_on=.true.
      case("toml"); native_mat_toml=.true.; native_mat_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_MAT must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_TEM",tem_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_TEM too long")
    if(ios==0 .and. len_trim(tem_env)>0) then
      select case(trim(tem_env))
      case("contract"); native_tem_receipt_on=.true.
      case("toml"); native_tem_toml=.true.; native_tem_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_TEM must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_MAN",man_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_MAN too long")
    if(ios==0 .and. len_trim(man_env)>0) then
      select case(trim(man_env))
      case("contract"); native_man_receipt_on=.true.
      case("toml"); native_man_toml=.true.; native_man_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_MAN must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_NRT",nrt_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_NRT too long")
    if(ios==0 .and. len_trim(nrt_env)>0) then
      select case(trim(nrt_env))
      case("contract"); native_nrt_receipt_on=.true.
      case("toml"); native_nrt_toml=.true.; native_nrt_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_NRT must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_GLB",glb_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_GLB too long")
    if(ios==0 .and. len_trim(glb_env)>0) then
      select case(trim(glb_env))
      case("contract"); native_glb_receipt_on=.true.
      case("toml"); native_glb_toml=.true.; native_glb_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_GLB must be contract or toml")
      end select
    endif
    call get_environment_variable("HSTAR_NATIVE_SOL",sol_env,status=ios)
    if(ios==-1) call native_error("HSTAR_NATIVE_SOL too long")
    if(ios==0 .and. len_trim(sol_env)>0) then
      select case(trim(sol_env))
      case("contract"); native_sol_receipt_on=.true.
      case("toml"); native_sol_toml=.true.; native_sol_receipt_on=.true.
      case default; call native_error("HSTAR_NATIVE_SOL must be contract or toml")
      end select
    endif
    call toml_read(trim(filename),doc)
    if(doc%failed) then
      write(error_unit,*) 'TOML line',doc%fail_line
      call native_error(trim(doc%message))
    endif
    allocate(seen(doc%n)); seen=.false.
    do i=1,doc%n
      do j=1,i-1
        if(doc%entry(i)%path == doc%entry(j)%path) call native_error('duplicate: '//trim(doc%entry(i)%path))
      enddo
    enddo
    if(native_glb_toml) then
      ! G2: case.name is only a name; the case skeleton comes from the TOML tables and the mesh.
      if(.not.(native_pre_toml .and. native_loa_toml .and. native_mat_toml .and. native_man_toml)) &
        call native_error("GLB: HSTAR_NATIVE_GLB=toml requires HSTAR_NATIVE_PRE, HSTAR_NATIVE_LOA, "// &
          "HSTAR_NATIVE_MAT and HSTAR_NATIVE_MAN = toml (their contract constants belong to the Cook/dam profiles)")
      if(native_tem_receipt_on .and. .not.native_tem_toml) &
        call native_error("GLB: HSTAR_NATIVE_GLB=toml requires HSTAR_NATIVE_TEM = toml or unset")
      call parse_glb()
      if(native_heat .and. .not.native_tem_toml) call native_error("GLB: analysis.type ""transient_heat"" "// &
        "requires HSTAR_NATIVE_TEM = toml (convection, pipes and ECWPIPE come from the TOML)")
    else
    if(native_nrt_toml) call native_error("NRT: HSTAR_NATIVE_NRT=toml requires HSTAR_NATIVE_GLB=toml")
    if(native_sol_toml) call native_error("SOL: HSTAR_NATIVE_SOL=toml requires HSTAR_NATIVE_GLB=toml")
    i=key_index('case.name')
    if(doc%entry(i)%kind/=TV_STR) call native_error('expected string: case.name')
    select case(trim(doc%entry(i)%svalue))
    case('cooks_membrane')
      include 'NativeContract.inc'
      native_bc(1:17,1)=native_nodes
      native_bc(1:17,2)=native_nodes
      native_elcounts(1)=256
    case('train01_gravdam_static')
      native_dam=.true.
      include 'NativeDamContract.inc'
      native_points=1705; native_elements=1600; native_materials=2; native_blocks=2
      native_elcounts(1:2)=[960,640]; native_appearance(1:2,1)=[1,0]; native_appearance(1:2,2)=[1,1]
      native_reset(1:2)=[0,1]; native_nbc=[50,41]
      native_kname(1)='Foundation'; native_kname(2)='Dam'; native_gid_f=1
      do i=1,25
        native_bc(2*i-1,1)=1+41*(i-1)
        native_bc(2*i,1)=41*i
      enddo
      do i=1,41
        native_bc(i,2)=i
      enddo
      do i=1,40
        native_edges(:,i)=[1009+17*(i-1),1026+17*(i-1),961+16*(i-1)]
      enddo
      native_edges(1,1)=997
      if(.not.native_loa_toml) native_pressure=number('step[2].load[2].distribution.scale')
      if(native_pressure<0 .or. native_pressure>1.e6_real64) call native_error('pressure scale out of range [0,1e6]')
    case default
      call native_error('unsupported case.name')
    end select
    endif
    if(native_mat_toml) call parse_mat()
    do i=1,native_materials
      if(native_mat_toml) exit
      native_e(i)=number(indexed('material',i)//'.E')
      native_nu(i)=number(indexed('material',i)//'.nu')
      native_density(i)=number(indexed('material',i)//'.density')
      if(native_e(i)<=0 .or. native_e(i)>1.e15_real64) call native_error('E out of range (0,1e15] Pa')
      if(native_nu(i)<=-1 .or. native_nu(i)>=0.5_real64) call native_error('nu out of range (-1,0.5)')
      if(native_density(i)<0 .or. native_density(i)>1.e6_real64) call native_error('density out of range [0,1e6]')   ! D3: 0 allowed (massless/weightless fill)
    enddo
    do i=1,native_blocks
      if(native_heat) exit                 ! transient_heat: increments come from [[step.increment]] (parse_man)
      if(.not.native_loa_toml) native_gravity(i)=number(indexed('step',i)//'.load[1].magnitude')
      native_substeps(i)=integer_value(indexed('step',i)//'.controls.substeps')
      native_iterations(i)=integer_value(indexed('step',i)//'.controls.max_iterations')
      if(.not.native_loa_toml .and. (native_gravity(i)<=0 .or. native_gravity(i)>1.e4_real64)) call native_error('gravity out of range (0,1e4]')
      if(native_substeps(i)<1 .or. native_substeps(i)>100) call native_error('substeps out of range [1,100]')
      if(native_iterations(i)<1 .or. native_iterations(i)>100) call native_error('max_iterations out of range [1,100]')
    enddo
    if(native_man_toml) call parse_man()
    if(native_loa_toml) call parse_loa()
    if(native_pre_toml) call parse_pre()
    if(native_heat) call heat_resolve()
    if(native_nrt_toml) call parse_nrt()
    if(sr_amp_on) native_mat_curve=amp_ref('analysis.strength_reduction.amplitude')   ! D3b: legacy MAT_DE curve
    if(native_stab_matde/=99999 .and. native_stab_matde>native_blocks) call native_error('GLB: analysis.strength_reduction.'// &
      'reset_state_each_increment_from_step = '//itoa(native_stab_matde)//' exceeds the number of steps')
    if((native_kstab/=0 .or. sr_amp_on) .and. .not.any_mcj) &
      call native_error('GLB: analysis.strength_reduction needs an mc_joint material (joint pilot D3b)')
    if(native_gid_ms+native_gid_ep+native_gid_y>0 .and. .not.(any_mcj .or. any(native_mat_contact))) &
      call native_error('GLB: output.field ms / ep / yield need an mc_joint or goodman_janbu material '// &
        '(joint pilot D3b / D3c)')
    if(native_type_nl==4 .and. .not.(any_mcj .or. any(native_mat_contact))) call native_error('GLB: stiffness_update '// &
      '"every_iteration" needs an mc_joint or goodman_janbu material (joint pilot D3b / D3c)')
    if(native_loa_toml) call loa_check_used()
    if(native_tem_toml .and. .not.native_heat) call parse_tem()
    native_pre_curve1 = (.not.native_pre_toml .and. native_dam) .or. pre_used_curve1
    if(native_sets_on) then
      if(len(sets_first_unreferenced())>0) call native_error("SETS: unreferenced "//sets_first_unreferenced())
      call sets_receipt("native-sets.txt")
    endif
    do i=1,doc%n
      if(.not.seen(i)) call native_error('unsupported key: '//trim(doc%entry(i)%path))
    enddo
    if(native_glb_toml) then
      call validate_mesh_general()
    else
      call validate_mesh()
    endif
    ! sol-pilot V9: with HSTAR_NATIVE_GLB=toml the inline set tables are no longer accepted (code removal later)
    if(native_glb_toml .and. .not.native_sets_on .and. &
       (doc%count_of("nset")>0 .or. doc%count_of("surface")>0 .or. doc%count_of("elset")>0)) &
      call native_error("SETS: inline [[nset]]/[[surface]]/[[elset]] tables are no longer accepted with "// &
        "HSTAR_NATIVE_GLB=toml; move them to 1.set ([mesh] sets_file; tools/fortran-interface/sets_tool.py convert)")
    native_mode=.true.
    print *, 'Native TOML input accepted (mesh relative to TOML file)'
    print *, 'TOML material E,nu,density:',native_e(1:native_materials),native_nu(1:native_materials),native_density(1:native_materials)
    if(.not.native_heat) print *, 'TOML gravity,substeps,max_iterations:',native_gravity(1:native_blocks), &
      native_substeps(1:native_blocks),native_iterations(1:native_blocks)
    if(native_heat) print *, 'TOML transient_heat: blocks, increments:',native_blocks,native_ninc(1:native_blocks)
  end subroutine

  function indexed(name,i) result(path)
    character(*), intent(in) :: name
    integer, intent(in) :: i
    character(:), allocatable :: path
    character(12) :: suffix
    write(suffix,'(i0)') i
    path=name//'['//trim(suffix)//']'
  end function

  subroutine native_material_receipt(i,e,nu,rho,alfa)
    integer, intent(in) :: i
    real(real64), intent(in) :: e,nu,rho,alfa
    integer :: u
    if(.not.native_dam) return
    open(newunit=u,file='native-material-'//achar(48+i)//'.txt',status='replace')
    write(u,'(a,es25.16)') 'E=',e
    write(u,'(a,es25.16)') 'nu=',nu
    write(u,'(a,es25.16)') 'density=',rho
    write(u,'(a,es25.16)') 'thermal_expansion=',alfa
    close(u)
  end subroutine

  subroutine native_transition(i,phase,norm)
    integer, intent(in) :: i
    character(*), intent(in) :: phase
    real(real64), intent(in) :: norm
    integer :: u
    if(.not.native_dam) return
    if(i==1 .and. phase=='before') then
      open(newunit=u,file='native-transitions.txt',status='replace')
    else
      open(newunit=u,file='native-transitions.txt',status='old',position='append')
    endif
    write(u,'(i0,1x,a,1x,es25.16)') i,phase,norm
    close(u)
  end subroutine

  ! Extra precision only for pressure-difference validation; GiD output is unchanged.
  subroutine native_force_begin()
    if(.not.native_dam) return
    force_frame=force_frame+1
    if(force_frame==1) then
      open(newunit=force_unit,file='native-force.txt',status='replace')
    else
      open(newunit=force_unit,file='native-force.txt',status='old',position='append')
    endif
  end subroutine

  subroutine native_force_row(node,values)
    integer, intent(in) :: node
    real(real64), intent(in) :: values(:)
    if(.not.native_dam) return
    write(force_unit,'(2(i0,1x),2(es25.16,1x))') force_frame,node,values
    if(node==native_points) close(force_unit)
  end subroutine

  subroutine native_stage_receipt(i,active,materials,reset,g,steps,iters,edges,curves)
    integer, intent(in) :: i,active(:),materials(:),reset,steps,iters,edges,curves(:)
    real(real64), intent(in) :: g
    integer :: u,j
    if(.not.native_dam) return
    open(newunit=u,file='native-stage-'//achar(48+i)//'.txt',status='replace')
    write(u,'(a,i0)') 'reset=',reset
    do j=1,size(active)
      write(u,'(a,i0,a,i0)') 'active_',j,'=',active(j)
      write(u,'(a,i0,a,i0)') 'material_',j,'=',materials(j)
      write(u,'(a,i0,a,i0)') 'gravity_curve_',j,'=',curves(j)
    enddo
    write(u,'(a,es25.16)') 'gravity=',g
    write(u,'(a,i0)') 'substeps=',steps
    write(u,'(a,i0)') 'max_iterations=',iters
    write(u,'(a,i0)') 'pressure_edges=',edges
    close(u)
  end subroutine

  subroutine native_pressure_receipt(c0,c1,p0,p1,scale)
    real(real64), intent(in) :: c0,c1,p0,p1,scale
    integer :: u
    open(newunit=u,file='native-pressure.txt',status='replace')
    write(u,'(a,es25.16)') 'at_1=',c0
    write(u,'(a,es25.16)') 'at_2=',c1
    write(u,'(a,es25.16)') 'value_1=',p0
    write(u,'(a,es25.16)') 'value_2=',p1
    write(u,'(a,es25.16)') 'scale=',scale
    close(u)
  end subroutine

  ! Called with actual solver state after material/load/step readers have consumed it.
  subroutine native_receipt(e,nu,rho,alfa,g,steps,iters)
    real(real64), intent(in) :: e,nu,rho,alfa,g
    integer, intent(in) :: steps,iters
    integer :: u
    open(newunit=u,file='native-consumed.txt',status='replace',action='write')
    write(u,'(a,es25.16)') 'E=',e
    write(u,'(a,es25.16)') 'nu=',nu
    write(u,'(a,es25.16)') 'density=',rho
    write(u,'(a,es25.16)') 'thermal_expansion=',alfa
    write(u,'(a,es25.16)') 'gravity=',g
    write(u,'(a,i0)') 'substeps=',steps
    write(u,'(a,i0)') 'max_iterations=',iters
    close(u)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! PRE pilot (HSTAR_NATIVE_PRE=toml): boundary data come from validated TOML
  ! instead of the profile constants. See pilot_design.md (decisions 1-10).
  ! ---------------------------------------------------------------------------
  function itoa(i) result(s)
    integer, intent(in) :: i
    character(:), allocatable :: s
    character(16) :: b
    write(b,"(i0)") i
    s=trim(b)
  end function

  logical function pre_owned(path)
    character(*), intent(in) :: path
    pre_owned=.false.
    if(.not.native_pre_toml) return
    pre_owned = index(path,"nset[")==1 .or. index(path,".boundary[")>0
  end function

  function string_value(path) result(s)
    character(*), intent(in) :: path
    character(:), allocatable :: s
    integer :: k
    k=key_index(path)
    if(doc%entry(k)%kind /= TV_STR) call native_error("expected string: "//path)
    s=trim(doc%entry(k)%svalue)
  end function

  subroutine parse_pre()
    type(pre_nset_t), allocatable :: ns(:)
    type(pre_cand_t), allocatable :: cand(:)
    type(pre_cand_t) :: t
    integer, allocatable :: owner(:,:)
    character(*), parameter :: dofname(10)=["ux","uy","uz","r4","r5","r6","p7","w8","a9","T "]
    character(:), allocatable :: p, q, sp, nm, dn, amp
    integer :: nn, i, j, k, n, v, s, nb, b, nd, d, idx, dof, kk, itc, outf, ncand, ndim_pre
    logical, allocatable :: nref(:)
    real(real64) :: val
    ndim_pre=native_ndim
    if(native_sets_on) then
      nn=0
      do i=1,sets_n
        if(sets(i)%kind==SK_NSET) nn=nn+1
      enddo
      ! D2 (conv step A): a 1.set without any nset is accepted (constraint-free deck, e.g. contact3d); a boundary
      ! entry that names an nset still fails the name lookup below.
      allocate(ns(nn)); j=0
      do i=1,sets_n
        if(sets(i)%kind/=SK_NSET) cycle
        j=j+1; ns(j)%name=sets(i)%name; ns(j)%nodes=sets(i)%m
        if(any(ns(j)%nodes>native_points)) call native_error("node out of range 1.."//itoa(native_points)//" in nset "//trim(ns(j)%name))
      enddo
    else
    nn=doc%count_of("nset")
    if(nn<1 .and. .not.native_heat) call native_error("at least one [[nset]] is required")
    allocate(ns(max(nn,1)))
    do i=1,nn
      p=indexed("nset",i)
      k=key_index(p//".@table")
      ns(i)%name=string_value(p//".name")
      if(len(ns(i)%name)==0) call native_error("empty nset name: "//p)
      do j=1,i-1
        if(trim(ns(j)%name)==trim(ns(i)%name)) call native_error("duplicate nset name: "//trim(ns(i)%name))
      enddo
      n=integer_value(p//".nodes.count")
      if(n<1) call native_error("nset has no nodes: "//trim(ns(i)%name))
      allocate(ns(i)%nodes(n))
      do k=1,n
        v=integer_value(indexed(p//".nodes",k))
        if(v<1 .or. v>native_points) call native_error("node out of range 1.."//itoa(native_points)//" in nset "//trim(ns(i)%name)//": "//itoa(v))
        do j=1,k-1
          if(ns(i)%nodes(j)==v) call native_error("duplicate node "//itoa(v)//" in nset "//trim(ns(i)%name))
        enddo
        ns(i)%nodes(k)=v
      enddo
    enddo
    endif
    allocate(nref(nn)); nref=.false.
    amp=""
    if(.not.native_loa_toml) then
      k=doc%find("amplitude[1].name")
      if(k==0) call native_error("missing: amplitude[1].name")
      amp=trim(doc%entry(k)%svalue)
    endif
    allocate(native_pre(native_blocks))
    do s=1,native_blocks
      sp=indexed("step",s)
      nb=doc%count_of(sp//".boundary")
      ! D2: a constraint-free deck (no nset at all) may have steps without boundary entries
      if(nb<1 .and. .not.native_heat .and. nn>0) call native_error("step "//itoa(s)//": at least one [[step.boundary]] is required")
      allocate(owner(10,native_points)); owner=0
      allocate(cand(max(3*nb,1))); ncand=0
      do b=1,nb
        q=indexed(sp//".boundary",b)
        k=key_index(q//".@table")
        nm=string_value(q//".nset")
        idx=0
        do i=1,nn
          if(trim(ns(i)%name)==nm) idx=i
        enddo
        if(idx==0) call native_error("step "//itoa(s)//" boundary "//itoa(b)//": undefined nset "//nm)
        nref(idx)=.true.
        if(native_sets_on) kk=sets_find(SK_NSET,nm)
        val=number(q//".value")
        itc=0
        kk=doc%find(q//".amplitude")
        if(kk>0) then
          if(native_loa_toml) then
            itc=amp_ref(q//".amplitude")
          else
            if(string_value(q//".amplitude")/=amp) call native_error("step "//itoa(s)//" boundary "//itoa(b)//": undefined amplitude")
            itc=1
            pre_used_curve1=.true.
          endif
        endif
        outf=1
        kk=doc%find(q//".output_reactions")
        if(kk>0) then
          kk=key_index(q//".output_reactions")
          if(doc%entry(kk)%kind /= TV_BOOL) call native_error("expected boolean: "//q//".output_reactions")
          outf=merge(1,0,doc%entry(kk)%lvalue)
        endif
        nd=integer_value(q//".dof.count")
        if(nd<1) call native_error("step "//itoa(s)//" boundary "//itoa(b)//": empty dof list")
        do d=1,nd
          dn=string_value(indexed(q//".dof",d))
          select case(dn)
          case("ux","uy","uz")
            if(native_heat) call native_error("dof "//dn//" not available in a transient_heat analysis (only T) "// &
              "(step "//itoa(s)//" boundary "//itoa(b)//")")
            dof=index("xyz",dn(2:2))
          case("T")
            if(.not.native_heat) call native_error("dof T needs analysis.type ""transient_heat"" (step "//itoa(s)// &
              " boundary "//itoa(b)//")")
            dof=10
          case("rx","ry","rz")          ! beam pilot C2: legacy slots 4/5/6 (3D), 2D rotation (rz) in slot 4
            if(.not.native_beam) call native_error("rotational dof "//dn//" needs a B2 (beam) section (step "//itoa(s)// &
              " boundary "//itoa(b)//")")
            if(native_ndim==2) then
              if(dn/="rz") call native_error("dof "//dn//" not available in 2D (only rz) (step "//itoa(s)//" boundary "//itoa(b)//")")
              dof=4
            else
              dof=3+index("xyz",dn(2:2))
            endif
          case default
            call native_error("unknown dof "//dn//" (expected ux, uy, uz, rx, ry, rz or T) in step "//itoa(s)//" boundary "//itoa(b))
          end select
          if(dof<=3 .and. dof>ndim_pre) call native_error("dof "//dn//" not available: mesh.dimension is "//itoa(ndim_pre))
          do k=1,size(ns(idx)%nodes)
            v=ns(idx)%nodes(k)
            if(owner(dof,v)/=0) call native_error("duplicate constraint: step "//itoa(s)//" node "//itoa(v)//" dof "//trim(dofname(dof))// &
                 " (boundary "//itoa(owner(dof,v))//" and boundary "//itoa(b)//")")
            owner(dof,v)=b
          enddo
          ncand=ncand+1
          cand(ncand)=pre_cand_t(dof,idx,b,itc,outf,val)
        enddo
      enddo
      ! canonical order (decision 7): dof ascending, then nset definition order
      do i=2,ncand
        t=cand(i); j=i-1
        do while(j>=1)
          if(cand(j)%dof<t%dof .or. (cand(j)%dof==t%dof .and. cand(j)%nset<t%nset)) exit
          cand(j+1)=cand(j); j=j-1
        enddo
        cand(j+1)=t
      enddo
      allocate(native_pre(s)%sets(ncand))
      do i=1,ncand
        native_pre(s)%sets(i)%ifixvar=cand(i)%dof
        native_pre(s)%sets(i)%itcurve=cand(i)%itcurve
        native_pre(s)%sets(i)%outfix=cand(i)%outfix
        native_pre(s)%sets(i)%value=cand(i)%val
        native_pre(s)%sets(i)%nodes=ns(cand(i)%nset)%nodes
      enddo
      deallocate(owner,cand)
    enddo
    do i=1,nn
      if(native_sets_on) exit             ! 1.set: one global unreferenced check at the end of native_start
      if(.not.nref(i) .and. .not.loa_nset_referenced(ns(i)%name)) call native_error("unreferenced nset: "//trim(ns(i)%name))
    enddo
  end subroutine

  integer function native_pre_nfix(iblk) result(n)
    integer, intent(in) :: iblk
    if(.not.native_pre_toml) then
      n=2
    else
      if(iblk<1 .or. iblk>size(native_pre)) call native_error("PRE: no boundary data for block "//itoa(iblk))
      n=size(native_pre(iblk)%sets)
    endif
  end function

  subroutine native_pre_head(iblk,iset,ifixvar,nfixnods,itcurve,tfixvar,outfix,jfixvar,gamawx,nextr)
    integer, intent(in) :: iblk,iset
    integer, intent(out) :: ifixvar,nfixnods,itcurve,tfixvar,outfix,jfixvar,nextr
    real(real64), intent(out) :: gamawx
    ifixvar=native_pre(iblk)%sets(iset)%ifixvar
    nfixnods=size(native_pre(iblk)%sets(iset)%nodes)
    itcurve=native_pre(iblk)%sets(iset)%itcurve
    outfix=native_pre(iblk)%sets(iset)%outfix
    tfixvar=0; jfixvar=0; gamawx=0._real64; nextr=0
  end subroutine

  subroutine native_pre_nodes(iblk,iset,list)
    integer, intent(in) :: iblk,iset
    integer, intent(out) :: list(:)
    list=native_pre(iblk)%sets(iset)%nodes
  end subroutine

  subroutine native_pre_values(iblk,iset,vals)
    integer, intent(in) :: iblk,iset
    real(real64), intent(out) :: vals(:)
    vals=native_pre(iblk)%sets(iset)%value
  end subroutine

  ! Receipt of what prescrib_set actually stored (read back from prescrib, not from TOML).
  subroutine native_pre_receipt(iblk,n,ifixset,ifixvar,node,ldofix,itcurve,outfix,jfixvar,vdofix)
    integer, intent(in) :: iblk,n,ifixset(n),ifixvar(n),node(n),ldofix(n),itcurve(n),outfix(n),jfixvar(n)
    real(real64), intent(in) :: vdofix(n)
    integer :: u,i
    if(iblk==1) then
      open(newunit=u,file="native-pre.txt",status="replace",action="write")
      write(u,"(a)") "native-pre v1 mode="//trim(merge("toml    ","contract",native_pre_toml))
    else
      open(newunit=u,file="native-pre.txt",status="unknown",position="append",action="write")
    endif
    write(u,"(a,i0,a,i0)") "block=",iblk," ndofix=",n
    write(u,"(a)") "ifixset ifixvar node ldofix itcurve outfix jfixvar value"
    do i=1,n
      write(u,"(7(i0,1x),es25.16)") ifixset(i),ifixvar(i),node(i),ldofix(i),itcurve(i),outfix(i),jfixvar(i),vdofix(i)
    enddo
    close(u)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! LOA pilot (HSTAR_NATIVE_LOA=toml): curves, point loads, edges, hydrostatic and
  ! gravity come from validated TOML instead of the profile constants.
  ! Design: docs/fortran-interface/loa-pilot/DESIGN.md (decisions D1-D9, confirmed).
  ! ---------------------------------------------------------------------------
  logical function loa_owned(path)
    character(*), intent(in) :: path
    loa_owned=.false.
    if(.not.native_loa_toml) return
    loa_owned = index(path,"amplitude[")==1 .or. index(path,"surface[")==1 .or. &
                index(path,"load[")==1 .or. index(path,".load[")>0
  end function

  integer function amp_by_name(nm) result(k)
    character(*), intent(in) :: nm
    integer :: i
    k=0
    do i=1,size(native_amp)
      if(trim(native_amp(i)%name)==trim(nm)) k=i
    enddo
  end function

  integer function amp_ref(path) result(k)
    character(*), intent(in) :: path
    character(:), allocatable :: nm
    nm=string_value(path)
    k=amp_by_name(nm)
    if(k==0) call native_error("undefined amplitude "//nm//" ("//path//")")
    native_amp(k)%used=.true.
  end function

  logical function present_key(path)
    character(*), intent(in) :: path
    present_key = doc%find(path)>0
  end function

  subroutine parse_loa()
    character(:), allocatable :: p, q, sp, nm, typ
    integer :: na, i, j, k, n, s, nl, l, ns, ie, total, ng, nh, ip, ii, nnset, idx
    integer :: e(3), wa
    real(real64) :: wl, bt, gw, g, d(3), fd(3)
    integer :: nd
    type(loa_hyd_t), allocatable :: tmp(:)
    logical, allocatable :: srf_step(:)
    ! ---- [[amplitude]]
    na=doc%count_of("amplitude")
    if(na<1) call native_error("LOA: at least one [[amplitude]] is required")
    allocate(native_amp(na))
    do i=1,na
      p=indexed("amplitude",i)
      k=key_index(p//".@table")
      nm=string_value(p//".name")
      if(len(nm)==0) call native_error("LOA: empty amplitude name: "//p)
      if(amp_by_name(nm)/=0) call native_error("LOA: duplicate amplitude name: "//nm)
      native_amp(i)%name=nm
      typ=string_value(p//".type")
      if(typ=="dabt") then
        ! T1: legacy DABT d(f)=a*b/(b+t)**2, theta(t)=a*t/(b+t); a=ttime_curve(1), b=dfact_curve(1)
        if(.not.native_heat) call native_error("LOA: amplitude "//nm//": type dabt needs analysis.type 'transient_heat'")
        native_amp(i)%typ="DABT"
        allocate(native_amp(i)%t(1),native_amp(i)%f(1))
        native_amp(i)%t(1)=number(p//".final_rise")
        native_amp(i)%f(1)=number(p//".half_time")
        if(native_amp(i)%t(1)<=0 .or. native_amp(i)%t(1)>1.e3_real64) &
          call native_error("LOA: amplitude "//nm//": final_rise out of range (0,1e3]")
        if(native_amp(i)%f(1)<=0 .or. native_amp(i)%f(1)>1.e6_real64) &
          call native_error("LOA: amplitude "//nm//": half_time out of range (0,1e6]")
        cycle
      endif
      if(typ/="linear") call native_error("LOA: amplitude "//nm//": type "//typ//" is not supported (only linear"// &
        trim(merge(", dabt",repeat(" ",6),native_heat))//")")
      if(.not.present_key(p//".before_first")) call native_error("LOA: amplitude "//nm//": before_first is required (only ""zero"")")
      if(string_value(p//".before_first")/="zero") call native_error("LOA: amplitude "//nm//": before_first must be ""zero""")
      n=integer_value(p//".points.count")
      if(n<2 .or. n>1000) call native_error("LOA: amplitude "//nm//": points must have 2..1000 entries")
      allocate(native_amp(i)%t(n),native_amp(i)%f(n))
      do j=1,n
        q=indexed(p//".points",j)
        if(integer_value(q//".count")/=2) call native_error("LOA: amplitude "//nm//": each point is [time, factor]")
        native_amp(i)%t(j)=number(q//"[1]")
        native_amp(i)%f(j)=number(q//"[2]")
        if(abs(native_amp(i)%f(j))>1.e6_real64 .or. abs(native_amp(i)%t(j))>1.e9_real64) &
          call native_error("LOA: amplitude "//nm//": value out of range")
        if(j>1) then
          if(native_amp(i)%t(j)<=native_amp(i)%t(j-1)) call native_error("LOA: amplitude "//nm//": time must be strictly increasing")
        endif
      enddo
    enddo
    ! ---- [[surface]] (edge sets, D5); with 1.set: the "surface ... edge" sets, in file order
    if(native_heat .and. .not.native_sets_on) call native_error("TEM: a transient_heat analysis needs [mesh] sets_file "// &
      "(convection faces are face sets in 1.set; T9)")
    if(native_sets_on .and. native_heat) then
      call heat_surfaces_from_sets(ns)
    elseif(native_sets_on) then
      call loa_surfaces_from_sets(ns)
    else
    ns=doc%count_of("surface")
    allocate(native_surf(ns))
    total=0
    do i=1,ns
      p=indexed("surface",i)
      total=total+integer_value(p//".edges.count")
    enddo
    allocate(native_loa_edges(5,max(total,1))); native_loa_edges=0
    native_loa_nedge=total
    total=0
    do i=1,ns
      p=indexed("surface",i)
      k=key_index(p//".@table")
      nm=string_value(p//".name")
      if(len(nm)==0) call native_error("LOA: empty surface name: "//p)
      do j=1,i-1
        if(trim(native_surf(j)%name)==nm) call native_error("LOA: duplicate surface name: "//nm)
      enddo
      native_surf(i)%name=nm
      if(string_value(p//".kind")/="edge2") call native_error("LOA: surface "//nm//": kind must be edge2")
      if(native_ndim/=2) call native_error("LOA: surface "//nm//": edge2 surfaces need mesh.dimension 2 (3D faces are deferred)")
      n=integer_value(p//".edges.count")
      if(n<1) call native_error("LOA: surface "//nm//" has no edges")
      native_surf(i)%first=total+1
      do j=1,n
        q=indexed(p//".edges",j)
        if(integer_value(q//".count")/=3) call native_error("LOA: surface "//nm//": edge is [node1, node2, element]")
        do k=1,3
          e(k)=integer_value(q//"["//itoa(k)//"]")
        enddo
        if(e(1)<1 .or. e(1)>native_points .or. e(2)<1 .or. e(2)>native_points) &
          call native_error("LOA: surface "//nm//": node out of range 1.."//itoa(native_points))
        if(e(1)==e(2)) call native_error("LOA: surface "//nm//": degenerate edge")
        if(e(3)<1 .or. e(3)>native_elements) call native_error("LOA: surface "//nm//": element out of range 1.."//itoa(native_elements))
        do ie=1,total
          if((native_loa_edges(1,ie)==e(1).and.native_loa_edges(2,ie)==e(2)) .or. &
             (native_loa_edges(1,ie)==e(2).and.native_loa_edges(2,ie)==e(1))) &
            call native_error("LOA: duplicate edge "//itoa(e(1))//"-"//itoa(e(2))//" in surface "//nm)
        enddo
        total=total+1
        native_loa_edges(1:2,total)=e(1:2); native_loa_edges(5,total)=e(3)
      enddo
      native_surf(i)%last=total
    enddo
    endif
    ! ---- [[step]] loads
    if(doc%count_of("step")/=native_blocks) call native_error("LOA: number of [[step]] must be "//itoa(native_blocks))
    allocate(native_loa(native_blocks))
    allocate(srf_step(max(ns,1)))
    do s=1,native_blocks
      sp=indexed("step",s)
      nl=doc%count_of(sp//".load")
      allocate(tmp(max(nl,1)))
      ng=0; nh=0; srf_step=.false.
      if(native_heat) then
        call heat_step_loads(s,sp,nl,ns,srf_step)
        native_loa(s)%hyd=tmp(1:0)
        deallocate(tmp)
        cycle
      endif
      do l=1,nl
        q=indexed(sp//".load",l)
        k=key_index(q//".@table")
        typ=string_value(q//".type")
        select case(typ)
        case("gravity")
          ng=ng+1
          if(ng>1) call native_error("LOA: step "//itoa(s)//": only one gravity load per step")
          g=number(q//".magnitude")
          if(g<=0 .or. g>1.e4_real64) call native_error('gravity out of range (0,1e4]')
          nd=native_ndim
          if(integer_value(q//".direction.count")/=nd) call native_error("LOA: step "//itoa(s)// &
            ": gravity direction needs "//itoa(nd)//" components (mesh.dimension)")
          d=0
          do j=1,nd
            d(j)=number(q//".direction["//itoa(j)//"]")
          enddo
          if(sum(d**2)==0) call native_error("LOA: step "//itoa(s)//": zero gravity direction")
          fd=d
          if(present_key(q//".fluid_direction.count")) then
            if(integer_value(q//".fluid_direction.count")/=nd) call native_error("LOA: step "//itoa(s)// &
              ": fluid_direction needs "//itoa(nd)//" components (mesh.dimension)")
            do j=1,nd
              fd(j)=number(q//".fluid_direction["//itoa(j)//"]")
            enddo
          endif
          if(.not.present_key(q//".amplitude")) call native_error("LOA: step "//itoa(s)//": gravity amplitude is required")
          native_loa(s)%gamp=amp_ref(q//".amplitude")
          if(string_value(q//".apply_to")/="all") call native_error("LOA: step "//itoa(s)//": gravity apply_to must be ""all""")
          native_loa(s)%g=g; native_loa(s)%dir=d; native_loa(s)%fdir=fd
          native_gravity(s)=g
        case("hydrostatic")
          nm=string_value(q//".surface")
          idx=0
          do i=1,ns
            if(trim(native_surf(i)%name)==nm) idx=i
          enddo
          if(idx==0) call native_error("LOA: step "//itoa(s)//": undefined surface "//nm)
          if(srf_step(idx)) call native_error("LOA: step "//itoa(s)//": surface "//nm//" loaded twice")
          srf_step(idx)=.true.; native_surf(idx)%used=.true.
          if(native_sets_on) i=sets_find(SK_SURFACE,nm)
          wl=number(q//".water_level"); bt=number(q//".bottom"); gw=number(q//".unit_weight")
          if(wl<=bt) call native_error("LOA: step "//itoa(s)//": water_level must be above bottom")
          if(gw<0 .or. gw>1.e6_real64) call native_error("LOA: step "//itoa(s)//": unit_weight out of range [0,1e6]")
          if(.not.present_key(q//".amplitude")) call native_error("LOA: step "//itoa(s)//": hydrostatic amplitude is required")
          ! vertical axis = legacy water (coordinate index): 2D always y; 3D faces need axis = "x" | "y" | "z"
          wa=2
          if(native_ndim==3) then
            if(.not.present_key(q//".axis")) call native_error("LOA: step "//itoa(s)//": hydrostatic on 3D faces needs "// &
              'axis = "x", "y" or "z" (the vertical coordinate)')
            typ=string_value(q//".axis")
            wa=index('xyz',typ)
            if(len(typ)/=1 .or. wa==0) call native_error("LOA: step "//itoa(s)//': hydrostatic axis must be "x", "y" or "z"')
          elseif(present_key(q//".axis")) then
            if(string_value(q//".axis")/="y") call native_error("LOA: step "//itoa(s)//': 2D hydrostatic axis must be "y"')
          endif
          nh=nh+1
          tmp(nh)=loa_hyd_t(native_surf(idx)%first,native_surf(idx)%last,amp_ref(q//".amplitude"),wl,bt,0._real64,wl-bt,gw,wa,0._real64)
        case("pressure")
          ! pressure pilot: uniform pressure (user 2026-10-16, clarified same day): face normal from node order
          ! (outward for the 1.set face/edge table forms); positive magnitude pushes INTO the element
          ! (Abaqus Dsload-P). Stored as magnitude; native_loa_press maps to legacy press = -magnitude
          nm=string_value(q//".surface")
          idx=0
          do i=1,ns
            if(trim(native_surf(i)%name)==nm) idx=i
          enddo
          if(idx==0) call native_error("LOA: step "//itoa(s)//": undefined surface "//nm)
          if(srf_step(idx)) call native_error("LOA: step "//itoa(s)//": surface "//nm//" loaded twice")
          srf_step(idx)=.true.; native_surf(idx)%used=.true.
          if(native_sets_on) i=sets_find(SK_SURFACE,nm)
          if(.not.present_key(q//".magnitude")) call native_error("LOA: step "//itoa(s)//": pressure magnitude is required")
          g=number(q//".magnitude")
          if(.not.ieee_is_finite(g) .or. abs(g)>1.e12_real64) call native_error("LOA: step "//itoa(s)//": pressure magnitude out of range")
          if(.not.present_key(q//".amplitude")) call native_error("LOA: step "//itoa(s)//": pressure amplitude is required")
          nh=nh+1
          tmp(nh)=loa_hyd_t(native_surf(idx)%first,native_surf(idx)%last,amp_ref(q//".amplitude"),0._real64,0._real64, &
            0._real64,0._real64,0._real64,0,g)
        case("concentrated")
          call native_error("LOA: step "//itoa(s)//": concentrated loads are global in HSTAR; define them in top-level [[load]]")
        case default
          call native_error("LOA: step "//itoa(s)//": unsupported load type "//typ)
        end select
      enddo
      if(ng==0) call native_error("LOA: step "//itoa(s)//": one gravity load is required")
      native_loa(s)%hyd=tmp(1:nh)
      deallocate(tmp)
    enddo
    do i=1,ns
      if(.not.native_surf(i)%used) call native_error("LOA: unreferenced surface "//trim(native_surf(i)%name))
    enddo
    ! ---- top-level [[load]]: concentrated (global over the whole analysis, legacy semantics)
    nnset=doc%count_of("nset")
    n=doc%count_of("load")
    allocate(native_pl(n), loa_nset_refs(n))
    do ip=1,n
      p=indexed("load",ip)
      k=key_index(p//".@table")
      if(string_value(p//".type")/="concentrated") call native_error("LOA: top-level [[load]] supports only type = ""concentrated""")
      nm=string_value(p//".nset")
      loa_nset_refs(ip)=nm
      if(native_sets_on) then
        idx=sets_find(SK_NSET,nm)
        if(idx==0) call native_error("LOA: load "//itoa(ip)//": undefined nset "//nm)
        native_pl(ip)%nodes=sets(idx)%m
        if(any(native_pl(ip)%nodes>native_points)) call native_error("LOA: node out of range in nset "//nm)
      else
      idx=0
      do i=1,nnset
        if(string_value(indexed("nset",i)//".name")==nm) idx=i
      enddo
      if(idx==0) call native_error("LOA: load "//itoa(ip)//": undefined nset "//nm)
      k=key_index(indexed("nset",idx)//".@table")
      k=integer_value(indexed("nset",idx)//".nodes.count")
      if(k<1) call native_error("LOA: nset "//nm//" has no nodes")
      allocate(native_pl(ip)%nodes(k))
      do j=1,k
        ii=integer_value(indexed(indexed("nset",idx)//".nodes",j))
        if(ii<1 .or. ii>native_points) call native_error("LOA: node out of range in nset "//nm)
        native_pl(ip)%nodes(j)=ii
      enddo
      endif
      if(integer_value(p//".components.count")/=native_ndim) call native_error("LOA: load "//itoa(ip)// &
        ": components must have mesh.dimension ("//itoa(native_ndim)//") entries")
      native_pl(ip)%comp=0
      do j=1,native_ndim
        native_pl(ip)%comp(j)=number(p//".components["//itoa(j)//"]")
      enddo
      native_pl(ip)%nudofn=2                        ! legacy native default (2D); 3D below
      if(native_ndim==3) native_pl(ip)%nudofn=3
      if(present_key(p//".moments.count")) then    ! beam pilot C2: nodal moments -> legacy pxyz rotation slots
        if(.not.native_beam) call native_error("LOA: load "//itoa(ip)//": moments need a B2 (beam) section")
        if(native_ndim==3) then
          if(integer_value(p//".moments.count")/=3) call native_error("LOA: load "//itoa(ip)//": moments needs 3 entries (mx, my, mz)")
          do j=1,3
            native_pl(ip)%comp(3+j)=number(p//".moments["//itoa(j)//"]")
          enddo
          native_pl(ip)%nudofn=6
        else
          if(integer_value(p//".moments.count")/=1) call native_error("LOA: load "//itoa(ip)//": moments needs 1 entry (mz) in 2D")
          native_pl(ip)%comp(4)=number(p//".moments[1]")
          native_pl(ip)%nudofn=4
        endif
      endif
      if(maxval(abs(native_pl(ip)%comp))>1.e12_real64) call native_error("LOA: load "//itoa(ip)//": component out of range")
      if(.not.present_key(p//".amplitude")) call native_error("LOA: load "//itoa(ip)//": amplitude is required")
      native_pl(ip)%itcurve=amp_ref(p//".amplitude")
    enddo
    ! ---- PRE contract mode uses curve 1 for every boundary amplitude (D9 plan b, static part)
    if(.not.native_pre_toml) then
      do s=1,native_blocks
        sp=indexed("step",s)
        do l=1,doc%count_of(sp//".boundary")
          q=indexed(sp//".boundary",l)
          if(present_key(q//".amplitude")) then
            k=amp_ref(q//".amplitude")
            if(k/=1) call native_error("LOA: boundary amplitude must be the first [[amplitude]] while HSTAR_NATIVE_PRE is contract")
            if(any(native_amp(1)%f/=1._real64) .or. native_amp(1)%t(1)>0) &
              call native_error("LOA: amplitude "//trim(native_amp(1)%name)//" used by contract-mode boundaries must be constant 1")
          endif
        enddo
      enddo
    endif
  end subroutine

  subroutine loa_check_used()
    integer :: i
    do i=1,size(native_amp)
      if(.not.native_amp(i)%used) call native_error("LOA: unreferenced amplitude "//trim(native_amp(i)%name))
    enddo
  end subroutine

  logical function loa_nset_referenced(nm)
    character(*), intent(in) :: nm
    integer :: i
    loa_nset_referenced=.false.
    if(.not.allocated(loa_nset_refs)) return
    do i=1,size(loa_nset_refs)
      if(trim(loa_nset_refs(i))==trim(nm)) loa_nset_referenced=.true.
    enddo
  end function

  ! ---- accessors used by Load.f90 (contract branch = the former profile constants)
  integer function native_loa_ncurves() result(n)
    n=1
    if(native_loa_toml) n=size(native_amp)
  end function

  subroutine native_loa_curve_head(ic,ntime,type_curve,nstoch,nline)
    integer, intent(in) :: ic
    integer, intent(out) :: ntime,nstoch,nline
    character(*), intent(out) :: type_curve
    ntime=2; type_curve='LINEAR'; nstoch=0; nline=2
    if(native_loa_toml) then
      ntime=size(native_amp(ic)%t)
      type_curve=native_amp(ic)%typ
    endif
  end subroutine

  subroutine native_loa_curve_times(ic,t)
    integer, intent(in) :: ic
    real(real64), intent(out) :: t(:)
    if(native_loa_toml) then
      t=native_amp(ic)%t
    else
      t=[0._real64,1._real64]
    endif
  end subroutine

  subroutine native_loa_curve_facts(ic,f)
    integer, intent(in) :: ic
    real(real64), intent(out) :: f(:)
    if(native_loa_toml) then
      f=native_amp(ic)%f
    else
      f=1._real64
    endif
  end subroutine

  ! D9 plan b (runtime part): when .pre data still use curve 1 implicitly, curve 1 must be constant 1.
  subroutine native_check_pre_curve1(type_curve,t,f)
    character(*), intent(in) :: type_curve
    real(real64), intent(in) :: t(:), f(:)
    if(.not.native_pre_curve1) return
    if(trim(type_curve)/='LINEAR' .or. any(f/=1._real64) .or. t(1)>0) &
      call native_error("PRE/LOA: curve 1 referenced by boundary data is not constant 1")
  end subroutine

  subroutine native_loa_pload_head(n,kpload)
    integer, intent(out) :: n,kpload
    n=0; kpload=0
    if(native_loa_toml) then
      n=size(native_pl)
      if(n>0) kpload=1
    endif
  end subroutine

  integer function native_loa_nedges() result(n)
    if(native_loa_toml) then
      n=native_loa_nedge
      if(native_heat) n=0                ! face4 surfaces belong to TEM (convection), not to .loa edges
    else
      n=merge(40,0,native_dam)
    endif
    loa_edge_group=0
  end function

  subroutine native_loa_edge_group(sedge,nnode,index,vdimn)
    integer, intent(out) :: sedge,nnode,index,vdimn
    sedge=40; nnode=2; index=1; vdimn=0
    if(native_loa_toml) then
      loa_edge_group=loa_edge_group+1
      sedge=native_surf(loa_edge_group)%last-native_surf(loa_edge_group)%first+1
      nnode=native_surf(loa_edge_group)%nnode; index=native_surf(loa_edge_group)%index
    endif
  end subroutine

  subroutine native_loa_edge(tedge,lnode,aelem)
    integer, intent(in) :: tedge
    integer, intent(out) :: lnode(:),aelem
    if(native_loa_toml) then
      lnode=native_loa_edges(1:size(lnode),tedge); aelem=native_loa_edges(5,tedge)
    else
      lnode=native_edges(1:2,tedge); aelem=native_edges(3,tedge)
    endif
  end subroutine

  subroutine native_loa_edge_load_head(iblk,nedgeload,ngroups)
    integer, intent(in) :: iblk
    integer, intent(out) :: nedgeload,ngroups
    integer :: i
    nedgeload=0; ngroups=0
    if(native_loa_toml) then
      ngroups=size(native_loa(iblk)%hyd)
      do i=1,ngroups
        nedgeload=nedgeload+native_loa(iblk)%hyd(i)%last-native_loa(iblk)%hyd(i)%first+1
      enddo
    elseif(native_dam .and. iblk==2) then
      nedgeload=40; ngroups=1
    endif
  end subroutine

  subroutine native_loa_hyd_head(iblk,ig,b,e,itc,water,code_load)
    integer, intent(in) :: iblk,ig
    integer, intent(out) :: b,e,itc,water,code_load
    b=1; e=40; itc=1; water=2; code_load=0
    if(native_loa_toml) then
      b=native_loa(iblk)%hyd(ig)%first; e=native_loa(iblk)%hyd(ig)%last; itc=native_loa(iblk)%hyd(ig)%itcurve
      water=native_loa(iblk)%hyd(ig)%water
    endif
  end subroutine

  ! water=0 rows (Load.f90 "read i0, press(edimn+1,1:nnode)"): TOML magnitude (positive = into the element)
  ! maps to legacy nodal press along the node-order normal as press = -magnitude
  subroutine native_loa_press(iblk,ig,pr)
    integer, intent(in) :: iblk,ig
    real(real64), intent(out) :: pr(:)
    if(.not.native_loa_toml) call native_error("LOA: per-face pressure rows need HSTAR_NATIVE_LOA=toml")
    pr=-native_loa(iblk)%hyd(ig)%press
  end subroutine

  subroutine native_loa_hyd_values(iblk,ig,c0,c1,p0,p1,fact)
    integer, intent(in) :: iblk,ig
    real(real64), intent(out) :: c0,c1,p0,p1,fact
    if(native_loa_toml) then
      c0=native_loa(iblk)%hyd(ig)%cor0; c1=native_loa(iblk)%hyd(ig)%cor1
      p0=native_loa(iblk)%hyd(ig)%p0; p1=native_loa(iblk)%hyd(ig)%p1; fact=native_loa(iblk)%hyd(ig)%fact
    else
      c0=50._real64; c1=0._real64; p0=0._real64; p1=50._real64; fact=native_pressure
    endif
  end subroutine

  subroutine native_loa_gravity_dirs(iblk,fg,ff)
    integer, intent(in) :: iblk
    real(real64), intent(out) :: fg(:),ff(:)
    if(native_loa_toml) then
      fg=native_loa(iblk)%dir(1:size(fg)); ff=native_loa(iblk)%fdir(1:size(ff))
    else
      fg=[0._real64,-1._real64]; ff=[0._real64,-1._real64]
    endif
  end subroutine

  subroutine native_loa_gravity_curves(iblk,tcg)
    integer, intent(in) :: iblk
    integer, intent(out) :: tcg(:)
    integer :: k
    k=1
    if(native_loa_toml) k=native_loa(iblk)%gamp
    tcg=k
    ! dam profile and GLB=toml: element sets not active in this step carry no gravity curve
    if(native_dam .or. native_glb_toml) tcg=merge(k,0,native_appearance(1:size(tcg),iblk)>0)
  end subroutine

  subroutine native_loa_deferred(what)
    character(*), intent(in) :: what
    call native_error("LOA: "//what//" is deferred in the LOA pilot and cannot be driven by TOML")
  end subroutine

  subroutine native_loa_open(u,iblk)
    integer, intent(out) :: u
    integer, intent(in) :: iblk
    if(iblk==0) then
      open(newunit=u,file="native-loa.txt",status="replace",action="write")
      write(u,"(a)") "native-loa v1 mode="//trim(merge("toml    ","contract",native_loa_toml))
    else
      open(newunit=u,file="native-loa.txt",status="old",position="append",action="write")
    endif
  end subroutine

  ! ---------------------------------------------------------------------------
  ! MAT pilot (HSTAR_NATIVE_MAT=toml): [[material]] (linear_elastic only) referenced by
  ! name from [[section]]. Design: docs/fortran-interface/mat-pilot/DESIGN.md (M1-M9).
  ! Material index (legacy imat) = [[material]] definition order (M3); thickness and
  ! gravity_factor (legacy ratio) live on [[section]] (M4); alfa optional, default 1e-5 (M5).
  ! ---------------------------------------------------------------------------
  logical function mat_owned(path)
    character(*), intent(in) :: path
    mat_owned=.false.
    if(.not.native_mat_toml) return
    mat_owned = index(path,"material[")==1 .or. index(path,"section[")==1 .or. index(path,"stress_strain_curve[")==1
  end function

  logical function in_list(seg,list)
    character(*), intent(in) :: seg, list
    in_list = index(" "//list//" "," "//seg//" ")>0
  end function

  ! Every key under one [[material]]/[[section]] table must be known; deferred features
  ! get a "deferred" message instead of the generic unsupported-key one.
  subroutine mat_check_keys(p,allowed,what)
    character(*), intent(in) :: p, allowed, what
    character(len=TOML_LEN_PATH) :: rest
    character(:), allocatable :: seg
    integer :: i, n, d
    character(*), parameter :: deferred = "plasticity plastic wetting liquefaction contact gap curve "// &
      "stress_strain hardening damage duncan_chang concrete cam_clay clay_pz sand_pz soil_pz "// &
      "fluid heat conductivity specific_heat geometry viscoelastic spring"
    n=len(p)+1
    do i=1,doc%n
      if(len_trim(doc%entry(i)%path)<=n) cycle
      if(doc%entry(i)%path(1:n)/=p//".") cycle
      rest=doc%entry(i)%path(n+1:)
      if(in_list(trim(rest),allowed)) then
        if(index(rest,"@table")>0) seen(i)=.true.
        cycle
      endif
      d=index(rest,".")
      seg=trim(rest); if(d>0) seg=rest(1:d-1)
      ! D3c: a table that is supported here (seg.@table allowed, e.g. section gap) reports the unknown key itself
      if(in_list(seg,deferred) .and. .not.in_list(seg//".@table",allowed)) call native_error("MAT: "//what//": "//seg// &
        " is deferred in the MAT pilot (only linear_elastic / goodman_fcm); use the legacy .mat")
      if(in_list(trim(rest),"E nu thermal_expansion")) call native_error("MAT: unsupported key: "// &
        trim(doc%entry(i)%path)//" (contract spelling; with HSTAR_NATIVE_MAT=toml use elastic.youngs_modulus,"// &
        " elastic.poisson_ratio, thermal.expansion_coefficient)")
      call native_error("MAT: unsupported key: "//trim(doc%entry(i)%path))
    enddo
  end subroutine

  subroutine parse_mat()
    logical :: cflag
    integer :: k2
    character(:), allocatable :: p, nm, s, el
    integer :: nmat, nsec, i, j, k, g, ng
    logical :: used(NATIVE_MAXG), covered(NATIVE_MAXG)
    character(len=TOML_LEN_TEXT) :: mname(NATIVE_MAXG), ename(NATIVE_MAXG), sname(NATIVE_MAXG)
    real(real64) :: v
    character(:), allocatable :: et
    ng=native_materials
    call parse_curves()
    nmat=doc%count_of("material")
    nsec=doc%count_of("section")
    if(nmat/=ng) call native_error("MAT: this profile needs exactly "//itoa(ng)// &
      " [[material]] table(s) (one per element set), found "//itoa(nmat))
    if(nsec/=ng) call native_error("MAT: this profile needs exactly "//itoa(ng)// &
      " [[section]] table(s) (one per element set), found "//itoa(nsec))
    do i=1,nmat
      p=indexed("material",i)
      k=key_index(p//".@table")
      nm=string_value(p//".name")
      if(len(nm)==0) call native_error("MAT: empty material name: "//p)
      do j=1,i-1
        if(trim(mname(j))==nm) call native_error("MAT: duplicate material name: "//nm)
      enddo
      mname(i)=nm
      if(native_heat) then
        call mat_check_keys(p,"@table name model diffusivity.count diffusivity[1] diffusivity[2] diffusivity[3] "// &
          "heat_source_curve placement_temperature_curve","material "//nm)
        if(string_value(p//".model")/="heat") call native_error("MAT: material "//nm// &
          ': a transient_heat analysis needs model "heat"')
        call parse_heat_mat(p,nm,i)
        cycle
      endif
      call mat_check_keys(p,"@table name model density elastic.@table elastic.youngs_modulus "// &
        "elastic.poisson_ratio thermal.@table thermal.expansion_coefficient creep.@table creep.model "// &
        "creep.kelvin_modulus creep.kelvin_viscosity creep.maxwell_viscosity goodman.@table goodman.softening "// &
        "goodman.fracture_energy goodman.tensile_strength mc_joint.@table mc_joint.cohesion mc_joint.hardening "// &
        "mc_joint.friction_angle mc_joint.dilation_angle mc_joint.tensile_strength mc_joint.cohesion_curve "// &
        "mc_joint.friction_curve mc_joint.dilation_curve mc_joint.tensile_curve janbu.@table janbu.k1 janbu.n janbu.kn "// &
        "janbu.ks_ref janbu.ks2_ref janbu.rf janbu.phi janbu.cohesion janbu.pa janbu.ft","material "//nm)
      s=string_value(p//".model")
      if(s=="goodman_janbu") native_mat_janbu(i)=.true.   ! D3c: CONTACT joint; the section must say contact = true
      if(s=="elastic_isotropic") call native_error("MAT: material "//nm// &
        ': model "elastic_isotropic" is the contract spelling; with HSTAR_NATIVE_MAT=toml use "linear_elastic"')
      if(s/="linear_elastic" .and. s/="goodman_fcm" .and. s/="mc_joint" .and. s/="goodman_janbu") &
        call native_error("MAT: material "//nm//": model "//s//" is deferred in the MAT pilot (only linear_elastic / "// &
        "goodman_fcm / goodman_janbu / mc_joint); use the legacy .mat")
      if(s/="goodman_janbu" .and. present_key(p//".janbu.@table")) call native_error("MAT: material "//nm// &
        ': [material.janbu] applies only to model "goodman_janbu"')
      if(s/="mc_joint" .and. present_key(p//".mc_joint.@table")) call native_error("MAT: material "//nm// &
        ': [material.mc_joint] applies only to model "mc_joint"')
      if(s/="goodman_fcm" .and. present_key(p//".goodman.@table")) call native_error("MAT: material "//nm// &
        ': [material.goodman] applies only to model "goodman_fcm"')
      native_density(i)=number(p//".density")
      if(native_density(i)<0 .or. native_density(i)>1.e6_real64) &   ! D3 (conv step A): density 0 allowed
        call native_error("MAT: material "//nm//": density out of range [0,1e6] kg/m^3")
      native_e(i)=number(p//".elastic.youngs_modulus")
      if(native_e(i)<=0 .or. native_e(i)>1.e15_real64) &
        call native_error("MAT: material "//nm//": youngs_modulus out of range (0,1e15] Pa")
      native_nu(i)=number(p//".elastic.poisson_ratio")
      if(native_nu(i)<=-1 .or. native_nu(i)>=0.5_real64) &
        call native_error("MAT: material "//nm//": poisson_ratio out of range (-1,0.5)")
      native_alfa(i)=1.e-5_real64
      if(present_key(p//".thermal.expansion_coefficient")) then
        native_alfa(i)=number(p//".thermal.expansion_coefficient")
        if(native_alfa(i)<0 .or. native_alfa(i)>1.e-2_real64) &
          call native_error("MAT: material "//nm//": expansion_coefficient out of range [0,1e-2] 1/K")
      endif
      call parse_creep(p,nm,i)
      if(s=="goodman_fcm") call parse_goodman_fcm(p,nm,i)
      if(s=="mc_joint") call parse_mc_joint(p,nm,i)
      if(s=="goodman_janbu") call parse_goodman_janbu(p,nm,i)
    enddo
    do k=1,native_nscurve
      if(.not.scurve_ref(k)) write(error_unit,"(a)") "MAT warning: stress_strain_curve "//trim(native_scurve(k)%name)// &
        " is not referenced by any material (kept as legacy curve "//itoa(k)//"; matline M2/M4)"
    enddo
    do j=1,ng
      if(native_sets_on) then
        ename(j)=native_kname(j)
      else
        ename(j)=string_value(indexed("elset",j)//".name")
      endif
    enddo
    used=.false.; covered=.false.
    do k=1,nsec
      p=indexed("section",k)
      i=key_index(p//".@table")
      nm=string_value(p//".name")
      if(len(nm)==0) call native_error("MAT: empty section name: "//p)
      do j=1,k-1
        if(trim(sname(j))==nm) call native_error("MAT: duplicate section name: "//nm)
      enddo
      sname(k)=nm
      call mat_check_keys(p,"@table name elset element formulation material thickness gravity_factor "// &
        "area torsion_constant inertia_y inertia_z contact layer_thickness gap.@table gap.initial gap.value "// &
        "orientation.@table orientation.node_pair.count orientation.node_pair[1] orientation.node_pair[2] "// &
        "orientation.vector.count orientation.vector[1] orientation.vector[2] orientation.vector[3] "// &
        "gap.tension_limit gap.use_tension_limit","section "//nm)
      el=string_value(p//".elset")
      g=0
      do j=1,ng
        if(trim(ename(j))==el) g=j
      enddo
      if(g==0) call native_error("MAT: section "//nm//": undefined elset "//el)
      if(native_sets_on) i=sets_find(SK_ELSET,el)
      if(covered(g)) call native_error("MAT: elset "//el//" is assigned by more than one section")
      covered(g)=.true.
      if(native_glb_toml) then
        ! element/formulation were checked against mesh.dimension by parse_glb (G7); L2 has no formulation
        s=string_value(p//".element")
        if(present_key(p//".formulation")) s=string_value(p//".formulation")
      else
      if(string_value(p//".element")/="Q4") call native_error("MAT: section "//nm//': only element "Q4" in this profile')
      if(string_value(p//".formulation")/="plane_strain") &
        call native_error("MAT: section "//nm//': only formulation "plane_strain" in this profile')
      endif
      s=string_value(p//".material")
      i=0
      do j=1,nmat
        if(trim(mname(j))==s) i=j
      enddo
      if(i==0) call native_error("MAT: section "//nm//": undefined material "//s)
      if(used(i)) call native_error("MAT: material "//s// &
        " is referenced by more than one section (one material per element set in the MAT pilot)")
      used(i)=.true.
      native_group_mat(g)=i
      et=string_value(p//".element")
      ! D3a/D3c: GOODMAN joints. goodman_janbu is the CONTACT usage (legacy: MECHANICAL CONTACT + igap0, Fem.f90:12261-12270,
      ! Stiff.f90:880); contact = true on its joint section -> legacy group name CONTACT + [section.gap] (igap0 card).
      cflag=.false.
      if(present_key(p//".contact")) then
        k2=key_index(p//".contact")
        if(doc%entry(k2)%kind/=TV_BOOL) call native_error("MAT: section "//nm//": contact must be true or false")
        cflag=doc%entry(k2)%lvalue
      endif
      if(native_mat_janbu(i)) then
        if(.not.cflag) call native_error("MAT: section "//nm//": material "//s// &
          ' is model "goodman_janbu", a contact joint: it needs contact = true on a joint section (legacy CONTACT group'// &
          ' with initial gap igap0); JANBU in a solid group is not valid')
        if(et/="joint") call native_error("MAT: section "//nm//': contact = true needs element "joint", '// &
          'formulation "goodman"')
        native_mat_contact(i)=.true.
        call parse_contact_gap(p,nm,i)
      elseif(cflag) then
        call native_error("MAT: section "//nm//': contact = true applies only to model "goodman_janbu" joints '// &
          '(MC joints in contact groups are not introduced yet)')
      elseif(present_key(p//".gap.@table")) then
        call native_error("MAT: section "//nm//": [section.gap] applies only to contact = true sections")
      endif
      if(et=="joint" .neqv. native_mat_goodman(i)) call native_error("MAT: section "//nm// &
        ': element "joint" and material model "goodman_fcm" / "goodman_janbu" go together (got element "'//et// &
        '", material '//s//')')
      ! D3b: mc_joint is an ordinary Q4 / B8 solid section; layer_thickness -> legacy group elcod_local (the
      ! thin-layer normal offset of the back nodes, Elements.f90 normal_local_inc)
      if(native_mat_mcj(i)) then
        if(.not.(et=="Q4" .or. et=="B8")) call native_error("MAT: section "//nm//": material "//s// &
          ' is model "mc_joint" (thin-layer solid): element must be "Q4" or "B8" (got "'//et//'")')
        if(.not.present_key(p//".layer_thickness")) call native_error("MAT: section "//nm// &
          ': model "mc_joint" needs layer_thickness (m, legacy elcod_local)')
        v=number(p//".layer_thickness")
        if(.not.(v>0 .and. v<=1.e3_real64)) call native_error("MAT: section "//nm//": layer_thickness out of range (0,1e3] m")
        native_elcod(g)=v
      elseif(present_key(p//".layer_thickness")) then
        ! user 2026-10-17: any Q4 / B8 solid section may be a thin layer (legacy elcod_local, normal_local_inc)
        if(.not.(et=="Q4" .or. et=="B8")) call native_error("MAT: section "//nm// &
          ': layer_thickness applies only to Q4 / B8 solid sections (got element "'//et//'")')
        v=number(p//".layer_thickness")
        if(.not.(v>0 .and. v<=1.e3_real64)) call native_error("MAT: section "//nm//": layer_thickness out of range (0,1e3] m")
        native_elcod(g)=v
      endif
      if(et=="L2" .or. et=="B2") then
        call section_geometry(p,nm,i,et)
      elseif(present_key(p//".area") .or. present_key(p//".torsion_constant") .or. present_key(p//".inertia_y") &
             .or. present_key(p//".inertia_z")) then
        call native_error("MAT: section "//nm//": area / torsion_constant / inertia_y / inertia_z apply only to "// &
          'line elements (element "L2" / "B2")')
      endif
      if(native_heat .and. (present_key(p//".thickness") .or. present_key(p//".gravity_factor"))) &
        call native_error("MAT: section "//nm//": thickness / gravity_factor do not apply to a transient_heat section")
      v=1._real64
      if(present_key(p//".thickness")) v=number(p//".thickness")
      if(present_key(p//".thickness") .and. native_ndim/=2) &
        call native_error("MAT: section "//nm//": thickness applies only to mesh.dimension 2")
      if(v<=0 .or. v>1.e3_real64) call native_error("MAT: section "//nm//": thickness out of range (0,1e3] m")
      native_thick(i)=v
      v=1._real64
      if(present_key(p//".gravity_factor")) v=number(p//".gravity_factor")
      if(v<0 .or. v>1.e3_real64) call native_error("MAT: section "//nm//": gravity_factor out of range [0,1e3]")
      native_ratio(i)=v
    enddo
    do i=1,nmat
      if(.not.used(i)) call native_error("MAT: unreferenced material "//trim(mname(i)))
    enddo
  end subroutine

  ! matline M5/M6: L2 section data. area (legacy aera) is required, J/Iy/Iz are optional (default 0, legacy
  ! decks write "A 0 0 0"); a bar takes its cross-section from area, so thickness is rejected.
  subroutine section_geometry(p,nm,i,et)
    character(*), intent(in) :: p, nm, et
    integer, intent(in) :: i
    real(real64) :: v
    if(.not.present_key(p//".area")) call native_error("MAT: section "//nm//': element "'//et//'" needs area (m^2)')
    if(present_key(p//".thickness")) call native_error("MAT: section "//nm// &
      ': thickness does not apply to element "'//et//'" (the line-element cross-section is area)')
    v=number(p//".area")
    if(.not.(v>0 .and. v<=1.e3_real64)) call native_error("MAT: section "//nm//": area out of range (0,1e3] m^2")
    native_area(i)=v
    native_tj(i)=geom_opt(p//".torsion_constant",nm)
    native_iy(i)=geom_opt(p//".inertia_y",nm)
    native_iz(i)=geom_opt(p//".inertia_z",nm)
    native_has_geom(i)=.true.
  end subroutine

  real(real64) function geom_opt(path,nm) result(v)
    character(*), intent(in) :: path, nm
    v=0._real64
    if(.not.present_key(path)) return
    v=number(path)
    if(.not.(v>=0 .and. v<=1.e3_real64)) call native_error("MAT: section "//nm//": "// &
      path(index(path,".",back=.true.)+1:)//" out of range [0,1e3] m^4")
  end function

  ! matline M1-M3: [[stress_strain_curve]] name, interpolation = "linear", strain = [...], stress = [...]
  ! -> legacy curve table scurves(1:nscurve) in definition order (parameter_find semantics unchanged, M3).
  ! M4: no TOML material model references a curve yet, so every curve is unreferenced -> warning (M2).
  subroutine parse_curves()
    character(:), allocatable :: p, nm
    integer :: nc, k, j, n, m, i
    real(real64) :: a
    nc=doc%count_of("stress_strain_curve")
    native_nscurve=nc
    allocate(native_scurve(nc))
    if(nc>0 .and. native_heat) call native_error("MAT: [[stress_strain_curve]] does not apply to a transient_heat analysis")
    do k=1,nc
      p=indexed("stress_strain_curve",k)
      i=key_index(p//".@table")
      call glb_check_keys(p,"@table name interpolation","strain stress")
      nm=string_value(p//".name")
      if(len(nm)==0) call native_error("MAT: empty stress_strain_curve name: "//p)
      do j=1,k-1
        if(trim(native_scurve(j)%name)==nm) call native_error("MAT: duplicate stress_strain_curve name: "//nm)
      enddo
      native_scurve(k)%name=nm
      if(string_value(p//".interpolation")/="linear") call native_error("MAT: stress_strain_curve "//nm// &
        ': interpolation must be "linear" (legacy type_curve LINEAR; other types are deferred)')
      n=integer_value(p//".strain.count"); m=integer_value(p//".stress.count")
      if(n<2) call native_error("MAT: stress_strain_curve "//nm//": needs at least 2 points")
      if(n>NATIVE_MAXC) call native_error("MAT: stress_strain_curve "//nm//": more than "//itoa(NATIVE_MAXC)//" points")
      if(m/=n) call native_error("MAT: stress_strain_curve "//nm//": strain and stress need the same number of points")
      allocate(native_scurve(k)%strain(n),native_scurve(k)%stress(n))
      do j=1,n
        a=number(indexed(p//".strain",j))
        if(.not.ieee_is_finite(a)) call native_error("MAT: stress_strain_curve "//nm//": non-finite strain")
        if(j>1) then
          if(a<=native_scurve(k)%strain(j-1)) call native_error("MAT: stress_strain_curve "//nm// &
            ": strain must be strictly increasing")
        endif
        native_scurve(k)%strain(j)=a
        a=number(indexed(p//".stress",j))
        if(.not.ieee_is_finite(a)) call native_error("MAT: stress_strain_curve "//nm//": non-finite stress")
        native_scurve(k)%stress(j)=a
      enddo
    enddo
    allocate(scurve_ref(nc)); scurve_ref=.false.
  end subroutine

  ! D3b: mc_joint curve keys name a [[stress_strain_curve]] (legacy scurves index, abscissa = plastic strain epstn)
  integer function scurve_ref_idx(path,nm) result(k)
    character(*), intent(in) :: path, nm
    character(:), allocatable :: c
    integer :: j
    k=0
    if(.not.present_key(path)) return
    c=string_value(path)
    do j=1,native_nscurve
      if(trim(native_scurve(j)%name)==c) k=j
    enddo
    if(k==0) call native_error("MAT: material "//nm//": undefined stress_strain_curve "//c//" ("//path//")")
    scurve_ref(k)=.true.
  end function

  ! Written from the stored scurves(i) (Material.f90).
  subroutine native_curve_receipt(i,typ,strain,stress)
    integer, intent(in) :: i
    character(*), intent(in) :: typ
    real(real64), intent(in) :: strain(:), stress(:)
    integer :: u, j
    if(.not.native_mat_receipt_on) return
    if(i==1) then
      open(newunit=u,file="native-mat-curves.txt",status="replace",action="write")
      write(u,"(a)") "native-mat-curves v1 mode="//trim(merge("toml    ","contract",native_mat_toml))
    else
      open(newunit=u,file="native-mat-curves.txt",status="old",position="append",action="write")
    endif
    write(u,"(a,i0,1x,a,1x,i0)") "curve ",i,trim(typ),size(strain)
    do j=1,size(strain)
      write(u,"(2es25.16)") strain(j),stress(j)
    enddo
    close(u)
  end subroutine

  ! Written from the stored props(imat)%geometry (Material.f90).
  subroutine native_geometry_receipt(i,aera,tj,iy,iz)
    integer, intent(in) :: i
    real(real64), intent(in) :: aera,tj,iy,iz
    integer :: u
    if(.not.native_mat_receipt_on) return
    open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    write(u,"(a,i0,a)") "geometry ",i," ROD_or_BEAM"
    write(u,"(a,4es25.16)") "aera,J,Iy,Iz=",aera,tj,iy,iz
    close(u)
  end subroutine

  ! [[material]] model = "heat" (T2): legacy 'HEAT' card alfa(1:ndimn), source_curve, place_curve,
  ! pipe_cooling=0 (equivalent pipe cooling deferred), ialfa=0. Curve names are resolved in heat_resolve.
  subroutine parse_heat_mat(p,nm,i)
    character(*), intent(in) :: p, nm
    integer, intent(in) :: i
    integer :: j
    native_mat_heat(i)=.true.
    if(integer_value(p//".diffusivity.count")/=native_ndim) call native_error("MAT: material "//nm// &
      ": diffusivity needs mesh.dimension ("//itoa(native_ndim)//") entries")
    do j=1,native_ndim
      native_diff(j,i)=number(indexed(p//".diffusivity",j))
      if(native_diff(j,i)<=0 .or. native_diff(j,i)>1.e3_real64) call native_error("MAT: material "//nm// &
        ": diffusivity out of range (0,1e3] m^2/time")
    enddo
    mat_src_name(i)=''
    if(present_key(p//".heat_source_curve")) mat_src_name(i)=string_value(p//".heat_source_curve")   ! D1: absent = no source (0)
    mat_place_name(i)=string_value(p//".placement_temperature_curve")
    native_density(i)=0; native_e(i)=0; native_nu(i)=0
  end subroutine

  ! Appended to native-mat.txt from the stored props(imat)%heat fields (Material.f90).
  subroutine native_heat_mat_receipt(i,alfa,src,place,pipe,ialfa)
    integer, intent(in) :: i,src,place,pipe,ialfa
    real(real64), intent(in) :: alfa(:)
    integer :: u
    if(.not.native_mat_receipt_on) return
    mat_receipt_n=mat_receipt_n+1
    if(mat_receipt_n==1) then
      open(newunit=u,file="native-mat.txt",status="replace",action="write")
      write(u,"(a)") "native-mat v1 mode="//trim(merge("toml    ","contract",native_mat_toml))
    else
      open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    endif
    write(u,"(a,i0,a)") "material ",i," HEAT"
    write(u,"(a,*(es25.16))") "diffusivity=",alfa
    write(u,"(a,4(1x,i0))") "source_curve,place_curve,pipe_cooling,ialfa=",src,place,pipe,ialfa
    close(u)
  end subroutine

  ! [[material]].creep (phase 1, C2): model "kelvin" | "burgers" -> legacy icreep=4
  ! (Ek, etak, etam). Kelvin stores etam = native_kelvin_etam (no Maxwell dashpot). Every
  ! other creep form (legacy icreep 1/2/3) is a hard error.
  ! D3a: [material.goodman] of model "goodman_fcm" -> legacy GOODMAN / FCM card (Material.f90:621-689):
  ! softening = legacy xlwmodel 1..4 (5 = user table, deferred), fracture_energy = Gf, tensile_strength = ft.
  ! kn0 = E*1e3, ks0 = E*1e3/(2(1+nu)) are derived by the legacy reader from youngs_modulus / poisson_ratio.
  subroutine parse_goodman_fcm(p,nm,i)
    character(*), intent(in) :: p, nm
    integer, intent(in) :: i
    real(real64) :: v
    if(.not.present_key(p//".goodman.@table")) call native_error("MAT: material "//nm// &
      ': model "goodman_fcm" needs [material.goodman] softening / fracture_energy / tensile_strength')
    if(present_key(p//".creep.@table")) call native_error("MAT: material "//nm//": creep does not apply to goodman_fcm")
    native_fcm_law(i)=integer_value(p//".goodman.softening")
    if(native_fcm_law(i)==5) call native_error("MAT: material "//nm//": goodman.softening = 5 (user table w0,w1,w2,ft,ft1)"// &
      " is deferred")
    if(native_fcm_law(i)<1 .or. native_fcm_law(i)>4) call native_error("MAT: material "//nm// &
      ": goodman.softening must be 1, 2, 3 or 4")
    ! D3a decision: law 3 uses w0 = 209e-6 regardless of Gf (Material.f90:678 @da3ebfc) -> accepted, with a warning
    if(native_fcm_law(i)==3) write(error_unit,"(a)") "MAT warning: material "//nm//": goodman.softening = 3 uses a fixed "// &
      "w0 = 209e-6 m (legacy); goodman.fracture_energy is read but has no effect"
    v=number(p//".goodman.fracture_energy")
    if(v<=0 .or. v>1.e9_real64) call native_error("MAT: material "//nm//": goodman.fracture_energy out of range (0,1e9] N/m")
    native_fcm_gf(i)=v
    v=number(p//".goodman.tensile_strength")
    if(v<=0 .or. v>1.e12_real64) call native_error("MAT: material "//nm//": goodman.tensile_strength out of range (0,1e12] Pa")
    native_fcm_ft(i)=v
    native_mat_goodman(i)=.true.
  end subroutine

  ! D3c: [material.janbu] of model "goodman_janbu" -> legacy GOODMAN / JANBU card (Material.f90:640-656):
  ! k1 = K1, n, kn = Kzz, ks_ref = Kzx, ks2_ref = Kzy (3D only, required there), rf = Rf, phi, cohesion = cohes,
  ! pa, ft = Ft. gamaw is not written (read but unused by the legacy PKPN, Stiff.f90:6089); point1 / point2 /
  ! uniax_cohes / frict_angle (JANBU header fields after the model name) are 0.
  subroutine parse_goodman_janbu(p,nm,i)
    character(*), intent(in) :: p, nm
    integer, intent(in) :: i
    integer :: i0
    if(.not.present_key(p//".janbu.@table")) call native_error("MAT: material "//nm// &
      ': model "goodman_janbu" needs [material.janbu] k1 / n / kn / ks_ref / rf / phi / cohesion / pa / ft')
    i0=key_index(p//".janbu.@table")
    if(present_key(p//".creep.@table")) call native_error("MAT: material "//nm//": creep does not apply to goodman_janbu")
    if(present_key(p//".goodman.@table")) call native_error("MAT: material "//nm// &
      ': [material.goodman] applies only to model "goodman_fcm"')
    native_janbu(1,i)=jv("k1",0._real64,1.e16_real64,.false.)
    native_janbu(2,i)=jv("n",-10._real64,10._real64,.true.)
    native_janbu(3,i)=jv("kn",0._real64,1.e16_real64,.false.)
    native_janbu(4,i)=jv("rf",0._real64,1._real64,.true.)
    native_janbu(5,i)=jv("phi",0._real64,90._real64,.true.)
    if(native_janbu(5,i)>=90) call native_error("MAT: material "//nm//": janbu.phi out of range [0,90) deg")
    native_janbu(6,i)=jv("ks_ref",0._real64,1.e16_real64,.true.)
    native_janbu(7,i)=jv("pa",0._real64,1.e9_real64,.false.)
    native_janbu(8,i)=0._real64
    native_janbu(9,i)=jv("cohesion",0._real64,1.e12_real64,.true.)
    native_janbu(10,i)=jv("ft",-1.e12_real64,1.e12_real64,.true.)
    if(native_ndim==3) then
      native_janbu(11,i)=jv("ks2_ref",0._real64,1.e16_real64,.true.)
    elseif(present_key(p//".janbu.ks2_ref")) then
      call native_error("MAT: material "//nm//": janbu.ks2_ref (legacy Kzy) applies only to mesh.dimension 3")
    endif
    native_mat_janbu(i)=.true.
    native_mat_goodman(i)=.true.
  contains
    real(real64) function jv(key,lo,hi,closed) result(v)
      character(*), intent(in) :: key
      real(real64), intent(in) :: lo, hi
      logical, intent(in) :: closed
      character(len=64) :: b
      if(.not.present_key(p//".janbu."//key)) call native_error("MAT: material "//nm//": missing janbu."//key)
      v=number(p//".janbu."//key)
      if(.not.(v<=hi .and. (v>lo .or. (closed .and. v>=lo)))) then
        write(b,"(a,es9.2,a,es9.2,a)") merge("[","(",closed),lo,",",hi,"]"
        call native_error("MAT: material "//nm//": janbu."//key//" out of range "//trim(b))
      endif
    end function
  end subroutine

  ! D3c: [section.gap] of a contact = true joint section -> legacy CONTACT card "igap0 gap0 ftcontact icft"
  ! (Material.f90:485-489). initial = "geometry" (default; igap0 99: gap = normal distance of the node pairs,
  ! natural_thickness = initial gap, Fem.f90:13113-13146) | "constant" (igap0 1, value = gap0). tension_limit = ftcontact,
  ! use_tension_limit = (icft /= 0), default true when tension_limit is given; without it the legacy open
  ! criterion uses 0.01 Pa (Fem.f90:13188-13192).
  subroutine parse_contact_gap(p,nm,i)
    character(*), intent(in) :: p, nm
    integer, intent(in) :: i
    character(:), allocatable :: s
    logical :: use
    integer :: k
    s="geometry"
    if(present_key(p//".gap.initial")) s=string_value(p//".gap.initial")
    select case(s)
    case("geometry")
      native_igap0(i)=99; native_gap(1,i)=0._real64
      if(present_key(p//".gap.value")) call native_error("MAT: section "//nm// &
        ': gap.value applies only to gap.initial = "constant"')
    case("constant")
      native_igap0(i)=1
      if(.not.present_key(p//".gap.value")) call native_error("MAT: section "//nm// &
        ': gap.initial = "constant" needs gap.value (m)')
      native_gap(1,i)=number(p//".gap.value")
      if(.not.(abs(native_gap(1,i))<=1._real64)) call native_error("MAT: section "//nm//": gap.value out of range [-1,1] m")
    case default
      call native_error("MAT: section "//nm//': gap.initial must be "geometry" or "constant" (legacy igap0 99 / 1; '// &
        'the circular gap table igap0 = 2 is not introduced yet)')
    end select
    use=present_key(p//".gap.tension_limit")
    if(present_key(p//".gap.use_tension_limit")) then
      k=key_index(p//".gap.use_tension_limit")
      if(doc%entry(k)%kind/=TV_BOOL) call native_error("MAT: section "//nm//": gap.use_tension_limit must be true or false")
      use=doc%entry(k)%lvalue
      if(use .and. .not.present_key(p//".gap.tension_limit")) call native_error("MAT: section "//nm// &
        ": gap.use_tension_limit = true needs gap.tension_limit (Pa)")
      if(.not.use .and. present_key(p//".gap.tension_limit")) call native_error("MAT: section "//nm// &
        ": gap.tension_limit is given but gap.use_tension_limit = false")
    endif
    native_icft(i)=merge(1,0,use)
    native_gap(2,i)=0._real64
    if(use) native_gap(2,i)=number(p//".gap.tension_limit")
    if(.not.(abs(native_gap(2,i))<=1.e12_real64)) call native_error("MAT: section "//nm// &
      ": gap.tension_limit out of range [-1e12,1e12] Pa")
  end subroutine

  ! D3c receipt: stored props(imat)%mechanical%solid CONTACT card + Goodman JANBU fields (Material.f90)
  subroutine native_janbu_receipt(i,igap0,gap0,ftc,icft,v)
    integer, intent(in) :: i, igap0, icft
    real(real64), intent(in) :: gap0, ftc, v(:)
    integer :: u
    if(.not.native_mat_receipt_on) return
    open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    write(u,"(a,i0,a)") "goodman ",i," JANBU name=CONTACT"
    write(u,"(a,i0,es25.16,es25.16,i3)") "igap0,gap0,ftcontact,icft= ",igap0,gap0,ftc,icft
    write(u,"(a,5es25.16)") "K1,n,Kzz,Rf,phi=",v(1:5)
    write(u,"(a,5es25.16)") "Kzx,pa,gamaw,cohes,Ft=",v(6:10)
    write(u,"(a,es25.16)") "Kzy=",v(11)
    close(u)
  end subroutine

  ! D3b: [material.mc_joint] of model "mc_joint" -> legacy CLASSICALEP / MCJOINT card (Material.f90 case CLASSICALEP):
  ! cohesion = sigma0, hardening (d c / d epstn, default 0), friction_angle / dilation_angle (deg), tensile_strength = ft
  ! (default 0), *_curve = csigma0 / cfrict / cdilan / cft. sigmat / csigmat are written 0 (inert in the legacy Q4/B8 path).
  subroutine parse_mc_joint(p,nm,i)
    character(*), intent(in) :: p, nm
    integer, intent(in) :: i
    real(real64) :: v
    integer :: i0
    if(.not.native_glb_toml) call native_error("MAT: material "//nm//': model "mc_joint" needs HSTAR_NATIVE_GLB=toml')
    if(.not.present_key(p//".mc_joint.@table")) call native_error("MAT: material "//nm// &
      ': model "mc_joint" needs [material.mc_joint] cohesion / friction_angle / dilation_angle')
    i0=key_index(p//".mc_joint.@table")
    if(present_key(p//".creep.@table")) call native_error("MAT: material "//nm//": creep does not apply to mc_joint")
    v=number(p//".mc_joint.cohesion")
    if(.not.(v>0 .and. v<=1.e12_real64)) call native_error("MAT: material "//nm//": mc_joint.cohesion out of range (0,1e12] Pa")
    native_mcj(1,i)=v
    v=0
    if(present_key(p//".mc_joint.hardening")) v=number(p//".mc_joint.hardening")
    if(.not.(abs(v)<=1.e15_real64)) call native_error("MAT: material "//nm//": mc_joint.hardening out of range [-1e15,1e15] Pa")
    native_mcj(2,i)=v
    v=number(p//".mc_joint.friction_angle")
    if(.not.(v>=0 .and. v<90)) call native_error("MAT: material "//nm//": mc_joint.friction_angle out of range [0,90) deg")
    native_mcj(3,i)=v
    v=number(p//".mc_joint.dilation_angle")
    if(.not.(v>=0 .and. v<90)) call native_error("MAT: material "//nm//": mc_joint.dilation_angle out of range [0,90) deg")
    native_mcj(4,i)=v
    v=0
    if(present_key(p//".mc_joint.tensile_strength")) v=number(p//".mc_joint.tensile_strength")
    if(.not.(v>=0 .and. v<=1.e12_real64)) call native_error("MAT: material "//nm// &
      ": mc_joint.tensile_strength out of range [0,1e12] Pa")
    native_mcj(5,i)=v
    native_mcj_curve(1,i)=scurve_ref_idx(p//".mc_joint.cohesion_curve",nm)
    ! D3b: curve-driven friction / dilation / tensile strength (legacy cfrict / cdilan / cft) are not introduced yet
    ! (user decision 2026-10-16); constant friction_angle / dilation_angle, cohesion via hardening or cohesion_curve
    if(present_key(p//".mc_joint.friction_curve") .or. present_key(p//".mc_joint.dilation_curve") .or. &
       present_key(p//".mc_joint.tensile_curve")) call native_error("MAT: material "//nm//": mc_joint.friction_curve / "// &
      "dilation_curve / tensile_curve are deferred (not introduced yet; use friction_angle / "// &
      "dilation_angle, cohesion_curve)")
    if(native_mcj(5,i)>0) write(error_unit,"(a)") "MAT warning: material "//nm// &
      ": mc_joint.tensile_strength has no effect (legacy Q4/B8 MCJOINT has no active tension"// &
      " cutoff; the Coulomb apex c/tan(phi) limits tension)"
    native_mat_mcj(i)=.true.
    any_mcj=.true.
  end subroutine

  ! D3b receipt: the stored props(imat)%mechanical%solid%ClassicalEP fields (after kstab reduction) and the globals
  subroutine native_mcj_receipt(i,sigma0,hard,frict,dilan,ft,cs,cf,cd,ct,kstab,matc,tload,tnl)
    integer, intent(in) :: i,cs,cf,cd,ct,matc,tnl
    real(real64), intent(in) :: sigma0,hard,frict,dilan,ft,kstab
    character(*), intent(in) :: tload
    integer :: u
    if(.not.native_mat_receipt_on) return
    open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    write(u,"(a,i0,a)") "classicalep ",i," MCJOINT"
    write(u,"(a,5es25.16)") "sigma0,hardening,frict,dilan,ft=",sigma0,hard,frict,dilan,ft
    write(u,"(a,4(1x,i0))") "csigma0,cfrict,cdilan,cft=",cs,cf,cd,ct
    write(u,"(a,es25.16,1x,i0,1x,a,1x,i0)") "kstab,mat_curve,type_load,type_nl=",kstab,matc,tload,tnl
    close(u)
  end subroutine

  ! D3a receipt: the stored props(imat)%mechanical%solid%Goodman%fcmp fields (Material.f90)
  subroutine native_goodman_receipt(i,law,gf,ft,w0,w1,ft1,kns0)
    integer, intent(in) :: i,law
    real(real64), intent(in) :: gf,ft,w0,w1,ft1,kns0(:)
    integer :: u
    if(.not.native_mat_receipt_on) return
    open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    write(u,"(a,i0,a,i0)") "goodman ",i," FCM softening=",law
    write(u,"(a,es25.16)") "Gf=",gf
    write(u,"(a,es25.16)") "ft=",ft
    write(u,"(a,es25.16)") "w0=",w0
    if(law==2) write(u,"(a,2es25.16)") "w1,ft1=",w1,ft1
    write(u,"(a,3es25.16)") "kns0=",kns0
    close(u)
  end subroutine

  subroutine parse_creep(p,nm,i)
    character(*), intent(in) :: p, nm
    integer, intent(in) :: i
    character(:), allocatable :: s
    native_icreep(i)=0
    if(.not.present_key(p//".creep.model")) then
      if(present_key(p//".creep.@table") .or. present_key(p//".creep.kelvin_modulus") .or. &
         present_key(p//".creep.kelvin_viscosity") .or. present_key(p//".creep.maxwell_viscosity")) &
        call native_error("MAT: material "//nm//": missing creep.model (kelvin or burgers)")
      return
    endif
    s=string_value(p//".creep.model")
    select case(s)
    case("kelvin","burgers")
    case default
      call native_error("MAT: material "//nm//": creep model "//s//" is deferred (creep pilot phase 1 "// &
        "supports only kelvin and burgers; legacy icreep 1/2/3 forms use the legacy .mat)")
    end select
    native_icreep(i)=4
    native_ek(i)=number(p//".creep.kelvin_modulus")
    if(native_ek(i)<=0 .or. native_ek(i)>1.e15_real64) &
      call native_error("MAT: material "//nm//": creep.kelvin_modulus out of range (0,1e15] Pa")
    native_etak(i)=number(p//".creep.kelvin_viscosity")
    if(native_etak(i)<=0 .or. native_etak(i)>1.e25_real64) &
      call native_error("MAT: material "//nm//": creep.kelvin_viscosity out of range (0,1e25] Pa*time")
    if(s=="kelvin") then
      if(present_key(p//".creep.maxwell_viscosity")) call native_error("MAT: material "//nm// &
        ": creep.maxwell_viscosity is not allowed for model kelvin (use model burgers)")
      native_etam(i)=native_kelvin_etam
    else
      native_etam(i)=number(p//".creep.maxwell_viscosity")
      if(native_etam(i)<=0 .or. native_etam(i)>1.e25_real64) &
        call native_error("MAT: material "//nm//": creep.maxwell_viscosity out of range (0,1e25] Pa*time")
    endif
  end subroutine

  ! Appended to native-mat.txt from the stored props(imat)%mechanical%solid%creep fields.
  subroutine native_creep_receipt(i,ek,etak,etam)
    integer, intent(in) :: i
    real(real64), intent(in) :: ek,etak,etam
    integer :: u
    if(.not.native_mat_receipt_on) return
    open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    write(u,"(a,es25.16)") "creep_Ek=",ek
    write(u,"(a,es25.16)") "creep_etak=",etak
    write(u,"(a,es25.16)") "creep_etam=",etam
    close(u)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! MAN seam (phase 1, C4): HSTAR_NATIVE_MAN=toml owns step[i].controls.time_increment,
  ! tolerance_force, tolerance_dof (static_U increment card: ditime, toler_force, toler_var).
  ! ---------------------------------------------------------------------------
  logical function man_owned(path)
    character(*), intent(in) :: path
    integer :: d
    man_owned=.false.
    if(.not.native_man_toml) return
    if(index(path,"step[")/=1) return
    d=index(path,"].controls.")
    if(d==0) return
    man_owned = in_list(trim(path(d+11:)),"time_increment tolerance_force tolerance_dof")
  end function

  ! transient_heat (T3): [[step.increment]] = one legacy time_dependent increment card.
  subroutine parse_man_heat()
    integer :: i, j, n, k
    character(:), allocatable :: p, q
    do i=1,native_blocks
      p=indexed("step",i)
      n=doc%count_of(p//".increment")
      if(n<1 .or. n>NATIVE_MAXI) call native_error("MAN: "//p//" needs 1.."//itoa(NATIVE_MAXI)//" [[step.increment]] tables")
      native_ninc(i)=n
      do j=1,n
        q=indexed(p//".increment",j)
        k=key_index(q//".@table")
        call glb_check_keys(q,"@table time_increment substeps max_iterations tolerance","")
        native_inc_dt(j,i)=number(q//".time_increment")
        if(native_inc_dt(j,i)<=0 .or. native_inc_dt(j,i)>1.e6_real64) call native_error("MAN: "//q//".time_increment out of range (0,1e6]")
        native_inc_nstep(j,i)=integer_value(q//".substeps")
        if(native_inc_nstep(j,i)<1 .or. native_inc_nstep(j,i)>100000) call native_error("MAN: "//q//".substeps out of range [1,100000]")
        native_inc_miter(j,i)=integer_value(q//".max_iterations")
        if(native_inc_miter(j,i)<1 .or. native_inc_miter(j,i)>100) call native_error("MAN: "//q//".max_iterations out of range [1,100]")
        native_inc_tol(j,i)=number(q//".tolerance")
        if(native_inc_tol(j,i)<=0 .or. native_inc_tol(j,i)>1._real64) call native_error("MAN: "//q//".tolerance out of range (0,1]")
      enddo
    enddo
  end subroutine

  subroutine parse_man()
    integer :: i
    character(:), allocatable :: p
    if(native_heat) then
      call parse_man_heat()
      return
    endif
    do i=1,native_blocks
      p=indexed("step",i)//".controls"
      native_ditime(i)=number(p//".time_increment")
      if(native_ditime(i)<=0 .or. native_ditime(i)>1.e6_real64) &
        call native_error("MAN: "//p//".time_increment out of range (0,1e6]")
      native_tol_force(i)=number(p//".tolerance_force")
      if(native_tol_force(i)<=0 .or. native_tol_force(i)>1._real64) &
        call native_error("MAN: "//p//".tolerance_force out of range (0,1]")
      native_tol_dof(i)=number(p//".tolerance_dof")
      if(native_tol_dof(i)<=0 .or. native_tol_dof(i)>1._real64) &
        call native_error("MAN: "//p//".tolerance_dof out of range (0,1]")
    enddo
  end subroutine

  subroutine native_man_deferred(what)
    character(*), intent(in) :: what
    call native_error("MAN: "//what//" is deferred in the creep pilot and cannot be driven by TOML")
  end subroutine

  ! Written from the stored static_U increment variables (Fem.f90), once per increment card.
  subroutine native_man_receipt(iblk,iincs,miter,ditime,nstep,inc_step,noutn,noutf,toler_force,toler_var)
    integer, intent(in) :: iblk,iincs,miter,nstep,inc_step,noutn,noutf
    real(real64), intent(in) :: ditime,toler_force,toler_var(:)
    integer :: u
    if(.not.native_man_receipt_on) return
    man_receipt_n=man_receipt_n+1
    if(man_receipt_n==1) then
      open(newunit=u,file="native-man.txt",status="replace",action="write")
      write(u,"(a)") "native-man v1 mode="//trim(merge("toml    ","contract",native_man_toml))
    else
      open(newunit=u,file="native-man.txt",status="old",position="append",action="write")
    endif
    write(u,"(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)") "block=",iblk," incr=",iincs," miter=",miter," nstep=",nstep, &
      " inc_step=",inc_step," noutn=",noutn," noutf=",noutf
    write(u,"(a,es25.16)") "ditime=",ditime
    write(u,"(a,es25.16)") "toler_force=",toler_force
    write(u,"(a,*(es25.16))") "toler_var=",toler_var
    close(u)
  end subroutine

  subroutine native_mat_deferred(what)
    character(*), intent(in) :: what
    call native_error("MAT: "//what//" is deferred in the MAT pilot and cannot be driven by TOML")
  end subroutine

  ! Written from the stored props(imat)%mechanical%solid fields (Material.f90).
  subroutine native_mat_receipt(i,material,rho,ratio,thick,e,nu,alfa,icreep,kind_wt,jliqu,ie,inu,rho_w)
    integer, intent(in) :: i,icreep,kind_wt,jliqu,ie,inu
    character(*), intent(in) :: material
    real(real64), intent(in) :: rho,ratio,thick,e,nu,alfa,rho_w
    integer :: u
    if(.not.native_mat_receipt_on) return
    mat_receipt_n=mat_receipt_n+1
    if(mat_receipt_n==1) then
      open(newunit=u,file="native-mat.txt",status="replace",action="write")
      write(u,"(a)") "native-mat v1 mode="//trim(merge("toml    ","contract",native_mat_toml))
    else
      open(newunit=u,file="native-mat.txt",status="old",position="append",action="write")
    endif
    write(u,"(a,i0,1x,a)") "material ",i,trim(material)
    write(u,"(a,es25.16)") "density=",rho
    write(u,"(a,es25.16)") "ratio=",ratio
    write(u,"(a,es25.16)") "thickness=",thick
    write(u,"(a,es25.16)") "E=",e
    write(u,"(a,es25.16)") "nu=",nu
    write(u,"(a,es25.16)") "alfa=",alfa
    write(u,"(a,5(1x,i0))") "icreep,kind_wt,jliqu,iE,iNu=",icreep,kind_wt,jliqu,ie,inu
    write(u,"(a,es25.16)") "density_w=",rho_w
    close(u)
  end subroutine

  ! Written from the stored group(igroup)%matno (Global.f90).
  subroutine native_mat_group_receipt(igroup,matno)
    integer, intent(in) :: igroup, matno
    integer :: u
    if(.not.native_mat_receipt_on) return
    if(igroup==1) then
      open(newunit=u,file="native-mat-groups.txt",status="replace",action="write")
    else
      open(newunit=u,file="native-mat-groups.txt",status="old",position="append",action="write")
    endif
    write(u,"(a,i0,a,i0)") "group ",igroup," matno=",matno
    close(u)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! TEM minimal step (HSTAR_NATIVE_TEM=toml, decisions T1/M8): a TOML case may not carry
  ! any temperature loading; the receipt records the counts boundt actually stored.
  ! ---------------------------------------------------------------------------
  subroutine parse_tem()
    integer :: i
    character(len=TOML_LEN_PATH) :: p
    do i=1,doc%n
      p=doc%entry(i)%path
      if(index(p,"material[")==1) cycle
      if(index(p,"thermal")==1 .or. index(p,"initial_temperature")==1 .or. index(p,"temperature_field")==1 &
         .or. index(p,".thermal")>0 .or. index(p,"temperature")>0) &
        call native_error("TEM: temperature loading ("//trim(p)//") is deferred; HSTAR_NATIVE_TEM=toml "// &
          "only accepts a case without temperature loading (use the legacy .tem)")
    enddo
  end subroutine

  subroutine native_tem_deferred(what)
    character(*), intent(in) :: what
    call native_error("TEM: "//what//" is deferred (TEM minimal step) and cannot be driven by TOML")
  end subroutine

  subroutine native_tem_receipt(nsurf,nedge,nelgroup,npipe,ecwpipe,outintr)
    integer, intent(in) :: nsurf,nedge,nelgroup,npipe,ecwpipe,outintr
    integer :: u
    if(.not.native_tem_receipt_on) return
    tem_receipt_n=tem_receipt_n+1
    if(tem_receipt_n==1) then
      open(newunit=u,file="native-tem.txt",status="replace",action="write")
      write(u,"(a)") "native-tem v1 mode="//trim(merge("toml    ","contract",native_tem_toml))
    else
      open(newunit=u,file="native-tem.txt",status="old",position="append",action="write")
    endif
    write(u,"(a,i0,6(a,i0))") "call=",tem_receipt_n," ntemp_surface=",nsurf," ntedge=",nedge, &
      " ntelgroup=",nelgroup," npipe=",npipe," ecwpipe=",ecwpipe," outintr=",outintr
    close(u)
  end subroutine

  ! ---------------------------------------------------------------------------
  ! GLB pilot (HSTAR_NATIVE_GLB=toml, glb-pilot/DESIGN.md G1-G12): global control from TOML.
  ! Counts (nodes, elements, groups, materials, steps, dofs) are derived, never keys (G4).
  ! ---------------------------------------------------------------------------
  ! beam pilot C2: section.orientation = { node_pair = [p1, p2] } -> legacy group point_direct (3D only; local z =
  ! coord(p2) - coord(p1), used as given, Elements.f90 direct_beam). Omitted -> 0 0 (legacy default rule).
  subroutine glb_beam_orientation(p,g)
    character(*), intent(in) :: p
    integer, intent(in) :: g
    integer :: j
    native_pdir(:,g)=0
    if(.not.present_key(p//'.orientation.@table')) return
    if(native_ndim/=3) call native_error('GLB: section '//p//': orientation applies only to mesh.dimension 3 '// &
      '(2D beams take the in-plane normal)')
    if(present_key(p//'.orientation.vector.count')) call native_error('GLB: section '//p// &
      ': orientation.vector is deferred (the solver orients a beam from two nodes, legacy point_direct; use node_pair)')
    if(.not.present_key(p//'.orientation.node_pair.count')) call native_error('GLB: section '//p// &
      ': orientation needs node_pair = [p1, p2]')
    if(integer_value(p//'.orientation.node_pair.count')/=2) call native_error('GLB: section '//p// &
      ': orientation.node_pair needs exactly 2 node ids')
    do j=1,2
      native_pdir(j,g)=integer_value(p//'.orientation.node_pair['//itoa(j)//']')
      if(native_pdir(j,g)<1) call native_error('GLB: section '//p//': orientation.node_pair node must be >= 1')
    enddo
    if(native_pdir(1,g)==native_pdir(2,g)) call native_error('GLB: section '//p// &
      ': orientation.node_pair needs two different nodes')
  end subroutine

  subroutine glb_check_keys(p,allowed,lists)
    character(*), intent(in) :: p, allowed, lists
    character(len=TOML_LEN_PATH) :: rest
    integer :: i, n, b
    n=len(p)+1
    do i=1,doc%n
      if(len_trim(doc%entry(i)%path)<=n) cycle
      if(doc%entry(i)%path(1:n)/=p//".") cycle
      rest=doc%entry(i)%path(n+1:)
      if(in_list(trim(rest),allowed)) cycle
      b=index(rest,"[")
      if(b==0) b=index(rest,".count")
      if(b>1) then
        if(in_list(rest(1:b-1),lists)) cycle
      endif
      call native_error("GLB: unsupported key: "//trim(doc%entry(i)%path)// &
        " (counts are derived from 1.cor/1.ele and the TOML tables; advanced switches are deferred)")
    enddo
  end subroutine

  integer function glb_count_rows(file) result(n)
    character(*), intent(in) :: file
    integer :: u, ios
    character(1024) :: row
    open(newunit=u,file=file,status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh '//file)
    n=0
    do
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) exit
      if(len_trim(row)==0) call native_error('GLB: blank row in mesh '//file)
      n=n+1
      if(n>NATIVE_MAXN) call native_error('GLB: mesh '//file//' exceeds '//itoa(NATIVE_MAXN)//' rows')
    enddo
    close(u)
  end function

  ! ---- SETS pilot helpers -------------------------------------------------
  subroutine sets_init()
    character(:), allocatable :: f, h, msg, dir
    integer :: j
    if(.not.native_glb_toml) call native_error("SETS: mesh.sets_file needs HSTAR_NATIVE_GLB=toml")
    if(doc%count_of("nset")>0) call native_error("SETS: inline [[nset]] is not allowed with mesh.sets_file (define it in 1.set)")
    if(doc%count_of("surface")>0) call native_error("SETS: inline [[surface]] is not allowed with mesh.sets_file (define it in 1.set)")
    if(doc%count_of("elset")>0) call native_error("SETS: inline [[elset]] is not allowed with mesh.sets_file (define it in 1.set)")
    f=string_value('mesh.sets_file')
    if(len(f)==0 .or. len(f)>64 .or. scan(f,'/'//achar(92))>0) &
      call native_error('SETS: mesh.sets_file must be a plain file name next to the TOML (e.g. "1.set")')
    if(present_key('mesh.sets_sha256')) then
      h=string_value('mesh.sets_sha256')          ! checked by the Python tools (S8)
      if(len(h)/=64 .or. verify(h,'0123456789abcdef')/=0) call native_error('SETS: mesh.sets_sha256 must be 64 lowercase hex digits')
    endif
    j=len_trim(native_mesh_prefix)-1
    dir=native_mesh_prefix(1:j)
    call sets_load(trim(native_mesh_prefix),dir//f,msg)
    if(len(msg)>0) call native_error(msg)
    native_sets_on=.true.
  end subroutine

  ! elsets of 1.set = element groups, in file order; each a contiguous range, together tiling 1..nelem
  subroutine glb_elsets_from_sets(ng,ename,total)
    integer, intent(out) :: ng, total
    character(len=TOML_LEN_TEXT), intent(inout) :: ename(:)
    integer :: i, j, k
    character(:), allocatable :: nm
    ng=0; total=0
    do i=1,sets_n
      if(sets(i)%kind/=SK_ELSET) cycle
      ng=ng+1
      if(ng>NATIVE_MAXG) call native_error('GLB: need 1..'//itoa(NATIVE_MAXG)//' elsets')
      nm=trim(sets(i)%name)
      if(len(nm)>20) call native_error('GLB: elset name longer than 20 characters: '//nm)
      do k=1,sets(i)%n
        if(sets(i)%m(k)/=total+k) call native_error('SETS: elset '//nm//' must be the contiguous element range '// &
          itoa(total+1)//'..'//itoa(total+sets(i)%n)//' (elsets follow the element order of 1.ele)')
      enddo
      ename(ng)=nm; native_kname(ng)=nm
      native_elcounts(ng)=sets(i)%n
      total=total+sets(i)%n
    enddo
    if(ng<1) call native_error('GLB: need 1..'//itoa(NATIVE_MAXG)//' elsets in 1.set')
    j=0
  end subroutine

  ! transient_heat: the "surface ... face" sets of 1.set (B8 faces, canonical node order of yl_sets) in file order
  subroutine heat_surfaces_from_sets(ns)
    integer, intent(out) :: ns
    integer :: i, j, k, total, f(4)
    ns=0; total=0
    do i=1,sets_n
      if(sets(i)%kind/=SK_SURFACE) cycle
      if(sets(i)%sub/=SS_FACE) call native_error("TEM: surface "//trim(sets(i)%name)// &
        ": a transient_heat analysis takes face surfaces only")
      ns=ns+1; total=total+sets(i)%n
    enddo
    if(native_ndim/=3) call native_error("TEM: face surfaces need mesh.dimension 3")
    allocate(native_surf(ns), native_faces(5,max(total,1)), native_loa_edges(5,1))
    native_loa_nedge=0; native_nface=0; k=0
    do i=1,sets_n
      if(sets(i)%kind/=SK_SURFACE) cycle
      k=k+1
      native_surf(k)%name=sets(i)%name; native_surf(k)%kind=4; native_surf(k)%first=native_nface+1
      do j=1,sets(i)%n
        call sets_face_nodes(sets(i)%m(2*j-1),sets(i)%m(2*j),f)
        native_nface=native_nface+1
        native_faces(1:4,native_nface)=f; native_faces(5,native_nface)=sets(i)%m(2*j-1)
      enddo
      native_surf(k)%last=native_nface
    enddo
  end subroutine

  subroutine loa_surfaces_from_sets(ns)
    ! 2D: "edge" (Q4 edge k as legacy (k+1, k)) and "edge_nodes" (as listed); 3D (pressure pilot): "face" (B8 face k in
    ! the PR_B8_FACE order) and "face_nodes" (as listed). Every surface = one legacy edge group, in file order.
    integer, intent(out) :: ns
    integer :: i, j, k, total, ie, e(3), st, el
    character(:), allocatable :: nm
    ns=0; total=0
    do i=1,sets_n
      if(sets(i)%kind==SK_SURFACE) then
        ns=ns+1; total=total+sets(i)%n
      endif
    enddo
    allocate(native_surf(ns), native_loa_edges(5,max(total,1))); native_loa_edges=0
    native_loa_nedge=total
    total=0; k=0
    do i=1,sets_n
      if(sets(i)%kind/=SK_SURFACE) cycle
      k=k+1; nm=trim(sets(i)%name); st=sets(i)%stride
      native_surf(k)%name=nm
      native_surf(k)%first=total+1
      select case(sets(i)%sub)
      case(SS_EDGE,SS_EDGE_NODES)
        if(native_ndim/=2) call native_error("LOA: surface "//nm//": edge surfaces need mesh.dimension 2")
        native_surf(k)%nnode=2; native_surf(k)%index=1
      case(SS_FACE,SS_FACE_NODES)
        if(native_ndim/=3) call native_error("LOA: surface "//nm//": face surfaces need mesh.dimension 3")
        native_surf(k)%nnode=4; native_surf(k)%index=5; native_surf(k)%kind=4
      end select
      do j=1,sets(i)%n
        el=sets(i)%m(st*(j-1)+1)
        total=total+1
        native_loa_edges(5,total)=el
        select case(sets(i)%sub)
        case(SS_EDGE)
          call sets_edge_nodes(el,sets(i)%m(2*j),e(1),e(2))
          native_loa_edges(1:2,total)=e(1:2)
        case(SS_EDGE_NODES)
          native_loa_edges(1:2,total)=sets(i)%m(st*(j-1)+2:st*j)
        case(SS_FACE)
          if(sets_elnn(el)/=8) call native_error("LOA: surface "//nm//": face loads need B8 elements (element "//itoa(el)//")")
          native_loa_edges(1:4,total)=sets_elnodes(PR_B8_FACE(:,sets(i)%m(2*j)),el)
        case(SS_FACE_NODES)
          native_loa_edges(1:4,total)=sets(i)%m(st*(j-1)+2:st*j)
        end select
        if(native_surf(k)%nnode==2) then
          do ie=native_surf(k)%first,total-1
            if((native_loa_edges(1,ie)==native_loa_edges(1,total).and.native_loa_edges(2,ie)==native_loa_edges(2,total)) .or. &
               (native_loa_edges(1,ie)==native_loa_edges(2,total).and.native_loa_edges(2,ie)==native_loa_edges(1,total))) &
              call native_error("LOA: duplicate edge "//itoa(native_loa_edges(1,total))//"-"//itoa(native_loa_edges(2,total))// &
                " in surface "//nm)
          enddo
        endif
      enddo
      native_surf(k)%last=total
    enddo
  end subroutine

  subroutine parse_glb()
    character(:), allocatable :: s, p, q, nm, el
    character(len=TOML_LEN_TEXT) :: ename(NATIVE_MAXG)
    integer :: i, j, k, n, ng, nb, nsec, nf, g, total
    logical :: has(6), act(NATIVE_MAXG), seen_el(NATIVE_MAXG)
    if(integer_value('version')/=1) call native_error('GLB: version must be 1')
    ! [case]
    k=key_index('case.@table')
    call glb_check_keys('case','@table name units description','')
    if(len(string_value('case.name'))==0) call native_error('GLB: empty case.name')
    if(string_value('case.units')/='SI') call native_error('GLB: case.units must be "SI"')
    if(present_key('case.description')) s=string_value('case.description')
    ! [analysis] (G3: static only)
    if(.not.present_key('analysis.type')) call native_error('GLB: [analysis] type = "static" is required with HSTAR_NATIVE_GLB=toml')
    k=key_index('analysis.@table')
    call glb_check_keys('analysis','@table type strength_reduction.@table strength_reduction.factor '// &
      'strength_reduction.amplitude strength_reduction.reset_state_each_increment_from_step','')
    if(present_key('analysis.strength_reduction.@table')) then        ! D3b: legacy kstab / MAT_DE + mat_curve
      k=key_index('analysis.strength_reduction.@table')
      if(.not.(present_key('analysis.strength_reduction.factor') .or. present_key('analysis.strength_reduction.amplitude') &
               .or. present_key('analysis.strength_reduction.reset_state_each_increment_from_step'))) &
        call native_error('GLB: [analysis.strength_reduction] needs factor, amplitude and/or '// &
          'reset_state_each_increment_from_step')
      if(present_key('analysis.strength_reduction.factor')) then
        native_kstab=number('analysis.strength_reduction.factor')
        if(.not.(native_kstab>0 .and. native_kstab<=100)) &
          call native_error('GLB: analysis.strength_reduction.factor out of range (0,100]')
      endif
      sr_amp_on=present_key('analysis.strength_reduction.amplitude')
      ! legacy stab_matde: from this block on, every increment starts from zero Gauss-point state (Fem.f90 stab_initialize)
      if(present_key('analysis.strength_reduction.reset_state_each_increment_from_step')) then
        native_stab_matde=integer_value('analysis.strength_reduction.reset_state_each_increment_from_step')
        if(native_stab_matde<1) call native_error('GLB: analysis.strength_reduction.reset_state_each_increment_from_step '// &
          'must be a step number >= 1')
      endif
    endif
    s=string_value('analysis.type')
    if(s=='transient_heat') then
      native_heat=.true.                  ! temperature pilot phase 2 (T1-T8)
    elseif(s/='static') then
      call native_error('GLB: analysis.type "'//s//'" is deferred (static, transient_heat)')
    endif
    ! [mesh]
    k=key_index('mesh.@table')
    call glb_check_keys('mesh','@table file format dimension sets_file sets_sha256','')
    if(string_value('mesh.file')/='1') call native_error('GLB: mesh.file must be "1" (1.cor / 1.ele next to the TOML)')
    if(string_value('mesh.format')/='hstar-legacy-cor-ele') call native_error('GLB: mesh.format must be "hstar-legacy-cor-ele"')
    native_ndim=integer_value('mesh.dimension')
    if(native_ndim/=2 .and. native_ndim/=3) call native_error('GLB: mesh.dimension must be 2 or 3')
    if(present_key('mesh.sets_file')) call sets_init()
    if(native_sets_on) then
      call glb_elsets_from_sets(ng,ename,total)
    else
    ! [[elset]] (group order = element order in 1.ele, legacy rule)
    ng=doc%count_of('elset')
    if(ng<1 .or. ng>NATIVE_MAXG) call native_error('GLB: need 1..'//itoa(NATIVE_MAXG)//' [[elset]] tables')
    total=0
    do i=1,ng
      p=indexed('elset',i)
      k=key_index(p//'.@table')
      call glb_check_keys(p,'@table name element_count','')
      nm=string_value(p//'.name')
      if(len(nm)==0) call native_error('GLB: empty elset name: '//p)
      if(len(nm)>20) call native_error('GLB: elset name longer than 20 characters: '//nm)
      do j=1,i-1
        if(trim(ename(j))==nm) call native_error('GLB: duplicate elset name: '//nm)
      enddo
      ename(i)=nm
      native_kname(i)=nm
      native_elcounts(i)=integer_value(p//'.element_count')
      if(native_elcounts(i)<1) call native_error('GLB: elset '//nm//': element_count must be >= 1')
      total=total+native_elcounts(i)
    enddo
    endif
    native_materials=ng
    ! element type per elset from [[section]] (G7); MAT checks the rest of the section
    nsec=doc%count_of('section')
    seen_el=.false.
    do k=1,nsec
      p=indexed('section',k)
      el=string_value(p//'.elset')
      g=0
      do j=1,ng
        if(trim(ename(j))==el) g=j
      enddo
      if(g==0) cycle                      ! MAT reports the undefined elset
      s=string_value(p//'.element')
      select case(s)
      case('Q4')
        if(native_ndim/=2) call native_error('GLB: section '//p//': element Q4 needs mesh.dimension 2')
        if(native_heat) call native_error('GLB: section '//p//': transient_heat with Q4 is deferred (only B8)')
        select case(string_value(p//'.formulation'))
        case('plane_strain'); native_gps(g)=.false.
        case('plane_stress'); native_gps(g)=.true.    ! legacy SPTYPE 'PS' (user 2026-10-17)
        case default
          call native_error('GLB: section '//p//': 2D Q4 formulation must be "plane_strain" or "plane_stress"')
        end select
        native_gname(g)='Q4'; native_gindex(g)=5; native_gnnode(g)=4
      case('B8')
        if(native_ndim/=3) call native_error('GLB: section '//p//': element B8 needs mesh.dimension 3')
        if(native_heat) then
          if(string_value(p//'.formulation')/='heat') &
            call native_error('GLB: section '//p//': a transient_heat analysis needs formulation "heat"')
        elseif(string_value(p//'.formulation')/='solid') then
          call native_error('GLB: section '//p//': 3D B8 supports only formulation "solid"')
        endif
        native_gname(g)='B8'; native_gindex(g)=9; native_gnnode(g)=8
      case('T4')                          ! joint pilot D1 (J2): 4-node tetrahedron, legacy index 7 (h4, Elements.f90:412)
        if(native_ndim/=3) call native_error('GLB: section '//p//': element T4 needs mesh.dimension 3')
        if(native_heat) then
          if(string_value(p//'.formulation')/='heat') &
            call native_error('GLB: section '//p//': a transient_heat analysis needs formulation "heat"')
        elseif(string_value(p//'.formulation')/='solid') then
          call native_error('GLB: section '//p//': 3D T4 supports only formulation "solid"')
        endif
        native_gname(g)='T4'; native_gindex(g)=7; native_gnnode(g)=4
      case('joint')                       ! joint pilot D3a: GOODMAN node-pair rows, legacy index 5 (2D Q4 row) / 9 (3D B8 row)
        if(native_heat) call native_error('GLB: section '//p//': element joint does not apply to transient_heat')
        if(string_value(p//'.formulation')/='goodman') &
          call native_error('GLB: section '//p//': element joint supports only formulation "goodman"')
        native_gjoint(g)=.true.
        if(native_ndim==2) then
          native_gname(g)='Q4'; native_gindex(g)=5; native_gnnode(g)=4
        else
          native_gname(g)='B8'; native_gindex(g)=9; native_gnnode(g)=8
        endif
      case('L2')                          ! matline L1-L4: 2-node bar (legacy index 1), 2D or 3D
        if(native_heat) call native_error('GLB: section '//p//': transient_heat with L2 is deferred (only B8)')
        if(present_key(p//'.formulation')) call native_error('GLB: section '//p// &
          ': formulation does not apply to element L2 (2-node bar)')
        native_gname(g)='L2'; native_gindex(g)=1; native_gnnode(g)=2
      case('B2')                          ! beam pilot C2: 2-node beam (legacy index 20, b2), 2D or 3D
        if(native_heat) call native_error('GLB: section '//p//': transient_heat with B2 is deferred')
        if(present_key(p//'.formulation')) call native_error('GLB: section '//p// &
          ': formulation does not apply to element B2 (2-node beam)')
        native_gname(g)='L2'; native_gindex(g)=20; native_gnnode(g)=2
        native_beam=.true.
        call glb_beam_orientation(p,g)
      case('bond_slip')                   ! steel25 pilot: legacy index 25 "steel", 2-node bond-slip spring (translations)
        if(native_heat) call native_error('GLB: section '//p//': transient_heat with bond_slip is deferred')
        if(present_key(p//'.formulation')) call native_error('GLB: section '//p// &
          ': formulation does not apply to element bond_slip (2-node bond-slip spring)')
        native_gname(g)='L2'; native_gindex(g)=25; native_gnnode(g)=2
        native_bond=.true.
      case('steel')
        call native_error('GLB: section '//p//': element "steel" is spelled "bond_slip" in TOML (legacy index 25)')
      case default
        call native_error('GLB: section '//p//': element "'//s//'" is deferred in the GLB pilot '// &
          '(only Q4 / B8 / T4 / L2 / B2 / joint / bond_slip)')
      end select
      seen_el(g)=.true.
    enddo
    call parse_bond_slip(ng,ename)
    ! [[step]] skeleton: active element sets, reset, increments (load/boundary/controls data: LOA/PRE/MAN)
    nb=doc%count_of('step')
    if(nb<1 .or. nb>NATIVE_MAXB) call native_error('GLB: need 1..'//itoa(NATIVE_MAXB)//' [[step]] tables')
    native_blocks=nb
    do i=1,nb
      p=indexed('step',i)
      k=key_index(p//'.@table')
      if(.not.native_heat) k=key_index(p//'.controls.@table')
      if(len(string_value(p//'.name'))==0) call native_error('GLB: empty step name: '//p)
      s=string_value(p//'.procedure')
      if(native_heat) then
        if(s/='transient_heat') call native_error('GLB: '//p//'.procedure must be "transient_heat" (analysis.type)')
      elseif(s/='static') then
        call native_error('GLB: '//p//'.procedure "'//s//'" is deferred (only static)')
      endif
      if(string_value(p//'.load_mode')/='load') call native_error('GLB: '//p//'.load_mode must be "load"')
      k=key_index(p//'.reset_state')
      if(doc%entry(k)%kind/=TV_BOOL) call native_error('expected boolean: '//p//'.reset_state')
      native_reset(i)=merge(1,0,doc%entry(k)%lvalue)
      if(native_heat .and. native_reset(i)/=0) call native_error('GLB: '//p//'.reset_state must be false in a transient_heat analysis')
      if(.not.native_heat) then
      if(integer_value(p//'.controls.increments')/=1) call native_error('GLB: '//p//'.controls.increments must be 1')
      s=string_value(p//'.controls.stiffness_update')
      if(s/='first_iteration' .and. s/='every_iteration') &
        call native_error('GLB: '//p//'.controls.stiffness_update must be "first_iteration" or "every_iteration"')
      k=merge(4,5,s=='every_iteration')             ! D3b: legacy type_nl 5 / 4 (full Newton), one value for all steps
      if(i>1 .and. k/=native_type_nl) call native_error('GLB: '//p//'.controls.stiffness_update must be the same in '// &
        'every step (legacy type_nl is global)')
      native_type_nl=k
      endif
      n=integer_value(p//'.active_elsets.count')
      if(n<1) call native_error('GLB: '//p//': active_elsets is empty')
      act=.false.
      do j=1,n
        q=string_value(indexed(p//'.active_elsets',j))
        g=0
        do k=1,ng
          if(trim(ename(k))==q) g=k
        enddo
        if(g==0) call native_error('GLB: '//p//': undefined elset '//q//' in active_elsets')
        if(native_sets_on) k=sets_find(SK_ELSET,q)
        if(act(g)) call native_error('GLB: '//p//': elset '//q//' listed twice in active_elsets')
        act(g)=.true.
      enddo
      native_appearance(1:ng,i)=merge(1,0,act(1:ng))
    enddo
    ! [solver] (G5: .sol is implied by linear = "profile")
    k=key_index('solver.@table')
    if(native_sol_toml) then
      call glb_check_keys('solver','@table linear pardiso.@table pardiso.matrix_type pardiso.threads pardiso.iparm','')
    else
      call glb_check_keys('solver','@table linear','')
    endif
    s=string_value('solver.linear')
    if(native_heat) then
      ! T6 option A: PARDISO, no environment variables; .sol (mtype 11 nonsymmetric) is implied
      if(s/='pardiso') call native_error('GLB: a transient_heat analysis needs solver.linear "pardiso" '// &
        '(nonsymmetric pipe advection; .sol mtype 11 is implied)')
    elseif(s=='pardiso' .and. native_sol_toml) then
      native_pardiso=.true.                 ! sol-pilot V1: static PARDISO from [solver.pardiso]
    elseif(s/='profile') then
      call native_error('GLB: solver.linear "'//s//'" is deferred for static (only profile; .sol 0 0 1 1 is implied; '// &
        'pardiso needs HSTAR_NATIVE_SOL=toml)')
    endif
    if(native_sol_toml) call parse_sol(s)
    ! [output] (G6)
    k=key_index('output.@table')
    if(native_heat) then
      call glb_check_keys('output','@table format','field')
      if(string_value('output.format')/='gid') call native_error('GLB: output.format must be "gid"')
      if(integer_value('output.field.count')/=1) call native_error('GLB: transient_heat output.field must be ["T"]')
      if(string_value('output.field[1]')/='T') call native_error('GLB: transient_heat output.field must be ["T"]')
      native_gid_f=0
    else
    call glb_check_keys('output','@table format stress_averaging','field stress_averaging_exclude')
    if(string_value('output.format')/='gid') call native_error('GLB: output.format must be "gid"')
    if(string_value('output.stress_averaging')/='direct') call native_error('GLB: output.stress_averaging must be "direct"')
    ! element sets excluded from nodal stress averaging (legacy .glb average_appear(group) = 0, Output.f90:5144-5149;
    ! legacy -1/-2 with 0 < elcod_local < 0.03 is the same exclusion). Line (L2) and joint sets are never averaged.
    if(present_key('output.stress_averaging_exclude.count')) then
      nf=integer_value('output.stress_averaging_exclude.count')
      if(nf<1) call native_error('GLB: output.stress_averaging_exclude must list at least one element set')
      do j=1,nf
        s=string_value(indexed('output.stress_averaging_exclude',j))
        k=0
        do g=1,ng
          if(trim(native_kname(g))==s) k=g
        enddo
        if(k==0) call native_error('GLB: output.stress_averaging_exclude: undefined element set "'//s//'"')
        if(native_avg_off(k)) call native_error('GLB: output.stress_averaging_exclude lists "'//s//'" twice')
        if(native_gnnode(k)==2 .or. native_gjoint(k)) call native_error('GLB: output.stress_averaging_exclude: "'//s// &
          '" is a line / joint element set, which is never nodally averaged')
        native_avg_off(k)=.true.
      enddo
    endif
    nf=integer_value('output.field.count')
    has=.false.
    do j=1,nf
      s=string_value(indexed('output.field',j))
      select case(s)
      case('u'); k=1
      case('s'); k=2
      case('f'); k=3
      case('ms'); k=4                  ! D3b: legacy gid_ms / gid_ep / gid_Y (mc_joint analyses only)
      case('ep'); k=5
      case('yield'); k=6
      case default; call native_error('GLB: output.field "'//s//'" is deferred (only u, s, f; ms, ep, yield with mc_joint)')
      end select
      if(has(k)) call native_error('GLB: output.field "'//s//'" listed twice')
      has(k)=.true.
    enddo
    if(.not.(has(1) .and. has(2))) call native_error('GLB: output.field must contain "u" and "s"')
    native_gid_f=merge(1,0,has(3))
    native_gid_ms=merge(1,0,has(4)); native_gid_ep=merge(1,0,has(5)); native_gid_y=merge(1,0,has(6))
    endif
    ! derived mesh counts (G4); full checks in validate_mesh_general
    native_points=glb_count_rows(trim(native_mesh_prefix)//'.cor')
    native_elements=glb_count_rows(trim(native_mesh_prefix)//'.ele')
    if(native_points<1) call native_error('GLB: 1.cor has no rows')
    ! orientation node ids: checked here, after 1.cor is counted (sections are parsed before the mesh counts)
    if(any(native_pdir(:,1:ng)>native_points)) call native_error('GLB: orientation.node_pair node out of range 1..'// &
      itoa(native_points))
    if(total/=native_elements) call native_error('GLB: sum of elset element_count ('//itoa(total)// &
      ') does not match the '//itoa(native_elements)//' element rows of 1.ele')
    do g=1,ng
      if(.not.seen_el(g)) native_gnnode(g)=4   ! MAT reports the missing section
    enddo
  end subroutine

  subroutine native_nrt_deferred(what)
    character(*), intent(in) :: what
    call native_error("NRT: "//what//" is deferred (temperature pilot: only TRAL ties from 1.set)")
  end subroutine

  ! steel25 pilot: top-level [crack] (ftcrack) and [bond_slip] (the legacy .glb global bond line). Values are passed
  ! through unchanged; only the representation is checked (DESIGN.md steel25-pilot section 5).
  subroutine parse_bond_slip(ng,ename)
    integer, intent(in) :: ng
    character(len=TOML_LEN_TEXT), intent(in) :: ename(:)
    character(:), allocatable :: s, q
    integer :: k, j, n, g, m, first, e
    real(real64) :: st, r
    logical :: has, taken(NATIVE_MAXG)
    if(present_key('crack.@table')) then
      k=key_index('crack.@table')
      call glb_check_keys('crack','@table tensile_strength','')
      if(.not.present_key('crack.tensile_strength')) call native_error('GLB: [crack] needs tensile_strength (legacy ftcrack)')
      native_ftcrack=number('crack.tensile_strength')
      if(.not.(native_ftcrack>0)) call native_error('GLB: crack.tensile_strength must be > 0')
    endif
    has=present_key('bond_slip.@table')
    if(native_bond .and. .not.has) call native_error('GLB: element "bond_slip" needs a top-level [bond_slip] table '// &
      '(legacy .glb bond line: law, unit_factor, normal)')
    if(has .and. .not.native_bond) call native_error('GLB: [bond_slip] is given but no [[section]] has element "bond_slip"')
    if(.not.has) return
    k=key_index('bond_slip.@table')
    call glb_check_keys('bond_slip','@table law stiffness unit_factor normal', &
      'normal_stiffness local_bar_sets axial_constraint_elsets')
    if(.not.present_key('bond_slip.unit_factor')) call native_error('GLB: [bond_slip] needs unit_factor (legacy coefMpa: '// &
      'model stress units per MPa, 1.0e6 for Pa, 1.0e3 for kPa)')
    native_coefmpa=number('bond_slip.unit_factor')
    if(.not.(native_coefmpa>0)) call native_error('GLB: bond_slip.unit_factor must be > 0')
    if(.not.present_key('bond_slip.law')) call native_error('GLB: [bond_slip] needs law (1..5 or "linear")')
    k=key_index('bond_slip.law')
    select case(doc%entry(k)%kind)
    case(TV_INT)
      native_ikindks=doc%entry(k)%ivalue
      if(native_ikindks<1 .or. native_ikindks>5) call native_error('GLB: bond_slip.law must be 1..5 (legacy ikindks) '// &
        'or "linear"')
      if(present_key('bond_slip.stiffness')) call native_error('GLB: bond_slip.stiffness applies only to law = "linear"')
    case(TV_STR)
      if(doc%entry(k)%svalue/='linear') call native_error('GLB: bond_slip.law must be 1..5 (legacy ikindks) or "linear"')
      if(.not.present_key('bond_slip.stiffness')) call native_error('GLB: law = "linear" needs bond_slip.stiffness '// &
        '(bond stress per unit slip, model stress units per m)')
      st=number('bond_slip.stiffness')
      if(.not.(st>0)) call native_error('GLB: bond_slip.stiffness must be > 0')
      ! legacy ikindks<0: ks=10**(-ikindks), ikindks>99: ks=ikindks, in MPa/mm; the solver scales by coefMpa*1000
      r=st/(native_coefmpa*1000._real64)
      native_ikindks=0
      n=nint(log10(r))
      if(n>=1 .and. n<=99) then
        if(abs(10._real64**n-r)<=1.e-12_real64*r) native_ikindks=-n
      endif
      if(native_ikindks==0 .and. r>99.5_real64 .and. r<2.e9_real64) then
        if(abs(anint(r)-r)<=1.e-9_real64*r) native_ikindks=nint(r)
      endif
      if(native_ikindks==0) call native_error('GLB: bond_slip.stiffness / (unit_factor*1000) must be 10^n (n = 1..99) '// &
        'or an integer > 99 (the legacy ikindks code of a linear bond law)')
    case default
      call native_error('GLB: bond_slip.law must be an integer 1..5 or "linear"')
    end select
    s=string_value('bond_slip.normal')
    select case(s)
    case('spring')
      native_doubsig=2
      if(.not.present_key('bond_slip.normal_stiffness.count')) call native_error('GLB: normal = "spring" needs '// &
        'bond_slip.normal_stiffness (legacy ktan1 [, ktan2])')
      n=integer_value('bond_slip.normal_stiffness.count')
      if(n/=native_ndim-1) call native_error('GLB: bond_slip.normal_stiffness needs '//itoa(native_ndim-1)// &
        ' value(s) for mesh.dimension '//itoa(native_ndim))
      native_ktan1=number('bond_slip.normal_stiffness[1]')
      native_ktan2=native_ktan1
      if(n==2) native_ktan2=number('bond_slip.normal_stiffness[2]')
      if(.not.(native_ktan1>0 .and. native_ktan2>0)) call native_error('GLB: bond_slip.normal_stiffness must be > 0')
      if(present_key('bond_slip.axial_constraint_elsets.count')) call native_error('GLB: bond_slip.axial_constraint_elsets '// &
        'applies only to normal = "constraint" (legacy lelenrt with doubsig 1)')
    case('constraint')
      native_doubsig=1
      if(present_key('bond_slip.normal_stiffness.count')) call native_error('GLB: bond_slip.normal_stiffness applies only '// &
        'to normal = "spring" (legacy doubsig 1 does not use ktan1/ktan2)')
    case default
      call native_error('GLB: bond_slip.normal must be "spring" (legacy doubsig 2) or "constraint" (doubsig 1)')
    end select
    native_nlocalbeam=0; taken=.false.
    if(present_key('bond_slip.local_bar_sets.count')) then
      n=integer_value('bond_slip.local_bar_sets.count')
      if(n<1) call native_error('GLB: bond_slip.local_bar_sets is empty')
      do j=1,n
        q=string_value(indexed('bond_slip.local_bar_sets',j))
        g=bond_elset(q,ng,ename,'local_bar_sets')
        if(native_gindex(g)/=1 .and. native_gindex(g)/=20) call native_error('GLB: bond_slip.local_bar_sets: "'//q// &
          '" is not an L2 / B2 element set')
        if(taken(g)) call native_error('GLB: bond_slip.local_bar_sets lists "'//q//'" twice')
        taken(g)=.true.
        native_nlocalbeam=native_nlocalbeam+1
        native_listglocbeam(native_nlocalbeam)=g
      enddo
    endif
    native_ndimnrt=0; taken=.false.
    if(present_key('bond_slip.axial_constraint_elsets.count')) then
      n=integer_value('bond_slip.axial_constraint_elsets.count')
      if(n<1) call native_error('GLB: bond_slip.axial_constraint_elsets is empty')
      m=0
      do j=1,n
        q=string_value(indexed('bond_slip.axial_constraint_elsets',j))
        g=bond_elset(q,ng,ename,'axial_constraint_elsets')
        if(native_gindex(g)/=25) call native_error('GLB: bond_slip.axial_constraint_elsets: "'//q// &
          '" is not a bond_slip element set')
        if(taken(g)) call native_error('GLB: bond_slip.axial_constraint_elsets lists "'//q//'" twice')
        taken(g)=.true.
        m=m+native_elcounts(g)
      enddo
      allocate(native_lelenrt(m))
      do j=1,n
        q=string_value(indexed('bond_slip.axial_constraint_elsets',j))
        g=bond_elset(q,ng,ename,'axial_constraint_elsets')
        first=sum(native_elcounts(1:g-1))           ! elsets follow the element order of 1.ele
        do e=1,native_elcounts(g)
          native_ndimnrt=native_ndimnrt+1
          native_lelenrt(native_ndimnrt)=first+e
        enddo
      enddo
    endif
  end subroutine

  integer function bond_elset(q,ng,ename,what) result(g)
    character(*), intent(in) :: q, what
    integer, intent(in) :: ng
    character(len=TOML_LEN_TEXT), intent(in) :: ename(:)
    integer :: j
    g=0
    do j=1,ng
      if(trim(ename(j))==q) g=j
    enddo
    if(g==0) call native_error('GLB: bond_slip.'//what//': undefined element set "'//q//'"')
  end function

  subroutine native_glb_deferred(what)
    character(*), intent(in) :: what
    call native_error("GLB: "//what//" is deferred (GLB pilot) and cannot be driven by TOML")
  end subroutine

  ! General Q4 / B8 mesh check for HSTAR_NATIVE_GLB=toml (DESIGN.md section 6 rules 2-3).
  subroutine validate_mesh_general()
    integer :: u, ios, k, id, g, last, nn, nd, j, mat, nt
    integer :: nodes(8)
    integer :: dn(8), ndn
    logical :: degen
    integer :: nmixed
    real(real64) :: jmax, jmin, xr(3,8)
    logical :: ok, frontrep
    character(4096) :: mixed_rows
    real(real64), allocatable :: xyz(:,:)
    logical, allocatable :: used(:)
    real(real64) :: c(3)
    character(1024) :: row
    logical, allocatable :: barnode(:)
    integer, allocatable :: bondrow(:,:)
    integer :: nbond
    nmixed=0; mixed_rows=''
    nd=native_ndim
    allocate(xyz(3,native_points), used(native_points)); xyz=0; used=.false.
    allocate(barnode(native_points), bondrow(3,native_elements)); barnode=.false.; nbond=0
    if(native_heat) then
      allocate(mesh_ele(8,native_elements)); mesh_ele=0
    endif
    open(newunit=u,file=trim(native_mesh_prefix)//'.cor',status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh 1.cor')
    do k=1,native_points
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) call native_error('missing mesh coordinate row')
      if(token_count(row) /= nd+1) call native_error('GLB: 1.cor row '//itoa(k)//' must have '//itoa(nd+1)// &
        ' columns (id + mesh.dimension coordinates)')
      c=0
      read(row,*,iostat=ios) id,c(1:nd)
      if(ios /= 0) call native_error('invalid mesh coordinates')
      if(id /= k) call native_error('GLB: 1.cor node numbers must be 1..n in order (row '//itoa(k)//')')
      if(.not.all(ieee_is_finite(c))) call native_error('non-finite mesh')
      xyz(:,k)=c
    enddo
    close(u)
    open(newunit=u,file=trim(native_mesh_prefix)//'.ele',status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh 1.ele')
    g=1; last=native_elcounts(1)
    do k=1,native_elements
      do while(k>last)
        g=g+1; last=last+native_elcounts(g)
      enddo
      nn=native_gnnode(g)
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) call native_error('missing mesh element row')
      ! legacy read_element reads "id nodes(1:nnode)"; a trailing group column is optional
      nt=token_count(row)
      if(nt /= nn+1 .and. nt /= nn+2) call native_error('GLB: 1.ele row '//itoa(k)//' must have '//itoa(nn+1)// &
        ' or '//itoa(nn+2)//' columns for element '//trim(native_gname(g))//' (elset '//trim(native_kname(g))//')')
      mat=g
      if(nt==nn+2) then
        read(row,*,iostat=ios) id,nodes(1:nn),mat
      else
        read(row,*,iostat=ios) id,nodes(1:nn)
      endif
      if(ios /= 0) call native_error('invalid mesh connectivity')
      if(id /= k) call native_error('GLB: 1.ele element numbers must be 1..n in order (row '//itoa(k)//')')
      if(mat /= g) call native_error('GLB: 1.ele row '//itoa(k)//': group column '//itoa(mat)//' does not match elset '// &
        itoa(g)//' (elsets follow the element order of 1.ele)')
      if(any(nodes(1:nn)<1) .or. any(nodes(1:nn)>native_points)) &
        call native_error('GLB: 1.ele row '//itoa(k)//': node out of range 1..'//itoa(native_points))
      if(native_gindex(g)==1 .or. native_gindex(g)==20) barnode(nodes(1:2))=.true.
      if(native_gindex(g)==25) then
        nbond=nbond+1; bondrow(1:2,nbond)=nodes(1:2); bondrow(3,nbond)=k
      endif
      degen=.false.
      ! user 2026-10-17: a B8 solid row may collapse in any node pattern (as legacy); the 2x2x2 Gauss-point Jacobian
      ! check below decides (a collapse that leaves any Gauss point with det J <= 0 is rejected)
      if(nn==8 .and. native_gindex(g)==9 .and. .not.native_gjoint(g)) then
        do j=2,8
          if(any(nodes(1:j-1)==nodes(j))) degen=.true.
        enddo
      endif
      ! D3b (train05): Q4 triangle collapse 3=4 (legacy Q4 row with a repeated last node); nodes 1-3 distinct
      if(nn==4 .and. native_gindex(g)==5 .and. .not.native_gjoint(g)) then
        degen=nodes(3)==nodes(4) .and. nodes(1)/=nodes(2) .and. nodes(1)/=nodes(3) .and. nodes(2)/=nodes(3)
      endif
      ! user 2026-10-17 (thin layer): with layer_thickness legacy never integrates the stored back nodes
      ! (Elements.f90 normal_local_inc rebuilds them), so repeated nodes are judged on the rebuilt element below
      if(native_elcod(g)>0 .and. .not.native_gjoint(g)) degen=.true.
      if(.not.degen) then
      do j=2,nn
        if(any(nodes(1:j-1)==nodes(j))) call native_error('GLB: 1.ele row '//itoa(k)//': repeated node'// &
          ' (only the Q4 triangle 3=4 collapse and B8 solid collapses are accepted)')
      enddo
      endif
      used(nodes(1:nn))=.true.
      if(allocated(mesh_ele)) mesh_ele(1:nn,k)=nodes(1:nn)
      if(native_elcod(g)>0 .and. .not.native_gjoint(g)) then
        ! D3b / user 2026-10-17: thin layer (section.layer_thickness = legacy elcod_local, any Q4/B8 solid): legacy
        ! rebuilds the back nodes as front + mean unit normal * layer_thickness (Elements.f90 normal_local_inc) and never
        ! integrates the stored ones. 2D: the rebuilt Q4 (1, 2, 2+n, 1+n, n = left normal of 1->2) is valid iff edge 1-2
        ! has length. 3D: the rebuilt B8 is checked at its 2x2x2 Gauss points (front face with repeated nodes: mixed
        ! sign read with the warning, all non-positive rejected; otherwise any non-positive rejected).
        if(native_ndim==2) then
          if(sum((xyz(1:2,nodes(2))-xyz(1:2,nodes(1)))**2)<=0) &
            call native_error('GLB: 1.ele row '//itoa(k)//': thin layer front edge (nodes 1-2) has zero length')
        else
          if(q4_face_area(xyz(:,nodes(1:4)))<=0) &
            call native_error('GLB: 1.ele row '//itoa(k)//': thin layer front face (nodes 1-4) has zero area')
          call thin_layer_rebuild(xyz(:,nodes(1:4)),native_elcod(g),xr,ok)
          if(.not.ok) call native_error('GLB: 1.ele row '//itoa(k)// &
            ': thin layer front face (nodes 1-4) has no normal at a Gauss point')
          jmin=b8_gauss_min_jacobian(xr,jmax)
          frontrep=nodes(1)==nodes(2) .or. nodes(1)==nodes(3) .or. nodes(1)==nodes(4) .or. nodes(2)==nodes(3) .or. &
                   nodes(2)==nodes(4) .or. nodes(3)==nodes(4)
          if(jmax<=0 .or. (jmin<=0 .and. .not.frontrep)) call native_error('GLB: 1.ele row '//itoa(k)// &
            ': thin layer (rebuilt back nodes) B8 Jacobian is not positive at the Gauss points')
          if(jmin<=0) then
            nmixed=nmixed+1
            if(len_trim(mixed_rows)<len(mixed_rows)-12) mixed_rows=trim(mixed_rows)//' '//itoa(k)
          endif
        endif
      elseif(native_gjoint(g)) then
        ! D3a: a joint row is two faces (node pairs 2D 1-4, 2-3; 3D i-(i+4)); the faces may coincide (zero thickness),
        ! so no volume Jacobian. The front face (integration face, Stiff.f90:123-131) must be non-degenerate.
        if(native_ndim==2) then
          if(sum((xyz(1:2,nodes(2))-xyz(1:2,nodes(1)))**2)<=0) &
            call native_error('GLB: 1.ele row '//itoa(k)//': joint front edge (nodes 1-2) has zero length')
        elseif(q4_face_area(xyz(:,nodes(1:4)))<=0) then
          call native_error('GLB: 1.ele row '//itoa(k)//': joint front face (nodes 1-4) has zero area')
        endif
      elseif(nn==2 .and. native_gindex(g)==25) then
        ! steel25 pilot: a bond-slip spring joins a bar node to a host node; zero length (coincident nodes) is normal
        if(nodes(1)==nodes(2)) call native_error('GLB: 1.ele row '//itoa(k)//': bond_slip element joins a node to itself')
      elseif(nn==2) then
        if(sum((xyz(:,nodes(2))-xyz(:,nodes(1)))**2)<=0) &
          call native_error('GLB: 1.ele row '//itoa(k)//': L2 element has zero length')
      elseif(native_gindex(g)==7) then
        if(t4_jacobian(xyz(:,nodes(1:4)))<=0) &
          call native_error('GLB: 1.ele row '//itoa(k)//': T4 Jacobian is not positive (node order / distortion)')
      elseif(nn==4 .and. degen) then
        ! collapsed Q4: the bilinear Jacobian is det = 2*area(1,2,3)*(1+eta)/4 > 0 at every Gauss point iff area > 0
        if((xyz(1,nodes(2))-xyz(1,nodes(1)))*(xyz(2,nodes(3))-xyz(2,nodes(1))) - &
           (xyz(2,nodes(2))-xyz(2,nodes(1)))*(xyz(1,nodes(3))-xyz(1,nodes(1)))<=0) &
          call native_error('GLB: 1.ele row '//itoa(k)//': collapsed Q4 (3=4) triangle 1-2-3 is not counter-clockwise')
      elseif(nn==4) then
        if(q4_min_jacobian(xyz(1:2,nodes(1:4)))<=0) &
          call native_error('GLB: 1.ele row '//itoa(k)//': Q4 Jacobian is not positive (node order / distortion)')
      elseif(degen) then
        ! user 2026-10-17: a collapsed B8 is rejected only if det J <= 0 at all 8 Gauss points; mixed sign is read
        ! (as legacy) with a warning listing the rows
        if(b8_gauss_min_jacobian(xyz(:,nodes(1:8)),jmax)<=0) then
          if(jmax<=0) call native_error('GLB: 1.ele row '//itoa(k)// &
            ': degenerate B8 Jacobian is not positive at any Gauss point')
          nmixed=nmixed+1
          if(len_trim(mixed_rows)<len(mixed_rows)-12) mixed_rows=trim(mixed_rows)//' '//itoa(k)
        endif
      else
        ! user 2026-10-17: regular B8 checked at the 2x2x2 Gauss points, where legacy evaluates the Jacobian (was corners)
        if(b8_gauss_min_jacobian(xyz(:,nodes(1:8)))<=0) &
          call native_error('GLB: 1.ele row '//itoa(k)//': B8 Jacobian is not positive at the Gauss points (node order / distortion)')
      endif
    enddo
    if(nmixed>0) write(*,'(a)') 'WARNING native mesh: '//itoa(nmixed)//' degenerate B8 row(s) with mixed-sign'// &
      ' Jacobian at the Gauss points (read as legacy): 1.ele rows'//trim(mixed_rows)
    read(u,'(a)',iostat=ios) row
    if(ios >= 0) call native_error('extra mesh element row')
    close(u)
    ! steel25 pilot: legacy steel_spring_parameter / Stiff.f90 index 25 need exactly one node of each spring on a
    ! bar (L2 / B2) element (icpspring), otherwise the solver stops ('stop here!')
    do j=1,nbond
      if(count(barnode(bondrow(1:2,j)))/=1) call native_error('GLB: 1.ele row '//itoa(bondrow(3,j))// &
        ': bond_slip element needs exactly one of its 2 nodes on an L2 / B2 bar element (legacy icpspring)')
    enddo
    do k=1,native_points
      if(.not.used(k)) call native_error('GLB: node '//itoa(k)//' is not used by any element')
    enddo
    if(native_heat) call heat_check_mesh_refs()
  end subroutine

  ! face4 rows must be a face of their B8 element; pipe segments must be an element edge of their elset.
  subroutine heat_check_mesh_refs()
    integer, parameter :: fc(4,6)=reshape([1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8],[4,6])
    integer, parameter :: ft(4,4)=reshape([1,2,3,3, 1,2,4,4, 2,3,4,4, 3,1,4,4],[4,4])
    integer :: t, e, k, j, ip, ie, g, first, last
    logical :: ok
    do t=1,native_nface
      e=native_faces(5,t); ok=.false.
      g=1; last=native_elcounts(1)
      do while(e>last)
        g=g+1; last=last+native_elcounts(g)
      enddo
      if(native_gindex(g)==7) then           ! D1: T4 face = degenerate 4-node row (yl_sets T4_FACE)
        do k=1,4
          if(all(native_faces(1:4,t)==mesh_ele(ft(:,k),e))) ok=.true.
        enddo
      elseif(native_gindex(g)==9) then
      do k=1,6
        if(all([(any(native_faces(1:4,t)==mesh_ele(fc(j,k),e)),j=1,4)])) ok=.true.
      enddo
      endif
      if(.not.ok) call native_error('TEM: face '//itoa(t)//' (nodes '//itoa(native_faces(1,t))//' '//itoa(native_faces(2,t))// &
        ' '//itoa(native_faces(3,t))//' '//itoa(native_faces(4,t))//') is not a face of element '//itoa(e))
    enddo
    if(.not.allocated(native_pipe)) return
    do ip=1,size(native_pipe)
      g=native_pipe(ip)%igroup
      first=sum(native_elcounts(1:g-1))+1; last=first+native_elcounts(g)-1
      do ie=1,size(native_pipe(ip)%seg,2)
        ok=.false.
        do e=first,last
          if(any(mesh_ele(1:8,e)==native_pipe(ip)%seg(1,ie)) .and. any(mesh_ele(1:8,e)==native_pipe(ip)%seg(2,ie))) then
            ok=.true.; exit
          endif
        enddo
        if(.not.ok) call native_error('TEM: pipe '//itoa(ip)//' segment '//itoa(ie)//': nodes do not share an element of elset '// &
          trim(native_kname(g)))
      enddo
    enddo
  end subroutine

  ! Smallest Jacobian determinant of a trilinear hexahedron (nodes 1-4 bottom counter-clockwise seen
  ! from the top face, 5-8 above them), sampled at the 8 corners and the centre.
  real(real64) function b8_min_jacobian(x) result(dmin)
    real(real64), intent(in) :: x(3,8)
    real(real64), parameter :: xi(8)=[-1,1,1,-1,-1,1,1,-1], et(8)=[-1,-1,1,1,-1,-1,1,1], ze(8)=[-1,-1,-1,-1,1,1,1,1]
    real(real64) :: p(3), j(3,3), dn(3)
    integer :: ip, a
    dmin=huge(1._real64)
    do ip=0,8
      if(ip==0) then
        p=0
      else
        p=[xi(ip),et(ip),ze(ip)]
      endif
      j=0
      do a=1,8
        dn(1)=xi(a)*(1+et(a)*p(2))*(1+ze(a)*p(3))/8
        dn(2)=et(a)*(1+xi(a)*p(1))*(1+ze(a)*p(3))/8
        dn(3)=ze(a)*(1+xi(a)*p(1))*(1+et(a)*p(2))/8
        j(:,1)=j(:,1)+x(:,a)*dn(1); j(:,2)=j(:,2)+x(:,a)*dn(2); j(:,3)=j(:,3)+x(:,a)*dn(3)
      enddo
      dmin=min(dmin, j(1,1)*(j(2,2)*j(3,3)-j(3,2)*j(2,3)) - j(1,2)*(j(2,1)*j(3,3)-j(3,1)*j(2,3)) &
                     + j(1,3)*(j(2,1)*j(3,2)-j(3,1)*j(2,2)))
    enddo
  end function

  ! D3a: area of a 3D quadrilateral face (two triangles 1-2-3, 1-3-4)
  ! Legacy thin-layer geometry (Elements.f90 normal_local_inc, 3D): unit normals of the bilinear front face at its 2x2
  ! Gauss points are averaged and normalised; back node i+4 = front node i + normal * h. ok=.false. if a normal vanishes.
  subroutine thin_layer_rebuild(f,h,x,ok)
    real(real64), intent(in) :: f(3,4), h
    real(real64), intent(out) :: x(3,8)
    logical, intent(out) :: ok
    real(real64), parameter :: xi(4)=[-1,1,1,-1], et(4)=[-1,-1,1,1]
    real(real64) :: g, s(3), t(3), n(3), r(3), a
    integer :: ip, i
    g=1._real64/sqrt(3._real64); r=0; ok=.true.
    do ip=1,4
      s=0; t=0
      do i=1,4
        s=s+f(:,i)*xi(i)*(1+et(i)*g*et(ip))/4
        t=t+f(:,i)*et(i)*(1+xi(i)*g*xi(ip))/4
      enddo
      n=[s(2)*t(3)-s(3)*t(2), s(3)*t(1)-s(1)*t(3), s(1)*t(2)-s(2)*t(1)]
      a=sqrt(sum(n**2))
      if(a<=0) then; ok=.false.; return; endif
      r=r+n/a
    enddo
    a=sqrt(sum(r**2))
    if(a<=0) then; ok=.false.; return; endif
    r=r/a
    x(:,1:4)=f
    do i=1,4
      x(:,i+4)=f(:,i)+r*h
    enddo
  end subroutine

  real(real64) function q4_face_area(x)
    real(real64), intent(in) :: x(3,4)
    real(real64) :: a(3), b(3), c(3)
    a=x(:,2)-x(:,1); b=x(:,3)-x(:,1); c=x(:,4)-x(:,1)
    q4_face_area=0.5_real64*(norm2([a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]) + &
                             norm2([b(2)*c(3)-b(3)*c(2), b(3)*c(1)-b(1)*c(3), b(1)*c(2)-b(2)*c(1)]))
  end function

  ! T4 (Elements.f90:2767 shape 1-s-t-u, s, t, u): constant Jacobian det[x2-x1; x3-x1; x4-x1]
  real(real64) function t4_jacobian(x) result(d)
    real(real64), intent(in) :: x(3,4)
    real(real64) :: a(3), b(3), c(3)
    a=x(:,2)-x(:,1); b=x(:,3)-x(:,1); c=x(:,4)-x(:,1)
    d=a(1)*(b(2)*c(3)-b(3)*c(2))-a(2)*(b(1)*c(3)-b(3)*c(1))+a(3)*(b(1)*c(2)-b(2)*c(1))
  end function

  ! D1: degenerate B8 (any collapse pattern, e.g. wedge 3=4,7=8 / pyramid 5=6=7=8): legacy integrates it with the
  ! 2x2x2 B8 rule, so its Jacobian must be positive at those Gauss points (it is 0 at the collapsed corners).
  real(real64) function b8_gauss_min_jacobian(x, dmax) result(dmin)
    real(real64), intent(in) :: x(3,8)
    real(real64), intent(out), optional :: dmax
    real(real64), parameter :: xi(8)=[-1,1,1,-1,-1,1,1,-1], et(8)=[-1,-1,1,1,-1,-1,1,1], ze(8)=[-1,-1,-1,-1,1,1,1,1]
    real(real64) :: p(3), j(3,3), dn(3), g, g8
    integer :: ip, a
    g=1._real64/sqrt(3._real64)
    dmin=huge(1._real64)
    if(present(dmax)) dmax=-huge(1._real64)
    do ip=1,8
      p=g*[xi(ip),et(ip),ze(ip)]
      j=0
      do a=1,8
        dn(1)=xi(a)*(1+et(a)*p(2))*(1+ze(a)*p(3))/8
        dn(2)=et(a)*(1+xi(a)*p(1))*(1+ze(a)*p(3))/8
        dn(3)=ze(a)*(1+xi(a)*p(1))*(1+et(a)*p(2))/8
        j(:,1)=j(:,1)+x(:,a)*dn(1); j(:,2)=j(:,2)+x(:,a)*dn(2); j(:,3)=j(:,3)+x(:,a)*dn(3)
      enddo
      g8=j(1,1)*(j(2,2)*j(3,3)-j(3,2)*j(2,3)) - j(1,2)*(j(2,1)*j(3,3)-j(3,1)*j(2,3)) &
         + j(1,3)*(j(2,1)*j(3,2)-j(3,1)*j(2,2))
      dmin=min(dmin,g8)
      if(present(dmax)) dmax=max(dmax,g8)
    enddo
  end function

  ! Smallest corner Jacobian determinant of a bilinear quad (all four corners > 0 <=> valid, counter-clockwise).
  real(real64) function q4_min_jacobian(x) result(dmin)
    real(real64), intent(in) :: x(2,4)
    integer :: i, ip, im
    real(real64) :: a(2), b(2)
    dmin=huge(1._real64)
    do i=1,4
      ip=mod(i,4)+1; im=mod(i+2,4)+1
      a=x(:,ip)-x(:,i); b=x(:,im)-x(:,i)
      dmin=min(dmin,a(1)*b(2)-a(2)*b(1))
    enddo
  end function

  ! ---------------------------------------------------------------------------
  ! TEMPERATURE pilot phase 2 (creep-tem-pilot/DESIGN.md section 9, T1-T8 confirmed)
  ! ---------------------------------------------------------------------------

  ! transient_heat step loads: only type = "convection" (T4); gravity is not part of a heat step. The
  ! legacy gravity record is implied (9.81 along -mesh axis ndim, curve 0 for every group = no gravity).
  subroutine heat_step_loads(s,sp,nl,ns,srf_step)
    integer, intent(in) :: s, nl, ns
    character(*), intent(in) :: sp
    logical, intent(inout) :: srf_step(:)
    character(:), allocatable :: q, typ, nm
    integer :: l, k, idx, i, nc
    type(tem_conv_t), allocatable :: cv(:)
    real(real64) :: b
    if(.not.allocated(native_tem_step)) allocate(native_tem_step(native_blocks))
    allocate(cv(max(nl,1))); nc=0
    do l=1,nl
      q=indexed(sp//".load",l)
      k=key_index(q//".@table")
      typ=string_value(q//".type")
      select case(typ)
      case("convection")
        nm=string_value(q//".surface")
        idx=0
        do i=1,ns
          if(trim(native_surf(i)%name)==nm) idx=i
        enddo
        if(idx==0) call native_error("TEM: step "//itoa(s)//": undefined surface "//nm)
        if(srf_step(idx)) call native_error("TEM: step "//itoa(s)//": surface "//nm//" loaded twice")
        srf_step(idx)=.true.; native_surf(idx)%used=.true.
        i=sets_find(SK_SURFACE,nm)
        b=number(q//".film_coefficient")
        if(b<=0 .or. b>1.e6_real64) call native_error("TEM: step "//itoa(s)//": film_coefficient out of range (0,1e6]")
        if(.not.present_key(q//".ambient_temperature")) &
          call native_error("TEM: step "//itoa(s)//": convection ambient_temperature (amplitude name) is required")
        nc=nc+1
        cv(nc)%surf=idx; cv(nc)%first=native_surf(idx)%first; cv(nc)%last=native_surf(idx)%last
        cv(nc)%beta=b
        cv(nc)%itcurve=amp_ref(q//".ambient_temperature")
        if(native_amp(cv(nc)%itcurve)%typ/="LINEAR") call native_error("TEM: step "//itoa(s)// &
          ": convection ambient_temperature must be a linear amplitude")
      case("gravity")
        call native_error("TEM: step "//itoa(s)//": a transient_heat step takes no gravity load")
      case default
        call native_error("TEM: step "//itoa(s)//": unsupported load type "//typ//" in a transient_heat step (only convection)")
      end select
    enddo
    native_tem_step(s)%conv=cv(1:nc)
    native_loa(s)%g=9.81_real64; native_loa(s)%dir=0; native_loa(s)%dir(native_ndim)=-1._real64
    native_loa(s)%fdir=native_loa(s)%dir; native_loa(s)%gamp=0
    native_gravity(s)=9.81_real64
  end subroutine

  ! After LOA/PRE: material curve references, [[pipe]] (T5) and the face4/active-set consistency.
  subroutine heat_resolve()
    integer :: i, k, np, j, n, g
    character(:), allocatable :: p, q, nm
    do i=1,native_materials
      if(.not.native_mat_heat(i)) cycle
      native_src_curve(i)=0
      if(len_trim(mat_src_name(i))>0) then
      k=amp_by_name(trim(mat_src_name(i)))
      if(k==0) call native_error("MAT: undefined amplitude "//trim(mat_src_name(i))//" (heat_source_curve)")
      if(native_amp(k)%typ/="DABT") call native_error("MAT: heat_source_curve "//trim(mat_src_name(i))//" must be a dabt amplitude")
      native_amp(k)%used=.true.; native_src_curve(i)=k
      endif
      k=amp_by_name(trim(mat_place_name(i)))
      if(k==0) call native_error("MAT: undefined amplitude "//trim(mat_place_name(i))//" (placement_temperature_curve)")
      if(native_amp(k)%typ/="LINEAR") call native_error("MAT: placement_temperature_curve "//trim(mat_place_name(i))// &
        " must be a linear amplitude")
      native_amp(k)%used=.true.; native_place_curve(i)=k
    enddo
    ! [[pipe]] (T5): one table per legacy pipe record; algo 3 only (advection along the pipe segments)
    np=doc%count_of("pipe")
    allocate(native_pipe(np))
    do i=1,np
      p=indexed("pipe",i)
      k=key_index(p//".@table")
      call glb_check_keys(p,"@table algorithm elset coefficient flow_rate water_conductivity water_density "// &
        "water_specific_heat begin_time end_time water_temperature","segments")
      if(integer_value(p//".algorithm")/=3) call native_error("TEM: pipe "//itoa(i)//": only algorithm 3 is supported "// &
        "(legacy algo_pipe 1/2/4 are deferred)")
      nm=string_value(p//".elset")
      g=0
      do j=1,native_materials
        if(trim(native_kname(j))==nm) g=j
      enddo
      if(g==0) call native_error("TEM: pipe "//itoa(i)//": undefined elset "//nm)
      if(native_sets_on) j=sets_find(SK_ELSET,nm)
      native_pipe(i)%igroup=g
      native_pipe(i)%coef=number(p//".coefficient")
      native_pipe(i)%qw=number(p//".flow_rate")
      native_pipe(i)%lamda=number(p//".water_conductivity")
      native_pipe(i)%dens=number(p//".water_density")
      native_pipe(i)%cw=number(p//".water_specific_heat")
      native_pipe(i)%t0=number(p//".begin_time")
      native_pipe(i)%t1=number(p//".end_time")
      if(native_pipe(i)%coef<=0 .or. native_pipe(i)%coef>1.e3_real64) call native_error("TEM: pipe "//itoa(i)//": coefficient out of range (0,1e3]")
      if(native_pipe(i)%qw<=0 .or. native_pipe(i)%qw>1.e9_real64) call native_error("TEM: pipe "//itoa(i)//": flow_rate out of range (0,1e9]")
      if(native_pipe(i)%lamda<=0 .or. native_pipe(i)%lamda>1.e12_real64) &
        call native_error("TEM: pipe "//itoa(i)//": water_conductivity out of range (0,1e12]")
      if(native_pipe(i)%dens<=0 .or. native_pipe(i)%dens>1.e5_real64) call native_error("TEM: pipe "//itoa(i)//": water_density out of range (0,1e5]")
      if(native_pipe(i)%cw<=0 .or. native_pipe(i)%cw>1.e7_real64) &
        call native_error("TEM: pipe "//itoa(i)//": water_specific_heat out of range (0,1e7]")
      if(native_pipe(i)%t1<native_pipe(i)%t0) call native_error("TEM: pipe "//itoa(i)//": end_time before begin_time")
      native_pipe(i)%twcurve=amp_ref(p//".water_temperature")
      n=integer_value(p//".segments.count")
      if(n<1 .or. n>100000) call native_error("TEM: pipe "//itoa(i)//": segments must have 1..100000 entries")
      allocate(native_pipe(i)%seg(2,n))
      do j=1,n
        q=indexed(p//".segments",j)
        if(integer_value(q//".count")/=2) call native_error("TEM: pipe "//itoa(i)//": a segment is [node1, node2]")
        native_pipe(i)%seg(1,j)=integer_value(q//"[1]")
        native_pipe(i)%seg(2,j)=integer_value(q//"[2]")
        if(any(native_pipe(i)%seg(:,j)<1) .or. any(native_pipe(i)%seg(:,j)>native_points)) &
          call native_error("TEM: pipe "//itoa(i)//": segment node out of range 1.."//itoa(native_points))
        if(native_pipe(i)%seg(1,j)==native_pipe(i)%seg(2,j)) call native_error("TEM: pipe "//itoa(i)//": degenerate segment "//itoa(j))
      enddo
    enddo
  end subroutine

  ! [[tie]] (T7, HSTAR_NATIVE_NRT=toml): legacy .nrt TRAL constraint T(slave) = sum(weights*T(masters)).
  subroutine parse_nrt()
    ! [[tie]] (T7, HSTAR_NATIVE_NRT=toml): set = name of a tie set in 1.set, dof = "T".
    ! Legacy .nrt TRAL: T(slave) = sum(weights*T(masters)); slave/chain/weight checks are done by yl_sets.
    integer :: n, i, j, k, r, nt, o1, o2, q
    character(:), allocatable :: p, nm
    integer, allocatable :: idx(:)
    logical :: rot
    n=doc%count_of("tie")
    if(n>0 .and. .not.native_sets_on) call native_error("NRT: [[tie]] needs [mesh] sets_file (ties are tie sets in 1.set)")
    allocate(idx(n)); nt=0
    do i=1,n
      p=indexed("tie",i)
      k=key_index(p//".@table")
      call glb_check_keys(p,"@table set dof","")
      rot=.false.
      select case(string_value(p//".dof"))
      case("T")
        if(.not.native_heat) call native_error("NRT: tie "//itoa(i)//": dof ""T"" needs analysis.type 'transient_heat'")
      case("U")                           ! joint pilot D1 (J7): legacy .nrt translg=88 TRAL, displacements 1..ndimn
        if(native_heat) call native_error("NRT: tie "//itoa(i)//": dof ""U"" in a transient_heat analysis "// &
          "(no displacement field; legacy translg=88 ties are a no-op there)")
        if(any(native_gindex(1:native_materials)==1)) call native_error("NRT: tie "//itoa(i)// &
          ": dof ""U"" with L2 elements is deferred (rotational dofs / ROTA rows)")
        native_nrt_dof_u=.true.
      case("UR")                          ! beam pilot C2: legacy .nrt 88 ROTAX/Y/Z rows from a 1.set rotation_tie set
        if(.not.native_beam) call native_error("NRT: tie "//itoa(i)//": dof ""UR"" needs a B2 (beam) section")
        native_nrt_dof_u=.true.; rot=.true.
      case default
        call native_error("NRT: tie "//itoa(i)//": dof must be ""T"" (transient_heat), ""U"" (displacements) or ""UR"" (rotations)")
      end select
      nm=string_value(p//".set")
      if(rot) then
        idx(i)=sets_find(SK_RTIE,nm)
        if(idx(i)==0) call native_error("NRT: tie "//itoa(i)//": undefined rotation_tie set "//nm)
      else
        idx(i)=sets_find(SK_TIE,nm)
        if(idx(i)==0) call native_error("NRT: tie "//itoa(i)//": undefined tie set "//nm)
      endif
      do j=1,i-1
        if(idx(j)==idx(i)) call native_error("NRT: tie set "//nm//" listed twice")
      enddo
      nt=nt+sets(idx(i))%n
    enddo
    allocate(native_tie(nt)); k=0
    do i=1,n
      do r=1,sets(idx(i))%n
        k=k+1
        o1=sets(idx(i))%off(r); o2=sets(idx(i))%off(r+1)-1
        if(sets(idx(i))%kind==SK_RTIE) then
          native_tie(k)%slave=sets(idx(i))%m(o1)
          native_tie(k)%axis=sets(idx(i))%m(o1+1)
          if(native_ndim==2 .and. native_tie(k)%axis/=3) call native_error("NRT: rotation_tie set "// &
            trim(sets(idx(i))%name)//": 2D rotation ties need axis 3 (z)")
          q=(o2-o1-1)/2
          allocate(native_tie(k)%masters(q),native_tie(k)%comps(q),native_tie(k)%weights(q))
          do j=1,q
            native_tie(k)%masters(j)=sets(idx(i))%m(o1+2*j)
            native_tie(k)%comps(j)=sets(idx(i))%m(o1+2*j+1)
            native_tie(k)%weights(j)=sets(idx(i))%w(o1+2*j)
          enddo
          if(native_ndim==2 .and. any(native_tie(k)%comps>2)) call native_error("NRT: rotation_tie set "// &
            trim(sets(idx(i))%name)//": 2D term component must be 1 or 2")
          if(native_tie(k)%slave>native_points .or. any(native_tie(k)%masters>native_points)) &
            call native_error("NRT: rotation_tie set "//trim(sets(idx(i))%name)//": node out of range")
          cycle
        endif
        native_tie(k)%slave=sets(idx(i))%m(o1)
        native_tie(k)%masters=sets(idx(i))%m(o1+1:o2)
        native_tie(k)%weights=sets(idx(i))%w(o1+1:o2)
        if(native_tie(k)%slave>native_points .or. any(native_tie(k)%masters>native_points)) &
          call native_error("NRT: tie set "//trim(sets(idx(i))%name)//": node out of range")
      enddo
    enddo
  end subroutine

  ! ---- accessors (Temper.f90 boundt, Fem.f90 time_dependent, Global.f90 read_nrt_constraints, Solver.f90 PARDISO)
  integer function native_tem_nedge() result(n)
    n=0
    if(native_heat) n=native_nface
    tem_face_blk=0
  end function

  subroutine native_tem_face_block(iblk,sedge,nnode,index,beta,ibeta)
    integer, intent(in) :: iblk
    integer, intent(out) :: sedge,nnode,index,ibeta
    real(real64), intent(out) :: beta
    integer :: i
    do
      tem_face_blk=tem_face_blk+1
      if(native_surf(tem_face_blk)%kind==4) exit
    enddo
    sedge=native_surf(tem_face_blk)%last-native_surf(tem_face_blk)%first+1
    nnode=4; index=5; ibeta=0; beta=0._real64
    do i=1,size(native_tem_step(iblk)%conv)
      if(native_tem_step(iblk)%conv(i)%surf==tem_face_blk) beta=native_tem_step(iblk)%conv(i)%beta
    enddo
  end subroutine

  subroutine native_tem_face(tedge,lnode,aelem)
    integer, intent(in) :: tedge
    integer, intent(out) :: lnode(:),aelem
    lnode=native_faces(1:4,tedge); aelem=native_faces(5,tedge)
  end subroutine

  integer function native_tem_ntelgroup(iblk) result(n)
    integer, intent(in) :: iblk
    integer :: i
    n=0
    tem_conv_grp=0
    if(.not.native_heat) return
    do i=1,size(native_tem_step(iblk)%conv)
      n=n+native_tem_step(iblk)%conv(i)%last-native_tem_step(iblk)%conv(i)%first+1
    enddo
  end function

  subroutine native_tem_conv_group(iblk,sedge,itcurve,nline,b,e)
    integer, intent(in) :: iblk
    integer, intent(out) :: sedge,itcurve,nline,b,e
    tem_conv_grp=tem_conv_grp+1
    b=native_tem_step(iblk)%conv(tem_conv_grp)%first
    e=native_tem_step(iblk)%conv(tem_conv_grp)%last
    sedge=e-b+1; itcurve=native_tem_step(iblk)%conv(tem_conv_grp)%itcurve; nline=0
  end subroutine

  subroutine native_tem_pipes(npipe,algo)
    integer, intent(out) :: npipe,algo
    npipe=0; algo=3
    if(native_heat) npipe=size(native_pipe)
  end subroutine

  subroutine native_tem_pipe_head(ip,nline,igroup,coef)
    integer, intent(in) :: ip
    integer, intent(out) :: nline,igroup
    real(real64), intent(out) :: coef
    nline=size(native_pipe(ip)%seg,2); igroup=native_pipe(ip)%igroup; coef=native_pipe(ip)%coef
  end subroutine

  subroutine native_tem_pipe_props(ip,qw,lamda,dens,cw,t0,t1,twc)
    integer, intent(in) :: ip
    real(real64), intent(out) :: qw,lamda,dens,cw,t0,t1
    integer, intent(out) :: twc
    qw=native_pipe(ip)%qw; lamda=native_pipe(ip)%lamda; dens=native_pipe(ip)%dens; cw=native_pipe(ip)%cw
    t0=native_pipe(ip)%t0; t1=native_pipe(ip)%t1; twc=native_pipe(ip)%twcurve
  end subroutine

  subroutine native_tem_pipe_seg(ip,ie,i0,n1,n2)
    integer, intent(in) :: ip,ie
    integer, intent(out) :: i0,n1,n2
    i0=ie; n1=native_pipe(ip)%seg(1,ie); n2=native_pipe(ip)%seg(2,ie)
  end subroutine

  ! Appended to native-tem.txt from what boundt stored (tedges, tedgeload, pipeinfo).
  subroutine native_tem_heat_receipt(iblk,nedge,beta_sum,nload,itc_list,npipe,nseg,coef_sum)
    integer, intent(in) :: iblk,nedge,nload,npipe,nseg,itc_list(:)
    real(real64), intent(in) :: beta_sum,coef_sum
    integer :: u
    if(.not.native_tem_receipt_on) return
    open(newunit=u,file="native-tem.txt",status="old",position="append",action="write")
    write(u,"(a,i0,a,i0,a,es25.16,a,i0,a,i0,a,i0)") "heat block=",iblk," faces=",nedge," beta_sum=",beta_sum, &
      " loaded_faces=",nload," npipe=",npipe," segments=",nseg
    write(u,"(a,es25.16)") "pipe_edstif_coef_sum=",coef_sum
    if(size(itc_list)>0) write(u,"(a,*(1x,i0))") "load_curves=",itc_list
    close(u)
  end subroutine

  integer function native_man_nincs(iblk) result(n)
    integer, intent(in) :: iblk
    n=native_ninc(iblk)
  end function

  subroutine native_man_inc(iblk,i,miter,ditime,noutn,noutf,nstep,inc_step,nresta,nmcon,nbasef,tolf,tolv)
    integer, intent(in) :: iblk,i
    integer, intent(out) :: miter,noutn,noutf,nstep,inc_step,nresta,nmcon,nbasef
    real(real64), intent(out) :: ditime,tolf,tolv(:)
    miter=native_inc_miter(i,iblk); ditime=native_inc_dt(i,iblk); nstep=native_inc_nstep(i,iblk)
    noutn=1; noutf=1; inc_step=1; nresta=1; nmcon=0; nbasef=0     ! validated implied values (T3)
    tolf=native_inc_tol(i,iblk); tolv=native_inc_tol(i,iblk)
  end subroutine

  integer function native_nrt_ntie() result(n)
    n=0
    if(native_nrt_toml) n=size(native_tie)
  end function

  subroutine native_nrt_tie(i,title,ipoin,nintf)
    integer, intent(in) :: i
    character(*), intent(out) :: title
    integer, intent(out) :: ipoin,nintf
    title='TRAL'; ipoin=native_tie(i)%slave; nintf=size(native_tie(i)%masters)
    if(native_tie(i)%axis>0) title='ROTA'//'XYZ'(native_tie(i)%axis:native_tie(i)%axis)
  end subroutine

  ! beam pilot C2: term k of ROTA tie i as the legacy row "jpoin jdimn f0"
  subroutine native_nrt_rota_row(i,k,jpoin,jdimn,f0)
    integer, intent(in) :: i,k
    integer, intent(out) :: jpoin,jdimn
    real(real64), intent(out) :: f0
    jpoin=native_tie(i)%masters(k); jdimn=native_tie(i)%comps(k); f0=native_tie(i)%weights(k)
  end subroutine

  subroutine native_nrt_tie_data(i,listf,rintf)
    integer, intent(in) :: i
    integer, intent(out) :: listf(:)
    real(real64), intent(out) :: rintf(:)
    listf=native_tie(i)%masters; rintf=native_tie(i)%weights
  end subroutine

  ! Written from the stored trans(:) table (Global.f90 read_nrt_constraints).
  subroutine native_nrt_receipt(nties,ndof,sum_w)
    integer, intent(in) :: nties,ndof
    real(real64), intent(in) :: sum_w
    integer :: u
    if(.not.native_nrt_receipt_on) return
    open(newunit=u,file="native-nrt.txt",status="replace",action="write")
    write(u,"(a)") "native-nrt v1 mode="//trim(merge("toml    ","contract",native_nrt_toml))
    write(u,"(a,i0,a,i0,a,es25.16)") "ties=",nties," constrained_dofs=",ndof," weight_sum=",sum_w
    close(u)
  end subroutine

  ! [solver.pardiso] (sol-pilot V1/V3/V4/V5), only with HSTAR_NATIVE_SOL=toml.
  subroutine parse_sol(lin)
    character(*), intent(in) :: lin
    integer :: k
    character(:), allocatable :: s
    if(lin/='pardiso') then
      if(doc%find('solver.pardiso.@table')>0) call native_error('SOL: [solver.pardiso] needs solver.linear = "pardiso"')
      return
    endif
    if(doc%find('solver.pardiso.@table')>0) k=key_index('solver.pardiso.@table')
    native_sol_mtype=merge(11,-2,native_heat)      ! V3: heat 11 (nonsymmetric), static -2 (symmetric indefinite)
    if(doc%find('solver.pardiso.matrix_type')>0) native_sol_mtype=integer_value('solver.pardiso.matrix_type')
    if(native_sol_mtype/=-2 .and. native_sol_mtype/=11) call native_error('SOL: solver.pardiso.matrix_type must be '// &
      '-2 (real symmetric indefinite) or 11 (real nonsymmetric); other PARDISO matrix types are deferred')
    if(native_heat .and. native_sol_mtype/=11) call native_error('SOL: transient_heat needs solver.pardiso.matrix_type 11 '// &
      '(nonsymmetric matrix storage)')
    if(.not.native_heat .and. native_sol_mtype/=-2) call native_error('SOL: static needs solver.pardiso.matrix_type -2 '// &
      '(symmetric upper-triangle storage)')
    native_sol_threads=1                           ! V5: 1 thread = reproducible (byte-exact acceptance)
    if(doc%find('solver.pardiso.threads')>0) native_sol_threads=integer_value('solver.pardiso.threads')
    if(native_sol_threads<1 .or. native_sol_threads>64) call native_error('SOL: solver.pardiso.threads out of range [1,64]')
    if(doc%find('solver.pardiso.iparm')>0) then
      s=string_value('solver.pardiso.iparm')
      if(s/='default') call native_error('SOL: solver.pardiso.iparm must be "default" (explicit iparm fields are deferred, V4)')
    endif
  end subroutine

  ! .sol for PARDISO (G5 analogue, T6 option A; sol-pilot V1-V5).
  ! contract: implied record (transient_heat only): mtype 11, ncpu 2, isdefault 1 + 8 fixed iparm values.
  ! toml: mtype / threads from [solver.pardiso]; static uses the MKL default iparm (isdefault 0), transient_heat keeps
  ! the internal 8-value record (V4). MKL threads are set to "threads" (V5), so threads = 1 is reproducible.
  subroutine native_sol_pardiso(mtype,ncpu,msglvl,isdefault,reducing_order,precgs,permutation,maxiter,out_of_core, &
                                eps_pivot,iparm11,iparm13)
    integer, intent(out) :: mtype,ncpu,msglvl,isdefault,reducing_order,precgs,permutation,maxiter,out_of_core, &
                            eps_pivot,iparm11,iparm13
    if(.not.native_heat .and. .not.native_pardiso) call native_glb_deferred('MAIN_PARDISO (.sol) outside transient_heat')
    mtype=11; ncpu=2; msglvl=0; isdefault=1
    reducing_order=2; precgs=0; permutation=0; maxiter=20; out_of_core=0; eps_pivot=20; iparm11=0; iparm13=0
    if(native_sol_toml) then
      mtype=native_sol_mtype; ncpu=native_sol_threads
      isdefault=merge(1,0,native_heat)
      call mkl_set_num_threads(ncpu)
    endif
    call native_sol_receipt('PARDISO',mtype,ncpu,msglvl,isdefault)
  end subroutine

  ! native-sol.txt: written at the first solver SET call (every block uses the same record, V7).
  subroutine native_sol_receipt(kind,a,b,c,d)
    character(*), intent(in) :: kind
    integer, intent(in) :: a,b,c,d
    integer :: u
    if(.not.native_sol_receipt_on) return
    sol_receipt_n=sol_receipt_n+1
    if(sol_receipt_n>1) return
    open(newunit=u,file="native-sol.txt",status="replace",action="write")
    write(u,"(a)") "native-sol v1 mode="//trim(merge("toml    ","contract",native_sol_toml))
    if(kind=='PROFILE') then
      write(u,"(a,4(a,i0))") "solver=PROFILE"," iafile=",a," icond=",b," ipdchk=",c," ising=",d
    else
      write(u,"(a,4(a,i0))") "solver=PARDISO"," mtype=",a," ncpu=",b," msglvl=",c," isdefault=",d
      if(.not.native_sol_toml) then
        write(u,"(a)") "threads=environment deterministic=unknown (contract record; MKL threads from the environment)"
      elseif(b==1) then
        write(u,"(a)") "threads=1 deterministic=yes (byte-exact acceptance)"
      else
        write(u,"(a,i0,a)") "threads=",b," deterministic=no (threads>1: PARDISO results are not byte-reproducible; "// &
          "acceptance by numerical tolerance only, sol-pilot V5)"
      endif
    endif
    close(u)
  end subroutine

  subroutine native_error(message)
    character(*), intent(in) :: message
    write(error_unit,'(a)') 'INVALID_INPUT TOML-03/04: '//message
    stop 2
  end subroutine

  integer function key_index(path) result(k)
    character(*), intent(in) :: path
    k=doc%find(path)
    if(k == 0) call native_error('missing: '//path)
    seen(k)=.true.
  end function

  real(real64) function number(path) result(v)
    character(*), intent(in) :: path
    integer :: k
    k=key_index(path)
    select case(doc%entry(k)%kind)
    case(TV_INT); v=real(doc%entry(k)%ivalue,real64)
    case(TV_REAL); v=doc%entry(k)%rvalue
    case default; call native_error('expected number: '//path)
    end select
    if(.not.ieee_is_finite(v)) call native_error('non-finite: '//path)
  end function

  integer function integer_value(path) result(v)
    character(*), intent(in) :: path
    integer :: k
    k=key_index(path)
    if(doc%entry(k)%kind /= TV_INT) call native_error('expected integer: '//path)
    v=doc%entry(k)%ivalue
  end function

  subroutine expect_int(path,v)
    character(*), intent(in) :: path
    integer, intent(in) :: v
    if(pre_owned(path) .or. loa_owned(path) .or. mat_owned(path) .or. man_owned(path)) return
    if(integer_value(path) /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_real(path,v)
    character(*), intent(in) :: path
    real(real64), intent(in) :: v
    if(pre_owned(path) .or. loa_owned(path) .or. mat_owned(path) .or. man_owned(path)) return
    if(number(path) /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_str(path,v)
    character(*), intent(in) :: path,v
    integer :: k
    if(pre_owned(path) .or. loa_owned(path) .or. mat_owned(path) .or. man_owned(path)) return
    k=key_index(path)
    if(doc%entry(k)%kind /= TV_STR) call native_error('expected string: '//path)
    if(doc%entry(k)%svalue /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_bool(path,v)
    character(*), intent(in) :: path
    logical, intent(in) :: v
    integer :: k
    if(pre_owned(path) .or. loa_owned(path) .or. mat_owned(path) .or. man_owned(path)) return
    k=key_index(path)
    if(doc%entry(k)%kind /= TV_BOOL) call native_error('expected boolean: '//path)
    if(doc%entry(k)%lvalue .neqv. v) call native_error('unsupported value: '//path)
  end subroutine

  integer function token_count(row) result(n)
    character(*), intent(in) :: row
    integer :: k
    logical :: in_token, space
    if(len_trim(row) >= len(row)) call native_error('mesh row too long')
    if(verify(trim(row),'0123456789.eE+- '//achar(9))/=0) call native_error('invalid mesh token')
    n=0; in_token=.false.
    do k=1,len_trim(row)
      space=row(k:k)==' ' .or. row(k:k)==achar(9)
      if(.not.space .and. .not.in_token) n=n+1
      in_token=.not.space
    enddo
  end function

  real(real64) function mesh_round6(x) result(y)
    real(real64), intent(in) :: x
    real(real64) :: scaled
    scaled=x*1.e6_real64
    ! Frozen generator uses Python decimal output: ties round to even, not away.
    if(scaled-floor(scaled)==0.5_real64) then
      y=real(2*nint(scaled/2),real64)/1.e6_real64
    else
      y=real(nint(scaled),real64)/1.e6_real64
    endif
  end function

  ! This milestone supports only the established 16x16 Cook mesh, including numbering.
  ! Check geometry/topology before any output is opened, not just row counts.
  subroutine validate_mesh()
    integer :: u,ios,k,id,ix,iy,mat,n(4),expected(4)
    real(real64) :: x,y,xx,yy
    character(1024) :: row
    open(newunit=u,file=trim(native_mesh_prefix)//'.cor',status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh 1.cor')
    do k=1,native_points
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) call native_error('missing mesh coordinate row')
      if(token_count(row) /= 3) call native_error('coordinate row must have 3 columns')
      read(row,*,iostat=ios) id,x,y
      if(ios /= 0) call native_error('invalid mesh coordinates')
      if(native_dam) then
        if(k<=1025) then
          ix=mod(k-1,41); iy=(k-1)/41
          if(ix<=12) then
            xx=80._real64*ix/12
          elseif(ix<=28) then
            xx=80._real64+30._real64*(ix-12)/16
          else
            xx=110._real64+80._real64*(ix-28)/12
          endif
          yy=-100._real64+100._real64*iy/24
        else
          ix=mod(k-1026,17); iy=(k-1026)/17+1
          yy=1.25_real64*iy
          xx=80._real64+(30._real64-0.5_real64*yy)*ix/16
        endif
        ! The frozen dam mesh was written with six decimal places.
        xx=mesh_round6(xx)
        yy=mesh_round6(yy)
      else
        ix=mod(k-1,17); iy=(k-1)/17
        xx=3._real64*ix
        yy=44._real64*ix/16 + (44._real64-28._real64*ix/16)*iy/16
      endif
      if(.not.ieee_is_finite(x) .or. .not.ieee_is_finite(y)) call native_error('non-finite mesh')
      if(id /= k .or. abs(x-xx)>1.e-7_real64 .or. abs(y-yy)>1.e-7_real64) &
        call native_error('unsupported mesh coordinates/numbering')
    enddo
    read(u,'(a)',iostat=ios) row
    if(ios >= 0) call native_error('extra mesh coordinate row')
    close(u)
    open(newunit=u,file=trim(native_mesh_prefix)//'.ele',status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh 1.ele')
    do k=1,native_elements
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) call native_error('missing mesh element row')
      if(token_count(row) /= 6) call native_error('element row must have 6 columns')
      read(row,*,iostat=ios) id,n,mat
      if(ios /= 0) call native_error('invalid mesh connectivity')
      if(native_dam) then
        if(k<=960) then
          ix=mod(k-1,40); iy=(k-1)/40
          expected=[1+ix+41*iy,2+ix+41*iy,43+ix+41*iy,42+ix+41*iy]
          if(mat/=1) call native_error('foundation mesh group must be 1')
        else
          ix=mod(k-961,16); iy=(k-961)/16
          if(iy==0) then
            expected=[997+ix,998+ix,1027+ix,1026+ix]
          else
            expected=[1009+ix+17*iy,1010+ix+17*iy,1027+ix+17*iy,1026+ix+17*iy]
          endif
          if(mat/=2) call native_error('dam mesh group must be 2')
        endif
      else
        ix=mod(k-1,16); iy=(k-1)/16
        expected=[1+ix+17*iy,2+ix+17*iy,19+ix+17*iy,18+ix+17*iy]
        if(mat/=1) call native_error('Cook mesh group must be 1')
      endif
      if(id/=k .or. any(n/=expected)) call native_error('unsupported mesh connectivity/numbering')
    enddo
    read(u,'(a)',iostat=ios) row
    if(ios >= 0) call native_error('extra mesh element row')
    close(u)
  end subroutine
end module native_input
