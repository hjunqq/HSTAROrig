! TOML-03/04: two bounded profiles. No solver state mirror or legacy-card emitter.
module native_input
  use iso_fortran_env, only: real64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use yl_authoring_toml
  implicit none
  character(1024) :: native_mesh_prefix='1'
  logical :: native_mode = .false.
  logical :: native_dam = .false.
  integer :: native_points=289, native_elements=256, native_materials=1, native_blocks=1
  integer :: native_elcounts(2)=[256,0], native_appearance(2,2)=1, native_reset(2)=0
  integer :: native_bc(50,2)=0, native_nbc(2)=[17,17], native_edges(3,40)
  real(real64) :: native_e(2), native_nu(2), native_density(2), native_gravity(2)
  real(real64) :: native_pressure=0.
  integer :: native_substeps(2), native_iterations(2)
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
    real(real64), allocatable :: t(:), f(:)
    logical :: used=.false.
  end type
  type(loa_amp_t), allocatable :: native_amp(:)
  type, private :: loa_surf_t
    character(len=TOML_LEN_TEXT) :: name=""
    integer :: first=0, last=0
    logical :: used=.false.
  end type
  type(loa_surf_t), allocatable, private :: native_surf(:)
  integer :: native_loa_nedge=0
  integer, allocatable :: native_loa_edges(:,:)
  type :: loa_hyd_t
    integer :: first, last, itcurve
    real(real64) :: cor0, cor1, p0, p1, fact
  end type
  type :: loa_step_t
    real(real64) :: g=0, dir(2)=0, fdir(2)=0
    integer :: gamp=0
    type(loa_hyd_t), allocatable :: hyd(:)
  end type
  type(loa_step_t), allocatable :: native_loa(:)
  type :: loa_pl_t
    integer :: itcurve=0
    real(real64) :: comp(2)=0
    integer, allocatable :: nodes(:)
  end type
  type(loa_pl_t), allocatable :: native_pl(:)
  integer, private :: loa_edge_group=0
  character(len=TOML_LEN_TEXT), allocatable, private :: loa_nset_refs(:)
contains
  subroutine native_start()
    character(1024) :: flag, filename
    character(256) :: pre_env, loa_env
    integer :: n, i, j, file_size, ios
    n = command_argument_count()
    if(n == 0) return
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
    i=key_index('case.name')
    if(doc%entry(i)%kind/=TV_STR) call native_error('expected string: case.name')
    select case(trim(doc%entry(i)%svalue))
    case('cooks_membrane')
      include 'NativeContract.inc'
      native_bc(1:17,1)=native_nodes
      native_bc(1:17,2)=native_nodes
    case('train01_gravdam_static')
      native_dam=.true.
      include 'NativeDamContract.inc'
      native_points=1705; native_elements=1600; native_materials=2; native_blocks=2
      native_elcounts=[960,640]; native_appearance(:,1)=[1,0]; native_appearance(:,2)=[1,1]
      native_reset=[0,1]; native_nbc=[50,41]
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
    do i=1,native_materials
      native_e(i)=number(indexed('material',i)//'.E')
      native_nu(i)=number(indexed('material',i)//'.nu')
      native_density(i)=number(indexed('material',i)//'.density')
      if(native_e(i)<=0 .or. native_e(i)>1.e15_real64) call native_error('E out of range (0,1e15] Pa')
      if(native_nu(i)<=-1 .or. native_nu(i)>=0.5_real64) call native_error('nu out of range (-1,0.5)')
      if(native_density(i)<=0 .or. native_density(i)>1.e6_real64) call native_error('density out of range (0,1e6]')
    enddo
    do i=1,native_blocks
      if(.not.native_loa_toml) native_gravity(i)=number(indexed('step',i)//'.load[1].magnitude')
      native_substeps(i)=integer_value(indexed('step',i)//'.controls.substeps')
      native_iterations(i)=integer_value(indexed('step',i)//'.controls.max_iterations')
      if(.not.native_loa_toml .and. (native_gravity(i)<=0 .or. native_gravity(i)>1.e4_real64)) call native_error('gravity out of range (0,1e4]')
      if(native_substeps(i)<1 .or. native_substeps(i)>100) call native_error('substeps out of range [1,100]')
      if(native_iterations(i)<1 .or. native_iterations(i)>100) call native_error('max_iterations out of range [1,100]')
    enddo
    if(native_loa_toml) call parse_loa()
    if(native_pre_toml) call parse_pre()
    if(native_loa_toml) call loa_check_used()
    native_pre_curve1 = (.not.native_pre_toml .and. native_dam) .or. pre_used_curve1
    do i=1,doc%n
      if(.not.seen(i)) call native_error('unsupported key: '//trim(doc%entry(i)%path))
    enddo
    call validate_mesh()
    native_mode=.true.
    print *, 'Native TOML input accepted (mesh relative to TOML file)'
    print *, 'TOML material E,nu,density:',native_e(1:native_materials),native_nu(1:native_materials),native_density(1:native_materials)
    print *, 'TOML gravity,substeps,max_iterations:',native_gravity(1:native_blocks),native_substeps(1:native_blocks),native_iterations(1:native_blocks)
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
    character(*), parameter :: dofname(3)=["ux","uy","uz"]
    character(:), allocatable :: p, q, sp, nm, dn, amp
    integer :: nn, i, j, k, n, v, s, nb, b, nd, d, idx, dof, kk, itc, outf, ncand, ndim_pre
    logical, allocatable :: nref(:)
    real(real64) :: val
    ndim_pre=2
    nn=doc%count_of("nset")
    if(nn<1) call native_error("at least one [[nset]] is required")
    allocate(ns(nn))
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
      if(nb<1) call native_error("step "//itoa(s)//": at least one [[step.boundary]] is required")
      allocate(owner(3,native_points)); owner=0
      allocate(cand(3*nb)); ncand=0
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
          case("ux"); dof=1
          case("uy"); dof=2
          case("uz"); dof=3
          case("rx","ry","rz")
            call native_error("rotational dof "//dn//" is not supported in the PRE pilot (step "//itoa(s)//" boundary "//itoa(b)//")")
          case default
            call native_error("unknown dof "//dn//" (expected ux, uy, uz) in step "//itoa(s)//" boundary "//itoa(b))
          end select
          if(dof>ndim_pre) call native_error("dof "//dn//" not available: mesh.dimension is "//itoa(ndim_pre))
          do k=1,size(ns(idx)%nodes)
            v=ns(idx)%nodes(k)
            if(owner(dof,v)/=0) call native_error("duplicate constraint: step "//itoa(s)//" node "//itoa(v)//" dof "//dofname(dof)// &
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
    integer :: e(3)
    real(real64) :: wl, bt, gw, g, d(2), fd(2)
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
      if(typ/="linear") call native_error("LOA: amplitude "//nm//": type "//typ//" is not supported in the LOA pilot (only linear)")
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
    ! ---- [[surface]] (edge sets, D5)
    ns=doc%count_of("surface")
    allocate(native_surf(ns))
    total=0
    do i=1,ns
      p=indexed("surface",i)
      total=total+integer_value(p//".edges.count")
    enddo
    allocate(native_loa_edges(3,max(total,1)))
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
        native_loa_edges(:,total)=e
      enddo
      native_surf(i)%last=total
    enddo
    ! ---- [[step]] loads
    if(doc%count_of("step")/=native_blocks) call native_error("LOA: number of [[step]] must be "//itoa(native_blocks))
    allocate(native_loa(native_blocks))
    allocate(srf_step(max(ns,1)))
    do s=1,native_blocks
      sp=indexed("step",s)
      nl=doc%count_of(sp//".load")
      allocate(tmp(max(nl,1)))
      ng=0; nh=0; srf_step=.false.
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
          if(integer_value(q//".direction.count")/=2) call native_error("LOA: step "//itoa(s)//": gravity direction needs 2 components")
          d(1)=number(q//".direction[1]"); d(2)=number(q//".direction[2]")
          if(sum(d**2)==0) call native_error("LOA: step "//itoa(s)//": zero gravity direction")
          fd=d
          if(present_key(q//".fluid_direction.count")) then
            if(integer_value(q//".fluid_direction.count")/=2) call native_error("LOA: step "//itoa(s)//": fluid_direction needs 2 components")
            fd(1)=number(q//".fluid_direction[1]"); fd(2)=number(q//".fluid_direction[2]")
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
          wl=number(q//".water_level"); bt=number(q//".bottom"); gw=number(q//".unit_weight")
          if(wl<=bt) call native_error("LOA: step "//itoa(s)//": water_level must be above bottom")
          if(gw<0 .or. gw>1.e6_real64) call native_error("LOA: step "//itoa(s)//": unit_weight out of range [0,1e6]")
          if(.not.present_key(q//".amplitude")) call native_error("LOA: step "//itoa(s)//": hydrostatic amplitude is required")
          nh=nh+1
          tmp(nh)=loa_hyd_t(native_surf(idx)%first,native_surf(idx)%last,amp_ref(q//".amplitude"),wl,bt,0._real64,wl-bt,gw)
        case("pressure")
          call native_error("LOA: step "//itoa(s)//": load type pressure is not accepted with HSTAR_NATIVE_LOA=toml; use hydrostatic")
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
      idx=0
      do i=1,nnset
        if(string_value(indexed("nset",i)//".name")==nm) idx=i
      enddo
      if(idx==0) call native_error("LOA: load "//itoa(ip)//": undefined nset "//nm)
      loa_nset_refs(ip)=nm
      k=key_index(indexed("nset",idx)//".@table")
      k=integer_value(indexed("nset",idx)//".nodes.count")
      if(k<1) call native_error("LOA: nset "//nm//" has no nodes")
      allocate(native_pl(ip)%nodes(k))
      do j=1,k
        ii=integer_value(indexed(indexed("nset",idx)//".nodes",j))
        if(ii<1 .or. ii>native_points) call native_error("LOA: node out of range in nset "//nm)
        native_pl(ip)%nodes(j)=ii
      enddo
      if(integer_value(p//".components.count")/=2) call native_error("LOA: load "//itoa(ip)//": components must have mesh.dimension (2) entries")
      native_pl(ip)%comp(1)=number(p//".components[1]"); native_pl(ip)%comp(2)=number(p//".components[2]")
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
    if(native_loa_toml) ntime=size(native_amp(ic)%t)
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
    endif
  end subroutine

  subroutine native_loa_edge(tedge,lnode,aelem)
    integer, intent(in) :: tedge
    integer, intent(out) :: lnode(:),aelem
    if(native_loa_toml) then
      lnode=native_loa_edges(1:2,tedge); aelem=native_loa_edges(3,tedge)
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
    endif
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
      fg=native_loa(iblk)%dir; ff=native_loa(iblk)%fdir
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
    if(native_dam) tcg=merge(k,0,native_appearance(1:size(tcg),iblk)>0)
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
    if(pre_owned(path) .or. loa_owned(path)) return
    if(integer_value(path) /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_real(path,v)
    character(*), intent(in) :: path
    real(real64), intent(in) :: v
    if(pre_owned(path) .or. loa_owned(path)) return
    if(number(path) /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_str(path,v)
    character(*), intent(in) :: path,v
    integer :: k
    if(pre_owned(path) .or. loa_owned(path)) return
    k=key_index(path)
    if(doc%entry(k)%kind /= TV_STR) call native_error('expected string: '//path)
    if(doc%entry(k)%svalue /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_bool(path,v)
    character(*), intent(in) :: path
    logical, intent(in) :: v
    integer :: k
    if(pre_owned(path) .or. loa_owned(path)) return
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
