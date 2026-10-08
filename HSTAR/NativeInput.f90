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
contains
  subroutine native_start()
    character(1024) :: flag, filename
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
      native_pressure=number('step[2].load[2].distribution.scale')
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
      native_gravity(i)=number(indexed('step',i)//'.load[1].magnitude')
      native_substeps(i)=integer_value(indexed('step',i)//'.controls.substeps')
      native_iterations(i)=integer_value(indexed('step',i)//'.controls.max_iterations')
      if(native_gravity(i)<=0 .or. native_gravity(i)>1.e4_real64) call native_error('gravity out of range (0,1e4]')
      if(native_substeps(i)<1 .or. native_substeps(i)>100) call native_error('substeps out of range [1,100]')
      if(native_iterations(i)<1 .or. native_iterations(i)>100) call native_error('max_iterations out of range [1,100]')
    enddo
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
    if(integer_value(path) /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_real(path,v)
    character(*), intent(in) :: path
    real(real64), intent(in) :: v
    if(number(path) /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_str(path,v)
    character(*), intent(in) :: path,v
    integer :: k
    k=key_index(path)
    if(doc%entry(k)%kind /= TV_STR) call native_error('expected string: '//path)
    if(doc%entry(k)%svalue /= v) call native_error('unsupported value: '//path)
  end subroutine
  subroutine expect_bool(path,v)
    character(*), intent(in) :: path
    logical, intent(in) :: v
    integer :: k
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
