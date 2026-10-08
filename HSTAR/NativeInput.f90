! TOML-03: single-step Cook profile. No solver state mirror or legacy-card emitter.
module native_input
  use iso_fortran_env, only: real64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use yl_authoring_toml
  implicit none
  character(1024) :: native_mesh_prefix='1'
  logical :: native_mode = .false.
  real(real64) :: native_e, native_nu, native_density, native_gravity
  integer :: native_substeps, native_iterations
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
    inquire(file=trim(filename),size=file_size,iostat=ios)
    if(ios /= 0 .or. file_size < 0) call native_error('cannot open TOML input')
    if(file_size > 65536) call native_error('TOML exceeds Cook profile limit (64 KiB)')
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
    include 'NativeContract.inc'
    native_e = number('material[1].E')
    native_nu = number('material[1].nu')
    native_density = number('material[1].density')
    native_gravity = number('step[1].load[1].magnitude')
    native_substeps = integer_value('step[1].controls.substeps')
    native_iterations = integer_value('step[1].controls.max_iterations')
    if(native_e <= 0 .or. native_e > 1.e15_real64) call native_error('E out of range (0,1e15] Pa')
    if(native_nu <= -1 .or. native_nu >= 0.5_real64) call native_error('nu out of range (-1,0.5)')
    if(native_density <= 0 .or. native_density > 1.e6_real64) call native_error('density out of range (0,1e6]')
    if(native_gravity <= 0 .or. native_gravity > 1.e4_real64) call native_error('gravity out of range (0,1e4]')
    if(native_substeps < 1 .or. native_substeps > 100) call native_error('substeps out of range [1,100]')
    if(native_iterations < 1 .or. native_iterations > 100) call native_error('max_iterations out of range [1,100]')
    do i=1,doc%n
      if(.not.seen(i)) call native_error('unsupported key: '//trim(doc%entry(i)%path))
    enddo
    call validate_mesh()
    native_mode=.true.
    print *, 'TOML-03 native Cook input accepted (mesh relative to TOML file)'
    print *, 'TOML material E,nu,density:',native_e,native_nu,native_density
    print *, 'TOML gravity,substeps,max_iterations:',native_gravity,native_substeps,native_iterations
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
    write(error_unit,'(a)') 'INVALID_INPUT TOML-03: '//message
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

  ! This milestone supports only the established 16x16 Cook mesh, including numbering.
  ! Check geometry/topology before any output is opened, not just row counts.
  subroutine validate_mesh()
    integer :: u,ios,k,id,ix,iy,mat,n(4),expected(4)
    real(real64) :: x,y,xx,yy
    character(1024) :: row
    open(newunit=u,file=trim(native_mesh_prefix)//'.cor',status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh 1.cor')
    do k=1,289
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) call native_error('missing mesh coordinate row')
      if(token_count(row) /= 3) call native_error('coordinate row must have 3 columns')
      read(row,*,iostat=ios) id,x,y
      if(ios /= 0) call native_error('invalid mesh coordinates')
      ix=mod(k-1,17); iy=(k-1)/17
      xx=3._real64*ix
      yy=44._real64*ix/16 + (44._real64-28._real64*ix/16)*iy/16
      if(.not.ieee_is_finite(x) .or. .not.ieee_is_finite(y)) call native_error('non-finite mesh')
      if(id /= k .or. abs(x-xx)>1.e-7_real64 .or. abs(y-yy)>1.e-7_real64) &
        call native_error('unsupported mesh coordinates/numbering')
    enddo
    read(u,'(a)',iostat=ios) row
    if(ios >= 0) call native_error('extra mesh coordinate row')
    close(u)
    open(newunit=u,file=trim(native_mesh_prefix)//'.ele',status='old',action='read',iostat=ios)
    if(ios /= 0) call native_error('cannot open mesh 1.ele')
    do k=1,256
      read(u,'(a)',iostat=ios) row
      if(ios /= 0) call native_error('missing mesh element row')
      if(token_count(row) /= 6) call native_error('element row must have 6 columns')
      read(row,*,iostat=ios) id,n,mat
      if(ios /= 0) call native_error('invalid mesh connectivity')
      ix=mod(k-1,16); iy=(k-1)/16
      expected=[1+ix+17*iy,2+ix+17*iy,19+ix+17*iy,18+ix+17*iy]
      if(id /= k .or. mat /= 1 .or. any(n /= expected)) call native_error('unsupported mesh connectivity/numbering')
    enddo
    read(u,'(a)',iostat=ios) row
    if(ios >= 0) call native_error('extra mesh element row')
    close(u)
  end subroutine
end module native_input
