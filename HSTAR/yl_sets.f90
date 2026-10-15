! 1.set reader (sets-pilot/DESIGN.md, S1-S10 confirmed 2026-10-14).
! Text, versioned ("hstar-sets 1"), "end" marker, per-set counts, generate ranges, full validation.
! Kinds: nset, elset (members; optional "generate" triples start end step), surface (pairs
! elem local-id; "face" = B8 faces 1..6, "edge" = Q4 edges 1..4), tie (slave nmaster masters weights).
! No dependency on native_input: errors come back in msg (caller stops with exit code 2).
module yl_sets
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: SK_NSET=1, SK_ELSET=2, SK_SURFACE=3, SK_TIE=4
  integer, parameter, public :: SS_PLAIN=0, SS_FACE=1, SS_EDGE=2
  integer, parameter, public :: SETS_LEN_NAME=32, SETS_MAX_LINE=4096
  character(7), parameter, public :: SETS_KIND_NAME(4)=['nset   ','elset  ','surface','tie    ']
  type, public :: set_t
    integer :: kind=0, sub=SS_PLAIN, n=0
    character(SETS_LEN_NAME) :: name=''
    logical :: generated=.false., referenced=.false.
    integer, allocatable :: m(:)          ! nset/elset members; surface: m(2*i-1)=elem, m(2*i)=local id
    integer, allocatable :: off(:)        ! tie: record i occupies m(off(i):off(i+1)-1) = slave, masters...
    real(real64), allocatable :: w(:)     ! tie weights, aligned with the masters
  end type
  type(set_t), allocatable, public :: sets(:)
  integer, public :: sets_n=0, sets_npoin=0, sets_nelem=0
  integer, allocatable, public :: sets_elnn(:), sets_elnodes(:,:)
  ! B8 local faces (verified on train_temp_creep: 1672/1672 unique, canonical order byte-identical)
  integer, parameter, public :: B8_FACE(4,6)=reshape([1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8],[4,6])
  ! tokenizer state
  integer :: u_=0, line_=0, pos_=1
  character(SETS_MAX_LINE) :: buf_=''
  character(64) :: back_=''
  logical :: has_back_=.false., eof_=.false.
  public :: sets_load, sets_find, sets_first_unreferenced, sets_edge_nodes, sets_face_nodes, sets_receipt, sets_reset
contains
  subroutine sets_reset()
    if(allocated(sets)) deallocate(sets)
    if(allocated(sets_elnn)) deallocate(sets_elnn)
    if(allocated(sets_elnodes)) deallocate(sets_elnodes)
    sets_n=0; sets_npoin=0; sets_nelem=0
  end subroutine

  ! prefix: mesh prefix ("dir/1"); setfile: path of 1.set. msg='' on success.
  subroutine sets_load(prefix,setfile,msg)
    character(*), intent(in) :: prefix, setfile
    character(:), allocatable, intent(out) :: msg
    character(64) :: tok
    integer :: ios, k
    type(set_t), allocatable :: tmp(:)
    msg=''
    call sets_reset()
    call load_mesh(prefix,msg); if(len(msg)>0) return
    open(newunit=u_,file=setfile,status='old',action='read',iostat=ios)
    if(ios/=0) then; msg='SETS: cannot open '//trim(setfile); return; endif
    line_=0; pos_=1; buf_=''; has_back_=.false.; eof_=.false.
    allocate(sets(16))
    call next(tok,msg); if(len(msg)>0) goto 9
    if(tok/='hstar-sets') then; msg='SETS: first line must be "hstar-sets 1"'; goto 9; endif
    call next(tok,msg); if(len(msg)>0) goto 9
    if(tok/='1') then; msg='SETS: unsupported format version '//trim(tok)//' (this program reads version 1)'; goto 9; endif
    do
      call next(tok,msg); if(len(msg)>0) goto 9
      if(eof_) then; msg='SETS: missing end marker (file truncated?)'; goto 9; endif
      if(tok=='end') then
        call next(tok,msg); if(len(msg)>0) goto 9
        if(.not.eof_) msg='SETS: data after end marker (line '//itoa(line_)//')'
        exit
      endif
      k=0
      select case(trim(tok))
      case('nset'); k=SK_NSET
      case('elset'); k=SK_ELSET
      case('surface'); k=SK_SURFACE
      case('tie'); k=SK_TIE
      case default
        if(is_int(tok) .and. sets_n>0) then
          msg='SETS: '//trim(SETS_KIND_NAME(sets(sets_n)%kind))//' '//trim(sets(sets_n)%name)// &
              ': more entries than its count '//itoa(sets(sets_n)%n)//' (line '//itoa(line_)//')'
        else
          msg='SETS: unknown set kind "'//trim(tok)//'" (line '//itoa(line_)//'; expected nset, elset, surface, tie or end)'
        endif
        goto 9
      end select
      if(sets_n==size(sets)) then
        allocate(tmp(2*size(sets))); tmp(1:sets_n)=sets(1:sets_n); call move_alloc(tmp,sets)
      endif
      sets_n=sets_n+1
      sets(sets_n)%kind=k
      call read_set(sets(sets_n),msg); if(len(msg)>0) goto 9
    enddo
    if(len(msg)==0) call check_ties(msg)
9   close(u_)
  end subroutine

  subroutine load_mesh(prefix,msg)
    character(*), intent(in) :: prefix
    character(:), allocatable, intent(inout) :: msg
    integer :: u, ios, nt, k, id
    character(1024) :: row
    sets_npoin=0
    open(newunit=u,file=trim(prefix)//'.cor',status='old',action='read',iostat=ios)
    if(ios/=0) then; msg='SETS: cannot open mesh 1.cor'; return; endif
    do
      read(u,'(a)',iostat=ios) row
      if(ios/=0) exit
      if(len_trim(row)>0) sets_npoin=sets_npoin+1
    enddo
    close(u)
    sets_nelem=0
    open(newunit=u,file=trim(prefix)//'.ele',status='old',action='read',iostat=ios)
    if(ios/=0) then; msg='SETS: cannot open mesh 1.ele'; return; endif
    do
      read(u,'(a)',iostat=ios) row
      if(ios/=0) exit
      if(len_trim(row)>0) sets_nelem=sets_nelem+1
    enddo
    rewind(u)
    allocate(sets_elnn(sets_nelem), sets_elnodes(8,sets_nelem)); sets_elnn=0; sets_elnodes=0
    k=0
    do
      read(u,'(a)',iostat=ios) row
      if(ios/=0) exit
      if(len_trim(row)==0) cycle
      k=k+1
      nt=ntok(row)
      select case(nt)
      case(3,4); sets_elnn(k)=2          ! matline L1: 2-node bar (L2)
      case(5,6); sets_elnn(k)=4
      case(9,10); sets_elnn(k)=8
      case default
        msg='SETS: 1.ele row '//itoa(k)//': only L2 (2 nodes), Q4 (4 nodes) and B8 (8 nodes) elements are supported'; close(u); return
      end select
      read(row,*,iostat=ios) id, sets_elnodes(1:sets_elnn(k),k)
      if(ios/=0 .or. id/=k) then; msg='SETS: 1.ele row '//itoa(k)//': invalid row'; close(u); return; endif
    enddo
    close(u)
  end subroutine

  subroutine read_set(s,msg)
    type(set_t), intent(inout) :: s
    character(:), allocatable, intent(inout) :: msg
    character(64) :: tok
    integer :: i, j, c, a, b, st, got, lim, nm, v, cap
    integer, allocatable :: mi(:)
    logical, allocatable :: mark(:), fmark(:,:)
    real(real64) :: wsum, x
    character(:), allocatable :: who
    call next(tok,msg); if(len(msg)>0) return
    if(eof_ .or. .not.valid_name(tok)) then
      msg='SETS: invalid set name "'//trim(tok)//'" (line '//itoa(line_)//'; letters, digits, _ ; at most 32)'; return
    endif
    s%name=tok
    who='SETS: '//trim(SETS_KIND_NAME(s%kind))//' '//trim(s%name)
    do i=1,sets_n-1
      if(sets(i)%kind==s%kind .and. sets(i)%name==s%name) then; msg=who//': duplicate name'; return; endif
    enddo
    call next(tok,msg); if(len(msg)>0) return
    if(.not.is_int(tok)) then; msg=who//': missing count'; return; endif
    read(tok,*) c
    if(c<1) then; msg=who//': count must be >= 1'; return; endif
    s%n=c
    call next(tok,msg); if(len(msg)>0) return
    select case(s%kind)
    case(SK_NSET,SK_ELSET)
      if(tok=='generate') then
        s%generated=.true.
      else
        call push(tok)
      endif
    case(SK_SURFACE)
      select case(trim(tok))
      case('face'); s%sub=SS_FACE
      case('edge'); s%sub=SS_EDGE
      case default; msg=who//': surface type must be face or edge'; return
      end select
    case(SK_TIE)
      call push(tok)
    end select
    lim=merge(sets_npoin,sets_nelem,s%kind/=SK_ELSET .and. s%kind/=SK_SURFACE)
    select case(s%kind)
    case(SK_NSET,SK_ELSET)
      allocate(s%m(c), mark(lim)); mark=.false.
      got=0
      do while(got<c)
        if(s%generated) then
          call num(a); if(len(msg)>0) return
          call num(b); if(len(msg)>0) return
          call num(st); if(len(msg)>0) return
          if(st<1 .or. b<a .or. mod(b-a,st)/=0) then
            msg=who//': generate needs start <= end, step >= 1 and end = start + k*step ('//itoa(a)//' '//itoa(b)//' '//itoa(st)//')'; return
          endif
          if(got+(b-a)/st+1>c) then; msg=who//': generate ranges give more members than the count '//itoa(c); return; endif
          do v=a,b,st
            call add(v)
            if(len(msg)>0) return
          enddo
        else
          call num(v); if(len(msg)>0) return
          call add(v); if(len(msg)>0) return
        endif
      enddo
    case(SK_SURFACE)
      allocate(s%m(2*c), fmark(6,sets_nelem)); fmark=.false.
      do i=1,c
        call num(a); if(len(msg)>0) return
        call num(b); if(len(msg)>0) return
        if(a<1 .or. a>sets_nelem) then; msg=who//': element out of range 1..'//itoa(sets_nelem)//': '//itoa(a); return; endif
        if(s%sub==SS_FACE) then
          if(sets_elnn(a)/=8) then; msg=who//': face surfaces need B8 elements (element '//itoa(a)//')'; return; endif
          if(b<1 .or. b>6) then; msg=who//': B8 face id must be 1..6 (element '//itoa(a)//': '//itoa(b)//')'; return; endif
        else
          if(sets_elnn(a)/=4) then; msg=who//': edge surfaces need Q4 elements (element '//itoa(a)//'; 3D edges are deferred)'; return; endif
          if(b<1 .or. b>4) then; msg=who//': Q4 edge id must be 1..4 (element '//itoa(a)//': '//itoa(b)//')'; return; endif
        endif
        if(fmark(b,a)) then; msg=who//': duplicate entry element '//itoa(a)//' id '//itoa(b); return; endif
        fmark(b,a)=.true.
        s%m(2*i-1)=a; s%m(2*i)=b
      enddo
    case(SK_TIE)
      cap=4*c
      allocate(mi(cap), s%off(c+1), s%w(cap))
      j=0
      do i=1,c
        s%off(i)=j+1
        call num(a); if(len(msg)>0) return
        call num(nm); if(len(msg)>0) return
        if(a<1 .or. a>sets_npoin) then; msg=who//': node out of range 1..'//itoa(sets_npoin)//': '//itoa(a); return; endif
        if(nm<1 .or. nm>8) then; msg=who//': tie '//itoa(i)//': number of masters must be 1..8'; return; endif
        do while(j+1+2*nm>cap)
          call grow()
        enddo
        j=j+1; mi(j)=a; s%w(j)=0
        do v=1,nm
          call num(b); if(len(msg)>0) return
          if(b<1 .or. b>sets_npoin) then; msg=who//': node out of range 1..'//itoa(sets_npoin)//': '//itoa(b); return; endif
          if(b==a) then; msg=who//': node '//itoa(a)//' is tied to itself'; return; endif
          j=j+1; mi(j)=b
        enddo
        wsum=0
        do v=1,nm
          call realnum(x); if(len(msg)>0) return
          if(.not.ieee_is_finite(x)) then; msg=who//': non-finite weight'; return; endif
          s%w(j-nm+v)=x; wsum=wsum+x
        enddo
        if(abs(wsum-1._real64)>1.e-6_real64) then; msg=who//': weights of slave '//itoa(a)//' must sum to 1'; return; endif
      enddo
      s%off(c+1)=j+1
      s%m=mi(1:j); s%w=s%w(1:j)
    end select
  contains
    subroutine add(v)
      integer, intent(in) :: v
      if(v<1 .or. v>lim) then
        msg=who//': '//trim(merge('node   ','element',s%kind==SK_NSET))//' out of range 1..'//itoa(lim)//': '//itoa(v); return
      endif
      if(mark(v)) then; msg=who//': duplicate member '//itoa(v); return; endif
      mark(v)=.true.; got=got+1; s%m(got)=v
    end subroutine
    subroutine num(v)
      integer, intent(out) :: v
      v=0
      call next(tok,msg); if(len(msg)>0) return
      if(eof_ .or. .not.is_int(tok)) then
        msg=who//': fewer entries than its count '//itoa(c)//' (line '//itoa(line_)//')'; return
      endif
      read(tok,*) v
    end subroutine
    subroutine realnum(x)
      real(real64), intent(out) :: x
      integer :: ios
      x=0
      call next(tok,msg); if(len(msg)>0) return
      read(tok,*,iostat=ios) x
      if(eof_ .or. ios/=0 .or. verify(trim(tok),'0123456789+-.eEdD')/=0) then
        msg=who//': expected a weight (line '//itoa(line_)//')'; return
      endif
    end subroutine
    subroutine grow()
      integer, allocatable :: t1(:)
      real(real64), allocatable :: t2(:)
      allocate(t1(2*cap), t2(2*cap)); t1(1:cap)=mi; t2(1:cap)=s%w
      call move_alloc(t1,mi); call move_alloc(t2,s%w); cap=2*cap
    end subroutine
  end subroutine

  ! ties over all tie sets: a slave appears once; no slave is a master (no chains)
  subroutine check_ties(msg)
    character(:), allocatable, intent(inout) :: msg
    logical, allocatable :: slave(:)
    integer :: i, r, k
    allocate(slave(max(sets_npoin,1))); slave=.false.
    do i=1,sets_n
      if(sets(i)%kind/=SK_TIE) cycle
      do r=1,sets(i)%n
        k=sets(i)%m(sets(i)%off(r))
        if(slave(k)) then; msg='SETS: tie '//trim(sets(i)%name)//': slave node '//itoa(k)//' is tied twice'; return; endif
        slave(k)=.true.
      enddo
    enddo
    do i=1,sets_n
      if(sets(i)%kind/=SK_TIE) cycle
      do r=1,sets(i)%n
        do k=sets(i)%off(r)+1,sets(i)%off(r+1)-1
          if(slave(sets(i)%m(k))) then
            msg='SETS: tie '//trim(sets(i)%name)//': master node '//itoa(sets(i)%m(k))//' is itself a slave (chained ties)'; return
          endif
        enddo
      enddo
    enddo
  end subroutine

  integer function sets_find(kind,name) result(k)
    integer, intent(in) :: kind
    character(*), intent(in) :: name
    integer :: i
    k=0
    do i=1,sets_n
      if(sets(i)%kind==kind .and. trim(sets(i)%name)==trim(name)) then
        k=i; sets(i)%referenced=.true.; return
      endif
    enddo
  end function

  function sets_first_unreferenced() result(s)
    character(:), allocatable :: s
    integer :: i
    s=''
    do i=1,sets_n
      if(.not.sets(i)%referenced) then
        s=trim(SETS_KIND_NAME(sets(i)%kind))//' '//trim(sets(i)%name); return
      endif
    enddo
  end function

  ! Q4 edge k joins local nodes k and k+1; HSTAR's legacy .loa lists it as (node k+1, node k)
  ! (verified on every LOA regression edge: 26188/26188).
  subroutine sets_edge_nodes(e,k,n1,n2)
    integer, intent(in) :: e,k
    integer, intent(out) :: n1,n2
    n1=sets_elnodes(mod(k,4)+1,e); n2=sets_elnodes(k,e)
  end subroutine

  subroutine sets_face_nodes(e,k,n)
    integer, intent(in) :: e,k
    integer, intent(out) :: n(4)
    n=sets_elnodes(B8_FACE(:,k),e)
  end subroutine

  subroutine sets_receipt(file)
    character(*), intent(in) :: file
    integer :: u,i
    open(newunit=u,file=file,status='replace',action='write')
    write(u,'(a)') 'native-sets v1'
    write(u,'(a,i0,a,i0,a,i0)') 'sets=',sets_n,' npoin=',sets_npoin,' nelem=',sets_nelem
    do i=1,sets_n
      write(u,'(a,1x,a,1x,i0,1x,i0)') trim(SETS_KIND_NAME(sets(i)%kind)),trim(sets(i)%name),sets(i)%sub,sets(i)%n
    enddo
    close(u)
  end subroutine

  ! ---- tokenizer: whitespace separated, "!" starts a comment
  subroutine push(tok)
    character(*), intent(in) :: tok
    back_=tok; has_back_=.true.
  end subroutine

  subroutine next(tok,msg)
    character(*), intent(out) :: tok
    character(:), allocatable, intent(inout) :: msg
    integer :: ios, a, c
    tok=''
    if(has_back_) then; tok=back_; has_back_=.false.; return; endif
    if(eof_) return
    do
      do while(pos_<=SETS_MAX_LINE)
        if(buf_(pos_:pos_)/=' ' .and. buf_(pos_:pos_)/=achar(9)) exit
        pos_=pos_+1
      enddo
      if(pos_<=SETS_MAX_LINE) then
        if(buf_(pos_:pos_)/='!') exit
      endif
      read(u_,'(a)',iostat=ios) buf_
      if(ios/=0) then; eof_=.true.; return; endif
      line_=line_+1; pos_=1
      if(buf_(SETS_MAX_LINE:SETS_MAX_LINE)/=' ') then
        msg='SETS: line '//itoa(line_)//' longer than '//itoa(SETS_MAX_LINE)//' characters'; return
      endif
      c=index(buf_,'!'); if(c>0) buf_(c:)=' '
      do c=1,len_trim(buf_)
        if(buf_(c:c)==achar(13)) buf_(c:c)=' '
      enddo
    enddo
    a=pos_
    do while(pos_<=SETS_MAX_LINE)
      if(buf_(pos_:pos_)==' ' .or. buf_(pos_:pos_)==achar(9)) exit
      pos_=pos_+1
    enddo
    if(pos_-a>len(tok)) then; msg='SETS: token too long (line '//itoa(line_)//')'; return; endif
    tok=buf_(a:pos_-1)
  end subroutine

  logical function is_int(t)
    character(*), intent(in) :: t
    is_int=len_trim(t)>0 .and. len_trim(t)<=9 .and. verify(trim(t),'0123456789')==0
  end function

  logical function valid_name(t)
    character(*), intent(in) :: t
    valid_name=len_trim(t)>0 .and. len_trim(t)<=SETS_LEN_NAME .and. &
      verify(trim(t),'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_')==0
  end function

  integer function ntok(row)
    character(*), intent(in) :: row
    integer :: i
    logical :: inw
    ntok=0; inw=.false.
    do i=1,len_trim(row)
      if(row(i:i)==' ' .or. row(i:i)==achar(9) .or. row(i:i)==achar(13)) then
        inw=.false.
      elseif(.not.inw) then
        inw=.true.; ntok=ntok+1
      endif
    enddo
  end function

  function itoa(i) result(s)
    integer, intent(in) :: i
    character(:), allocatable :: s
    character(16) :: b
    write(b,'(i0)') i
    s=trim(b)
  end function
end module
