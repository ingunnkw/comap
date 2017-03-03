! A new acceptlist module, which emphasizes on simplicity of manipulating the
! accept lists. Therefore each accept list is an int(ndi, nmod, nces).
! This makes the acceptlists that can be represented dependent on the
! ces module and the module module. runs and segments outside the
! current runlist cannot be represented. This makes the module less
! reusable than the previous one. The previous one was problematic because
! it was in C++, which means HK can't maintain it, and because it
! would requrie flattening before synchronizing with others.

module quiet_acceptlist_mod
  use quiet_utils
  use quiet_ces_mod
  use quiet_module_mod
  implicit none

  interface is_accepted_scan_any
     module procedure is_accepted
  end interface

  type acceptlist
     integer(i4b) :: ndi, nmod, nces
     integer(i4b), dimension(:,:,:), allocatable :: status
  end type acceptlist

  logical(lgt), private :: initialized = .false.
  integer(i4b), parameter :: REJECTED_NONE = 0, REJECTED_ALIST = 1, REJECTED_NOREASON = 2, REJECTED_TARGET = 3, REJECTED_SPLIT = 4

contains

  subroutine allocate_acceptlist(alist)
    implicit none
    type(acceptlist)            :: alist
    call deallocate_acceptlist(alist)
    alist%ndi  = get_num_diodes()
    alist%nmod = get_num_modules()
    alist%nces = get_num_ces()
    allocate(alist%status(0:alist%ndi-1,0:alist%nmod-1,alist%nces))
  end subroutine allocate_acceptlist

  subroutine deallocate_acceptlist(alist)
    implicit none
    type(acceptlist) :: alist
    if (allocated(alist%status)) deallocate(alist%status)
  end subroutine deallocate_acceptlist

  subroutine get_accepted_ceses(alist, cnums)
    implicit none
    type(acceptlist)                        :: alist
    integer(i4b), dimension(:), allocatable :: cnums
    integer(i4b)                            :: i, n
    n = 0
    do i = 1, size(alist%status,3)
       if(any(alist%status(:,:,i) == REJECTED_NONE)) n = n+1
    end do
    allocate(cnums(n))
    n = 0
    do i = 1, size(alist%status,3)
       if(any(alist%status(:,:,i) == REJECTED_NONE)) then
          n = n+1
          cnums(n) = i
       end if
    end do
  end subroutine

  ! Like above, but is seen as accepted if any of the alists accept it
  subroutine get_accepted_ceses_multi(alists, cnums)
    implicit none
    type(acceptlist)                        :: alists(:)
    integer(i4b), dimension(:), allocatable :: cnums
    integer(i4b)                            :: i, j, n
    n = 0
    do i = 1, size(alists(1)%status,3)
       do j = 1, size(alists)
          if(any(alists(j)%status(:,:,i) == REJECTED_NONE)) then
             n = n+1
             exit
          end if
       end do
    end do
    allocate(cnums(n))
    n = 0
    do i = 1, size(alists(1)%status,3)
       do j = 1, size(alists)
          if(any(alists(j)%status(:,:,i) == REJECTED_NONE)) then
             n = n+1
             cnums(n) = i
             exit
          end if
       end do
    end do
  end subroutine

  function is_accepted(alist, cid, mod, di) result(res)
    implicit none
    type(acceptlist)       :: alist
    integer(i4b)           :: cid, cnum
    integer(i4b), optional :: mod, di
    logical(lgt)           :: res
    cnum = lookup_ces(cid)
    if (present(di)) then
       res  = alist%status(di,mod,cnum) == REJECTED_NONE
    else if (present(mod)) then
       res  = any(alist%status(:,mod,cnum) == REJECTED_NONE)
    else 
       res  = any(alist%status(:,:,cnum) == REJECTED_NONE)
    end if
  end function is_accepted

  subroutine accept(alist, cid, mod, di)
    implicit none
    type(acceptlist)  :: alist
    integer(i4b), optional :: cid, mod, di
    integer(i4b)      :: cnum
    if (.not. present(cid)) then
       alist%status = REJECTED_NONE
       return
    end if
    cnum = lookup_ces(cid)
    if (cnum < 1 .or. cnum > alist%nces) return
    if (present(di)) then
       alist%status(di,mod,cnum) = REJECTED_NONE
    else if (present(mod)) then
       alist%status(:,mod,cnum)  = REJECTED_NONE
    else 
       alist%status(:,:,cnum)    = REJECTED_NONE
    end if
  end subroutine accept

  subroutine reject(alist, cid, mod, di, reason)
    implicit none
    type(acceptlist)  :: alist
    integer(i4b), optional :: cid, mod, di, reason
    integer(i4b)      :: r, cnum
    cnum = lookup_ces(cid)
    if (cnum < 1 .or. cnum > alist%nces) return
    r = REJECTED_NOREASON; if(present(reason)) r = reason
    if (present(di)) then
       alist%status(di,mod,cnum) = r
    else if (present(mod)) then
       alist%status(:,mod,cnum)  = r
    else if (present(cid)) then
       alist%status(:,:,cnum)    = r
    else
       alist%status              = r
    end if
  end subroutine reject

  subroutine initialize_accept_list(filename, alist, default)
    implicit none
    character(len=*)    :: filename
    type(acceptlist)    :: alist
    character(len=1024) :: line
    integer(i4b) :: i, unit, cid, mod, di, cnum, status, defval
    integer(i4b), optional :: default
    integer(i4b), dimension(:,:), allocatable :: rejects
    ! Is this a dummy wildcard file? If so, don't add anything, but
    ! set the default value to be accept instead of reject.
    defval = REJECTED_ALIST; if(present(default)) defval = default
    call allocate_acceptlist(alist)
    if(filename == "*") then
       call accept(alist)
       return
    end if
    alist%status = defval
    ! Actually read data
    unit = getlun()
    open(unit,file=filename,action="read")
    do
       read(unit,fmt="(a)",end=2) line
       if(line(1:1) == '#') cycle
       read(line,*) cid, mod, di, status
       if (status == 1) then
          call accept(alist, cid, mod, di)
       else
          call reject(alist, cid, mod, di, REJECTED_ALIST)
       end if
    end do
2   close(unit)
  end subroutine initialize_accept_list

  subroutine append_accept_list(filename, cid, alist)
    implicit none
    character(len=*)    :: filename
    integer(i4b)        :: cid
    type(acceptlist)    :: alist

    character(len=4096) :: line, diode
    integer(i4b) :: i, j, unit, cnum

    cnum = lookup_ces(cid)
    if (cnum < 1 .or. cnum > alist%nces) return

    ! Write line to file
    unit = getlun()
    open(unit,file=filename,action="write",position="append",recl=4096)
    do i = 0, alist%nmod-1
       do j = 0, alist%ndi-1
          if (alist%status(j,i,cnum) == REJECTED_NONE) then
             write(*,fmt='(i7,i5,i3,i3)') cid, i, j, 1
          else
             write(*,fmt='(i7,i5,i3,i3)') cid, i, j, 0
          end if
       end do
    end do
    close(unit)

  end subroutine append_accept_list

end module quiet_acceptlist_mod
