! A new acceptlist module, which emphasizes on simplicity of manipulating the
! accept lists. Therefore each accept list is an int(ndet, nces).
! This makes the acceptlists that can be represented dependent on the
! sacn module and the detector module. runs and segments outside the
! current runlist cannot be represented. 
! This has been copied from quiet_acceptlist_mod.

module comap_acceptlist_mod
  use quiet_utils
  use comap_scan_mod
  use comap_detector_mod
  use comap_frequency_mod
  implicit none

  interface is_accepted_scan_any
     module procedure is_accepted
  end interface

  type acceptlist
     integer(i4b) :: nfreq, ndet, nscan
     integer(i4b), dimension(:,:,:), allocatable :: status !(nfreq, ndet, nscan)
  end type acceptlist

  logical(lgt), private :: initialized = .false.
  integer(i4b), parameter :: REJECTED_NONE = 0, REJECTED_ALIST = 1, REJECTED_NOREASON = 2, REJECTED_TARGET = 3, REJECTED_SPLIT = 4

contains

  subroutine allocate_acceptlist(alist)
    implicit none
    type(acceptlist)            :: alist
    call deallocate_acceptlist(alist)
    alist%nfreq = get_num_freqs()
    alist%ndet  = get_num_detectors()
    alist%nscan = get_num_scans()
    allocate(alist%status(1:alist%nfreq,1:alist%ndet,1:alist%nscan))
  end subroutine allocate_acceptlist

  subroutine deallocate_acceptlist(alist)
    implicit none
    type(acceptlist) :: alist
    if (allocated(alist%status)) deallocate(alist%status)
  end subroutine deallocate_acceptlist

  subroutine get_accepted_scans(alist, snums)
    implicit none
    type(acceptlist)                        :: alist
    integer(i4b), dimension(:), allocatable :: snums
    integer(i4b)                            :: i, n
    n = 0
    do i = 1, size(alist%status,3)
       if(any(alist%status(:,:,i) == REJECTED_NONE)) n = n+1
    end do
    allocate(snums(n))
    n = 0
    do i = 1, size(alist%status,3)
       if(any(alist%status(:,:,i) == REJECTED_NONE)) then
          n = n+1
          snums(n) = i
       end if
    end do
  end subroutine

  ! Like above, but is seen as accepted if any of the alists accept it
  subroutine get_accepted_scans_multi(alists, snums)
    implicit none
    type(acceptlist)                        :: alists(:)
    integer(i4b), dimension(:), allocatable :: snums
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
    allocate(snums(n))
    n = 0
    do i = 1, size(alists(1)%status,3)
       do j = 1, size(alists)
          if(any(alists(j)%status(:,:,i) == REJECTED_NONE)) then
             n = n+1
             snums(n) = i
             exit
          end if
       end do
    end do
  end subroutine

  function is_accepted(alist, sid, freq, det) result(res)
    implicit none
    type(acceptlist)       :: alist
    integer(i4b)           :: sid, snum
    integer(i4b), optional :: det, freq
    logical(lgt)           :: res
    snum = lookup_scan(sid)
    if (present(freq) .and. present(det)) then
       res  = alist%status(freq,det,snum) == REJECTED_NONE
    else if (present(freq)) then
       res  = any(alist%status(freq,:,snum)   == REJECTED_NONE)
    else if (present(det)) then
       res  = any(alist%status(:,det,snum)    == REJECTED_NONE)
    else 
       res  = any(alist%status(:,:,snum)  == REJECTED_NONE)
    end if
  end function is_accepted

  subroutine accept(alist, sid, freq, det)
    implicit none
    type(acceptlist)       :: alist
    integer(i4b), optional :: sid, freq, det
    integer(i4b)           :: snum
    if (.not. present(sid)) then
       alist%status = REJECTED_NONE
       return
    end if
    snum = lookup_scan(sid)
    if (snum < 1 .or. snum > alist%nscan) return
    if (present(freq) .and. present(det)) then
       alist%status(freq,det,snum) = REJECTED_NONE
    else if (present(freq)) then
       alist%status(freq,:,snum)   = REJECTED_NONE
    else if (present(det)) then
       alist%status(:,det,snum)    = REJECTED_NONE
    else
       alist%status(:,:,snum)      = REJECTED_NONE
    end if
  end subroutine accept

  subroutine reject(alist, sid, freq, det, reason)
    implicit none
    type(acceptlist)       :: alist
    integer(i4b), optional :: sid, freq, det, reason
    integer(i4b)           :: r, snum
    snum = lookup_scan(sid)
    if (snum < 1 .or. snum > alist%nscan) return
    r = REJECTED_NOREASON; if(present(reason)) r = reason
    if (present(freq) .and. present(det)) then
       alist%status(freq,det,snum) = r
    else if (present(freq)) then
       alist%status(freq,:,snum)   = r
    else if (present(det)) then
       alist%status(:,det,snum)    = r
    else if (present(sid)) then
       alist%status(:,:,snum)      = r
    else
       alist%status                = r
    end if
  end subroutine reject

  subroutine initialize_accept_list(filename, alist, default)
    implicit none
    character(len=*)       :: filename
    type(acceptlist)       :: alist
    character(len=1024)    :: line
    integer(i4b)           :: i, j, unit, sid, det, snum, status, defval
    integer(i4b)           :: numscans, numdets, numrej, freq 
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
       if(line(1:1) == '#') then 
          cycle
       else
          read(line,*) numscans
          exit
       end if
    end do

    do i = 1, numscans
       read(unit,*) sid, numdets
       do j =1, numdets
          read(unit,fmt="(a)",end=2) line
          ! Then chunk line
          ! read(line,*) det, numrej, ++
          if (numrej == 0) then
             call accept(alist, sid, det)
          else
             call reject(alist, sid, det, REJECTED_ALIST) ! add which freqs
          end if
        end do
     end do
2    close(unit)
   end subroutine initialize_accept_list

end module comap_acceptlist_mod
