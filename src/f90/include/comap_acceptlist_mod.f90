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

  type rejectfreqs
     integer(i4b) :: numranges       ! Number of rejected frequency intervals
     integer(i4b) :: numrf           ! Number of rejected frequencies
     integer(i4b) :: numaf           ! Number of accepted frequencies
     integer(i4b), dimension(:,:), allocatable :: ranges !(3(beg,end,stat),numranges)
  end type rejectfreqs

  type acceptlist
     integer(i4b) :: nfreq, ndet, nscan
     integer(i4b),      dimension(:),   allocatable :: scans
     integer(i4b),      dimension(:,:), allocatable :: status  !(ndet, nscan)
     type(rejectfreqs), dimension(:,:), allocatable :: rf      !(ndet, nscan)
  end type acceptlist

  type afreq
     integer(i4b) :: n
     integer(i4b), allocatable, dimension(:) :: rejected
  end type afreq

  type adet
     type(afreq), allocatable, dimension(:,:) :: adet_sb
  end type adet

  logical(lgt), private :: initialized = .false.
  integer(i4b), parameter :: REJECTED_NONE = 0, REJECTED_ALIST = 1, REJECTED_NOREASON = 2, REJECTED_TARGET = 3, REJECTED_SPLIT = 4

contains

  subroutine allocate_acceptlist(alist)
    implicit none
    type(acceptlist)            :: alist
    call deallocate_acceptlist(alist)
    alist%nfreq = get_num_freqs()
    alist%ndet  = get_num_dets()
    alist%nscan = get_num_scans()
    allocate(alist%status( 1:alist%ndet, 1:alist%nscan))
    allocate(alist%rf(     1:alist%ndet, 1:alist%nscan))
  end subroutine allocate_acceptlist

  subroutine deallocate_acceptlist(alist)
    implicit none
    type(acceptlist) :: alist
    if (allocated(alist%status)) deallocate(alist%status)
    if (allocated(alist%rf))     deallocate(alist%rf)
  end subroutine deallocate_acceptlist

  subroutine initialize_accept_list(filename, alist, default)
    implicit none
    character(len=*)    :: filename
    type(acceptlist)    :: alist
    character(len=8)    :: sid
    character(len=1024) :: line
    integer(i4b)        :: i, j, k, l, unit, numscans, numdet, numsb, numrej, rfreq

    !call allocate_acceptlist(alist)
    
    unit = getlun()
    open(unit, file=filename, action="read")
    ! Read number of scans
    do
       read(unit, fmt="(a)", end=2) line
       if (line(1:1) == '#') then
          cycle
       else
          read(line,*) numscans
          !alist%nscan = numscans
          allocate(alist%scans(numscans))
          exit
       end if
    end do
    ! Read rejections per scan
    do i = 1, numscans
       read(unit, fmt="a", end=2)
       read(line,*) sid, numdet
       do j = 1, numdet
          read(unit, fmt="a", end=2)
          read(line,*) det, numsb
          do k = 1, numsb
             read(unit, fmt="a", end=2)
             read(line,*) sb, numrej
             allocate(adet_sb(det,sb)%rejected(numrej))
             do l = 1, numrej
                read(unit, fmt="a", end=2)
                read(line,*) rfreq
                adet_sb(det,sb)%rejected(l) = rfreq
             end do
          end do
       end do
    end do

  end subroutine initialize_accept_list


!  subroutine initialize_accept_list(filename, alist, default)
!     implicit none
!     character(len=*)       :: filename
!     type(acceptlist)       :: alist
!     character(len=9)       :: sid
!     character(len=1024)    :: line
!     integer(i4b)           :: i, j, k, unit, det, snum, status, defval
!     integer(i4b)           :: numscans, numdets, numrej, rfreqs 
!     integer(i4b), optional :: default
!     integer(i4b), dimension(1:2*alist%nfreq) :: rejects
!     ! Is this a dummy wildcard file? If so, don't add anything, but
!     ! set the default value to be accept instead of reject.
!     !defval = REJECTED_ALIST; if(present(default)) defval = default
!     !call allocate_acceptlist(alist)
!     !if(filename == "*") then
!     !   call accept(alist)
!     !   return
!     !end if
!     alist%status = 0 ! Everything accepted as default
!     ! Read number of scans
!     unit = getlun()
!     open(unit,file=filename,action="read")
!     do
!        read(unit,fmt="(a)",end=2) line
!        if(line(1:1) == '#') then 
!           cycle
!        else
!           read(line,*) numscans
!           exit
!        end if
!     end do
!     ! Read rejections per scan
!     do i = 1, numscans
!        read(unit,*) sid, numdets
!        snum = lookup_scan(sid)
!        do j = 1, numdets
!           ! Read number of rejected frequency ranges
!           read(unit,fmt="(a)",end=2) line
!           read(line,*) det, numrej
!           rfreqs=0
!           if (numrej /= 0) then
!              allocate(alist%rf(det,snum)%ranges(2,numrej))
!              alist%rf(det,snum)%numranges = numrej
!              read(line,*) det, numrej, rejects
!              do k = 1, numrej
!                 alist%rf(det,snum)%ranges(1:2,k) = rejects(2*k-1:2*k) 
!              end do
!              ! Count number of rejected frequencies
!              do k = 1, numrej
!                 rfreqs = rfreqs + alist%rf(det,snum)%ranges(2,k)-alist%rf(det,snum)%ranges(1,k)+1
!              end do
!           end if
!           alist%rf(det,snum)%numrf = rfreqs
!           alist%status = rfreqs
!           alist%rf(det,snum)%ranges(3,k) = REJECTED_ALIST 
!        end do
!     end do
! 2   close(unit)
!   end subroutine initialize_accept_list

  subroutine get_accepted_scans(alist, snums)
    implicit none
    type(acceptlist)                        :: alist
    integer(i4b), dimension(:), allocatable :: snums
    integer(i4b)                            :: i, n
    n = 0
    do i = 1, size(alist%status,2)
       if(any(alist%status(:,i) /= alist%nfreq)) n = n+1
    end do
    allocate(snums(n))
    n = 0
    do i = 1, size(alist%status,2)
       if(any(alist%status(:,i) /= alist%nfreq)) then
          n = n+1
          snums(n) = i
       end if
    end do
  end subroutine

  subroutine get_amatrix_per_scan(alist, sid, amat)
    implicit none
    type(acceptlist)                          :: alist
    character(len=9)                          :: sid
    integer(i4b), dimension(:,:), allocatable :: amat
    integer(i4b)                              :: det
    integer(i4b), dimension(:), allocatable   :: avec

    ! allocate(amat(alist%nfreq,alist%ndet))
    amat = 0;
    allocate(avec(alist%nfreq))
    do det = 1, alist%ndet
       call get_avector_per_detscan(alist, det, sid, avec)
       amat(:,det) = avec
    end do

  end subroutine get_amatrix_per_scan

  subroutine get_avector_per_detscan(alist, det, sid, avec)
    implicit none
    type(acceptlist)                          :: alist
    integer(i4b)                              :: det
    character(len=9)                          :: sid
    integer(i4b), dimension(:), allocatable   :: avec
    integer(i4b)                              :: i, snum
    
    !allocate(avec(alist%nfreq))
    avec = 0;
    snum = lookup_scan(sid)    
    do i = 1, alist%rf(det,snum)%numranges
       avec(alist%rf(det,snum)%ranges(1:2,i)) = avec(alist%rf(det,snum)%ranges(3,i))
    end do
  end subroutine get_avector_per_detscan

  function freq_is_accepted(alist, sid, det, freq) result(res)
    implicit none
    type(acceptlist)       :: alist
    integer(i4b)           :: det, freq
    character(len=9)       :: sid
    logical(lgt)           :: res
    integer(i4b)           :: snum
    snum = lookup_scan(sid)
  end function freq_is_accepted




  function is_accepted(alist, sid, det, freq) result(res)
    implicit none
    type(acceptlist)       :: alist
    character(len=9)       :: sid
    integer(i4b), optional :: det, freq
    logical(lgt)           :: res
    integer(i4b)           :: snum

    snum = lookup_scan(sid)

    if (present(freq)) then
       !res  = alist%status(freq,det,snum) == REJECTED_NONE
    else
       

    end if

  end function is_accepted

end module comap_acceptlist_mod


!!$  subroutine accept(alist, sid, det, freqs)
!!$    implicit none
!!$    type(acceptlist)            :: alist
!!$    integer(i4b),      optional :: sid, det
!!$    type(rejectfreqs), optional :: freqs
!!$    integer(i4b)                :: snum
!!$    if (.not. present(sid)) then
!!$       alist%status = REJECTED_NONE
!!$       return
!!$    end if
!!$    snum = lookup_scan(sid)
!!$    if (snum < 1 .or. snum > alist%nscan) return
!!$    if (.not. present(det)) then
!!$       alist%status(:,snum) = REJECTED_NONE
!!$    else
!!$       alist%status(det,snum) = REJECTED_NONE
!!$    end if
!!$    if (present(freqs)) then
!!$
!!$   alist%status(:,det,snum)    = REJECTED_NONE
!!$    else
!!$       alist%status(:,:,snum)      = REJECTED_NONE
!!$    end if
!!$  end subroutine accept
!!$
!!$  subroutine reject(alist, sid, freq, det, reason)
!!$    implicit none
!!$    type(acceptlist)       :: alist
!!$    integer(i4b), optional :: sid, freq, det, reason
!!$    integer(i4b)           :: r, snum
!!$    snum = lookup_scan(sid)
!!$    if (snum < 1 .or. snum > alist%nscan) return
!!$    r = REJECTED_NOREASON; if(present(reason)) r = reason
!!$    if (present(freq) .and. present(det)) then
!!$       alist%status(freq,det,snum) = r
!!$    else if (present(freq)) then
!!$       alist%status(freq,:,snum)   = r
!!$    else if (present(det)) then
!!$       alist%status(:,det,snum)    = r
!!$    else if (present(sid)) then
!!$       alist%status(:,:,snum)      = r
!!$    else
!!$       alist%status                = r
!!$    end if
!!$  end subroutine reject
!!$
!!$  ! Like above, but is seen as accepted if any of the alists accept it
!!$  subroutine get_accepted_scans_multi(alists, snums)
!!$    implicit none
!!$    type(acceptlist)                        :: alists(:)
!!$    integer(i4b), dimension(:), allocatable :: snums
!!$    integer(i4b)                            :: i, j, n
!!$    n = 0
!!$    do i = 1, size(alists(1)%status,3)
!!$       do j = 1, size(alists)
!!$          if(any(alists(j)%status(:,:,i) == REJECTED_NONE)) then
!!$             n = n+1
!!$             exit
!!$          end if
!!$       end do
!!$    end do
!!$    allocate(snums(n))
!!$    n = 0
!!$    do i = 1, size(alists(1)%status,3)
!!$       do j = 1, size(alists)
!!$          if(any(alists(j)%status(:,:,i) == REJECTED_NONE)) then
!!$             n = n+1
!!$             snums(n) = i
!!$             exit
!!$          end if
!!$       end do
!!$    end do
!!$  end subroutine
