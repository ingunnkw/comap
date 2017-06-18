module comap_scan_mod
  use healpix_types
  use quiet_utils
  use sort_utils
  implicit none

  ! Information about a comap_scan_info. Only the information that can be
  ! extracted directly from the level2-runlist. Category and object
  ! are indices into name arrays.
  type comap_scan_info
    integer(i4b)       :: sid
    real(dp)           :: mjd(2), az, el, dk, lat, lon, time_of_day
    character(len=512) :: l2file, l3file, object, scanmode
    character(len=512), allocatable, dimension(:) :: l1files
  end type comap_scan_info

  type comap_runlist
     integer(i4b)                                     :: n
     type(comap_scan_info), dimension(:), allocatable :: scans
     integer(i4b),          dimension(:), allocatable :: sidmap
  end type comap_runlist

  integer(i4b), parameter :: PROP_TOD = 1, PROP_AZ = 2, PROP_EL = 3, &
    & PROP_DECK = 4, PROP_DUR = 5, PROP_MJD = 6

  type(comap_runlist)          :: scan_db
  integer(i4b),    allocatable :: sid_list(:), sid_sort(:)
  logical(lgt),        private :: initialized = .false.
  character(len=512),  private :: l1dir, l2dir, l3dir

contains

  ! Public interface: These are meant to be called from outside.

  subroutine initialize_scan_mod(parfile)
    implicit none
    character(len=*)       :: parfile
    character(len=512)     :: runlist
    if(initialized) return
    call get_parameter(0, parfile, 'RUNLIST',     par_string=runlist)
    call get_parameter(0, parfile, 'LEVEL1_DIR',  par_string=l1dir)
    call get_parameter(0, parfile, 'LEVEL2_DIR',  par_string=l2dir)
    call get_parameter(0, parfile, 'LEVEL3_DIR',  par_string=l3dir)

    call read_runlist(runlist, l1dir, l2dir, l3dir, scan_db)
    allocate(sid_list(scan_db%n), sid_sort(scan_db%n))
    sid_list = scan_db%scans%sid
    sid_sort = sid_list
    call quicksort_int(sid_sort)
    initialized = .true.
  end subroutine

  function get_num_scans() result(res)
    implicit none
    integer(i4b) :: res
    res = size(scan_db%scans)
  end function

  function lookup_scan(sid) result(res)
    implicit none
    integer(i4b) :: sid, res
    res = 0
    if(sid < 0 .or. sid > size(scan_db%sidmap)) return
    res = scan_db%sidmap(sid)
  end function

  subroutine get_scan_info(ind, scan)
    implicit none
    type(comap_scan_info) :: scan
    integer(i4b)    :: ind
    call copy_scan_info(scan_db%scans(ind), scan)
  end subroutine

  subroutine get_l2_param(sid, index, mjd, az, el, dk, time, target)
    implicit none
    integer(i4b),                   intent(in),  optional :: sid, index
    real(dp),                       intent(out), optional :: az, el, dk, time
    real(dp),         dimension(2), intent(out), optional :: mjd
    character(len=*),               intent(out), optional :: target
    type(comap_scan_info) :: info

    if(present(index)) then
       call get_scan_info(index, info)
    else
       if(.not. present(sid)) then
          write(stderr,*) "Either index or sid must be given to get_l2_param!"
          stop
       end if
       call assert(sid > 0, "Non-positive sid in get_l2_param!")
       call get_scan_info(scan_db%sidmap(sid), info)
    end if

    if (present(az))   then; az   = info%az; end if
    if (present(el))   then; el   = info%el; end if
    if (present(dk))   then; dk   = info%dk; end if
    if (present(time)) then; time = info%time_of_day; end if
    if (present(mjd))  then; mjd  = info%mjd; end if
    if (present(target)) then; target = info%object; end if
    call free_scan_info(info)
  end subroutine get_l2_param

  !-----------------------------------------!
  !--------- The helper routines -----------!
  !-----------------------------------------!

  subroutine read_runlist(file, l1dir, l2dir, l3dir, runlist)
    implicit none
    character(len=*)   :: file
    character(len=512) :: name, l1dir, l2dir, l3dir, line
    integer(i4b)       :: unit, nobj, nscan, nfile, sid, i, j, k, n, cnum, cmax
    real(dp)           :: mjd(2)
    type(comap_runlist):: runlist
    type(comap_scan_info)    :: scan
    n = count_num_scan(file)
    call free_runlist(runlist)
    allocate(runlist%scans(n))
    cnum = 0
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) nobj
    do i = 1, nobj
       read(unit,*) name, nscan
       do j = 1, nscan
          cnum = cnum+1
          read(unit,*) scan%sid, scan%mjd, nfile, scan%az, scan%el, scan%dk, scan%lon, scan%lat
          scan%object      = name
          scan%time_of_day = mjd2chile_time(mean(scan%mjd))
          allocate(scan%l1files(nfile))
          do k = 1, nfile
             read(unit,fmt="(a)") line
             scan%l1files(k) = trim(l1dir) // "/" // trim(get_token(line, " ", 3))
          end do
          scan%l2file = trim(l2dir) // "/" // trim(name) // "/" // trim(name) // "_" // &
           & trim(itoa(scan%sid)) // ".h5"
          scan%l3file = trim(l3dir) // "/" // trim(name) // "/" // trim(name) // "_" // &
           & trim(itoa(scan%sid)) // ".h5"
          call copy_scan_info(scan, runlist%scans(cnum))
          call free_scan_info(scan)
       end do
    end do
    close(unit)
    ! Set up the sid mapping
    cmax = maxval(runlist%scans%sid)
    allocate(runlist%sidmap(cmax))
    runlist%sidmap = 0
    do i = 1, size(runlist%scans)
       runlist%sidmap(runlist%scans(i)%sid) = i
    end do
    runlist%n = size(runlist%scans)
  end subroutine read_runlist

  function count_num_scan(file) result(n)
    implicit none
    character(len=*)   :: file
    character(len=512) :: name
    integer(i4b)       :: unit, nobj, nscan, nfile, sid, i, j, k, n
    real(dp)           :: mjd(2)
    n = 0
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) nobj
    do i = 1, nobj
       read(unit,*) name, nscan
       n = n+nscan
       do j = 1, nscan
          read(unit,*) sid, mjd, nfile
          do k = 1, nfile
             read(unit,*)
          end do
       end do
    end do
    close(unit)
  end function

  subroutine free_scan_info(a)
    implicit none
    type(comap_scan_info) :: a
    if(allocated(a%l1files)) deallocate(a%l1files)
  end subroutine

  subroutine copy_scan_info(a, b)
    implicit none
    type(comap_scan_info) :: a, b
    call free_scan_info(b)
    allocate(b%l1files(size(a%l1files)))
    b = a
  end subroutine

  subroutine free_runlist(runlist)
    implicit none
    type(comap_runlist) :: runlist
    integer(i4b) :: i
    if(allocated(runlist%scans)) then
       do i = 1, size(runlist%scans)
          call free_scan_info(runlist%scans(i))
       end do
       deallocate(runlist%scans)
    end if
    if(allocated(runlist%sidmap)) deallocate(runlist%sidmap)
  end subroutine

  subroutine copy_runlist(a, b)
    implicit none
    type(comap_runlist) :: a
    type(comap_runlist) :: b
    integer(i4b)        :: i
    call free_runlist(b)
    allocate(b%scans(size(a%scans)))
    allocate(b%sidmap(size(a%sidmap)))
    do i = 1, size(a%scans)
       call copy_scan_info(a%scans(i),b%scans(i))
    end do
    b%sidmap = a%sidmap
  end subroutine

  function mjd2chile_time(mjd) result(res)
    implicit none
    real(dp) :: mjd, res
    res = modulo((mjd-int(mjd))*24-3,24d0)
  end function

end module
