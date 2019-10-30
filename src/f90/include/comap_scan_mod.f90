module comap_scan_mod
  use healpix_types
  use quiet_utils
  use sort_utils
  implicit none

  ! Information about a comap_scan_info. Only the information that can be
  ! extracted directly from the level2-runlist. Category and object
  ! are indices into name arrays.
  type comap_subscan
    integer(i4b)       :: id, feature
    real(dp)           :: mjd(2), az, el, dk, lat, lon, time_of_day
    character(len=512) :: l2file, l3file, scanmode
  end type comap_subscan

  type comap_scan_info
    integer(i4b)         :: id, nsub, objectnum
    real(dp)             :: mjd(2)
    character(len=512)   :: l1file, object
    type(comap_subscan), allocatable, dimension(:) :: ss
  end type comap_scan_info

  type comap_runlist
     integer(i4b)                                     :: n, nsub
     type(comap_scan_info), dimension(:), allocatable :: scans
  end type comap_runlist

  integer(i4b), parameter :: PROP_TOD = 1, PROP_AZ = 2, PROP_EL = 3, &
    & PROP_DECK = 4, PROP_DUR = 5, PROP_MJD = 6

  type(comap_runlist)          :: scan_db
  logical(lgt),        private :: initialized = .false.
  character(len=512),  private :: l1dir, l2dir, l3dir

contains

  ! Public interface: These are meant to be called from outside.

  subroutine initialize_scan_mod(parfile, object)
    implicit none
    character(len=*)       :: parfile
    character(len=512)     :: runlist
    character(len=512), optional :: object
    if(initialized) return
    call get_parameter(0, parfile, 'RUNLIST',     par_string=runlist)
    call get_parameter(0, parfile, 'LEVEL1_DIR',  par_string=l1dir)
    call get_parameter(0, parfile, 'LEVEL2_DIR',  par_string=l2dir)
    call get_parameter(0, parfile, 'LEVEL3_DIR',  par_string=l3dir)

    call read_runlist(runlist, l1dir, l2dir, l3dir, scan_db, object)
    initialized = .true.
  end subroutine

  function get_num_scans() result(res)
    implicit none
    integer(i4b) :: res
    res = scan_db%n
  end function

  function get_num_subscans() result(res)
    implicit none
    integer(i4b) :: res
    res = scan_db%nsub
  end function

!!$  function lookup_scan(sid) result(res)
!!$    implicit none
!!$    character(len=9) :: sid
!!$    integer(i4b) :: res, i
!!$    res = 0
!!$    do i = 1, scan_db%n
!!$       if (scan_db%scans(i)%id == sid) then
!!$          res = i
!!$          exit
!!$       end if
!!$    end do
!!$  end function

  subroutine get_scan_info(ind, scan)
    implicit none
    type(comap_scan_info) :: scan
    integer(i4b)    :: ind
    call copy_scan_info(scan_db%scans(ind), scan)
  end subroutine

!!$  subroutine get_l2_param(sid, index, mjd, az, el, dk, time, target)
!!$    implicit none
!!$    integer(i4b),                   intent(in),  optional :: sid, index
!!$    real(dp),                       intent(out), optional :: az, el, dk, time
!!$    real(dp),         dimension(2), intent(out), optional :: mjd
!!$    character(len=*),               intent(out), optional :: target
!!$    type(comap_scan_info) :: info
!!$
!!$    if(present(index)) then
!!$       call get_scan_info(index, info)
!!$    else
!!$       if(.not. present(sid)) then
!!$          write(stderr,*) "Either index or sid must be given to get_l2_param!"
!!$          stop
!!$       end if
!!$       call assert(sid > 0, "Non-positive sid in get_l2_param!")
!!$       call get_scan_info(scan_db%sidmap(sid), info)
!!$    end if
!!$
!!$    if (present(az))   then; az   = info%az; end if
!!$    if (present(el))   then; el   = info%el; end if
!!$    if (present(dk))   then; dk   = info%dk; end if
!!$    if (present(time)) then; time = info%time_of_day; end if
!!$    if (present(mjd))  then; mjd  = info%mjd; end if
!!$    if (present(target)) then; target = info%object; end if
!!$    call free_scan_info(info)
!!$  end subroutine get_l2_param

  !-----------------------------------------!
  !--------- The helper routines -----------!
  !-----------------------------------------!

  subroutine read_runlist(file, l1dir, l2dir, l3dir, runlist, object)
    implicit none
    character(len=*)   :: file
    character(len=512) :: name, l1dir, l2dir, l3dir, line, l1file, objlist(29)
    character(len=512), optional :: object
    character(len=9)   :: subsid
    integer(i4b)       :: unit, nobj, nscan, nfile, sid, i, j, k, n, cnum, cmax, nsub, feature, a, b, c, d
    real(dp)           :: mjd(2), mjd_file(2)
    type(comap_runlist):: runlist
    type(comap_scan_info)    :: scan
    n = count_num_scan(file)
    call free_runlist(runlist)
    allocate(runlist%scans(n))
    objlist = (/'jupiter', 'venus', 'TauA', 'shela', 'hetdex', 'patch1', &
         'patch2', 'co1', 'co2', 'co3', 'co4', 'co5', 'co6', 'co7', 'mars', &
         'fg1', 'fg2', 'fg3', 'fg4', 'fg5', 'fg6', 'fg7', 'ambient_load', &
         'ground_scan', 'stationary', 'sky_dip', 'CasA', 'CygA', 'other' /)
    cnum = 0
    runlist%nsub = 0
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) nobj
    do i = 1, nobj
       read(unit,*) name, nscan
       do j = 1, nscan
          read(unit,'(a)') line
          read(line,*) scan%id, scan%mjd, scan%nsub
          scan%l1file = trim(l1dir) // "/" // trim(get_token(line, " ", 5))
          scan%object = name
          scan%objectnum = findloc(objlist, trim(name), 1)
          allocate(scan%ss(scan%nsub))
          cnum         = cnum+1
          runlist%nsub = runlist%nsub + scan%nsub
          do k = 1, scan%nsub
             read(unit,*) scan%ss(k)%id, scan%ss(k)%mjd, feature, &
                  & scan%ss(k)%az, scan%ss(k)%el, scan%ss(k)%dk, a, b, c, d
             call int2string(scan%ss(k)%id, subsid)
             if (present(object)) then
                if (trim(object) /= trim(name)) cycle
             end if
             scan%ss(k)%feature = feature
             if (feature == 16) then
                scan%ss(k)%scanmode = 'circ'
             else if (feature == 32) then
                scan%ss(k)%scanmode = 'ces'
             else if (feature == 512) then
                scan%ss(k)%scanmode = 'raster'
             else if (feature == 8192) then
                scan%ss(k)%scanmode = 'tsys'
             else if (feature == 32768) then
                scan%ss(k)%scanmode = 'liss'
             else
                scan%ss(k)%scanmode = 'none'
             end if
             !scan%ss(k)%id          = k
             scan%ss(k)%time_of_day = mjd2chile_time(0.5d0*(scan%ss(k)%mjd(1)+scan%ss(k)%mjd(2)))
             scan%ss(k)%l2file = trim(l2dir) // "/" // trim(name) // "/" // trim(name) // "_" // &
                  & trim(subsid) // ".h5"
             scan%ss(k)%l3file = trim(l3dir) // "/" // trim(name) // "/" // trim(name) // "_" // &
                  & trim(subsid) // ".h5"
          end do
          call copy_scan_info(scan, runlist%scans(cnum))
          call free_scan_info(scan)
       end do
    end do
    close(unit)
    ! Set up the sid mapping
    runlist%n = size(runlist%scans)
  end subroutine read_runlist

  function count_num_scan(file, object) result(n)
    implicit none
    character(len=*)   :: file
    character(len=512), optional :: object
    character(len=512) :: name, line
    integer(i4b)       :: unit, nobj, nscan, nfile, sid, i, j, k, n
    real(dp)           :: mjd(2)
    n = 0
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) nobj
    do i = 1, nobj
       read(unit,*) name, nscan
       do j = 1, nscan
          read(unit,*) sid, mjd, nfile
          if (present(object)) then
             if (trim(object) /= trim(name)) cycle
          end if
          n = n + 1
          do k = 1, nfile
             read(unit,'(a)') line
          end do
       end do
    end do
    close(unit)
  end function

  subroutine free_scan_info(a)
    implicit none
    type(comap_scan_info) :: a
    if (allocated(a%ss)) deallocate(a%ss)
  end subroutine

  subroutine copy_scan_info(a, b)
    implicit none
    type(comap_scan_info) :: a, b
    integer(i4b) :: i

    call free_scan_info(b)
    b%id = a%id
    b%nsub = a%nsub
    b%mjd  = a%mjd
    b%l1file = a%l1file
    b%object = a%object
    b%objectnum = a%objectnum
    allocate(b%ss(b%nsub))
    do i = 1, b%nsub
       b%ss(i)%id          = a%ss(i)%id
       b%ss(i)%mjd         = a%ss(i)%mjd
       b%ss(i)%az          = a%ss(i)%az
       b%ss(i)%el          = a%ss(i)%el
       b%ss(i)%dk          = a%ss(i)%dk
       b%ss(i)%lat         = a%ss(i)%lat
       b%ss(i)%lon         = a%ss(i)%lon
       b%ss(i)%time_of_day = a%ss(i)%time_of_day
       b%ss(i)%l2file      = a%ss(i)%l2file
       b%ss(i)%l3file      = a%ss(i)%l3file
       b%ss(i)%scanmode    = a%ss(i)%scanmode
       b%ss(i)%feature     = a%ss(i)%feature
    end do
  end subroutine

  subroutine free_runlist(runlist)
    implicit none
    type(comap_runlist) :: runlist
    integer(i4b) :: i, j
    if(allocated(runlist%scans)) then
       do i = 1, size(runlist%scans)
          call free_scan_info(runlist%scans(i))
       end do
       deallocate(runlist%scans)
    end if
  end subroutine

  subroutine copy_runlist(a, b)
    implicit none
    type(comap_runlist) :: a
    type(comap_runlist) :: b
    integer(i4b)        :: i
    call free_runlist(b)
    allocate(b%scans(size(a%scans)))
    do i = 1, size(a%scans)
       call copy_scan_info(a%scans(i),b%scans(i))
    end do
  end subroutine

  function mjd2chile_time(mjd) result(res)
    implicit none
    real(dp) :: mjd, res
    res = modulo((mjd-int(mjd))*24-3,24d0)
  end function

end module
