module quiet_ces_mod
  use healpix_types
  use quiet_utils
  use sort_utils
  implicit none

  ! Information about a quiet_ces_info. Only the information that can be
  ! extracted directly from the level2-runlist. Category and object
  ! are indices into name arrays.
  type quiet_ces_info
    integer(i4b)       :: cid
    real(dp)           :: mjd(2), az, el, dk, lat, lon, time_of_day
    character(len=512) :: l2file, l3file, object
    character(len=512), allocatable, dimension(:) :: l1files
  end type quiet_ces_info

  type quiet_runlist
     integer(i4b)                               :: n
     type(quiet_ces_info), dimension(:), allocatable :: ceses
     integer(i4b),         dimension(:), allocatable :: cidmap
  end type quiet_runlist

  integer(i4b), parameter :: PROP_TOD = 1, PROP_AZ = 2, PROP_EL = 3, &
    & PROP_DECK = 4, PROP_DUR = 5, PROP_MJD = 6

  type(quiet_runlist)          :: ces_db
  integer(i4b),    allocatable :: cid_list(:), cid_sort(:)
  logical(lgt),        private :: initialized = .false.
  character(len=512),  private :: l1dir, l2dir, l3dir

contains

  ! Public interface: These are meant to be called from outside.

  subroutine initialize_ces_mod(parfile, planck)
    implicit none
    character(len=*)       :: parfile
    character(len=512)     :: runlist
    logical(lgt), optional :: planck
    logical(lgt)           :: planck_
    if(initialized) return
    call get_parameter(0, parfile, 'RUNLIST',     par_string=runlist)
    call get_parameter(0, parfile, 'LEVEL1_DIR',  par_string=l1dir)
    call get_parameter(0, parfile, 'LEVEL2_DIR',  par_string=l2dir)
    call get_parameter(0, parfile, 'LEVEL3_DIR',  par_string=l3dir)

    planck_ = .false.; if(present(planck)) planck_ = planck
    if (planck_) then
       call read_runlist_planck(runlist, l2dir, l3dir, ces_db)
    else
       call read_runlist(runlist, l1dir, l2dir, l3dir, ces_db)
    end if
    allocate(cid_list(ces_db%n), cid_sort(ces_db%n))
    cid_list = ces_db%ceses%cid
    cid_sort = cid_list
    call quicksort_int(cid_sort)
    initialized = .true.
  end subroutine

  function get_num_ces() result(res)
    implicit none
    integer(i4b) :: res
    res = size(ces_db%ceses)
  end function

  function lookup_ces(cid) result(res)
    implicit none
    integer(i4b) :: cid, res
    res = 0
    if(cid < 0 .or. cid > size(ces_db%cidmap)) return
    res = ces_db%cidmap(cid)
  end function

  subroutine get_ces_info(ind, ces)
    implicit none
    type(quiet_ces_info) :: ces
    integer(i4b)    :: ind
    call copy_ces_info(ces_db%ceses(ind), ces)
  end subroutine

  subroutine get_l2_param(cid, index, mjd, az, el, dk, time, target)
    implicit none
    integer(i4b),                   intent(in),  optional :: cid, index
    real(dp),                       intent(out), optional :: az, el, dk, time
    real(dp),         dimension(2), intent(out), optional :: mjd
    character(len=*),               intent(out), optional :: target
    type(quiet_ces_info) :: info

    if(present(index)) then
       call get_ces_info(index, info)
    else
       if(.not. present(cid)) then
          write(stderr,*) "Either index or cid must be given to get_l2_param!"
          stop
       end if
       call assert(cid > 0, "Non-positive cid in get_l2_param!")
       call get_ces_info(ces_db%cidmap(cid), info)
    end if

    if (present(az))   then; az   = info%az; end if
    if (present(el))   then; el   = info%el; end if
    if (present(dk))   then; dk   = info%dk; end if
    if (present(time)) then; time = info%time_of_day; end if
    if (present(mjd))  then; mjd  = info%mjd; end if
    if (present(target)) then; target = info%object; end if
    call free_ces_info(info)
  end subroutine get_l2_param

  !-----------------------------------------!
  !--------- The helper routines -----------!
  !-----------------------------------------!

  subroutine read_runlist(file, l1dir, l2dir, l3dir, runlist)
    implicit none
    character(len=*)   :: file
    character(len=512) :: name, l1dir, l2dir, l3dir, line
    integer(i4b)       :: unit, nobj, nces, nfile, cid, i, j, k, n, cnum, cmax
    real(dp)           :: mjd(2)
    type(quiet_runlist):: runlist
    type(quiet_ces_info)    :: ces
    n = count_num_ces(file)
    call free_runlist(runlist)
    allocate(runlist%ceses(n))
    cnum = 0
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) nobj
    do i = 1, nobj
       read(unit,*) name, nces
       do j = 1, nces
          cnum = cnum+1
          read(unit,*) ces%cid, ces%mjd, nfile, ces%az, ces%el, ces%dk, ces%lon, ces%lat
          ces%object      = name
          ces%time_of_day = mjd2chile_time(mean(ces%mjd))
          allocate(ces%l1files(nfile))
          do k = 1, nfile
             read(unit,fmt="(a)") line
             ces%l1files(k) = trim(l1dir) // "/" // trim(get_token(line, " ", 3))
          end do
          ces%l2file = trim(l2dir) // "/" // trim(name) // "/" // trim(name) // "_" // &
           & trim(itoa(ces%cid)) // ".hdf"
          ces%l3file = trim(l3dir) // "/" // trim(name) // "/" // trim(name) // "_" // &
           & trim(itoa(ces%cid)) // ".hdf"
          call copy_ces_info(ces, runlist%ceses(cnum))
          call free_ces_info(ces)
       end do
    end do
    close(unit)
    ! Set up the cid mapping
    cmax = maxval(runlist%ceses%cid)
    allocate(runlist%cidmap(cmax))
    runlist%cidmap = 0
    do i = 1, size(runlist%ceses)
       runlist%cidmap(runlist%ceses(i)%cid) = i
    end do
    runlist%n = size(runlist%ceses)
  end subroutine read_runlist

  subroutine read_runlist_planck(file, l2dir, l3dir, runlist)
    implicit none
    character(len=*)   :: file
    character(len=512) :: name, l2dir, l3dir
    integer(i4b)       :: unit, nobj, nces, nfile, cid, i, cmax
    real(dp)           :: mjd(2)
    type(quiet_runlist):: runlist
    type(quiet_ces_info)    :: ces
    call free_runlist(runlist)
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) name, nces
    allocate(runlist%ceses(nces))
    do i = 1, nces
       read(unit,*) ces%cid
       ces%object      = name
       ces%l2file = trim(l2dir) // "/" // trim(name) // "/" // trim(itoa(ces%cid)) // ".h5"
       ces%l3file = trim(l3dir) // "/" // trim(name) // "/" // trim(itoa(ces%cid)) // ".hdf"
       call copy_ces_info(ces, runlist%ceses(i))
       call free_ces_info(ces)
    end do
    close(unit)
    ! Set up the cid mapping
    cmax = maxval(runlist%ceses%cid)
    allocate(runlist%cidmap(cmax))
    runlist%cidmap = 0
    do i = 1, size(runlist%ceses)
       runlist%cidmap(runlist%ceses(i)%cid) = i
    end do
    runlist%n = size(runlist%ceses)
  end subroutine read_runlist_planck

  function count_num_ces(file) result(n)
    implicit none
    character(len=*)   :: file
    character(len=512) :: name
    integer(i4b)       :: unit, nobj, nces, nfile, cid, i, j, k, n
    real(dp)           :: mjd(2)
    n = 0
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    read(unit,*) nobj
    do i = 1, nobj
       read(unit,*) name, nces
       n = n+nces
       do j = 1, nces
          read(unit,*) cid, mjd, nfile
          do k = 1, nfile
             read(unit,*)
          end do
       end do
    end do
    close(unit)
  end function

  subroutine free_ces_info(a)
    implicit none
    type(quiet_ces_info) :: a
    if(allocated(a%l1files)) deallocate(a%l1files)
  end subroutine

  subroutine copy_ces_info(a, b)
    implicit none
    type(quiet_ces_info) :: a, b
    call free_ces_info(b)
    allocate(b%l1files(size(a%l1files)))
    b = a
  end subroutine

  subroutine free_runlist(runlist)
    implicit none
    type(quiet_runlist) :: runlist
    integer(i4b) :: i
    if(allocated(runlist%ceses)) then
       do i = 1, size(runlist%ceses)
          call free_ces_info(runlist%ceses(i))
       end do
       deallocate(runlist%ceses)
    end if
    if(allocated(runlist%cidmap)) deallocate(runlist%cidmap)
  end subroutine

  subroutine copy_runlist(a, b)
    implicit none
    type(quiet_runlist) :: a
    type(quiet_runlist) :: b
    integer(i4b)        :: i
    call free_runlist(b)
    allocate(b%ceses(size(a%ceses)))
    allocate(b%cidmap(size(a%cidmap)))
    do i = 1, size(a%ceses)
       call copy_ces_info(a%ceses(i),b%ceses(i))
    end do
    b%cidmap = a%cidmap
  end subroutine

  function mjd2chile_time(mjd) result(res)
    implicit none
    real(dp) :: mjd, res
    res = modulo((mjd-int(mjd))*24-3,24d0)
  end function

end module
