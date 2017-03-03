module quiet_apex_mod
  use healpix_types
  use spline_1d_mod
  use quiet_utils
  implicit none

  real(dp),     parameter, private :: one_minute = 1.d0 / (24.d0*60.d0)
  integer(i4b),            private :: nsamp

  integer(i4b), parameter :: APEX_PWV            = 1
  integer(i4b), parameter :: APEX_DEW_POINT      = 2
  integer(i4b), parameter :: APEX_HUMIDITY       = 3
  integer(i4b), parameter :: APEX_PRESSURE       = 4
  integer(i4b), parameter :: APEX_TEMPERATURE    = 5
  integer(i4b), parameter :: APEX_WIND_DIRECTION = 6
  integer(i4b), parameter :: APEX_WIND_SPEED     = 7
  integer(i4b), parameter :: APEX_NUM_TYPE       = 7

  character(len=24), dimension(APEX_NUM_TYPE) :: apex_names = [&
     'APEX_PWV                ',&
     'APEX_DEW_POINT          ',&
     'APEX_HUMIDITY           ',&
     'APEX_PRESSURE           ',&
     'APEX_TEMPERATURE        ',&
     'APEX_WIND_DIRECTION     ',&
     'APEX_WIND_SPEED         ']

  real(dp), allocatable, dimension(:),   private :: time
  real(dp), allocatable, dimension(:,:), private :: data

contains

  subroutine initialize_quiet_apex_mod(parfile)
    implicit none

    character(len=*), intent(in) :: parfile

    integer(i4b)       :: unit, i
    character(len=512) :: apexfile

    unit = getlun()

    call get_parameter(unit, parfile, 'APEXFILE', par_string=apexfile)

    ! Find number of time samples
    nsamp = 0
    open(unit, file=trim(apexfile))
    do while (.true.)
       read(unit,*,end=1)
       nsamp = nsamp+1
    end do
1   close(unit)

    ! Read data file
    allocate(time(nsamp), data(nsamp,APEX_NUM_TYPE))
    open(unit, file=trim(apexfile))
    do i = 1, nsamp
       read(unit,*) time(i), data(i,:)
    end do
    close(unit)

  end subroutine initialize_quiet_apex_mod

  function get_apex_data(mjd, apex_type)
    implicit none
    
    real(dp),     intent(in) :: mjd
    integer(i4b), intent(in) :: apex_type
    real(dp)                 :: get_apex_data

    integer(i4b) :: i
    real(dp)     :: x

    i = locate(time, mjd)

    ! Return 0 if we're outside valid range
    if (i <= 1 .or. i > nsamp) then
       get_apex_data = 0.d0
       return
    end if

    ! Return 0 if APEX has stopped taking data for the current time
    if (abs(mjd-time(i)) > 60.d0*one_minute .or. abs(time(i+1)-mjd) > 60.d0*one_minute) then
       get_apex_data = 0.d0
       return
    end if

    ! Do a linear interpolation if things look good
    x = (mjd-time(i)) / (time(i+1)-time(i))
    get_apex_data = (1.d0-x) * data(i,apex_type) + x * data(i+1,apex_type)

  end function get_apex_data

  subroutine get_apex_range(mjd, apex_type, vals)
    implicit none
    
    real(dp),     dimension(:), intent(in)  :: mjd
    integer(i4b),               intent(in)  :: apex_type
    real(dp),     dimension(:), intent(out) :: vals

    integer(i4b) :: i

    do i = 1, size(mjd)
       vals(i) = get_apex_data(mjd(i), apex_type)
    end do

  end subroutine get_apex_range


end module quiet_apex_mod
