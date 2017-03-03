! This module provides access to pmac info as a function of time.
! Some of this could be read directly from the level3-files, but one
! sometimes need the information globally too. A more general
! solution might be to have l3gen or ces_validate output the
! whole-season condensed information in a special hdf file.
module quiet_pmac_mod
  use quiet_utils
  use locate_mod
  use quiet_fileutils
  implicit none

  real(dp), allocatable, private, dimension(:,:) :: az_currents
  real(dp), private :: mean_curr

contains

  subroutine initialize_pmac_mod(parfile, az_current_file)
    implicit none
    character(len=*),  intent(in)   :: parfile
    character(len=512)              :: azfile
    character(len=*),  optional     :: az_current_file
    logical(lgt),      save         :: initialized = .false.
    if(initialized) return
    if(present(az_current_file)) then
       azfile = az_current_file
    else
       call get_parameter(0, parfile, "AZ_CURRENT_FILE", par_string=azfile)
    end if
    if(azfile /= "") then
       call read_table(azfile, az_currents)
       mean_curr = mean(az_currents(:,2))
    end if
    initialized = .true.
  end subroutine

  ! Return the nearest az current to the specified date. It
  ! caches the last left side of the range (prev), to optimize
  ! the most common cases.
  function get_az_current(mjd) result(current)
    implicit none
    real(dp),     intent(in) :: mjd
    real(dp)                 :: current
    integer(i4b), save       :: prev = 0
    if(.not. allocated(az_currents)) then
       current = NaN
    else
       current = get_nearest_dp(az_currents, mjd, prev) - mean_curr
    end if
  end function

  function get_nearest_dp(table, key, hint) result(res)
    implicit none
    real(dp),     intent(in)    :: table(:,:), key
    integer(i4b), intent(inout) :: hint
    real(dp)                    :: res
    integer(i4b)                :: i, j, n
    res = NaN
    n = size(table,1)
    if(n == 0) return
    res = table(1,2)
    if(n == 1 .or. key <= table(1,1)) return
    res = table(n,2)
    if(key >= table(n,1)) return
    ! Ok, we now know that we have at least 2 points,
    ! and the boundary conditions have been taken care of
    if(hint <= 0 .or. hint >= n) then
       hint = min(n-1,locate(table(:,1), key))
    elseif(key < table(hint,1) .or. key > table(hint+1,1)) then
       hint = min(n-1,locate(table(:,1), key))
    end if
    ! Ok, key is between hint and hint+1. Find the closest one
    if(key-table(hint,1) < table(hint+1,1) - key) then
       res = table(hint,2)
    else
       res = table(hint+1,2)
    end if
  end function

end module
