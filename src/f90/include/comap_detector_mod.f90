module comap_detector_mod
  use quiet_utils
  implicit none

  type comap_detector
     integer(i4b) :: id, telescope, horn
     real(dp)     :: phi, theta, psi, freq, fwhm
     logical(lgt) :: ok
  end type

  integer(i4b), private :: num_sideband
  type(comap_detector),  dimension(:), allocatable, public :: comap_detectors  ! (1:ndet)

contains

  subroutine initialize_detector_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    character(len=512)           :: dfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    call get_parameter(0, parfile, "DETECTOR_FILE", par_string=dfile)
    call get_parameter(0, parfile, "NUM_SIDEBAND", par_int=num_sideband)
    call read_detectors(dfile, comap_detectors)
    initialized = .true.
  end subroutine

  ! Helper functions below
  subroutine read_detectors(filename, detectors)
    implicit none
    character(len=*),                                intent(in)    :: filename
    type(comap_detector), dimension(:), allocatable, intent(inout) :: detectors

    character(len=512) :: line
    integer(i4b)       :: i, j, k, m, n, id, unit, telescope, horn
    real(dp)           :: theta, phi, psi, fwhm, freq
    logical(lgt)       :: ok
    if (allocated(detectors)) deallocate(detectors)
    unit = getlun()
    open(unit,file=filename,status="old",action="read")
    ! Get the maximum detector number
    n = -1
    do
       read(unit,'(a)',end=1) line
       if(line(1:1) == "#" .or. line == "") cycle
       read(line,*) k
       n = max(n,k)
    end do
    1 rewind(unit)
    call assert(n>0, "Error! No valid detectors in "//trim(filename))
    allocate(detectors(n))
    do
       read(unit,'(a)',end=2) line
       if(line(1:1) == "#" .or. line == "") cycle
       read(line,*) id, telescope, horn, theta, phi, psi, fwhm, freq, ok
       detectors(id)%id        = id
       detectors(id)%telescope = telescope
       detectors(id)%horn      = horn
       detectors(id)%theta     = theta * DEG2RAD          ! Input in deg
       detectors(id)%phi       = phi   * DEG2RAD          ! Input in deg
       detectors(id)%psi       = psi   * DEG2RAD          ! Input in deg
       detectors(id)%fwhm      = fwhm  * DEG2RAD / 60     ! Input in arcmin
       detectors(id)%freq      = freq  * 1d9              ! Input in GHz
       detectors(id)%ok        = ok 
    end do
    2 close(unit)
  end subroutine

  function get_num_dets() result(res)
    implicit none
    integer(i4b) :: res
    res = size(comap_detectors)
  end function

  function get_num_sideband() result(res)
    implicit none
    integer(i4b) :: res
    res = num_sideband
  end function

  function detector_valid(id) result(res)
    implicit none
    integer(i4b) :: id
    logical(lgt) :: res
    if(id < 1 .or. id > size(comap_detectors)) then
       res = .false.
    else
       res = comap_detectors(id)%ok
    end if
  end function

  subroutine get_detector_pos(id, theta, phi)
    implicit none
    integer(i4b) :: id
    real(dp)     :: theta, phi
    theta = comap_detectors(id)%theta
    phi   = comap_detectors(id)%phi
  end subroutine

  function get_detector_fwhm(id) result(fwhm)
    implicit none
    integer(i4b) :: id
    real(dp)     :: fwhm
    fwhm = comap_detectors(id)%fwhm
  end function

  function get_max_fwhm() result(fwhm)
    implicit none
    real(dp)     :: fwhm
    fwhm = maxval(comap_detectors%fwhm)
  end function

  function get_detector_angle(id) result(res)
    implicit none
    integer(i4b) :: id
    real(dp)     :: res
    res = comap_detectors(id)%psi
  end function

  function get_detector_freq(id) result(res)
    implicit none
    integer(i4b) :: id
    real(dp)     :: res
    res = comap_detectors(id)%freq
  end function

  function is_alive(id) result(res)
    implicit none
    integer(i4b) :: id
    logical(lgt) :: res
    res = comap_detectors(id)%ok
  end function

  function get_focalplane_radius() result(res)
    implicit none
    real(dp) :: res
    res = maxval(comap_detectors%theta)
  end function

  function get_center_detector() result(res)
    implicit none
    integer(i4b) :: res
    res = minloc(comap_detectors%theta,1)
  end function

  ! Warning: This routine is slow, so use sparingly
  subroutine detector_abs2rel(telescope, horn, id)
    implicit none
    integer(i4b), intent(in)  :: telescope, horn
    integer(i4b), intent(out) :: id
    integer(i4b) :: i
    do i = 1, size(comap_detectors)
       if (comap_detectors(i)%telescope == telescope .and. comap_detectors(i)%horn == horn) then
          id = i
          return
       end if
    end do
    id = -1
  end subroutine

  subroutine detector_rel2abs(id, telescope, horn)
    implicit none
    integer(i4b) :: id, telescope, horn
    telescope = comap_detectors(id)%telescope
    horn      = comap_detectors(id)%horn
  end subroutine detector_rel2abs

end module
