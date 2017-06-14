! This module replaces quiet_module_mod. 
!
! Horn format:
! id telescope horn theta phi psi fwhm freq ok
!
! The horn and diode arrays are made public instead of having a horde
! of accessor functions, and this makes array operations simpler.
! For example, to get the board corresponding to each diode, perhaps for
! jackknife purposes, you would simply do quiet_horns(quiet_diodes%horn)%board.
! Diode indices are absolute, but the relative diode number is available as "sub".
!
! To translate from (mod,sub) to di, use: quiet_horns(mod)%diodes(sub)

module comap_detector_mod
  use quiet_utils
  implicit none

  type comap_detector
     integer(i4b) :: id, telescope, horn
     real(dp)     :: phi, theta, psi, freq, fwhm
     logical(lgt) :: ok
  end type

  type(comap_detector),  dimension(:), allocatable, public :: comap_detectors  ! (1:ndet)

contains

  subroutine init_detector_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    character(len=512)           :: dfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    call get_parameter(0, parfile, "DETECTOR_FILE", par_string=dfile)
    call read_detectors(dfile, comap_detectors)
    call setup_diodes(quiet_horns, quiet_diodes)
    initialized = .true.
  end subroutine

  ! Helper functions below
  subroutine read_detectors(filename, detectors)
    implicit none
    character(len=*),                            intent(in)    :: filename
    type(quiet_horn), dimension(:), allocatable, intent(inout) :: detectors
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
       horns(id)%id        = id
       horns(id)%telescope = telescope
       horns(id)%horn      = horn
       horns(id)%theta     = theta * DEG2RAD          ! Input in deg
       horns(id)%phi       = phi   * DEG2RAD          ! Input in deg
       horns(id)%psi       = psi   * DEG2RAD          ! Input in deg
       horns(id)%fwhm      = fwhm  * DEG2RAD / 60     ! Input in arcmin
       horns(id)%freq      = freq  * 1d9              ! Input in GHz
       horns(id)%ok        = ok 
    end do
    2 close(unit)
  end subroutine

  function get_num_detectors() result(res)
    implicit none
    integer(i4b) :: res
    res = size(quiet_detector)
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
    fwhm = comap_detectors(mod)%fwhm
  end function

  function get_max_fwhm() result(fwhm)
    implicit none
    real(dp)     :: fwhm
    fwhm = maxval(comap_detectors%fwhm)
  end function

  function get_detector_angle(id) result(res)
    implicit none
    integer(i4b) :: mod, di
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
    integer(i4b) :: diode, mod, di
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
