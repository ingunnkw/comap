module comap_pointing_mod
  use healpix_types
  use pix_tools
  use quiet_utils
  use math_tools
  use comap_detector_mod
  use powell_mod
  use locate_mod
  implicit none

  ! Note on coordinate systems /-conversions:
  ! Tele: Telescope pointing in horizontal coords according to encoder
  ! Tele + mount model + focal plane -> hor
  ! Hor: Individual detector pointing in horizontal coords
  ! Then we do apparent corrections and conversion to celestial and galactic coords

  real(dp) :: COMAP_GEODETIC_LONGITUDE  = -118.28194444d0  ! 118 degrees, 16 minutes, 55 seconds WEST
  real(dp) :: COMAP_GEODETIC_LATITUDE   =   37.23388888d0  ! 37 degrees, 14 minute, 2 seconds NORTH
  real(dp) :: COMAP_ALTITUDE            =   1222.d0        ! meters

  real(dp), parameter :: COMAP_GEOCENTRIC_LATITUDE = 37.23388888d0   ! degrees; same as GEODETIC for now
  real(dp), parameter :: COMAP_RADIUS              = 1.000277634d0   ! units of earth radii

  real(dp), parameter :: J2000_RA_NGP = 192.859498564d0       !< Right ascension of the North Galactic Pole, J2000 epoch
  real(dp), parameter :: J2000_DEC_NGP = 27.128335953d0       !< Declination of the North Galactic Pole, J2000 epoch
  real(dp), parameter :: J2000_RA_GC = 266.405090965d0        !< Right ascension of the galactic center, J2000 epoch
  real(dp), parameter :: J2000_DEC_GC = -28.936174393d0       !< Declination of the galactic center, J2000 epoch
  real(dp), parameter :: J2000_EULER_ALPHA = -57.068351386d0  !< Precomputed Euler angle for conversion between equitorial and galactic coordinates, J2000 epoch
  real(dp), parameter :: J2000_EULER_BETA = -62.871663896d0   !< Precomputed Euler angle for conversion between equitorial and galactic coordinates, J2000 epoch
  real(dp), parameter :: J2000_EULER_GAMMA = -192.859498564d0 !< Precomputed Euler angle for conversion between equitorial and galactic coordinates, J2000 epoch

  real(dp), parameter :: UT2MJD = 40587.d0 ! Offset, in days, between Unix Time and Modified Julian Day
  real(dp), parameter :: DTOR  = 1.7453292519943295769d-2 ! Conversion from degrees to radians
  real(dp), parameter :: RTOD  = 1.d0 / 1.7453292519943295769d-2 ! Conversion from radians to degrees
  real(dp), parameter :: LEAPS = 34
  real(dp), parameter :: DAYSECS = 24d0*60*60

  ! EQU and CEL are synonyms
  integer(i4b), parameter :: COORD_TELE = 0, COORD_HOR = 1, &
    & COORD_EQU = 2, COORD_CEL = 2, COORD_GAL = 3, COORD_CENTER = 4, &
    & NUM_COORD = 5
  ! The next step of the shortest path from one system to another
  integer(i4b)            :: next_step(0:NUM_COORD-1,0:NUM_COORD-1)

  integer(i4b),                                private :: ndet
  real(dp),     dimension(3,3),                private :: M_cel2gal, M_gal2cel
  real(dp),     dimension(:,:,:), allocatable, private :: diode_euler

  real(dp),     dimension(:),     allocatable, private :: mjds
  integer(i4b),                                private :: prev_cnum

  type inv_rot_params
     integer(i4b) :: det_id
     real(dp)     :: mjd, M(3,3)
  end type

  type(inv_rot_params), private :: inv_rot

  logical(lgt), private :: initialized = .false.

contains

  subroutine initialize_simple_pointing
    implicit none
    logical(lgt), save :: initialized = .false.
    if(initialized) return
    M_cel2gal = get_matrix_equ2gal()
    M_gal2cel = transpose(M_cel2gal)
    initialized = .true.
  end subroutine

  subroutine initialize_comap_pointing_mod(paramfile)
    implicit none
    character(len=*), intent(in) :: paramfile

    integer(i4b)       :: i, nparam
    character(len=512) :: focalplane
    real(dp) :: phi1, theta1, psi1, phi2, theta2, psi2, mjd
    if(initialized) return
    call initialize_simple_pointing
    call initialize_detector_mod(paramfile)
    call setup_detector_euler_matrices(diode_euler)

    next_step = &
      & transpose(reshape((/ &
      &  COORD_TELE, COORD_HOR, COORD_HOR, COORD_HOR, COORD_HOR, &
      &  COORD_HOR,  COORD_HOR, COORD_CEL, COORD_CEL, COORD_CEL, &
      &  COORD_CEL,  COORD_CEL, COORD_CEL, COORD_GAL, COORD_GAL, &
      &  COORD_CEL,  COORD_CEL, COORD_CEL, COORD_GAL, COORD_CENTER, &
      &  COORD_GAL,  COORD_GAL, COORD_GAL, COORD_GAL, COORD_CENTER  &
      & /), (/NUM_COORD, NUM_COORD/)))

    initialized = .true.
  end subroutine initialize_comap_pointing_mod

  ! NOTE: All angles here are in the healpix convention!

  ! Driver routine for switching between coordinate systems.
  ! Note: mjd is needed for TELE and HOR to
  ! the others. mod and diode is needed to and from TELE, unless you
  ! want the borepoint.
  !
  ! The coordinate systems are connected like this:
  !
  ! center-gal-cel-hor-tele
  !
  ! We transform stepwise twoards the target. The setup does support
  ! branches in the tree, but we don't have any yet. The current
  ! implementation is about the same length as the old one,
  ! but this one scales as O(N), while the previous one scaled
  ! as O(N^2), where N is the number of coordinate systems.
  subroutine coord_convert(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, det_id, phic, thetac, euler)
    implicit none
    integer(i4b)           :: sys1, sys2, sys, nsys
    integer(i4b), optional :: det_id
    real(dp)               :: phi1, theta1, psi1, phi2, theta2, psi2
    real(dp),     optional :: mjd, thetac, phic, euler(3,3)
    real(dp)               :: mat(3,3)

    if(sys1 == sys2) then
       phi2 = phi1; theta2 = theta1; psi2 = psi1
       return
    end if

    ! Step 1: Go to matrix representation
    sys = sys1
    if(sys1 == COORD_TELE) then
       mat = angles2hor(mjd, phi1, theta1, psi1, det_id)
       sys = COORD_HOR
    else
       mat = angles2rot(phi1, theta1, psi1)
    end if

    ! Step 2: Transform mat until we are done
    do
       nsys = next_step(sys, sys2)
       if(nsys == sys) exit
       select case(sys)
          case(COORD_TELE)
             call abort("Bug in COORD_TELE conversion!")
          case(COORD_HOR)
             select case(nsys)
                case(COORD_CEL);    mat = matmul(rot_hor2equ(mjd), mat)
                case(COORD_TELE);   continue
                case default; call abort("Bug in COORD_HOR conversion!")
             end select
          case(COORD_CEL)
             select case(nsys)
                case(COORD_GAL);    mat = matmul(rot_equ2gal(), mat)
                case default; call abort("Bug in COORD_CEL conversion!")
             end select
          case(COORD_GAL)
             select case(nsys)
                case(COORD_CEL);    mat = matmul(rot_gal2equ(), mat)
                case(COORD_CENTER); mat = matmul(rot_gal2center(phic,thetac), mat)
                case default; call abort("Bug in COORD_GAL conversion!")
             end select
          case(COORD_CENTER)
             select case(nsys)
                case(COORD_GAL);    mat = matmul(rot_center2gal(phic,thetac), mat)
                case default; call abort("Bug in COORD_CENTER conversion!")
             end select
          case default; call abort("Unhandled coordinate system " // trim(itoa(sys)) // "!")
       end select
       sys = nsys
    end do

    ! Step 3: Read out result
    if(sys2 == COORD_TELE) then
       call hor2angles(mat, mjd, phi2, theta2, psi2, det_id)
    else
       call rot2angles(mat, phi2, theta2, psi2)
    end if
    if(present(euler)) euler = mat
  end subroutine

  ! Shortcuts to go in and out of matrix representation
  function angles2rot(phi, theta, psi) result(mat)
    implicit none
    real(dp) :: phi, theta, psi, mat(3,3)
    call compute_euler_matrix_zyz(phi, theta, psi, mat)
  end function

  subroutine rot2angles(mat, phi, theta, psi)
    implicit none
    real(dp) :: phi, theta, psi, mat(3,3)
    call convert_euler_matrix_to_angles_zyz(mat, phi, theta, psi)
  end subroutine

  ! Composite utility functions for going from reported boresight angles
  ! to corrected diode pointing in horizontal coordinates, and back.
  ! The goal of this rewrite was to get things as modular as possible,
  ! but I think it is natural to have these two shortcut functions.
  function angles2hor(mjd, phi, theta, psi, det_id) result(mat)
    implicit none
    real(dp),     intent(in)  :: phi, theta, psi, mjd
    real(dp), dimension(3,3)  :: mat
    integer(i4b)              :: det
    integer(i4b), optional    :: det_id
    det = -1; if(present(det_id)) det = det_id
    mat = rot_boresight2hor(mjd, phi, theta, psi, det)
  end function

  subroutine hor2angles(mat, mjd, phi, theta, psi, det_id)
    implicit none
    real(dp)                  :: mjd, phi, theta, psi, mat(3,3)
    integer(i4b)              :: det
    integer(i4b), optional    :: det_id
    det = -1; if(present(det_id)) det = det_id
    call rot_hor2boresight(mjd, mat, det, phi, theta, psi)
  end subroutine

  ! The non-complicated rotations
  function rot_hor2equ(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = get_matrix_hor2equ(mjd2lst(mjd, COMAP_GEODETIC_LONGITUDE), &
      & COMAP_GEODETIC_LATITUDE)
  end function

  function rot_hor2gal(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = matmul(rot_equ2gal(),rot_hor2equ(mjd))
  end function

  function rot_equ2hor(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = get_matrix_equ2hor(mjd2lst(mjd, COMAP_GEODETIC_LONGITUDE), &
      & COMAP_GEODETIC_LATITUDE)
  end function

  function rot_equ2gal() result(mat)
    implicit none
    real(dp) :: mat(3,3)
    mat = M_cel2gal
  end function

  function rot_gal2equ() result(mat)
    implicit none
    real(dp) :: mat(3,3)
    mat = M_gal2cel
  end function

  function rot_gal2hor(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = transpose(rot_hor2gal(mjd))
  end function

  ! Rotate from galactic coordinates to a coordinate
  ! system where coordinates phi, theta become
  ! 0, pi/2. That is, the object at phi, theta
  ! is rotated to be where the galactic center would be
  function rot_gal2center(phi, theta) result(mat)
    implicit none
    real(dp) :: mat(3,3), phi, theta
    call compute_euler_matrix_zyz(0d0, pi/2-theta, -phi, mat)
  end function

  function rot_center2gal(phi, theta) result(mat)
    implicit none
    real(dp) :: mat(3,3), phi, theta
    mat = transpose(rot_gal2center(phi, theta))
  end function

  !!!!!!!!!!!!!!!!!!!!
  ! Helper functions !
  !!!!!!!!!!!!!!!!!!!!

  ! Rotate from module coordinates to boresight. If diode is not included,
  ! the diode angle is not part of the rotation.
  function rot_detector2boresight(det_id) result(mat)
    implicit none
    integer(i4b), intent(in)           :: det_id
    real(dp), dimension(3,3)           :: mat
    integer(i4b)                       :: di
    mat = diode_euler(:,:,det_id)
  end function

  function rot_boresight2module(det_id) result(mat)
    implicit none
    integer(i4b), intent(in)           :: det_id
    real(dp), dimension(3,3)           :: mat
    mat = transpose(rot_detector2boresight(det_id))
  end function

  ! Perform the full rotation from naive boresight coordinates
  ! to corrected horizontal coordinates.
  function rot_boresight2hor(mjd, phi, theta, psi, det_id) result(mat)
    implicit none
    real(dp),     intent(in)  :: phi, theta, psi, mjd
    integer(i4b), intent(in)  :: det_id
    real(dp)                  :: p(3), az, el, dk, coaz, ab(2), ab2(2), p1(2), p2(2)
    real(dp)                  :: del, coeff, dev, dpsi
    real(dp), dimension(3,3)  :: mat, M_bore, M_az, M_el, M_fp, fp1, fp2
    real(dp), dimension(3,3)  :: M_test1, M_test2, M_id
    p   = [phi,theta,psi]
    az  = -p(1); el = pi/2-p(2); dk = p(3)
    mat = get_identity(3)
    ! Boresight pointing
    call compute_euler_matrix_zyz(p(1), p(2), p(3), M_bore)
    mat = matmul(M_bore, mat)

    mat = matmul(mat, diode_euler(:,:,det_id))
  end function

  ! Go from horizontal boresight coordinates to telescope coordinates.
  ! We use powell. This is a pretty slow function. It may be sped up
  ! by using a better starting point by using the fiducial focalplane
  ! layout as an approximation.
  subroutine rot_hor2boresight(mjd, mat_hor, det_id, phi, theta, psi)
    implicit none
    real(dp) :: mjd, mat_hor(3,3), phi, theta, psi, p(3)
    integer(i4b) :: i, det_id, err
    inv_rot%det_id = det_id
    inv_rot%mjd    = mjd
    inv_rot%M      = mat_hor
    call convert_euler_matrix_to_angles_zyz(mat_hor, p(1), p(2), p(3))
    call powell(p, powell_inv_rot, err)
    phi = p(1); theta = p(2); psi = p(3)
  end subroutine

  function powell_inv_rot(p) result(res)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: res
    res = sum(abs(rot_boresight2hor(inv_rot%mjd, p(1),p(2),p(3), inv_rot%det_id)-inv_rot%M))
  end function

  function parse_coord_name(name) result(res)
    implicit none
    character(len=*) :: name
    integer(i4b)     :: res
    if(name == 'galactic' .or. name == 'gal') then; res = COORD_GAL
    elseif(name == 'celestial' .or. name == 'cel') then; res = COORD_CEL
    elseif(name == 'equatorial' .or. name == 'equ') then; res = COORD_EQU
    elseif(name == 'horizontal' .or. name == 'hor') then; res = COORD_HOR
    elseif(name == 'telesope' .or. name == 'tele') then; res = COORD_TELE
    elseif(name == 'centered' .or. name == 'center') then; res = COORD_CENTER
    else
       write(*,*) "Unrecognized coordinate system '" // trim(name) // "'!"
       stop
    end if
  end function

  ! What the normal convention is depends on the system. In all cases, we
  ! swap between healpix and the norm.
  subroutine swap_coordinate_convention(phi, theta, psi, system)
    implicit none
    real(dp)     :: phi, theta, psi
    integer(i4b), optional :: system
    integer(i4b) :: sys
    sys = COORD_HOR; if(present(system)) sys = system
    if(sys == COORD_HOR .or. sys == COORD_TELE) then
       phi = -phi
       theta = 0.5d0*pi - theta
    else
       theta = 0.5d0*pi - theta
    end if
  end subroutine

  ! Convert from Unix time to modified Julian day.
  !
  ! Unix time is defined as seconds since midnight UTC on January 1, 1970. For this function, 
  ! it is given as two integer values: seconds (current value is roughly \f$1.2e+09\f$) and 
  ! microseconds (which has values between 0 and 999,999). Assuming 4 byte integers, 
  ! the \a sec argument is large enough to be valid until approximately the year 2106.
  !
  ! The returned value is the decimal number of days since midnight of November 17, 1858, UTC. 
  ! Current values of modified Julian day are roughly \f$2.45e+06\f$.
  function unixtime2mjd(sec, usec) 
    implicit none

    real(dp), intent(in) :: sec, usec
    real(dp)             :: unixtime2mjd

    ! convert seconds and microseconds to days (divide by 86400)
    ! then add the offset (40587.0 days) between Unix time (Jan 1, 1970) and MJD (Nov 17, 1858)
    unixtime2mjd = (sec + (usec / 1.d6)) / 86400.d0 + UT2MJD
  end function unixtime2mjd

  ! Convert from modified Julian day to local sidereal time.
  !
  ! All QPoint functions that deal with equitorial coordinates take \a lst as an argument. 
  ! The qpoint_mjd2lst function takes as an argument the modified Julian day (\a mjd), 
  ! which is the decimal number of days since November 17, 1958, UTC.
  !
  ! The value returned is the local sidereal time, calculated at the specified longitude. 
  ! The result is in units of degrees, in the range [0,360).
  function mjd2lst(mjd, longitude)
    implicit none
    real(dp), intent(in) :: mjd, longitude
    real(dp)             :: mjd2lst

    real(dp) :: s, du, tu, gmst, lst

    s  = 86400.d0 * (mjd - floor(mjd))   ! seconds since midnight UTC
    du = mjd - 51544.5d0                 ! days since J2000 epoch
    tu = du / 36525.d0                   ! convert du to centuries

    ! Greenwich Mean Sidereal Time
    ! Formula from Astrophysical Formulae, Volume 2 (3Ed) by Kenneth Lang, 1999.
    ! Conforms to IAU 1976 System of Astronomical Constants, 1980 IAU Theory of Nutation, FK5 catalog
    gmst = s + 24110.54841d0 + tu * (8640184.812866d0 + tu * (0.093104d0 - tu * 0.0000062d0)); 

    gmst = gmst * (360.d0 / 86400.d0);
  
    ! add telescope longitude, fix to the range [0,360)
    lst = gmst + longitude
    do while (lst < 0.d0)    
       lst = lst + 360.d0
    end do
    do while (lst >= 360.d0) 
       lst = lst - 360.d0
    end do
    mjd2lst = lst
  end function mjd2lst

  function get_matrix_hor2equ(lst, lat) result(matrix)
    implicit none

    real(dp), intent(in)     :: lst, lat
    real(dp), dimension(3,3) :: matrix

    call compute_euler_matrix_zyz(DTOR * (lst - 180.d0), DTOR * (lat - 90.d0), 0.d0, matrix)
  end function get_matrix_hor2equ

  function get_matrix_equ2hor(lst, lat) result (matrix)
    implicit none
    real(dp), intent(in)     :: lst, lat
    real(dp), dimension(3,3) :: matrix
    matrix = get_matrix_hor2equ(lst, lat)
    matrix = transpose(matrix)
  end function get_matrix_equ2hor

  function get_matrix_equ2gal() result(matrix)
    implicit none
    real(dp), dimension(3,3) :: matrix
    call compute_euler_matrix_zyz(DTOR * J2000_EULER_ALPHA, DTOR * J2000_EULER_BETA, &
         & DTOR * J2000_EULER_GAMMA, matrix)
  end function get_matrix_equ2gal

  function get_matrix_gal2equ() result(matrix)
    implicit none
    real(dp), dimension(3,3) :: matrix
    matrix = get_matrix_equ2gal()
    matrix = transpose(matrix)
  end function get_matrix_gal2equ

  subroutine setup_detector_euler_matrices(mats)
    implicit none
    real(dp), dimension(:,:,:), allocatable :: mats
    integer(i4b)                            :: ndet, i
    real(dp)                                :: theta, phi
    ndet = get_num_dets()
    allocate(mats(3,3,0:ndet))
    call compute_euler_matrix_zyz(-pi, 0d0, pi, mats(:,:,0))
    do i = 1, ndet
       call get_detector_pos(i, theta, phi)
       call compute_euler_matrix_zyz(phi-pi, theta, -phi+pi, mats(:,:,i))
    end do
  end subroutine setup_detector_euler_matrices

  subroutine verify_detector(det_id)
    implicit none
    integer(i4b), intent(in) :: det_id
    if (det_id < 0 .or. det_id > get_num_dets()) then
       write(*,*) 'comap_pointing_mod: ERROR -- requested module does not exist:', det_id
       stop
    end if
  end subroutine verify_detector

  ! Hack: Let pointing mod know that pointing information might have changed
  subroutine update_pointing_mod
    implicit none
    deallocate(diode_euler)
    call setup_detector_euler_matrices(diode_euler)
  end subroutine

end module comap_pointing_mod
