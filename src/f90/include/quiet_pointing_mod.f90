module quiet_pointing_mod
  use healpix_types
  use pix_tools
  use quiet_utils
  use math_tools
  use l1_read_mod
  use quiet_module_mod
  use quiet_pmac_mod
  !use quiet_ces_mod
  use quiet_hdf_mod
  use powell_mod
  use locate_mod
  implicit none

  ! Note on coordinate systems /-conversions:
  ! Tele: Telescope pointing in horizontal coords according to encoder
  ! Tele + mount model + focal plane -> hor
  ! Hor: Individual detector pointing in horizontal coords
  ! Then we do apparent corrections and conversion to celestial and galactic coords

  real(dp) :: QUIET_GEODETIC_LONGITUDE  = -67.76166667d0  ! 67 degrees, 45 minutes, 42.0 seconds WEST
  real(dp) :: QUIET_GEODETIC_LATITUDE   = -23.02822222d0 ! 23 degrees, 1 minute, 41.6 seconds SOUTH
  real(dp) :: QUIET_ALTITUDE            =  5020.d0        ! meters

!  real(dp), parameter :: QUIET_GEODETIC_LONGITUDE  = -67.7667d0      ! Keith's longitude
!  real(dp), parameter :: QUIET_GEODETIC_LATITUDE   = -23.0333d0      ! Keith's latitude
!  real(dp), parameter :: QUIET_ALTITUDE            =  5080.d0        ! Keith's altitude

  real(dp), parameter :: QUIET_GEOCENTRIC_LATITUDE = -22.890103588d0 ! degrees
  real(dp), parameter :: QUIET_RADIUS              = 1.000277634d0   ! units of earth radii

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
  integer(i4b), parameter :: COORD_TELE = 0, COORD_HOR = 1, COORD_APP = 2, &
    & COORD_EQU = 3, COORD_CEL = 3, COORD_GAL = 4, COORD_CENTER = 5, &
    & NUM_COORD = 6
  ! The next step of the shortest path from one system to another
  integer(i4b)            :: next_step(0:NUM_COORD-1,0:NUM_COORD-1)

  type quiet_mount
     real(dp)     :: mjd_start
     real(dp)     :: enc_offset(3), kf, omega, theta, theta_E
     real(dp)     :: theta_c, psi_c, qsag, theta_o, psi_o, azcurr
     real(dp)     :: tilt_scale(2), fp_flex(2), ellcol(2), encflex(4), azcomp(4)
     real(dp)     :: akito_el(5), akito_dk(2), akito_dir
  end type quiet_mount

  integer(i4b),       parameter                            :: num_mount_par = 8
  logical(lgt),                                    private :: apparent_to_ast
  integer(i4b),                                    private :: nmod, num_mount
  type(quiet_mount),  allocatable, target,         private :: mount_data(:)
  type(quiet_mount),               target,         private :: mount_over(1)
  type(quiet_mount),               target                  :: mount_none
  type(quiet_mount),  pointer,                     private :: mount(:)
  real(dp),                        dimension(3,3), private :: M_cel2gal, M_gal2cel
  real(dp), dimension(:,:,:,:), allocatable, private       :: diode_euler

  real(dp), dimension(:,:),     allocatable, private       :: tilts ! (omega/theta,cnum)
  real(dp), dimension(:),       allocatable, private       :: mjds
  integer(i4b),                              private       :: prev_cnum

  type inv_rot_params
     integer(i4b) :: mod, di
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

  subroutine initialize_quiet_pointing_mod(paramfile, apply_mount_model, apparent_correction)
    implicit none
    character(len=*), intent(in) :: paramfile
    logical(lgt),     intent(in), optional :: apply_mount_model
    logical(lgt),     intent(in), optional :: apparent_correction

    integer(i4b)       :: i, nparam
    logical(lgt)       :: apply_mount
    type(hdf_file)     :: hfile
    character(len=512) :: focalplane, mountfile, tilt_file, mountformat
    real(dp) :: phi1, theta1, psi1, phi2, theta2, psi2, mjd
    if(initialized) return
    call initialize_simple_pointing
    call initialize_module_mod(paramfile)
    call initialize_pmac_mod(paramfile)
    !call initialize_ces_mod(paramfile)
    call setup_diode_euler_matrices(diode_euler)

    call flatten_mount(mount_none, [0d0], 0, nparam)
    call flatten_mount(mount_none, dzeroes(nparam), -1)

    !call get_parameter(0, paramfile, 'TILT_FILE', par_string=tilt_file, desc=&
    ! & "HDF file with per ces tilt values.")
    if(present(apply_mount_model)) then
       apply_mount = apply_mount_model
    else
       call get_parameter(0, paramfile, 'APPLY_MOUNT_MODEL', par_lgt=apply_mount)
    end if
    if (apply_mount) then
       call get_parameter(0, paramfile, 'MOUNT_MODEL_FILE',    par_string=mountfile)
       call get_parameter(0, paramfile, 'MOUNT_MODEL_FORMAT',  par_string=mountformat)
    end if
    if(present(apparent_correction)) then
       apparent_to_ast = apparent_correction
    else
       call get_parameter(0, paramfile, 'APPARENT_POINTING_CORRECTION', par_lgt=apparent_to_ast)
    end if

    if(apply_mount) then
       call set_mount_params(mountfile, mountformat, mount_data)
    else
       allocate(mount_data(1))
       mount_data = mount_none
    end if
    mount => mount_data

    next_step = &
      & transpose(reshape((/ &
      &  COORD_TELE, COORD_HOR, COORD_HOR, COORD_HOR, COORD_HOR, COORD_HOR, &
      &  COORD_TELE, COORD_HOR, COORD_APP, COORD_APP, COORD_APP, COORD_APP, &
      &  COORD_HOR,  COORD_HOR, COORD_APP, COORD_CEL, COORD_CEL, COORD_CEL, &
      &  COORD_APP,  COORD_APP, COORD_APP, COORD_CEL, COORD_GAL, COORD_GAL, &
      &  COORD_CEL,  COORD_CEL, COORD_CEL, COORD_CEL, COORD_GAL, COORD_CENTER, &
      &  COORD_GAL,  COORD_GAL, COORD_GAL, COORD_GAL, COORD_GAL, COORD_CENTER  &
      & /), (/NUM_COORD, NUM_COORD/)))

    ! Setup NOVAS
    call loacc
    call setdt(65.8d0) ! Should be set properly

    !prev_cnum = 0
    !call open_hdf_file(tilt_file, hfile, "r")
    !call read_alloc_hdf(hfile, "tilt", tilts)
    !call close_hdf_file(hfile)
    !allocate(mjds(size(ces_db%ceses)))
    !mjds = ces_db%ceses(ces_db%cidmap(cid_sort))%mjd(1)
    initialized = .true.
  end subroutine initialize_quiet_pointing_mod

  ! NOTE: All angles here are in the healpix convention!

  ! Driver routine for switching between coordinate systems.
  ! Note: mjd is needed for TELE and HOR to
  ! the others. mod and diode is needed to and from TELE, unless you
  ! want the borepoint.
  !
  ! The coordinate systems are connected like this:
  !
  ! center-gal-cel-app-hor-tele
  !
  ! We transform stepwise twoards the target. The setup does support
  ! branches in the tree, but we don't have any yet. The current
  ! implementation is about the same length as the old one,
  ! but this one scales as O(N), while the previous one scaled
  ! as O(N^2), where N is the number of coordinate systems.
  subroutine coord_convert(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, mod, diode, phic, thetac, euler)
    implicit none
    integer(i4b)           :: sys1, sys2, sys, nsys
    integer(i4b), optional :: mod, diode
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
       mat = angles2hor(mjd, phi1, theta1, psi1, mod, diode)
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
                case(COORD_APP);    mat = matmul(rot_hor2equ(mjd), mat)
                case(COORD_TELE);   continue
                case default; call abort("Bug in COORD_HOR conversion!")
             end select
          case(COORD_APP)
             select case(nsys)
                case(COORD_CEL);    if(apparent_to_ast) mat = app2cel(mat, mjd)
                case(COORD_HOR);    mat = matmul(rot_equ2hor(mjd), mat)
                case default; call abort("Bug in COORD_APP conversion!")
             end select
          case(COORD_CEL)
             select case(nsys)
                case(COORD_GAL);    mat = matmul(rot_equ2gal(), mat)
                case(COORD_APP);    if(apparent_to_ast) mat = cel2app(mat, mjd)
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
       call hor2angles(mat, mjd, phi2, theta2, psi2, mod, diode)
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
  function angles2hor(mjd, phi, theta, psi, module, diode) result(mat)
    implicit none
    real(dp),     intent(in)  :: phi, theta, psi, mjd
    real(dp), dimension(3,3)  :: mat
    integer(i4b)              :: mod, di
    integer(i4b), optional    :: module, diode
    mod = -1; if(present(module)) mod = module
    di  = -1; if(present(diode))  di  = diode
    mat = rot_boresight2hor(mjd, phi, theta, psi, mod, di)
  end function

  subroutine hor2angles(mat, mjd, phi, theta, psi, module, diode)
    implicit none
    real(dp)                  :: mjd, phi, theta, psi, mat(3,3)
    integer(i4b)              :: mod, di
    integer(i4b), optional    :: module, diode
    mod = -1; if(present(module)) mod = module
    di  = -1; if(present(diode))  di  = diode
    call rot_hor2boresight(mjd, mat, mod, di, phi, theta, psi)
  end subroutine

  ! The non-complicated rotations
  function rot_hor2equ(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = get_matrix_hor2equ(mjd2lst(mjd, QUIET_GEODETIC_LONGITUDE), &
      & QUIET_GEODETIC_LATITUDE)
  end function

  function rot_hor2gal(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = matmul(rot_equ2gal(),rot_hor2equ(mjd))
  end function

  function rot_equ2hor(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd
    mat = get_matrix_equ2hor(mjd2lst(mjd, QUIET_GEODETIC_LONGITUDE), &
      & QUIET_GEODETIC_LATITUDE)
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

  function app2cel(mat1, mjd) result(mat2)
    implicit none
    real(dp) :: mat1(3,3), mat2(3,3), ra1, ra2, dec1, dec2, psi, mjd, tjd
    call rot2angles(mat1, ra1, dec1, psi)
    ra1 = 12/pi*ra1; dec1 = 90-180/pi*dec1
    tjd = mjd + (LEAPS + 32.184d0)/DAYSECS + 2400000.5d0
    call mpstar(tjd, 3, ra1, dec1, ra2, dec2)
    ra2 = pi/12*ra2; dec2 = pi/2 - pi/180*dec2
    mat2 = angles2rot(ra2, dec2, psi)
  end function

  function cel2app(mat1, mjd) result(mat2)
    implicit none
    real(dp) :: mat1(3,3), mat2(3,3), ra1, ra2, dec1, dec2, psi, mjd, tjd
    call rot2angles(mat1, ra1, dec1, psi)
    ra1 = 12/pi*ra1; dec1 = 90-180/pi*dec1
    tjd = mjd + (LEAPS + 32.184d0)/DAYSECS + 2400000.5d0
    call apstar(tjd, 3, ra1, dec1, 0d0, 0d0, 0d0, 0d0, ra2, dec2)
    ra2 = pi/12*ra2; dec2 = pi/2 - pi/180*dec2
    mat2 = angles2rot(ra2, dec2, psi)
  end function

  !!!!!!!!!!!!!!!!!!!!
  ! Helper functions !
  !!!!!!!!!!!!!!!!!!!!

  ! Rotate from module coordinates to boresight. If diode is not included,
  ! the diode angle is not part of the rotation.
  function rot_module2boresight(mod, diode) result(mat)
    implicit none
    integer(i4b), intent(in)           :: mod
    integer(i4b), intent(in), optional :: diode
    real(dp), dimension(3,3)           :: mat
    integer(i4b)                       :: di
    call verify_module(mod)
    di = -1; if(present(diode)) di = diode
    mat = diode_euler(:,:,di,mod)
  end function

  function rot_boresight2module(mod, diode) result(mat)
    implicit none
    integer(i4b), intent(in)           :: mod
    integer(i4b), intent(in), optional :: diode
    real(dp), dimension(3,3)           :: mat
    mat = transpose(rot_module2boresight(mod, diode))
  end function

  ! Perform the full rotation from naive boresight coordinates
  ! to corrected horizontal coordinates.
  function rot_boresight2hor(mjd, phi, theta, psi, mod, di) result(mat)
    implicit none
    real(dp),     intent(in)  :: phi, theta, psi, mjd
    integer(i4b), intent(in)  :: mod, di
    real(dp)                  :: p(3), az, el, dk, coaz, ab(2), ab2(2), p1(2), p2(2)
    real(dp)                  :: del, coeff, dev, dpsi
    real(dp), dimension(3,3)  :: mat, M_col, M_bore, M_flex, M_az, M_el, M_fp, fp1, fp2, M_azcomp
    real(dp), dimension(3,3)  :: M_test1, M_test2, M_detector
    type(quiet_mount), pointer:: mnt
    mnt => mount(get_mount_idx(mjd))
    p   = [phi,theta,psi]-[mnt%enc_offset(1)+mnt%encflex(2)*cos(phi+mnt%encflex(1)),&
         & mnt%enc_offset(2)+mnt%encflex(3)*cos(phi+mnt%encflex(1)),&
         & mnt%enc_offset(3)+mnt%encflex(4)*cos(phi+mnt%encflex(1))]
    az  = -p(1); el = pi/2-p(2); dk = p(3)
    mat = get_identity(3)
    ! Collimation, including ellipticity
    call compute_euler_matrix_zyz(mnt%psi_c+pi, mnt%theta_c*(1+mnt%ellcol(1)*cos(dk+mnt%ellcol(2))), -mnt%psi_c+pi, M_col)
    mat = matmul(M_col, mat)
    ! Fp flex
    !coaz = quiet_horns(mod)%theta*sin(quiet_horns(mod)%phi+dk)
    !call compute_euler_matrix_zyz(-mnt%fp_flex(1), mnt%fp_flex(2)*coaz*sin(el), mnt%fp_flex(1), M_fp)
    !mat = matmul(M_fp, mat)
    ! Boresight pointing
    call compute_euler_matrix_zyz(p(1), p(2), p(3), M_bore)
    mat = matmul(M_bore, mat)
    ! Flexure
    call compute_euler_matrix_zyz(-az, mnt%kf*cos(el), az, M_flex)
    mat = matmul(M_flex, mat)
    ! Az tilt
    call compute_euler_matrix_zyz(-mnt%omega, mnt%theta, mnt%omega, M_az)
    mat = matmul(M_az, mat)
    ! El tilt
    call compute_euler_matrix_zyz(-(0.5d0*pi+az), mnt%theta_E, (0.5d0*pi+az), M_el)
    mat = matmul(M_el, mat)
    ! Az compression along collimation circle
    call compute_euler_matrix_zyz(dk-mnt%azcomp(1), mnt%azcomp(2)*cos(az+mnt%azcomp(3))*&
         & sin(0.25*(dk+mnt%azcomp(4))), -dk+mnt%azcomp(1), M_azcomp)
    mat = matmul(M_azcomp, mat)

    ! Akito stuff
!!$    if(mod < 0 .or. all(mnt%akito_el(:4) == 0) .and. mnt%akito_dk(1) == 0) then
!!$       M_detector = diode_euler(:,:,di,mod)
!!$    else
!!$       del   = (el - mnt%akito_el(5))*RAD2DEG
!!$       coeff = sum(mnt%akito_el(:4)*[1d0,del,del**2,del**3]) + &
!!$        & mnt%akito_dk(1)*sin(dk-mnt%akito_dk(2))
!!$       p1    = [quiet_horns(mod)%theta, quiet_horns(mod)%phi]
!!$       coaz  = p1(1)*sin(p1(2)+dk)
!!$       dev   = coaz*coeff*RAD2DEG
!!$       ab    = [-p1(1)*cos(p1(2)+pi/2),p1(1)*sin(p1(2)+pi/2)]
!!$       ab2   = ab - dev/60*DEG2RAD*[cos(mnt%akito_dir),-sin(mnt%akito_dir)]
!!$       p2    = [sqrt(sum(ab2**2)),atan2(ab2(2),-ab2(1))-pi/2]
!!$       dpsi  = 0
!!$       if(di >= 0) dpsi = get_diode_angle(mod,di)
!!$       call compute_euler_matrix_zyz(p2(2)-pi, p2(1), &
!!$        & dpsi-p2(2)+pi, M_detector)
!!$    end if
!!$    mat = matmul(mat, M_detector)
    mat = matmul(mat, diode_euler(:,:,di,mod))
  end function

  ! Go from horizontal boresight coordinates to telescope coordinates.
  ! We use powell. This is a pretty slow function. It may be sped up
  ! by using a better starting point by using the fiducial focalplane
  ! layout as an approximation.
  subroutine rot_hor2boresight(mjd, mat_hor, mod, di, phi, theta, psi)
    implicit none
    real(dp) :: mjd, mat_hor(3,3), phi, theta, psi, p(3)
    integer(i4b) :: i, mod, di, err
    inv_rot%mod = mod
    inv_rot%di  = di
    inv_rot%mjd = mjd
    inv_rot%M   = mat_hor
    call convert_euler_matrix_to_angles_zyz(mat_hor, p(1), p(2), p(3))
    call powell(p, powell_inv_rot, err)
    phi = p(1); theta = p(2); psi = p(3)
  end subroutine

  function powell_inv_rot(p) result(res)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: res
    res = sum(abs(rot_boresight2hor(inv_rot%mjd, p(1),p(2),p(3), inv_rot%mod, inv_rot%di)-inv_rot%M))
  end function

  function parse_coord_name(name) result(res)
    implicit none
    character(len=*) :: name
    integer(i4b)     :: res
    if(name == 'galactic' .or. name == 'gal') then; res = COORD_GAL
    elseif(name == 'celestial' .or. name == 'cel') then; res = COORD_CEL
    elseif(name == 'equatorial' .or. name == 'equ') then; res = COORD_EQU
    elseif(name == 'apparent' .or. name == 'app') then; res = COORD_APP
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

  function get_mount_idx(mjd) result (idx)
    implicit none
    real(dp),         intent(in) :: mjd
    integer(i4b)                 :: idx
    integer(i4b) :: i
    if (size(mount) == 1) then
       idx = 1
       return
    end if

    i = 1
    do while (mount(i+1)%mjd_start < mjd)
       i = i+1
       if (i == num_mount) exit
    end do
    idx = i
  end function

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

  subroutine set_mount_params(mount_file, mount_format, mount)
    implicit none
    character(len=*), intent(in)    :: mount_file, mount_format
    type(quiet_mount), dimension(:), allocatable, intent(inout) :: mount
    integer(i4b)       :: i, unit
    character(len=512) :: line

    ! Find number of valid lines
    unit = getlun()
    num_mount = 0
    open(unit, file=trim(mount_file),status='old')
    do while (.true.)
       read(unit,'(a)',end=200) line
       if (line(1:1) == '#') then
          cycle 
       else
          read(line,*,end=200) 
          num_mount = num_mount+1
       end if
    end do
200 close(unit)
    if(allocated(mount)) deallocate(mount)
    allocate(mount(num_mount))
    i = 1
    open(unit, file=trim(mount_file), status='old')
    do while (.true.)
       read(unit,'(a)',end=210) line
       if (line(1:1) == '#') then
          cycle 
       else
          ! This is the format in the mount model file. If we assume that
          ! the files always are kept updated when new parameters are added,
          ! this could be done much more simply. But to avoid making the older
          ! files unreadable whenever a new parameter is added.
          ! Parameters to add: az_offset, el_offset, tilt_scaling?
          mount(i) = mount_none
          if (trim(mount_format) == 'W') then
             read(line,*) &
                  & mount(i)%mjd_start, mount(i)%enc_offset, mount(i)%omega, mount(i)%theta, &
                  & mount(i)%theta_e,   mount(i)%theta_c,    mount(i)%psi_c, mount(i)%ellcol, &
                  & mount(i)%encflex
          else 
             read(line,*) &
                  & mount(i)%mjd_start, mount(i)%enc_offset(3), mount(i)%qsag,  &
                  & mount(i)%kf,        mount(i)%omega,         mount(i)%theta, &
                  & mount(i)%theta_e,   mount(i)%theta_c,       mount(i)%psi_c, &
                  & mount(i)%theta_o,   mount(i)%psi_o
          end if
          mount(i)%kf            = mount(i)%kf            * DTOR
          mount(i)%omega         = mount(i)%omega         * DTOR
          mount(i)%theta         = mount(i)%theta         * DTOR
          mount(i)%theta_e       = mount(i)%theta_e       * DTOR
          mount(i)%theta_c       = mount(i)%theta_c       * DTOR
          mount(i)%psi_c         = mount(i)%psi_c         * DTOR
          mount(i)%qsag          = mount(i)%qsag          * DTOR
          mount(i)%theta_o       = mount(i)%theta_o       * DTOR
          mount(i)%psi_o         = mount(i)%psi_o         * DTOR
          mount(i)%azcurr        = mount(i)%azcurr        * DTOR
          mount(i)%enc_offset    = mount(i)%enc_offset    * DTOR
          mount(i)%ellcol        = mount(i)%ellcol        * DTOR
          mount(i)%encflex       = mount(i)%encflex       * DTOR
          i = i+1
       end if
    end do
210 close(unit)
  end subroutine set_mount_params

  subroutine set_mount_override(status, mount_params, mount_model)
    implicit none
    logical(lgt)                :: status
    real(dp),          optional :: mount_params(:)
    type(quiet_mount), optional :: mount_model
    if(present(mount_params)) then
       call flatten_mount(mount_over(1), mount_params, -1)
    elseif(present(mount_model)) then
       mount_over(1) = mount_model
    end if
    if(status) then
       mount => mount_over
    else
       mount => mount_data
    end if
  end subroutine set_mount_override

  subroutine flathelp(a,b,i,dir)
    implicit none
    real(dp)     :: a, b(:)
    integer(i4b) :: dir, i
    if(dir > 0) then; b(i) = a; elseif(dir < 0) then; a = b(i); end if
  end subroutine

  subroutine flatten_mount(mount, flat, dir, n)
    implicit none
    type(quiet_mount) :: mount
    integer(i4b)      :: i, dir
    real(dp)          :: flat(:)
    integer(i4b), optional :: n
    i = 0
    i=i+1; call flathelp(mount%enc_offset(1), flat, i, dir)
    i=i+1; call flathelp(mount%enc_offset(2), flat, i, dir)
    i=i+1; call flathelp(mount%enc_offset(3), flat, i, dir)
    i=i+1; call flathelp(mount%kf,            flat, i, dir)
    i=i+1; call flathelp(mount%omega,         flat, i, dir)
    i=i+1; call flathelp(mount%theta,         flat, i, dir)
    i=i+1; call flathelp(mount%theta_e,       flat, i, dir)
    i=i+1; call flathelp(mount%theta_c,       flat, i, dir)
    i=i+1; call flathelp(mount%psi_c,         flat, i, dir)
    i=i+1; call flathelp(mount%tilt_scale(1), flat, i, dir)
    i=i+1; call flathelp(mount%tilt_scale(2), flat, i, dir)
    i=i+1; call flathelp(mount%fp_flex(1),    flat, i, dir)
    i=i+1; call flathelp(mount%fp_flex(2),    flat, i, dir)
    i=i+1; call flathelp(mount%azcurr,        flat, i, dir)
    i=i+1; call flathelp(mount%azcomp(1),     flat, i, dir)
    i=i+1; call flathelp(mount%azcomp(2),     flat, i, dir)
    i=i+1; call flathelp(mount%azcomp(3),     flat, i, dir)
    i=i+1; call flathelp(mount%azcomp(4),     flat, i, dir)
    i=i+1; call flathelp(mount%ellcol(1),     flat, i, dir)
    i=i+1; call flathelp(mount%ellcol(2),     flat, i, dir)
    i=i+1; call flathelp(mount%encflex(4),    flat, i, dir)
    i=i+1; call flathelp(mount%akito_el(1),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_el(2),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_el(3),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_el(4),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_el(5),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_dk(1),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_dk(2),   flat, i, dir)
    i=i+1; call flathelp(mount%akito_dir,     flat, i, dir)
    if(present(n)) n = i
  end subroutine

  subroutine setup_diode_euler_matrices(mats)
    implicit none
    real(dp), dimension(:,:,:,:), allocatable :: mats
    integer(i4b)                              :: nmod, ndi, m, d
    real(dp)                                  :: theta, phi
    nmod = get_num_modules()
    ndi  = get_num_diodes()
    allocate(mats(3,3,-1:ndi-1,-1:nmod-1))
    call compute_euler_matrix_zyz(-pi, 0d0, pi, mats(:,:,-1,-1))
    do m = 0, nmod-1
       call get_module_pos(m, theta, phi)
       call compute_euler_matrix_zyz(phi-pi, theta, -phi+pi, mats(:,:,-1,m))
       do d = 0, ndi-1
          call compute_euler_matrix_zyz(phi-pi, theta, &
            & get_diode_angle(m,d)-phi+pi, mats(:,:,d,m))
       end do
    end do
  end subroutine

  subroutine verify_module(mod)
    implicit none
    integer(i4b), intent(in) :: mod
    if (mod < -1 .or. mod >= get_num_modules()) then
       write(*,*) 'quiet_pointing_mod: ERROR -- requested module does not exist:', mod
       stop
    end if
  end subroutine verify_module

  ! Hack: Let pointing mod know that pointing information might have changed
  subroutine update_pointing_mod
    implicit none
    deallocate(diode_euler)
    call setup_diode_euler_matrices(diode_euler)
  end subroutine

  !function get_tilt(mjd) result(res)
  !  implicit none
  !  real(dp)       :: mjd, res(2)
  !  type(hdf_file) :: hfile
  !  if(prev_cnum <= 0 .or. prev_cnum > size(ces_db%ceses)) then
  !     prev_cnum = ces_db%cidmap(cid_sort(locate(mjds,mjd)))
  !  elseif(ces_db%ceses(prev_cnum)%mjd(1) > mjd .or. &
  !   & ces_db%ceses(prev_cnum)%mjd(2) < mjd) then
  !     ! This will trigger all the time if we are outside a ces range. May be slow
  !     prev_cnum = ces_db%cidmap(cid_sort(locate(mjds,mjd)))
  !  end if
  !  res = tilts(:,prev_cnum)
  !end function

end module quiet_pointing_mod
