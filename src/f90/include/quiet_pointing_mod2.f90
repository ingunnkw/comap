! The goal of this module is to provide the euler matrices making up
! the chain of transformations from detector coordinates all the way
! to galactic coordinates.
!
! * Detector coordinates: The local coordinate system of each detector.
! * Focalplane coordinates: The coordinates the detector positions are
!   specified in. If you want a map of detector position, you would use
!   these coordinates.
! * ??? coordinates: As above, but with the deck rotation axis as
!   center.
! * Boresight-local coordinates: Coordinate system pointing straight up relative to
!   the elevation drive top and with north direction given by local ideal azimuth.
!   As above, but with deck rotation (enc + offset) applied.
! * El coordinates: Coordinate system pointing straigt up relative to
!   the elevation drive bottom (i.e. as above, but with 90 - (el enc + offset) applied
!   around y axis).
! * Etilt coordinates: Coordinate system pointing straight up relative to
!   the azimuth plane.
! * Az coordinates: coordinate system with az = 0 corresponding to north. I.e. as
!   above, but with -az (enc + offset) rotation applied around z axis.
! * Atilt coordinates1: Coordinate system pointing straight up relative to
!   the local direction of gravity.
! * Flex coordinates: Coordinates where the drooping of the telescope due to
!   gravity has been taken into account.
! * Atilt coordinates2: Coordinate system pointing straight up relative to
!   the direction from the center of the earth.
! * Horizontal coordinates: Coordinates where the refraction of the atmosphere
!   has been taken into account.
! * Apparent equatorial coordinates: Equatorial coordinates before taking into
!   account the earth's motion and orbital changes.
! * Equatorial coordinates: Epoch 2000 equatorial equatorial coordinates.
! * Galactic coordinates: Epoch 2000 galactic coordinates.
!
! Each of these can be implemented as an Euler rotation, though some of them
! could be more efficiently implemented than that. Several of them are
! nonlinear in the sense that the magnitude of the rotation depends on the
! current position. Of course, a rotation applies to the whole sky, so what
! do I mean by "current position"? Up until the point where the light leaves
! the telescope, the direction needed is that of the telescope tube itself,
! i.e. the coordinates of the boresight at the current step in the chain.
! After the light leaves the telescope, the "current position" is the
! actual direction of each point on the sky - these are lensing effects that
! distort rather than just rotate.
!
! How to do this in pratice? We may wish to go from, say, Atilt coordinates 1
! to Boresight-local coordinates for example. If so, the actual boresight az
! and el are not known, only that for the detector we look at. The first 4
! systems are detector-specific. The next 6 operate on the boresight. And the
! rest operate on the combination. We can divide these into groups D, B and C,
! such that C*B*D = galactic detector pointing, C*B = galactic boresight.
!
! So how do we move backwards from a given galactic pointing to
! detector coordinates, for example? If we know that the coordinates we
! have correspond to a given detector, we can calculate C, and then do:
! gal = C*B*D*I => B = C'*gal*D', and knowing B we can extract the boresight
! az/el and thus find anything we want.
!
! For general coordinates, though, the I is replaced with a general position in
! detector-local coordinates, and it is impossible to extract the az/el.
! This is not a big problem, though, as we are not that interested in
! coordinates that are not hit by detectors.
!
! To be exact, when calculating detector pointings in Rigid coordinates,
! we must first start from boresight pointing, pass it through the rotations
! to determine the parameters for the B rotation, and then use that for
! each of the individual detectors. This might not be so bad, actually,
! as the composite B + part of D can then be stored and reused for all the
! detectors. For each detector, the only nonlinear part will then be
! the apparent pointing correction.

module quiet_pointing_mod2
  use quiet_utils
  use quiet_detector_mod
  use quiet_ces_mod
  use quiet_hdf_mod
  implicit none

  real(dp), parameter :: QUIET_GEODETIC_LONGITUDE  = -67.76166667d0 ! -67 45'42.0"
  real(dp), parameter :: QUIET_GEODETIC_LATITUDE   = -23.02822222d0 ! -23 01'41.6"
  real(dp), parameter :: QUIET_ALTITUDE            =  5020.d0       ! m
  ! Precomputed euler angles for equ->gal, J2000
  real(dp), parameter :: J2000_EULER_ALPHA         =  -57.068351386d0
  real(dp), parameter :: J2000_EULER_BETA          =  -62.871663896d0
  real(dp), parameter :: J2000_EULER_GAMMA         = -192.859498564d0
  real(dp), parameter :: LEAPS                     = 34 ! leap seconds since ??
  real(dp), parameter :: DAYSECS                   = 24d0*60*60

  integer(i4b), parameter :: mpar_denc(3) = [1,2,3], mpar_col(2) = [4,5], &
   & mpar_atilt(2) = [6,7], mpar_etilt = 8, mpar_flex = 9, mpar_max = 9

  integer(i4b), parameter :: nrdef_max = 32, rchain_max = 16
  integer(i4b), parameter :: &
   & rdef_di      =  1, rdef_horn    =  2, rdef_col     =  3, rdef_dk      =  4, &
   & rdef_el      =  5, rdef_etilt   =  6, rdef_az      =  7, rdef_gtilt   =  8, &
   & rdef_flex    =  9, rdef_atilt   = 10, rdef_atm     = 11, rdef_hor2app = 12, &
   & rdef_app2equ = 13, rdef_equ2gal = 14, nrdef_raw    = 14
  integer(i4b) :: rdef_det, rdef_bore, rdef_hor2equ, rdef_hor2gal

  real(dp), dimension(3,3), parameter, private :: identity3 = &
   & reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[3,3])

  ! Time-dependent mount model
  type mount_model
     real(dp),           allocatable :: range(:,:), params(:,:)
  end type

  type rot_chain
     integer(i4b) :: n = 0, parts(rchain_max)
  end type

  ! Everything needed for calculating a full rotation matrix.
  type rotinfo
     type(mount_model),  pointer :: model
     type(mount_model)           :: simple
     real(dp)                    :: mjd, enc(3)
     integer(i4b)                :: mod, di
     type(rot_chain)             :: rots(nrdef_max)
     real(dp)                    :: mats(3,3,nrdef_max)
     ! Handle the time-dependent external lookup
     integer(i4b)                :: mount_idx, cnum
     real(dp), allocatable       :: tilts(:,:), tilt_mjds(:)
  end type

  type(mount_model), target:: mount_null, mount_default

  type(rotinfo), private :: info

contains

  subroutine init_quiet_pointing_mod2(parfile)
    implicit none
    character(len=*) :: parfile
    call initialize_ces_mod(parfile)
    call init_detector_mod(parfile)
    allocate(mount_null%range(2,1), mount_null%params(mpar_max,1))
    mount_null%range(:,1) = [-infinity,infinity]
    mount_null%params = 0
    call init_rotinfo(info)
  end subroutine

  subroutine init_rotinfo(info)
    implicit none
    type(rotinfo), target :: info
    call free_rotinfo(info)
    rdef_det = new_rot(info, [rdef_di, rdef_horn, rdef_col])
    rdef_bore= new_rot(info, [rdef_dk, 0, rdef_el, rdef_etilt, rdef_az, rdef_gtilt, &
     & rdef_flex, rdef_atilt])
    rdef_hor2equ = new_rot(info, [rdef_atm,rdef_hor2app,rdef_app2equ])
    rdef_hor2gal = new_rot(info, [rdef_atm,rdef_hor2app,rdef_app2equ,rdef_equ2gal])
    allocate(info%simple%range(2,1),info%simple%params(mpar_max,1))
    info%simple%range(:,1) = [-infinity,infinity]
    info%simple%params = 0
    info%model => info%simple
    ! Tilt stuff
    info%cnum = 0
    allocate(info%tilts(2,size(ces_db%ceses)), info%tilt_mjds(size(ces_db%ceses)))
    info%tilts     = nan ! nan indicates that this nees to be read from file
    info%tilt_mjds = ces_db%ceses(ces_db%cidmap(cid_sort))%mjd(1)
  end subroutine

  function calc_rot(info, rots) result(mat)
    implicit none
    type(rotinfo), intent(inout) :: info
    integer(i4b),  intent(in)    :: rots(:)
    real(dp)                     :: mat(3,3)
    call calc_chain(info, rots, identity3, mat)
  end function

  subroutine free_rotinfo(info)
    implicit none
    type(rotinfo) :: info
    nullify(info%model) ! We do not clean it up - it belongs to the caller
    if(allocated(info%tilts))     deallocate(info%tilts)
    if(allocated(info%tilt_mjds)) deallocate(info%tilt_mjds)
    call free_mount_model(info%simple)
    info%rots%n = 0
  end subroutine

  subroutine free_mount_model(model)
    implicit none
    type(mount_model) :: model
    if(allocated(model%range)) deallocate(model%range)
    if(allocated(model%params))deallocate(model%params)
  end subroutine

  function get_mount_model(info) result(res)
    implicit none
    type(rotinfo),     target  :: info
    type(mount_model), pointer :: res
    res => info%model
  end function

  subroutine set_mount_model(info, model)
    implicit none
    type(rotinfo)             :: info
    type(mount_model), target :: model
    info%model => model
  end subroutine

  subroutine set_mount_simple(info, params)
    implicit none
    type(rotinfo),     target :: info
    real(dp)                  :: params(:)
    info%simple%params(:,1) = params
    info%model => info%simple
  end subroutine

  ! In these functions 'maz' and 'mel' are actually zentih angle
  ! and -az in order to make the euler matrices nicer.

  ! The following functions are the parts of the encoder-independent
  ! matrix "D" discussed above. The full D = dk*col*horn*di.

  function rot_di(di) result(mat)
    implicit none
    integer(i4b) :: di
    real(dp)     :: mat(3,3)
    call compute_euler_matrix_zyz(0d0, 0d0, quiet_diodes(di)%psi, mat)
  end function

  function rot_horn(mod) result(mat)
    implicit none
    integer(i4b) :: mod
    real(dp)     :: mat(3,3), phi, theta
    ! theta,phi are strange, but I think it has to do with whether
    ! they are defined looking up or down on the focalplane.
    phi = quiet_horns(mod)%phi; theta = quiet_horns(mod)%theta
    call compute_euler_matrix_zyz(phi, -theta, -phi, mat)
  end function

  function rot_col(phi_c, theta_c) result(mat)
    implicit none
    real(dp) :: phi_c, theta_c, mat(3,3)
    ! Not sure why we rotate in the opposite direction here.
    ! Are we sure we want to do this?
    call compute_euler_matrix_zyz(phi_c, -theta_c, -phi_c, mat)
  end function

  function rot_dk(dk,ddk) result(mat)
    implicit none
    real(dp) :: dk, ddk, mat(3,3)
    call compute_euler_matrix_zyz(0d0, 0d0, dk-ddk, mat)
  end function

  ! The following rotations implement the boresight-related
  ! pointing of the telescope, called "B" above.
  ! B = atilt(static)*flex*atilt(dyn)*az*etilt*el

  function rot_el(mel,dmel) result(mat)
    implicit none
    real(dp) :: mel, dmel, mat(3,3)
    call compute_euler_matrix_zyz(0d0, mel-dmel, 0d0, mat)
  end function

  function rot_etilt(etilt) result(mat)
    implicit none
    real(dp) :: etilt, mat(3,3)
    call compute_euler_matrix_zyz(-pi/2, etilt, pi/2, mat)
  end function

  function rot_az(maz,dmaz) result(mat)
    implicit none
    real(dp) :: maz, dmaz, mat(3,3)
    call compute_euler_matrix_zyz(0d0, 0d0, maz-dmaz, mat)
  end function

  function rot_atilt(omega, theta, maz) result(mat)
    implicit none
    real(dp) :: omega, theta, maz, mat(3,3)
    call compute_euler_matrix_zyz(-omega+maz, theta, omega-maz, mat)
  end function

  function rot_flex(maz_grav, mel_grav, kf) result(mat)
    implicit none
    real(dp) :: maz_grav, mel_grav, kf, mat(3,3)
    call compute_euler_matrix_zyz(maz_grav, kf*sin(mel_grav), -maz_grav, mat)
  end function

  ! Implementation of what happens outside the telescope, called
  ! matrix "C" above. When going all the way to galactic coordinates,
  ! C = equ2gal * app2equ * hor2app * atm

  function rot_atm(maz_atm, mel_atm, A, B) result(mat)
    implicit none
    real(dp) :: maz_atm, mel_atm, A, B, mat(3,3)
    call compute_euler_matrix_zyz(maz_atm, A*tan(mel_atm) + B*tan(mel_atm)**3, -maz_atm, mat)
  end function

  function rot_hor2app(mjd) result(mat)
    implicit none
    real(dp) :: mat(3,3), mjd, lst
    lst = mjd2lst(mjd, QUIET_GEODETIC_LONGITUDE)
    call compute_euler_matrix_zyz(DEG2RAD * (lst - 180.d0), &
     & DEG2RAD * (QUIET_GEODETIC_LATITUDE - 90.d0), 0.d0, mat)
  end function

  function rot_app2equ(mra_vac, mdec_vac, mjd) result(mat)
    implicit none
    real(dp) :: mra_vac, mdec_vac, mjd, mat(3,3), ra1, dec1, tjd, ra2, dec2, phi, theta
    ! This one does not deal properly with the psi angle. It assumes that
    ! the rotation keep shte "up" direction unchanged.
    ra1 = 12/pi*mra_vac; dec1 = 90-180/pi*mdec_vac
    tjd = mjd + (LEAPS + 32.184d0)/DAYSECS + 2400000.5d0
    call mpstar(tjd, 3, ra1, dec1, ra2, dec2)
    phi = pi/12*ra2; theta = pi/2 - pi/180*dec2
    call compute_euler_matrix_zyz(phi, theta-mdec_vac, -mra_vac, mat)
  end function

  function rot_equ2gal() result(mat)
    implicit none
    real(dp) :: mat(3,3)
    call compute_euler_matrix_zyz(DEG2RAD * J2000_EULER_ALPHA, &
     & DEG2RAD * J2000_EULER_BETA, DEG2RAD * J2000_EULER_GAMMA, mat)
  end function


  ! ======= Support routines ========

  ! Defines a new composite rotation
  function new_rot(info, parts) result(id)
    implicit none
    type(rotinfo), intent(inout) :: info
    integer(i4b),  intent(in)    :: parts(:)
    integer(i4b)                 :: id
    do id = nrdef_raw+1, nrdef_max
       if(info%rots(id)%n == 0) exit
    end do
    call assert(id < nrdef_max, "Ran out of rotation definition slots!")
    info%rots(id)%n = size(parts)
    info%rots(id)%parts(1:size(parts)) = parts
  end function

  subroutine free_rot(info, id)
    implicit none
    type(rotinfo), intent(inout) :: info
    integer(i4b),  intent(in)    :: id
    info%rots(id)%n = 0
  end subroutine

  ! Calculate the values of the rotations given. Does not try
  ! to be smart - no caching is going on here (though it is
  ! very tempting). Recursively invokes itself for user-defined chains.
  ! Use negative ids to reuse the previous value.
  recursive subroutine calc_chain(info, rots, imat, omat)
    implicit none
    type(rotinfo), intent(inout), target :: info
    integer(i4b),  intent(in)    :: rots(:)
    real(dp),      intent(in)    :: imat(3,3)
    real(dp),      intent(out)   :: omat(3,3)
    real(dp)                     :: mat(3,3), cmat(3,3), phi, theta, psi, tilt(2), refr(2)
    integer(i4b)                 :: i, mi
    real(dp),           pointer  :: model(:)
    cmat  = imat
    omat  = identity3
    model => get_mount_params(info)

    do i = 1, size(rots)
       if(rots(i) < 0) then
          cmat = matmul(info%mats(:,:,-rots(i)),cmat)
          omat = matmul(info%mats(:,:,-rots(i)),omat)
       elseif(rots(i) == 0) then
          cmat = identity3
       else
          select case(rots(i))
             case(rdef_di);   mat  = rot_di(info%di)
             case(rdef_horn); mat  = rot_horn(info%mod)
             case(rdef_col);  mat  = rot_col(model(mpar_col(1)), model(mpar_col(2)))
             case(rdef_dk);   mat  = rot_dk(info%enc(3), model(mpar_denc(3)))
             case(rdef_el);   mat  = rot_el(info%enc(2), model(mpar_denc(2)))
             case(rdef_etilt);mat  = rot_etilt(model(mpar_etilt))
             case(rdef_az);   mat  = rot_az(info%enc(1), model(mpar_denc(1)))
             case(rdef_gtilt)
                tilt = get_tilt_params(info)
                call convert_euler_matrix_to_angles_zyz(cmat, phi, theta, psi)
                mat = rot_atilt(tilt(1), tilt(2), phi)
             case(rdef_flex)
                call convert_euler_matrix_to_angles_zyz(cmat, phi, theta, psi)
                mat = rot_flex(phi, theta, model(mpar_flex))
             case(rdef_atilt)
                mat = rot_atilt(model(mpar_atilt(1)), model(mpar_atilt(2)),0d0)
             case(rdef_atm)
                call convert_euler_matrix_to_angles_zyz(cmat, phi, theta, psi)
                call get_refraction_params(info, refr)
                mat = rot_atm(phi, theta, refr(1), refr(2))
             case(rdef_hor2app); mat = rot_hor2app(info%mjd)
             case(rdef_app2equ)
                call convert_euler_matrix_to_angles_zyz(cmat, phi, theta, psi)
                mat = rot_app2equ(phi, theta, info%mjd)
             case(rdef_equ2gal); mat = rot_equ2gal()
             case default
                call calc_chain(info, info%rots(rots(i))%parts(1:info%rots(rots(i))%n), cmat, mat)
          end select
          info%mats(:,:,rots(i)) = mat
          cmat = matmul(mat, cmat)
          omat = matmul(mat, omat)
       end if
    end do
  end subroutine

  function get_mount_params(info) result(res)
    implicit none
    type(rotinfo),      intent(inout), target  :: info
    real(dp),                          pointer :: res(:)
    integer(i4b)                               :: i
    i = info%mount_idx
    if(i < 1 .or. i > size(info%model%range,2)) i = 1
    info%mount_idx = i
    if(info%mjd<info%model%range(1,i) .or. info%mjd>info%model%range(2,i)) then
       info%mount_idx = locate(info%model%range(2,:), info%mjd)
    end if
    res => info%model%params(:,info%mount_idx)
  end function

  function get_tilt_params(info) result(res)
    implicit none
    type(rotinfo)  :: info
    real(dp)       :: res(2)
    type(hdf_file) :: hfile
    if(info%cnum <= 0 .or. info%cnum > size(ces_db%ceses)) then
       info%cnum = ces_db%cidmap(cid_sort(locate(info%tilt_mjds, info%mjd)))
    elseif(ces_db%ceses(info%cnum)%mjd(1) > info%mjd .or. &
     & ces_db%ceses(info%cnum)%mjd(2) < info%mjd) then
       ! This will trigger all the time if we are outside a ces range. May be slow
       info%cnum = ces_db%cidmap(cid_sort(locate(info%tilt_mjds, info%mjd)))
    end if
    if(isnan(info%tilts(1,info%cnum))) then
       call open_hdf_file(ces_db%ceses(info%cnum)%l3file, hfile, "r")
       call read_hdf(hfile, "tilt", slice([1,2],[1]), info%tilts(:,info%cnum))
       call close_hdf_file(hfile)
    end if
    res = info%tilts(:,info%cnum)
  end function

  subroutine get_refraction_params(info, refr)
    implicit none
    type(rotinfo), intent(inout) :: info
    real(dp),      intent(out)   :: refr(2)
    refr = 0
  end subroutine

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

end module
