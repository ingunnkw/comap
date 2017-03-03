! This module provides the angular two-point correlation function
! based on an angular power spectrum, including polarization.
! The power spectrum must already be smoothed with the appropriate
! beams. The main usage of this module is:
!
!  call angcorr_init(fun, cls)
!  call angcorr_get(fun, v1, v2, C) as many times as necessary
!  call angcorr_free(fun)
!
! Here cls(0:lmax,3,3) is an angular power spectrum in matrix form,
! and C(3,3) is the covariance between the points on the sky with
! direction unit vectors v1 and v2.
!
! angcorr_get uses interpolation by default. Ideally, the spline
! nodes would have a density profile adapted to the local feature
! size, but I haven't found a satisfactory formula for that yet
! (though 1/(x+a) wasn't that bad).
module quiet_angcorr_mod
  use spline_1D_mod
  implicit none

  type angcorr_type
    real(dp), allocatable :: cls(:,:,:), coeff(:,:,:,:,:), x(:,:)
    real(dp)              :: zmax, zmin, zmid
    integer(i4b)          :: nspline, lmin
  end type

  interface angcorr_get
    module procedure angcorr_get_norot, angcorr_get_rot
  end interface

contains

  ! Set up the angular correlation function. The default is to
  ! spline for angles > 2 degrees, with 10000 spline knots.
  ! This can be overridden using the optional parameters.
  ! Setting nspline = 0 disables splining.
  subroutine angcorr_init(fun, cls, nspline, lmin, zmax, zmin, rmin, rmax)
    implicit none
    type(angcorr_type), intent(inout) :: fun
    real(dp),           intent(in)    :: cls(0:,:,:)
    real(dp),           optional      :: rmin, zmax, zmin, rmax
    integer(i4b),       optional      :: nspline, lmin
    integer(i4b)                      :: lmax, i, j, k
    real(dp)                          :: ol
    call angcorr_free(fun)
    fun%nspline = 10000; fun%lmin = 0; fun%zmin = -1; fun%zmax = cos(0.01)
    fun%zmid    = cos(17d0*pi/180)
    if(present(nspline)) fun%nspline = nspline
    if(present(lmin))    fun%lmin    = lmin
    if(present(zmax))    fun%zmax    = zmax
    if(present(zmin))    fun%zmin    = zmin
    if(present(rmin))    fun%zmax    = cos(rmin)
    if(present(rmax))    fun%zmin    = cos(rmax)
    ol = 0.05 ! Overlap between spline regions to handle bad behavoir there

    lmax = ubound(cls,1)
    allocate(fun%cls(0:lmax,3,3), fun%coeff(fun%nspline,3,3,2,2), fun%x(fun%nspline,2))
    fun%cls = cls
    if(fun%nspline > 0) then
       do i = 1, fun%nspline
          fun%x(i,1) = fun%zmid + (min(1d0,fun%zmax+ol)-fun%zmid) * real(i-1,dp) / (fun%nspline-1)
          fun%x(i,2) = fun%zmin + (fun%zmid+ol-fun%zmin) * real(i-1,dp) / (fun%nspline-1)
          call angcorr_get_raw(fun, fun%x(i,1), fun%coeff(i,:,:,1,1))
          call angcorr_get_raw(fun, fun%x(i,2), fun%coeff(i,:,:,1,2))
       end do
       do k = 1, 2
          do j = 1, size(fun%coeff,3)
             do i = j, size(fun%coeff,2)
                call spline(fun%x(:,k), fun%coeff(:,i,j,1,k), 1d30, 1d30, fun%coeff(:,i,j,2,k))
                if(i /= j) fun%coeff(:,j,i,2,k) = fun%coeff(:,i,j,2,k)
             end do
          end do
       end do
    end if
  end subroutine

  subroutine angcorr_get_raw(fun, z, C)
    implicit none
    type(angcorr_type) :: fun
    real(dp)           :: z, C(3,3)
    call compute_corrfuncs(fun%cls, z, C, fun%lmin)
  end subroutine

  subroutine angcorr_get_norot(fun, z, C)
    implicit none
    type(angcorr_type) :: fun
    real(dp)           :: z, C(3,3)
    integer(i4b)       :: i, j, k
    if(fun%nspline > 0 .and. z < fun%zmax .and. z > fun%zmin) then
       k = 1; if (z < fun%zmid) k = 2
       do j = 1, size(C,2)
          do i = j, size(C,1)
             C(i,j) = splint_uniform_grid(fun%x(:,k), fun%coeff(:,i,j,1,k), fun%coeff(:,i,j,2,k), z)
             if(i /= j) C(j,i) = C(i,j)
          end do
       end do
    else
       call angcorr_get_raw(fun, z, C)
    end if
  end subroutine

  subroutine angcorr_get_rot(fun, v1, v2, C)
    implicit none
    type(angcorr_type) :: fun
    real(dp)           :: v1(3), v2(3), C(3,3), t2a_1(3,3), t2a_2(3,3), norot(3,3)
    call angcorr_get_norot(fun, dot_product(v1,v2), norot)
    call compute_rotation_angle(v1, v2, t2a_1, t2a_2)
    C = matmul(t2a_1, matmul(norot, transpose(t2a_2)))
  end subroutine

  subroutine angcorr_free(fun)
    implicit none
    type(angcorr_type) :: fun
    if(allocated(fun%cls))   deallocate(fun%cls)
    if(allocated(fun%coeff)) deallocate(fun%coeff)
    if(allocated(fun%x))     deallocate(fun%x)
  end subroutine

  !----------------------------------!
  !-------- Helper functions --------!
  !----------------------------------!

  ! Given z = cos(theta) and a beam-smoothed power spectrum
  ! cls(0:lmax,ncomp,ncomp), compute the rotationally invariant
  ! correlation between points separated by an angle theta
  subroutine compute_corrfuncs(cls, z, C, lmin)
    implicit none
    real(dp),                            intent(in)           :: z
    real(dp),     dimension(0:,:,:),     intent(in)           :: cls
    real(dp),     dimension(:,:),        intent(out)          :: C
    integer(i4b), optional,              intent(in)           :: lmin
    integer(i4b) :: l, lmax, lmin_
    real(dp)     :: Q, R
    real(dp), allocatable, dimension(:,:,:) :: pls
    lmax = ubound(cls,1)
    lmin_= 0; if(present(lmin)) lmin_ = lmin
    allocate(pls(0:lmax,-1:1,-1:1))
    call compute_wigner_d2(z, pls)
    C = 0
    do l = lmin_, lmax
       Q = 0.5d0 * (pls(l,1,1) + pls(l,1,-1))
       R = 0.5d0 * (pls(l,1,1) - pls(l,1,-1))
       C(1,1) = C(1,1) + (2*l+1) * cls(l,1,1) * pls(l,0,0)
       C(1,2) = C(1,2) - (2*l+1) * cls(l,1,2) * pls(l,0,1)
       C(1,3) = C(1,3) - (2*l+1) * cls(l,1,3) * pls(l,0,1)
       C(2,2) = C(2,2) + (2*l+1) * (cls(l,2,2) * Q + cls(l,3,3) * R)
       C(2,3) = C(2,3) + (2*l+1) * cls(l,2,3) * pls(l,1,-1)
       C(3,3) = C(3,3) + (2*l+1) * (cls(l,2,2) * R + cls(l,3,3) * Q)
       C(2,1) = C(1,2)
       C(3,1) = C(1,3)
       C(3,2) = C(2,3)
    end do
    C = C / (4.d0*pi)
    deallocate(pls)
  end subroutine compute_corrfuncs

  ! Compute the IQU-equivalents of the
  ! legendre polynomials, which are given in terms
  ! of the Wigner d-matrices. Here z = cos(theta)
  subroutine compute_wigner_d2(z, pls)
    implicit none
    real(dp),                      intent(in)  :: z
    real(dp), dimension(0:,1:,1:), intent(out) :: pls
    integer(i4b) :: i, s_i1, s_i2, s1, s2, l, lmax
    real(dp)     :: rho
    real(dp), dimension(-1:1,-1:1) :: pl_m1, pl_00, pl_p1

    lmax = size(pls(:,1,1))-1

    pl_m1 = 0.d0; pl_00 = 0.d0; pl_p1 = 0.d0
    pls   = 0.d0

    ! Initialize recursions
    pl_m1(0,  0) = 1.d0
    pl_00(0,  0) = z
    pl_p1(0,  0) = 0.5d0 * (3.d0*z**2 - 1.d0)

    pl_p1( 1,  0) = 0.25d0 * sqrt(6.d0) * (1.d0+z) * (1.d0-z)
    pl_p1( 0,  1) = pl_p1(1, 0)
    pl_p1(-1,  0) = -pl_p1(1, 0)
    pl_p1( 0, -1) = -pl_p1(1, 0)

    pl_p1( 1,  1) = 0.25d0 * (1.d0+z)**2 
    pl_p1(-1, -1) = pl_p1(1, 1)
    pl_p1( 1, -1) = 0.25d0 * (1.d0-z)**2 
    pl_p1(-1,  1) = pl_p1(1,-1)

    pls(0,:,:) = pl_m1
    pls(1,:,:) = pl_00
    pls(2,:,:) = pl_p1

    pl_m1 = pl_00
    pl_00 = pl_p1

    ! Do the recursions
    do l = 2, lmax-1
       do s_i1 = -1, 1
          do s_i2 = -1, 1
             s1 = 2*s_i1
             s2 = 2*s_i2
             rho = sqrt(real(l**2 - s1**2,dp) * real(l**2 - s2**2,dp)) / real(l,dp)
             pl_p1(s_i1,s_i2) = real(2*l+1,dp) * (z - real(s1*s2,dp)/real(l*(l+1),dp)) * &
                  & pl_00(s_i1, s_i2) - rho * pl_m1(s_i1, s_i2)
             rho = sqrt(real((l+1)**2 - s1**2,dp) * real((l+1)**2 - s2**2,dp)) / &
                  & real(l+1,dp)
             pl_p1(s_i1,s_i2) = pl_p1(s_i1,s_i2) / rho
          end do
       end do
       pls(l+1,:,:) = pl_p1
       pl_m1 = pl_00
       pl_00 = pl_p1
    end do
  end subroutine

  ! Computes the IQU rotation matrix for parallel transport from
  ! vec1 to vec2 and vec2 to vec1. This is needed to compute
  ! the full correlation between two points.
  subroutine compute_rotation_angle(vec1, vec2, t2a_1, t2a_2)
    implicit none
    real(dp), dimension(3),   intent(in)  :: vec1, vec2
    real(dp), dimension(:,:), intent(out) :: t2a_1, t2a_2
    integer(i4b) :: i, j, nfield
    real(dp) :: len_u, len_v, z, cos_theta, sgn
    real(dp), dimension(3) :: u, v
    nfield = size(t2a_1(:,1))

    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       do i = 1, nfield
          do j = 1, nfield
             if (i == j) then
                t2a_1(i,j) = 1.d0
                t2a_2(i,j) = 1.d0
             else
                t2a_1(i,j) = 0.d0
                t2a_2(i,j) = 0.d0
             end if
          end do
       end do
       return
    end if

    sgn    = 1.d0
    if (vec1(1)*vec2(2)-vec1(2)*vec2(1) < 0.d0) sgn = -1.d0

    ! Rotation from vec1 to vec 2
    u         = vec1(3) * vec1 
    u(3)      = u(3) - 1.d0
    v         = vec2 - z * vec1
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    if (nfield == 1) then
       t2a_1(1,1) = 1.d0
    else if (nfield == 2) then
       t2a_1(1,1)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_1(2,2)  = t2a_1(1,1)             ! UU
       t2a_1(2,1)  = sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_1(1,2)  = -t2a_1(2,1)                                    ! QU
    else if (nfield == 3) then
       t2a_1(1,1)   = 1.d0 ! TT
       t2a_1(1,2)   = 0.d0 ! TQ
       t2a_1(1,3)   = 0.d0 ! TU
       t2a_1(2,1)   = 0.d0 ! QT
       t2a_1(3,1)   = 0.d0 ! UT
       t2a_1(2,2)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_1(3,3)  = t2a_1(2,2)             ! UU
       t2a_1(3,2)  = sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_1(2,3)  = -t2a_1(3,2)                                    ! QU
    end if

    ! Rotation from vec2 to vec 1; sgn is opposite from 1->2
    u         = vec2(3) * vec2
    u(3)      = u(3) - 1.d0
    v         = vec1 - z * vec2
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    cos_theta = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    if (nfield == 1) then
       t2a_2(1,1) = 1.d0
    else if (nfield == 2) then
       t2a_2(1,1)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_2(2,2)  = t2a_2(1,1)             ! UU
       t2a_2(2,1)  = -sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_2(1,2)  = -t2a_2(2,1)                                    ! QU
    else if (nfield == 3) then
       t2a_2(1,1)   = 1.d0 ! TT
       t2a_2(1,2)   = 0.d0 ! TQ
       t2a_2(1,3)   = 0.d0 ! TU
       t2a_2(2,1)   = 0.d0 ! QT
       t2a_2(3,1)   = 0.d0 ! UT
       t2a_2(2,2)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_2(3,3)  = t2a_2(2,2)             ! UU
       t2a_2(3,2)  = -sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_2(2,3)  = -t2a_2(3,2)                                    ! QU
    end if
  end subroutine

end module
