! This module contains routines simplifying healpix use, for example
! when accessing pixel windows and ring weights.
module quiet_healpix_mod
  use rngmod
  use quiet_system_mod
  use quiet_utils
  use math_tools
  implicit none

  type :: arrwrap
     real(dp), dimension(:,:), allocatable :: arr
  end type

  type(arrwrap), dimension(:), allocatable, target :: ringweights, pixwins

contains

  function get_hpix_pixwin(nside) result(win)
    implicit none
    integer(i4b) :: nside, order
    real(dp), dimension(:,:), pointer :: win
    if(.not. allocated(pixwins)) allocate(pixwins(0:255))
    order = nint(log(real(nside,dp))/log(2d0))
    if(.not. allocated(pixwins(order)%arr)) call read_hpix_pixwin(nside, pixwins(order)%arr)
    win => pixwins(order)%arr
  end function

  function get_hpix_ringweights(nside) result(weights)
    implicit none
    integer(i4b) :: nside, order
    real(dp), dimension(:,:), pointer :: weights
    if(.not. allocated(ringweights)) allocate(ringweights(0:255))
    order = nint(log(real(nside,dp))/log(2d0))
    if(.not. allocated(ringweights(order)%arr)) call read_hpix_ringweights(nside, ringweights(order)%arr)
    weights => ringweights(order)%arr
  end function

  subroutine read_hpix_aux(path, arr)
    implicit none
    real(dp), dimension(:,:), allocatable :: arr
    integer(i4b)       :: unit, bz, status, hdutype, nrow, ncol, i, repeat, width, datacode
    logical(lgt)       :: anyf
    character(len=512) :: path
    unit = getlun()
    status = 0
    call ftopen(unit, trim(path), 0, bz, status)
    call assert(status == 0, "Error opening " // trim(path) // "! Is $HEALPIX set correctly?")
    call ftmahd(unit, 2, hdutype, status)
    call ftgnrw(unit, nrow, status)
    call ftgncl(unit, ncol, status)
    call ftgtcl(unit, 1, datacode,repeat,width,status)
    allocate(arr(nrow*repeat, ncol))
    ! Inefficient due to fits format stupidity, but I don't know how to do it
    ! more efficiently
    do i = 1, ncol
       call ftgcvd(unit, i, 1, 1, nrow*repeat, 0, arr(:,i), anyf, status)
    end do
    call ftclos(unit, status)
  end subroutine

  subroutine read_hpix_ringweights(nside, weights)
    implicit none
    real(dp), dimension(:,:), allocatable :: weights
    integer(i4b)       :: nside
    character(len=512) :: path
    call getenv("HEALPIX", path)
    if(path == "") then
       path = "weight_ring_n" // trim(itoa(nside, 5)) // ".fits"
    else
       path = trim(path) // "/data/weight_ring_n" // trim(itoa(nside, 5)) // ".fits"
    end if
    call read_hpix_aux(path, weights)
    weights = weights + 1
  end subroutine

  subroutine read_hpix_pixwin(nside, window)
    implicit none
    real(dp), dimension(:,:), allocatable :: window
    integer(i4b)       :: nside
    character(len=512) :: path
    call getenv("HEALPIX", path)
    path = trim(path) // "/data/pixel_window_n" // trim(itoa(nside, 4)) // ".fits"
    call read_hpix_aux(path, window)
  end subroutine

  subroutine cl2alm(ps, rng, alm)
    implicit none
    real(dp)         :: ps(0:,:)
    complex(dpc)     :: alm(:,0:,0:)
    type(planck_rng) :: rng
    real(dp),     dimension(:,:), allocatable :: ps_mat, chol
    real(dp),  dimension(:,:), allocatable :: cchol, r
    integer(i4b)     :: lmax, l, i, n, m
    lmax = ubound(ps,1)
    n    = 3
    alm  = 0
    allocate(ps_mat(n,n), chol(n,n), cchol(n,n), r(n,2))
    do l = 0, lmax
       if(all(ps(l,:) == 0)) cycle
       call vec2smat(ps(l,:), ps_mat(:,:))
       call cholesky_decompose_with_mask(ps_mat(:,:), chol)
       do i = 1, n; r(i,1) = rand_gauss(rng); end do
       alm(:,l,0) = reshape(dcmplx(matmul(chol,r(:,1)),0d0),[n])
       do m = 1, l
          do i = 1, n; r(i,:) = [ rand_gauss(rng), rand_gauss(rng) ]/sqrt(2d0); end do
          alm(:,l,m) = reshape(dcmplx(matmul(chol,r(:,1:1)),matmul(chol,r(:,2:2))),[n])
       end do
    end do
    deallocate(ps_mat, chol, cchol, r)
  end subroutine

  ! Given an angular power spectrum cl, evaluate the
  ! two-point correlation function at regular intervals
  ! from 0 to pi. If x is provided, evaluate it at these
  ! distances instead.
  subroutine cl2corr(cl, corr, x)
    implicit none
    real(dp),     intent(in)           :: cl(0:)
    real(dp),     intent(out)          :: corr(:)
    real(dp),     intent(in), optional :: x(:)
    real(dp)                           :: dists(size(corr)), pl(0:size(cl)-1)
    integer(i4b)                       :: i, l
    if(present(x)) then
       dists = x
    else
       dists = pi*(irange(size(corr))-1)/(size(corr)-1)
    end if

    ! The two point correlation function is
    !  <f(x)*f(x+r)'> = 1/4pi * sum (2l+1)*Cl*Pl(cos(r))
    ! Our function for calculating Pl is normalized to 1/2pi, not
    ! 2/(2l+1), so Pl = sqrt(4pi/(2l+1)) pl, meaning that we get
    !  <f(x)*f(x+r)'> = sum sqrt((2l+1)/4pi))*Cl*pl(cos(r))
    corr = 0
    do i = 1, size(corr)
       call comp_normalised_plm(size(pl)-1, 0, dists(i), pl)
       do l = 0, size(pl)-1
          corr(i) = corr(i) + (2*l+1)**0.5*cl(l)*pl(l)
       end do
       corr(i) = corr(i) / sqrt(4*pi)
    end do
  end subroutine

  ! Healpix frequently works with vectors representing the uniqe elements of
  ! a symmetric matrix, for example TT, EE, BB, TE; TB, EB. These routines
  ! provide translation between the full matrix representation and the
  ! vetor version.
  subroutine vec2smat_map(map, n)
    implicit none
    integer(i4b) :: map(:,:), i, j, k, n
    do i = 1, n; map(i,:) = [ i, i]; end do
    k = n
    do i = 1, n
       do j = i+1, n
          k = k+1
          map(k,:) = [ i, j ]
       end do
    end do
  end subroutine

  subroutine smat2vec_map(map, n)
    implicit none
    integer(i4b) :: map(:,:), i, j, k, n
    do i = 1, n; map(i,i) = i; end do
    k = n
    do i = 1, n
       do j = i+1, n
          k = k+1
          map(i,j) = k
          map(j,i) = k
       end do
    end do
  end subroutine

  subroutine smat2vec(mat, flat)
    implicit none
    real(dp)     :: mat(:,:), flat(:)
    integer(i4b) :: n, i, j, k
    n = size(mat,2)
    k = 0
    do i = 1, n
       k = k+1
       flat(k) = mat(i,i)
    end do
    do i = 1, n
       do j = i+1, n
          k = k+1
          flat(k) = mat(j,i)
       end do
    end do
  end subroutine

  subroutine vec2smat(flat, mat)
    implicit none
    real(dp)     :: mat(:,:), flat(:)
    integer(i4b) :: n, i, j, k
    n = size(mat,2)
    k = 0
    do i = 1, n
       k = k+1
       mat(i,i) = flat(k)
    end do
    do i = 1, n
       do j = i+1, n
          k = k+1
          mat(j,i) = flat(k)
          mat(i,j) = flat(k)
       end do
    end do
  end subroutine

end module
