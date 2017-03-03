! Read a covariance matrix in wmap fits format, and produce one in
! our format.
program fits2cov
  use quiet_utils
  use quiet_fileutils
  implicit none
  character(len=512) :: arg, ifile, ofile
  character(len=80)  :: comment, order_str
  logical(lgt)       :: anyf, inv
  integer(i4b)       :: i, j, k, m, n, unit, bs, status = 0, naxes(256)
  integer(i4b)       :: naxis, bitpix, pcount, gcount, extended, simple
  integer(i4b)       :: nside, npix, ncomp, order
  integer(i4b), dimension(:),       allocatable :: pix
  real(dp),     dimension(:,:,:,:), allocatable :: cov
  call getarg(1, ifile)
  call getarg(2, ofile)
  inv = .true. ! wmap matrix is inverse.

  call ftopen(unit, trim(ifile), 0, bs, status)
  call assert(status == 0, "Error opening file: " // trim(itoa(status)))
  call ftghpr(unit, size(naxes), simple, bitpix, naxis, naxes, pcount, &
   & gcount, extended, status)
  call assert(status == 0, "Error reading header: " // trim(itoa(status)))
  call assert(naxis == 2, "Got " // trim(itoa(naxis)) // " dims, expected 2!")
  call assert(naxes(1) == naxes(2), "Non-square matrix!")
  call ftgkyj(unit, "NSIDE", nside, comment, status)
  call assert(status == 0, "Error reading nside: " // trim(itoa(status)))
  npix  = 12*nside**2
  n     = naxes(1)
  ncomp = n/npix
  call assert(ncomp*npix == n, "ncomp*npix /= n!")
  call ftgkys(unit, "ORDERING", order_str, comment, status)
  call assert(status == 0, "Error reading ORDERING!")
  if(order_str == "RING") then; order = RING; else; order = NEST; end if

  allocate(cov(npix,ncomp,npix,ncomp))
  allocate(pix(npix))

  ! set up pixels: Full sky coverage
  do i = 1, npix; pix(i) = i-1; end do
  ! Read the data into the covariance matrix.
  call ftgpvd(unit, 1, 1, size(cov), 0, cov, anyf, status)
  call assert(status == 0, "Error raeding matrix from file: " // trim(itoa(status)))

  ! We have all we need; time to output
  call write_covmat(cov, pix, nside, order, ofile, inv)
  deallocate(pix,cov)
  call ftclos(unit, status)
end program
