program cov_unf2hdf
  use healpix_types
  use quiet_fileutils
  implicit none

  integer(i4b)       :: i, n, ordering, unit, pol_flag
  integer(i8b)       :: m
  logical(lgt)       :: inv
  character(len=128) :: icovfile, ocovfile, imapfile
  real(dp),     allocatable :: covmat(:,:), map(:,:)
  integer(i4b), allocatable :: pixels(:)
  integer(i4b)       :: nside, order, ncomp, order2

  if(iargc() == 2) then
     call getarg(1, icovfile)
     call getarg(2, ocovfile)
     call read_covmat(covmat, pixels, nside, order, ncomp, icovfile, inv, verbose=.true.)
     call write_covmat(covmat, pixels, nside, order, ncomp, ocovfile, inv, verbose=.true.)
  elseif(iargc() == 3) then
     call getarg(1, icovfile)
     call getarg(2, imapfile)
     call getarg(3, ocovfile)

     ! Read map to get pixels
     call read_map(map, pixels, nside, order, imapfile)
     ! Read cov manually
     unit = getlun()
     open(unit, file=icovfile, form='unformatted', action="read", status="old")
     read(unit) n
     read(unit) order2
     call assert(order == order2, "Map and covmat do not agree on pixel ordering")
     read(unit) ncomp
     m = n
     allocate(covmat(m,m))
     do i = 1, n
        read(unit) covmat(:,i)
     end do
     read(unit) inv
     close(unit)
     ! Write to new format
     call write_covmat(covmat, pixels, nside, order, ncomp, ocovfile, inv, verbose=.true.)
  else
     write(*,*) "Usage input_cov.unf [input_map.fits] output_cov.hdf"
     stop
  end if
end program
