! Given a sigma map in focalplane coordiantes, produce a mask by
! hierarchically thresholding.
program sunmask
  use quiet_utils
  use quiet_fileutils
  implicit none
  character(len=512) :: mapname, maskname, arg
  integer(i4b)       :: i, j, k, m, n
  integer(i4b)       :: nside, ordering, npix, nstep
  real(dp)           :: threshold
  real(dp), dimension(:,:), allocatable :: map, mask
  real(dp), dimension(:),   allocatable :: tmap, tmask, fmask

  if(iargc() /= 4) then
     write(stderr,*) "Usage: sunmask sigma.fits threshold nstep mask.fits"
     stop
  end if

  call getarg(1, mapname)
  call getarg(2, arg); read(arg,*) threshold
  call getarg(3, arg); read(arg,*) nstep
  call getarg(4, maskname)

  call read_map(map, ordering, mapname, nside=nside)
  npix = 12*nside**2
  if(ordering == RING) call convert_ring2nest(nside, map)
  allocate(mask(0:npix-1,1),fmask(0:npix-1))
  ! Missing pixels will be accepted, so set their sigma 0 zero
  where(healnot(map(:,1))) map(:,1) = 0

  ! 1 accept, 0 reject for mask. For tmask and fmask it is opposite
  mask = 1

  ! Inefficient, but bah!
  do i = 0, nstep-1
     allocate(tmap(0:npix/4**i-1),tmask(0:npix/4**i-1))
     call udgrade_nest(map(:,1), nside, tmap, nside/2**i)
     tmap = tmap*2**i ! that is, sum divided by sqrt(n)
     where(tmap > threshold)
        tmask = 1
     elsewhere
        tmask = 0
     end where
     call udgrade_nest(tmask, nside/2**i, fmask, nside)
     where(fmask /= 0) mask(:,1) = 0
     deallocate(tmap, tmask)
  end do
  if(ordering == RING) call convert_nest2ring(nside, mask)
  call write_map(mask, ordering, maskname)
  deallocate(map, mask, fmask)

end program
