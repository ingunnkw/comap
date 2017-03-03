! Given a 'multimask' a mask file with 0 for nothing, and regions with
! 1, 2, ... n, produce folders for each region, with a sparsity for
! each. Only the first component of the mask is considered.

program region_split
  use quiet_utils
  use quiet_fileutils
  implicit none
  character(len=512)                        :: arg, infile, oprefix, oname
  real(dp),     dimension(:,:), allocatable :: map
  integer(i4b), dimension(:),   allocatable :: pixels, vals, opix
  integer(i4b) :: nside, ordering, i, j, k, n, nreg
  call getarg(1, infile)
  call getarg(2, oprefix)

  call read_map(map, pixels, nside, ordering, infile)
  allocate(vals(size(pixels)))
  vals = nint(map(:,1))
  deallocate(map)
  ! First find out which regions are defined. We do not support
  ! gaps in regions. If region 7 is defined, 1-6 must also exist.
  nreg = maxval(vals)
  ! Then, for each region, get the corresponding pixels
  do i = 1, nreg
     n = count(vals == i)
     allocate(opix(n))
     k = 0
     do j = 1, size(pixels)
        if(vals(j) /= i) cycle
        k = k+1
        opix(k) = pixels(j)
     end do
     allocate(map(n,3))
     map = 1
     oname = trim(oprefix) // trim(itoa(i)) // "/sparsity.fits"
     call mkdirs(trim(oname), .true.)
     call write_map(map, opix, nside, ordering, oname)
     deallocate(map, opix)
  end do
end program
