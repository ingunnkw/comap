! Simply adds pixel info at the end of a covariance matrix, thus updating to the
! new format.

program addpix_cov
  use quiet_utils
  use quiet_fileutils
  implicit none
  character(len=512) :: ifile, mapfile, ofile, type
  integer(i4b)       :: ntot, ncomp, nside, order, morder, n, input, output, i
  logical(lgt)       :: inverse
  integer(i4b), dimension(:),   allocatable :: pixels
  real(dp),     dimension(:),   allocatable :: col
  real(dp),     dimension(:,:), allocatable :: map

  if(iargc() < 4) then
     write(*,*) "Usage: addpix_cov type input_cov.unf sparse_map.fits output_cov.unf"
     write(*,*) "type can be 'map', for a sparse map, or map2mask for a map2mask."
     write(*,*) "This affects how the third argument is interpreted."
     stop
  end if

  call getarg(1, type)
  call getarg(2, ifile)
  call getarg(3, mapfile)
  call getarg(4, ofile)

  if(type == "map2mask") then
     call read_map(map, morder, mapfile)
     allocate(pixels(nint(maxval(map(:,1)))))
     do i = 0, size(map,1)-1
        if(map(i,1) <= 0) cycle
        pixels(nint(map(i,1))) = i
     end do
     nside = npix2nside(size(map,1))
  else
     call read_map(map, pixels, nside, morder, mapfile)
  end if

  input  = getlun()
  open(input, file=ifile,form="unformatted",action="read",status="old")
  output = getlun()
  open(output,file=ofile,form="unformatted",action="write")

  read(input) ntot;  write(output) ntot
  read(input) order; write(output) order
  read(input) ncomp; write(output) ncomp

  n = ntot/ncomp
  if(n /= size(pixels)) call error("Different number of pixels in map and cov: " // &
   & trim(itoa(size(pixels))) // " vs. " // trim(itoa(n)) // "!")

  allocate(col(ntot))
  do i = 1, ntot
     read(input) col; write(output) col
  end do
  inverse = .false.
  read(input,end=1) inverse
1 write(output)     inverse
  close(input)

  do i = 1, n
     call set_ordering(nside, morder, order, pixels(i), pixels(i))
  end do

  write(output) nside
  write(output) pixels
  close(output)

  deallocate(pixels, map, col)

contains

  subroutine error(msg)
    implicit none
    character(len=*) msg
    write(*,*) msg
    stop
  end subroutine

end program
