! Given an equation set and a compatible raw rhs, adds the raw rhs at the end of the
! existing rhs part, producing a new file.

program addrhs_eqn
  use quiet_utils
  use quiet_fileutils
  implicit none
  character(len=512) :: ifile, rhsfile, ofile
  integer(i4b)       :: ntot(2), ncomp(2), nside, order(2), input, output, rhs, i, inmap, rnmap
  integer(i4b)       :: npix
  logical(lgt)       :: inverse
  real(dp),     dimension(:),   allocatable :: col, map
  real(sp),     dimension(:),   allocatable :: fmap
  integer(i4b), dimension(:),   allocatable :: map2mask

  if(iargc() < 3) call error("Usage: addrhs_eqn eqn_in.unf rawrhs_in.unf eqn_out.unf")

  call getarg(1, ifile)
  call getarg(2, rhsfile)
  call getarg(3, ofile)
  call dset(level=1)
  call dmem("init")

  input  = getlun()
  open(input, file=ifile,  form="unformatted",action="read",status="old")
  call dmem("opened " // trim(ifile))
  rhs    = getlun()
  open(rhs,   file=rhsfile,form="unformatted",action="read",status="old")
  call dmem("opened " // trim(rhsfile))
  output = getlun()
  open(output,file=ofile,  form="unformatted",action="write")
  call dmem("opened " // trim(ofile))

  ! Read basic parameters
  read(input) ntot(1);  write(output) ntot(1)
  read(input) order(1); write(output) order(1)
  read(input) ncomp(1); write(output) ncomp(1)
  call dmem("eqn info: ntot: " // trim(itoa(ntot(1))) // ", order: " // &
   & trim(itoa(order(1))) // ", ncomp: " // trim(itoa(ncomp(1))))

  read(rhs)   ntot(2)
  read(rhs)   order(2)
  read(rhs)   ncomp(2)
  read(rhs)   rnmap
  call dmem("rhs info: ntot: " // trim(itoa(ntot(2))) // ", order: " // &
   &trim(itoa(order(2))) // ", ncomp: " // trim(itoa(ncomp(2))) // ", nmap: " // &
   & trim(itoa(rnmap)))

  if(any([ntot(1),order(1),ncomp(1)]/=[ntot(2),order(2),ncomp(2)])) &
   & call error("Incompatible parameters in eqn file and rawrhs file!")

  allocate(col(ntot(1)+1),map(ntot(1)),fmap(ntot(1)))
  call dmem("allocate")
  do i = 1, ntot(1)
     read(input) col; write(output) col
  end do
  call dmem("cov done")
  read(input) inmap;  write(output) inmap+rnmap
  do i = 1, inmap; read(input) map; write(output) map; end do
  call dmem("old rhs")
  do i = 1, rnmap; read(rhs)  fmap; write(output) real(fmap,dp); end do
  call dmem("new rhs")

  read(input) npix;      write(output) npix
  call dmem("npix: " // trim(itoa(npix)))
  allocate(map2mask(0:npix-1))
  call dmem("allocate map2mask")
  read(input) map2mask;  write(output) map2mask
  call dmem("done")

  close(input)
  close(output)
  close(rhs)

  deallocate(col, map, fmap, map2mask)

contains

  subroutine error(msg)
    implicit none
    character(len=*) msg
    write(*,*) msg
    stop
  end subroutine

end program
