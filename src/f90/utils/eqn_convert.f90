! Convert from the old to the new equation set forms
program eqn_convert
  use quiet_utils
  implicit none
  character(len=512) :: str, ifile, ofile
  integer(i4b)       :: iunit, ounit, i, j, n, order, pol, nrhs, npix
  logical(lgt)       :: reverse
  real(dp),     dimension(:), allocatable :: col
  integer(i4b), dimension(:), allocatable :: map2mask

  reverse = .false.
  i = 1
  j = 0
  do while(i <= iargc())
     call getarg(i, str)
     if(str == "-r") then; reverse = .true.
     else
        if(j == 0) then;     ifile = str
        elseif(j == 1) then; ofile = str
        end if
        j = j+1
     end if
     i = i+1
  end do
  if(j < 2) call help

  iunit = getlun()
  open(iunit,file=ifile,form="unformatted",action="read")
  ounit = getlun()
  open(ounit,file=ofile,form="unformatted",action="write")

  read (iunit) n
  write(ounit) n
  read (iunit) order
  write(ounit) order
  read (iunit) pol
  write(ounit) pol
  allocate(col(n))
  do i = 1, n
     read (iunit) col
     write(ounit) col
  end do
  if(reverse) then
     read (iunit) nrhs
  else
     nrhs = 1
     write(ounit) nrhs
  end if
  do i = 1, nrhs
     read (iunit) col
     write(ounit) col
  end do
  read (iunit) npix
  write(ounit) npix
  allocate(map2mask(0:npix-1))
  read (iunit) map2mask
  write(ounit) map2mask

  deallocate(col, map2mask)
  close(iunit)
  close(ounit)
contains
  subroutine help
     write(*,*) "Synatex: eqn_convert [options] infile outfile"
     write(*,*) "Options: -r: Convert from new to old instead of old to new"
     stop
  end subroutine
end program
