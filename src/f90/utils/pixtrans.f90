! Simple program for transforming between pixel orderings
program pixtrans
  use quiet_utils
  implicit none
  character(len=512) :: arg
  integer(i4b)       :: nside, iorder, oorder, ip, op
  call getarg(1, arg); read(arg,*) nside
  call getarg(2, arg); if(arg == "ring") then; iorder = ring; else; iorder = nest; endif
  call getarg(3, arg); if(arg == "ring") then; oorder = ring; else; oorder = nest; endif
  do
     read(*,*,end=1) ip
     call set_ordering_pix(nside, iorder, oorder, ip, op)
     write(*,*) op
  end do
  1 continue
end program
