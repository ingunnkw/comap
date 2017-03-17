program time_test
  use comap_Lx_mod
  implicit none

  include "mpif.h"

  character(len=1024) :: filename
  type(lx_struct)  :: data
  real(dp)         :: start_time, stop_time, total_time
  integer(i8b)     :: filesize
  integer(i4b)     :: myid, numprocs, ierr, root

  filename = "/mn/stornext/d5/comap/testdata/data_20h_1046_lvl1.h5"

  inquire(file=trim(filename), size=filesize)

  call wall_time(start_time)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)
  root = 0

  call read_l1_file(trim(filename), data)

  call wall_time(stop_time)

  total_time = stop_time - start_time
  write(*,*) myid, total_time, " s, ", filesize/1024.**2, "Mb"
  write(*,*) myid, filesize/total_time/1024.**2, " Mb/s"

  call mpi_finalize(ierr)

end program time_test
