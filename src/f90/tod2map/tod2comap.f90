program tod2comap
  use comap_lx_mod
  use tod2comap_utils
  use tod2comap_mapmaker
  use tod2comap_cl_mod
  implicit none

!  include "mpif.h"

!  integer(i4b)  :: myid, numprocs, ierr, root

  type map_type
     integer(i4b) :: n_x, n_y, numfreq, n_k
     real(dp)     :: dtheta
     real(dp), allocatable, dimension(:)     :: x, y, f, k ! (n_x or n_y or numfreq or n_k)
     real(dp), allocatable, dimension(:,:,:) :: m, rms, dsum, nhit, div ! (n_x, n_y, numfreq)
  end type map_type

!  call mpi_init(ierr)
!  call mpi_comm_rank(mpi_comm_world, myid, ierr)
!  call mpi_comm_size(mpi_comm_world, numprocs, ierr)
!  root = 0

  type(map_type) :: map














!  call mpi_finalize(ierr)

end program tod2comap
