module quiet_mpi_mod
  use healpix_types
  use rngmod
  implicit none

  include "mpif.h"
contains

  subroutine initialize_random_seeds(comm, base_seed, rng_handle)
    implicit none

    integer(i4b),     intent(in)  :: comm, base_seed
    type(planck_rng), intent(out) :: rng_handle

    integer(i4b) :: i, myid, root, ierr, numprocs, seed
    integer(i4b), dimension(MPI_STATUS_SIZE) :: status

    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, numprocs, ierr)
    root = 0

    ! Initialize random number generator
    if (myid == root) then
       
       seed = base_seed
       call rand_init(rng_handle, seed)
       do i = 1, numprocs-1
          seed = nint(rand_uni(rng_handle)*1000000.d0)
          call mpi_send(seed, 1, MPI_INTEGER, i, 98, comm, ierr)
       end do
     
    else 
       
       call mpi_recv(seed, 1, MPI_INTEGER, root, 98, comm, status, ierr)
       call rand_init(rng_handle, seed)
       
    end if

  end subroutine initialize_random_seeds

  subroutine mpi_write_shortcut(unit, string, newline)
    implicit none
    character(len=*) string
    character(len=len_trim(string)+1) buffer
    integer(i4b) :: k, unit
    logical(lgt), optional :: newline
    logical(lgt) :: nl
    nl = .true.; if(present(newline)) nl = newline
    if(nl) then
       buffer = string
       buffer(len(buffer):len(buffer)) = achar(10)
       call mpi_file_write_shared(unit, buffer, len(buffer), MPI_CHARACTER, MPI_STATUS_IGNORE, k)
       !call mpi_file_sync(unit, k)
     else
       call mpi_file_write_shared(unit, string, len_trim(string), MPI_CHARACTER, MPI_STATUS_IGNORE, k)
     end if
  end subroutine

end module quiet_mpi_mod
  
