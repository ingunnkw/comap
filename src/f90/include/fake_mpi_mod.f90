module fake_mpi_mod
  use healpix_types
  implicit none
  integer(i4b), private   :: nproc = 1, id = 0
  logical(lgt), private   :: initialized = .false.
  integer(i4b), parameter :: MPI_COMM_WORLD = 0

contains

  subroutine mpi_init(ierr)
    implicit none
    integer(i4b) :: ierr
    character(len=512) :: str
    ierr = 0
    ! First try OMPI stuff
    call getenv("OMPI_COMM_WORLD_SIZE", str)
    read(str, *, err=1, end=1) nproc
    call getenv("OMPI_COMM_WORLD_RANK", str)
    read(str, *, err=1, end=1) id
    return
1   call getenv("SLURM_NPROCS", str)
    read(str, *, err=1, end=2) nproc
    call getenv("SLURM_PROCID", str)
    read(str, *, err=1, end=2) id
    return
2   id    = 0
    nproc = 1
  end subroutine

  subroutine mpi_comm_rank(comm, myid, ierr)
    implicit none
    integer(i4b) :: comm, myid, ierr
    if(comm == MPI_COMM_WORLD) then
        myid = id
        ierr = 0
    else
        ierr = 1
    end if
  end subroutine

  subroutine mpi_comm_size(comm, numprocs, ierr)
    implicit none
    integer(i4b) :: comm, numprocs, ierr
    if(comm == MPI_COMM_WORLD) then
        numprocs = nproc
        ierr = 0
    else
        ierr = 1
    end if
  end subroutine

  subroutine mpi_finalize(ierr)
    implicit none
    integer(i4b) :: ierr
    ierr = 0
  end subroutine

end module fake_mpi_mod
  
