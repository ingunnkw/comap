program tod2comap
  use comap_lx_mod
  implicit none

!  include "mpif.h"

!  integer(i4b)  :: myid, numprocs, ierr, root

  type tod_type ! ???
     real(dp)     :: samprate, Tsys
     integer(i4b) :: numsamp, numdet, numfreq
     real(dp)     :: fmin, fmax, df
     real(dp), allocatable, dimension(:)     :: t, f             ! (time or freq)
     real(dp), allocatable, dimension(:,:,:) :: d, d_raw, g, rms ! (time, freq, det) 
     real(dp), allocatable, dimension(:,:,:) :: point            ! (time, 3, numdet)
  end type tod_type

  type map_type
     integer(i4b) :: n_x, n_y, numfreq, n_k
     real(dp)     :: dtheta
     real(dp), allocatable, dimension(:)     :: x, y, f, k ! (n_x or n_y or numfreq or n_k)
     real(dp), allocatable, dimension(:,:,:) :: m, rms, dsum, nhit, div ! (n_x, n_y, numfreq)
  end type map_type

  type(tod_type)  :: tod
  type(map_type)  :: map
  type(lx_struct) :: data

  character(len=512) :: filename
  integer(i4b) :: i, j, k

!  call mpi_init(ierr)
!  call mpi_comm_rank(mpi_comm_world, myid, ierr)
!  call mpi_comm_size(mpi_comm_world, numprocs, ierr)
!  root = 0
!  call mpi_finalize(ierr)

  ! Read data
  call read_l3_file(filename, data)


contains

  subroutine compute_maps(data, map)
    implicit none


  end subroutine compute_maps

 
end program tod2comap
