module l1_read_comap_mod
  use healpix_types
  use spline_1D_mod
  use quiet_utils
  use quiet_l1_defs
  use quiet_apex_mod
  use quiet_hdf_mod
!  use comap_Lx_mod
  implicit none

  !integer(i4b), private :: num_det = 19
  type hk_type
     integer(i4b)                                    :: n, n_t
     character(len=128), allocatable, dimension(:)   :: name
     real(dp),           allocatable, dimension(:)   :: time
     real(dp),           allocatable, dimension(:,:) :: value
  end type hk_type

  type point_type
     real(dp)                                  :: start_mjd
     integer(i2b), allocatable, dimension(:)   :: mode
     real(dp),     allocatable, dimension(:)   :: time
     real(dp),     allocatable, dimension(:,:) :: orig_point
  end type point_type
contains
  ! ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ! ;; MAIN l1_read PROCEDURE;;
  ! ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine l1_read_point(filename, housekeeping, pointing, detectors, status_code)
    implicit none

    character(len=*),                   intent(in)              :: filename
    integer(i4b)                                                :: ext(7), npoint, nsamp
    type(hdf_file)                                              :: file
    integer(i4b),        dimension(1:), intent(inout), optional :: detectors
    !type(Lx_struct),                    intent(inout), optional :: data
    type(hk_type),                    intent(inout), optional :: housekeeping
    type(point_type),                 intent(inout), optional :: pointing
    integer(i4b),                       intent(out),   optional :: status_code
    !call free_lx_struct(data)
    
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "point", ext)
    npoint = ext(2); nsamp = ext(1)
    call allocate_point_type(nsamp, npoint, pointing)
    
    call read_hdf(file, "time",  pointing%time)
    
    call read_hdf(file, "point", pointing%orig_point)
                                                                                                          
    call close_hdf_file(file)

  end subroutine L1_read_point
  
  subroutine allocate_point_type(num_samples, num_pointings, pointing)
    implicit none

    integer(i4b),       intent(in)  :: num_samples, num_pointings
    type(point_type), intent(out) :: pointing

    allocate(pointing%mode(num_samples))
    allocate(pointing%time(num_samples))
    allocate(pointing%orig_point(num_pointings,num_samples))
  end subroutine allocate_point_type

  subroutine deallocate_point_type(pointing)
    implicit none

    type(point_type), intent(inout) :: pointing

    pointing%start_mjd = 0.d0
    if (allocated(pointing%mode))              deallocate(pointing%mode)
    if (allocated(pointing%time))              deallocate(pointing%time)
    if (allocated(pointing%orig_point))        deallocate(pointing%orig_point)
    
  end subroutine deallocate_point_type


end module l1_read_comap_mod
