module comap_map_mod
  use healpix_types
  use quiet_utils
  use quiet_hdf_mod
  implicit none

  real(dp), parameter :: MAP_BASE_PIXSIZE = 1.d0 ! Arcmin

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, n_k, ntheta ! 2^ntheta
     real(dp)     :: x0, y0, f0, df
     real(dp)     :: dtheta
     real(dp), allocatable, dimension(:)     :: x, y, f, k ! (n_x or n_y or nfreq or n_k)
     real(dp), allocatable, dimension(:,:,:) :: m, rms, dsum, nhit, div ! (n_x, n_y, nfreq)
  end type map_type


contains

  subroutine output_map_h5(filename,map)
    implicit none
    character(len=*), intent(in)    :: filename
    type(map_type),   intent(inout) :: map

    type(hdf_file)     :: file
    
    call open_hdf_file(trim(filename), file, "w")
    call write_hdf(file, "n_x", map%n_x)
    call write_hdf(file, "n_y", map%n_y)
    call write_hdf(file, "x",   map%x)
    call write_hdf(file, "y",   map%y)
    call write_hdf(file, "map", map%m)
    call write_hdf(file, "rms", map%rms)
    call write_hdf(file, "nhit", map%nhit)
    call close_hdf_file(file)

  end subroutine output_map_h5


  subroutine read_map_h5(filename,map)
    implicit none
    character(len=*), intent(in)  :: filename
    type(map_type),   intent(out) :: map

    type(hdf_file) :: file
    integer(i4b)   :: nx, ny, nfreq, ext(7)

    call open_hdf_file(trim(filename), file, "r")

    call get_size_hdf(file, "map", ext)
    nx = ext(1); ny = ext(2); nfreq = ext(3)

    allocate(map%x(nx), map%y(ny))
    allocate(map%m(nx,ny,nfreq), map%rms(nx,ny,nfreq), map%nhit(nx,ny,nfreq))

    call read_hdf(file, "n_x", map%n_x)
    call read_hdf(file, "n_y", map%n_y)
    call read_hdf(file, "x",   map%x)
    call read_hdf(file, "y",   map%y)
    call read_hdf(file, "map", map%m)
    call read_hdf(file, "rms", map%rms)
    call read_hdf(file, "nhit", map%nhit)
    call read_hdf_file(file)

  end subroutine read_map_h5

end module comap_map_mod
