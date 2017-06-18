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

  ! Writes an h5 file with maps/rms/nhit
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

    call free_map_type(map)

  end subroutine output_map_h5


  ! Reads an h5 file
  subroutine read_map_h5(filename,map)
    implicit none
    character(len=*), intent(in)  :: filename
    type(map_type),   intent(out) :: map

    type(hdf_file) :: file
    integer(i4b)   :: nx, ny, nfreq, ext(7)

    call free_map_type(map)

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
    call close_hdf_file(file)

  end subroutine read_map_h5


  ! Creates a .dat file with the maps/rms/nhit for each frequency
  subroutine output_maps(prefix, map)
    implicit none
    character(len=*), intent(in)    :: prefix
    type(map_type),   intent(inout) :: map

    integer(i4b)       :: i, j, k, unit
    character(len=4)   :: itext
    character(len=512) :: filename

    unit = getlun()
    do i = 1, map%nfreq
       call int2string(i,itext)
       filename = trim(prefix)//'_freq'//itext//'_map.dat'
       open(unit, file=trim(filename), recl=100000)
       write(unit,*) '# n_x = ', map%n_x
       write(unit,*) '# n_y = ', map%n_y
       write(unit,*) '# x   = ', real(map%x,sp)
       write(unit,*) '# y   = ', real(map%y,sp)
       do j = 1, map%n_x
          do k = 1, map%n_y
             write(unit,fmt='(e16.8)',advance='no') map%m(j,k,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

    unit = getlun()
    do i = 1, map%nfreq
       call int2string(i,itext)
       filename = trim(prefix)//'_freq'//itext//'_rms.dat'
       open(unit, file=trim(filename), recl=100000)
       write(unit,*) '# n_x = ', map%n_x
       write(unit,*) '# n_y = ', map%n_y
       write(unit,*) '# x   = ', real(map%x,sp)
       write(unit,*) '# y   = ', real(map%y,sp)
       do j = 1, map%n_x
          do k = 1, map%n_y
             write(unit,fmt='(e16.8)',advance='no') map%rms(j,k,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

    unit = getlun()
    do i = 1, map%nfreq
       call int2string(i,itext)
       filename = trim(prefix)//'_freq'//itext//'_nhit.dat'
       open(unit, file=trim(filename), recl=100000)
       write(unit,*) '# n_x = ', map%n_x
       write(unit,*) '# n_y = ', map%n_y
       write(unit,*) '# x   = ', real(map%x,sp)
       write(unit,*) '# y   = ', real(map%y,sp)
       do j = 1, map%n_x
          do k = 1, map%n_y
             write(unit,fmt='(e16.8)',advance='no') map%nhit(j,k,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

    call free_map_type(map)

  end subroutine output_maps


  subroutine free_map_type(map)
    implicit none
    type(map_type), intent(inout) :: map

    if (allocated(map%x))    deallocate(map%x) 
    if (allocated(map%y))    deallocate(map%y)
    if (allocated(map%f))    deallocate(map%f)
    if (allocated(map%k))    deallocate(map%k)
    if (allocated(map%m))    deallocate(map%m)
    if (allocated(map%rms))  deallocate(map%rms)
    if (allocated(map%dsum)) deallocate(map%dsum)
    if (allocated(map%nhit)) deallocate(map%nhit)
    if (allocated(map%div))  deallocate(map%div)

  end subroutine free_map_type

end module comap_map_mod
