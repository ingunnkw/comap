module comap_map_mod
  use healpix_types
  use quiet_utils
  use quiet_hdf_mod
  implicit none

  real(dp), parameter :: MAP_BASE_PIXSIZE = 1.d0 ! Arcmin

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, nsb, ndet, ndet_tot,  n_k, ntheta ! 2^ntheta
     !real(dp)     :: x0, y0, f0, 
     real(dp)     :: dthetax, dthetay, df
     real(dp)     :: mean_az, mean_el, time(2)
     character(len=512) :: name
     integer(i4b), allocatable, dimension(:)       :: feeds
     real(dp),     allocatable, dimension(:)       :: x, y, k ! (n_x or n_y or n_k)
     real(dp),     allocatable, dimension(:,:)     :: freq    ! (nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:,:) :: m, rms, dsum, div ! (n_x, n_y, nfreq, nsb, ndet)
     integer(i4b), allocatable, dimension(:,:,:,:,:) :: nhit
  end type map_type


contains

  ! Writes an h5 file with maps/rms/nhit
  subroutine output_map_h5(prefix, map, det, sb)
    implicit none
    character(len=*), intent(in)    :: prefix
    type(map_type),   intent(inout) :: map
    integer(i4b), optional :: det, sb

    type(hdf_file)     :: file
    character(len=512) :: filename
    
    filename = trim(prefix)//'_'//trim(map%name)//'.h5'
    call open_hdf_file(trim(filename), file, "w")
    call write_hdf(file, "n_x", map%n_x)
    call write_hdf(file, "n_y", map%n_y)
    call write_hdf(file, "x",   map%x)
    call write_hdf(file, "y",   map%y)
    if (present(det)) then
       call write_hdf(file, "map", map%m(:,:,:,sb:sb,det:det))
       call write_hdf(file, "rms", map%rms(:,:,:,sb:sb,det:det))
       call write_hdf(file, "nhit", map%nhit(:,:,:,sb:sb,det:det))
    else
       call write_hdf(file, "map", map%m)
       call write_hdf(file, "rms", map%rms)
       call write_hdf(file, "nhit", map%nhit)
       call write_hdf(file, "map_beam", sum(map%m,dim=5))
       call write_hdf(file, "rms_beam", sum(map%rms,dim=5))
       call write_hdf(file, "nhit_beam", sum(map%nhit,dim=5))
    end if
    call write_hdf(file, "freq", map%freq)
    call write_hdf(file, "mean_az", map%mean_az)
    call write_hdf(file, "mean_el", map%mean_el)
    call write_hdf(file, "time", map%time)
    call write_hdf(file, "feeds", map%feeds)
    call close_hdf_file(file)

    !call free_map_type(map)

  end subroutine output_map_h5


  subroutine output_submap_h5(prefix, map)
    implicit none
    character(len=*), intent(in) :: prefix
    type(map_type),   intent(in) :: map

    type(hdf_file)     :: file
    character(len=512) :: filename

    filename = trim(prefix)//'_'//trim(map%name)//'.h5'
    call open_hdf_file(trim(filename), file, "w")
    call write_hdf(file, "n_x", map%n_x)
    call write_hdf(file, "n_y", map%n_y)
    call write_hdf(file, "x",   map%x)
    call write_hdf(file, "y",   map%y)
    call write_hdf(file, "map", map%m)!(:,:,1,1,1))
    call write_hdf(file, "rms", map%rms)!(:,:,1,1,1))
    call write_hdf(file, "nhit", map%nhit)!(:,:,1,1,1))
    call write_hdf(file, "freq", map%freq)
    call write_hdf(file, "mean_az", map%mean_az)
    call write_hdf(file, "mean_el", map%mean_el)
    call write_hdf(file, "time", map%time)
    call write_hdf(file, "feeds", map%feeds)
    call close_hdf_file(file)

  end subroutine output_submap_h5


  ! Reads an h5 file
  subroutine read_map_h5(filename,map)
    implicit none
    character(len=*), intent(in)  :: filename
    type(map_type),   intent(out) :: map

    type(hdf_file) :: file
    integer(i4b)   :: nx, ny, nfreq, nsb, ndet, ext(7)

    call free_map_type(map)

    call open_hdf_file(trim(filename), file, "r")

    call get_size_hdf(file, "map", ext)
    nx = ext(1); ny = ext(2); nfreq = ext(3); nsb = ext(4); ndet = ext(5)

    allocate(map%x(nx), map%y(ny))
    allocate(map%m(nx,ny,nfreq,nsb,ndet), map%rms(nx,ny,nfreq,nsb,ndet), map%nhit(nx,ny,nfreq, nsb,ndet))
    allocate(map%freq(nfreq,nsb))

    call read_hdf(file, "n_x", map%n_x)
    call read_hdf(file, "n_y", map%n_y)
    call read_hdf(file, "x",   map%x)
    call read_hdf(file, "y",   map%y)
    call read_hdf(file, "map", map%m)
    call read_hdf(file, "rms", map%rms)
    call read_hdf(file, "nhit", map%nhit)
    call read_hdf(file, "freq", map%freq)
    call read_hdf(file, "time", map%time)
    call read_hdf(file, "mean_az", map%mean_az)
    call read_hdf(file, "mean_el", map%mean_el)
    call read_hdf(file, "feeds", map%feeds)
    call close_hdf_file(file)

  end subroutine read_map_h5


!   ! Creates a .dat file with the maps/rms/nhit for each frequency
!   subroutine output_maps(prefix, map)
!     implicit none
!     character(len=*), intent(in)    :: prefix
!     type(map_type),   intent(inout) :: map

!     integer(i4b)       :: i, j, k, unit
!     character(len=4)   :: itext
!     character(len=512) :: filename

!     unit = getlun()
!     do i = 6, 6
!     !do i = 1, map%nfreq
!        call int2string(i,itext)
!        filename = trim(prefix)//'_freq'//itext//'_map.dat'
!        open(unit, file=trim(filename), recl=100000)
!        write(unit,*) '# n_x = ', map%n_x
!        write(unit,*) '# n_y = ', map%n_y
!        write(unit,*) '# x   = ', real(map%x,sp)
!        write(unit,*) '# y   = ', real(map%y,sp)
!        do j = 1, map%n_x
!           do k = 1, map%n_y
!              write(unit,fmt='(e16.8)',advance='no') map%m(j,k,i,:)
!           end do
!           write(unit,*)
!        end do
!        close(unit)
!     end do

!     unit = getlun()
!     do i = 6, 6
! !    do i = 1, map%nfreq
!        call int2string(i,itext)
!        filename = trim(prefix)//'_freq'//itext//'_rms.dat'
!        open(unit, file=trim(filename), recl=100000)
!        write(unit,*) '# n_x = ', map%n_x
!        write(unit,*) '# n_y = ', map%n_y
!        write(unit,*) '# x   = ', real(map%x,sp)
!        write(unit,*) '# y   = ', real(map%y,sp)
!        do j = 1, map%n_x
!           do k = 1, map%n_y
!              write(unit,fmt='(e16.8)',advance='no') map%rms(j,k,i,:)
!           end do
!           write(unit,*)
!        end do
!        close(unit)
!     end do

!     unit = getlun()
!     do i = 6, 6
! !    do i = 1, map%nfreq
!        call int2string(i,itext)
!        filename = trim(prefix)//'_freq'//itext//'_nhit.dat'
!        open(unit, file=trim(filename), recl=100000)
!        write(unit,*) '# n_x = ', map%n_x
!        write(unit,*) '# n_y = ', map%n_y
!        write(unit,*) '# x   = ', real(map%x,sp)
!        write(unit,*) '# y   = ', real(map%y,sp)
!        do j = 1, map%n_x
!           do k = 1, map%n_y
!              write(unit,fmt='(e16.8)',advance='no') map%nhit(j,k,i,:)
!           end do
!           write(unit,*)
!        end do
!        close(unit)
!     end do

!     call free_map_type(map)

!   end subroutine output_maps

  subroutine nullify_map_type(map)
    implicit none
    type(map_type), intent(inout) :: map

    map%m    = 0.d0
    map%rms  = 0.d0
    map%dsum = 0.d0
    map%div  = 0.d0
    map%nhit = 0.d0

  end subroutine nullify_map_type


  subroutine free_map_type(map)
    implicit none
    type(map_type), intent(inout) :: map

    if (allocated(map%x))     deallocate(map%x) 
    if (allocated(map%y))     deallocate(map%y)
    if (allocated(map%freq))  deallocate(map%freq)
    if (allocated(map%k))     deallocate(map%k)
    if (allocated(map%m))     deallocate(map%m)
    if (allocated(map%rms))   deallocate(map%rms)
    if (allocated(map%dsum))  deallocate(map%dsum)
    if (allocated(map%nhit))  deallocate(map%nhit)
    if (allocated(map%div))   deallocate(map%div)
    if (allocated(map%feeds)) deallocate(map%feeds)

  end subroutine free_map_type

end module comap_map_mod
