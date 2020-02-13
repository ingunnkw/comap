module comap_map_mod
  use healpix_types
  use quiet_utils
  use quiet_hdf_mod
  implicit none

  real(dp), parameter :: MAP_BASE_PIXSIZE = 1.d0 ! Arcmin

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, nsb, ndet, ndet_tot, n_k, ntheta, nside ! 2^ntheta
     !real(dp)     :: x0, y0, f0, 
     real(dp)     :: dthetax, dthetay, df
     real(dp)     :: mean_az, mean_el, time(2), center(2)
     character(len=512) :: name

     integer(i4b), allocatable, dimension(:)          :: feeds
     real(dp),     allocatable, dimension(:)          :: x, y, k                           ! (n_x or n_y or n_k)
     real(dp),     allocatable, dimension(:,:)        :: freq                              ! (nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:,:)  :: m, rms, dsum, div                 ! (n_x, n_y, nfreq, nsb, ndet)
     real(sp),     allocatable, dimension(:,:,:,:)    :: m_co, rms_co, dsum_co, div_co     ! (n_x, n_y, nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:)    :: m_sim, rms_sim, dsum_sim, div_sim ! (n_x, n_y, nfreq, nsb, ndet)
     integer(i4b), allocatable, dimension(:,:,:,:,:)  :: nhit                              ! (n_x, n_y, nfreq, nsb, ndet)
     integer(i4b), allocatable, dimension(:,:,:,:)    :: nhit_co                           ! (n_x, n_y, nfreq, nsb)

  end type map_type


contains

  ! Writes an h5 file with maps/rms/nhit
  subroutine output_map_h5(filename, map, det, sb)
    implicit none
    character(len=*), intent(in)    :: filename
    type(map_type),   intent(inout) :: map
    integer(i4b), optional :: det, sb

    type(hdf_file)     :: file
    
    call open_hdf_file(trim(filename), file, "w")
    call write_hdf(file, "n_x", map%n_x)
    call write_hdf(file, "n_y", map%n_y)
    call write_hdf(file, "x",   map%x)
    call write_hdf(file, "y",   map%y)
    call write_hdf(file, "patch_center", map%center)
    call write_hdf(file, "nside", map%nside)
    if (present(det)) then
       call write_hdf(file, "map", map%m(:,:,:,sb:sb,det:det))
       call write_hdf(file, "rms", map%rms(:,:,:,sb:sb,det:det))
       call write_hdf(file, "nhit", map%nhit(:,:,:,sb:sb,det:det))
    else
       call write_hdf(file, "map", map%m)
       call write_hdf(file, "rms", map%rms)
       call write_hdf(file, "nhit", map%nhit)
       call write_hdf(file, "map_beam", map%m_co)
       call write_hdf(file, "rms_beam", map%rms_co)
       call write_hdf(file, "nhit_beam", map%nhit_co)

    end if
    call write_hdf(file, "freq", map%freq)
    call write_hdf(file, "mean_az", map%mean_az)
    call write_hdf(file, "mean_el", map%mean_el)
    call write_hdf(file, "time", map%time)
    call write_hdf(file, "feeds", map%feeds)
    call close_hdf_file(file)

    !call free_map_type(map)

  end subroutine output_map_h5


  subroutine output_submap_sim_h5(filename, map, sim)
    implicit none
    character(len=*), intent(in) :: filename
    type(map_type),   intent(in) :: map
    integer(i4b),     intent(in) :: sim

    type(hdf_file)     :: file
 

    call open_hdf_file(trim(filename), file, "w")

    ! For simulated data 
    call write_hdf(file, "map_sim", map%m_sim) 
    call write_hdf(file, "rms_sim", map%rms_sim)
    call write_hdf(file, "sim", sim)  
    ! call write_hdf(file, "map_sim_beam", sum(map%m_sim, dim=5))
    ! call write_hdf(file, "rms_sim_beam", sum(map%rms_sim, dim=5))

    call close_hdf_file(file)

  end subroutine output_submap_sim_h5


  subroutine output_submap_h5(filename, map)
    implicit none
    character(len=*), intent(in) :: filename
    type(map_type),   intent(in) :: map

    type(hdf_file)     :: file
    
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
    call write_hdf(file, "patch_center", map%center)
    call write_hdf(file, "nside", map%nside)
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
    allocate(map%freq(nfreq,nsb), map%feeds(ndet))
    

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
    call read_hdf(file, "patch_center", map%center)
    call read_hdf(file, "nside", map%nside)

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

    map%m       = 0.d0
    map%rms     = 0.d0
    map%dsum    = 0.d0
    map%div     = 0.d0
    map%nhit    = 0.d0
    map%m_co    = 0.d0
    map%rms_co  = 0.d0
    map%dsum_co = 0.d0
    map%div_co  = 0.d0
    map%nhit_co = 0.d0

    ! Simulated data
    map%m_sim    = 0.d0
    map%rms_sim  = 0.d0
    map%dsum_sim = 0.d0
    map%div_sim  = 0.d0 

  end subroutine nullify_map_type


  subroutine free_map_type(map)
    implicit none
    type(map_type), intent(inout) :: map

    if (allocated(map%x))       deallocate(map%x) 
    if (allocated(map%y))       deallocate(map%y)
    if (allocated(map%freq))    deallocate(map%freq)
    if (allocated(map%k))       deallocate(map%k)
    if (allocated(map%m))       deallocate(map%m)
    if (allocated(map%rms))     deallocate(map%rms)
    if (allocated(map%dsum))    deallocate(map%dsum)
    if (allocated(map%nhit))    deallocate(map%nhit)
    if (allocated(map%div))     deallocate(map%div)
    if (allocated(map%feeds))   deallocate(map%feeds)
    if (allocated(map%m_co))    deallocate(map%m_co)
    if (allocated(map%rms_co))  deallocate(map%rms_co)
    if (allocated(map%nhit_co)) deallocate(map%nhit_co)
    if (allocated(map%div_co))  deallocate(map%div_co)
    if (allocated(map%dsum_co)) deallocate(map%dsum_co)

    ! Simulated data 
    if (allocated(map%m_sim))    deallocate(map%m_sim)
    if (allocated(map%rms_sim))  deallocate(map%rms_sim) 
    if (allocated(map%dsum_sim)) deallocate(map%dsum_sim) 
    if (allocated(map%div_sim))  deallocate(map%div_sim) 
 
  end subroutine free_map_type

end module comap_map_mod
