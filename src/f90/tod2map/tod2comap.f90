program tod2comap
  use comap_lx_mod
  use comap_scan_mod
  use comap_acceptlist_mod
  use quiet_fft_mod
  implicit none

!  include "mpif.h"

  real(dp), parameter :: MAP_BASE_PIXSIZE = 1.d0 ! Arcmin

  type tod_type
     real(dp)     :: samprate, Tsys
     integer(i4b) :: nsamp, ndet, nfreq
     real(dp)     :: fmin, fmax, df
     real(dp), allocatable, dimension(:)     :: t, f             ! (time or freq)
     real(dp), allocatable, dimension(:,:,:) :: d, d_raw, g, rms ! (time, freq, det) 
     real(dp), allocatable, dimension(:,:)   :: point            ! (3, time)
  end type tod_type

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, n_k, ntheta ! 2^ntheta
     real(dp)     :: x0, y0, f0, df
     real(dp)     :: dtheta
     real(dp), allocatable, dimension(:)     :: x, y, f, k ! (n_x or n_y or nfreq or n_k)
     real(dp), allocatable, dimension(:,:,:) :: m, rms, dsum, nhit, div ! (n_x, n_y, nfreq)
  end type map_type

  type(tod_type)        :: tod
  type(map_type)        :: map
  type(lx_struct)       :: data
  type(comap_scan_info) :: scan
  type(acceptlist)      :: alist

  character(len=512)    :: filename, parfile, acceptfile, prefix
  integer(i4b)          :: nscan, i, j, k
  integer(i4b)          :: myid, numprocs, ierr, root

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)
  root = 0

  parfile = 'param_test.txt'

  call initialize_scan_mod(parfile)
  nscan = get_num_scans()

  do i = 1, nscan 
     if (myid == 0) write(*,*) i, 'of', nscan
     if (allocated(alist%status)) deallocate(alist%status)
     call get_scan_info(i,scan)
     filename = scan%l3file
     prefix = 'files/'//trim(scan%object)//'_'//trim(itoa(scan%cid)) ! patchID_scanID
     ! Get TOD / read level 3 file
     call get_tod(trim(filename), data, tod)
     allocate(alist%status(tod%nfreq,tod%ndet,nscan))
     alist%status = 0 ! TODO: check if any parts of the scan has been rejected
     ! Write TOD to file
     !if (myid == 0) call output_tod(trim(prefix), 1, tod, alist, i)
     ! Compute maps
     call compute_maps(data, tod, map, alist, i)
     ! Write maps to file (includes rms and hit count)
     if (myid == 0) call output_maps(trim(prefix), map)
  end do

  if (myid == 0) write(*,*) 'Done'
  call mpi_finalize(ierr)

contains

  subroutine hp_filter(nu_cut, tod_d)
    implicit none

    integer(i4b),                intent(in)    :: nu_cut
    real(dp),     dimension(1:), intent(inout) :: tod_d

    integer(i4b) :: n
    complex(dpc), allocatable, dimension(:) :: tod_fft

    n = size(tod_d)/2 + 1
    allocate(tod_fft(n))

    call fft(tod_d, tod_fft, 1)
    tod_fft(1:nu_cut) = 0.d0
    call fft(tod_d, tod_fft, -1)
    
    deallocate(tod_fft)

  end subroutine hp_filter


  subroutine get_tod(l3file, data, tod)
    implicit none
    character(len=*), intent(in)    :: l3file
    type(lx_struct),  intent(inout) :: data
    type(tod_type),   intent(inout) :: tod

    integer(i4b) :: i, j, k, nu_cut
   
    nu_cut = 10000

    ! Read data
    call read_l3_file(l3file, data)
    call free_tod_type(tod)

    tod%samprate = data%samprate
    tod%nsamp = size(data%time)
    tod%nfreq = size(data%nu)
    tod%ndet  = size(data%tod,3)

    allocate( tod%t(tod%nsamp), tod%f(tod%nfreq), &
         & tod%point(3,tod%nsamp), &
         & tod%d_raw(tod%nsamp, tod%nfreq, tod%ndet), &
         & tod%d(tod%nsamp, tod%nfreq, tod%ndet), &
         & tod%g(0,0,0), & !(tod%nsamp, tod%nfreq, tod%ndet), &               ! ??
         & tod%rms(tod%nsamp, tod%nfreq, tod%ndet))

    tod%t = data%time; tod%f = data%nu
    tod%point = data%point ! call make_angles_safe(tod%point(1,:),maxang)
    tod%g     = data%gain

    !write(*,*) shape(data%point), tod%nsamp

    do k = 1, tod%ndet
       do j = 1, tod%nfreq
          do i = 1, tod%nsamp
            ! tod%d_raw(i,j,k) = mean
             tod%d_raw(i,j,k) = data%tod(i,j,k) ! ??
          end do
          
          ! Apply high pass filter
          tod%d(:,j,k) = tod%d_raw(:,j,k)
          call hp_filter(nu_cut, tod%d(:,j,k))

          ! Estimate RMS
          tod%rms(:,j,k) = sqrt(variance(tod%d(:,j,k)))
       end do
    end do

  end subroutine get_tod


  subroutine output_tod(prefix, det, tod, alist, scan_nr)
    implicit none
    character(len=*), intent(in) :: prefix
    type(tod_type),   intent(in) :: tod
    type(acceptlist), intent(in) :: alist
    integer(i4b),     intent(in) :: det, scan_nr

    character(len=512) :: filename
    character(len=4)   :: jtext
    integer(i4b) :: unit, i, j

    unit = getlun()
    do j = 1, tod%nfreq
       if (alist%status(j,det,scan_nr) == 0) then
          call int2string(j,jtext)
          filename = trim(prefix) // '_freq' // jtext // '_tod.dat'
          open(unit, file=trim(filename), recl=4096)
          do i = 1, tod%nsamp
             write(unit, fmt='(f16.8)', advance='no') tod%t(i)
             write(unit, fmt='(2f24.8)', advance='no') tod%d(i,j,det), tod%d_raw(i,j,det)
             write(unit,*)
          end do
          close(unit)
       end if
    end do

  end subroutine output_tod
  

  subroutine compute_maps(data, tod, map, alist, scan_nr)
    implicit none
    type(lx_struct),  intent(in)    :: data
    type(tod_type),   intent(in)    :: tod
    type(map_type),   intent(inout) :: map
    type(acceptlist), intent(in)    :: alist
    integer(i4b),     intent(in)    :: scan_nr
    
    integer(i4b) :: i, j, k, p, q, fs
    real(dp)     :: x_min, x_max, y_min, y_max, pad, gain_hc

    call free_map_type(map)

    ! Set up map grid
    fs = 4
    pad = 0.d0!0.3d0 ! degrees
    map%nfreq = tod%nfreq
    map%dtheta = 5.d0/60.d0 ! Arcmin
    !write(*,*) data%point_lim

    x_min = minval(tod%point(1,fs:)) - pad; x_max =  maxval(tod%point(1,fs:)) + pad
    y_min = minval(tod%point(2,fs:)) - pad; y_max =  maxval(tod%point(2,fs:)) + pad
    !x_min = data%point_lim(1) - pad; x_max = data%point_lim(2) + pad
    !y_min = data%point_lim(3) - pad; y_max = data%point_lim(4) + pad
    map%n_x = (x_max-x_min)/map%dtheta+1; map%n_y = (y_max-y_min)/map%dtheta+1
    allocate(map%x(map%n_x), map%y(map%n_y))
    do i = 1, map%n_x
       map%x(i) = x_min + (i-1)*map%dtheta
    end do
    do i = 1, map%n_y
       map%y(i) = y_min + (i-1)*map%dtheta
    end do

    
    ! Set up map structures
    if (.not. allocated(map%dsum)) then
       allocate(map%m(map%n_x, map%n_y, map%nfreq), map%dsum(map%n_x, map%n_y, map%nfreq), &
            & map%nhit(map%n_x, map%n_y, map%nfreq), map%div(map%n_x, map%n_y, map%nfreq), &
            & map%rms(map%n_x, map%n_y, map%nfreq))
       map%dsum = 0.d0
       map%nhit = 0.d0
       map%div  = 0.d0
    end if

    gain_hc = 1.0 ! Replaces tod%g(i,j,k) in the code below
    
    ! Co-add into maps
    do k = 1, tod%ndet
       do i = fs, tod%nsamp
          p = min(max(int((tod%point(1,i)-x_min)/map%dtheta),1),map%n_x)
          q = min(max(int((tod%point(2,i)-y_min)/map%dtheta),1),map%n_y)
          do j = 1, tod%nfreq
             if (alist%status(j,k,scan_nr) == 0) then
                map%dsum(p,q,j) = map%dsum(p,q,j) + gain_hc    / tod%rms(i,j,k)**2 * tod%d(i,j,k)
                map%div(p,q,j)  = map%div(p,q,j)  + gain_hc**2 / tod%rms(i,j,k)**2
                map%nhit(p,q,j) = map%nhit(p,q,j) + 1.d0
             end if
          end do
       end do
    end do
    where(map%nhit > 0)
       map%m   = map%dsum / map%div
       map%rms = 1.d0 / sqrt(map%div)
    elsewhere
       map%m   = 0.d0
       map%rms = 0.d0
    end where
    
!!$    ! Report reduced chisquares
!!$    do i = 1, map%nfreq
!!$       nu    = count(map%nhit(:,:,i)>0)
!!$       chisq = sum((map%m(:,:,i)/map%rms(:,:,i))**2,map%nhit(:,:,i)>0)
!!$       write(*,fmt='(a,i5,a,f8.3,a,f8.3)') 'Freq = ', i, ', red chisq = ', chisq/nu, ', sigma = ', (chisq-nu)/sqrt(2.d0*nu)
!!$    end do

  end subroutine compute_maps


  subroutine output_maps(prefix, map)
    implicit none
    character(len=*), intent(in) :: prefix
    type(map_type),   intent(in) :: map

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

  end subroutine output_maps


  subroutine free_tod_type(tod)
    implicit none
    type(tod_type), intent(inout) :: tod 

    if (allocated(tod%t))     deallocate(tod%t)
    if (allocated(tod%f))     deallocate(tod%f)
    if (allocated(tod%d))     deallocate(tod%d)
    if (allocated(tod%d_raw)) deallocate(tod%d_raw)
    if (allocated(tod%g))     deallocate(tod%g)
    if (allocated(tod%rms))   deallocate(tod%rms)
    if (allocated(tod%point)) deallocate(tod%point)

  end subroutine free_tod_type

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
   
end program tod2comap
