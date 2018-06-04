program tod2comap
  use comap_lx_mod
  use comap_map_mod
  use comap_scan_mod
  use comap_acceptlist_mod
  use quiet_fft_mod
  use quiet_hdf_mod
  use tod2comap_mapmaker
  !use tod2comap_utils
  implicit none

!  include "mpif.h"

  !type tod_type
  !   real(dp)     :: samprate, Tsys
  !   integer(i4b) :: nsamp, ndet, nfreq, nsb
  !   real(dp)     :: fmin, fmax, df
  !   real(dp), allocatable, dimension(:)       :: t                ! (time or freq)
  !   real(dp), allocatable, dimension(:,:,:,:) :: d, d_raw, g, rms ! (time, freq,  sb, det) 
  !   real(dp), allocatable, dimension(:,:)     :: point, f         ! (3, time) or (sb, freq)
  !end type tod_type


  type(tod_type), allocatable, dimension(:) :: tod
  type(map_type)        :: map
  type(comap_scan_info) :: scan
  type(acceptlist)      :: alist

  integer(i4b), allocatable, dimension(:,:) :: pixels
  character(len=512)    :: filename, parfile, acceptfile, prefix, pre
  integer(i4b)          :: nscan, i, j, k, det, sb, freq
  integer(i4b)          :: myid, numprocs, ierr, root

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)
  root = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! mpi per scan, all detectors
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1, parfile)

  call initialize_scan_mod(parfile)
  nscan = get_num_scans()
  allocate(tod(nscan))
  !nscan = 1
!  call free_map_type(map)

  ! This loop currently requiers that all scans are of the same patch

  do i = 1, nscan 
     !if (myid == 0) write(*,*) i, 'of', nscan
     if (allocated(alist%status)) deallocate(alist%status)
     call get_scan_info(i,scan)

     ! Get TOD / read level 3 file
     filename = scan%l3file
     call get_tod(trim(filename), tod(i))


     allocate(alist%status(tod(i)%nfreq, tod(i)%ndet))
     alist%status = 0 ! TODO: check if any parts of the scan has been rejected

     !pre = 'files/'
     !call output_tod(trim(pre), 1, tod(i), alist)

     if (i == 1) call initialize_mapmaker(map, tod(1))

     call time2pix(tod(i), map)

  end do



  ! Compute and co-add maps
  !call binning(map, tod(i), alist)
  if (myid == 0) write(*,*) "CG mapmaker ..."
  det = 1
  do sb = 1, tod(1)%nsb
     if (myid == 0) write(*,*) "sb", sb
     !do freq = tod(1)%nfreq/4, tod(1)%nfreq/4
     do freq = 1, tod(1)%nfreq
        call pcg_mapmaker(tod, map, alist, det, sb, freq)
     end do
  end do
  

  !end do
  if (myid == 0) write(*,*) "Writing to file ..."
  !write(*,*) trim(itoa(scan%sid))
  prefix = '/mn/stornext/u3/mariekf/work/comap/src/f90/tod2map/files/'//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID

  !if (binning) call finalize_binning(map)

  if (myid == 0) call output_map_h5(trim(prefix), map)! call output_maps(trim(prefix), map)

  !end do

  if (myid == 0) write(*,*) 'Done'
  do i = 1, nscan
     call free_tod_type(tod(i))
  end do
  deallocate(tod)
  call free_map_type(map)
  call mpi_finalize(ierr)

contains

  ! subroutine hp_filter(nu_cut, tod_d, samp_rate)
  !   implicit none

  !   real(dp),                    intent(in)    :: nu_cut, samp_rate
  !   real(dp),     dimension(1:), intent(inout) :: tod_d

  !   integer(i4b) :: n, ind_cut
  !   complex(dpc), allocatable, dimension(:) :: tod_fft

  !   n = size(tod_d)/2 + 1
  !   allocate(tod_fft(n))
  !   ind_cut = freq2ind(nu_cut, samp_rate, n) ! n or nsamps??????
  !   call fft(tod_d, tod_fft, 1)
  !   tod_fft(1:ind_cut) = 0.d0
  !   call fft(tod_d, tod_fft, -1)
    
  !   deallocate(tod_fft)

  ! end subroutine hp_filter


  ! subroutine get_tod(l3file, data, tod)
  !   implicit none
  !   character(len=*), intent(in)    :: l3file
  !   type(lx_struct),  intent(inout) :: data
  !   type(tod_type),   intent(inout) :: tod

  !   integer(i4b) :: i, j, k, l
  !   real(dp)     :: nu_cut
   
  !   nu_cut = 0.1d0

  !   ! Read data
  !   call read_l3_file(l3file, data)
  !   call free_tod_type(tod)

  !   tod%samprate = data%samprate
  !   tod%nsamp = size(data%time)
  !   tod%nfreq = size(data%nu,1)
  !   tod%nsb   = size(data%tod,3)
  !   tod%ndet  = size(data%tod,4)

  !   !write(*,*) tod%nsamp, tod%nfreq, tod%nsb, tod%ndet

  !   allocate( tod%t(tod%nsamp), tod%f(tod%nsb, tod%nfreq), &
  !        & tod%point(3,tod%nsamp), &
  !        & tod%d_raw(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
  !        & tod%d(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
  !        & tod%g(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), & 
  !        & tod%rms(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet))

  !   tod%t = data%time; tod%f = data%nu
  !   tod%point = data%point!_cel ! call make_angles_safe(tod%point(1,:),maxang)
  !   tod%g     = data%gain

  !   !write(*,*) shape(data%point), tod%nsamp

  !   do k = 1, tod%ndet
  !      do l = 1, tod%nsb
  !         do j = 1, tod%nfreq
  !            !do j = 6, 6
  !            do i = 1, tod%nsamp
  !               tod%d_raw(i,j,l,k) = data%tod(i,j,l,k)
  !            end do
          
  !            ! Apply high pass filter
  !            tod%d(:,j,l,k) = tod%d_raw(:,j,l,k)
  !            tod%d(:,j,l,k) = tod%d(:,j,l,k) - mean(tod%d(:,j,l,k))
  !            call hp_filter(nu_cut, tod%d(:,j,l,k),tod%samprate)

  !            ! Estimate RMS
  !            tod%rms(:,j,l,k) = sqrt(variance(tod%d(:,j,l,k)))
  !         end do
  !      end do
  !   end do

  ! end subroutine get_tod


  subroutine output_tod(prefix, det, tod, alist)
    implicit none
    character(len=*), intent(in) :: prefix
    type(tod_type),   intent(in) :: tod
    type(acceptlist), intent(in) :: alist
    integer(i4b),     intent(in) :: det

    character(len=512) :: filename
    character(len=4)   :: jtext
    integer(i4b) :: unit, i, j, l

    unit = getlun()
    j=1!do j = 1, tod%nfreq
       if (alist%status(j,det) == 0) then
          call int2string(j,jtext)
          filename = trim(prefix) // 'tod.dat' !trim(prefix) // '_freq' // jtext // '_tod.dat'
          open(unit, file=trim(filename), recl=4096)
          do i = 1, tod%nsamp-100000
             write(unit, fmt='(f16.8)', advance='no') tod%t(i)
             write(unit, fmt='(2f24.8)', advance='no') tod%d(i,j,:,det)!, tod%d_raw(i,j,det)
             write(unit,*)
          end do
          close(unit)
       end if
       !end do
       
  end subroutine output_tod

  
  ! subroutine compute_scan_maps(data, tod, map, alist, det)
  !   implicit none
  !   type(lx_struct),  intent(in)    :: data
  !   type(tod_type),   intent(in)    :: tod
  !   type(map_type),   intent(inout) :: map
  !   type(acceptlist), intent(in)    :: alist
  !   integer(i4b),     intent(in)    :: det
    
  !   integer(i4b) :: i, j, k, l, p, q, fs, st
  !   real(dp)     :: x_min, x_max, y_min, y_max, pad, gain_hc
  !   real(8), parameter :: PI = 4*atan(1.d0)

  !   if (myid == 0) write(*,*) 'Beginning coadding'

    ! ! Set up map grid
    ! fs = 200 !85291
    ! st = tod%nsamp-200 !93579!tod%namp! - 100000
    ! pad = 0.3d0 ! degrees
    ! map%nfreq = tod%nfreq
    ! map%nsb = tod%nsb
    ! map%dthetay = 1.5d0/60.d0 ! resolution
    ! map%dthetax = map%dthetay/abs(cos(mean(data%point_tel(2,:))*PI/180.d0))

    ! x_min = minval(tod%point(1,fs:st)) - pad; x_max =  maxval(tod%point(1,fs:st)) + pad
    ! y_min = minval(tod%point(2,fs:st)) - pad; y_max =  maxval(tod%point(2,fs:st)) + pad
    ! !write(*,*) 'RA', x_min, x_max
    ! !x_min = data%point_lim(1); x_max = data%point_lim(2)
    ! !y_min = data%point_lim(3); y_max = data%point_lim(4)
    ! !write(*,*) x_min, x_max, y_min, y_max
    ! !write(*,*) data%point_lim
    ! if (.not. allocated(map%x)) then 
    !    map%n_x = (x_max-x_min)/map%dthetax+1; map%n_y = (y_max-y_min)/map%dthetay+1
    !    if (myid == 0) write(*,*) map%n_x, map%n_y
    !    allocate(map%x(map%n_x), map%y(map%n_y))
    !    do i = 1, map%n_x
    !       map%x(i) = x_min + (i-1)*map%dthetax
    !    end do
    !    do i = 1, map%n_y
    !       map%y(i) = y_min + (i-1)*map%dthetay
    !    end do
    ! end if
    ! !write(*,*) map%x(8), map%y(9)

    ! ! Set up map structures
    ! if (.not. allocated(map%dsum)) then
    !    allocate(map%m(map%n_x, map%n_y, map%nfreq, map%nsb), &
    !         & map%dsum(map%n_x, map%n_y, map%nfreq, map%nsb), &
    !         & map%nhit(map%n_x, map%n_y, map%nfreq, map%nsb), &
    !         & map%div(map%n_x, map%n_y, map%nfreq, map%nsb), &
    !         & map%rms(map%n_x, map%n_y, map%nfreq, map%nsb))
    !    map%dsum = 0.d0
    !    map%nhit = 0.d0
    !    map%div  = 0.d0
    ! end if

    ! gain_hc = 1.0 ! Replaces tod%g(i,j,k) in the code below

    
    ! Co-add into maps

    !do k = 1, 1
    !   l = 1
    !   open(58,file='tod.dat', recl=1024)
    !   do i = fs, st
    !      write(58,*) i, tod%d(i,1,l,k), data%point_tel(1,i), data%point_tel(2,i)
    !   end do
    !   close(58)
    !end do

    ! do i = fs, st!tod%nsamp
    !    p = min(max(int((tod%point(1,i)-x_min)/map%dthetax),1),map%n_x)
    !    q = min(max(int((tod%point(2,i)-y_min)/map%dthetay),1),map%n_y)
    !    do l = 1, tod%nsb
    !       !do j = 6, 6
    !       do j = 1, tod%nfreq
    !          if (alist%status(j,det) == 0) then
    !             if (tod%g(1,j,l,det) .ne. 0.d0) then

    !                map%dsum(p,q,j,l) = map%dsum(p,q,j,l) + tod%g(1,j,l,det)    / tod%rms(i,j,l,det)**2 * tod%d(i,j,l,det)
    !                map%div(p,q,j,l)  = map%div(p,q,j,l)  + tod%g(1,j,l,det)**2 / tod%rms(i,j,l,det)**2
    !                map%nhit(p,q,j,l) = map%nhit(p,q,j,l) + 1.d0

    !             end if
    !          end if
    !       end do
    !    end do
    ! end do

    !if (myid == 0) write(*,*) 'Ending coadding'
    
!!$    ! Report reduced chisquares
!!$    do i = 1, map%nfreq
!!$       nu    = count(map%nhit(:,:,i)>0)
!!$       chisq = sum((map%m(:,:,i)/map%rms(:,:,i))**2,map%nhit(:,:,i)>0)
!!$       write(*,fmt='(a,i5,a,f8.3,a,f8.3)') 'Freq = ', i, ', red chisq = ', chisq/nu, ', sigma = ', (chisq-nu)/sqrt(2.d0*nu)
!!$    end do

  !end subroutine compute_scan_maps



  ! subroutine finalize_mapmaking(map)
  !   implicit none
  !   type(map_type), intent(inout) :: map

  !   where(map%nhit > 0)
  !      map%m   = map%dsum / map%div
  !      map%rms = 1.d0 / sqrt(map%div)
  !   elsewhere
  !      map%m   = 0.d0
  !      map%rms = 0.d0
  !   end where    

  ! end subroutine finalize_mapmaking


  ! subroutine free_tod_type(tod)
  !   implicit none
  !   type(tod_type), intent(inout) :: tod 

  !   if (allocated(tod%t))     deallocate(tod%t)
  !   if (allocated(tod%f))     deallocate(tod%f)
  !   if (allocated(tod%d))     deallocate(tod%d)
  !   if (allocated(tod%d_raw)) deallocate(tod%d_raw)
  !   if (allocated(tod%g))     deallocate(tod%g)
  !   if (allocated(tod%rms))   deallocate(tod%rms)
  !   if (allocated(tod%point)) deallocate(tod%point)

  ! end subroutine free_tod_type
   
end program tod2comap
