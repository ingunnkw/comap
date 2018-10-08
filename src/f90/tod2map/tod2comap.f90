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
  character(len=512)    :: filename, parfile, acceptfile, prefix, pre, map_name
  integer(i4b)          :: nscan, i, j, k, det, sb, freq
  integer(i4b)          :: myid, nproc, ierr, root
  logical               :: binning_split

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  root = 0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! mpi per scan, all detectors
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1, parfile)

  call initialize_scan_mod(parfile)
  nscan = get_num_scans()
  write(*,*) nscan
  stop
  allocate(tod(nscan))
  !nscan = 1
!  call free_map_type(map)

  call get_parameter(0, parfile, 'MAP_DIR', par_string=pre)
  !call get_parameter()

  binning_split = .false.


  ! This loop currently requiers that all scans are of the same patch

  do i = 1, nscan 
     if (myid == 0) write(*,*) i, 'of', nscan
     if (allocated(alist%status)) deallocate(alist%status)
     call get_scan_info(i,scan)

     ! Get TOD / read level 3 file
     !filename = scan%l3file
     filename = '/mn/stornext/d5/comap/protodir/level3_new/Ka/'//trim(scan%object)//'/'//trim(scan%object)//'_00'//trim(itoa(scan%sid))//'_01.h5'
     if (i == 1) write(*,*) trim(filename)
     call get_tod(trim(filename), tod(i))

     allocate(alist%status(tod(i)%nfreq, tod(i)%ndet))
     alist%status = 0 ! TODO: check if any parts of the scan has been rejected

     !prefix = pre//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID

     if (binning_split) then
        call initialize_mapmaker(map, tod(i:i), parfile)
        call time2pix(tod(i:i), map)
        do det = 3, 3
           do sb = 1, tod(i)%nsb
              do freq = 1, tod(i)%nfreq
                 call binning(map, tod(i:i), alist, det, sb, freq)
              end do
           end do
        end do
        call finalize_binning(map)
        prefix = trim(pre)//trim(scan%object)//'_'//itoa(scan%sid)
        call output_map(trim(prefix), map)
     end if

  end do

  !if (myid == 0) write(*,*) trim(scan%object), trim(itoa(scan%sid))
  !prefix = trim(pre)//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID
  !if (myid == 0) write(*,*) prefix
  !call mpi_finalize(ierr)
  !stop

  call initialize_mapmaker(map, tod, parfile)
  call time2pix(tod, map)

  ! Compute and co-add maps
  !call binning(map, tod(i), alist)
  if (myid == 0) write(*,*) "CG mapmaker ..."
  do det = 3, 3
     do sb = 3,3!1, tod(1)%nsb
        if (myid == 0) write(*,*) "sb", sb
        !do freq = tod(1)%nfreq/4, tod(1)%nfreq/4
        do freq = 30,30
        !do freq = 1, tod(1)%nfreq
           if (myid == 0 .and. modulo(freq, 10) == 0) write(*,*) 'freq', freq, 'of', tod(1)%nfreq
           call pcg_mapmaker(tod, map, alist, det, sb, freq, parfile)
           !if (binning) call binning(map, tod, alist, det, sb, freq)
        end do
     end do
  end do
  !if (binning) call finalize_binning(map)

  !end do
  if (myid == 0) write(*,*) "Writing to file ..."
  !write(*,*) trim(itoa(scan%sid))
  prefix = trim(pre)//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID

  if (myid == 0) call output_map_h5(trim(prefix), map)

  !end do

  if (myid == 0) write(*,*) 'Done'
  do i = 1, nscan
     call free_tod_type(tod(i))
  end do
  deallocate(tod)
  call free_map_type(map)
  call mpi_finalize(ierr)

contains



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

     
end program tod2comap
