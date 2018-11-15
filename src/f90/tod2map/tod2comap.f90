program tod2comap
  use comap_lx_mod
  use comap_map_mod
  use comap_scan_mod
  use comap_acceptlist_mod
  use comap_patch_mod
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
  type(map_type)        :: map_tot, map_scan, buffer
  type(comap_scan_info) :: scan
  type(acceptlist)      :: alist
  type(patch_info)      :: pinfo

  integer(i4b), allocatable, dimension(:,:) :: pixels
  character(len=512)    :: filename, parfile, acceptfile, prefix, pre, map_name, object
  integer(i4b)          :: nscan, i, j, k, det, sb, freq
  integer(i4b)          :: myid, nproc, ierr, root
  logical               :: binning_split, found
  real(dp), allocatable, dimension(:,:) :: offsets
  real(dp)              :: my_x_max, my_y_max

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  root = 0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! mpi per scan, all detectors
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1, parfile)
  call get_parameter(0, parfile, 'MAP_DIR', par_string=pre)
  call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)
  call get_parameter(0, parfile, 'ACCEPTLIST', par_string=acceptfile)
  binning_split = .true.
  !call get_parameter(0, parfile, 'BIN_SPLIT', par_)
  !call get_parameter()
  call initialize_scan_mod(parfile, object)
  call initialize_comap_patch_mod(parfile)
  call initialize_detector_mod(parfile)
  call initialize_accept_list(trim(acceptfile), alist)
  found = get_patch_info(object, pinfo)
  if (.not. found) then
     write(*,*) "Error: patch not found"
     call mpi_finalize(ierr)
     stop
  end if
  nscan = get_num_scans()
  write(*,*) nscan
  !stop
  allocate(tod(nscan))
  !nscan = 1
  !call free_map_type(map)


  !allocate(alist%status(tod(i)%nfreq, tod(i)%ndet))
  !alist%status = 0 ! TODO: check if any parts of the scan has been rejected


  ! This loop currently requiers that all scans are of the same patch
  do i = 1+myid, nscan, nproc 
     write(*,*) myid, i, 'of', nscan
     !if (allocated(alist%status)) deallocate(alist%status)
     call get_scan_info(i,scan)

     ! Get TOD / read level 3 file
     filename = scan%l3file
     !filename = '/mn/stornext/d5/comap/protodir/level3_split/Ka/'//trim(scan%object)//'/'//trim(scan%object)//'_00'//trim(itoa(scan%sid))//'_01.h5'
     !if (i == 1) write(*,*) trim(filename)
     call get_tod(trim(filename), tod(i), parfile)

  end do


!!$    open(58,file='/mn/stornext/d5/comap/protodir/auxiliary/Receiver_offset.dat')
!!$    allocate(offsets(19,5))
!!$    read(58,*)
!!$    read(58,*)
!!$    read(58,*)
!!$    do k = 1, 19
!!$       read(58,*) offsets(k,:)
!!$    end do
!!$    close(58)
!!$    offsets = offsets / 60.
!!$
!!$    my_x_max = maxval(tod(1)%point(1,:,1))
!!$    my_y_max = maxval(tod(1)%point(2,:,1))
!!$    do k = 1, size(tod)
!!$       do j = 1, tod(k)%ndet
!!$          if (j == 9 .or. j == 10) cycle
!!$          do i = 1, size(tod(k)%point(1,:,j))
!!$             if (sqrt((tod(k)%point(1,i,j)-237.9)**2 + (tod(k)%point(2,i,j)+19.44)**2) > 0.1) then
!!$                tod(k)%d(i,:,:,j) = 0.d0
!!$                tod(k)%point(1,i,j) = my_x_max
!!$                tod(k)%point(2,i,j) = my_y_max
!!$             end if
!!$          end do
!!$
!!$          tod(k)%point(1,:,j) = tod(k)%point(1,:,j) + offsets(j,4)
!!$          tod(k)%point(2,:,j) = tod(k)%point(2,:,j) + offsets(j,5)
!!$       end do
!!$    end do


  if (myid==0) write(*,*) "Initialising mapmaker"
  call initialize_mapmaker(map_scan, tod, parfile, pinfo)
  call initialize_mapmaker(map_tot,  tod, parfile, pinfo)
  call initialize_mapmaker(buffer,   tod, parfile, pinfo)

  do i = 1+myid, nscan, nproc
     !prefix = pre//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID

     !if (binning_split) then

        call time2pix(tod(i:i), map_scan)
        call time2pix(tod(i:i), map_tot)
        call nullify_map_type(map_scan)
        write(*,*) myid, "making maps, scan", i
        call binning(map_tot, map_scan, tod(i), alist, i)
        call finalize_scan_binning(map_scan)
        prefix = trim(pre)//trim(scan%object)//'_'//scan%id
        call output_map_h5(trim(prefix), map_scan)
        !call free_map_type(map_scan)
     !end if



  end do

  call mpi_reduce(map_tot%div, buffer%div, size(map_tot%div), MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)
  call mpi_reduce(map_tot%dsum, buffer%dsum, size(map_tot%dsum), MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)


  if (myid == 0) then
     map_tot%div = buffer%div
     map_tot%dsum = buffer%dsum
     write(*,*) 'sum', sum(abs(map_tot%dsum)), sum(abs(map_tot%div))
     write(*,*) "Finalising"
     call finalize_binning(map_tot)
     prefix = trim(pre)//trim(scan%object)
     call output_map_h5(trim(prefix), map_tot)
     call free_map_type(map_tot)
  end if
  call mpi_finalize(ierr)
  stop

  !if (myid == 0) write(*,*) trim(scan%object), trim(itoa(scan%sid))
  !prefix = trim(pre)//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID
  !if (myid == 0) write(*,*) prefix
  !call mpi_finalize(ierr)
  !stop

  !if (.not. binning_split) then
  !call initialize_mapmaker(map_tot, tod, parfile)
  !call time2pix(tod, map)

  ! Compute and co-add maps
  !call binning(map, tod(i), alist)
  ! if (myid == 0) write(*,*) "mapmaker ..."
  ! do det = 3, 3
  !    do sb = 1, tod(1)%nsb
  !       if (myid == 0) write(*,*) "sb", sb
  !       !do freq = tod(1)%nfreq/4, tod(1)%nfreq/4
  !       !do freq = 30,30
  !       do freq = 1, tod(1)%nfreq
  !          if (myid == 0 .and. modulo(freq, 10) == 0) write(*,*) 'freq', freq, 'of', tod(1)%nfreq
  !          !call pcg_mapmaker(tod, map, alist, det, sb, freq, parfile)
  !          !call binning(map_tot, map_scan, tod, alist)
  !       end do
  !    end do
  ! end do
  ! call finalize_binning(map)

  ! !end do
  ! if (myid == 0) write(*,*) "Writing to file ..."
  ! !write(*,*) trim(itoa(scan%sid))
  ! prefix = trim(pre)//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID

  ! if (myid == 0) call output_map_h5(trim(prefix), map)

  !end if

  !end do

  if (myid == 0) write(*,*) 'Done'
  do i = 1, nscan
     call free_tod_type(tod(i))
  end do
  deallocate(tod)
  !call free_map_type(map)
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
