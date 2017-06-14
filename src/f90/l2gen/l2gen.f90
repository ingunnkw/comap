! Given a runlist, extract CESes from level1 files, producing level2-files
! in hdf format.
program l2gen
  !use fake_mpi_mod
  !use l1_read_mod
  use quiet_utils
  use quiet_fileutils
  !use quiet_detector_mod
  use comap_scan_mod
  use comap_detector_mod
  use comap_Lx_mod
  implicit none

  character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile
  integer(i4b)         :: i, j, k, l, m, n, snum, nscan, unit, myid, nproc, ierr, ndet
  integer(i4b)         :: mstep, i2, phase_offset, decimation, mod, di, nsamp, status
  integer(i4b)         :: debug, num_l1_files
  logical(lgt)         :: exist, reprocess, check_existing, gonext
  real(dp)             :: timing_offset, mjd_tol, mjd(2), dt_error, samprate_in, samprate, scanfreq
  integer(i4b), dimension(:),   allocatable :: detectors
  real(sp),     dimension(:,:), allocatable :: buf
  type(comap_scan_info) :: scan
  !type(point_struct)   :: point_in, point
  !type(data_struct)    :: data_in, data
  !type(hk_struct)      :: hkeep
  !type(hdf_file)       :: file
  type(Lx_struct), allocatable, dimension(:) :: data_l1
  type(Lx_struct)                            :: data_l2

  call getarg(1, parfile)
  call get_parameter(unit, parfile, 'L2_SAMPRATE',         par_dp=samprate)
  call get_parameter(unit, parfile, 'REPROCESS_ALL_FILES', par_lgt=reprocess)
  call get_parameter(unit, parfile, 'DEBUG',               par_int=debug)

  check_existing = .true.
  mjd_tol  = 10d0/60/60/24

  call initialize_scan_mod(parfile)
  call init_detector_mod(parfile)
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call dset(id=myid,level=debug)

  ndet = get_num_detectors()
  allocate(detectors(ndet))
  do i = 1, ndet; detectors(i) = i; end do
  mstep = 10

  nscan = get_num_scans()
  do snum = 1+myid, nscan, nproc
     call get_scan_info(snum, scan)
     inquire(file=scan%l2file,exist=exist)
     if(exist .and. .not. reprocess) then
        if(check_existing) then
           ! Check if the file is as long as it should be
           call get_l2_time_stats(scan%l2file, mjd, dt_error)
           gonext = dt_error < 0.25 .and. (mjd(2) - scan%mjd(2) < mjd_tol .or. scan%mjd(1) - mjd(1) < mjd_tol)
        else
           gonext = .true.
        end if
        if(gonext) then
           write(*,fmt="(i3,a,2i5,a)") myid, " skipping already finished scan:", snum, scan%cid
           cycle
        end if
     end if
     write(*,fmt="(i3,a,i4,a)") myid, " processing scan ", scan%cid, " (" // trim(itoa(snum)) // "/" // trim(itoa(nscan)) // ")"
     call dmem("scan start")

     ! Read in Level 1 data
     num_l1_files = size(scan%l1files)
     allocate(data_l1(num_l1_files))
     do i = 1, num_l1_files
        call read_l1_file(scan%l1files(i), data_l1(i))
     end do

     ! Reformat L1 data into L2 format, and truncate
     !call merge_l1_into_l2_files(scan%mjd, data_l1, data_l2)
     
     ! Write L2 file to disk
     call write_l2_file(scan%l2file, data_l2)

     do i = 1, num_l1_files
        call free_lx_struct(data_l1(i))
     end do
     call free_lx_struct(data_l2)

  end do
  call mpi_finalize(ierr)

contains

  subroutine get_l2_time_stats(filename, mjd, dt_error)
    implicit none
    character(len=*) :: filename
    real(dp)         :: mjd(2), dt_ideal, dt_error, srate
    integer(i4b)     :: n(7)
    type(hdf_file)   :: file
    real(dp),     dimension(:), allocatable :: time
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "time", n)
    allocate(time(n(1)))
    call read_hdf(file, "time", time)
    call read_hdf(file, "samprate", srate)
    dt_ideal = 1000/srate
    dt_error = maxval(abs((time(2:)-time(1:size(time)-1))*24*60*60*1000 - dt_ideal))/dt_ideal
    mjd = [ time(1), time(size(time)) ]
    call close_hdf_file(file)
    deallocate(time)
  end subroutine

end program
