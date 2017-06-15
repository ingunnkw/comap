! Given a runlist, extract CESes from level1 files, producing level2-files
! in hdf format.
program l2gen
  use quiet_utils
  use quiet_fileutils
  use comap_scan_mod
  use comap_detector_mod
  use comap_Lx_mod
  implicit none

  character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile
  integer(i4b)         :: i, j, k, l, m, n, snum, nscan, unit, myid, nproc, ierr, ndet
  integer(i4b)         :: mstep, i2, decimation, mod, di, nsamp, status, numfreq
  integer(i4b)         :: debug, num_l1_files
  logical(lgt)         :: exist, reprocess, check_existing, gonext
  real(dp)             :: timing_offset, mjd_tol, mjd(2), dt_error, samprate_in, samprate, scanfreq
  integer(i4b), dimension(:),   allocatable :: detectors
  type(comap_scan_info) :: scan
  type(Lx_struct), allocatable, dimension(:) :: data_l1
  type(Lx_struct)                            :: data_l2_fullres, data_l2_decimated

  call getarg(1, parfile)
  call get_parameter(unit, parfile, 'L2_SAMPRATE',         par_dp=samprate)
  call get_parameter(unit, parfile, 'NUMFREQ',             par_int=numfreq)
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
        gonext = .true.
!!$        if(check_existing) then
!!$           ! Check if the file is as long as it should be
!!$           call get_l2_time_stats(scan%l2file, mjd, dt_error)
!!$           gonext = dt_error < 0.25 .and. (mjd(2) - scan%mjd(2) < mjd_tol .or. scan%mjd(1) - mjd(1) < mjd_tol)
!!$        else
!!$           gonext = .true.
!!$        end if
        if(gonext) then
           write(*,fmt="(i3,a,2i5,a)") myid, " skipping already finished scan:", snum, scan%cid
           cycle
        end if
     else if (exist .and. reprocess) then
        call rm(scan%l2file)
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
     call merge_l1_into_l2_files(scan%mjd, data_l1, data_l2_fullres)

     ! If necessary, decimate L2 file in both time and frequency
     call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)

     ! Write L2 file to disk
     call write_l2_file(scan%l2file, data_l2_decimated)

     ! Clean up data structures
     do i = 1, num_l1_files
        call free_lx_struct(data_l1(i))
     end do
     deallocate(data_l1)
     call free_lx_struct(data_l2_decimated)
     call free_lx_struct(data_l2_fullres)

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


  subroutine merge_l1_into_l2_files(mjd, data_l1, data_l2)
    implicit none
    real(dp),                                   intent(in)  :: mjd(2)
    type(Lx_struct), allocatable, dimension(:), intent(in)  :: data_l1
    type(Lx_struct),                            intent(out) :: data_l2
    
    integer(i4b) :: i, j, k, m, n, nsamp, nsamp_tot, num_l1_files, nfreq, nsb, ndet
    real(dp)     :: samprate
    integer(i4b), allocatable, dimension(:,:) :: ind
    type(Lx_struct) :: data_l2_fullres

    ! Find number of samples
    num_l1_files = size(data_l1)
    allocate(ind(num_l1_files,2))   ! First and last accepted index of each L1 file
    nsamp_tot = 0
    do i = 1, num_l1_files
       nsamp = size(data_l1(i)%tod_l1,1)
       if (i == 1) then
          nfreq    = size(data_l1(i)%tod_l1,2)
          nsb      = size(data_l1(i)%tod_l1,3)
          ndet     = size(data_l1(i)%tod_l1,4)
          samprate = data_l1(i)%samprate
       else
          call assert(size(data_l1(i)%nu_l1,1)  == nfreq, 'Different number of frequencies in L1 files')
          call assert(size(data_l1(i)%tod_l1,3) == nsb, 'Different number of sidebands in L1 files')
          call assert(size(data_l1(i)%tod_l1,4) == ndet, 'Different number of detectors in L1 files')
          call assert(data_l1(i)%samprate == samprate, 'Different sample rates in L1 files')
       end if

       ! Check that there are some valid samples inside current L1 file
       if (mjd(1) > data_l1(i)%time(nsamp) .or. mjd(2) < data_l1(i)%time(1)) then
          ! No acceptable samples
          ind(i,:) = -1
          cycle
       end if

       ! Find start position
       ind(i,1) = 1
       do while (data_l1(i)%time(ind(i,1)) < mjd(1) .and. ind(i,1) <= nsamp)
          ind(i,1) = ind(i,1) + 1
       end do

       ! Find end position
       ind(i,2) = nsamp
       do while (data_l1(i)%time(ind(i,2)) > mjd(2) .and. ind(i,2) >= 1)
          ind(i,2) = ind(i,2) - 1
       end do

       nsamp_tot = nsamp_tot + ind(i,2)-ind(i,1)+1
    end do

    call assert(any(ind(:,1) /= -1) .and. any(ind(:,2) /= -1), 'No valid ranges in L1 files')

    ! Allocate full-resolution L2 structure
    allocate(data_l2%time(nsamp_tot))
    allocate(data_l2%nu(nsb*nfreq))
    allocate(data_l2%tod(nsamp_tot, nsb*nfreq, ndet))
    allocate(data_l2%point_tel(3,nsamp_tot))
    allocate(data_l2%point_cel(3,nsamp_tot))
    allocate(data_l2%flag(nsamp_tot))

    ! Merge L1 data
    data_l2%decimation_time = 1
    data_l2%decimation_nu   = 1
    data_l2%samprate        = samprate
    data_l2%scanmode        = data_l1(1)%scanmode_l1(1)
    j                               = 1
    do i = 1, num_l1_files
       nsamp = ind(i,2)-ind(i,1)+1
       data_l2%time(j:j+nsamp-1)        = data_l1(i)%time(ind(i,1):ind(i,2))
       data_l2%flag(j:j+nsamp-1)        = data_l1(i)%flag(ind(i,1):ind(i,2))
       data_l2%point_tel(:,j:j+nsamp-1) = data_l1(i)%point_tel(:,ind(i,1):ind(i,2))
       data_l2%point_cel(:,j:j+nsamp-1) = data_l1(i)%point_cel(:,ind(i,1):ind(i,2))
       k = 1
       do m = 1, nsb
          do n = 1, nfreq
             if (i == 1) data_l2%nu(k)    = data_l1(i)%nu_l1(n,m)
             data_l2%tod(j:j+nsamp-1,k,:) = data_l1(i)%tod_l1(ind(i,1):ind(i,2),n,m,:)
             k                                    = k+1
          end do
       end do
       call assert(all(data_l1(i)%scanmode_l1 == data_l2%scanmode), 'Varying scanmode within L1 files!')
       j = j + nsamp
    end do

    deallocate(ind)

  end subroutine merge_l1_into_l2_files

  subroutine decimate_L2_data(samprate_out, numfreq_out, data_in, data_out)
    implicit none
    real(dp),        intent(in) :: samprate_out
    integer(i4b),    intent(in) :: numfreq_out
    type(Lx_struct), intent(in)  :: data_in
    type(Lx_struct), intent(out) :: data_out

    integer(i4b) :: i, j, k, l, nsamp_in, nsamp_out, ndet, dt, dnu

    ndet                     = size(data_in%tod,3)
    data_out%samprate        = samprate_out
    dt                       = nint(samprate_out/data_in%samprate)
    data_out%decimation_time = dt
    call assert(data_out%decimation_time >= 1, 'Cannot ask for higher output sample rate than input')

    dnu                    = size(data_in%nu)/numfreq_out
    data_out%decimation_nu = dnu
    call assert(data_out%decimation_nu >= 1, 'Cannot ask for more frequencies than in input files')

    nsamp_out = int(size(data_in%time)/data_out%decimation_time)

    allocate(data_out%time(nsamp_out))
    allocate(data_out%nu(numfreq_out))
    allocate(data_out%tod(nsamp_out, numfreq_out, ndet))
    allocate(data_out%point_tel(3,nsamp_out))
    allocate(data_out%point_cel(3,nsamp_out))
    allocate(data_out%flag(nsamp_out))

    ! Make angles safe for averaging
    call make_angles_safe(data_in%point_tel(:,1), real(2.d0*pi,sp)) ! Phi
    call make_angles_safe(data_in%point_tel(:,3), real(2.d0*pi,sp)) ! Psi
    call make_angles_safe(data_in%point_cel(:,1), real(2.d0*pi,sp)) ! Phi
    call make_angles_safe(data_in%point_cel(:,3), real(2.d0*pi,sp)) ! Psi

    do i = 1, nsamp_out
       data_out%time(i) = mean(data_in%time((i-1)*dt+1:i*dt))  ! Time
       data_out%flag(i) = data_in%time((i-1)*dt+1)             ! Pick first flag in segment
       data_out%point_tel(1,i) = mean(data_in%point_tel(1,(i-1)*dt+1:i*dt)) ! Phi
       data_out%point_tel(2,i) = mean(data_in%point_tel(2,(i-1)*dt+1:i*dt)) ! Theta
       data_out%point_tel(3,i) = mean(data_in%point_tel(3,(i-1)*dt+1:i*dt)) ! Psi
       data_out%point_cel(1,i) = mean(data_in%point_cel(1,(i-1)*dt+1:i*dt)) ! Phi
       data_out%point_cel(2,i) = mean(data_in%point_cel(2,(i-1)*dt+1:i*dt)) ! Theta
       data_out%point_cel(3,i) = mean(data_in%point_cel(3,(i-1)*dt+1:i*dt)) ! Psi

       do k = 1, numfreq_out
          if (i == 1) data_out%nu(k) = mean(data_in%nu((k-1)*dnu+1:k*dnu)) ! Frequency
          do l = 1, ndet           ! Time-ordered data
             data_out%tod(i,k,l) = mean(reshape(data_in%tod((i-1)*dt+1:i*dt,(k-1)*dnu+1:k*dnu,l), &
                  & (/ dt*dnu/)))
          end do
          
       end do
    end do

  end subroutine decimate_L2_data


end program
