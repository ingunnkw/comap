! Given a runlist, extract CESes from level1 files, producing level2-files
! in hdf format.
program l2gen
  use quiet_utils
  use quiet_fileutils
  use comap_scan_mod
  use comap_detector_mod
  use comap_Lx_mod
  use rngmod
  use quiet_fft_mod
  use spline_1D_mod
  implicit none

  character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile
  integer(i4b)         :: i, j, k, l, m, n, snum, nscan, unit, myid, nproc, ierr, ndet
  integer(i4b)         :: mstep, i2, decimation, mod, di, nsamp, status, numfreq
  integer(i4b)         :: debug, num_l1_files, seed
  logical(lgt)         :: exist, reprocess, check_existing, gonext
  real(dp)             :: timing_offset, mjd_tol, mjd(2), dt_error, samprate_in, samprate, scanfreq
  integer(i4b), dimension(:),   allocatable :: detectors
  type(comap_scan_info) :: scan
  type(Lx_struct), allocatable, dimension(:) :: data_l1
  type(Lx_struct)                            :: data_l2_fullres, data_l2_decimated
  type(planck_rng)     :: rng_handle

  call getarg(1, parfile)
  call get_parameter(unit, parfile, 'L2_SAMPRATE',         par_dp=samprate)
  call get_parameter(unit, parfile, 'NUMFREQ',             par_int=numfreq)
  call get_parameter(unit, parfile, 'REPROCESS_ALL_FILES', par_lgt=reprocess)
  call get_parameter(unit, parfile, 'DEBUG',               par_int=debug)
  call get_parameter(unit, parfile, 'SEED',                par_int=seed)

  check_existing = .true.
  mjd_tol  = 10d0/60/60/24

  call initialize_scan_mod(parfile)
  call initialize_detector_mod(parfile)
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call dset(id=myid,level=debug)

  call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)

  ndet = get_num_dets()
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
           write(*,fmt="(i3,a,2i5,a)") myid, " skipping already finished scan:", snum, scan%sid
           cycle
        end if
     else if (exist .and. reprocess) then
        call rm(scan%l2file)
     end if
     write(*,fmt="(i3,a,i4,a)") myid, " processing scan ", scan%sid, " (" // trim(itoa(snum)) // "/" // trim(itoa(nscan)) // ")"
     call dmem("scan start")

     ! Read in Level 1 data
     num_l1_files = size(scan%l1files)
     allocate(data_l1(num_l1_files))
     do i = 1, num_l1_files
        call read_l1_file(scan%l1files(i), data_l1(i))
     end do

     ! Reformat L1 data into L2 format, and truncate
     call merge_l1_into_l2_files(scan%mjd, data_l1, data_l2_fullres)

     ! Fourier transform frequency direction
     !call convert_GHz_to_k(data_l2_fullres(i))

     ! If necessary, decimate L2 file in both time and frequency
     call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)

     ! Replace TOD with simulated data
     if (.false.) call simulate_gain_data(rng_handle, data_l2_decimated)

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
    
    integer(i4b) :: i, j, k, l, m, n, nsamp, nsamp_tot, num_l1_files, nfreq, nsb, ndet, nsamp_point, buffer
    real(dp)     :: samprate, offset_mjd, buffer_width
    integer(i4b), allocatable, dimension(:,:) :: ind, ind_point
    type(Lx_struct)   :: data_l2_fullres
    type(spline_type), allocatable, dimension(:) :: point_tel_spline, point_cel_spline

    offset_mjd   = 0.00045d0+400.d0*0.02/3600.d0/24.d0 ! 1303, Jupiter
    buffer_width = 0.03

    ! Find number of samples
    num_l1_files = size(data_l1)
    allocate(ind(num_l1_files,2))         ! First and last accepted index of each L1 file
    allocate(ind_point(num_l1_files,2))   ! First and last accepted index of each L1 file in time_point
    nsamp_tot = 0
    do i = 1, num_l1_files
       nsamp       = size(data_l1(i)%tod_l1,1)
       nsamp_point = size(data_l1(i)%time_point)
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

       ! Find start position for pointing
       ind_point(i,1) = 1
       do while (data_l1(i)%time_point(ind_point(i,1)) < mjd(1) .and. ind_point(i,1) <= nsamp_point)
          ind_point(i,1) = ind_point(i,1) + 1
       end do

       ! Find end position for pointing
       ind_point(i,2) = nsamp_point
       do while (data_l1(i)%time_point(ind_point(i,2)) > mjd(2) .and. ind_point(i,2) >= 1)
          ind_point(i,2) = ind_point(i,2) - 1
       end do

       nsamp_tot = nsamp_tot + ind(i,2)-ind(i,1)+1

       buffer         = buffer_width*(ind_point(i,2)-ind_point(i,1))
       ind_point(i,1) = max(1, ind_point(i,1)-buffer)
       ind_point(i,2) = min(ind_point(i,2)+buffer, nsamp_point)
    
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

    allocate(point_tel_spline(3), point_cel_spline(3))
    do i = 1, num_l1_files

       ! Spline pointing
       do l = 1, 3
          open(58,file='time.dat')
          do k = ind_point(i,1), ind_point(i,2)
             write(58,*) data_l1(i)%time_point(k)
          end do
          close(58)
          call spline(point_tel_spline(l), data_l1(i)%time_point(ind_point(i,1):ind_point(i,2)), &
               & real(data_l1(i)%point_tel(l,ind_point(i,1):ind_point(i,2)),dp))
          call spline(point_cel_spline(l), data_l1(i)%time_point(ind_point(i,1):ind_point(i,2)), &
               & real(data_l1(i)%point_cel(l,ind_point(i,1):ind_point(i,2)),dp))
       end do

       nsamp = ind(i,2)-ind(i,1)+1
       data_l2%time(j:j+nsamp-1)        = data_l1(i)%time(ind(i,1):ind(i,2))
       data_l2%flag(j:j+nsamp-1)        = data_l1(i)%flag(ind(i,1):ind(i,2))
       !data_l2%point_tel(:,j:j+nsamp-1) = data_l1(i)%point_tel(:,ind(i,1):ind(i,2))
       !data_l2%point_cel(:,j:j+nsamp-1) = data_l1(i)%point_cel(:,ind(i,1):ind(i,2))
       do l = 1, 3
          do k = j, j+nsamp-1
             data_l2%point_tel(l,k) = splint(point_tel_spline(l), data_l2%time(k)+offset_mjd)
             data_l2%point_cel(l,k) = splint(point_cel_spline(l), data_l2%time(k)+offset_mjd)
          end do
       end do
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

    deallocate(ind, ind_point, point_cel_spline, point_tel_spline)

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
    call make_angles_safe(data_in%point_tel(1,:), real(2.d0*pi,sp)) ! Phi
    call make_angles_safe(data_in%point_tel(3,:), real(2.d0*pi,sp)) ! Psi
    call make_angles_safe(data_in%point_cel(1,:), real(2.d0*pi,sp)) ! Phi
    call make_angles_safe(data_in%point_cel(3,:), real(2.d0*pi,sp)) ! Psi

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

  subroutine simulate_gain_data(rng_handle, data)
    implicit none
    type(planck_rng), intent(inout) :: rng_handle
    type(Lx_struct),  intent(inout) :: data

    real(dp)     :: T_0, Tsys, bw, gain, sigma, rms_drift, max_drift_freq
    integer(i4b) :: i, j, k, l, n, ind_cut
    real(sp),     dimension(:), allocatable :: drift
    complex(spc), dimension(:), allocatable :: ffts

    n         = size(data%time)
    T_0       = 8.d0   ! Atmosphere brightness temperature at zenith in K
    Tsys      = 35.d0  ! System temperature in K
    bw        = 8.d9 / 1024.d0   ! Assume 8 GHz bandwidth and 1024 channels
    gain      = 5d8    ! Gain in ADU/K
    sigma     = Tsys / sqrt(bw/data%samprate)
    rms_drift = 5d-6
    max_drift_freq = 1.d0/3.d0 ! High-frequency cutoff in Hz
    ind_cut   = freq2ind(max_drift_freq, data%samprate, n)

    ! Replace TOD with simulated data
    allocate(ffts(0:n/2), drift(1:n))
    do k = 1, size(data%tod,3)        ! Detector  
       do j = 1, size(data%tod,2)     ! Frequency
          ! Start with large-scale drifts
          !data%tod(:,j,k) = 36.d0 ! Monopole
          data%tod(:,j,k) = 0.d0 ! Monopole
          ffts = 0.d0
          do l = 1, ind_cut
             ffts(l) = sqrt(rms_drift * (real(l,dp)/real(ind_cut,dp))**(-3.d0)) * &
                  & cmplx(rand_gauss(rng_handle),rand_gauss(rng_handle))
          end do
          call fft(drift, ffts, -1)
          !data%tod(:,j,k) = data%tod(:,j,k) + drift
          
          do i = 1, size(data%tod,1)  ! Time
             data%tod(i,j,k) = data%tod(i,j,k) + T_0 / sin(data%point_tel(2,i)*DEG2RAD)  ! Co-secant model
             !data%tod(i,j,k) = data%tod(i,j,k) + sigma * rand_gauss(rng_handle) ! Radiometer equation noise
             data%tod(i,j,k) = data%tod(i,j,k) * gain  ! Apply gain
          end do
       end do
    end do
    deallocate(ffts, drift)

    if (myid == 0) then
       open(58,file='gain_sim.dat')
       do i = 1, size(data%tod,1)
          write(58,*) data%time(i), data%tod(i,1,1), data%point_tel(2,i)
          !write(58,*) i, data%tod(i,1,1), data%point_tel(2,i)
       end do
       close(58)
    end if

  end subroutine simulate_gain_data

  subroutine convert_GHz_to_k(data)
    implicit none
    type(Lx_struct), intent(inout)  :: data


  end subroutine convert_GHz_to_k

end program
