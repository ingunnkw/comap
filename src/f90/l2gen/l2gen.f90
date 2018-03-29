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
        call correct_missing_time_steps(data_l1(i)%time_point)
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

  subroutine compute_time_offset(time_el, el, time, tod, offset)
    implicit none
    real(dp), dimension(1:), intent(in)  :: time_el, el, time, tod
    real(dp),                intent(out) :: offset 

    integer(i4b) :: i, j, p, q, n_tod, n_el, start_tod, end_tod, start_el, end_el
    integer(i4b) :: off(1), left, center, right, n_val
    real(dp)     :: mjd_min, mjd_max, dmjd, threshold, max_slew_length, delta_tod, delta_el, delta_cross
    real(dp), allocatable, dimension(:) :: el_spline, tod_int, del, del2, numsteps, time_tod, tod_resamp
    integer(i4b), allocatable, dimension(:) :: slew_sample_id
    type(spline_type)  :: point_spline
    complex(dpc), dimension(:), allocatable :: fft_el, fft_tod

    mjd_min = max(minval(time_el), minval(time))
    mjd_max = min(maxval(time_el), maxval(time))
    dmjd    = mjd_max - mjd_min
    mjd_min = mjd_min + 0.01d0*dmjd
    mjd_max = mjd_max - 0.01d0*dmjd
    start_tod = locate(time,    mjd_min)
    end_tod   = locate(time,    mjd_max)
    start_el  = locate(time_el, mjd_min)
    end_el    = locate(time_el, mjd_max)
    n_el      = end_el  - start_el  + 1
    n_tod     = end_tod - start_tod + 1
    threshold = 0.3d0
    max_slew_length = 2000 ! Maximum number of samples between CES's

    ! Resample pointing on radiometer time grid
    allocate(el_spline(n_tod))
    call spline(point_spline, time_el(start_el:end_el), el(start_el:end_el))
    do i = start_tod, end_tod
       el_spline(i-start_tod+1) = splint(point_spline, time(i))
    end do

    ! Search for jumps, using the derivative of the elevation
    allocate(del(n_tod), del2(n_tod))
    del(1) = el_spline(2) - el_spline(1)
    do i = 2, n_tod-1
       del(i) = 0.5d0*(el_spline(i+1)-el_spline(i-1))
    end do
    del(n_tod) = el_spline(n_tod)-el_spline(n_tod-1)
    del        = abs(del)

    ! Threshold derivative
    where (del > threshold*maxval(del))
       del = 1.d0
    elsewhere
       del = 0.d0
    end where
    call median_filter(del, del2, 10)
    del = del2

    ! Identify valid slews
    allocate(slew_sample_id(10000))
    j = 0
    i = 1
    do while (i < n_tod)
       if (del(i) == 0.d0) then
          i = i+1
       else
          k = 0
          do while (del(i) == 1.d0)
             k = k+1
             i = i+1
             if (i == n_tod) exit
          end do
          if (k <= max_slew_length) then
             if (j == 0) then
                j                 = j+1
                slew_sample_id(j) = i - k/2
             else if (i-k/2-slew_sample_id(j) > 2*max_slew_length) then
                j                 = j+1
                slew_sample_id(j) = i - k/2
             end if
          end if
       end if
    end do

    ! Allocate shifted time and tod arrays 
    allocate(time_tod(n_tod), tod_resamp(n_tod))
    time_tod = time(start_tod:end_tod)
    tod_resamp = tod(start_tod:end_tod)

    ! Search for crude offset, giving equal weight to all points
    allocate(numsteps(-n_tod:n_tod))
    numsteps = 0.d0
    off          = -10000000
    do i = -n_tod, n_tod
       n_val    = 0
       do k = 1, j
          left   = slew_sample_id(k) + i - 0.5*max_slew_length
          center = slew_sample_id(k) + i 
          right  = slew_sample_id(k) + i + 0.5*max_slew_length
          if (left < 1 .or. right > n_tod)     cycle
          if (left-i < 1 .or. right-i > n_tod) cycle
          n_val = n_val+1
          delta_tod   = (tod_resamp(right)-tod_resamp(center))*(tod_resamp(center)-tod_resamp(left))
          delta_el    = (el_spline(right-i)-el_spline(center-i))*(el_spline(center-i)-el_spline(left-i))
          delta_cross = (tod_resamp(right)-tod_resamp(left))*(el_spline(right-i)-el_spline(left-i))
          if (delta_tod > 0.d0 .and. delta_el > 0.d0 .and. delta_cross < 0.d0) then
             numsteps(i) = numsteps(i) + 1.d0
          end if
       end do
       if (n_val > j/2) then 
          numsteps(i) = numsteps(i) / n_val
       else
          numsteps(i) = 0.d0
       end if
    end do

    where (numsteps == maxval(numsteps))
       numsteps = 1.d0
    elsewhere
       numsteps = 0.d0
    end where
    deallocate(del2)
    allocate(del2(-n_tod:n_tod))
    do i = -n_tod, n_tod
       p = max(i-100,-n_tod)
       q = min(i+100,n_tod)
       del2(i) = maxval(numsteps(p:q))
    end do
    numsteps = del2

    p = -n_tod
    do while (numsteps(p) == 0.d0)
       p = p+1
       if (p == n_tod) exit
    end do
       
    q = p
    do while (numsteps(q) == 1.d0)
       q = q+1
       if (q == n_tod) exit
    end do

    off = (q+p)/2

    ! Return offset
    offset = off(1) * (time(2)-time(1))

!!$    open(58,file='test.dat')
!!$    do i = 1, n_tod
!!$       write(58,*) time_tod(i), tod_resamp(i)
!!$    end do
!!$    write(58,*)
!!$    do i = 1, n_tod
!!$       if (i+off(1)>1 .and. i+off(1) < n_tod) then
!!$          write(58,*) time_tod(i), el_spline(i+off(1))
!!$       end if
!!$    end do
!!$    close(58)
!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(el_spline, del, del2, numsteps, slew_sample_id, time_tod, tod_resamp)
    call free_spline(point_spline)

  end subroutine compute_time_offset


  subroutine merge_l1_into_l2_files(mjd, data_l1, data_l2)
    implicit none
    real(dp),                                   intent(in)  :: mjd(2)
    type(Lx_struct), allocatable, dimension(:), intent(in)  :: data_l1
    type(Lx_struct),                            intent(out) :: data_l2
    
    integer(i4b) :: i, j, k, l, m, n, nsamp, nsamp_tot, num_l1_files, nfreq, nsb, ndet, nsamp_point, buffer, ndet0
    real(dp)     :: samprate, offset_mjd, mjd_min, mjd_max
    integer(i4b), allocatable, dimension(:,:) :: ind, ind_point
    type(Lx_struct)   :: data_l2_fullres
    type(spline_type), allocatable, dimension(:) :: point_tel_spline, point_cel_spline

    ! Find number of samples
    num_l1_files = size(data_l1)
    allocate(ind(num_l1_files,2))         ! First and last accepted index of each L1 file
    allocate(ind_point(num_l1_files,2))   ! First and last accepted index of each L1 file in time_point
    nsamp_tot = 0
    do i = 1, num_l1_files

       ! Match pointing time with radiometer time
       !call compute_time_offset(data_l1(i)%time_point, real(data_l1(i)%point_tel(2,:),dp), data_l1(i)%time, &
       !     & real(data_l1(i)%tod_l1(:,1,1,1),dp), offset_mjd)
       offset_mjd = -1.d0 / 3600.d0 / 24.d0

       ! Find basic information
       nsamp       = size(data_l1(i)%tod,1)
       nsamp_point = size(data_l1(i)%time_point)

       if (i == 1) then
          nfreq    = size(data_l1(i)%tod,2)
          nsb      = size(data_l1(i)%tod,3)
          ndet     = size(data_l1(i)%tod,4)
          samprate = data_l1(i)%samprate
       else
          call assert(size(data_l1(i)%nu,1)  == nfreq, 'Different number of frequencies in L1 files')
          call assert(size(data_l1(i)%tod,3) == nsb, 'Different number of sidebands in L1 files')
          call assert(size(data_l1(i)%tod,4) == ndet, 'Different number of detectors in L1 files')
          call assert(data_l1(i)%samprate == samprate, 'Different sample rates in L1 files')
       end if

       mjd_min = max(mjd(1), max(minval(data_l1(i)%time), minval(data_l1(i)%time_point+offset_mjd)))
       mjd_max = min(mjd(2), min(maxval(data_l1(i)%time), maxval(data_l1(i)%time_point+offset_mjd)))

       ! Check that there are some valid samples inside current L1 file
       if (mjd_min > data_l1(i)%time(nsamp) .or. mjd_max < data_l1(i)%time(1)) then
          ! No acceptable samples
          ind(i,:) = -1
          cycle
       end if

       ! Find start position
       ind(i,1) = 1
       do while (data_l1(i)%time(ind(i,1)) < mjd_min .and. ind(i,1) <= nsamp)
          ind(i,1) = ind(i,1) + 1
       end do

       ! Find end position
       ind(i,2) = nsamp
       do while (data_l1(i)%time(ind(i,2)) > mjd_max .and. ind(i,2) >= 1)
          ind(i,2) = ind(i,2) - 1
       end do

       ! Find start position for pointing
       ind_point(i,1) = 1
       do while (data_l1(i)%time_point(ind_point(i,1))+offset_mjd < mjd_min .and. ind_point(i,1) <= nsamp_point)
          ind_point(i,1) = ind_point(i,1) + 1
       end do

       ! Find end position for pointing
       ind_point(i,2) = nsamp_point
       do while (data_l1(i)%time_point(ind_point(i,2))+offset_mjd > mjd_max .and. ind_point(i,2) >= 1)
          ind_point(i,2) = ind_point(i,2) - 1
       end do

       nsamp_tot = nsamp_tot + ind(i,2)-ind(i,1)+1
    
    end do


    call assert(any(ind(:,1) /= -1) .and. any(ind(:,2) /= -1), 'No valid ranges in L1 files')

    ! Allocate full-resolution L2 structure
    allocate(data_l2%time(nsamp_tot))
    allocate(data_l2%nu(nfreq,nsb))
    allocate(data_l2%tod(nsamp_tot, nfreq, nsb, ndet))
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
          call spline(point_tel_spline(l), data_l1(i)%time_point(ind_point(i,1):ind_point(i,2))+offset_mjd, &
               & real(data_l1(i)%point_tel(l,ind_point(i,1):ind_point(i,2)),dp))
          call spline(point_cel_spline(l), data_l1(i)%time_point(ind_point(i,1):ind_point(i,2))+offset_mjd, &
               & real(data_l1(i)%point_cel(l,ind_point(i,1):ind_point(i,2)),dp))
       end do

       nsamp = ind(i,2)-ind(i,1)+1
       data_l2%time(j:j+nsamp-1)        = data_l1(i)%time(ind(i,1):ind(i,2))
       data_l2%flag(j:j+nsamp-1)        = data_l1(i)%flag(ind(i,1):ind(i,2))
       !data_l2%point_tel(:,j:j+nsamp-1) = data_l1(i)%point_tel(:,ind(i,1):ind(i,2))
       !data_l2%point_cel(:,j:j+nsamp-1) = data_l1(i)%point_cel(:,ind(i,1):ind(i,2))
       do l = 1, 3
          do k = j, j+nsamp-1
             data_l2%point_tel(l,k) = splint(point_tel_spline(l), data_l2%time(k))
             data_l2%point_cel(l,k) = splint(point_cel_spline(l), data_l2%time(k))
          end do
       end do

       do m = 1, nsb
          do n = 1, nfreq
             if (i == 1) data_l2%nu(n,m)  = data_l1(i)%nu(n,m)
             data_l2%tod(j:j+nsamp-1,n,m,:) = data_l1(i)%tod(ind(i,1):ind(i,2),n,m,:)
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

    integer(i4b) :: i, j, k, l, nsamp_in, nsamp_out, ndet, dt, dnu, nsb

    nsb                      = size(data_in%tod,3)
    ndet                     = size(data_in%tod,4)
    data_out%samprate        = samprate_out
    dt                       = nint(samprate_out/data_in%samprate)
    data_out%decimation_time = dt
    call assert(data_out%decimation_time >= 1, 'Cannot ask for higher output sample rate than input')

    dnu                    = size(data_in%nu,1)/numfreq_out
    data_out%decimation_nu = dnu
    call assert(data_out%decimation_nu >= 1, 'Cannot ask for more frequencies than in input files')

    nsamp_out = int(size(data_in%time)/data_out%decimation_time)

    allocate(data_out%time(nsamp_out))
    allocate(data_out%nu(numfreq_out,nsb))
    allocate(data_out%tod(nsamp_out, numfreq_out, nsb, ndet))
    allocate(data_out%point_tel(3,nsamp_out))
    allocate(data_out%point_cel(3,nsamp_out))
    allocate(data_out%flag(nsamp_out))

    ! Make angles safe for averaging
    call make_angles_safe(data_in%point_tel(1,:), real(360.d0,sp)) ! Phi
    call make_angles_safe(data_in%point_tel(3,:), real(360.d0,sp)) ! Psi
    call make_angles_safe(data_in%point_cel(1,:), real(360.d0,sp)) ! Phi
    call make_angles_safe(data_in%point_cel(3,:), real(360.d0,sp)) ! Psi

    do i = 1, nsamp_out
       data_out%time(i) = mean(data_in%time((i-1)*dt+1:i*dt))  ! Time
       data_out%flag(i) = data_in%time((i-1)*dt+1)             ! Pick first flag in segment
       data_out%point_tel(1,i) = mean(data_in%point_tel(1,(i-1)*dt+1:i*dt)) ! Phi
       data_out%point_tel(2,i) = mean(data_in%point_tel(2,(i-1)*dt+1:i*dt)) ! Theta
       data_out%point_tel(3,i) = mean(data_in%point_tel(3,(i-1)*dt+1:i*dt)) ! Psi
       data_out%point_cel(1,i) = mean(data_in%point_cel(1,(i-1)*dt+1:i*dt)) ! Phi
       data_out%point_cel(2,i) = mean(data_in%point_cel(2,(i-1)*dt+1:i*dt)) ! Theta
       data_out%point_cel(3,i) = mean(data_in%point_cel(3,(i-1)*dt+1:i*dt)) ! Psi

       do j = 1, nsb
          do k = 1, numfreq_out
             if (i == 1) data_out%nu(k,j) = mean(data_in%nu((k-1)*dnu+1:k*dnu,j)) ! Frequency
             do l = 1, ndet           ! Time-ordered data
                data_out%tod(i,k,j,l) = mean(reshape(data_in%tod((i-1)*dt+1:i*dt,(k-1)*dnu+1:k*dnu,j,l), &
                     & (/ dt*dnu/)))
             end do
          end do
       end do
    end do

  end subroutine decimate_L2_data

  subroutine simulate_gain_data(rng_handle, data)
    implicit none
    type(planck_rng), intent(inout) :: rng_handle
    type(Lx_struct),  intent(inout) :: data

    real(dp)     :: T_0, Tsys, bw, gain, sigma, rms_drift, max_drift_freq
    integer(i4b) :: i, j, k, l, s, n, ind_cut
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
    do k = 1, size(data%tod,4)        ! Detector  
       do s = 1, size(data%tod,3)     ! Sideband
          do j = 1, size(data%tod,2)     ! Frequency
             ! Start with large-scale drifts
             !data%tod(:,j,k) = 36.d0 ! Monopole
             data%tod(:,j,s,k) = 0.d0 ! Monopole
             ffts = 0.d0
             do l = 1, ind_cut
                ffts(l) = sqrt(rms_drift * (real(l,dp)/real(ind_cut,dp))**(-3.d0)) * &
                     & cmplx(rand_gauss(rng_handle),rand_gauss(rng_handle))
             end do
             call fft(drift, ffts, -1)
             !data%tod(:,j,k) = data%tod(:,j,k) + drift
             
             do i = 1, size(data%tod,1)  ! Time
                data%tod(i,j,s,k) = data%tod(i,j,s,k) + T_0 / sin(data%point_tel(2,i)*DEG2RAD)  ! Co-secant model
                !data%tod(i,j,k) = data%tod(i,j,k) + sigma * rand_gauss(rng_handle) ! Radiometer equation noise
                data%tod(i,j,s,k) = data%tod(i,j,s,k) * gain  ! Apply gain
             end do
          end do
       end do
    end do
    deallocate(ffts, drift)

    if (myid == 0) then
       open(58,file='gain_sim.dat')
       do i = 1, size(data%tod,1)
          write(58,*) data%time(i), data%tod(i,1,1,1), data%point_tel(2,i)
          !write(58,*) i, data%tod(i,1,1), data%point_tel(2,i)
       end do
       close(58)
    end if

  end subroutine simulate_gain_data

  subroutine convert_GHz_to_k(data)
    implicit none
    type(Lx_struct), intent(inout)  :: data


  end subroutine convert_GHz_to_k

  subroutine correct_missing_time_steps(time)
    implicit none
    real(dp), dimension(:), intent(inout) :: time

    integer(i4b) :: i, n, step
    real(dp)     :: dt

    ! Check time array
    n  = size(time)
    dt = time(2)-time(1)
    if (dt == 0) then
       write(*,*) 'Error: First time sample is buggy'
       call mpi_finalize(ierr)
       stop
    end if

    ! Check each sample
    do i = 2, n
       if (abs(time(i)-(time(i-1)+dt))/dt > 0.01d0) time(i) = time(i-1) + dt
    end do

  end subroutine correct_missing_time_steps


end program
