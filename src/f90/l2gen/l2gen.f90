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
  use quiet_status_mod
  use comap_gain_mod
  implicit none

  character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile, freqmaskfile, monitor_file_name
  integer(i4b)         :: i, j, k, l, m, n, snum, nscan, unit, myid, nproc, ierr, ndet
  integer(i4b)         :: mstep, i2, decimation, nsamp, numfreq
  integer(i4b)         :: debug, num_l1_files, seed, bp_filter
  logical(lgt)         :: exist, reprocess, check_existing, gonext, norm_with_tp
  real(dp)             :: timing_offset, mjd(2), dt_error, samprate_in, samprate, scanfreq, nu_gain, alpha_gain
  type(comap_scan_info) :: scan
  type(Lx_struct)                            :: data_l1, data_l2_fullres, data_l2_decimated
  type(planck_rng)     :: rng_handle
  type(status_file)    :: status

  call getarg(1, parfile)
  call get_parameter(unit, parfile, 'L2_SAMPRATE',              par_dp=samprate)
  call get_parameter(unit, parfile, 'NUMFREQ',                  par_int=numfreq)
  call get_parameter(unit, parfile, 'REPROCESS_ALL_FILES',      par_lgt=reprocess)
  call get_parameter(unit, parfile, 'DEBUG',                    par_int=debug)
  call get_parameter(unit, parfile, 'SEED',                     par_int=seed)
  call get_parameter(unit, parfile, 'FREQUENCY_MASK',           par_string=freqmaskfile)
  call get_parameter(unit, parfile, 'GAIN_NORMALIZATION_NU',    par_dp=nu_gain)
  call get_parameter(unit, parfile, 'GAIN_NORMALIZATION_ALPHA', par_dp=alpha_gain)
  call get_parameter(unit, parfile, 'BANDPASS_FILTER_ORDER',    par_int=bp_filter)
  call get_parameter(unit, parfile, 'NORMALIZE_GAIN_WITH_TOTAL_POWER', par_lgt=norm_with_tp)

  if (.not. norm_with_tp .and. bp_filter >= 0) then
     write(*,*) 'Error: Not possible to poly-filter data without first normalize with total power!'
     stop
  end if

  check_existing = .true.
  call initialize_scan_mod(parfile)
  call initialize_detector_mod(parfile)
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call init_status(status, 'l2gen_mon.txt'); 
  call update_status(status, 'init')
  !call dset(id=myid,level=debug)


  call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)

  nscan = get_num_scans()
  do snum = 1+myid, nscan, nproc
     call get_scan_info(snum, scan)
     inquire(file=scan%l2file,exist=exist)
     if(exist .and. .not. reprocess) then
        gonext = .true.
        if(gonext) then
           write(*,fmt="(i3,a,2i5,a)") myid, " skipping already finished scan:", snum, scan%sid
           cycle
        end if
     else if (exist .and. reprocess) then
        call rm(scan%l2file)
     end if
     write(*,fmt="(i3,a,i4,a)") myid, " processing scan ", scan%sid, " (" // trim(itoa(snum)) // "/" // trim(itoa(nscan)) // ")"
     call update_status(status, 'scan_start')


     ! Read in Level 1 data
     num_l1_files = size(scan%l1files)
     call read_l1_file(scan%l1files(1), data_l1); call update_status(status, 'read_l1')
     call correct_missing_time_steps(data_l1%time_point)

     ! Initialize frequency mask
     call initialize_frequency_mask(freqmaskfile, numfreq, data_l2_fullres)
     call update_status(status, 'freq_mask')

     ! Reformat L1 data into L2 format, and truncate
     call merge_l1_into_l2_files(scan%mjd, data_l1, data_l2_fullres)
     call update_status(status, 'merge')

     ! Elevation-gain renormalization for circular scans
!     if (circular) call remove_elevation_gain(data_l2_fullres) 
     call remove_elevation_gain(data_l2_fullres) 

     ! Normalize gain
     if (norm_with_tp) then
        call normalize_gain(data_l2_fullres, nu_gain, alpha_gain)
        call update_status(status, 'gain_norm')
     end if

     ! Poly-filter if requested
     call polyfilter_TOD(data_l2_fullres, bp_filter)
     call update_status(status, 'polyfilter')

     ! Fourier transform frequency direction
     !call convert_GHz_to_k(data_l2_fullres(i))

     ! If necessary, decimate L2 file in both time and frequency
     call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
     call update_status(status, 'decimate')
     !write(*,*) 'c'

     ! Replace TOD with simulated data
     if (.false.) call simulate_gain_data(rng_handle, data_l2_decimated)

     ! Write L2 file to disk
     write(*,*) 'Writing ', trim(scan%l2file)
     call write_l2_file(scan%l2file, data_l2_decimated)
     call update_status(status, 'write_l2')

     ! Clean up data structures
     call free_lx_struct(data_l1)
     call free_lx_struct(data_l2_decimated)
     call free_lx_struct(data_l2_fullres)

  end do
  call mpi_finalize(ierr)

  call free_status(status)

contains

  subroutine normalize_gain(data_l2, nu_gain, alpha_gain)
    implicit none
    type(Lx_struct),                            intent(inout) :: data_l2
    real(dp),                                   intent(in)    :: nu_gain, alpha_gain

    integer(i4b) :: i, j, k, l, n, nsamp, nfreq, nsb, ndet
    real(dp)     :: samprate, nu
    real(dp),     allocatable, dimension(:) :: dt
    complex(dpc), allocatable, dimension(:) :: dv

    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    samprate    = data_l2%samprate    
    n           = nsamp+1

    allocate(dt(2*nsamp), dv(0:n-1))
    allocate(data_l2%mean_tp(nfreq,nsb,ndet))
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       write(*,*) 'Normalizing gains for det = ', i
       do j = 1, nsb
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle

             dt(1:nsamp)            = data_l2%tod(:,k,j,i)
             dt(2*nsamp:nsamp+1:-1) = dt(1:nsamp)
             call fft(dt, dv, 1)
             ! Apply lowpass filter
             do l = 0, n-1
                nu = ind2freq(l+1, samprate, n)
                dv(l) = dv(l) * 1.d0/(1.d0 + (nu/nu_gain)**alpha_gain)
             end do
             call fft(dt, dv, -1)
             data_l2%tod(:,k,j,i)     = data_l2%tod(:,k,j,i) / dt(1:nsamp)
             data_l2%mean_tp(k,j,i)   = mean(dt(1:nsamp))
          end do
       end do
    end do
    deallocate(dt, dv)

  end subroutine normalize_gain


  subroutine polyfilter_TOD(data_l2, bp_filter)
    implicit none
    type(Lx_struct),                            intent(inout) :: data_l2
    integer(i4b),                               intent(in)    :: bp_filter

    integer(i4b) :: i, j, k, l, n, nsamp, nfreq, nsb, ndet, p, stat
    real(dp)     :: samprate, nu, mu
    real(dp),     allocatable, dimension(:)   :: x
    real(dp),     allocatable, dimension(:,:) :: T, A

    data_l2%polyorder = bp_filter
    if (bp_filter < 0) return

    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    p           = bp_filter

    allocate(T(nfreq,0:p), A(0:p,0:p), x(0:p))
    allocate(data_l2%tod_poly(nsamp,0:p,nsb,ndet))

    ! Precompute polynomial basis
    do k = 1, nfreq
       mu = max(min(2.d0*real(k-1,dp)/real(nfreq-1,dp)-1.d0,1.d0),-1.d0)
       call get_legendre_polynomials(mu,T(k,:))
    end do

    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          !write(*,*) 'Polyfiltering det, sb = ', i, j

          ! Pre-compute Cholesky factor of coupling matrix for current sideband
          do m = 0, p
             do n = 0, m
                A(m,n) = sum(T(:,m)*T(:,n)*data_l2%freqmask_full(:,j,i))
             end do
          end do
          call dpotrf('L', p+1, A, p+1, stat )

          ! Solve for polynomial coefficients
          !$OMP PARALLEL PRIVATE(k,m,x,l,stat)
          !$OMP DO SCHEDULE(guided)
          do k = 1, nsamp
             do m = 0, p
                x(m) = sum(data_l2%tod(k,:,j,i)*T(:,m)*data_l2%freqmask_full(:,j,i))
             end do
             call dpotrs( 'L', p+1, 1, A, p+1, x, p+1, stat)

             ! Store poly coeffs as separate TOD
             data_l2%tod_poly(k,:,j,i) = x

             ! Subtract fitted polynomial from main TOD
             do l = 1, nfreq
                if (data_l2%freqmask_full(l,j,i) == 0.d0) cycle
                data_l2%tod(k,l,j,i) = data_l2%tod(k,l,j,i) - sum(T(l,:)*x)
             end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
       end do
    end do
    deallocate(T, A, x)

  end subroutine polyfilter_TOD


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
    real(dp),                                   intent(in)    :: mjd(2)
    type(Lx_struct),                            intent(in)    :: data_l1
    type(Lx_struct),                            intent(inout) :: data_l2
    
    integer(i4b) :: i, j, k, l, m, n, nsamp, nsamp_tot, num_l1_files, nfreq, nsb, ndet, nsamp_point, buffer, ndet0, err
    integer(i4b) :: ind(2), ind_point(2)
    real(dp)     :: samprate, mjd_min, mjd_max
    type(Lx_struct)   :: data_l2_fullres
    type(spline_type), allocatable, dimension(:) :: point_tel_spline, point_cel_spline

    ! Find basic information
    nsamp       = size(data_l1%tod,1)
    nsamp_point = size(data_l1%time_point)
    nfreq       = size(data_l1%tod,2)
    nsb         = size(data_l1%tod,3)
    ndet        = size(data_l1%tod,4)
    samprate    = data_l1%samprate
    mjd_min     = max(mjd(1), minval(data_l1%time))
    mjd_max     = min(mjd(2), maxval(data_l1%time))
    
    ! Find start position
    ind(1) = 1
    do while (data_l1%time(ind(1)) < mjd_min .and. ind(1) <= nsamp)
       ind(1) = ind(1) + 1
    end do
    
    ! Find end position
    ind(2) = nsamp
    do while (data_l1%time(ind(2)) > mjd_max .and. ind(2) >= 1)
       ind(2) = ind(2) - 1
    end do
    
    ! Find start position for pointing
    ind_point(1) = 1
    do while (data_l1%time_point(ind_point(1)) < mjd_min .and. ind_point(1) <= nsamp_point)
       ind_point(1) = ind_point(1) + 1
    end do
    
    ! Find end position for pointing
    ind_point(2) = nsamp_point
    do while (data_l1%time_point(ind_point(2)) > mjd_max .and. ind_point(2) >= 1)
       ind_point(2) = ind_point(2) - 1
    end do
    
    nsamp_tot = ind(2)-ind(1)+1

    ! Allocate full-resolution L2 structure
    allocate(data_l2%time(nsamp_tot), stat=err)
    allocate(data_l2%nu(nfreq,nsb,ndet))
    allocate(data_l2%tod(nsamp_tot, nfreq, nsb, ndet))
    allocate(data_l2%point_tel(3,nsamp_tot,ndet))
    allocate(data_l2%point_cel(3,nsamp_tot,ndet))
    !allocate(data_l2%flag(nsamp_tot))

    ! Merge L1 data
    data_l2%decimation_time = 1
    data_l2%decimation_nu   = 1
    data_l2%samprate        = samprate
    data_l2%scanmode        = data_l1%scanmode_l1(1)
    data_l2%nu              = data_l1%nu
    j                       = 1

    allocate(point_tel_spline(3), point_cel_spline(3))

    ! Spline pointing
    nsamp = ind(2)-ind(1)+1
    data_l2%time(j:j+nsamp-1)        = data_l1%time(ind(1):ind(2))
!    data_l2%flag(j:j+nsamp-1)        = data_l1%flag(ind(1):ind(2))
    do i = 1, ndet
       if (.not. is_alive(i)) then
          data_l2%point_tel(:,:,i) = 0.d0
          data_l2%point_cel(:,:,i) = 0.d0
          data_l2%tod(:,:,:,i)     = 0.d0
          cycle
       end if
       do l = 1, 3
          call spline(point_tel_spline(l), data_l1%time_point(ind_point(1):ind_point(2)), &
               & real(data_l1%point_tel(l,ind_point(1):ind_point(2),i),dp))
          call spline(point_cel_spline(l), data_l1%time_point(ind_point(1):ind_point(2)), &
               & real(data_l1%point_cel(l,ind_point(1):ind_point(2),i),dp))
       end do

       do l = 1, 3
          do k = j, j+nsamp-1
             data_l2%point_tel(l,k,i) = splint(point_tel_spline(l), data_l2%time(k))
             data_l2%point_cel(l,k,i) = splint(point_cel_spline(l), data_l2%time(k))
          end do
       end do
    end do

    do j = 1, ndet
       do m = 1, nsb
          do n = 1, nfreq
             data_l2%tod(1:nsamp,n,m,j) = data_l1%tod(ind(1):ind(2),n,m,j)
          end do
       end do
    end do
    deallocate(point_cel_spline, point_tel_spline)

  end subroutine merge_l1_into_l2_files

  subroutine decimate_L2_data(samprate_out, numfreq_out, data_in, data_out)
    implicit none
    real(dp),                          intent(in)  :: samprate_out
    integer(i4b),                      intent(in)  :: numfreq_out
    type(Lx_struct),                   intent(in)  :: data_in
    type(Lx_struct),                   intent(out) :: data_out

    integer(i4b) :: i, j, k, l, m, n, nsamp_in, nsamp_out, ndet, dt, dnu, nsb
    real(dp)     :: w, weight
    real(dp), allocatable, dimension(:,:,:) :: sigmasq

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
    allocate(data_out%nu(numfreq_out,nsb,ndet))
    allocate(data_out%tod(nsamp_out, numfreq_out, nsb, ndet))
    allocate(data_out%point_tel(3,nsamp_out,ndet))
    allocate(data_out%point_cel(3,nsamp_out,ndet))
    allocate(data_out%flag(nsamp_out))
    allocate(data_out%freqmask(numfreq_out,nsb,ndet))
    allocate(data_out%freqmask_full(size(data_in%nu,1,1),nsb,ndet))
    allocate(data_out%mean_tp(size(data_in%nu,1),nsb,ndet))
    allocate(data_out%var_fullres(size(data_in%nu,1),nsb,ndet))
    if (allocated(data_in%mean_tp)) then
       data_out%mean_tp = data_in%mean_tp
    else
       data_out%mean_tp = 0.d0
    end if
    data_out%freqmask      = data_in%freqmask
    data_out%freqmask_full = data_in%freqmask_full

    ! Make angles safe for averaging
    do j = 1, ndet
       if (.not. is_alive(j)) cycle
       call make_angles_safe(data_in%point_tel(1,:,j), real(360.d0,sp)) ! Phi
       call make_angles_safe(data_in%point_tel(3,:,j), real(360.d0,sp)) ! Psi
       call make_angles_safe(data_in%point_cel(1,:,j), real(360.d0,sp)) ! Phi
       call make_angles_safe(data_in%point_cel(3,:,j), real(360.d0,sp)) ! Psi
    end do

    ! Compute variance per frequency channel
    !open(58,file='variance.dat')
    do k = 1, ndet
       if (.not. is_alive(k)) cycle
       do j = 1, nsb
          do i = 1, size(data_in%nu,1)
             data_out%var_fullres(i,j,k) = variance(data_in%tod(:,i,j,k))
             !write(58,*) i, data_out%var_fullres(i,j,k)
          end do
          !write(58,*)
       end do
    end do
    !close(58)
!!$    call mpi_finalize(ierr)
!!$    stop


    
    
    !$OMP PARALLEL PRIVATE(i,j,k,l,n,m,weight,w)
    !$OMP DO SCHEDULE(guided)    
    do i = 1, nsamp_out
       data_out%time(i) = mean(data_in%time((i-1)*dt+1:i*dt))  ! Time
       data_out%flag(i) = data_in%time((i-1)*dt+1)             ! Pick first flag in segment
       do j = 1, ndet
          data_out%point_tel(1,i,j) = mean(data_in%point_tel(1,(i-1)*dt+1:i*dt,j)) ! Phi
          data_out%point_tel(2,i,j) = mean(data_in%point_tel(2,(i-1)*dt+1:i*dt,j)) ! Theta
          data_out%point_tel(3,i,j) = mean(data_in%point_tel(3,(i-1)*dt+1:i*dt,j)) ! Psi
          data_out%point_cel(1,i,j) = mean(data_in%point_cel(1,(i-1)*dt+1:i*dt,j)) ! Phi
          data_out%point_cel(2,i,j) = mean(data_in%point_cel(2,(i-1)*dt+1:i*dt,j)) ! Theta
          data_out%point_cel(3,i,j) = mean(data_in%point_cel(3,(i-1)*dt+1:i*dt,j)) ! Psi
       end do

       do j = 1, nsb
          do k = 1, numfreq_out
             do l = 1, ndet           ! Time-ordered data
                if (.not. is_alive(l)) cycle
                data_out%nu(k,j,l)    = mean(data_in%nu((k-1)*dnu+1:k*dnu,j,l)) ! Frequency
                data_out%tod(i,k,j,l) = 0.d0
                weight                = 0.d0
                do n = (k-1)*dnu+1, k*dnu
                   if (data_out%var_fullres(n,j,l) <= 0) then
                      w = 0.d0
                   else
                      w      = 1.d0 / data_out%var_fullres(n,j,l) * data_in%freqmask_full(n,j,l)
                   end if
                   weight = weight + w !data_in%freqmask_full(n,j,l)
                   do m = (i-1)*dt+1, i*dt
                      data_out%tod(i,k,j,l) = data_out%tod(i,k,j,l) + w * data_in%tod(m,n,j,l) !* data_in%freqmask_full(n,j,l)
                   end do
                end do
                if (weight > 0.d0) then
                   data_out%tod(i,k,j,l) = data_out%tod(i,k,j,l) / weight
                else
                   data_out%tod(i,k,j,l) = 0.d0
                end if
             end do
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Polyfiltered TOD
    data_out%polyorder = data_in%polyorder
    if (data_out%polyorder >= 0) then
       allocate(data_out%tod_poly(nsamp_out, 0:data_out%polyorder, nsb, ndet))
       do l = 1, ndet           
          if (.not. is_alive(l)) then
             data_out%tod_poly(:,:,:,l) = 0.d0
             cycle
          end if
          do j = 1, nsb
             do k = 0, data_out%polyorder
                do i = 1, nsamp_out
                   data_out%tod_poly(i,k,j,l) = mean(data_in%tod_poly((i-1)*dt+1:i*dt,k,j,l))
                end do
             end do
          end do
       end do
    end if

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
                data%tod(i,j,s,k) = data%tod(i,j,s,k) + T_0 / sin(data%point_tel(2,i,k)*DEG2RAD)  ! Co-secant model
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
          write(58,*) data%time(i), data%tod(i,1,1,1), data%point_tel(2,i,1)
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

  subroutine initialize_frequency_mask(freqmaskfile, nfreq, data)
    implicit none
    character(len=*),                                intent(in)    :: freqmaskfile
    integer(i4b),                                    intent(in)    :: nfreq
    type(Lx_struct),                                 intent(inout) :: data
    

    integer(i4b) :: i, j, k, nfreq_full, nsb, ndet, unit, det, sb, freq, dfreq, ierr
    logical(lgt) :: first
    character(len=1024) :: line, val, equal

    ndet = get_num_dets()
    nsb  = get_num_sideband()
    unit = getlun()
    open(unit, file=trim(freqmaskfile), recl=1024)
    first = .true.
    do while (.true.)
       read(unit,'(a)', end=99) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (first) then
          read(line,*) nfreq_full
          allocate(data%freqmask(nfreq,nsb,ndet))
          allocate(data%freqmask_full(nfreq_full,nsb,ndet))
          data%freqmask_full = 1.d0
          first = .false.
       else
          read(line,*) freq, sb, det
          if (det == 0 .and. sb == 0 .and. freq == 0) then
             write(*,*) 'ERROR: All frequencies removed by freqmask!'
             stop
          else if (det == 0 .and. sb == 0) then
             data%freqmask_full(freq,:,:) = 0.d0
          else if (det == 0 .and. freq == 0) then
             data%freqmask_full(:,sb,:) = 0.d0
          else if (sb == 0 .and. freq == 0) then
             data%freqmask_full(:,:,det) = 0.d0
          else if (det == 0) then
             data%freqmask_full(freq,sb,:) = 0.d0
          else if (sb == 0) then
             data%freqmask_full(freq,:,det) = 0.d0
          else if (freq == 0) then
             data%freqmask_full(:,sb,det) = 0.d0
          else 
             data%freqmask_full(freq,sb,det) = 0.d0
          end if
          if (all(data%freqmask_full == 0)) then
             write(*,*) 'ERROR: All frequencies removed by freqmask!'
             stop
          end if
       end if
    end do
99  close(unit)

    dfreq = nfreq_full/nfreq
    if (dfreq == 1) then
       data%freqmask = data%freqmask_full
    else
       data%freqmask = 1.
       do k = 1, ndet
          do j = 1, nsb
             do i = 1, nfreq
                if (all(data%freqmask_full((i-1)*dfreq+1:i*dfreq,j,k) == 0.d0)) data%freqmask(i,j,k) = 0.0
             end do
          end do
       end do
    end if

  end subroutine initialize_frequency_mask

  subroutine remove_elevation_gain(data) 
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: n, ndet, nsb, nfreq, i, j, k, l, m, tmin, tmax, parfile_time
    real(dp)        :: g, a, sigma0, chisq, tsys, tau, dnu, const
    real(dp), dimension(:), allocatable :: el, dat

!    m=a*data%samprate/data%scanfreq(2)      ! Number of tod-samples per gain estimate
    m   = (size(data%time))
    write(*,*) 'Number of samples per gain estimate    =', m
!    n   = (size(data%time)+m-1)/m           ! Number of gain samples
    n   = (size(data%time))/m               ! Number of gain samples
    
    write(*,*) 'Number of gain estimates for this scan =', n
    write(*,*) 'n*m                                    =', n*m
    write(*,*) 'Total number of samples                =', size(data%time)
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)
    write(*,*) nfreq, '= nfreq', nsb, '= nsb', ndet, '= ndet'
    write(*,*) '---------------------------------------------------------'
    !allocate(data%time_gain(n), data%gain(n,nfreq,nsb,ndet))
    allocate(el(m), dat(m))
    !data%time_gain = data%time(::m) ! OBS may reallocate length of time_gain!!!
    !open(13,file='gain.dat')
    !open(14,file='chisq.dat')
    do k = 1, ndet
       do l = 1, nsb
          do j = 1, nfreq
          do i = 1, n
             tmin = (i-1)*m+1
             tmax = i*m
             el  = data%point_tel(2,tmin:tmax,k)
             dat = data%tod(tmin:tmax,j,l,k)
             write(*,*) tmin, tmax, 'min max'
             call estimate_gain(el,dat,g,sigma0,chisq)
             data%tod(tmin:tmax,j,l,k) = data%tod(tmin:tmax,j,l,k) - g*1/(sin(el*pi/180.))
!             write(13,*) g
!             write(13,*) (g-5.6d10)/5.6d10*100
!             write(14,*) chisq
          end do
          end do
       end do
    end do
    !close(13)
    !close(14)
    deallocate(el, dat)
  end subroutine remove_elevation_gain

end program
