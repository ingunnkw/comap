! l3gen: Produce proper L3 files from L2 files
program l3gen
  use quiet_utils
  use comap_defs
  use quiet_fileutils
  use quiet_system_mod
  use comap_scan_mod
  use comap_detector_mod
  use quiet_mpi_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use comap_Lx_mod
  !use quiet_fft_mod
  use comap_pointing_mod
  !use quiet_noise_estimation_mod
  use powell_mod
  !use quiet_gain_mod
  !use quiet_stat_mod
  !use quiet_patch_detect_mod
  !use quiet_sidelobe_mod
  use quiet_nr_mod
  use quiet_filter_mod
  use quiet_hdf_mod
  use spline_1d_mod
  implicit none

  type info_struct
     integer(i4b)       :: id, nproc
  end type info_struct

  character(len=512)    :: parfile, odir, outfile
  character(len=512)    :: lockfile, tmpfile, point_objs, tilt_file, offset_mask_file, coord_out
  integer(i4b)          :: ierr, snum, nmod, i, j, isys, osys, mod, nside_l3, debug
  integer(i4b)          :: num_corr_bins
  real(dp)              :: fix_highpass_freq_scan, fix_lowpass_freq, t1, t2
  logical(lgt)          :: reprocess, exist, scanmask, inter_module, no_filters, use_templates
  type(task_list)       :: tasks
  type(info_struct)     :: info
  type(comap_scan_info) :: scan
  type(lx_struct)       :: data
  real(sp),     dimension(:,:), allocatable :: powspecs, templates
  complex(spc), dimension(:,:), allocatable :: ffts
  logical(lgt), dimension(:,:), allocatable :: mask

  call getarg(1, parfile)
  call get_parameter(0, parfile, 'OUTPUT_DIR',           par_string=odir)
  call get_parameter(0, parfile, 'REPROCESS_ALL_FILES',  par_lgt=reprocess)
  call get_parameter(0, parfile, 'DEBUG',                par_int=debug)
  call get_parameter(0, parfile, 'L3_FAST',              par_lgt=no_filters, desc=&
       & "Set this to .false.")
  call get_parameter(0, parfile, 'FIX_HIGHPASS_FREQ',    par_dp=fix_highpass_freq_scan, desc=&
       & "Set this to negative to fit")
  call get_parameter(0, parfile, 'FIX_LOWPASS_FREQ',     par_dp=fix_lowpass_freq, desc=&
       & "Set this to negative to fit")
  call get_parameter(0, parfile, 'L2_COORDINATE_SYSTEM', par_string=coord_out, desc=&
       & "Coordinate system for L2 and L3 files; galactic or celestial")

  isys = COORD_CEL
  if (trim(coord_out) == 'galactic') then
     osys = COORD_GAL
  else
     osys = COORD_CEL
  end if

  call wall_time(t1)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, info%id, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, info%nproc, ierr)
  call dset(id=info%id,level=debug)
  !call print_host_mapping

  lockfile = trim(odir) // "/lock.dat"
  call mkdirs(trim(lockfile), .true.)

  call initialize_scan_mod(parfile);             call dmem("scan mod")
  !call initialize_noise_estimation_mod(parfile); call dmem("noise mod")
  call initialize_comap_pointing_mod(parfile);   call dmem("pointing mod")
  !call initialize_gain_mod(parfile);             call dmem("gain mod")
  !call initialize_patch_detect_mod(parfile);     call dmem("patch detect mod")
  !call initialize_filter_mod(parfile);           call dmem("filter mod")

  ! Process all CES's
  call init_task_list(tasks, lockfile, get_num_scans(), MPI_COMM_WORLD)
  do while(get_next_task(tasks, snum))
     call get_scan_info(snum, scan)
     tmpfile = trim(scan%l3file) // ".part"
     inquire(file=tmpfile,exist=exist)
     if(exist) then
        write(*,fmt="(i3,a,2i5)") info%id, " found incomplete run:", snum, scan%cid
        call rm(tmpfile)
     end if
     inquire(file=scan%l3file,exist=exist)
     if(exist .and. .not. reprocess) then
        write(*,fmt="(i3,a,2i5,a)") info%id, " skipping already finished scan:", snum, scan%cid
        cycle
     end if
     inquire(file=scan%l2file,exist=exist)
     if(.not. exist) then
        write(*,fmt="(i3,a,2i5,a)") info%id, " data missing for scan:", snum, scan%cid, " " // trim(scan%l2file)
        cycle
     end if
     write(*,fmt="(i3,a,i4,a)") info%id, " processing scan ", scan%cid, " (" // trim(itoa(snum)) // "/" // &
          & trim(itoa(get_num_scans())) // ")"

     call dmem("scan start")

     ! Read L2 data
     call read_L2_file(scan%l2file, data)     ; call dmem("read l2")

     ! Process it
     call calc_point (data,  isys, osys)      ; call dmem("point")
     !call calc_objrel(data,  point_objs)      ; call dmem("objrel")
     !call calc_pixels(data, nside_l3)         ; call dmem("pixels")
     call calc_scanfreq(data);                ; call dmem("scanfreq")
     !call subtract_offset(data, offset_mask_file) ; call dmem("offset")
     !call calc_mask(data, osys, mask);        ; call dmem("mask")
     !call calc_noise(data, mask, ffts);       ; call dmem("noise")
     !call calc_fourier(data, ffts, powspecs)  ; call dmem("fourier")
 
     !call fit_noise(data, powspecs, cnum)     ; call dmem("noise")
     !call calc_corr(data, ffts)               ; call dmem("corr")
     !call calc_apex(data)                     ; call dmem("apex")
     !call calc_gain(data)                     ; call dmem("gain")
     !call calc_diode_stats(data, powspecs)    ; call dmem("diode_stats")
     !call calc_stats(data)                    ; call dmem("stats")
     !allocate(data%filter_par(size(data%tod,2),NUM_FILTER_PAR))
     !data%filter_par = -1
     !call calc_az_filter_par(data)            ; call dmem("az_par")
     !call calc_bandpass_filter_par(data, powspecs)     ; call dmem("bandpass")

     ! Output L3 data
     call mkdirs(tmpfile, .true.)
     call write_L3_file(tmpfile, data)        ; call dmem("write l3")

     call mkdirs(trim(scan%l3file), .true.)
     call mv(tmpfile, scan%l3file)

     !deallocate(ffts, powspecs)
     call free_lx_struct(data)
  end do

  call free_task_list(tasks)
  write(*,fmt="(i3,a)") info%id, " finished"
  
  call mpi_finalize(ierr)
  call wall_time(t2)

  write(*,*) 'Time elapsed: ', t2-t1
contains

!!$  subroutine subtract_offset(data, maskfile)
!!$    implicit none
!!$    type(lx_struct)                             :: data
!!$    character(len=512)                          :: maskfile
!!$    integer(i4b)                                :: nside, ndi, nsamp, di, nmod, nhit, nmasked, i, pix, m, numdi, ordering
!!$    logical(lgt)                                :: exist
!!$    real(dp), dimension(:),   allocatable   :: mask_samp
!!$    real(dp), dimension(:,:), allocatable   :: mask
!!$    real(dp), dimension(:),   allocatable   :: offset
!!$
!!$    ndi = size(data%tod,2)
!!$    nsamp = size(data%tod,1)
!!$    nmod = get_num_modules()
!!$    nhit = size(data%pixels,1)
!!$    numdi = get_num_diodes()
!!$
!!$    inquire(file=trim(maskfile), exist=exist)
!!$    if (exist) then
!!$       call read_map(mask, ordering, trim(maskfile))
!!$       nside = npix2nside(size(mask,1))
!!$    else
!!$       nside = 1
!!$       allocate(mask(0:12*nside**2-1,3))
!!$       mask = 1.d0
!!$       ordering = 1
!!$    end if
!!$
!!$    allocate(offset(ndi))
!!$    allocate(mask_samp(nsamp))
!!$
!!$    do m = 0,nmod-1
!!$       do i = 1,nsamp
!!$          pix = ang2pix(nside, ordering, real(data%point(2,i,m),dp), real(data%point(1,i,m),dp))
!!$          mask_samp(i) = mask(pix,1)
!!$       end do
!!$       
!!$!       write(*,*) 'ratio = ', sum(mask_samp)/nsamp
!!$       if (sum(mask_samp) > 0.1d0 * nsamp) then
!!$          do i = 1, numdi
!!$             !di = diode_rel2abs(m,i)
!!$             di = m*numdi + i
!!$             offset(di) = sum(data%tod(:,di)*mask_samp)/real(sum(mask_samp),dp)
!!$          end do
!!$       else
!!$          write(*,*) 'Warning: Too few unmasked samples to subtract offset, module', trim(itoa(m))
!!$          offset(m*numdi+1:(m+1)*numdi) = 0.d0
!!$       end if
!!$    end do
!!$    
!!$    do di=1,ndi
!!$       data%tod(:,di) = data%tod(:,di) - offset(di)
!!$    end do
!!$
!!$    if(allocated(mask)) deallocate(mask)
!!$    deallocate(offset, mask_samp)
!!$  end subroutine subtract_offset


  subroutine calc_point(data, isys, osys)
    implicit none
    type(lx_struct) :: data
    integer(i4b) :: i, nsamp, mod, nmod, isys, osys
    real(dp)     :: op(3), np(3), mat(3,3)
    nsamp = size(data%tod,1)
    allocate(data%point(3,nsamp))
    do i = 1, nsamp
      op = data%point_cel(:,i)
      call swap_coordinate_convention(op(1), op(2), op(3), isys)
      call coord_convert(isys, op(1), op(2), op(3), osys, np(1), np(2), np(3), &
       & mjd=data%time(i), euler=mat)
      data%point(:,i) = np
    end do
    data%coord_sys = osys

    ! Find map extent
    call make_angles_safe(data%point(:,1),real(2.d0*pi,sp))
    data%point_lim(1) = minval(data%point(:,1))
    data%point_lim(2) = maxval(data%point(:,1))
    data%point_lim(3) = minval(data%point(:,2))
    data%point_lim(4) = maxval(data%point(:,2))
    data%point(:,1)   = mod(data%point(:,1),2.d0*pi)
  end subroutine

  subroutine calc_scanfreq(data)
    implicit none
    type(lx_struct) :: data
    call get_scanfreq(data%point_tel(1,:), data%samprate, data%scanfreq(1))
    call get_scanfreq(data%point_tel(2,:), data%samprate, data%scanfreq(2))
  end subroutine calc_scanfreq

!!$  subroutine fit_noise(data, powspecs, cnum)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    integer(i4b) :: ndi, i, cnum
!!$    character(len=4) :: myid_text
!!$    real(sp)     :: powspecs(:,:)
!!$    real(dp)     :: chisq
!!$    type(quiet_ces_info) :: ces
!!$    ndi   = size(data%tod,2)
!!$    allocate(data%sigma0(ndi), data%alpha(ndi), data%fknee(ndi))
!!$    call get_ces_info(cnum, ces)
!!$    call int2string(info%id, myid_text)
!!$!TMR debug outputting to figure out the bias
!!$!    open(42,file='chisq_andstuff_first.dat')
!!$    do i = 1, ndi
!!$       if (is_alive(i)) then
!!$          if(no_filters) then
!!$             data%sigma0(i) = 1e-5; data%alpha(i) = -1; data%fknee(i) = 0.02
!!$          else
!!$!             call fit_1overf_profile_old(data%samprate, data%scanfreq, 1d-2, data%sigma0(i), &
!!$!               & data%alpha(i), data%fknee(i), tod_ps=real(powspecs(:,i),dp), &
!!$!               & cnum=ces%cid, diode=i, chisq_out=chisq, apply_scanmask=scanmask)
!!$!             write(*,*) i, data%sigma0(i), data%alpha(i), data%fknee(i)
!!$!             write(*,*) 'chisq1 = ', chisq
!!$             call fit_1overf_profile(data%samprate, data%scanfreq, 1d-2, data%sigma0(i), &
!!$               & data%alpha(i), data%fknee(i), tod_ps=real(powspecs(:,i),dp), &
!!$               & cnum=ces%cid, diode=i, chisq_out=chisq, apply_scanmask=scanmask)
!!$!             write(*,*) i, data%sigma0(i), data%alpha(i), data%fknee(i)
!!$!write(42,"(i4,f14.5,g14.5,2f14.5)") i, chisq, data%sigma0(i), data%fknee(i), data%alpha(i)
!!$
!!$!             write(*,*) 'chisq2 = ', chisq
!!$!             stop
!!$          end if
!!$       else
!!$          data%sigma0(i) = 0.d0
!!$          data%alpha(i)  = 0.d0
!!$          data%fknee(i)  = 0.d0
!!$       end if
!!$    end do
!!$!    close(42)
!!$!    stop
!!$  end subroutine fit_noise



  subroutine get_scanfreq(point, samprate, scanfreq)
    implicit none
    real(sp)                                :: point(:)
    real(dp)                                :: scanfreq, samprate
    real(sp),     dimension(:), allocatable :: tod, pow
    complex(spc), dimension(:), allocatable :: ft
    integer(i4b)        :: n, m, i
    n = size(point) ! min(10000, size(az))
    m = n/2+1
    allocate(tod(n), ft(m), pow(m))
    tod = point(1:n)
    call fft(tod, ft, 1)
    pow = ft*conjg(ft)
    i = maxloc(pow(2:),1)+1
    scanfreq = ind2freq(i, samprate, m)
    deallocate(tod, ft, pow)
  end subroutine

!!$  subroutine calc_gain(data)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    integer(i4b)    :: n, ndi, i, j, m
!!$    real(dp), dimension(:), allocatable :: tmp
!!$    m   = 100
!!$    n   = (size(data%time)+m-1)/m
!!$    ndi = size(data%tod,2)
!!$    allocate(data%time_gain(n), data%gain(n,ndi), tmp(n))
!!$    data%time_gain = data%time(::m)
!!$    do i = 1, ndi
!!$       call get_gains(data%time_gain, i, tmp)
!!$       data%gain(:,i) = tmp
!!$    end do
!!$    deallocate(tmp)
!!$  end subroutine
!!$
!!$  subroutine calc_stats(data)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    integer(i4b)    :: i, j, k, m, n
!!$    real(dp)        :: mjd
!!$    data%stats = NaN
!!$    mjd        = (data%time(1)+data%time(size(data%time)))/2
!!$    data%stats(STAT_MJD)             = mjd
!!$    data%stats(STAT_LST)             = mjd2lst(mjd, QUIET_GEODETIC_LONGITUDE)
!!$    data%stats(STAT_AZ)              = average_ang(real(data%orig_point(1,:),dp))
!!$    data%stats(STAT_EL)              = average_ang(real(data%orig_point(2,:),dp))
!!$    data%stats(STAT_DK)              = average_ang(real(data%orig_point(3,:),dp))
!!$    data%stats(STAT_PWV)             = mean(real(data%hk%apex%value(:,APEX_PWV),dp))
!!$    data%stats(STAT_PWV_CHANGE)      = sqrt(variance(real(data%hk%apex%value(:,APEX_PWV),dp)))
!!$    data%stats(STAT_HUMIDITY)        = mean(real(data%hk%apex%value(:,APEX_HUMIDITY),dp))
!!$    data%stats(STAT_HUMIDITY_CHANGE) = sqrt(variance(real(data%hk%apex%value(:,APEX_HUMIDITY),dp)))
!!$    data%stats(STAT_WIND)            = mean(real(data%hk%apex%value(:,APEX_WIND_SPEED),dp))
!!$    data%stats(STAT_T_AMBIENT)       = mean(real(data%hk%apex%value(:,APEX_TEMPERATURE),dp))
!!$    data%stats(STAT_TENC)            = mean(real(data%hk%encl%value(:,P3T09),dp))
!!$    data%stats(STAT_TENC_CHANGE)     = sqrt(variance(real(data%hk%encl%value(:,P3T09),dp)))
!!$    data%stats(STAT_CRYO)            = mean(real(data%hk%cryo%value(:,P2T19),dp))
!!$    data%stats(STAT_CRYO_CHANGE)     = sqrt(variance(real(data%hk%cryo%value(:,P2T19),dp)))
!!$    !data%stats(STAT_BIAS_ELECTRONICS_TEMERATURE and CHANGE) = NaN
!!$    !data%stats(STAT_SUN_SIDELOBE and MOON) = NaN
!!$    !data%stats(STAT_SUN_DISTANCE and MOON) = NaN
!!$  end subroutine calc_stats
!!$
!!$  subroutine calc_diode_stats(data, powspecs)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    real(sp)        :: powspecs(:,:)
!!$    integer(i4b)    :: i, j, k, m, n, ndi, obj, pix
!!$    logical(lgt)    :: alive
!!$    real(dp)        :: mjd, avgpt(3), phi, theta, psi, va(3), v(3), b2h(3,3), mat(3,3)
!!$    real(dp)        :: t1, p1
!!$    logical(lgt), dimension(:), allocatable :: sidelobe_hit
!!$
!!$    ndi = size(data%tod(1,:))
!!$    allocate(data%diode_stats(ndi,NUM_DIODE_STATS))
!!$    data%diode_stats = NaN
!!$    do i = 1, ndi
!!$       data%diode_stats(i,:) = 0.d0
!!$       if (is_alive(i)) then
!!$          !print *, "In here is ", i
!!$          call compute_typeB_chisq(data%samprate, data%sigma0(i), data%tod(:,i), data%tp(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_TYPEB))
!!$          call compute_weather_stat(data%samprate, 10.d0, data%tp(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_WEATHER1))  ! 10 sec
!!$          call compute_weather_stat(data%samprate, 30.d0, data%tp(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_WEATHER2))  ! 30 sec
!!$          call jump_finder(1001, i, data%tod(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_JUMP)) ! check for jumps
!!$          !call tp_rms_finder(1000, i, data%tp(:,i), data%diode_stats(i,DIODE_STAT_TP_RMS))
!!$          call compute_single_fft_chisq(data%samprate, 9.9d0, 10.1d0, data%sigma0(i), powspecs(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_10HZ))
!!$          call compute_single_fft_chisq(data%samprate, 1.19d0, 1.21d0, &
!!$               & data%sigma0(i), powspecs(:,i), data%diode_stats(i,DIODE_STAT_1_2HZ))
!!$          call compute_single_fft_chisq(data%samprate, data%scanfreq-0.001d0, data%scanfreq+0.001d0, &
!!$               & data%sigma0(i), powspecs(:,i), data%diode_stats(i,DIODE_STAT_SSS))
!!$          ! Check biases for non-gaussianity/outliers
!!$          !call compute_bias_stat(data%hk%bias%value, data%diode_stats(i,DIODE_STAT_BIAS)) 
!!$          data%diode_stats(i,DIODE_STAT_SIGMA0)   = data%sigma0(i)
!!$          data%diode_stats(i,DIODE_STAT_ALPHA)    = data%alpha(i)
!!$          data%diode_stats(i,DIODE_STAT_FKNEE)    = data%fknee(i)
!!$          data%diode_stats(i,DIODE_STAT_GAIN)     = real(sum(data%gain(:,i),1)/size(data%gain(:,i),1),dp)
!!$       end if
!!$    end do
!!$    ! We also want some sidelobe stats: Sidelobe elevation and
!!$    ! sidelobe sun hit. The sun is assumed to be the first object
!!$    ! in the object-relative coordinates.
!!$
!!$    ! Find the sidelobe hits
!!$    allocate(sidelobe_hit(0:size(quiet_horns)-1))
!!$    sidelobe_hit = .false.
!!$    data%diode_stats(:,DIODE_STAT_SIDELOBE_HIT) = 0
!!$    obj = 1
!!$    do j = 1, size(data%point_objrel,2)
!!$       pix = ang2pix(sidelobes%nside, sidelobes%order, &
!!$        & real(data%point_objrel(2,j,obj),dp), real(data%point_objrel(1,j,obj),dp))
!!$       k = lookup_sidelobe_range(sidelobes, data%time(j))
!!$       sidelobe_hit = sidelobe_hit .or. sidelobes%masks(pix,:,k)
!!$    end do
!!$    where(sidelobe_hit(quiet_diodes%horn)) data%diode_stats(:,DIODE_STAT_SIDELOBE_HIT) = 1
!!$    deallocate(sidelobe_hit)
!!$
!!$    ! Find sidelobe horizontal pointing
!!$    call calc_average_pointing(data, avgpt)
!!$    call swap_coordinate_convention(avgpt(1), avgpt(2), avgpt(3))
!!$    k   = lookup_sidelobe_range(sidelobes, data%time(1))
!!$    b2h = rot_boresight2hor(data%time(1), avgpt(1), avgpt(2), avgpt(3), -1, -1)
!!$    data%diode_stats(:,DIODE_STAT_SIDELOBE_AZ) = 0
!!$    data%diode_stats(:,DIODE_STAT_SIDELOBE_EL) = 0
!!$    do mod = 0, size(quiet_horns)-1
!!$       ! For each pixel in the sidelobe, transform to telescope coordinates
!!$       if(isnan(sidelobes%centers(1,mod,k))) then
!!$          phi = nan; theta = nan
!!$       else
!!$          call vec2ang(sidelobes%centers(:,mod,k), theta, phi)
!!$          mat = matmul(b2h,angles2rot(phi, theta, 0d0))
!!$          call rot2angles(mat, phi, theta, psi)
!!$          call swap_coordinate_convention(phi, theta, psi)
!!$       end if
!!$       data%diode_stats(quiet_horns(mod)%diodes,DIODE_STAT_SIDELOBE_EL) = theta
!!$       data%diode_stats(quiet_horns(mod)%diodes,DIODE_STAT_SIDELOBE_AZ) = phi
!!$!call vec2ang(sidelobes%centers(:,mod,k), t1, p1)
!!$!write(*,'(i4,6f9.2)') mod, avgpt(1:2)*RTOD, p1*RTOD, t1*RTOD, phi*RTOD, theta*RTOD
!!$    end do
!!$  end subroutine calc_diode_stats
!!$
!!$
!!$  ! Find the parameters for the bandpass filter. At the lower end,
!!$  ! we need to worry about the scanning frequency and its multiples.
!!$  ! We want to go as low as possible while maintaining a good chisquare,
!!$  ! and also want to ensure that the filter has finished falling by the
!!$  ! time we reach a harmonic with nasty stuff in it. We will therefore
!!$  ! go in steps 0.5f, 1.5f, 2.5f, 3.5f, ...
!!$  ! We take the spikefilter into account, so that we do not cut because
!!$  ! of spikes which will be filtered anyway.
!!$  subroutine calc_bandpass_filter_par(data, powspecs)
!!$    implicit none
!!$    type(lx_struct)              :: data
!!$    real(sp)                     :: powspecs(:,:)
!!$    integer(i4b)                 :: i, j, k, numbin, ndi, numsamp, nfreq, ind1, ind2, n
!!$    integer(i4b)                 :: delta, n_smooth, delta_high, nd, delta_harm, nharm, indh
!!$    integer(i4b)                 :: nlim
!!$    integer(i4b),allocatable     :: pos(:,:), mask(:)
!!$    real(dp)                     :: dnu, chisq0, nu_min, nu_max, nu, cum, sigma_lim
!!$    real(dp)                     :: acceptlim, highlow_boundary
!!$    real(dp), allocatable, dimension(:) :: chisq, N_fft, chisq_cumdev, chisq_harm, filter, chisq_scan, spikefilter
!!$    type(filter_params) :: filter_opts
!!$!real(dp), allocatable :: cumdevs(:,:), harms(:,:)
!!$!type(hdf_file)        :: hfile
!!$
!!$    if(no_filters) return
!!$    ndi        = size(data%tod,2)
!!$    numsamp    = size(data%tod,1)
!!$    nfreq      = size(powspecs,1)
!!$    dnu        = ind2freq(2, data%samprate, nfreq)
!!$    delta_harm = 19
!!$    sigma_lim  = 4
!!$    acceptlim  = 0.75
!!$    delta_high = 100 ! Smoothing scale for high-frequency chi-squares. Should be even
!!$    ! We do not want the filters to overlap, so we restrict highpass to
!!$    ! be below highlow_boundary and lowpass to be above this.
!!$    highlow_boundary = 2.0
!!$    nlim             = freq2ind(highlow_boundary, data%samprate, nfreq)
!!$    nharm            = nlim/freq2ind(data%scanfreq, data%samprate, nfreq)-1
!!$    allocate(chisq_cumdev(nharm), chisq(nfreq), N_fft(nfreq), chisq_harm(nharm))
!!$    allocate(pos(ndi,2))
!!$    pos = -1
!!$!allocate(harms(nharm,ndi),cumdevs(nharm,ndi))
!!$    
!!$    ! Load the filter to check for spike filtering
!!$    call get_default_filter(filter_opts)
!!$    allocate(mask(nfreq))
!!$    mask = 1
!!$
!!$    if(filter_opts%apply_spike) then
!!$       allocate(spikefilter(nfreq))
!!$       spikefilter = 1.d0
!!$       call spike_filter_fft(.false., filter_opts%spikefreqs, filter_opts%num_spikes, numsamp, data%samprate, nfreq, spikefilter)
!!$       where(spikefilter < 1.d0) mask = 0 ! Mask skips all frequencies with spikes, if any
!!$    end if
!!$
!!$    do i = 1, ndi
!!$       if (data%sigma0(i) == 0) cycle
!!$       call get_inv_noise_filter_fft(.false., numsamp, data%samprate, &
!!$        & data%sigma0(i), data%fknee(i), data%alpha(i), N_fft, filter_opts%apply_highpass)
!!$       ! This is based on complex numbers, so each entry is a chisquare
!!$       ! with 2 degrees of freedom.
!!$
!!$       ! Just mask out the spikefiltered frequencies everywhere.
!!$       chisq(2:nfreq) = powspecs(2:nfreq,i) / N_fft(2:nfreq)
!!$       chisq = chisq*mask ! Masking spikefilter frequencies - if apply_spike=.false., this doesn't do anything
!!$
!!$       cum = 0
!!$       do j = nharm, 1, -1
!!$          ind1 = freq2ind((j-0.5)*data%scanfreq, data%samprate, nfreq)
!!$          ind2 = freq2ind((j+0.5)*data%scanfreq, data%samprate, nfreq)
!!$          indh = freq2ind(      j*data%scanfreq, data%samprate, nfreq)
!!$!          n   = size(chisq(ind1:nlim)) ! Version without spikefilter
!!$          n= sum(mask(ind1:nlim)) ! If there are spikefiltered frequencies, they are not counted.
!!$
!!$          if(j == nharm) then
!!$             cum = sum(chisq(ind1:nlim))   ! And I set them to zero in the chisq array, so they do not count here either
!!$          else
!!$             cum = cum + sum(chisq(ind1:ind2-1))
!!$          end if
!!$          ! Number of sigmas away from expected chisquare in gaussian approx.
!!$          ! DET ER SÅNN SOM DETTE HER ALTSÅ!
!!$          chisq_cumdev(j) = (cum - n)/sqrt(1d0*n)
!!$          ! Find the same thing just around the current harmonic
!!$          ind1  = max(1,min(nfreq,indh-delta_harm))
!!$          ind2  = max(1,min(nfreq,indh+delta_harm))
!!$
!!$          !n     = ind2-ind1+1     ! Version without spikefilter
!!$          n = sum(mask(ind1:ind2)) ! taking mask into account again
!!$          chisq_harm  (j) = (sum(chisq(ind1:ind2))-n)/sqrt(1d0*n)
!!$       end do
!!$
!!$       ! Accept the first position where the cum chisq is ok and all
!!$       ! following harm chisqs are also ok
!!$       do j = nharm, 2, -1
!!$          if(chisq_harm(j) > sigma_lim) exit
!!$       end do
!!$       do k = j, nharm-1
!!$          if(chisq_cumdev(k) < sigma_lim) exit
!!$       end do
!!$       pos(i,1) = k 
!!$!harms(:,i) = chisq_harm
!!$!cumdevs(:,i) = chisq_cumdev
!!$
!!$       ! Then handle the high frequency part of the filter
!!$       ind1   = nlim
!!$       ind2   = nfreq ! nint(12.5 / dnu)
!!$       
!!$       ! Search first for thin spikes
!!$       j = ind1
!!$       do while (chisq(j) < 15.d0 .and. j < ind2)
!!$          j = j+1
!!$       end do
!!$
!!$       ! Search for extended features
!!$       do k = max(ind1,delta_high/2+1), min(ind2,nfreq-delta_high/2)
!!$          chisq0 = sum(chisq(k-delta_high/2:k+delta_high/2))
!!$!          chisq0 = (chisq0 - (delta_high+1)) / sqrt(delta_high+1.d0) ! Version wthout spikefilter
!!$          n = sum(mask(k-delta_high/2:k+delta_high/2)) 
!!$          chisq0 = (chisq0 - n) / sqrt(1.d0*n)         
!!$          if (chisq0 > 5.d0) exit
!!$       end do
!!$       pos(i,2) = min(j,k)
!!$
!!$    end do
!!$    ! The first index (0.5*scanfreq) is extra suspicious. We only
!!$    ! accept it if most of the others also accept it.
!!$    if(count(pos(:,1)==1) < ndi*acceptlim) where(pos(:,1)==1) pos(:,1) = 2
!!$!call open_hdf_file("foo.hdf", hfile, "w")
!!$!call write_hdf(hfile, "harm", harms)
!!$!call write_hdf(hfile, "cum",  cumdevs)
!!$!call close_hdf_file(hfile)
!!$!deallocate(harms, cumdevs)
!!$!stop
!!$
!!$    ! Translate pos into filter parameters
!!$    do i = 1, ndi
!!$       if(data%sigma0(i) == 0) cycle
!!$       if (fix_highpass_freq_scan > 0.d0) then
!!$          data%filter_par(i,FILTER_HIGH_NU_SCAN) = fix_highpass_freq_scan
!!$       else
!!$          data%filter_par(i,FILTER_HIGH_NU_SCAN)= (pos(i,1)-0.5)
!!$       end if
!!$       data%filter_par(i,FILTER_HIGH_ALPHA) = -30
!!$       if (fix_lowpass_freq > 0.d0) then
!!$          data%filter_par(i,FILTER_LOW_NU)     = fix_lowpass_freq
!!$       else
!!$          data%filter_par(i,FILTER_LOW_NU)     = ind2freq(pos(i,2), data%samprate, nfreq) - 0.2d0
!!$       end if
!!$       data%filter_par(i,FILTER_LOW_ALPHA)  = -300
!!$    end do
!!$
!!$    deallocate(chisq_cumdev, chisq, N_fft, chisq_harm, pos, mask)
!!$    if(allocated(spikefilter)) deallocate(spikefilter)
!!$
!!$  end subroutine calc_bandpass_filter_par
!!$
!!$  ! This routine is currently not in use.
!!$  subroutine compute_fft_chisq(N_fft, powspec, chisq)
!!$    implicit none
!!$
!!$    real(dp),     dimension(0:), intent(in)  :: N_fft
!!$    real(sp),     dimension(0:), intent(in)  :: powspec
!!$    real(dp),                    intent(out) :: chisq
!!$
!!$    integer(i4b) :: i
!!$    real(dp)     :: n, w_sq, f
!!$
!!$    ! Compute effective chi-square and expected variance            
!!$    chisq = 0.d0
!!$    n     = 0.d0
!!$    do i = 0, size(powspec)-1
!!$       chisq = chisq + powspec(i) / N_fft(i)
!!$       n     = n     + 1.d0
!!$    end do
!!$    chisq = (chisq - n) / sqrt(n)
!!$  end subroutine compute_fft_chisq
!!$
!!$  ! Find which samples of each module hit something nasty
!!$  subroutine calc_mask(data, osys, mask)
!!$    implicit none
!!$    type(lx_struct)                            :: data
!!$    integer(i4b)                               :: osys, i, j, k, m, n, nmod, cmod, red
!!$    real(dp)                                   :: frad, brad
!!$    logical(lgt), dimension(:,:),  allocatable :: mask
!!$    integer(i4b), dimension(:),    allocatable :: strong, objs, uobjs, tmp, relevant
!!$    real(dp),     dimension(:,:,:),allocatable :: pos
!!$    n    = size(data%time)
!!$    nmod = get_num_modules()
!!$    cmod = get_center_module()
!!$    frad = get_focalplane_rad()
!!$    red  = 100
!!$    brad = get_max_fwhm()*fwhm2sigma*10;
!!$    allocate(mask(n,0:nmod-1))
!!$    mask = .true.
!!$    call wherei(patches%priority >= patch_strong_lim, strong)
!!$    allocate(pos(3,(n-1)/red+1,size(strong)),objs(n))
!!$    ! Which patches hit anything at all?
!!$    call get_patch_pos_multi(patches(strong), data%time(1:n:red), osys, pos)
!!$    call get_patch_hits(patches(strong), pos, real(data%point(1,:,cmod),dp), &
!!$     & real(data%point(2,:,cmod),dp), frad + brad, objs)
!!$    call uniqi(objs, uobjs)
!!$    call wherei(uobjs /= 0, tmp)
!!$    allocate(relevant(size(tmp)))
!!$    relevant = strong(uobjs(tmp))
!!$    ! Ok, we now have the list of patches that actually matter
!!$    do mod = 0, nmod-1
!!$       call get_patch_hits(patches(relevant), pos(:,:,uobjs(tmp)), &
!!$        & real(data%point(1,:,mod),dp), real(data%point(2,:,mod),dp), brad, objs)
!!$       where(objs /= 0) mask(:,mod) = .false.
!!$    end do
!!$    deallocate(uobjs, tmp, objs, strong, pos, relevant)
!!$  end subroutine
!!$
!!$  subroutine calc_fourier(data, ffts, powspecs)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    integer(i4b)    :: ndi, i, j, k, nf
!!$    complex(spc), dimension(:,:), allocatable :: ffts
!!$    real(sp),     dimension(:,:), allocatable :: powspecs
!!$    ndi   = size(data%tod,2)
!!$    nf    = size(data%tod,1)/2+1
!!$    allocate(ffts(nf,ndi), powspecs(nf, ndi))
!!$    call fft_multi(data%tod, ffts, 1)
!!$    do i = 1, ndi
!!$       call extract_powspec(ffts(:,i), powspecs(:,i))
!!$    end do
!!$  end subroutine
!!$
!!$  function ptrans_rot(phi1, theta1, phi2, theta2) result(angle)
!!$    implicit none
!!$    real(dp)           :: phi1, phi2, theta1, theta2, c2a, s2a, r1(3), r2(3), angle
!!$    call ang2vec(theta1,phi1,r1)
!!$    call ang2vec(theta2,phi2,r2)
!!$    call qu_transport_rot(r1,r2,c2a,s2a)
!!$    angle = atan2(s2a,c2a)/2
!!$  end function
!!$
!!$  subroutine calc_average_pointing(data, avgpt)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    real(dp)        :: avgpt(:), v(3), va(3), pt(2), foo
!!$    integer(i4b)    :: i
!!$    va  = 0
!!$    foo = 0
!!$    do i = 1, size(data%orig_point,2)
!!$       pt = data%orig_point(:,i)
!!$       call swap_coordinate_convention(pt(1), pt(2), foo)
!!$       call ang2vec(pt(2), pt(1), v)
!!$       va = va + v
!!$    end do
!!$    call vec2ang(va, avgpt(2), avgpt(1))
!!$    call swap_coordinate_convention(avgpt(1), avgpt(2), foo)
!!$    avgpt(3) = average_ang(real(data%orig_point(3,:),dp))
!!$  end subroutine

end program l3gen
