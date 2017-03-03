module ces_validate_mod
  use healpix_types
  use pix_tools
  use quiet_utils
  use quiet_fft_mod
  use quiet_ces_mod
  use quiet_Lx_mod
  use math_tools
  use quiet_module_mod
  use quiet_acceptlist_mod
  use quiet_mpi_mod
  implicit none

  integer(i4b), parameter :: st     = sp
  integer(i4b), parameter :: mpi_st = mpi_real

  type validation_struct
     integer(i8b) :: status
     real(dp)     :: full_fft_chisq, scanfreq_fft_chisq, tod_chisq
     real(dp)     :: sigma0, alpha, fknee
     real(dp)     :: azorder, nu_low, nu_high
     real(dp)     :: map_absmax, tod_absmax, tod_chi_az
     real(dp)     :: ces_accept_ratio
     real(dp)     :: typeb, wthr1, wthr2, jump
     real(dp), dimension(6)   :: band_fft_chisq
  end type validation_struct

  real(dp)     :: fft_chisq_threshold, map_absmax_threshold, map_min_mean_chisq, map_max_mean_chisq
  real(dp)     :: tod_chisq_threshold, tod_absmax_threshold, tod_az_max_chisq, tod_az_binsize
  real(dp)     :: min_diode_accept_ratio, fft_chisq_1_2Hz_threshold
  real(dp)     :: max_noise_alpha, min_noise_alpha, fbias_threshold
  logical(lgt) :: validate_highfreq_chisq
  logical(lgt), private :: apply_ces_diode_cuts
  real(dp)     :: fft_chisq_high_threshold, fft_chisq_scan_threshold
  real(dp)     :: fft_chisq_outlier_threshold, fft_chisq_low_threshold
  real(dp)     :: weather_threshold, weather_ces_fraction, typeb_threshold, fknee_abs_threshold_T, fknee_abs_threshold_P
  real(dp)     :: apex_max_pwv, sigma0_threshold, fknee_threshold, jump_finder_fractional_jump_threshold
  real(dp)     :: nu_high_threshold, nu_low_threshold
  real(dp)     :: fft_spike_min_freq, fft_spike_max_freq, fft_spike_threshold, fft_spike_dnu
  real(dp), allocatable, dimension(:,:,:), private :: cut_values

  logical(lgt), private :: initialized = .false.

  ! Status bit descriptors
  integer(i4b),       parameter             :: max_status = 64
  integer(i8b),       parameter             :: status_ok  = 0
  integer(i4b)                              :: num_status = 0
  integer(i4b)                              :: mask_danger  = 0, mask_bright = 0
  integer(i8b),       dimension(max_status) :: status_codes
  character(len=256), dimension(max_status) :: status_desc

  integer(i8b) :: status_baddi, status_nofile, status_static, status_elementary, &
   & status_corr, status_apex, status_typeb, status_weather, status_sigma0,      &
   & status_jump, status_tod, status_fft, status_sun, status_acceptlim, status_gain, &
   & status_fknee, status_filter, status_fbias
  integer(i8b) :: status_tod_chisq, status_tod_absmax, status_tod_az
  integer(i8b) :: status_fft_sfreq, status_fft_ulow, status_fft_full, &
   & status_fft_uhigh, status_fft_1_2, status_fft_spike, status_fft_1_0

contains

  subroutine initialize_ces_validate_mod(parfile)
    implicit none

    character(len=*),                   intent(in) :: parfile

    integer(i4b) :: unit, ndi, nces, i
    logical(lgt) :: havevalidate_highfreq_chisq, havefft_chisq_high_threshold

    if(initialized) return

    unit = getlun()

    ! Read parameters
    call get_parameter(unit, parfile,  'APPLY_CES_DIODE_CUTS',       par_lgt=apply_ces_diode_cuts)
    call get_parameter(unit, parfile,  'WEATHER_CHISQ_THRESHOLD',    par_dp=weather_threshold)
    call get_parameter(unit, parfile,  'WEATHER_CES_FRACTION',       par_dp=weather_ces_fraction)
    !call get_parameter(unit, parfile,  'JUMP_FINDER_CONSECUTIVE_SLOPE_THRESHOLD',    par_dp=jump_finder_consecutive_slope_threshold)
    !call get_parameter(unit, parfile,  'JUMP_FINDER_MAX_SLOPE_THRESHOLD',       par_dp=jump_finder_max_slope_threshold)
    call get_parameter(unit, parfile,  'JUMP_FINDER_FRACTIONAL_JUMP_THRESHOLD',       par_dp=jump_finder_fractional_jump_threshold)
    call get_parameter(unit, parfile,  'SIGMA0_THRESHOLD',           par_dp=sigma0_threshold)
    call get_parameter(unit, parfile,  'FKNEE_THRESHOLD',            par_dp=fknee_threshold)
    call get_parameter(unit, parfile,  'FKNEE_ABS_THRESHOLD_T',      par_dp=fknee_abs_threshold_T)
    call get_parameter(unit, parfile,  'FKNEE_ABS_THRESHOLD_P',      par_dp=fknee_abs_threshold_P)
    call get_parameter(unit, parfile,  'NU_LOWPASS_THRESHOLD',       par_dp=nu_low_threshold)
    call get_parameter(unit, parfile,  'NU_HIGHPASS_SCAN_THRESHOLD', par_dp=nu_high_threshold)
    call get_parameter(unit, parfile,  'TYPE_B_CHISQ_THRESHOLD',     par_dp=typeb_threshold)
    call get_parameter(unit, parfile,  'FFT_CHISQ_THRESHOLD',        par_dp=fft_chisq_threshold)
    call get_parameter(unit, parfile,  'FFT_CHISQ_SCAN_THRESHOLD',   par_dp=fft_chisq_scan_threshold)
    call get_parameter(unit, parfile,  'FFT_CHISQ_1.2HZ_THRESHOLD',  par_dp=fft_chisq_1_2Hz_threshold)
    call get_parameter(unit, parfile,  'FFT_CHISQ_HIGH_THRESHOLD',   par_dp=fft_chisq_high_threshold)
    call get_parameter(unit, parfile,  'FFT_CHISQ_LOW_THRESHOLD',    par_dp=fft_chisq_low_threshold)
    call get_parameter(unit, parfile,  'TOD_CHISQ_THRESHOLD',        par_dp=tod_chisq_threshold)
    call get_parameter(unit, parfile,  'TOD_ABSMAX_THRESHOLD',       par_dp=tod_absmax_threshold)
    call get_parameter(unit, parfile,  'TOD_AZ_BINSIZE',             par_dp=tod_az_binsize)
    call get_parameter(unit, parfile,  'TOD_AZ_MAX_CHISQ',           par_dp=tod_az_max_chisq)
    call get_parameter(unit, parfile,  'MAP_MIN_MEAN_CHISQ',         par_dp=map_min_mean_chisq)
    call get_parameter(unit, parfile,  'MAP_MAX_MEAN_CHISQ',         par_dp=map_max_mean_chisq)
    call get_parameter(unit, parfile,  'MAP_ABSMAX_THRESHOLD',       par_dp=map_absmax_threshold)
    call get_parameter(unit, parfile,  'CES_MIN_DIODE_ACCEPT_RATIO', par_dp=min_diode_accept_ratio)
    call get_parameter(unit, parfile,  'MAX_NOISE_ALPHA',            par_dp=max_noise_alpha)
    call get_parameter(unit, parfile,  'MIN_NOISE_ALPHA',            par_dp=min_noise_alpha)
    call get_parameter(unit, parfile,  'APEX_MAX_PWV',               par_dp=apex_max_pwv)
    call get_parameter(unit, parfile,  'FFT_CHISQ_SPIKE_MIN_FREQ',   par_dp=fft_spike_min_freq)
    call get_parameter(unit, parfile,  'FFT_CHISQ_SPIKE_MAX_FREQ',   par_dp=fft_spike_max_freq)
    call get_parameter(unit, parfile,  'FFT_CHISQ_SPIKE_DELTA_FREQ', par_dp=fft_spike_dnu)
    call get_parameter(unit, parfile,  'FFT_CHISQ_SPIKE_THRESHOLD',  par_dp=fft_spike_threshold)
    call get_parameter(unit, parfile,  'ONEOVERF_BIAS_THRESHOLD',    par_dp=fbias_threshold)

    status_baddi       = add_status("removing dead diodes")
    status_nofile      = add_status("removing bad files  ")
    status_static      = add_status("static cuts         ")
    status_gain        = add_status("gain cut            ")
    status_elementary  = add_status("elementary cuts     ")
    status_corr        = add_status("correlation cuts    ")
    status_apex        = add_status("APEX cuts           ", danger=.true.)
    status_typeb       = add_status("Type-B cuts         ")
    status_weather     = add_status("weather cuts        ", danger=.true.)
    status_sigma0      = add_status("sigma0 cut          ")
    status_fknee       = add_status("fknee cut           ")
    status_filter      = add_status("filter cut          ")
    status_jump        = add_status("jump cuts           ", danger=.true., bright=.true.)
    status_tod_chisq   = add_status("tod chisq cut       ", danger=.true., bright=.true.)
    status_tod_absmax  = add_status("tod absmax cut      ", danger=.true., bright=.true.)
    status_tod_az      = add_status("tod azimuth cut     ", danger=.true., bright=.true.)
    !status_tod         = add_status("tod cut             ", danger=.true., bright=.true.)
    !status_tod_chisq   = status_tod
    !status_tod_absmax  = status_tod
    !status_tod_az      = status_tod
    status_fft_sfreq   = add_status("scan freq chisq     ", danger=.true., bright=.true.)
    status_fft_ulow    = add_status("low  freq chisq     ", danger=.true., bright=.true.)
    status_fft_full    = add_status("all  freq chisq     ", danger=.true., bright=.true.)
    status_fft_uhigh   = add_status("high freq chisq     ", danger=.true., bright=.true.)
    status_fft_1_2     = add_status("1.2 Hz spike        ", danger=.true., bright=.true.)
    status_fft_1_0     = add_status("1.0 Hz spike        ", danger=.true., bright=.true.)
    status_fft_spike   = add_status("general spike       ", danger=.true., bright=.true.)
    !status_fft         = add_status("fft cut             ", danger=.true., bright=.true.)
    !status_fft_sfreq   = status_fft
    !status_fft_ulow    = status_fft
    !status_fft_full    = status_fft
    !status_fft_uhigh   = status_fft
    !status_fft_1_2     = status_fft
    !status_fft_1_0     = status_fft
    !status_fft_spike   = status_fft
    status_sun         = add_status("sun cut             ", danger=.true.)
    status_acceptlim   = add_status("accept ratio cut    ")
    status_fbias       = add_status("1/f sigma0 bias")

    initialized = .true.
  end subroutine initialize_ces_validate_mod

  subroutine initialize_diode_ces_cuts(diode_stats)
    implicit none
    real(st),         dimension(:,:,:), intent(in) :: diode_stats
    integer(i4b) :: i, ndi

    ! Set up diode specific cut values
    ndi = size(diode_stats(:,1,1))
    allocate(cut_values(ndi,DIODE_STAT_MAX,2))
    do i = 1, ndi
       call compute_conf_limit(real(diode_stats(i,DIODE_STAT_WEATHER1,:),dp), &
        & weather_threshold, 'two_sided', cut_values(i,DIODE_STAT_WEATHER1,:))
       call compute_conf_limit(real(diode_stats(i,DIODE_STAT_WEATHER2,:),dp), &
        & weather_threshold, 'two_sided', cut_values(i,DIODE_STAT_WEATHER2,:))
       call compute_conf_limit(real(log(diode_stats(i,DIODE_STAT_SIGMA0,:)),dp),&
        & sigma0_threshold, 'two_sided', cut_values(i,DIODE_STAT_SIGMA0,:))
       cut_values(i,DIODE_STAT_SIGMA0,:) = exp(cut_values(i,DIODE_STAT_SIGMA0,:))
       call compute_conf_limit(real(diode_stats(i,DIODE_STAT_FKNEE,:),dp),&
        & fknee_threshold, 'upper', cut_values(i,DIODE_STAT_FKNEE,:))
    end do

  end subroutine initialize_diode_ces_cuts


  ! Driver routines
  subroutine validate_tod(polarization, sigma_sq, az, tod, stat)
    implicit none

    logical(lgt),                intent(in)    :: polarization
    real(dp),                    intent(in)    :: sigma_sq
    real(sp),     dimension(1:), intent(in)    :: tod, az
    type(validation_struct),     intent(inout) :: stat

    integer(i4b) :: i
    logical(lgt) :: status

    ! Note: We thin the TOD by a given factor to avoid correlations between samples.
    !       Such correlations make it difficult to estimate the effective chi-square

    ! Check TOD chisquare
    status       = validate_tod_chisq(tod(::10), sigma_sq, stat%tod_chisq)
    if(.not. status) stat%status = ior(stat%status, status_tod_chisq)
    
    ! Check for strong outliers
    status       = validate_tod_absmax(tod, sigma_sq, stat%tod_absmax)
    if(.not. status) stat%status = ior(stat%status, status_tod_absmax)
    
    ! Check for azimuth structure
    status       = validate_tod_chi_az(tod(::10), az(::10), stat%tod_chi_az)
    if(.not. status) stat%status = ior(stat%status, status_tod_az)
  end subroutine

  subroutine validate_fft(polarization, samprate, scanfreq, N_filter, N_fft, powspec, stat)
    implicit none

    logical(lgt),                intent(in)    :: polarization
    real(dp),                    intent(in)    :: samprate, scanfreq
    real(dp),     dimension(0:), intent(in)    :: N_filter, N_fft
    real(sp),     dimension(0:), intent(in)    :: powspec
    type(validation_struct),     intent(inout) :: stat

    real(dp)     :: dnu, dnu_spike, c
    integer(i4b) :: i, ind_min, ind_max, ind, nu
    logical(lgt) :: status
    real(dp), allocatable, dimension(:) :: mask

    ! Check Fourier-space chi-square
    dnu          = ind2freq(2, samprate, size(powspec))


    ! Check full-range spectrum; mostly for verifying the noise model; check, but don't cut
    ind_min = 1
    ind_max = size(powspec)-1
    status       = validate_fft_chisq(.false., polarization, N_filter(ind_min:ind_max), &
         & N_fft(ind_min:ind_max), powspec(ind_min:ind_max), stat%full_fft_chisq)

    ! Check scan frequency; include a 10 mHz band around the scan frequency
    ind_min = nint((scanfreq-0.001d0) / dnu)
    ind_max = nint((scanfreq+0.001d0) / dnu)
    status       = validate_fft_chisq(.false., polarization, N_filter(ind_min:ind_max), N_fft(ind_min:ind_max), &
         & powspec(ind_min:ind_max), stat%scanfreq_fft_chisq, chi2_thresh=fft_chisq_scan_threshold)
    if(.not. status) stat%status = ior(stat%status, status_fft_sfreq)

    ! Unfiltered low frequencies; 0 to 0.2 Hz
    ind_min = 1
    ind_max = nint(0.2d0 / dnu)
    status       = validate_fft_chisq(.false., polarization, N_filter(ind_min:ind_max), N_fft(ind_min:ind_max), &
         & powspec(ind_min:ind_max), stat%band_fft_chisq(1), chi2_thresh=fft_chisq_low_threshold)
    if(.not. status) stat%status = ior(stat%status, status_fft_ulow)

    ! CMB data; filtered full range
    ind_min = 1
    ind_max = size(powspec)-1
    status       = validate_fft_chisq(.true., polarization, N_filter(ind_min:ind_max), N_fft(ind_min:ind_max), &
         & powspec(ind_min:ind_max), stat%band_fft_chisq(2))
    if(.not. status) stat%status = ior(stat%status, status_fft_full)

    ! High frequencies; 10 to 25 Hz; check, but don't cut necessarily. Pass in an explicit Chi^2 threshold
    ind_min = nint(10.d0 / dnu)
    ind_max = size(powspec)-1
    status       = validate_fft_chisq(.false., polarization, N_filter(ind_min:ind_max), N_fft(ind_min:ind_max), &
         & powspec(ind_min:ind_max), stat%band_fft_chisq(3),chi2_thresh=fft_chisq_high_threshold)
    if(.not. status) stat%status = ior(stat%status, status_fft_uhigh)

    ! Check for 1.2Hz spikes and harmonics
    allocate(mask(0:size(powspec)-1))
    mask = 0.d0
    do nu = 1, 7
       ind = nint(1.2d0*nu / dnu)
       mask(ind-2:min(ind+2,size(powspec)-1)) = 1.d0
    end do
    ind_min = 1
    ind_max = size(powspec)-1
    status = validate_fft_chisq(.true., polarization, mask(ind_min:ind_max), N_fft(ind_min:ind_max), &
         & powspec(ind_min:ind_max), stat%band_fft_chisq(4), chi2_thresh=fft_chisq_1_2Hz_threshold)
    deallocate(mask)
    if(.not. status) stat%status = ior(stat%status, status_fft_1_2)


    ! Check for FFT spikes 
    dnu_spike = nint(fft_spike_dnu      / dnu)
    ind_min   = nint(fft_spike_min_freq / dnu)
    ind_max   = ind_min
    stat%band_fft_chisq(5) = 0.d0
    do while (ind_max < nint(fft_spike_max_freq / dnu))
       ind_max   = ind_min + dnu_spike
       status       = validate_fft_chisq(.true., polarization, N_filter(ind_min:ind_max), N_fft(ind_min:ind_max), &
            & powspec(ind_min:ind_max), c, chi2_thresh=fft_spike_threshold)
       stat%band_fft_chisq(5) = max(stat%band_fft_chisq(5), c)
       if(.not. status) stat%status = ior(stat%status, status_fft_spike)
       ind_min   = ind_max
    end do

    ! Check for 1Hz spikes and harmonics, including 5Hz
    allocate(mask(0:size(powspec)-1))
    mask = 0.d0
    do nu = 1, 7
       ind = nint(nu / dnu)
       mask(ind-1:min(ind+1,size(powspec)-1)) = 1.d0
    end do
    ind_min = 1
    ind_max = size(powspec)-1
    status = validate_fft_chisq(.true., polarization, mask(ind_min:ind_max), N_fft(ind_min:ind_max), &
         & powspec(ind_min:ind_max), stat%band_fft_chisq(6), chi2_thresh=fft_chisq_1_2Hz_threshold)
    deallocate(mask)
    if(.not. status) stat%status = ior(stat%status, status_fft_1_0)
  end subroutine validate_fft

  subroutine validate_type_B_chisq(chisq, flag, stats)
    implicit none

    real(st), dimension(:,:),                intent(in)    :: chisq
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats

    integer(i4b) :: i, j, ndi, nces

    ndi  = size(chisq,1)
    nces = size(chisq,2)

    do i = 1, nces
       do j = 1, ndi
          if (chisq(j,i) > typeb_threshold) then
             stats(j,i)%status = ior(stats(j,i)%status, flag)
          end if
       end do
    end do

  end subroutine validate_type_B_chisq

  subroutine validate_weather(vals, weather_type, flag, stats)
    implicit none

    real(st), dimension(:,:),                intent(in)    :: vals
    integer(i4b),                            intent(in)    :: weather_type
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats

    integer(i4b) :: i, j, ndi, nces, num_accept, mod, di, num_reject
    real(dp)     :: chisq

    ndi  = size(vals,1)
    nces = size(vals,2)

    do i = 1, nces
       num_accept = 0
       num_reject = 0
       do j = 1, ndi
          if (.not. is_active(stats(j,i)%status)) cycle
          if (vals(j,i) >= cut_values(j,weather_type,1) .and. vals(j,i) <= cut_values(j,weather_type,2)) then
             num_accept = num_accept + 1
          else
             num_reject = num_reject + 1
             stats(j,i)%status = ior(stats(j,i)%status, flag)
          end if
       end do

       if (num_accept < weather_ces_fraction * (num_accept+num_reject)) then
          ! Reject all diodes due to bad weather overall
          do j = 1, ndi
             stats(j,i)%status = ior(stats(j,i)%status, flag)
          end do
       end if

    end do

  end subroutine validate_weather

  subroutine validate_jump(vals, flag, stats)
    implicit none

    real(dp), dimension(:,:),                intent(in)    :: vals
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats

    integer(i4b) :: i, j, ndi, nces, num_accept, mod, di, num_reject

    ndi  = size(vals,1)
    nces = size(vals,2)
    do i = 1, nces
       do j = 1, ndi
          if (.not. is_active(stats(j,i)%status)) cycle
          print *, "The diode is ", i, " ", j, " and it's stat value is ", vals(j,i), &
               & " and jump_finder fractional is ", jump_finder_fractional_jump_threshold
          if (vals(j,i).gt.jump_finder_fractional_jump_threshold) then
             stats(j,i)%status = ior(stats(j,i)%status, flag)
             print *, "This will be flagged"
          end if
       end do
     end do

  end subroutine validate_jump

  subroutine validate_filter_params(nu_low, nu_high, scanf, flag, stats)
    implicit none
    
    real(sp), dimension(:,:),                intent(in)    :: nu_low, nu_high
    real(dp), dimension(:),                  intent(in)    :: scanf
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats
    
    integer(i4b) :: i, j, ndi, nces, mod, di
    
    ndi  = size(nu_low,1)
    nces = size(nu_low,2)
    do i = 1, nces
       do j = 1, ndi
          if (.not. is_active(stats(j,i)%status)) cycle
          if (nu_low(j,i) < 0 .or. nu_high(j,i) < 0 .or. nu_low(j,i) < nu_low_threshold .or. nu_high(j,i) > nu_high_threshold) then
             stats(j,i)%status = ior(stats(j,i)%status, flag)
          end if

          ! TMR: Added interval validation - almost replaces the absolute thresholds, but not quite
          ! Note hardcoded value! Consider making another paramete.
          if (nu_low(j,i) - nu_high(j,i)*scanf(i) < 2.0) then
             stats(j,i)%status = ior(stats(j,i)%status, flag)
          end if

       end do
    end do
    
  end subroutine validate_filter_params

  subroutine validate_sigma0(vals, flag, stats)
    implicit none

    real(st), dimension(:,:),                intent(in)    :: vals
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats

    integer(i4b) :: i, j, ndi, nces, mod, di
    real(dp)     :: chisq

    ndi  = size(vals,1)
    nces = size(vals,2)
    do i = 1, nces
       do j = 1, ndi
          if (.not. is_active(stats(j,i)%status)) cycle
          if (vals(j,i) < cut_values(j,DIODE_STAT_SIGMA0,1) .or. vals(j,i) > cut_values(j,DIODE_STAT_SIGMA0,2)) then
             stats(j,i)%status = ior(stats(j,i)%status, flag)
          end if
       end do
     end do

  end subroutine validate_sigma0

  subroutine validate_fknee(vals, flag, stats)
    implicit none

    real(st), dimension(:,:),                intent(in)    :: vals
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats

    integer(i4b) :: i, j, ndi, nces, mod, di, ndiode
    real(dp)     :: chisq

    ndiode = get_num_diodes()
    ndi  = size(vals,1)
    nces = size(vals,2)
    do i = 1, nces
       do j = 1, ndi
          if (.not. is_active(stats(j,i)%status)) cycle
          if (vals(j,i) > cut_values(j,DIODE_STAT_FKNEE,2))                                 &
               & stats(j,i)%status = ior(stats(j,i)%status, flag)
          if (is_polarization_module((j-1)/ndiode) .and. vals(j,i) > fknee_abs_threshold_P) &
               & stats(j,i)%status = ior(stats(j,i)%status, flag)
          if (is_temperature_module((j-1)/ndiode)  .and. vals(j,i) > fknee_abs_threshold_T) &
               & stats(j,i)%status = ior(stats(j,i)%status, flag)
       end do
     end do
  end subroutine validate_fknee

  subroutine validate_fbias(vals, flag, stats)
    implicit none
    real(st), dimension(:,:,:),  intent(in)    :: vals
    integer(i8b),                intent(in)    :: flag
    type(validation_struct), dimension(:,:)    :: stats
    integer(i4b) :: i, j, ndi, nces, mod, di
    real(dp)     :: bias
    ndi  = size(vals,1)
    nces = size(vals,3)
    do i = 1, nces
       do j = 1, ndi
          if (.not. is_active(stats(j,i)%status)) cycle
          bias = (1/vals(j,DIODE_STAT_FKNEE,i))**vals(j,DIODE_STAT_ALPHA,i)
          if(bias > fbias_threshold) stats(j,i)%status = ior(stats(j,i)%status, flag)
       end do
    end do
   end subroutine

  ! Remove things that are obviously useless, without need of parameters.
  ! Things that go here typically break assumptions of tod2map.
  subroutine validate_elementary(dstats, flag, stats)
    implicit none
    real(st), dimension(:,:,:),              intent(in)    :: dstats
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats
    integer(i4b) :: i
!    where(dstats(:,DIODE_STAT_GAIN  ,:) ==   0) stats%status = ior(stats%status, flag)
    where(dstats(:,DIODE_STAT_SIGMA0,:) <=   0) stats%status = ior(stats%status, flag)
    ! 100 is the default value, which means that the fit failed.
    where(dstats(:,DIODE_STAT_FKNEE,:)  >= 100) stats%status = ior(stats%status, flag)
    where(dstats(:,DIODE_STAT_FKNEE,:)  <=   0) stats%status = ior(stats%status, flag)

! TMR debug trial
!    where(dstats(:,DIODE_STAT_ALPHA,:)  < -7.0 ) stats%status = ior(stats%status, flag)
  end subroutine

  ! Remove things that are obviously useless, without need of parameters.
  ! Things that go here typically break assumptions of tod2map.
  subroutine validate_gain(dstats, flag, stats)
    implicit none
    real(st), dimension(:,:,:),              intent(in)    :: dstats
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats
    integer(i4b) :: i
    where(dstats(:,DIODE_STAT_GAIN  ,:) ==   0) stats%status = ior(stats%status, flag)
  end subroutine

  subroutine validate_sun(vals, flag, stats)
    implicit none
    real(st), dimension(:,:),                intent(in)    :: vals
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:)                :: stats
    where(vals > 0) stats%status = ior(stats%status, flag)
  end subroutine validate_sun

  subroutine validate_static(flag, stats)
    implicit none
    integer(i8b),                           intent(in)    :: flag
    type(validation_struct),dimension(:,:), intent(inout) :: stats
    integer(i4b)                                          :: i

    ! Runs in July 2009 use different, poorly calibrated hardware.
    do i = 1, size(ces_db%ceses)
       if(ces_db%ceses(i)%mjd(2) < 55044 .and. ces_db%ceses(i)%mjd(1) > 55013) &
        & stats(:,i)%status = ior(stats(:,i)%status, flag)

    ! Q-band observations before Oct 24 2008 were calibration/testing only
       if(ces_db%ceses(i)%mjd(2) < 54763) stats(:,i)%status = ior(stats(:,i)%status, flag)
    end do
    
  end subroutine

  !subroutine validate_corr(corr, flag, stats)
  !  implicit none
  !  real(dp),                intent(in)    :: corr(:,:)
  !  integer(i8b),            intent(in)    :: flag
  !  type(validation_struct), intent(inout) :: stats(:)
  !  integer(i4b)                           :: status, i
  !  integer(i4b), allocatable              :: inds(:)
  !  real(dp),     allocatable              :: mat(:,:)
  !  ! Find a valid diode
  !  do i = 1, size(corr,2)
  !     if(any(corr(:i-1,i) == corr(:i-1,i)) .or. any(corr(i+1:,i) == corr(i+1:,i))) exit
  !  end do
  !  ! Use this to find all valid diodes
  !  call wherei(corr(:,i) == corr(:,i), inds)
  !  allocate(mat(size(inds),size(inds)))
  !  mat = corr(inds,inds)
  !  call invert_matrix(mat, cholesky=.true., status=status)
  !  if(status /= 0) stats%status = ior(stats%status, flag)
  !  deallocate(inds, mat)
  !end subroutine

  subroutine validate_APEX(ces_stats, flag, stats)
    implicit none

    real(st),                dimension(:,:), intent(in)    :: ces_stats
    integer(i8b),                            intent(in)    :: flag
    type(validation_struct), dimension(:,:), intent(inout) :: stats    

    integer(i4b) :: i, j, ndi, nces, mod, di
    type(quiet_ces_info) :: ces

    ndi  = size(stats,1)
    nces = size(stats,2)

    do i = 1, nces
       call get_ces_info(i, ces)
       ! Check PWV
       if (ces_stats(STAT_PWV,i) > apex_max_pwv) then
          do j = 1, ndi
             if (is_active(stats(j,i)%status)) stats(j,i)%status = ior(stats(j,i)%status, flag)
          end do
       end if
    end do

  end subroutine validate_APEX

  ! If more than X% of the diodes in this ces were cut for 'dangerous'
  ! reasons, cut them all.
  subroutine validate_global_ces(data, stats)
    implicit none
    type(Lx_struct),                       intent(in)    :: data
    type(validation_struct), dimension(:), intent(inout) :: stats
    integer(i4b) :: nbad
    real(dp)     :: ces_accept_ratio

    ! Check ratio of accepted diodes; if lower than threshold, reject all
    nbad             = count(iand(stats%status,mask_danger) /= 0)
    ces_accept_ratio = 1-real(nbad,dp)/size(stats)

    if(ces_accept_ratio < min_diode_accept_ratio) &
     & stats%status = ior(stats%status,status_acceptlim)
    stats%ces_accept_ratio = ces_accept_ratio
  end subroutine validate_global_ces

  ! Individual test routines
  function validate_tod_chisq(tod, sigma_sq, chisq_out)
    implicit none

    real(sp),     dimension(1:), intent(in)  :: tod
    real(dp),                    intent(in)  :: sigma_sq
    logical(lgt)                             :: validate_tod_chisq
    real(dp),                    intent(out) :: chisq_out
    
    integer(i4b) :: i
    real(dp)     :: chisq, w, w_sq, rms, mu, n

    ! Compute effective chi-square and expected variance
    n   = real(size(tod),dp)
!    rms = sqrt(sum(tod**2) / n)
    chisq = 0.d0
    do i = 1, size(tod)
       chisq = chisq + tod(i)**2 / sigma_sq
    end do

    chisq_out = (chisq - n) / sqrt(2.d0*n)
    validate_tod_chisq = abs(chisq_out) < tod_chisq_threshold

  end function validate_tod_chisq

  function validate_tod_absmax(tod, sigma_sq, absmax)
    implicit none

    real(sp),     dimension(1:), intent(in)  :: tod
    real(dp),                    intent(in)  :: sigma_sq
    logical(lgt)                             :: validate_tod_absmax
    real(dp),                    intent(out) :: absmax
    
    integer(i4b) :: i
    real(dp)     :: chisq, w, w_sq, rms, mu, n

    rms = sqrt(sigma_sq)
    absmax = maxval(abs(tod/rms))

    validate_tod_absmax = abs(absmax) < tod_absmax_threshold

  end function validate_tod_absmax

  function validate_tod_chi_az(tod, az, chi_az)
    implicit none

    real(sp),     dimension(1:), intent(in)  :: tod, az
    logical(lgt)                             :: validate_tod_chi_az
    real(dp),                    intent(out) :: chi_az

    integer(i4b) :: i, numbin
    real(dp)     :: sigma, chisq, sigma0, mu, n
    integer(i4b), allocatable, dimension(:) :: counts
    real(sp),     allocatable, dimension(:) :: binned, rms
    
    mu     = sum(tod) / size(tod)
    sigma0 = sqrt(sum((tod-mu)**2) / size(tod))

    numbin = int(360.d0 / tod_az_binsize)+1
    allocate(counts(numbin), binned(numbin), rms(numbin))
    call azbin_tod(tod, az, binned, counts)

    chisq = 0.d0
    n     = 0.d0
    do i = 1, numbin
       if (counts(i) > 0) then
          sigma = sigma0 / sqrt(real(counts(i),dp))
          chisq = chisq + binned(i)**2 / sigma**2
          n     = n + 1.d0
       end if
    end do

    chi_az = (chisq - n) / sqrt(2.d0*n)
    validate_tod_chi_az = abs(chi_az) < tod_az_max_chisq

    deallocate(counts, binned, rms)

  end function validate_tod_chi_az


  function validate_fft_chisq(apply_filter, polarization, N_filter, N_fft, powspec, chisq_out, chi2_thresh)
    implicit none

    logical(lgt),                intent(in)  :: apply_filter, polarization
    real(dp),     dimension(0:), intent(in)  :: N_filter, N_fft
    real(sp),     dimension(0:), intent(in)  :: powspec
    logical(lgt)                             :: validate_fft_chisq
    real(dp),                    intent(out) :: chisq_out
    real(dp),     optional,      intent(in)  :: chi2_thresh
    real(dp)                                 :: chi2_thresh_

    integer(i4b) :: i
    real(dp)     :: chisq, w, w_sq, f

    ! Compute effective chi-square and expected variance
    chisq = 0.d0
    w     = 0.d0
    w_sq  = 0.d0
    do i = 0, size(powspec)-1
       if (apply_filter) then
          f = N_filter(i)
       else
          f = 1.d0
       end if
       chisq = chisq + f * powspec(i) / N_fft(i)
       w     = w     + f
       w_sq  = w_sq  + f**2
    end do

    chisq_out = (chisq - w) / sqrt(w_sq)

    chi2_thresh_=fft_chisq_threshold
    if(present(chi2_thresh)) chi2_thresh_ = chi2_thresh
    validate_fft_chisq = abs(chisq_out) < chi2_thresh_

  end function validate_fft_chisq


  function validate_fft_outliers(polarization, samprate, N_filter, N_fft, powspec, outlier_stat)
    implicit none

    logical(lgt),                intent(in)  :: polarization
    real(dp),                    intent(in)  :: samprate
    real(dp),     dimension(0:), intent(in)  :: N_filter, N_fft, powspec
    logical(lgt)                             :: validate_fft_outliers
    real(dp),                    intent(out) :: outlier_stat

    integer(i4b) :: i, num_outlier, n, w, upper, lower, ind_min, ind_max
    real(dp)     :: nu_min, nu_max, nu_mean, nu_stddev, wn_level, threshold
    real(dp), allocatable, dimension(:) :: p_smooth

    nu_min  = 5.d0 ! Hz
    nu_max  = 8.d0 ! Hz
    w       = 50 ! Smoothing window
    ind_min = freq2ind(nu_min, samprate, size(powspec))
    ind_max = freq2ind(nu_max, samprate, size(powspec))

    if (ind_min < w .or. ind_max > size(powspec)-w) then
       write(*,*) 'Error in validate_fft_outliers: Too few points in spectrum'
    end if

    ! Smooth power spectrum
    allocate(p_smooth(ind_min:ind_max))
    do i = ind_min, ind_max
       p_smooth(i) = sum(powspec(i-w:i+w)) / real(2*w+1,dp)
    end do

    ! Estimate white noise level
    wn_level  = minval(p_smooth)
    threshold = fft_chisq_outlier_threshold * wn_level
    
    outlier_stat = maxval(p_smooth) / wn_level
    validate_fft_outliers = .not. any(p_smooth > threshold)

    deallocate(p_smooth)

  end function validate_fft_outliers

  function accepted_map_absmax(polarization, rms, map)
    implicit none

    logical(lgt),                intent(in) :: polarization
    real(dp),     dimension(1:), intent(in) :: rms, map
    logical(lgt)                            :: accepted_map_absmax

    real(dp) :: mu, n, r

    n  = real(size(map),dp)
    mu = sum(map) / n

    accepted_map_absmax = maxval(abs((map-mu)/rms)) < map_absmax_threshold

  end function accepted_map_absmax

  function accepted_map_chisq(polarization, rms, map)
    implicit none

    logical(lgt),                intent(in) :: polarization
    real(dp),     dimension(1:), intent(in) :: rms, map
    logical(lgt)                            :: accepted_map_chisq
    
    integer(i4b) :: i, n, m
    real(dp) :: mu, r

    n  = size(map)
    m  = 0
    mu = 0.d0
    do i = 1, n
       if (map(i) /= 0.d0) then
          mu = mu + map(i)
          m  = m  + 1
       end if
    end do
    mu = mu / real(m,dp)

    r = 0.d0
    do i = 1, n
       if (map(i) /= 0.d0) then
          r  = r + ((map(i)-mu)/rms(i))**2
       end if
    end do
    r = r / real(m,dp)

    write(*,*) 'map_mean_chisq = ', r
    accepted_map_chisq = r > map_min_mean_chisq .and. r < map_max_mean_chisq

  end function accepted_map_chisq

  subroutine initialize_validation_struct(stat)
    implicit none
    type(validation_struct), intent(out) :: stat
    stat%status = status_ok
    stat%full_fft_chisq     = 0
    stat%scanfreq_fft_chisq = 0
    stat%tod_chisq          = 0
    stat%tod_chi_az         = 0
    stat%tod_absmax         = 0
    stat%map_absmax         = 0
    stat%ces_accept_ratio   = 0
    stat%band_fft_chisq     = 0
    stat%sigma0             = 0
    stat%alpha              = 0
    stat%fknee              = 0
    stat%azorder            = 0
    stat%nu_high            = 0
    stat%nu_low             = 0
    stat%band_fft_chisq     = 0
    stat%jump               = 0
  end subroutine initialize_validation_struct

  subroutine initialize_validation_structs(stats)
    implicit none
    type(validation_struct), dimension(:,:), intent(out) :: stats
    integer(i4b) :: i, j
    do i = 1, size(stats,1)
       do j = 1, size(stats,2)
          call initialize_validation_struct(stats(i,j))
       end do
    end do
  end subroutine initialize_validation_structs

  subroutine update_current_stat(mod, di, stat, data, diode_stats)
    implicit none
    integer(i4b), intent(in) :: mod, di
    type(validation_struct)  :: stat
    type(Lx_struct)          :: data
    real(st), intent(in)     :: diode_stats(:,:)
    integer(i4b) :: diode
    diode = mod*get_num_diodes()+di+1
    stat%sigma0             = data%sigma0(diode)
    stat%alpha              = data%alpha(diode)
    stat%fknee              = data%fknee(diode)
    stat%azorder            = diode_stats(diode,DIODE_STAT_AZORDER)
    stat%nu_low             = diode_stats(diode,DIODE_STAT_NU_LOW)
    stat%nu_high            = diode_stats(diode,DIODE_STAT_NU_HIGH)
    stat%typeb              = diode_stats(diode,DIODE_STAT_TYPEB)
    stat%wthr1              = diode_stats(diode,DIODE_STAT_WEATHER1)
    stat%wthr2              = diode_stats(diode,DIODE_STAT_WEATHER2)
    stat%jump               = diode_stats(diode,DIODE_STAT_JUMP)
  end subroutine update_current_stat

  ! Find the average pointing over all diodes in assembly, by looking at every
  ! 100 samples.
  function average_position(pixels, nside) result(pos)
    implicit none
    integer(i4b) :: pixels(:), nside, n, i
    real(dp)     :: pos(2), avg(3), vec(3), theta, phi
    n   = size(pixels)
    avg = 0
    do i = 1, n, max(1, n/1000)
       call pix2vec_nest(nside, pixels(i), vec)
       avg = avg+vec
    end do
    avg = avg/sqrt(sum(avg**2))
    call vec2ang(avg, pos(1), pos(2))
  end function

  function angular_distance(pos1, pos2) result(r)
    implicit none
    real(dp) :: pos1(2), pos2(2), r, v1(3), v2(3)
    call ang2vec(pos1(1),pos1(2), v1)
    call ang2vec(pos2(1),pos2(2), v2)
    call angdist(v1, v2, r)
  end function

  subroutine status2cuts(status, cuts, level)
    implicit none

    integer(i4b),                intent(in)  :: status
    logical(lgt), dimension(1:), intent(out) :: cuts
    integer(i4b),                intent(out) :: level

    integer(i4b) :: i

    cuts  = .false.
    level = 0
    do i = 1, size(cuts)
       cuts(i) = btest(status,i-1)
       if (cuts(i) .and. level == 0) level = i
    end do

  end subroutine status2cuts

  function is_active(status)
    implicit none

    integer(i8b), intent(in) :: status
    logical(lgt)             :: is_active
    is_active = .not. btest(status,0) .and. .not. btest(status,1)
  end function is_active

  function add_status(desc, danger, bright) result(bit)
    implicit none
    character(len=*), intent(in)           :: desc
    logical(lgt),     intent(in), optional :: danger, bright
    integer(i8b)                           :: bit
    logical(lgt)                           :: danger_, bright_
    danger_ = .false.; if(present(danger)) danger_ = danger
    bright_ = .false.; if(present(bright)) bright_ = bright
    bit                      = ishft(int(1,i8b),num_status)
    num_status               = num_status + 1
    status_codes(num_status) = bit
    status_desc (num_status) = desc
    if(danger_) mask_danger = ior(mask_danger, bit)
    if(bright_) mask_bright = ior(mask_bright, bit)
  end function

end module ces_validate_mod
