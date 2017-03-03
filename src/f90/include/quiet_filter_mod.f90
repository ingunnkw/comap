module quiet_filter_mod
  use healpix_types
  use quiet_assembly_mod
  use quiet_fft_mod
  use math_tools
  use spline_1D_mod
  implicit none

  type filter_params
     logical(lgt)       :: apply_lowpass, apply_highpass, apply_spike
     logical(lgt)       :: ignore_oneoverf, ignore_diode_corr, ignore_diode_corr_freq
     logical(lgt)       :: use_precomp_filter, use_woodbury_spikes=.false.
     real(dp)           :: nu_highpass, alpha_highpass, nu_lowpass, alpha_lowpass
     real(dp), allocatable, dimension(:) :: spikefreqs
     integer(i4b)       :: num_spikes
  end type filter_params

  type filter_struct
     logical(lgt) :: apply_lowpass_filter, apply_highpass_filter, apply_spikefilter
     integer(i4b) :: numsamp, n
     real(dp)     :: fft_low_freq, fft_low_alpha
     real(dp)     :: fft_high_freq, fft_high_alpha
     real(dp)     :: sigma0, fknee, alpha
     real(dp),     pointer,     dimension(:) :: N_corr_real, N_corr_F_sq_real
     real(dp),     allocatable, dimension(:) :: N_corr_fft, N_corr_F_sq_fft
  end type filter_struct

  type noise_corr
     integer(i4b)                        :: n
     real(dp), allocatable, dimension(:) :: N_corr
  end type noise_corr

!TMR
  type noise_template
     logical(lgt)                              :: use_templates
     integer(i4b)                              :: ndi
     integer(i4b), allocatable, dimension(:)   :: nbins
     real(dp),     allocatable, dimension(:,:) :: nu, temps, derivs
  end type noise_template

  type(filter_params),  private :: default_filter
  type(noise_template), private :: templates    ! TMR
  real(dp),             private :: corr_func_cutoff = 1d-5
  integer(i4b),         private :: corrlen_cap = -1


contains

  subroutine initialize_filter_mod(paramfile)
    implicit none
    character(len=*), intent(in) :: paramfile
    call read_filters_parfile(paramfile, default_filter, templates)
    if (default_filter%apply_spike) call read_spikes(paramfile, default_filter)
    if (templates%use_templates)    call read_templates(paramfile, templates) ! TMR
  end subroutine initialize_filter_mod

  ! Reads list of frequencies to be cut by spikefilter
  subroutine read_spikes(paramfile,filter)
    implicit none
    type(filter_params), intent(inout) :: filter
    character(len=*)   , intent(in)    :: paramfile
    character(len=256)                 :: spikefile
    integer(i4b)                       :: unit, i, nspikes
    unit=getlun()
    call get_parameter(unit, paramfile, 'SPIKELIST', par_string=spikefile)    

    open(unit,file=trim(spikefile), action='read',status='old')
    read(unit,*) nspikes
    filter%num_spikes=nspikes
    allocate(filter%spikefreqs(nspikes))
    do i=1,nspikes
       read(unit,*) filter%spikefreqs(i)
    end do
    close(unit)
  end subroutine read_spikes

  ! Read weather template functions - one file per diode - and preparing for spline
  subroutine read_templates(paramfile,t)
    implicit none
    type(noise_template), intent(inout)   :: t
    character(len=*),     intent(in)      :: paramfile
    character(len=256)                    :: temp_dir, file
    integer(i4b)                          :: unit, i, di, ndi, maxbins, lastbin
    unit = getlun()
    call get_parameter(unit, paramfile, 'TEMPLATES_DIR', par_string=temp_dir)
    file = trim(temp_dir) // '/info.dat'
    open(unit, file=file, action='read', status='old')
    read(unit,*) ndi, maxbins
    close(unit)
    t%ndi = ndi
    allocate(t%nbins(ndi))
    allocate(t%temps(maxbins,ndi))
    allocate(t%nu(maxbins,ndi))
    allocate(t%derivs(maxbins,ndi))
    ! Temps are not of equal length (small differences)
    
    ! Reading template functions from file
    do di=1,ndi
       file = trim(temp_dir) // '/diode_' // trim(itoa(di,3)) // '.dat'
       open(unit, file=file, action='read', status='old')
       read(unit,*) t%nbins(di)
       do i=1,t%nbins(di)
          read(unit,*) t%nu(i,di), t%temps(i,di)
       end do
       close(unit)
    end do

    ! Doing first step of spline procedure (so we won't have to do the same stuff for each ces)
    do di = 1,ndi
       lastbin = t%nbins(di)
       call spline(t%nu(1:lastbin,di),t%temps(1:lastbin,di),1.d30,1.d30,t%derivs(1:lastbin,di))    
    end do

  end subroutine read_templates

! Retrieving template functions from noise_template structure, splined to correct length.
  subroutine get_templates(t_full, nf, srate, diodes_in)
    implicit none
    real(sp), dimension(:,:), allocatable, intent(out) :: t_full
    real(dp),                              intent(in)  :: srate
    integer(i4b),                          intent(in)  :: nf
    integer(i4b), optional, dimension(:),  intent(in)  :: diodes_in
    integer(i4b)                                       :: a, b, di, i, lastbin
    integer(i4b), dimension(:), allocatable            :: diodes

    ! Use full set of diodes (from local structure) unless otherwise specified
    if(present(diodes_in)) then
       allocate(diodes(size(diodes_in)))
       diodes = diodes_in
    else
       allocate(diodes(templates%ndi))
       do i=1,templates%ndi
          diodes(i) = i
       end do
    end if

    allocate(t_full(nf,size(diodes)))
    do di = 1,size(diodes)
       lastbin = templates%nbins(diodes(di))
       do i = 1,nf
          t_full(i,di) = splint(templates%nu(1:lastbin,diodes(di)),&
               & templates%temps(1:lastbin,diodes(di)),templates%derivs(1:lastbin,diodes(di)),&
               & ind2freq(i,srate,nf))
       end do
    end do

!!$open(42,file='splined_temps.dat',action='write')
!!$do i=1,nf
!!$   write(42,'(77g14.6)') ind2freq(i,srate,nf), t_full(i,:)
!!$end do
!!$close(42)

  end subroutine get_templates


  function get_inv_N_covar_matrix(get_F_sq, i, j, diode1, diode2, assembly)
    implicit none

    logical(lgt),         intent(in) :: get_F_sq
    integer(i4b),         intent(in) :: i, j, diode1, diode2
    type(quiet_assembly), intent(in) :: assembly
    real(dp)                         :: get_inv_N_covar_matrix

    integer(i4b) :: r
    real(dp)     :: res

    r = abs(i-j)
    if (r < assembly%noise_corr_length) then
       if (get_F_sq) then
          res = assembly%N_corr_F_sq_real(abs(i-j),diode1,diode2)
       else
          res = assembly%N_corr_real(abs(i-j),diode1,diode2)
       end if
    else
       res = 0.d0
    end if

    get_inv_N_covar_matrix = res

  end function get_inv_N_covar_matrix


  subroutine get_default_filter(filteropts)
    implicit none
    type(filter_params), intent(out) :: filteropts
    filteropts = default_filter
  end subroutine get_default_filter

  subroutine get_filter_info(cid, mod, diode, filter, filtopts)
    implicit none

    integer(i4b),            intent(in)  :: cid, mod, diode
    type(filter_struct),     intent(out) :: filter
    type(filter_params),     optional    :: filtopts

    integer(i4b)        :: status
    real(dp)            :: nu
    type(filter_params) :: fopts

    fopts = default_filter; if(present(filtopts)) fopts = filtopts

    if(fopts%ignore_oneoverf) then
       filter%fknee = 1e-6
       filter%alpha  = -10
    end if

    filter%apply_highpass_filter = fopts%apply_highpass
    filter%apply_lowpass_filter  = fopts%apply_lowpass
    filter%apply_spikefilter     = fopts%apply_spike
    
    ! Set up apodized highpass filter
    if (fopts%apply_highpass) then
       filter%fft_high_freq  = fopts%nu_highpass
       filter%fft_high_alpha = fopts%alpha_highpass
    end if

    ! Set up apodized lowpass filter
    if (fopts%apply_lowpass) then
       filter%fft_low_freq  = fopts%nu_lowpass
       filter%fft_low_alpha = fopts%alpha_lowpass
    end if

  end subroutine get_filter_info

  subroutine deallocate_filter(filter)
    implicit none

    type(filter_struct), intent(inout) :: filter

    if (associated(filter%N_corr_real))      deallocate(filter%N_corr_real)
    if (allocated(filter%N_corr_fft))        deallocate(filter%N_corr_fft)
    if (associated(filter%N_corr_F_sq_real)) deallocate(filter%N_corr_F_sq_real)
    if (allocated(filter%N_corr_F_sq_fft))   deallocate(filter%N_corr_F_sq_fft)

  end subroutine deallocate_filter


  ! Generate inverse power spectrum of 1/f noise
  subroutine get_inv_noise_filter_fft(inverse, numsamp, samprate, sigma_0, fknee, alpha, N_filter, highpass)
    implicit none

    logical(lgt),                         intent(in)  :: inverse, highpass
    integer(i4b),                         intent(in)  :: numsamp
    real(dp),                             intent(in)  :: samprate, sigma_0, fknee, alpha
    real(dp),     dimension(0:numsamp/2), intent(out) :: N_filter

    integer(i4b) :: i
    real(dp)     :: f

    do i = 1, size(N_filter)-1
       f = ind2freq(i+1, samprate, size(N_filter))
       if (inverse) then
          N_filter(i) = 1d0/(sigma_0**2*(1+(f/fknee)**alpha))
       else
          N_filter(i) = sigma_0**2*(1+(f/fknee)**alpha)
       end if
    end do
    N_filter(0) = N_filter(1)
    if(highpass) N_filter(0) = 0.d0

  end subroutine get_inv_noise_filter_fft

  ! Generate inverse power spectrum of a 1/f-type highpass/lowpass filter.
  subroutine apodize_filter_fft(inverse, numsamp, samprate, nu_cut, alpha, highpass, N_filter)
    implicit none

    integer(i4b),                         intent(in)  :: numsamp
    real(dp),                             intent(in)  :: samprate, nu_cut, alpha
    logical(lgt),                         intent(in)  :: inverse, highpass
    real(dp),     dimension(0:numsamp/2), intent(inout) :: N_filter

    real(dp)     :: f, a, fact
    integer(i4b) :: i

    a = alpha
    if (.not. highpass) a = -a
    N_filter(0) = 1.d0
    if (highpass) N_filter(0) = 0.d0
    ! Hm. This means this sub needs to be run first for lowpass, then for highpass, if both filters are on.

    do i = 1, numsamp/2
       f = ind2freq(i+1, samprate, size(N_filter))
       fact = (1+(f/nu_cut)**a)
       fact = max(fact, 1.d-100)
       fact = min(fact, 1.d100)
       if (inverse) then
          N_filter(i) = N_filter(i)/fact
       else
          N_filter(i) = N_filter(i)*fact
       end if
    end do

  end subroutine apodize_filter_fft

  ! Filter frequencies from a priveded function using the list of frequencies
  ! read into the default filter in initialize_filter_mod. 
  subroutine spike_filter_fft(inverse, spikes, nspikes, nsamp, srate, nf, fltfunc)
    implicit none
    logical(lgt)          , intent(in)    :: inverse
    integer(i4b)          , intent(in)    :: nspikes, nsamp, nf
    real(dp)              , intent(in)    :: srate, spikes(:)
    real(dp), dimension(:), intent(inout) :: fltfunc
    integer(i4b)                          :: unit, i, spikeind
    do i=1,nspikes
       spikeind = freq2ind(spikes(i), srate, nf)
       ! TMR: Testing if there is a difference between 1 and 3 freqs per spike
!       fltfunc(spikeind-1:spikeind+1) = 1.d+100 
       if (inverse) then
          fltfunc(spikeind) = 1.d+100 ! The filter is inverse
       else
          fltfunc(spikeind) = 1.d-100
       end if
    end do

  end subroutine spike_filter_fft

  ! Get the inverse time correlation function N_filter_real corresponding to N_filter_fft
  ! by using a fourier transform. The resulting length n is that needed to contain
  ! the parts of the output with non-negligible values.
  subroutine get_inv_noise_filter_real(numsamp, N_filter_fft, n, N_filter_real)
    implicit none

    integer(i4b),                                  intent(in)  :: numsamp
    real(dp),              dimension(0:numsamp/2), intent(in)  :: N_filter_fft
    integer(i4b),                                  intent(out) :: n
    real(dp),     allocatable, dimension(:),           intent(out) :: N_filter_real

    integer(i4b) :: i, n_tot
    real(dp)     :: area, cum, offset
    real(dp),     allocatable, dimension(:) :: d_real
    complex(dpc), allocatable, dimension(:) :: d_fft

    ! The null-filter
    if(all(N_filter_fft == 1)) then
       n = 0
       allocate(N_filter_real(-0:0))
       N_filter_real = 1
       return
    end if

    if (sum(abs(N_filter_fft)) == 0.d0) then
       n = 0
       allocate(N_filter_real(-numsamp+1:numsamp-1))
       N_filter_real = 0.d0
       return
    end if

    ! Go from power spectrum to time domain. Need to pad to go from real to
    ! 'complex'.
    n_tot = numsamp/2+1
    allocate(d_real(0:numsamp-1))
    allocate(d_fft(0:n_tot-1))
    d_fft = N_filter_fft
    call fft(d_real, d_fft, -1)

    ! Find area with noticable values. We do so by integrating from
    ! the start and out until we have 1-corr_func_cutoff of the
    ! total absolute area under the function.
    offset = d_real((numsamp-1)/2)
    area = sum(abs(d_real(0:(numsamp-1)/2)-offset))
    cum  = 0
    do i = 0, numsamp-1
       cum = cum + abs(d_real(i)-offset)
!write(*,'(i6,5e15.7)') i, d_real(i), d_real(i)-offset, cum, 1-cum/area, corr_func_cutoff
       if(1-cum/area < corr_func_cutoff) exit
    end do
    n = min(i,numsamp-1)

    ! Rescale to get the correlation function
    d_real = d_real/sqrt(real(numsamp,dp))    

    allocate(N_filter_real(-numsamp+1:numsamp-1))
    N_filter_real = 0
    do i = -numsamp+1, numsamp-1
       N_filter_real(i) = d_real(abs(i))
    end do
    deallocate(d_real)
    deallocate(d_fft)

    ! DO THIS SUBTRACTION ELSEWHERE!!!
    !N_filter_real = N_filter_real - sum(N_filter_real(-n:n))/size(N_filter_real(-n:n))

  end subroutine get_inv_noise_filter_real

  ! Given the power spectrum N_filter_fft of a filter, apply that to
  ! tod by using the convolution theorem. Tod should have length numsamp.
  subroutine apply_filter_fft(N_filter_fft, tod)
    implicit none

    real(dp),     dimension(0:), intent(in)    :: N_filter_fft
    real(dp),     dimension(1:), intent(inout) :: tod

    integer(i4b) :: n
    complex(dpc), allocatable, dimension(:) :: tod_fft

    n = size(tod)/2+1
    allocate(tod_fft(n))

    call fft(tod, tod_fft, 1)
    tod_fft = tod_fft*N_filter_fft
    call fft(tod, tod_fft, -1)
    deallocate(tod_fft)
  end subroutine apply_filter_fft

  ! Normal PROPERTY = value pairs are separated by semicolons and whitespace
  subroutine read_filters_str(line, f)
    implicit none
    character(len=*)    :: line
    type(filter_params) :: f
    integer(i4b)        :: i, j, n
    character(len=512), dimension(:), allocatable :: p
    n = num_tokens(line, ";", group="''" // '""')
    allocate(p(n))
    call get_tokens(line, ";", p, group="''" // '""')

    call get_parameter_arr(p, 'NU_HIGHPASS_IN_SCANFREQ',       par_dp=f%nu_highpass)
    call get_parameter_arr(p, 'ALPHA_HIGHPASS',                par_dp=f%alpha_highpass)
    call get_parameter_arr(p, 'NU_LOWPASS',                    par_dp=f%nu_lowpass)
    call get_parameter_arr(p, 'ALPHA_LOWPASS',                 par_dp=f%alpha_lowpass)
    call get_parameter_arr(p, 'APPLY_HIGHPASS_FILTER',         par_lgt=f%apply_highpass)
    call get_parameter_arr(p, 'APPLY_LOWPASS_FILTER',          par_lgt=f%apply_lowpass)
    call get_parameter_arr(p, 'IGNORE_ONEOVERF',               par_lgt=f%ignore_oneoverf)
    call get_parameter_arr(p, 'IGNORE_DIODE_CORRELATIONS',     par_lgt=f%ignore_diode_corr)
    call get_parameter_arr(p, 'IGNORE_DIODE_CORRFREQ',         par_lgt=f%ignore_diode_corr_freq)
    call get_parameter_arr(p, 'USE_PRECOMPUTED_FILTER',        par_lgt=f%use_precomp_filter)
    deallocate(p)
  end subroutine

  subroutine read_filters_parfile(paramfile, f, t)
    implicit none
    character(len=*)    :: paramfile
    type(filter_params) :: f
    type(noise_template):: t
    integer(i4b)        :: unit
    unit = getlun()
    call get_parameter(unit, paramfile, 'NU_HIGHPASS_IN_SCANFREQ',       par_dp=f%nu_highpass)
    call get_parameter(unit, paramfile, 'ALPHA_HIGHPASS',                par_dp=f%alpha_highpass)
    call get_parameter(unit, paramfile, 'NU_LOWPASS',                    par_dp=f%nu_lowpass)
    call get_parameter(unit, paramfile, 'ALPHA_LOWPASS',                 par_dp=f%alpha_lowpass)
    call get_parameter(unit, paramfile, 'APPLY_HIGHPASS_FILTER',         par_lgt=f%apply_highpass)
    call get_parameter(unit, paramfile, 'APPLY_LOWPASS_FILTER',          par_lgt=f%apply_lowpass)
    call get_parameter(unit, paramfile, 'APPLY_SPIKEFILTER',             par_lgt=f%apply_spike)    
    if (f%apply_spike) then
       call get_parameter(unit, paramfile, 'USE_WOODBURY_SPIKEFILTER',      par_lgt=f%use_woodbury_spikes)
    end if
    call get_parameter(unit, paramfile, 'IGNORE_ONEOVERF',               par_lgt=f%ignore_oneoverf)
    call get_parameter(unit, paramfile, 'IGNORE_DIODE_CORRELATIONS',     par_lgt=f%ignore_diode_corr)
    call get_parameter(unit, paramfile, 'IGNORE_DIODE_CORRFREQ',         par_lgt=f%ignore_diode_corr_freq)
    call get_parameter(unit, paramfile, 'NOISE_CORR_FUNC_CUTOFF',        par_dp=corr_func_cutoff)    
    call get_parameter(unit, paramfile, 'USE_PRECOMPUTED_FILTER',        par_lgt=f%use_precomp_filter)
    call get_parameter(unit, paramfile, 'CORRELATION_LENGTH_CAP',        &
         & par_int=corrlen_cap, desc="The maximum length allowed for the correlation function. "// &
         & "Negative values disable the cap. If enabled, the resulting covariance matrix will "// &
         & "be incorrect. Thus, for production runs, this option should always be -1.")
    call get_parameter(unit, paramfile, 'USE_TEMPLATES',                 par_lgt=t%use_templates)  
  end subroutine read_filters_parfile

  subroutine get_F_filter(pow, f_scan, samprate, F, filter_opts)
    implicit none

    real(dp),                intent(in)  :: pow, f_scan, samprate
    real(dp), dimension(1:), intent(out) :: F
    type(filter_params),     optional    :: filter_opts
    
    real(dp)            :: nu, a, fact
    integer(i4b)        :: i, nfreq, nsamp
    type(filter_params) :: opts
    logical(lgt)        :: inverse=.false.

    opts = default_filter; if (present(filter_opts)) opts = filter_opts

    nfreq = size(F)
    nsamp = 2*(nfreq-1)
    F     = 1.d0

    if (opts%apply_highpass) then
       F(1) = 0.d0
       do i = 2, nfreq
          nu = ind2freq(i, samprate, nfreq)
          fact = min(max((1+(nu/(f_scan*opts%nu_highpass))**opts%alpha_highpass), 1.d-100), 1.d100)
          F(i) = F(i) * fact**(-pow)
       end do
    end if

    if (opts%apply_lowpass) then
       do i = 1, nfreq
          nu = ind2freq(i, samprate, nfreq)
          fact = min(max((1+(nu/opts%nu_lowpass)**(-opts%alpha_lowpass)), 1.d-100), 1.d100)
          F(i) = F(i) * fact**(-pow)
       end do
    end if

! TMR: testing adding spikefilter to this function to include it in various chisq tests in cesvalidate. 
! This routine is not used by any other program.
! Update: No longer using this in cesval either since this uses the one-size-fits-all filter
    if (pow < 0.d0) inverse=.true.
    if (opts%apply_spike) call spike_filter_fft(inverse, opts%spikefreqs, opts%num_spikes, nsamp, samprate, nfreq, F)

  end subroutine get_F_filter

  subroutine get_N_filter(samprate, sigma_0, fknee, alpha, pow, N_filter, highpass)
    implicit none

    logical(lgt),                intent(in)  :: highpass
    real(dp),                    intent(in)  :: samprate, sigma_0, fknee, alpha, pow
    real(dp),     dimension(0:), intent(out) :: N_filter
    integer(i4b) :: i
    real(dp)     :: f

    do i = 1, size(N_filter)-1
       f = ind2freq(i+1, samprate, size(N_filter))
       N_filter(i) = (sigma_0**2*(1+(f/fknee)**alpha))**pow
    end do

    N_filter(0) = N_filter(1)
    if(highpass) N_filter(0) = 0.d0

  end subroutine get_N_filter

  subroutine apply_filter_fft_matrix(N_filter_fft, tod, N_d)
    implicit none
    real(dp),         dimension(:,:,0:), intent(in)   :: N_filter_fft
    real(dp),         dimension(:,:)                  :: tod
    real(dp),         dimension(:,:),    intent(out)  :: N_d
    complex(dpc),     dimension(:,:),    allocatable  :: tod_fft
    integer(i4b)                                      :: i, n

    n = (size(tod,1)/2+1)
    allocate(tod_fft(n,size(tod,2)))
    call fft_multi(tod, tod_fft, 1)
    ! Multiply with inverse noise matrix in Fourier space
    do i = 1, n
       tod_fft(i,:)   = matmul(N_filter_fft(:,:,i-1), tod_fft(i,:))
    end do
    call fft_multi(N_d, tod_fft, -1)
    deallocate(tod_fft)
  end subroutine apply_filter_fft_matrix

  subroutine apply_filter_real_matrix(tcorr, itod, otod)
    implicit none
    real(dp),   intent(in) :: tcorr(:,:,:), itod(:,:)
    real(dp),   intent(out):: otod(:,:)
    integer(i4b) :: i, j, k, n, ndi, nc
    nc = (size(tcorr,1)-1)/2
    n  = size(itod,1)
    ndi= size(itod,2)
    do i = 1, n
       otod(i,:) = 0
       do j = max(1, i-nc), min(n, i+nc)
          k = j-i + nc+1
          otod(i,:) = otod(i,:) + matmul(tcorr(k,:,:), itod(j,:))
       end do
    end do
  end subroutine

  subroutine apply_filter_fft_multi(N_filter_fft, tod, N_d)
    implicit none
    real(dp),         dimension(:,:,0:), intent(in)   :: N_filter_fft
    real(dp),         dimension(:,:,:)                :: tod
    real(dp),         dimension(:,:,:),  intent(out)  :: N_d
    complex(dpc),     dimension(:,:,:),  allocatable  :: tod_fft
    integer(i4b)                                      :: i, j, n

    n = (size(tod,1)/2+1)
    allocate(tod_fft(n,size(tod,2),size(tod,3)))
    call fft_multi(tod, tod_fft, 1)
    ! Multiply with inverse noise matrix in Fourier space
    !$OMP PARALLEL PRIVATE(j,i)
    do j = 1, size(tod,3)
       !$OMP DO
       do i = 1, n
          tod_fft(i,:,j)   = matmul(N_filter_fft(:,:,i-1), tod_fft(i,:,j))
       end do
       !$OMP END DO
    end do
    !$OMP END PARALLEL
    call fft_multi(N_d, tod_fft, -1)
    deallocate(tod_fft)
  end subroutine apply_filter_fft_multi

  subroutine apply_filter_real_multi(tcorr, itod, otod)
    implicit none
    real(dp),         dimension(:,:,0:), intent(in)   :: tcorr
    real(dp),         dimension(:,:,:)                :: itod
    real(dp),         dimension(:,:,:),  intent(out)  :: otod
    integer(i4b)                                      :: i
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i = 1, size(itod,3)
       call apply_filter_real_matrix(tcorr, itod(:,:,i), otod(:,:,i))
    end do
    !$OMP END PARALLEL
  end subroutine

  subroutine initialize_filter_assembly(assembly, filtopts)
    implicit none
    type(quiet_assembly),  intent(inout) :: assembly
    type(filter_params),   optional      :: filtopts

    integer(i4b) :: i, j, k, ndi, nsamp, n, n_max, n_max_F_sq, nf, nbin
    real(dp)     :: nu, fft_low_freq, fft_low_alpha, fft_high_freq, fft_high_alpha, scanfreq, f_scan
    real(dp)     :: srate
    type(noise_corr), allocatable, dimension(:,:) :: inv_N_corr, inv_N_corr_F_sq
    real(dp),         allocatable, dimension(:)   :: flt
    real(sp),         allocatable, dimension(:,:) :: temps_full
    type(filter_params)                           :: filter
    real(dp)     :: t1, t2

    filter = default_filter; if(present(filtopts)) filter = filtopts

    ndi   = assembly%num_diodes
    nsamp = assembly%numsamp
    srate = assembly%samprate
    nf    = nsamp/2+1

    ! Set up apodized highpass filter parameters
    if (filter%use_precomp_filter) then
       fft_high_freq  = assembly%fft_high_freq(1)
!fixme ingunn       fft_high_freq  = assembly%fft_high_freq(1)
!       fft_high_freq  = max(assembly%fft_high_freq(1), 1.5d0*assembly%scanfreq)
       fft_high_alpha = assembly%fft_high_alpha(1)
       fft_low_freq   = assembly%fft_low_freq(1)
       fft_low_alpha  = assembly%fft_low_alpha(1)
    else
       fft_high_freq  = filter%nu_highpass * assembly%scanfreq
       fft_high_alpha = filter%alpha_highpass
       fft_low_freq   = filter%nu_lowpass
       fft_low_alpha  = filter%alpha_lowpass
    end if

    ! Initialize Fourier space noise correlation function
    if(allocated(assembly%N_corr_fft))      deallocate(assembly%N_corr_fft)
    if(allocated(assembly%N_corr_F_sq_fft)) deallocate(assembly%N_corr_F_sq_fft)
    if(allocated(assembly%inv_N_filter))    deallocate(assembly%inv_N_filter)
    if(allocated(assembly%inv_N_fft))       deallocate(assembly%inv_N_fft)
    allocate(assembly%N_corr_fft     (ndi,ndi, 0:nf-1))
    allocate(assembly%N_corr_F_sq_fft(ndi, ndi, 0:nf-1))
    allocate(assembly%inv_N_filter   (0:nf-1,ndi))
    allocate(assembly%inv_N_fft      (0:nf-1,ndi))
    assembly%N_corr_fft      = 0.d0
    assembly%N_corr_F_sq_fft = 0.d0
    assembly%inv_N_filter    = 1.d0 ! Leftovers?
    assembly%inv_N_fft       = 0.d0

!    call wall_time(t2)
!    write(*,*) 'init = ', t2-t1

    ! Set up the basic noise filters
!    call wall_time(t1)

    ! Set up correlations
    if(filter%ignore_diode_corr) then
       assembly%N_corr_fft = spread(get_identity(ndi),3,nf)
    elseif(filter%ignore_diode_corr_freq) then
       assembly%N_corr_fft = spread(assembly%corr(size(assembly%corr,1),:,:),3,nf)
    else
       do i = 1, size(assembly%corr,1)
          assembly%N_corr_fft(:,:,i-1) = assembly%corr(i,:,:)
       end do
    end if

    do i = 1, ndi
       if(.not. filter%ignore_oneoverf) then
          ! We fetch the non-inverse noisemod here
          call get_inv_noise_filter_fft(.false., nsamp, srate, assembly%sigma0(i), &
            & assembly%fknee(i), assembly%alpha(i), assembly%inv_N_fft(:,i), filter%apply_highpass)
       else
          assembly%inv_N_fft(:,i) = assembly%sigma0(i)**2
       end if
    end do

    if(templates%use_templates) then
       call get_templates(temps_full,nf,srate,assembly%diodes_abs) ! ass%diodes_abs is an array of 3-4 elements with the diode index of the assembly diodes
       assembly%inv_N_fft = assembly%inv_N_fft * real(temps_full,dp)
       deallocate(temps_full)
    end if

    do i = 1, ndi
       do j = 1, ndi
          assembly%N_corr_fft(i,j,:) = assembly%N_corr_fft(i,j,:)*sqrt(assembly%inv_N_fft(:,i)*assembly%inv_N_fft(:,j))
       end do
    end do

!    call wall_time(t1)
    assembly%N_corr_F_sq_fft = assembly%N_corr_fft

    ! Apply extra filters
    allocate(flt(0:nsamp/2))
    flt = 1

    if(filter%apply_lowpass)  call apodize_filter_fft(.false., nsamp, &
       & srate, fft_low_freq,  fft_low_alpha,  .false., flt)
    if(filter%apply_highpass) call apodize_filter_fft(.false., nsamp, &
       & srate, fft_high_freq, fft_high_alpha, .true., flt)

    ! flt is now
    ! \               /
    !  \             /
    !   \-----------/
    ! i.e. it is the inverse weighting filter.

    if(filter%apply_spike .and. (.not. filter%use_woodbury_spikes)) &
         & call spike_filter_fft(.true., filter%spikefreqs, filter%num_spikes, nsamp, srate, nf, flt)
    
    do i = 1, ndi
       assembly%inv_N_fft(:,i)       = (assembly%inv_N_fft(:,i)  / flt)
       do j = 1, ndi
          assembly%N_corr_F_sq_fft(j,i,:) = assembly%N_corr_F_sq_fft(j,i,:) * flt**2
          assembly%N_corr_fft     (j,i,:) = assembly%N_corr_fft     (j,i,:) * flt   
       end do
    end do

    ! Invert Fourier space correlation function
    do i = 0, nsamp/2

! TMR: If we want to run with offset subtraction we need to kill this thing, but for the moment the results are raving.
      if(.not. filter%ignore_oneoverf .and. i == 0) then
          assembly%N_corr_fft     (:,:,0) = 0.d0    ! Remove offset by hand
          assembly%N_corr_F_sq_fft(:,:,0) = 0.d0    ! Remove offset by hand
       else
          call eigen_pow(assembly%N_corr_fft(:,:,i),      -1d0, assembly%N_corr_fft(:,:,i))
          call eigen_pow(assembly%N_corr_F_sq_fft(:,:,i), -1d0, assembly%N_corr_F_sq_fft(:,:,i))
       end if
    end do

    ! These are used when calculating div. Div ignores pixel correlations
    ! and thus cannot represent time correlations properly. But it should
    ! be as close as possible to the operation rhs is exposed to, so that
    ! inv(div)*rhs is as close as possible to the map. rhs = P*F*iN*d,
    ! so its amplitude is reduced by mean(F*iN). We will use mean(F*iN)
    ! as a representative value of the noise+filter.
    allocate(assembly%N_eff(ndi, ndi), assembly%inv_N_eff(ndi,ndi))

    assembly%inv_N_eff = sum(assembly%N_corr_fft(:,:,1:),3)/(size(assembly%N_corr_fft,3)-1.d0)
    assembly%N_eff = assembly%inv_N_eff
    call eigen_pow(assembly%N_eff, -1d0, assembly%N_eff)
    !call invert_matrix(assembly%N_eff, cholesky=.true.)

    ! Compute real-space correlation function
!    call wall_time(t1)
    allocate(inv_N_corr(ndi,ndi))
    allocate(inv_N_corr_F_sq(ndi,ndi))
    n_max      = 0

    do i = 1, ndi
       do j = i, ndi
          call get_inv_noise_filter_real(nsamp, assembly%N_corr_fft(i,j,:), &
               & inv_N_corr(i,j)%n, inv_N_corr(i,j)%N_corr)
!call dump_matrix(assembly%N_corr_fft(1,1,:),"foo.txt")
!call dump_matrix(inv_N_corr(1,1)%N_corr,"bar.txt")
!stop
          call get_inv_noise_filter_real(nsamp, assembly%N_corr_F_sq_fft(i,j,:), &
               & inv_N_corr_F_sq(i,j)%n, inv_N_corr_F_sq(i,j)%N_corr)
          n_max = max(n_max, inv_N_corr(i,j)%n)
          n_max = max(n_max, inv_N_corr_F_sq(i,j)%n)
!write(*,*) i, j, inv_N_corr(i,j)%n, inv_N_corr_F_sq(i,j)%n, n_max
       end do
    end do
    if(corrlen_cap >= 0) n_max = min(n_max,corrlen_cap)
    if (n_max > 0) then
       ! We subtract the mean from the correlation function in order to ensure that absolute offsets in the 
       ! time domain are nullified, and to be consistent with setting the zeroth frequency to zero.
       ! Need to do it here, after n_max has been finally settled.
       do i = 1, ndi
          do j = i, ndi
             inv_N_corr(i,j)%N_corr      = inv_N_corr(i,j)%N_corr      - &
                  & sum(inv_N_corr(i,j)%N_corr(-n_max:n_max)) / (2*n_max+1)
             inv_N_corr_F_sq(i,j)%N_corr = inv_N_corr_F_sq(i,j)%N_corr - &
                  & sum(inv_N_corr_F_sq(i,j)%N_corr(-n_max:n_max)) / (2*n_max+1)
          end do
       end do
    end if

!    call wall_time(t1)
    if (allocated(assembly%N_corr_real))      deallocate(assembly%N_corr_real)
    if (allocated(assembly%N_corr_F_sq_real)) deallocate(assembly%N_corr_F_sq_real)
    allocate(assembly%N_corr_real     (-n_max:n_max, ndi, ndi))
    allocate(assembly%N_corr_F_sq_real(-n_max:n_max, ndi, ndi))
    assembly%noise_corr_length = n_max
    assembly%N_corr_real       = 0.d0
    assembly%N_corr_F_sq_real  = 0.d0
    do i = 1, ndi
       do j = i, ndi
          assembly%N_corr_real(-n_max:n_max,i,j) = inv_N_corr(i,j)%N_corr(-n_max:n_max)
          assembly%N_corr_real(:,j,i)            = assembly%N_corr_real(:,i,j)
          assembly%N_corr_F_sq_real(-n_max:n_max,i,j) = inv_N_corr_F_sq(i,j)%N_corr(-n_max:n_max)
          assembly%N_corr_F_sq_real(:,j,i)    = assembly%N_corr_F_sq_real(:,i,j)
          deallocate(inv_N_corr(i,j)%N_corr)
          deallocate(inv_N_corr_F_sq(i,j)%N_corr)
       end do
    end do
    deallocate(inv_N_corr)
    deallocate(inv_N_corr_F_sq)
  end subroutine initialize_filter_assembly

end module quiet_filter_mod
