module quiet_noise_estimation_mod
  use healpix_types
  use quiet_fft_mod
  use math_tools
  use quiet_utils
  use quiet_mpi_mod
  use quasi_newton_mod
  implicit none

  real(dp),                            private :: samprate_int
  integer(i4b),                        private :: ind_max
  real(dp),                            private :: sigma1, srange(2)
  real(dp), allocatable, dimension(:), private :: f, freqs
  real(dp), allocatable, dimension(:), private :: mask

  integer(i4b), private :: FIT_ALPHA  = 1
  integer(i4b), private :: FIT_F_KNEE = 2

  integer(i4b), parameter, private :: npar = 2
  real(dp),                private :: nu_low_wn, nu_high_wn
  real(dp),                private :: prior(npar,2), f_scan_int, df_scan_int
  logical(lgt),            private :: initialized = .false.

contains

  subroutine initialize_noise_estimation_mod(parfile)
    implicit none
    character(len=*) :: parfile

    if(initialized) return
    call get_parameter(0, parfile, 'NOISE_EST_LOW_NU', par_dp=nu_low_wn)
    call get_parameter(0, parfile, 'NOISE_EST_HIGH_NU', par_dp=nu_high_wn)

    ! Set up priors
    prior(FIT_ALPHA,1)   = -10.d0
    prior(FIT_ALPHA,2)   = -0.1d0
    prior(FIT_F_KNEE,1)  = 1.d-5
    prior(FIT_F_KNEE,2)  = 1.d0
    ! We're not using these anymore

    initialized = .true.
  end subroutine initialize_noise_estimation_mod

  subroutine fit_1overf_profile(samprate, f_scan, df_scan, sigma0, alpha, f_knee, tod, tod_ps, &
       & cnum, diode, chisq_out, apply_scanmask, refit, limits)
    implicit none

    real(dp),                intent(in)            :: samprate, f_scan, df_scan
    integer(i4b),            intent(in),  optional :: cnum, diode, limits(2) ! TMR - optional index limits
    real(dp), dimension(1:),              optional :: tod, tod_ps
    real(dp),                intent(inout)         :: sigma0, alpha, f_knee ! Changed the intent since these are also input in the case of refit=.true.
    real(dp),                intent(out), optional :: chisq_out
    logical(lgt),            intent(in),  optional :: apply_scanmask, refit   ! optional for second noise model fit!

    integer(i4b) :: i, j, n, m, numsamp, numbin, ind1, ind2, ierr, iter, numiter
    logical(lgt) :: accept, est_sigma0_from_low_freq, use_grid, apply_scanmask_, refit_ ! TMR
    real(dp)     :: alpha_min, alpha_max, f_knee_min, f_knee_max, dalpha, df_knee, nu, P_nu, &
         & dnu_wn, fret, gtol
    real(dp)     :: ratio_median_to_mean, chisq, t1, t2, dx, p_full(0:2)
    real(dp),    allocatable, dimension(:)   :: params_lnL, params_mean, rms, alphas, f_knees, &
         & p, f_sort, a
    complex(dp), allocatable, dimension(:)   :: ft
    real(dp),    allocatable, dimension(:)   :: N_filter, N_ft, lnL_s, x
    real(dp),    allocatable, dimension(:,:) :: lnL
    integer(i4b), dimension(2) :: max_pos

    apply_scanmask_ = .true.; if(present(apply_scanmask)) apply_scanmask_ = apply_scanmask
    refit_ = .false.;         if(present(refit)) refit_ = refit 

    accept = .false.
    ! Check that the module isn't dead
    if (present(tod)) then
       if (is_nan(sum(abs(tod)))) then
          sigma0 = 0.d0
          alpha  = 0.d0
          f_knee = 0.d0
          return
       end if
    else if (present(tod_ps)) then
       if (is_nan(sum(abs(tod_ps)))) then
          sigma0 = 0.d0
          alpha  = 0.d0
          f_knee = 0.d0
          return
       end if       
    end if

    f_scan_int   = f_scan
    df_scan_int  = df_scan
    samprate_int = samprate

    ! Prepare power spectrum
    if (present(tod_ps)) then
       n = size(tod_ps)
       allocate(f(0:n-1))
       f = tod_ps
    else if (present(tod)) then
       n = (size(tod)/2+1)
       allocate(ft(0:n-1), f(0:n-1))
       call fft(tod, ft, 1)
       call extract_powspec(ft, f)       
       deallocate(ft)
    else
       write(*,*) 'quiet_noise_estimation_mod: Error -- either tod or fft must be present'
       stop
    end if

    ! If all entries are zero, return with zero variance
    if (all(f < 1.d-30)) then
       sigma0 = 0.d0
       alpha  = -1.d0
       f_knee =  0.1d0
       deallocate(f)
       return
    end if

    ! Set up scan frequency mask
    allocate(mask(0:n-1), freqs(0:n-1))
    do i = 0, n-1
       freqs(i) = ind2freq(i+1, samprate, n)
    end do
    dnu_wn = ind2freq(2, samprate, n)
    mask   = 1.d0
    if(apply_scanmask_) then
       ind1   = 0
       ind2   = 0
       nu     = 0.d0
       do while (ind2 < n-1)
          nu     = nu + f_scan_int
          ind1   = max(nint((nu-df_scan_int) / dnu_wn),1)
          ind2   = min(nint((nu+df_scan_int) / dnu_wn),n-1)
          !if (nu == f_scan_int) mask(ind1:ind2) = 0.d0
          mask(ind1:ind2) = 0.d0
       end do
    end if
    mask(0)   = 0.d0
    mask(n-1) = 0.d0

    if (present(limits)) then
       mask(0:limits(1)) = 0.d0
       mask(limits(2):n-1) = 0.d0
    end if

    allocate(p(0:npar-1))

    if (.not. refit_) then

       ! Compute white noise level
       ! est_sigma0_from_low_freq = .true.                      ! What's this? obsolete?
       dnu_wn = ind2freq(2, samprate, n)                        ! This has already been done.
       ind1   = max(nint(nu_low_wn / dnu_wn),2)
       ind2   = min(nint(nu_high_wn / dnu_wn),n-2)
       !100 continue                                               ! obsolete?
              
       ! This is the measured sigma0 in the cmb region under the assumption
       ! that we have no 1/f component. To get the true sigma0, we
       ! will have to compensate for the bias this introduces.
       sigma1 = sqrt(sum(f(ind1:ind2)*mask(ind1:ind2)) / sum(mask(ind1:ind2)))
       srange = [ ind1, ind2 ] * dnu_wn       ! This is no longer in use since introducing the joint fit of all 3 parameters. No bias correction needed.
           
!       ind_max = freq2ind(1d0, samprate, n)
       ind_max = freq2ind(4.5d0, samprate, n) ! Moving the limit since we are using joint fit of all 3 parameters

       ! Fit 1/f component, using data up to 1 Hz -- this sucks up weather to some extent. 
       ! Should come back to this, and implement a noise model that actually fits for both instrumental noise
       ! and weather

       p_full(0) = log(sigma1**2)    ! variance
       p_full(1) = log(0.001d0) ! f_knee
       p_full(2) = log(2.d0) ! alpha

    else
       ! Refit: Do the noise fit on full filtered spectrum.
       if(present(limits)) then
          ind_max = limits(2)
       else
          ind_max = n-1 
       end if
       ! Using results from first fit as initial values for refit
       p_full(0) = log(sigma0**2)    ! variance
       p_full(1) = log(f_knee) ! f_knee
       p_full(2) = log(-alpha) ! alpha
    end if

    gtol = 1.d-5
    numiter = 5
    do i = 1, numiter
       !call dfpmin(p,gtol,iter,fret,lnL_noise_powell,dlnL_noise_powell, ierr=ierr)
       !write(*,*) sigma1, real(p,sp), iter
       !write(*,*) 'fret1 = ', fret

       !p_full(1) = log(0.001d0) ! f_knee
       !p_full(2) = log(2.d0) ! alpha
       call dfpmin(p_full,gtol,iter,fret,lnL_noise_powell_full,dlnL_noise_powell_full, ierr=ierr)
!       write(*,*) i, real(sqrt(abs(p_full(0))),sp), real(p_full(1:2),sp), iter, ierr
       !write(*,*) 'fret2 = ', fret

!       stop
       if (ierr == 1) then
          if (present(cnum)) then
             write(*,*) 'Error: NaNs in noise fit -- this should be fixed, cnum and diode= ', cnum, diode
!             stop
          else
             write(*,*) 'Error: NaNs in noise fit -- this should be fixed'
          end if

          sigma0 = 0; alpha = 0; f_knee = 0
          goto 101                                    ! This aborts the noise estimation
          stop
       else if (ierr == 2) then
          sigma0 = sqrt(exp(p_full(0)))
          f_knee = exp(p_full(1))
          alpha  = -exp(p_full(2))

          ! Check goodness-of-fit
          ind1   = nint(0.2d0 / dnu_wn)
          ind2   = nint(nu_high_wn / dnu_wn)
          accept = check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, f, chisq)
          if(dtest(2) .and. abs(chisq) > 10.d0 .and. i == numiter) &
               & write(*,*) 'Warning: Noise fit did not converge properly. Estimates may be OK anyway'
       else
          sigma0 = sqrt(exp(p_full(0)))
          f_knee = exp(p_full(1))
          alpha  = -exp(p_full(2))
          exit
       end if
    end do
    deallocate(p)

    ! Add a safety for extremely low f_knee's, to avoid NaN's in later programs
    if (f_knee < 1.d-10 .or. alpha > -0.01d0) then
       ! Return white noise
       f_knee = 1.d-10
       alpha  = -10.d0
    end if

!    open(57,file='pow.dat')
!    do i = 1, n-1
!       write(57,*) freqs(i), f(i), sigma0**2 * (1 + (freqs(i)/f_knee)**alpha)
!    end do
!    close(57)

    if (present(chisq_out)) then
       ! Check goodness-of-fit
       ind1   = nint(0.2d0 / dnu_wn)
       ind2   = nint(nu_high_wn / dnu_wn)
       accept = check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, f, chisq)
       chisq_out = chisq
    end if

101 deallocate(f)
    deallocate(mask, freqs)

  end subroutine fit_1overf_profile

  function lnL_noise_powell(p)
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp)             :: lnL_noise_powell

    integer(i4b) :: i
    real(dp)     :: P_nu, alpha, f_knee, x, current_sigma0

    f_knee =  exp(p(1))
    alpha  = -exp(p(2))

    if (f_knee > 3.d0 .or. alpha < -10.d0) then
       lnL_noise_powell = 1.d30
       return
    end if

    lnL_noise_powell = 0.d0
    current_sigma0 = fix_sigma(sigma1, f_knee, alpha, srange)
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       P_nu      = current_sigma0**2 * (1.d0 + (freqs(i)/f_knee)**alpha)
       lnL_noise_powell = lnL_noise_powell + f(i) / P_nu + log(P_nu)
    end do

    ! Add prior on alpha
    lnL_noise_powell = lnL_noise_powell / sum(mask(1:ind_max))
  end function lnL_noise_powell

  function dlnL_noise_powell(p) 
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp), dimension(size(p))        :: dlnL_noise_powell

    integer(i4b) :: i
    real(dp)     :: nu, P_nu, alpha, f_knee, x, dPda, dPdf, current_sigma0, dS2da, dS2df

    f_knee = exp(p(1))
    alpha  = -exp(p(2))

    dlnL_noise_powell = 0.d0
    current_sigma0 = fix_sigma(sigma1, f_knee, alpha, srange)
    call deriv_sigma2(sigma1, f_knee, alpha, srange, dS2da, dS2df)
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       nu        = freqs(i)
       P_nu      = current_sigma0**2 * (1.d0 + (nu/f_knee)**alpha)
       dPdf      = -alpha * current_sigma0**2 * nu**alpha * f_knee**(-alpha-1) + &
                 & P_nu/current_sigma0**2 * dS2df
       dPda      = current_sigma0**2 * (nu/f_knee)**alpha * log(nu/f_knee) + &
                 & P_nu/current_sigma0**2 * dS2da
       dlnL_noise_powell(1)     = dlnL_noise_powell(1) + dPdf/P_nu * (-f(i)/P_nu + 1.d0)
       dlnL_noise_powell(2)     = dlnL_noise_powell(2) + dPda/P_nu * (-f(i)/P_nu + 1.d0)
    end do
    dlnL_noise_powell(1) =  dlnL_noise_powell(1) * f_knee
    dlnL_noise_powell(2) =  dlnL_noise_powell(2) * alpha

    dlnL_noise_powell = dlnL_noise_powell / sum(mask(1:ind_max))
  end function dlnL_noise_powell


  function lnL_noise_powell_full(p)
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp)             :: lnL_noise_powell_full

    integer(i4b) :: i
    real(dp)     :: P_nu, alpha, f_knee, x, sigma_sq

    sigma_sq =  exp(p(1))
    f_knee   =  exp(p(2))
    alpha    = -exp(p(3))

    !write(*,*) 'par1 = ', sigma_sq, f_knee, alpha

! TMR: changed prior on alpha from -10 to -6
    if (f_knee > 3.d0 .or. alpha < -10.d0 .or. sigma_sq == 0.d0) then
       lnL_noise_powell_full = 1.d30
       return
    end if

    lnL_noise_powell_full = 0.d0
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       P_nu      = sigma_sq * (1.d0 + (freqs(i)/f_knee)**alpha)
       lnL_noise_powell_full = lnL_noise_powell_full + f(i) / P_nu + log(P_nu)
    end do

    ! Add prior on alpha
    lnL_noise_powell_full = lnL_noise_powell_full / sum(mask(1:ind_max))
    !write(*,*) 'lnL = ', lnL_noise_powell_full
  end function lnL_noise_powell_full

  function dlnL_noise_powell_full(p) 
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp), dimension(size(p))        :: dlnL_noise_powell_full

    integer(i4b) :: i
    real(dp)     :: nu, P_nu, alpha, f_knee, x, dPda, dPdf, current_sigma0, &
         & dS2da, dS2df, sigma_sq, dLdP

    sigma_sq = exp(p(1))
    f_knee   = exp(p(2))
    alpha    = -exp(p(3))

    !write(*,*) 'par2 = ', sigma_sq, f_knee, alpha

    dlnL_noise_powell_full = 0.d0
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       nu        = freqs(i)
       P_nu      = sigma_sq * (1.d0 + (nu/f_knee)**alpha)
       dLdP      = -f(i)/P_nu**2 + 1.d0 / P_nu
       dlnL_noise_powell_full(1)     = dlnL_noise_powell_full(1) + dLdP * (1.d0 + (nu/f_knee)**alpha) * sigma_sq
       dlnL_noise_powell_full(2)     = dlnL_noise_powell_full(2) - dLdP * sigma_sq * alpha * (nu/f_knee)**alpha
       dlnL_noise_powell_full(3)     = dlnL_noise_powell_full(3) + dLdP * sigma_sq * alpha * (nu/f_knee)**alpha * log(nu/f_knee)
    end do

    dlnL_noise_powell_full = dlnL_noise_powell_full / sum(mask(1:ind_max))
    !write(*,*) 'dlnL = ', dlnL_noise_powell_full
  end function dlnL_noise_powell_full

  ! Compensate for bias in sigma due to 1/f profile,
  ! when sigma1 was measured in range srange.
  !  int{a:b}(sigma0**2*(1+(x/fknee)**alpha)) =
  !  sigma0**2*(b-a + (b**(alpha+1)-a**(alpha+1))/(fknee*(alpha+1))) = 
  !  sigma1**2*(b-a)
  ! The case alpha = -1 should be handled as an exception here.
  function fix_sigma(sigma1, fknee, alpha, srange) result(sigma0)
    implicit none
    real(dp)     :: sigma1, sigma0, fknee,alpha, srange(2), a, b
    a = srange(1)
    b = srange(2)
    sigma0 = sigma1/sqrt(1 + (b**(alpha+1)-a**(alpha+1))/(fknee**alpha*(alpha+1)*(b-a)))
  end function
  subroutine deriv_sigma2(sigma1, fknee, alpha, srange, dS2df, dS2da)
    implicit none
    real(dp)     :: sigma1, sigma0, fknee,alpha, srange(2), a, b, dS2df, dS2da
    real(dp)     :: d, da, a1, fa
    a = srange(1)
    b = srange(2)
    d = b-a; a1 = alpha+1; da = b**a1 - a**a1; fa = fknee**(-alpha)
    sigma0 = sigma1/sqrt(1+da/(d*a1)*fa)
    dS2df = alpha/a1 * sigma0**4/sigma1**2*da/d*fknee**(-alpha-1)
    dS2da = -sigma0**4/sigma1**2/d*(1/a1*fa*(log(b/fknee)*b**a1-log(a/fknee)*a**a1)-&
     & 1/a1**2*(b**a1-a**a1)/fknee**alpha)
  end subroutine

  function check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, spec, chisq)
    implicit none

    integer(i4b),                intent(in)  :: ind1, ind2
    real(dp),                    intent(in)  :: sigma0, alpha, f_knee, samprate
    real(dp),     dimension(0:), intent(in)  :: spec
    real(dp),                    intent(out) :: chisq
    logical(lgt)                             :: check_noise_fit

    integer(i4b) :: i, n
    real(dp)     :: nu, N_fft, f

    n = size(spec)

    ! Compute effective chi-square and expected variance
    chisq = 0.d0
    nu    = 0.d0
    do i = ind1, ind2
       if (mask(i) == 1.d0) then
          f     = ind2freq(i, samprate, n)
          N_fft = sigma0**2 * (1.d0 + (f/f_knee)**alpha)
          chisq = chisq + spec(i) / N_fft
          nu    = nu + 1
       end if
    end do

    chisq = (chisq - nu) / sqrt(nu)

    check_noise_fit = abs(chisq) < 4.d0

  end function check_noise_fit

!  subroutine debug_fit(samprate, f_scan, df_scan, sigma0, f_knee, alpha, ps, y)
!    implicit none
!    real(dp)            :: samprate, f_scan, df_scan, ps(:), dnu_wn
!    real(dp)            :: sigma0, alpha, f_knee, p(2), df, da, y(:,:), v(2)
!    integer(i4b)        :: f_scan_int, df_scan_int, ind1, ind2
!    integer(i4b)        :: i, n, m
!
!    n = size(y,1)
!    m = size(ps)
!    allocate(f(0:m-1), mask(0:m-1), freqs(0:m-1))
!    f = ps
!    mask = 1
!    do i = 0, m-1
!       freqs(i) = ind2freq(i+1, samprate, m)
!    end do
!
!    f_scan_int   = f_scan
!    df_scan_int  = df_scan
!    samprate_int = samprate
!
!    dnu_wn = ind2freq(2, samprate, m)
!    ind1   = nint(nu_low_wn / dnu_wn)
!    ind2   = nint(nu_high_wn / dnu_wn)
!
!    sigma1 = sqrt(sum(f(ind1:ind2)*mask(ind1:ind2)) / sum(mask(ind1:ind2)))
!    srange = [ ind1, ind2 ]
!
!    dnu     = ind2freq(2, samprate, m)
!    ind_max = nint(1.d0 / dnu)
!
!    p = log([ f_knee, -alpha ])
!
!    do i = 1, n
!       df = (i-n/2)*0.01/(n/2)
!       p(1) = log(f_knee+df)
!       y(i,1) = lnL_noise_powell(p)
!v = dlnL_noise_powell(p)
!write(*,'(i5,4e15.7)') i, y(i,1), v, p(1)
!    end do
!    p = log([ f_knee, -alpha ])
!    do i = 1, n
!       da = (i-n/2)*0.1/(n/2)
!       p(2) = log(-alpha-da)
!       y(i,2) = lnL_noise_powell(p)
!v = dlnL_noise_powell(p)
!write(*,'(i5,4e15.7)') i, y(i,2), v, p(2)
!    end do
!    deallocate(f,mask,freqs)
!  end subroutine

!  subroutine fit_1overf_profile(samprate, f_scan, df_scan, sigma0, alpha, f_knee, tod, tod_ps, cnum, diode, chisq_out, apply_scanmask)
!    implicit none
!
!    real(dp),                intent(in)            :: samprate, f_scan, df_scan
!    integer(i4b),            intent(in), optional  :: cnum, diode
!    real(dp), dimension(1:),              optional :: tod, tod_ps
!    real(dp),                intent(out)           :: sigma0, alpha, f_knee
!    real(dp),                intent(out), optional :: chisq_out
!    logical(lgt),            intent(in),  optional :: apply_scanmask
!
!    integer(i4b) :: i, j, n, m, numsamp, numbin, ind1, ind2, ierr, iter, numiter
!    logical(lgt) :: accept, est_sigma0_from_low_freq, use_grid, apply_scanmask_
!    real(dp)     :: alpha_min, alpha_max, f_knee_min, f_knee_max, dalpha, df_knee, nu, P_nu, dnu_wn, fret, gtol
!    real(dp)     :: ratio_median_to_mean, chisq, t1, t2, dx
!    real(dp),    allocatable, dimension(:)   :: params_lnL, params_mean, rms, alphas, f_knees, p, f_sort, a
!    complex(dp), allocatable, dimension(:)   :: ft
!    real(dp),    allocatable, dimension(:)   :: N_filter, N_ft, lnL_s, x
!    real(dp),    allocatable, dimension(:,:) :: lnL
!    integer(i4b), dimension(2) :: max_pos
!
!    apply_scanmask_ = .true.; if(present(apply_scanmask)) apply_scanmask_ = apply_scanmask
!
!!    call dmem('a')
!
!    accept = .false.
!    ! Check that the module isn't dead
!    if (present(tod)) then
!       if (is_nan(sum(abs(tod)))) then
!          sigma0 = 0.d0
!          alpha  = 0.d0
!          f_knee = 0.d0
!          return
!       end if
!    else if (present(tod_ps)) then
!       if (is_nan(sum(abs(tod_ps)))) then
!          sigma0 = 0.d0
!          alpha  = 0.d0
!          f_knee = 0.d0
!          return
!       end if       
!    end if
!
!    f_scan_int   = f_scan
!    df_scan_int  = df_scan
!    samprate_int = samprate
!
!!    call dmem('b')
!
!    ! Prepare power spectrum
!    if (present(tod_ps)) then
!       n = size(tod_ps)
!       allocate(f(0:n-1))
!       f = tod_ps
!    else if (present(tod)) then
!       n = (size(tod)/2+1)
!       allocate(ft(0:n-1), f(0:n-1))
!       call fft(tod, ft, 1)
!       call extract_powspec(ft, f)       
!       deallocate(ft)
!    else
!       write(*,*) 'quiet_noise_estimation_mod: Error -- either tod or fft must be present'
!       stop
!    end if
!
!    ! If all entries are zero, return with zero variance
!    if (all(f < 1.d-30)) then
!       sigma0 = 0.d0
!       alpha  = -1.d0
!       f_knee =  0.1d0
!       deallocate(f)
!       return
!    end if
!
!!    call dmem('c')
!
!    ! Set up scan frequency mask
!    allocate(mask(0:n-1), freqs(0:n-1))
!    do i = 0, n-1
!       freqs(i) = ind2freq(i+1, samprate, n)
!    end do
!    dnu_wn = ind2freq(2, samprate, n)
!    mask   = 1.d0
!    if(apply_scanmask_) then
!       ind1   = 0
!       ind2   = 0
!       nu     = 0.d0
!       do while (ind2 < n-1)
!          nu     = nu + f_scan_int
!          ind1   = max(nint((nu-df_scan_int) / dnu_wn),1)
!          ind2   = min(nint((nu+df_scan_int) / dnu_wn),n-1)
!          mask(ind1:ind2) = 0.d0
!       end do
!    end if
!    mask(0)   = 0.d0
!    mask(n-1) = 0.d0
!!    call dmem('d')
!
!    ! Compute white noise level
!    est_sigma0_from_low_freq = .true.
!    dnu_wn = ind2freq(2, samprate, n)
!    ind1   = nint(nu_low_wn / dnu_wn)
!    ind2   = nint(nu_high_wn / dnu_wn)
!100 continue
!
!    ! This is the measured sigma0 in the cmb region under the assumption
!    ! that we have no 1/f component. To get the true sigma0, we
!    ! will have to compensate for the bias this introduces.
!    sigma1 = sqrt(sum(f(ind1:ind2)*mask(ind1:ind2)) / sum(mask(ind1:ind2)))
!    srange = [ ind1, ind2 ]
!    
!!    call dmem('e')
!
!    ! Fit 1/f component, using data up to 1 Hz -- this sucks up weather to some extent. 
!    ! Should come back to this, and implement a noise model that actually fits for both instrumental noise
!    ! and weather
!    dnu     = ind2freq(2, samprate, n)
!    ind_max = nint(1.d0 / dnu)
!
!    allocate(p(0:npar-1))
!    p(0) = log(0.001d0) ! f_knee
!    p(1) = log(2.d0) ! alpha
!    gtol = 1.d-5
!!    call dmem('f')
!    numiter = 5
!    do i = 1, numiter
!       call dfpmin(p,gtol,iter,fret,lnL_noise_powell,dlnL_noise_powell, ierr=ierr)
!       !write(*,*) ierr, real(p,sp), iter
!       if (ierr == 1) then
!          if (present(cnum)) then
!             write(*,*) 'Error: NaNs in noise fit -- this should be fixed, cnum = ', cnum
!          else
!             write(*,*) 'Error: NaNs in noise fit -- this should be fixed'
!          end if
!          sigma0 = 0; alpha = 0; f_knee = 0
!          goto 101
!          stop
!       else if (ierr == 2) then
!          f_knee = exp(p(0))
!          alpha  = -exp(p(1))
!          sigma0 = sigma1/sqrt(get_bias(f_knee, alpha))
!
!          ! Check goodness-of-fit
!          ind1   = nint(0.2d0 / dnu_wn)
!          ind2   = nint(nu_high_wn / dnu_wn)
!          accept = check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, f, chisq)
!          if(dtest(2) .and. abs(chisq) > 10.d0 .and. i == numiter) write(*,*) 'Warning: Noise fit did not converge properly. Estimates may be OK anyway'
!       else
!          f_knee = exp(p(0))
!          alpha  = -exp(p(1))
!          sigma0 = sigma1/sqrt(get_bias(f_knee, alpha))
!          exit
!       end if
!    end do
!    deallocate(p)
!
!    ! Add a safety for extremely low f_knee's, to avoid NaN's in later programs
!    if (f_knee < 1.d-10 .or. alpha > -0.01d0) then
!       ! Return white noise
!       f_knee = 1.d-10
!       alpha  = -10.d0
!    end if
!
!    if (present(chisq_out)) then
!       ! Check goodness-of-fit
!       ind1   = nint(0.2d0 / dnu_wn)
!       ind2   = nint(nu_high_wn / dnu_wn)
!       accept = check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, f, chisq)
!       chisq_out = chisq
!    end if
!
!101 deallocate(f)
!    deallocate(mask, freqs)
!
!  end subroutine fit_1overf_profile
!
!  function lnL_noise_powell(p) result(lnl)
!    use healpix_types
!    implicit none
!    real(dp), dimension(:), intent(in) :: p
!    real(dp)  :: alpha, f_knee, sigma0, lnl
!    real(dp)  :: x, xa, xa1, m_tot, m_fxa1, m_lxa1, s_tot, s_xa, bias
!    integer(i4b) :: i
!    lnl    = 1d30
!    f_knee =  exp(p(1))
!    alpha  = -exp(p(2))
!    if (f_knee > 3.d0 .or. alpha < -10.d0) return
!    m_tot = 0; m_fxa1 = 0; m_lxa1 = 0; s_tot = 0; s_xa = 0
!    do i = 1, max(ind_max, srange(2))
!       if(mask(i) == 0) cycle
!       x      = freqs(i)/f_knee
!       xa     = x**alpha
!       xa1    = xa+1
!       if(i <= ind_max) then
!          m_tot  = m_tot   + mask(i)
!          m_fxa1 = m_fxa1  + mask(i) * f(i)/xa1
!          m_lxa1 = m_lxa1  + mask(i) * log(xa1)
!       end if
!       if(i >= srange(1) .and. i <= srange(2)) then
!          s_tot   = s_tot   + mask(i)
!          s_xa    = s_xa    + mask(i) * xa
!       end if
!    end do
!    bias   = 1/(1+s_xa/s_tot)
!    sigma0 = sigma1*sqrt(bias)
!    lnl    = (m_fxa1/sigma0**2 + m_lxa1)/m_tot + log(sigma0**2)
!  end function
!
!  function dlnL_noise_powell(p) result(dlnl)
!    use healpix_types
!    implicit none
!    real(dp), dimension(:), intent(in) :: p
!    real(dp)  ::alpha, f_knee, sigma0, dlnl(size(p))
!    real(dp)  :: x, xa, xa1, dxa1(2), m_tot, m_fxa1, m_lxa1, s_tot, s_xa, s_xalx
!    real(dp)  :: bias, dbias(2), dm_fxa1(2), dm_lxa1(2)
!    integer(i4b) :: i
!
!    f_knee =  exp(p(1))
!    alpha  = -exp(p(2))
!    ! Precompute the terms needed. The variable names here look pretty
!    ! nasty, but they try to describe what goes into that variable.
!    m_tot = 0; m_fxa1 = 0; m_lxa1 = 0; s_tot = 0; s_xa = 0; s_xalx = 0
!    dm_fxa1 = 0; dm_lxa1 = 0
!    do i = 1, max(ind_max, srange(2))
!       if(mask(i) == 0) cycle
!       x      = freqs(i)/f_knee
!       xa     = x**alpha
!       xa1    = xa+1
!       if(i <= ind_max) then
!          dxa1   = [ -alpha/f_knee, log(x) ] * xa
!          ! To build the bias and likelihood, we need various weighted sums.
!          ! m_* indicate those using the full mask while s_* use the sigma
!          ! mask.
!          m_tot   = m_tot   + mask(i)
!          m_fxa1  = m_fxa1  + mask(i) * f(i)/xa1
!          m_lxa1  = m_lxa1  + mask(i) * log(xa1)
!          dm_fxa1 = dm_fxa1 - mask(i)*f(i)/xa1**2 * dxa1
!          dm_lxa1 = dm_lxa1 + mask(i)/xa1         * dxa1
!       end if
!       if(i >= srange(1) .and. i <= srange(2)) then
!          s_tot   = s_tot   + mask(i)
!          s_xa    = s_xa    + mask(i) * xa
!          s_xalx  = s_xalx  + mask(i) * xa * log(x)
!       end if
!    end do
!    bias   = 1/(1+s_xa/s_tot)
!    dbias  = bias**2/s_tot * [ alpha/f_knee*s_xa, -s_xalx ]
!    sigma0 = sigma1*sqrt(bias)
!
!    ! Ok, with all those nasty symbols defined, we can finally
!    ! calculate the change in the likelihood. The likelihood is
!    ! L = m_fxa1/m_tot/sigma0**2 + m_lxa1/m_tot + log(sigma0**2),
!    ! which gives
!    dlnl = ((dm_fxa1 - m_fxa1*dbias)/sigma0**2 + dm_lxa1)/m_tot + dbias
!    ! And finally transform back to the original coordinates
!    dlnl   = dlnl * [ f_knee, alpha ]
!!write(*,*) "A", sigma1, f_knee, alpha, sigma0
!!write(*,*) "B", dlnl
!  end function
!
!  function get_bias(f_knee, alpha) result(bias)
!    use healpix_types
!    implicit none
!    real(dp)  :: alpha, f_knee, xa, s_tot, s_xa, bias
!    integer(i4b) :: i
!    s_tot = 0; s_xa = 0
!    do i = srange(1), srange(2)
!       if(mask(i) == 0) cycle
!       xa     = (freqs(i)/f_knee)**alpha
!       s_tot   = s_tot   + mask(i)
!       s_xa    = s_xa    + mask(i) * xa
!    end do
!    bias   = 1/(1+s_xa/s_tot)
!  end function


end module quiet_noise_estimation_mod
