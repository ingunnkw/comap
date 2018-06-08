module comap_noise_estimation_mod
  use healpix_types
  use quiet_fft_mod
  use math_tools
  use quiet_utils
  use quiet_mpi_mod
  use quasi_newton_mod
  implicit none

  real(dp),                            private :: samprate_int
  integer(i4b),                        private :: ind_max
  real(dp), allocatable, dimension(:), private :: f, freqs
  real(dp), allocatable, dimension(:), private :: mask

  integer(i4b), private :: FIT_ALPHA  = 1
  integer(i4b), private :: FIT_F_KNEE = 2

  integer(i4b), parameter, private :: npar = 2
  real(dp),                private :: nu_low_wn, nu_high_wn
  real(dp),                private :: prior(npar,2), f_scan_int(2), df_scan_int
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
       & snum, frequency, detector, chisq_out, apply_scanmask, refit, limits)
    implicit none

    real(dp),                intent(in)            :: samprate, f_scan(2), df_scan
    integer(i4b),            intent(in),  optional :: snum, frequency, detector, limits(2) 
    real(dp), dimension(1:),              optional :: tod, tod_ps
    real(dp),                intent(inout)         :: sigma0, alpha, f_knee 
    real(dp),                intent(out), optional :: chisq_out
    logical(lgt),            intent(in),  optional :: apply_scanmask, refit 

    integer(i4b) :: i, j, n, m, numsamp, numbin, ind1, ind2, ierr, iter, numiter
    logical(lgt) :: accept, est_sigma0_from_low_freq, use_grid, apply_scanmask_, refit_ ! TMR
    real(dp)     :: alpha_min, alpha_max, f_knee_min, f_knee_max, dalpha, df_knee, nu, P_nu, &
         & dnu_wn, fret, gtol, sigma1
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
       write(*,*) 'comap_noise_estimation_mod: Error -- either tod or fft must be present'
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
       if (f_scan_int(1) > 0.d0) then
          ind1   = 0
          ind2   = 0
          nu     = 0.d0
          do while (ind2 < n-1)
             nu     = nu + f_scan_int(1)
             ind1   = max(nint((nu-df_scan_int) / dnu_wn),1)
             ind2   = min(nint((nu+df_scan_int) / dnu_wn),n-1)
             mask(ind1:ind2) = 0.d0
          end do
       end if
       if (f_scan_int(2) > 0.d0) then
          ind1   = 0
          ind2   = 0
          nu     = 0.d0
          do while (ind2 < n-1)
             nu     = nu + f_scan_int(2)
             ind1   = max(nint((nu-df_scan_int) / dnu_wn),1)
             ind2   = min(nint((nu+df_scan_int) / dnu_wn),n-1)
             mask(ind1:ind2) = 0.d0
          end do
       end if
    end if
    mask(0)   = 0.d0
    mask(n-1) = 0.d0

    if (present(limits)) then
       mask(0:limits(1)) = 0.d0
       mask(limits(2):n-1) = 0.d0
    end if

!!$    open(58,file='powspec_noise.dat')
!!$    do i = 1, n
!!$       write(58,*) ind2freq(i, samprate, n), f(i-1), mask(i-1)
!!$    End do
!!$    close(58)
!!$    call mpi_finalize(i)
!!$    stop


    allocate(p(0:npar-1))

    if (.not. refit_) then

       ! This is the measured sigma0 in the signal region under the assumption
       ! that we have no 1/f component. 
       ind1    = max(nint(nu_low_wn / dnu_wn),2)
       ind2    = min(nint(nu_high_wn / dnu_wn),n-2)
       sigma1  = sqrt(sum(f(ind1:ind2)*mask(ind1:ind2)) / sum(mask(ind1:ind2)))
       ind_max = n-1

       p_full(0) = log(sigma1**2)    ! variance
       p_full(1) = log(1.d0) ! f_knee
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
       call dfpmin(p_full,gtol,iter,fret,lnL_noise_powell_full,dlnL_noise_powell_full, ierr=ierr)

       if (ierr == 1) then
          if (present(snum)) then
             write(*,*) 'Error: NaNs in noise fit -- {snum,freq,det} = ', snum, frequency, detector
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

!!$    open(57,file='pow.dat')
!!$    do i = 1, n-1
!!$       write(57,*) freqs(i), f(i), sigma0**2 * (1 + (freqs(i)/f_knee)**alpha)
!!$    end do
!!$    close(57)
!!$    call mpi_finalize(i)
!!$    stop

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

  function lnL_noise_powell_full(p)
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp)             :: lnL_noise_powell_full

    integer(i4b) :: i
    real(dp)     :: P_nu, alpha, f_knee, x, sigma_sq

    sigma_sq =  exp(min(p(1),200.d0))
    if (p(2) < -10.d0) then
       f_knee = exp(-10.d0)
    else if (p(2) > 3.d0) then
       f_knee = exp(3.d0)
    else
       f_knee   =  exp(p(2))
    end if
    if (p(3) > 6.d0) then ! Avoid overflow
       alpha = -exp(6.d0)
    else
       alpha    = -exp(p(3))
    end if

    if (f_knee > 10.d0 .or. alpha < -10.d0 .or. sigma_sq == 0.d0) then
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

  end function lnL_noise_powell_full

  function dlnL_noise_powell_full(p) 
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp), dimension(size(p))        :: dlnL_noise_powell_full

    integer(i4b) :: i
    real(dp)     :: nu, P_nu, alpha, f_knee, x, dPda, dPdf, current_sigma0, &
         & dS2da, dS2df, sigma_sq, dLdP

    sigma_sq =  exp(min(p(1),200.d0))
    if (p(2) < -10.d0) then
       f_knee = exp(-10.d0)
    else if (p(2) > 3.d0) then
       f_knee = exp(3.d0)
    else
       f_knee   =  exp(p(2))
    end if
    if (p(3) > 6.d0) then ! Avoid overflow
       alpha = -exp(6.d0)
    else
       alpha    = -exp(p(3))
    end if

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
  end function dlnL_noise_powell_full


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


end module comap_noise_estimation_mod
