module comap_noise_estimation_mod
  use healpix_types
  use quiet_fft_mod
  use math_tools
  use quiet_utils
  use quiet_mpi_mod
  use quasi_newton_mod
  implicit none

  real(dp),                            private :: samprate_int
  integer(i4b),                        private :: ind_min, ind_max
  real(dp), allocatable, dimension(:), private :: f, freqs
  real(dp), allocatable, dimension(:), private :: mask

  integer(i4b), private :: FIT_ALPHA  = 1
  integer(i4b), private :: FIT_F_KNEE = 2
  real(dp),     private :: DEFAULT_ALPHA = -1.d0
  real(dp),     private :: DEFAULT_FKNEE =  1.d0

  logical(lgt) :: fit_param(2)

  real(dp),                private :: nu_low_wn, nu_high_wn
  real(dp),                private :: f_scan_int(2), df_scan_int
  logical(lgt),            private :: initialized = .false.

contains

  subroutine initialize_noise_estimation_mod(parfile)
    implicit none
    character(len=*) :: parfile

    if(initialized) return
    call get_parameter(0, parfile, 'NOISE_EST_LOW_NU', par_dp=nu_low_wn)
    call get_parameter(0, parfile, 'NOISE_EST_HIGH_NU', par_dp=nu_high_wn)

    initialized = .true.
  end subroutine initialize_noise_estimation_mod

  subroutine fit_1overf_profile(samprate, f_scan, df_scan, sigma0, alpha, f_knee, tod, tod_ps, &
       & snum, frequency, detector, chisq_out, apply_scanmask, refit, limits, fit_par)
    implicit none

    real(dp),                intent(in)            :: samprate, f_scan(2), df_scan
    real(dp),                intent(in),  optional :: limits(2)
    integer(i4b),            intent(in),  optional :: snum, frequency, detector
    logical(lgt),            intent(in),  optional :: fit_par(2)
    real(dp), dimension(1:),              optional :: tod, tod_ps
    real(dp),                intent(inout)         :: sigma0, alpha, f_knee 
    real(dp),                intent(out), optional :: chisq_out
    logical(lgt),            intent(in),  optional :: apply_scanmask, refit 

    integer(i4b) :: i, j, n, m, numsamp, numbin, ind1, ind2, ierr, iter, numiter, npar
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
    fit_param = .true.;       if(present(fit_par)) fit_param = fit_par

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
       if (.false. .and. f_scan_int(2) > 0.d0) then
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
       ind_min             = freq2ind(limits(1), samprate, n)
       ind_max             = freq2ind(limits(2), samprate, n)
       mask(0:ind_min-1)   = 0.d0
       mask(ind_max+1:n-1) = 0.d0
    else
       ind_min = 1
       ind_max = n-1 
    end if

!!$    open(57,file='pow.dat')
!!$    do i = 1, n-1
!!$       write(57,*) freqs(i), f(i)!, sigma0**2 * (1 + (freqs(i)/f_knee)**alpha)
!!$    end do
!!$    close(57)
!!$    call mpi_finalize(i)
!!$    stop

!!$    open(58,file='powspec_noise.dat')
!!$    do i = 1, n
!!$       write(58,*) ind2freq(i, samprate, n), f(i-1), mask(i-1)
!!$    End do
!!$    close(58)
!!$    call mpi_finalize(i)
!!$    stop

    npar = count(fit_param)
    allocate(p(npar))

    ! Compute white noise solution
    ind1    = max(nint(nu_low_wn / dnu_wn),2)
    ind2    = min(nint(nu_high_wn / dnu_wn),n-2)
    if (all(.not. fit_param)) then
       ! Return white noise only
       alpha   = -2.d0
       f_knee  = 1d-6
       sigma0  = sqrt(compute_sigma_sq(alpha, f_knee))
       deallocate(p,f,mask,freqs)
       return
    end if

    ! Initialize search parameters
    i = 1
    if (fit_param(1)) then
       p(i) = log(DEFAULT_FKNEE) ! f_knee
       i    = i+1
    end if
    if (fit_param(2)) then
       p(i) = log(-DEFAULT_ALPHA) ! alpha
       i    = i+1
    end if

    gtol = 1.d-5
    numiter = 5
    do i = 1, numiter
       call dfpmin(p,gtol,iter,fret,lnL_noise_powell_full,dlnL_noise_powell_full, ierr=ierr)
       
       if (ierr == 1) then
          if (present(snum)) then
             write(*,*) 'Error: NaNs in noise fit -- {snum,freq,det} = ', snum, frequency, detector
          else
             write(*,*) 'Error: NaNs in noise fit -- this should be fixed'
          end if

          sigma0 = 0; alpha = 0; f_knee = 0
          goto 101                                    ! This aborts the noise estimation
          stop
       else 
          j = 1
          if (fit_param(1)) then
             f_knee = exp(p(j))
             j      = j+1
          else
             f_knee = 1e-6
          end if
          if (fit_param(2)) then
             alpha  = -exp(p(j))
          else
             alpha = DEFAULT_ALPHA
          end if
          sigma0 = sqrt(compute_sigma_sq(alpha, f_knee))

          ! Check goodness-of-fit between 0.03 and 3 Hz
          ind1   = max(nint(0.03d0 / dnu_wn),2)
          ind2   = nint(3.00d0 / dnu_wn)
          accept = check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, f, chisq)
          if(dtest(2) .and. abs(chisq) > 10.d0 .and. i == numiter) &
               & write(*,*) 'Warning: Noise fit did not converge properly. Estimates may be OK anyway'
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
       ind1   = max(nint(0.03d0 / dnu_wn),2)
       ind2   = nint(3.00d0 / dnu_wn)
       accept = check_noise_fit(sigma0, alpha, f_knee, ind1, ind2, samprate, f, chisq)
       chisq_out = chisq
    end if

101 deallocate(f)
    deallocate(mask, freqs)

  end subroutine fit_1overf_profile

  function compute_sigma_sq(alpha, f_knee)
    implicit none
    real(dp), intent(in) :: alpha, f_knee
    real(dp)             :: compute_sigma_sq
    
    integer(i4b) :: i
    real(dp)     :: A, b, P_nu

    A = 0.d0
    b = 0.d0
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       P_nu  = 1.d0 + (freqs(i)/f_knee)**alpha
       A     = A + P_nu * f(i)
       b     = b + P_nu * P_nu
    end do

    compute_sigma_sq = A/b

  end function compute_sigma_sq

  function compute_dsigma0_dalpha(alpha, f_knee)
    implicit none
    real(dp), intent(in) :: alpha, f_knee
    real(dp)             :: compute_dsigma0_dalpha
    
    integer(i4b) :: i
    real(dp)     :: A, b, P_nu, dAda, dbda

    A    = 0.d0
    b    = 0.d0
    dAda = 0.d0
    dbda = 0.d0
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       P_nu  = 1.d0 + (freqs(i)/f_knee)**alpha
       A     = A + P_nu * f(i)
       b     = b + P_nu * P_nu

       dAda = dAda + f(i) * (freqs(i)/f_knee)**alpha * log(freqs(i)/f_knee)
       dbda = dbda + 2.d0 * P_nu * (freqs(i)/f_knee)**alpha * log(freqs(i)/f_knee)
    end do

    compute_dsigma0_dalpha = (dAda * b - A * dbda) / b**2

  end function compute_dsigma0_dalpha

  function compute_dsigma0_dfknee(alpha, f_knee)
    implicit none
    real(dp), intent(in) :: alpha, f_knee
    real(dp)             :: compute_dsigma0_dfknee
    
    integer(i4b) :: i
    real(dp)     :: A, b, P_nu, dAda, dbda

    A    = 0.d0
    b    = 0.d0
    dAda = 0.d0
    dbda = 0.d0
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       P_nu  = 1.d0 + (freqs(i)/f_knee)**alpha
       A     = A + P_nu * f(i)
       b     = b + P_nu * P_nu

       dAda = dAda + f(i) * alpha * (freqs(i)/f_knee)**(alpha-1) * (-freqs(i)/f_knee**2)
       dbda = dbda + 2.d0 * P_nu * alpha * (freqs(i)/f_knee)**(alpha-1) * (-freqs(i)/f_knee**2)
    end do

    compute_dsigma0_dfknee = (dAda * b - A * dbda) / b**2

  end function compute_dsigma0_dfknee

  function lnL_noise_powell_full(p)
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp)             :: lnL_noise_powell_full

    integer(i4b) :: i
    real(dp)     :: P_nu, alpha, f_knee, x, sigma_sq, ssum

    i = 1
    if (fit_param(1)) then
       if (p(i) < -10d0) then
          f_knee = exp(-10.d0)
       else if (p(i) > 3.d0) then
          f_knee = exp(3.d0)
       else
          f_knee = exp(p(i))
       end if
       i      = i+1
    else
       f_knee = 1e-6
    end if
    if (fit_param(2)) then
       if (p(i) > 6.d0) then
          alpha = -exp(6.d0)
       else
          alpha  = -exp(p(i))
       end if
       i      = i+1
    else
       alpha = DEFAULT_ALPHA
    end if

    if (f_knee > 10.d0 .or. alpha < -10.d0 .or. alpha > -0.01d0) then
       lnL_noise_powell_full = 1.d30
       return
    end if
    sigma_sq = compute_sigma_sq(alpha, f_knee)

    lnL_noise_powell_full = 0.d0
    !!$OMP PARALLEL PRIVATE(i,ssum,P_nu)
    ssum = 0.d0
    !!$OMP DO SCHEDULE(guided)
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       P_nu  = sigma_sq * (1.d0 + (freqs(i)/f_knee)**alpha)
       ssum  = ssum + f(i) / P_nu + log(P_nu)
    end do
    !!$OMP END DO
    !!$OMP ATOMIC
    lnL_noise_powell_full = lnL_noise_powell_full  + ssum
    !!$OMP END PARALLEL

    ! Normalize for better numerical precision
    lnL_noise_powell_full = lnL_noise_powell_full / sum(mask(1:ind_max))

    !write(*,*) real(sigma_sq,sp), real(alpha,sp), real(f_knee,sp), real(lnL_noise_powell_full,sp)

  end function lnL_noise_powell_full

  function dlnL_noise_powell_full(p) 
    use healpix_types
    implicit none

    real(dp), dimension(:), intent(in) :: p
    real(dp), dimension(size(p))        :: dlnL_noise_powell_full

    integer(i4b) :: i, j
    real(dp)     :: nu, P_nu, alpha, f_knee, x, dPda, dPdf, current_sigma0, &
         & dS2da, dS2df, sigma_sq, dLdP, ssum(3), dsda, dsdf

    i = 1
    if (fit_param(1)) then
       if (p(i) < -10d0) then
          f_knee = exp(-10.d0)
       else if (p(i) > 3.d0) then
          f_knee = exp(3.d0)
       else
          f_knee = exp(p(i))
       end if
       i      = i+1
    else
       f_knee = 1d-6
    end if
    if (fit_param(2)) then
       if (p(i) > 6.d0) then
          alpha = -exp(6.d0)
       else
          alpha  = -exp(p(i))
       end if
       i      = i+1
    else
       alpha = DEFAULT_ALPHA
    end if
    sigma_sq = compute_sigma_sq(alpha, f_knee)

    dsda = compute_dsigma0_dalpha(alpha, f_knee)
    dsdf = compute_dsigma0_dfknee(alpha, f_knee)
          

    dlnL_noise_powell_full = 0.d0
    !!$OMP PARALLEL PRIVATE(i,ssum,nu,P_nu,dLdP)
    ssum = 0.d0
    !!$OMP DO SCHEDULE(guided)
    do i = 1, ind_max
       if(mask(i) == 0) cycle
       nu        = freqs(i)
       P_nu      = sigma_sq * (1.d0 + (nu/f_knee)**alpha)
       dLdP      = -f(i)/P_nu**2 + 1.d0 / P_nu

       j = 1
       if (fit_param(1)) then
          ! dL/dtheta = dL/dP * dP/da * da/dtheta
          ssum(j) = ssum(j) + dLdP * (sigma_sq * alpha * (nu/f_knee)**(alpha-1) * (-nu/f_knee**2) + dsdf * (1.d0+(nu/f_knee)**alpha)) * f_knee
          j = j+1
       end if

       if (fit_param(2)) then
          ssum(j) = ssum(j) + dLdP * (sigma_sq * (nu/f_knee)**alpha * log(nu/f_knee) + dsda * (1+(nu/f_knee)**alpha)) * alpha
          !ssum(j) = ssum(j) - dLdP * sigma_sq * alpha * (nu/f_knee)**alpha 
          !ssum(j) = ssum(j) - dLdP * dsda * (nu/f_knee)**alpha !* alpha
       end if
    end do
    !!$OMP END DO
    do j = 1, size(p)
       !!$OMP ATOMIC
       dlnL_noise_powell_full(j) = dlnL_noise_powell_full(j)  + ssum(j)
    end do
    !!$OMP END PARALLEL

    dlnL_noise_powell_full = dlnL_noise_powell_full / sum(mask(1:ind_max))

!    write(*,*) real(alpha,sp), real(p(1),sp), '-0.0035964393', real(dlnL_noise_powell_full,sp)
!    stop


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
