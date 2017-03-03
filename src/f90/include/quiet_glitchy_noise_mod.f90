module quiet_glitchy_noise_mod
  use healpix_types
  use rngmod
  use quiet_utils
  use quiet_fft_mod
  use quiet_constrained_mod
  use quasi_newton_mod
  implicit none

  type glitch
     real(dp),     dimension(:),   allocatable :: shape
     integer(i4b) :: begtime, length
  end type glitch

  type(planck_rng), private :: rng
  integer(i4b),     private :: cg_seed
  integer(i4b),     private :: myid_int
  logical(lgt),     private :: dosample
  integer(i4b),             parameter, private :: SIGMA0_IND = 1
  integer(i4b),             parameter, private :: ALPHA_IND  = 2
  integer(i4b),             parameter, private :: FKNEE_IND  = 3
  real(dp),                            private :: samprate 
  real(sp), allocatable, dimension(:), private :: powspec, tod

contains

  subroutine initialize_glitchy_noise_mod(parfile, seed_in)
    implicit none

    character(len=*), intent(in), optional :: parfile
    integer(i4b),     intent(in), optional :: seed_in
    integer(i4b) :: seed

    if (present(seed_in)) then
       seed = seed_in
    else
       call get_parameter(0, parfile, 'BASE_SEED', par_int=seed)
    end if
    call rand_init(rng, seed)
  end subroutine initialize_glitchy_noise_mod

  ! Driver routine -- this is what external users should call
  subroutine fit_1overf_profile_with_glitches(samprate_in, tod_in, mask, sigma0, alpha, fknee, &
       & numit, myid, todid, dosample_in, var, skew, templates, constrained_realization)
    implicit none
    real(dp),                   intent(in)            :: samprate_in
    real(sp),     dimension(:), intent(in)            :: tod_in
    logical(lgt), dimension(:), intent(in)            :: mask
    real(dp),                   intent(inout)         :: sigma0, alpha, fknee
    integer(i4b),               intent(out)           :: numit
    integer(i4b),               intent(in)            :: myid, todid
    logical(lgt),               intent(in)            :: dosample_in
    real(dp),                   intent(out), optional :: var(3), skew(3)
    real(sp),     dimension(:), intent(in),  optional :: templates
    real(sp),     dimension(:), intent(out), optional :: constrained_realization
    integer(i4b) :: i, j, ntod, nfreq, iter, numgib, burnin
    logical(lgt) :: converged
    real(dp)     :: lim, tamp
    character(len=3)            :: filnum
    real(dp),     allocatable, dimension(:,:) :: samples
    type(glitch), allocatable, dimension(:)   :: glitches

    dosample = dosample_in
    samprate = samprate_in
    lim = 0.0001d0
    cg_seed  = nint(rand_uni(rng)*1000000.d0)

    ! Check that the module isn't dead
    if (is_nan(real(sum(abs(tod_in)),dp))) then
       sigma0 = 0.d0
       alpha  = 0.d0
       fknee = 0.d0
       return
    end if

  ! Set up spectrum structure
    ntod = size(tod_in)
    nfreq = (ntod/2)+1
    allocate(powspec(nfreq), tod(ntod))
    tod = tod_in
  
!if (todid==538) then
!open(13, file='tod_before.dat')
!do i = 1, ntod
!   write(13,*) tod(i)
!end do
!close(13)
!end if

    ! Read glitch shape info
!    allocate(glitches(1))
!    glitches(1)%begtime = 10000
!    glitches(1)%length  = 200
!    allocate(glitches(1)%shape(glitches(1)%length))
!    do i = 1, glitches(1)%length
!       glitches(1)%shape(i) = exp(-real(i, dp)/25.d0)*1.1d-4
!    end do
!open(13, file='glitch.dat')
!do i = 1, glitches(1)%length
!   write(13,*) glitches(1)%shape(i)
!end do
!close(13)

!tod(glitches(1)%begtime:glitches(1)%begtime + glitches(1)%length -1) = &
!     & tod(glitches(1)%begtime:glitches(1)%begtime + glitches(1)%length -1) + glitches(1)%shape
!do i = 1, glitches(1)%length
!   glitches(1)%shape(i) = exp(-real(i, dp)/25.d0)*1d-4
!end do
!open(13, file='tod_after.dat')
!do i = 1, ntod
!   write(13,*) tod(i)
!end do
!close(13)
!stop

    ! Initialize chain at input values
    numgib = 100
    burnin = 10
    allocate(samples(3,numgib))
    samples(SIGMA0_IND,1) = sigma0
    samples(ALPHA_IND, 1) = alpha
    samples(FKNEE_IND, 1) = fknee
    numit = -10
    tamp = 0.d0
    ! Do the Gibbs sampling
    do i = 1, numgib-1
!       call sample_glitch_amplitudes(samples(SIGMA0_IND,i), samples(ALPHA_IND,i), samples(FKNEE_IND,i))
       call sample_noise_spectrum_given_params(mask, samples(SIGMA0_IND,i), &
        & samples(ALPHA_IND,i), samples(FKNEE_IND,i), i==1)
       call sample_template_amplitude (templates, mask, samples(SIGMA0_IND,i), samples(ALPHA_IND,i), samples(FKNEE_IND,i), tamp)
       samples(:,i+1) = samples(:,i)
       call sample_noise_params_given_spectrum(samples(SIGMA0_IND,i+1), samples(ALPHA_IND,i+1), samples(FKNEE_IND,i+1), & 
            & i, myid, todid)
       if (.not. dosample .and. i>1) then
          if (all(abs(samples(:,i+1)-samples(:,i))/abs(samples(:,i)) < lim) .and. &
               & all(abs(samples(:,i)-samples(:,i-1))/abs(samples(:,i)) < lim)) exit
       end if
    end do
    numgib = i
    numit = numgib

    ! Return mean of second half of chain as estimate
    if (dosample) then
!       sigma0 = sum(samples(SIGMA0_IND,numgib/2+1:numgib)) / real(numgib/2,dp)
!       alpha  = sum(samples(ALPHA_IND, numgib/2+1:numgib)) / real(numgib/2,dp)
!       fknee  = sum(samples(FKNEE_IND, numgib/2+1:numgib)) / real(numgib/2,dp)
!       if (present(var)) then
!          var(1) = sqrt(sum((samples(SIGMA0_IND,numgib/2+1:numgib)-sigma0)**2) / real(numgib/2-1,dp))
!          var(2) = sqrt(sum((samples(ALPHA_IND,numgib/2+1:numgib)-alpha)**2) / real(numgib/2-1,dp))
!          var(3) = sqrt(sum((samples(FKNEE_IND,numgib/2+1:numgib)-fknee)**2) / real(numgib/2-1,dp))
!       end if
       sigma0 = sum(samples(SIGMA0_IND,burnin+1:numgib)) / real(numgib-burnin,dp)
       alpha  = sum(samples(ALPHA_IND, burnin+1:numgib)) / real(numgib-burnin,dp)
       fknee  = sum(samples(FKNEE_IND, burnin+1:numgib)) / real(numgib-burnin,dp)
       if (present(var)) then
          var(1) = sqrt(sum((samples(SIGMA0_IND,burnin+1:numgib)-sigma0)**2) / real(numgib-burnin-1,dp))
          var(2) = sqrt(sum((samples(ALPHA_IND,burnin+1:numgib)-alpha)**2) / real(numgib-burnin-1,dp))
          var(3) = sqrt(sum((samples(FKNEE_IND,burnin+1:numgib)-fknee)**2) / real(numgib-burnin-1,dp))
       end if
       if (present(skew)) then
          skew(1) = (sum((samples(SIGMA0_IND,burnin+1:numgib)-sigma0)**3) / real(numgib-burnin-1,dp))/var(1)**3
          skew(2) = (sum((samples(ALPHA_IND,burnin+1:numgib)-alpha)**3) / real(numgib-burnin-1,dp))/var(2)**3
          skew(3) = (sum((samples(FKNEE_IND,burnin+1:numgib)-fknee)**3) / real(numgib-burnin-1,dp))/var(3)**3
       end if
    else
       sigma0 = samples(SIGMA0_IND, numgib)
       alpha  = samples(ALPHA_IND, numgib)
       fknee  = samples(FKNEE_IND, numgib)
    end if

    if (present(constrained_realization)) then
       ! Generate a constrained realization based on best-fit parameters
    end if

    deallocate(samples, powspec, tod)!glitches

  end subroutine fit_1overf_profile_with_glitches

  ! Sigurd's routine
  subroutine sample_noise_spectrum_given_params(mask, sigma0, alpha, fknee, firstcall)
    implicit none
    logical(lgt),               intent(in)  :: firstcall
    logical(lgt), dimension(:), intent(in)  :: mask
    real(dp),                   intent(in)  :: sigma0, alpha, fknee
    real(dp),     dimension(:), allocatable :: tod2
    integer(i4b)                            :: i, j, m, n
    type(planck_rng)                        :: rng2
    integer(i4b), save                      :: myseed

    call rand_init(rng2, cg_seed)
    n = size(tod); m = n/2+1
    powspec(1) = 0
    do i = 2, m
       powspec(i) = sigma0**2*(1+(ind2freq(i, samprate, m)/fknee)**alpha)
    end do
    allocate(tod2(n))!,ft(m))!,tmp(n))
    if(dosample) then
       call constrained_realization(real(powspec,dp), real(tod,dp), mask, rng, tod2)
    else
       call constrained_realization(real(powspec,dp), real(tod,dp), mask, rng2, tod2)
    end if
    tod = tod2
    deallocate(tod2)
  end subroutine sample_noise_spectrum_given_params

  ! Ingunn
  subroutine sample_template_amplitude (templates, mask, sigma0, alpha, fknee, tamplitude)
    implicit none
    real(sp),     dimension(:)              :: templates
    logical(lgt), dimension(:), intent(in)  :: mask
    real(dp),                   intent(in)  :: sigma0, alpha, fknee
    real(dp)                                :: tamplitude
    complex(spc), allocatable :: ftemp(:), ftod(:)
    real(sp),     allocatable :: cov(:)
    real(dp)                  :: teller, nevner
    integer(i4b)              :: i, n, m, j
    n = size(tod); m = n/2+1
    ! Add old template estimate to tod
    tod = tod + tamplitude * templates
    ! Estimate new template amplitude
    allocate(ftemp(m), ftod(m), cov(m))
    cov(1) = 0
    do i = 2, m
       cov(i) = sigma0**2*(1+(ind2freq(i, samprate, m)/fknee)**alpha)
    end do
    call fft(tod,ftod,1)
    call fft(templates,ftemp,1)
    teller = 0; nevner = 0
    do i = 2, m
       teller = teller + conjg(ftemp(i))*ftod(i)/cov(i)
       nevner = nevner + conjg(ftemp(i))*ftemp(i)/cov(i)
    end do
    if (dosample) then
       tamplitude = teller/nevner + rand_gauss(rng)/sqrt(nevner)
    else
       tamplitude = teller/nevner
    end if
    ! Subtract new template estimat from tod
    tod = tod - tamplitude * templates
  end subroutine sample_template_amplitude

  ! Ingunn's routine
  subroutine sample_noise_params_given_spectrum(sigma0, alpha, fknee, stepid, myid, todid)
    implicit none
    real(dp),                 intent(inout) :: sigma0, alpha, fknee
    integer(i4b),             intent(in)    :: stepid, myid, todid
    real(dp),   dimension(:,:), allocatable :: params, corr
    complex(spc), dimension(:), allocatable :: ft
    real(dp)                    :: t, tvec(3), loglike1, loglike2, fret, mean(3), var(3), corrlimit, bpar(3)
    integer(i4b)                :: n, i, j, s, iter, ierr, accept, maxcorr, mark 
    character(len=3)            :: filnum
    real(dp), save              :: smatrix(3,3)
    integer(i4b), save          :: maxit
    logical(lgt), save          :: firstcall = .true.  

    myid_int = myid
    mark = 0
    allocate(ft(size(powspec)))
    call fft(tod,ft,1)
    call extract_powspec(ft, powspec)
    
    ! Set up arrays of param points
    accept = 0
    if (firstcall) then
       maxit = 10000
    end if
    allocate(params(3,0:maxit))
    params(1,0) = sigma0
    params(2,0) = alpha
    params(3,0) = fknee
    ! If not gibbssampling (or for first for gibbs step) find max-like point
    if (.not. dosample .or. stepid==1) then
       ! Find maxlike point
       bpar = params(:,0)
       call dfpmin(params(:,0), 1d-6, iter,fret,minusloglike,dml,ierr)
       if(ierr/=0 .or. iter==1) then
          do j = 1, 20
             params(:,0)=bpar/2*(1+real(j,dp)/10)
             call dfpmin(params(:,0), 1d-6, iter,fret,minusloglike,dml,ierr)
!             write(*,fmt='(i4,a,i3,3e16.6,2i3)') myid, ' ml', j, params(:,0), ierr, iter
             if (ierr ==0 .and. iter>1 )exit
          end do
       end if
       if (ierr /= 0 .or. iter==1) then
          write(*,fmt='(i4,a,i3,4e16.6,i2)') myid, ' A1', j, params(:,0), fret, ierr
          write(*,fmt='(i4,a,i3,a,3e12.4)') myid, ' A2', iter, ' Best params maxlike search:', params(:,0)
          write(*,fmt='(i4,a)') myid, ' A3 Quitting'
          call mpi_finalize
          stop
       end if
    end if

    ! If not gibbs sampling then finished
    if (.not. dosample) then
       ! Output params
       sigma0 = abs(params(1,0))
       alpha  = params(2,0)
       fknee  = params(3,0)
    ! If gibbs sampling then continue   
    else 
       ! Set optimal step length
       if (stepid==1) call calc_step_matrix(params(1,0), params(2,0), params(3,0), smatrix, myid, todid)
       ! Calculate likelihood for starting point
       loglike1 = loglike(params(:,0))
       ! Loop over all iterations
       do s = 1, maxit
          ! Choose new point
          do i = 1, 3
             tvec(i) = rand_gauss(rng)
          end do
          params(:, s) = params(:, s-1) + matmul(smatrix,tvec)
          ! Calculate new likelihood
          loglike2 = loglike(params(:,s))
          ! Metropolis-Hastings desicion
          t = rand_uni(rng)
          if (exp(loglike2-loglike1) > t) then 
             loglike1 = loglike2
             accept = accept + 1
!             write(*,*) s
          else
             params(:, s) = params(:, s-1)
          end if
       end do
       ! Output params
       sigma0 = abs(params(1,maxit))
       alpha  = params(2,maxit)
       fknee  = params(3,maxit)
       ! Ouput chains
!       open(42, file="parchains.txt")
!       do i = 0, maxit
!          write(42,*) params(1,i), params(2,i), params(3,i)
!       end do
!       close(42)
       if (real(accept,dp)/real(maxit,dp)<0.05 .or. real(accept,dp)/real(maxit,dp)> 0.95) &
            & write(*,*) myid, stepid, 'Accept ratio =', real(accept,dp)/real(maxit,dp), maxit
       if (firstcall) then
          ! Calculate correlation length
          maxcorr = 100
          allocate(corr(1:3, 0:maxcorr))
          corr = 0.d0
          mean = 0.d0
          var = 0.d0
          do s = 1, maxit
             mean(:) = mean(:) + params(:, s)
          end do
          mean = mean/real(maxit,dp)
          do s = 1, maxit
             var(:) = var(:) + (params(:,s)-mean(:))**2
          end do
          var = var/real(maxit,dp)
          do i = 0, maxcorr
             do s = 1, maxit-i
                corr(:,i) = corr(:,i) + (params(:,s)-mean(:))*(params(:,s+i)-mean(:))
             end do
             corr(:,i) = corr(:,i)/real(maxit-i,dp)
          end do
          do j = 1, 3
             corr(j,:) = corr(j,:)/var(j)
          end do
          ! Search for no more than 10% correlation
          corrlimit = 0.2d0
          do i = 1, maxcorr
             if (all(corr(:,i) < corrlimit)) exit
          end do
          if (i == maxcorr) then
             write(*,*) 'Too much correlations. quiting'
             stop
          else
             write(*,fmt='(i4,i4,a,f5.2)') myid,  i, ' iterations needed to get correlations below', corrlimit
             maxit = i*2
          end if
          ! Ouput correlation length
          call int2string(myid, filnum)
          open(42, file="corrlength_sigma"//filnum//".txt")
          do i = 0, maxcorr
             write(42,*) corr(1,i)
          end do
          close(42)
          open(42, file="corrlength_alpha"//filnum//".txt")
          do i = 0, maxcorr
             write(42,*) corr(2,i)
          end do
          close(42)
          open(42, file="corrlength_fknee"//filnum//".txt")
          do i = 0, maxcorr
             write(42,*) corr(3,i)
          end do
          close(42)
          deallocate(corr)
          firstcall = .false.
       end if
    end if
    deallocate(params)
  end subroutine sample_noise_params_given_spectrum

  function loglike(params)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in) :: params
    real(dp)                           :: loglike
    integer(i4b)                       :: n, i 
    real(dp)                           :: sigma0, alpha, fknee, nu, pnu
    sigma0 = params(1)
    alpha  = params(2)
    fknee  = params(3)
    if (fknee < 0.d0 .or. sigma0 < 0.d0) then
       loglike = -1.d30
    else
       n = size(powspec)
       loglike = 0.d0
       do i = 2, n
          nu   = ind2freq(i, samprate, n)
          pnu  = sigma0**2 * (1+(nu/fknee)**alpha)
          loglike = loglike - powspec(i)/pnu - log(pnu)
       end do
    end if
  end function loglike
  
  function minusloglike(params)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in) :: params
    real(dp)                           :: minusloglike
    minusloglike = -loglike(params)
  end function minusloglike
  
  function dml(params)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in) :: params
    real(dp), dimension(size(params))  :: dml
    integer(i4b)                       :: n, i 
    real(dp)                           :: sigma0, alpha, fknee, nu, p, b, dpp(3)
    sigma0 = params(1)
    alpha  = params(2)
    fknee  = params(3)
    n = size(powspec)
    dml = 0.d0
    do i = 2, n
       nu   = ind2freq(i, samprate, n)
       p    = sigma0**2 * (1+(nu/fknee)**alpha)
       b    = (1-powspec(i)/p)/p
       dpp(1) = 2*sigma0*(1+(nu/fknee)**alpha)
       dpp(2) = sigma0**2*(nu/fknee)**alpha*log(nu/fknee)
       dpp(3) = -sigma0**2*alpha/fknee*(nu/fknee)**(alpha)
       dml = dml +b*dpp
    end do
  end function dml

  subroutine calc_step_matrix(sigma0, alpha, fknee, matrix, myid, todid)
    implicit none

    real(dp),               intent(in)    :: sigma0, alpha, fknee
    real(dp),               intent(out)   :: matrix(3,3)
    integer(i4b), optional, intent(in)    :: myid, todid

    integer(i4b)    :: n, i, status
    real(dp)        :: nu, p, a, b, dpp(3), mat(3,3), eigenvectors(3,3), eigenvals(3)

    n = size(powspec)
    matrix = 0.d0
    do i = 2, n
       nu    = ind2freq(i, samprate, n)
       p     = sigma0**2 * (1+(nu/fknee)**alpha)
       a     = (2*powspec(i)/p-1)/(p**2)
       b     = (1-powspec(i)/p)/p
       dpp(1) = 2*sigma0*(1+(nu/fknee)**alpha)
       dpp(2) = sigma0**2*(nu/fknee)**alpha*log(nu/fknee)
       dpp(3) = -sigma0**2*alpha/fknee*(nu/fknee)**alpha
       matrix(1,1) = matrix(1,1) + a*dpp(1)*dpp(1) +b*2*(1+(nu/fknee)**alpha)
       matrix(1,2) = matrix(1,2) + a*dpp(1)*dpp(2) +b*2*sigma0*(nu/fknee)**alpha*log(nu/fknee)
       matrix(1,3) = matrix(1,3) + a*dpp(1)*dpp(3) -b*2*sigma0*alpha/fknee*(nu/fknee)**alpha
       matrix(2,1) = matrix(1,2)
       matrix(2,2) = matrix(2,2) + a*dpp(2)*dpp(2) +b*sigma0**2*(nu/fknee)**alpha*(log(nu/fknee))**2
       matrix(2,3) = matrix(2,3) + a*dpp(2)*dpp(3) -b*sigma0**2/fknee*(nu/fknee)**alpha*(1+alpha*log(nu/fknee))
       matrix(3,1) = matrix(1,3)
       matrix(3,2) = matrix(2,3)
       matrix(3,3) = matrix(3,3) + a*dpp(3)*dpp(3) +b*sigma0**2*alpha*(alpha+1)/fknee**2*(nu/fknee)**alpha
    end do

    call invert_singular_matrix(matrix, 0.d0)
    call cholesky_decompose(matrix, mat, status)
    matrix=mat
    if (status /= 0) then
       write(*,fmt='(i4,a,i3,3e16.6)') myid, ' stat =', status, sigma0, alpha, fknee
       !call mpi_finalize
       stop
    end if

!    call get_eigen_decomposition(myid, matrix, eigenvals, eigenvectors) 
!    !call invert_eigenvals(myid, eigenvals, unit, outprefix, outtext)
!    mat = transpose(eigenvectors)
!    do i = 1, 3
!       mat(i,:) = mat(i,:) / sqrt(eigenvals(i))
!    end do
!    call DGEMM('N', 'N', 3, 3, 3, 1.d0, eigenvectors, 3, mat, 3, 0.d0, matrix, 3)

  end subroutine calc_step_matrix

end module quiet_glitchy_noise_mod
