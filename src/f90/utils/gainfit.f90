! Fit the gain model m(t,d) = A(d)*(1+B(t-t0)+C*cos(E(t-t0))+D*sin(E(t-t0)))
! to a set of (d,t,v,sigma) values.
program gainfit
  use quiet_utils
  use math_tools
  use powell_mod
  use quiet_hdf_mod
  use rngmod
  implicit none

  type gain_model
     real(dp) :: t0, slope, harm_amp, freq, phase, chisq
     real(dp), allocatable :: diode_amps(:)
  end type
  type(gain_model)              :: gmod_powell
  real(dp),         allocatable :: data_powell(:,:)
  integer(i4b),     allocatable :: diode_powell(:)

  character(len=512)        :: fname, modfile, resfile, line, wfile
  real(dp),     allocatable :: data(:,:), template(:), t(:), w(:,:)
  integer(i4b), allocatable :: dimap(:)
  integer(i4b)              :: i, j, k, m, n, unit, err, nsim
  logical(lgt)              :: mcmc
  real(dp)                  :: v, gain, gdev
  type(gain_model)          :: model
  type(hdf_file)            :: hfile

  if(iargc() < 3) then
     write(*,*) "Usage: Gainfit measurements.txt model.txt residuals.txt"
     write(*,*) "The measurements file should contain columns of"
     write(*,*) " [abs diode] [mjd] [gain] [stddev]"
     write(*,*) "The resulting gain model will be output in model.txt,"
     write(*,*) "and the residuals in residuals.txt."
     stop
  end if

  call getarg(1, fname)
  call getarg(2, modfile)
  call getarg(3, resfile)
  mcmc = iargc() > 3
  nsim = 1000

  unit = getlun()
  open(unit, file=fname, action="read", status="old")
  n = 0
  do
     read(unit,'(a)',end=1) line
     n = n+1
  end do
  1 rewind(unit)
  allocate(data(n,3),dimap(n))
  do i = 1, n
     read(unit,*) dimap(i), data(i,:)
  end do
  close(unit)
  ! Diodes count from 1 internally, but 0 externally
  dimap = dimap+1

  ! Normalize the data and go to mV/K instead of V/K
  data(:,1) = data(:,1)
  data(:,2:3) = data(:,2:3)!*1e3

  ! Find the best fit model
  call fit_gain_model(data, dimap, model)

  ! Possibly do a thurough analysis of the uncertainty in
  ! the absolute gain
  if(mcmc) then
     call getarg(4, wfile)
     call open_hdf_file(wfile, hfile, "r")
     call read_alloc_hdf(hfile, "time", t)
     call read_alloc_hdf(hfile, "weight", w)
     call close_hdf_file(hfile)
     call analyse_mean_gain(data, dimap, model, w, t, nsim, gain, gdev)
     write(*,'(a,f9.5,a,f9.5)') 'Season average: ', gain, ' ± ', gdev
 end if

  ! Now output the model
  open(unit, file=modfile)
  call output_model(model, unit)
  close(unit)

  ! Output residuals
  allocate(template(size(data,1)))
  call calc_shape(data(:,1), model, template)
  open(unit, file=resfile)
  do i = 1, size(data,1)
     v = model%diode_amps(dimap(i))*template(i)
     write(unit,'(i3,f15.7,4e15.7)') dimap(i)-1, data(i,1), data(i,2:3), v, (v-data(i,2))/data(i,3)
  end do
  close(unit)

contains

  subroutine fit_gain_model(data, dilist, model)
    implicit none
    real(dp)               :: data(:,:), p(4), chisq
    integer(i4b)           :: dilist(:)
    type(gain_model)       :: model
    integer(i4b)           :: ndi, err
    ndi = maxval(dilist)
    if(allocated(model%diode_amps))       deallocate(model%diode_amps)
    if(allocated(gmod_powell%diode_amps)) deallocate(gmod_powell%diode_amps)
    allocate(model%diode_amps(ndi), gmod_powell%diode_amps(ndi))
    allocate(data_powell(size(data,1),size(data,2)))
    allocate(diode_powell(size(dilist)))
    data_powell    = data
    diode_powell   = dilist
    gmod_powell%t0 = minval(data(:,1))
    p = [ -0.00124774/4.15424, 0.295267/4.15424, -0.0202857, 4.31926 ]
    call powell(p, powell_chisq, err)
    model%slope    = p(1)
    model%harm_amp = p(2)
    model%freq     = p(3)
    model%phase    = p(4)
    model%t0       = gmod_powell%t0
    call fit_linear(data, dilist, model, model%chisq)
    deallocate(data_powell, diode_powell, gmod_powell%diode_amps)
  end subroutine

  function powell_chisq(p) result(chisq)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp) :: chisq
    ! Expand to gain model format
    gmod_powell%slope    = p(1)
    gmod_powell%harm_amp = p(2)
    gmod_powell%freq     = p(3)
    gmod_powell%phase    = p(4)
    call fit_linear(data_powell, diode_powell, gmod_powell, chisq)
  end function

  subroutine fit_linear(data, dilist, model, chisq)
    implicit none
    real(dp)              :: data(:,:), chisq, t0
    real(dp), allocatable :: A(:), b(:)
    integer(i4b)          :: dilist(:), d, i, ndi
    type(gain_model)      :: model
    real(dp), allocatable :: template(:)
    ! Each diode is independent, so we can solve for one at a time
    ndi = size(model%diode_amps)
    allocate(template(size(data,1)), A(ndi), b(ndi))
    call calc_shape(data(:,1), model, template)
    ! Fit all the diode amplitudes. This includes the missing
    ! diodes, which will end up with NaN, which must be checked
    ! for later
    b = 0; A = 0
    do i = 1, size(data,1)
       d = dilist(i)
       b(d) = b(d) + template(i)*data(i,2)/data(i,3)**2
       A(d) = A(d) + template(i)**2/data(i,3)**2
    end do
    model%diode_amps = b/A
    chisq = 0
    do i = 1, size(data,1)
       chisq = chisq + (data(i,2)-model%diode_amps(dilist(i))*template(i))**2/data(i,3)**2
    end do
    chisq = chisq/size(data,1)
    deallocate(template)
  end subroutine

  subroutine calc_shape(time, model, template)
    implicit none
    real(dp)         :: time(:), template(:), t0
    type(gain_model) :: model
    t0 = model%t0
    template = 1+model%slope*(time-t0)+model%harm_amp*sin(model%freq*(time-t0)+model%phase)
  end subroutine

  subroutine output_model(model, unit)
    implicit none
    type(gain_model) :: model
    integer(i4b)     :: unit, i, ndi, d
    real(dp)         :: off, slope, harm_amp, freq, phase
    ndi = size(model%diode_amps)
    write(unit,'(a,e15.7)') "# Chisq: ", model%chisq
    do i = 1, ndi
       if(model%diode_amps(i) /= model%diode_amps(i)) cycle
       d        = i-1
       off      = 1             *model%diode_amps(i)
       slope    = model%slope   *model%diode_amps(i)
       harm_amp = model%harm_amp*model%diode_amps(i)
       freq     = model%freq
       phase    = model%phase
       write(unit,'(i4,i2,a2,f10.3,5f10.6)') d/4, modulo(d,4), "H", model%t0, off, slope, harm_amp, freq, phase
    end do
  end subroutine

  subroutine analyse_mean_gain(data, dilist, model, w, t, nsim, gain, gdev)
    implicit none
    real(dp)              :: data(:,:), w(:,:), t(:), nextra, gain, gdev
    real(dp)              :: gsum, gvar
    real(dp), allocatable :: data2(:,:)
    integer(i4b)          :: dilist(:), i, j, k, n, nsim
    type(gain_model)      :: model, model2
    type(planck_rng)      :: rng
    call rand_init(rng, 1)
    ! Based on the chisquare, find out how much larger the real
    ! scatter is than the expected scatter. This is taken to
    ! be a result of an extra unmodeled noise component, which
    ! is added as extra noise in this simulation.
    nextra = sqrt(model%chisq)
    allocate(data2(size(data,1),size(data,2)))
    data2 = data
    gsum = 0
    gvar = 0
    do i = 1, nsim
       ! Generate a new noise realization and fit a model to it
       do j = 1, size(data,1)
          data2(j,2) = data(j,2) + rand_gauss(rng)*data(j,3)*nextra
       end do
       call fit_gain_model(data2, dilist, model2)
       ! Find the season average gain for this model
       call season_gain(model2, w, t, gain)
       gsum = gsum + gain
       gvar = gvar + gain**2
       write(*,'(i5,3e15.7)') i, gain, gsum/i, sqrt(gvar/i-(gsum/i)**2)
    end do
    deallocate(data2)
    gsum = gsum/nsim
    gvar = gvar/nsim
    gvar = gvar-gsum**2
    gain = gsum
    gdev = sqrt(gvar)
  end subroutine

  ! Given weights and their time, calculate the average
  ! gain in the same units as the input gains.
  ! The gains are weighted by gain/var.
  subroutine season_gain(model, w, t, gain)
    implicit none
    type(gain_model)       :: model
    real(dp)               :: w(:,:), t(:), gain, v, denom
    real(dp),  allocatable :: template(:)
    integer(i4b)           :: d
    allocate(template(size(t)))
    call calc_shape(t, model, template)
    gain  = 0
    denom = 0
    do d = 1, size(model%diode_amps)
       if(model%diode_amps(d) /= model%diode_amps(d)) cycle
       gain  = gain  + sum((template*model%diode_amps(d))**2*w(:,d))
       denom = denom + sum(template*model%diode_amps(d)*w(:,d))
    end do
    gain = gain/denom
    deallocate(template)
  end subroutine

end program
