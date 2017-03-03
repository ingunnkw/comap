! Fits one amplitude per diode. Inputs are gain measurements:
!   mjd di gain dgain
! And template points
!   mjd di val
! The template points are binned in mjd, and a nonlinear search
! is performed per diode to find the best amplitude.
program gainfit2
  use quiet_utils
  use powell_mod
  use quiet_fileutils
  implicit none
  character(len=512)         :: arg, gfile, tfile, otfile, sfile, mfile, resfile
  integer(i4b)               :: unit, di, n, m, ndi, nbin, i,j,k, b, ni,nj,ndof
  real(dp)                   :: mjd, val, dval, mjd_range(2), bsize, amp, llik, chisq, chisq2, tmp(3)
  real(dp),    allocatable   :: gains(:,:,:), template(:,:,:), amps(:)
  real(dp),    allocatable   :: bcorr(:,:,:)
  integer(i4b),allocatable   :: templens(:,:), gainlens(:)

  ! For powell
  real(dp),     pointer :: pgain(:,:), ptemp(:,:), pmjd(:)
  integer(i4b), pointer :: plen(:)

  bsize = 5
  call getarg(1, gfile)
  call getarg(2, tfile)
  call getarg(3, otfile)
  call getarg(4, sfile)
  call getarg(5, mfile)
  call getarg(6, resfile)
  ! Read the templates. First scan the file to find the data properties...
  unit = getlun()
  n     = 0
  ndi   = 0
  open(unit,file=tfile,action="read",status="old")
  do
     read(unit,*,end=1) mjd, j, di, val
     n = n+1
     if(n==1) mjd_range = [mjd,mjd]
     mjd_range(1) = min(mjd_range(1),mjd)
     mjd_range(2) = max(mjd_range(2),mjd)
     di = di+1
     ndi = max(ndi,di)
  end do
  1 rewind(unit)
  nbin = ceiling((mjd_range(2)-mjd_range(1))/bsize)+1
  allocate(template(n,nbin,ndi),templens(nbin,ndi))
  ! Then read it in properly
  templens = 0
  do i = 1, n
     read(unit,*) mjd, j, di, val
     di = di+1
     b = min(nbin,floor((mjd-mjd_range(1))/(mjd_range(2)-mjd_range(1))*nbin)+1)
     templens(b,di) = templens(b,di) + 1
     template(templens(b,di),b,di) = val
  end do
  close(unit)

  ! Then get the actual measurements
  allocate(gainlens(ndi))
  gainlens = 0
  m        = 0
  open(unit,file=gfile,action="read",status="old")
  do
     read(unit,*,end=2) di, mjd, val, dval
     di = di+1
     gainlens(di) = gainlens(di)+1
     m  = m+1
  end do
  2 rewind(unit)
  allocate(gains(3,maxval(gainlens),ndi))
  gainlens = 0
  do i = 1, m
     read(unit,*) di, mjd, val, dval
     di = di+1
     gainlens(di) = gainlens(di)+1
     gains(:,gainlens(di),di) = [mjd,val,dval]
  end do
  close(unit)

  ! Phew! Reading in files is so cumbersome in Fortran.
  allocate(amps(ndi))
  chisq = 0; ndof = 0; chisq2 = 0
  do di = 1, ndi
     call fit_template(gains(:,:gainlens(di),di), template(:,:,di), templens(:,di), mjd_range, amps(di), llik)
     ndof  = ndof  + gainlens(di)
     chisq = chisq + 2*llik
     chisq2= chisq2+ 2*loglik_mean(gains(:,:gainlens(di),di), template(:,:,di), templens(:,di), mjd_range, amps(di))
  end do
  write(*,*) chisq, chisq2, ndof

  ! Output the model
  open(unit,file=mfile)
  do di = 1, ndi
     write(unit,'(i3,i2,a2,e15.7)') (di-1)/4, modulo(di-1,4), "N", amps(di)
  end do
  close(unit)

  ! Output scaled template
  open(unit,file=otfile)
  do di = 1, ndi
     do i = 1, nbin
        do j = 1, templens(i,di)
           write(unit,'(i3,3e15.7)') di-1, (i-1)*(mjd_range(2)-mjd_range(1))/nbin+mjd_range(1), template(j,i,di)*amps(di), amps(di)
        end do
     end do
  end do
  close(unit)

  ! Output the residuals
  open(unit,file=resfile)
  do di = 1, ndi
     do i = 1, gainlens(di)
        mjd = gains(1,i,di)
        val = gains(2,i,di)
        dval= gains(3,i,di)
        b   = floor((mjd-mjd_range(1))/(mjd_range(2)-mjd_range(1))*nbin)+1
        if(b <= 0) b = 1
        if(b > nbin) b = nbin
        tmp(1) = sum(template(:templens(b,di),b,di))/templens(b,di)*amps(di)
        tmp(2)  = sqrt(dval**2+amps(di)**2*variance(template(:templens(b,di),b,di)/nbin))
        write(unit,'(i4,f15.7,4e15.7)') di-1, mjd, val, tmp(2), &
         & tmp(1), (val-tmp(1))/tmp(2)
     end do
  end do
  close(unit)

  ! Slices
  open(unit,file=sfile)
  do i = 1, 1000
     amp = 1e6*(i-1)/1000
     write(unit,advance="no",fmt="(e15.7)") amp
     do di = 1, ndi
        write(unit,advance="no",fmt="(e15.7)") loglik(gains(:,:gainlens(di),di), template(:,:,di), templens(:,di), mjd_range, amp)
     end do
     write(unit,*)
  end do
  close(unit)

contains

  subroutine fit_template(gains, template, lens, mjd_range, amp, lik)
    implicit none
    real(dp),     target :: gains(:,:), template(:,:), mjd_range(:)
    real(dp)             :: amp, meantemp(2), p(1), lik
    integer(i4b), target :: lens(:)
    integer(i4b)         :: err
    pgain => gains
    ptemp => template
    plen  => lens
    pmjd  => mjd_range
    ! Find starting guess
    meantemp = 0
    do i = 1, size(template,2)
       do j = 1, lens(i)
          meantemp(1) = meantemp(1) + template(j,i)
          meantemp(2) = meantemp(2) + 1
       end do
    end do
    p(1) = mean(gains(2,:))/(meantemp(1)/meantemp(2))
    call powell(p, powell_lik, err)
    amp = p(1)
    lik = powell_lik(p)
  end subroutine

  function powell_lik(p) result(lik)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp) :: lik
    lik = loglik_mean(pgain, ptemp, plen, pmjd, p(1))
  end function

  function loglik(gains, template, lens, mjd_range, amp)
    implicit none
    real(dp)        :: gains(:,:), template(:,:), mjd_range(:), amp
    real(dp)        :: mjd, val, dval, loglik, lik
    integer(i4b)    :: lens(:), i, b
    loglik = 0
    do i = 1, size(gains,2)
       mjd = gains(1,i)
       val = gains(2,i)
       dval= gains(3,i)
       b   = floor((mjd-mjd_range(1))/(mjd_range(2)-mjd_range(1))*size(lens))+1
       if(b <= 0) b = 1
       if(b > size(lens)) b = size(lens)
       if(lens(b)==0) cycle
       lik = 0
       do j = 1, lens(b)
          lik = lik + exp(-0.5*((val-template(j,b)*amp)/dval)**2)
       end do
       loglik = loglik - log(lik/sqrt(2*pi*dval**2)/lens(b))
    end do
  end function

  function loglik_mean(gains, template, lens, mjd_range, amp) result(loglik)
    implicit none
    real(dp)        :: gains(:,:), template(:,:), mjd_range(:), amp
    real(dp)        :: mjd, val, dval, loglik, lik, mean, dev
    integer(i4b)    :: lens(:), i, b, n
    loglik = 0
    do i = 1, size(gains,2)
       mjd = gains(1,i)
       val = gains(2,i)
       dval= gains(3,i)
       b   = floor((mjd-mjd_range(1))/(mjd_range(2)-mjd_range(1))*size(lens))+1
       if(b <= 0) b = 1
       if(b > size(lens)) b = size(lens)
       if(lens(b) == 0) cycle
       mean = sum(template(:lens(b),b))/lens(b)
       dev  = sqrt(variance(template(:lens(b),b)))/lens(b)**0.5
       if(dev/=dev)dev = 0
!write(*,'(5e15.7)') val, mean*amp, dval, dev*amp, 0.5*(val-mean*amp)**2/(dval**2+dev**2*amp**2)
       loglik = loglik + 0.5*(val-mean*amp)**2/(dval**2+dev**2*amp**2*0)
    end do
  end function

end program
