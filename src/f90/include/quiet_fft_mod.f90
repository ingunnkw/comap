! ***********************************************************************
!
!      Module for performing fourier transformations and bandpass
!      filtering.
!
!
!
! ***********************************************************************

module quiet_fft_mod
  use healpix_types
  implicit none
  include 'fftw3.f'

  interface fft
     module procedure fft_outplace_cplx_sp, fft_outplace_cplx_dp
  end interface

  interface fft_multi
     module procedure fft_multi_cplx_sp, fft_multi_cplx_dp, fft_multi2_cplx_dp, fft_multi2_cplx_sp
  end interface

  interface extract_powspec
     module procedure extract_powspec_tod_sp, extract_powspec_tod_dp, extract_powspec_cplx_sp, extract_powspec_cplx_dp
  end interface

contains

  subroutine fft_thread_init
    implicit none
    logical(lgt), save :: initialized = .false.
    integer(i4b)       :: ret, omp_get_max_threads, n
    if(initialized) return
    !$OMP PARALLEL NUM_THREADS(1)
    n = omp_get_max_threads()
    !$OMP END PARALLEL
    call dfftw_init_threads(ret)
    call sfftw_init_threads(ret)
    call dfftw_plan_with_nthreads(n)
    call sfftw_plan_with_nthreads(n)
    initialized = .true.
  end subroutine

  subroutine fft_outplace_cplx_dp(tod_real, tod_fft, direction)
    implicit none
    real(dp),    dimension(0:), intent(inout) :: tod_real
    complex(dpc), dimension(0:), intent(inout) :: tod_fft
    integer :: direction
    integer*8 :: plan

    call fft_thread_init
    if (size(tod_real)/2+1 /= size(tod_fft)) then
       write(*,*) 'quiet_fft_mod -- Incompatible array lengths in fft_outplace'
       stop
    end if

    if(direction > 0) then
       call dfftw_plan_dft_r2c_1d(plan, size(tod_real), tod_real, tod_fft, &
            & fftw_estimate + fftw_unaligned)
       call dfftw_execute_dft_r2c(plan, tod_real, tod_fft)
       call dfftw_destroy_plan(plan)
       tod_fft = tod_fft/sqrt(real(size(tod_real),dp))
    else
       call dfftw_plan_dft_c2r_1d(plan, size(tod_real), tod_fft, tod_real, &
            & fftw_estimate + fftw_unaligned)
       call dfftw_execute_dft_c2r(plan, tod_fft, tod_real)
       call dfftw_destroy_plan(plan)
       tod_real = tod_real/sqrt(real(size(tod_real),dp))
    end if
  end subroutine fft_outplace_cplx_dp

  subroutine fft_outplace_cplx_sp(tod_real, tod_fft, direction)
    implicit none
    real(sp),     dimension(0:), intent(inout) :: tod_real
    complex(spc), dimension(0:), intent(inout) :: tod_fft
    integer :: direction
    integer*8 :: plan

    call fft_thread_init
    if (size(tod_real)/2+1 /= size(tod_fft)) then
       write(*,*) 'quiet_fft_mod -- Incompatible array lengths in fft_outplace'
       stop
    end if

    if(direction > 0) then
       call sfftw_plan_dft_r2c_1d(plan, size(tod_real), tod_real, tod_fft, &
            & fftw_estimate + fftw_unaligned)
       call sfftw_execute_dft_r2c(plan, tod_real, tod_fft)
       call sfftw_destroy_plan(plan)
       tod_fft = tod_fft/sqrt(real(size(tod_real),sp))
    else
       call sfftw_plan_dft_c2r_1d(plan, size(tod_real), tod_fft, tod_real, &
            & fftw_estimate + fftw_unaligned)
       call sfftw_execute_dft_c2r(plan, tod_fft, tod_real)
       call sfftw_destroy_plan(plan)
       tod_real = tod_real/sqrt(real(size(tod_real),sp))
    end if
  end subroutine fft_outplace_cplx_sp

  subroutine fft_multi_cplx_dp(tod_real, tod_fft, direction)
    implicit none
    real(dp), dimension(:,:), intent(inout) :: tod_real
    complex(dpc), dimension(:,:), intent(inout) :: tod_fft
    real(dp)                                :: norm
    integer(i8b)                            :: plan
    integer(i4b)                            :: n, nf, ndi, direction, i
    n   = size(tod_real, 1)
    nf  = size(tod_fft,  1)
    ndi = size(tod_real, 2)
    norm = 1d0/sqrt(real(n,dp))
    call fft_thread_init

    if(direction > 0) then
       call dfftw_plan_many_dft_r2c(plan, 1, n, ndi, tod_real, &
        & n, 1, n, tod_fft, nf, 1, nf, fftw_unaligned + fftw_estimate)
       call dfftw_execute_dft_r2c(plan, tod_real, tod_fft)
       call dfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_fft = tod_fft * norm
       !$OMP END PARALLEL WORKSHARE
    else
       call dfftw_plan_many_dft_c2r(plan, 1, n, ndi, tod_fft, &
        & nf, 1, nf, tod_real, n, 1, n, fftw_unaligned + fftw_estimate)
       call dfftw_execute_dft_c2r(plan, tod_fft, tod_real)
       call dfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_real = tod_real * norm
       !$OMP END PARALLEL WORKSHARE
    end if
  end subroutine

  subroutine fft_multi_cplx_sp(tod_real, tod_fft, direction)
    implicit none
    real(sp),    dimension(:,:), intent(inout) :: tod_real
    complex(spc), dimension(:,:), intent(inout) :: tod_fft
    real(sp)                                :: norm
    integer(i8b)                            :: plan
    integer(i4b)                            :: n, nf, ndi, direction, i
    n   = size(tod_real, 1)
    nf  = size(tod_fft,  1)
    ndi = size(tod_real, 2)
    norm = 1.0/sqrt(real(n,sp))
    call fft_thread_init

    if(direction > 0) then
       call sfftw_plan_many_dft_r2c(plan, 1, n, ndi, tod_real, &
        & n, 1, n, tod_fft, nf, 1, nf, fftw_unaligned + fftw_estimate)
       call sfftw_execute_dft_r2c(plan, tod_real, tod_fft)
       call sfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_fft = tod_fft * norm
       !$OMP END PARALLEL WORKSHARE
    else
       call sfftw_plan_many_dft_c2r(plan, 1, n, ndi, tod_fft, &
        & nf, 1, nf, tod_real, n, 1, n, fftw_unaligned + fftw_estimate)
       call sfftw_execute_dft_c2r(plan, tod_fft, tod_real)
       call sfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_real = tod_real * norm
       !$OMP END PARALLEL WORKSHARE
    end if
  end subroutine


  ! O, woe is me. Why do I need to write the same function
  ! so many times?
  subroutine fft_multi2_cplx_dp(tod_real, tod_fft, direction)
    implicit none
    real(dp),     dimension(:,:,:), intent(inout) :: tod_real
    complex(dpc), dimension(:,:,:), intent(inout) :: tod_fft
    real(dp)                                :: norm
    integer(i8b)                            :: plan
    integer(i4b)                            :: n, nf, ndi, ns, direction, i, j
    nf  = size(tod_fft,  1)
    n   = size(tod_real, 1)
    ndi = size(tod_real, 2)
    ns  = size(tod_real, 3)
    norm = 1d0/sqrt(real(n,dp))
    call fft_thread_init

    if(direction > 0) then
       call dfftw_plan_many_dft_r2c(plan, 1, n, ndi*ns, tod_real, &
        & n, 1, n, tod_fft, nf, 1, nf, fftw_unaligned + fftw_estimate)
       call dfftw_execute_dft_r2c(plan, tod_real, tod_fft)
       call dfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_fft = tod_fft * norm
       !$OMP END PARALLEL WORKSHARE
    else
       call dfftw_plan_many_dft_c2r(plan, 1, n, ndi*ns, tod_fft, &
        & nf, 1, nf, tod_real, n, 1, n, fftw_unaligned + fftw_estimate)
       call dfftw_execute_dft_c2r(plan, tod_fft, tod_real)
       call dfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_real = tod_real * norm
       !$OMP END PARALLEL WORKSHARE
    end if
  end subroutine

  subroutine fft_multi2_cplx_sp(tod_real, tod_fft, direction)
    implicit none
    real(sp),     dimension(:,:,:), intent(inout) :: tod_real
    complex(spc), dimension(:,:,:), intent(inout) :: tod_fft
    real(sp)                                :: norm
    integer(i8b)                            :: plan
    integer(i4b)                            :: n, nf, ndi, ns, direction, i, j
    nf  = size(tod_fft,  1)
    n   = size(tod_real, 1)
    ndi = size(tod_real, 2)
    ns  = size(tod_real, 3)
    norm = 1d0/sqrt(real(n,sp))
    call fft_thread_init

    if(direction > 0) then
       call sfftw_plan_many_dft_r2c(plan, 1, n, ndi*ns, tod_real, &
        & n, 1, n, tod_fft, nf, 1, nf, fftw_unaligned + fftw_estimate)
       call sfftw_execute_dft_r2c(plan, tod_real, tod_fft)
       call sfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_fft = tod_fft * norm
       !$OMP END PARALLEL WORKSHARE
    else
       call sfftw_plan_many_dft_c2r(plan, 1, n, ndi*ns, tod_fft, &
        & nf, 1, nf, tod_real, n, 1, n, fftw_unaligned + fftw_estimate)
       call sfftw_execute_dft_c2r(plan, tod_fft, tod_real)
       call sfftw_destroy_plan(plan)
       !$OMP PARALLEL WORKSHARE
       tod_real = tod_real * norm
       !$OMP END PARALLEL WORKSHARE
    end if
  end subroutine

  function freq2ind(freq, samp_rate, n) result(ind)
    implicit none
    real(dp) :: freq, samp_rate
    integer(i4b) :: n, ind
    ind = int(freq/(samp_rate/2)*(n-1))+1
  end function freq2ind

  ! Rescale from an integer index into a power spectrum
  ! or similar frequency representation to the actual
  ! frequency. The conversion is rather trivial, but
  ! having it explicitly in a function makes what's happening
  ! more explicit. Rescales the interval 1:n -> 0:srate/2.
  function ind2freq(ind, samp_rate, n) result(freq)
    implicit none
    real(dp) :: freq, samp_rate
    integer(i4b) :: n, ind
    freq = (ind-1)*(samp_rate/2)/(n-1)
  end function ind2freq

  subroutine extract_powspec_tod_dp(tod,powspec)
    implicit none
    real(dp)                  :: tod(:), powspec(:)
    complex(dpc), allocatable :: ft(:)
    allocate(ft(size(powspec)/2))
    call fft(tod, ft, 1)
    call extract_powspec(ft, powspec)
    deallocate(ft)
  end subroutine

  subroutine extract_powspec_tod_sp(tod,powspec)
    implicit none
    real(sp)                  :: tod(:), powspec(:)
    complex(spc), allocatable :: ft(:)
    allocate(ft(size(powspec)/2))
    call fft(tod, ft, 1)
    call extract_powspec(ft, powspec)
    deallocate(ft)
  end subroutine

  subroutine extract_powspec_cplx_dp(ft,powspec)
    implicit none
    complex(dpc), dimension(:) :: ft
    real(dp),     dimension(:) :: powspec

    integer(i4b) :: i,j,k,n,m
    if(size(powspec) < size(ft)) then
       write(*,*) "Second argument is too small to hold power spectrum!"
       stop
    end if
    powspec = real(ft * conjg(ft),dp)
  end subroutine

  subroutine extract_powspec_cplx_sp(ft,powspec)
    implicit none
    complex(spc), dimension(:) :: ft
    real(sp),     dimension(:) :: powspec

    integer(i4b) :: i,j,k,n,m
    if(size(powspec) < size(ft)) then
       write(*,*) "Second argument is too small to hold power spectrum!"
       stop
    end if
    powspec = real(ft * conjg(ft),sp)
  end subroutine

end module quiet_fft_mod
