module quiet_constrained_mod
   use quiet_utils
   use quiet_fft_mod
   use rngmod
   use math_tools
   implicit none

contains

   subroutine constrained_realization(ps, data, mask, rng, out, mode)
     implicit none
     real(dp),    intent(in)         :: ps(:), data(:)
     real(dp),    intent(out)        :: out(:)
     logical(lgt),intent(in)         :: mask(:)
     type(planck_rng) :: rng
     integer(i4b)     :: status
     logical(lgt), optional :: mode
     call constrained_realization_bcg(ps, data, mask, rng, out, status,mode=mode)
     if(status /= 0) then
        write(*,*)'OBS constrained_realization_bcg status =', status
        call constrained_realization_mat(ps, data, mask, rng, out, mode=mode)
     end if
   end subroutine constrained_realization

   subroutine constrained_realization_mat(ps, data, mask, rng, out, mode)
     implicit none
     real(dp),    intent(in)         :: ps(:), data(:)
     real(dp),    intent(out)        :: out(:)
     logical(lgt),intent(in)         :: mask(:)
     logical(lgt)     :: find_mode
     type(planck_rng) :: rng
     logical(lgt),                    optional :: mode
     real(dp),     dimension(:),   allocatable :: afull, cafull, icorr, ca, y, yfull
     real(dp),     dimension(:),   allocatable :: x, xfull
     complex(dpc), dimension(:),   allocatable :: ft,  fafull, fcafull
     integer(i4b), dimension(:),   allocatable :: inds
     real(dp),     dimension(:,:), allocatable :: bmat, bchol
     integer(i4b) :: i, j, k, m, n, p

!write(*,*) "using slow method"
     find_mode = .false.; if(present(mode)) find_mode = mode
     if(all(mask)) then; out = data; return; end if
     call iwhere(.not. mask, inds)
     n = size(data); m = size(ps); p = size(inds)

     ! First find the expectation value in the gap.
     ! First calculate Ca by matmul(N^-1, mask*data)
     allocate(afull(n),cafull(n),fafull(m),fcafull(m),ca(p))
     afull       = data
     afull(inds) = 0
     call fft(afull, fafull,   1)
     fcafull = 0
     where(ps/=0) fcafull = fafull/ps
     call fft(cafull,  fcafull, -1)
     ca = cafull(inds)
     deallocate(afull, cafull, fafull, fcafull)

     ! We must now solve By = Ca. For this, we need the
     ! pizel-space inverse correlation.
     allocate(ft(m), icorr(0:n-1), bmat(p,p))
     ft = 0
     where(ps/=0) ft = 1/ps
     call fft(icorr, ft, -1)
     icorr = icorr / sqrt(real(n))
     do i = 1, p
        do j = 1, p
           k = inds(i)-inds(j)
           bmat(j,i) = icorr(abs(k))
        end do
     end do
     deallocate(ft, icorr)
     call invert_matrix(bmat)

     ! Finally, solve the equation
     allocate(y(p))
     y = -matmul(bmat, ca)

     ! And generate a new realization
     allocate(bchol(p,p),x(p))
     if(find_mode) then
        x = y
     else
        call cholesky_decompose(bmat, bchol)
        do i = 1, p
           x(i) = rand_gauss(rng)
        end do
        x = matmul(bchol, x)+y
     end if

     out = data
     out(inds) = x
     deallocate(x, y, ca, bmat, bchol)
   end subroutine

   subroutine constrained_realization_bcg(ps, data, mask, rng, out, status, mode)
     implicit none
     real(dp),    intent(in)         :: ps(:), data(:)
     real(dp),    intent(out)        :: out(:)
     logical(lgt),intent(in)         :: mask(:)
     logical(lgt)     :: find_mode
     real(dp)         :: sc
     type(planck_rng) :: rng
     type(bcg_search) :: bcg
     integer(i4b)     :: i, j, k, m, n, maxit
     integer(i4b), optional                    :: status
     logical(lgt), optional                    :: mode
     integer(i4b), dimension(:),   allocatable :: inds, rest
     complex(dpc), dimension(:),   allocatable :: sf
     real(dp),     dimension(:),   allocatable :: rhs, tmp

     find_mode = .false.; if(present(mode)) find_mode = mode
     if(all(mask)) then; out = data; return; end if
     if(present(status)) status = 0
     call iwhere(mask, inds)
     call iwhere(.not. mask, rest)
     n = size(data); m = size(ps)
     sc = mean(ps(m/2:))
     maxit = 1000

     ! Will solve [ s_1, (iSs)_2 ] = [ d_1, (icSr)_2 ] by
     ! conjgate gradients for s. r is a random array,
     ! while iS is 1/ps. If iS is far from 1, we need to
     ! precondition by multiplying with [1, mean(S)]

     allocate(sf(m),rhs(n), tmp(n))
     do i = 1, m
        sf(i) = cmplx(rand_gauss(rng), rand_gauss(rng))/sqrt(2d0)
     end do
     where(ps/=0) sf = sf/ps**0.5
     call fft(rhs,sf,-1)
     if(find_mode) rhs = rhs*1e-20
     rhs(inds) = data(inds)

     ! Now start the CG search
     call bcg_init(bcg, rhs)
     do i = 1, maxit
        tmp = bcg%p; tmp(rest) = tmp(rest)*sc
        call bcg1(bcg, tmp)

        call fft(bcg%y, sf, 1)
        where(ps/=0) sf = sf / ps
        call fft(tmp,   sf,-1)
        tmp(inds) = bcg%y(inds)
        call bcg2(bcg, tmp)

        tmp = bcg%s; tmp(rest) = tmp(rest)*sc
        call bcg3(bcg, tmp)

        call fft(bcg%z, sf, 1)
        where(ps/=0) sf = sf / ps
        call fft(tmp,   sf,-1)
        tmp(inds) = bcg%z(inds)
        call bcg4(bcg, tmp)

        tmp = bcg%t; tmp(rest) = tmp(rest)*sc
        call bcg5(bcg, tmp)

        ! Test
        call fft(bcg%x, sf, 1)
        where(ps/=0) sf = sf / ps
        call fft(tmp,   sf,-1)
        tmp(inds) = bcg%z(inds)

        if(bcg%err < 1e-6) exit
     end do

     out = bcg%x
     call bcg_free(bcg)
     if(present(status) .and. i > maxit) status = 1
     deallocate(sf, rhs, tmp, inds, rest)
   end subroutine

   subroutine iwhere(mask, w)
     implicit none
     logical(lgt) :: mask(:)
     integer(i4b) :: i, n
     integer(i4b), dimension(:), allocatable :: w
     n = count(mask)
     allocate(w(n))
     n = 0
     do i = 1, size(mask)
        if(.not. mask(i)) cycle
        n = n+1
        w(n) = i
     end do
   end subroutine

end module
