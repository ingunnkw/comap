program freq2corr
  use quiet_utils
  use quiet_fft_mod
  implicit none
  character(len=512)        :: cmd
  real(dp),     allocatable :: x(:), r(:), x2(:)
  complex(dpc), allocatable :: c(:)
  integer(i4b)              :: n, i, nt
  real(dp)                  :: v(2), df, dt
  call getarg(1, cmd)
  allocate(r(16),x(16))
  n = 0
  do
     read(*,*,end=2) v
     n = n+1
     if(n > size(r)) call resize(r, 2*size(r))
     if(n > size(x)) call resize(x, 2*size(x))
     x(n) = v(1)
     r(n) = v(2)
  end do
2 call resize(r, n)
  call resize(x, n)

  select case(cmd)
     case("identity")
        allocate(x2(n))
        x2 = x
     case("f2c")
        allocate(c(n),x2((n-1)*2))
        c = r
        call resize(r, (n-1)*2)
        call fft(r, c, -1)
        deallocate(c)
        r = r/size(r)**0.5
        dt=2/x(n)
        x2=(irange(size(x2))-1)*dt
     case("c2f")
        allocate(c(n/2+1),x2(n/2+1))
        r = r*size(r)**0.5
        call fft(r, c,  1)
        call resize(r, n/2+1)
        r = real(c, dp)
        deallocate(c)
        df=2/(x(2)-x(1))/(size(x2)-1)
        x2=(irange(size(x2))-1)*df
  end select

  do i = 1, size(x2)
     write(*,'(2e25.16)') x2(i), r(i)
  end do

contains

  subroutine resize(a, n)
    implicit none
    real(dp), allocatable :: a(:), b(:)
    integer(i4b)          :: n
    allocate(b(size(a)))
    b = a
    deallocate(a)
    allocate(a(n))
    a(:min(size(a),size(b))) = b(:min(size(a),size(b)))
    deallocate(b)
  end subroutine

end program
