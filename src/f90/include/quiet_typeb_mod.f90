module quiet_typeb_mod
   use quiet_utils
   use quiet_fft_mod
   use quiet_detector_mod
   !use l1_read_mod
   implicit none

   type binfo
      real(dp)                              :: preamp, sigma0
      integer(i4b)                          :: pattern
      real(dp),    dimension(:),allocatable :: x0, y0, A, cum
   end type

   type(binfo), dimension(:), allocatable          :: bs
   real(dp),                  parameter            :: bit2volt = 2d0**(-16)
   real(dp),         private, parameter            :: nsig = 5
   logical(lgt),     private                       :: active = .false.

!   private f, df, g

contains

   subroutine initialize_typeb_mod(parfile)
     implicit none
     character(len=*)   :: parfile
     character(len=512) :: fname
     logical(lgt), save :: initialized = .false.
     if(initialized) return
     call init_detector_mod(parfile)
     call get_parameter(0, parfile, 'TYPEB_CORRECT', par_lgt=active)
     if(active) then
        call get_parameter(0, parfile, 'TYPEB_FILE', par_string=fname)
        call read_binfos(fname, bs)
     end if
     initialized = .true.
   end subroutine

   ! Srate = scanfreq/samprate. Typical value: 0.001
   ! This is used to avoid being biased by strong signals when
   ! doing the internal noise measurement.
!!$   subroutine typeb_correct(data)
!!$     implicit none
!!$     type(data_struct) :: data
!!$     integer(i4b)      :: i, m, d, n, ndi
!!$     if(.not. active) return
!!$     do i = 1, size(quiet_diodes)
!!$        m = quiet_diodes(i)%horn
!!$        d = quiet_diodes(i)%sub
!!$        call typeb_correct_single(data%RQ(m+1)%demod(d,:), data%RQ(m+1)%avg(d,:), bs(i))
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine typeb_correct_single(tod, tp, b)
!!$     implicit none
!!$     real(dp)     :: tod(:), tp(:)
!!$     type(binfo)  :: b
!!$     integer(i4b) :: i, n
!!$     real(dp), dimension(:), allocatable :: wtod, wtp
!!$     if(size(b%x0) == 0) return
!!$     allocate(wtod(size(tod)), wtp(size(tp)))
!!$     wtod = tod; wtp = tp
!!$     do i = 1, 2
!!$        call measure_noise(wtod, b)
!!$        wtod = tod; wtp = tp
!!$        call apply_typeb(wtod, wtp, b)
!!$     end do
!!$     tod = wtod; tp = wtp
!!$     deallocate(wtod, wtp)
!!$   end subroutine
!!$
!!$   ! Use the upper half of the frequencies to measure
!!$   ! the noise. Extrapolate to 800 kHz. Due to intermediate
!!$   ! states, there are only 6880 samples for each 100 Hz sample.
!!$   subroutine measure_noise(tod, b)
!!$     implicit none
!!$     real(dp)                                :: tod(:)
!!$     type(binfo)                             :: b
!!$     integer(i4b)                            :: i, j, k, m, n, s, r(2), ns
!!$     real(dp),     dimension(:), allocatable :: dd
!!$     complex(dpc), dimension(:), allocatable :: ft
!!$     ! We must double demodulate to get something with the right noise
!!$     ! properties, since we would be affected by all the demod jumping otherwise.
!!$     n     = size(tod)/2
!!$     m     = n/2+1
!!$     allocate(ft(m),dd(n))
!!$     dd = (tod(2:2*n:2)-tod(1:2*n:2))/2
!!$     call fft(dd, ft, 1)
!!$     b%sigma0 = sqrt(6880*100/50*mean(real(ft(m/2:)*conjg(ft(m/2:)),dp)))*b%preamp
!!$     deallocate(ft,dd)
!!$   end subroutine
!!$
!!$   ! We know the noise, so just correct
!!$   subroutine apply_typeb(tod, tp, b)
!!$     implicit none
!!$     real(dp)    :: tod(:), tp(:), x(2), y(2)
!!$     type(binfo) :: b
!!$     integer(i4b):: i, r(2)
!!$     ! For each sample, we need y1 = tp(i)+tod(i), y2 = tp(i)-tod(i).
!!$     ! We then typeb-correct each of these, getting x1 and x2, and
!!$     ! combine these to find the true tp(i) and tod(i).
!!$     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(r,i,y,x)
!!$     r = [ 1, size(b%y0) ]
!!$     !$OMP DO
!!$     do i = 1, size(tod)
!!$        y = [ tp(i)+tod(i),       tp(i)-tod(i)       ]
!!$        x = [ deglitch(y(1),b,r), deglitch(y(2),b,r) ]
!!$        tp (i) = (x(1)+x(2))/2
!!$        tod(i) = (x(1)-x(2))/2
!!$     end do
!!$     !$OMP END DO
!!$     !$OMP END PARALLEL
!!$   end subroutine
!!$
!!$   function deglitch(y, b, r) result(x)
!!$     implicit none
!!$     real(dp)     :: y, x
!!$     type(binfo)  :: b
!!$     integer(i4b) :: r(:), n
!!$     call find_closest(y, b%y0, b%sigma0, r)
!!$     if(r(2)<r(1)) then ! empty range
!!$        x = y - b%cum(r(1))
!!$     else
!!$        x = g(y, b%sigma0, b%x0(r(1):r(2)), b%A(r(1):r(2)), b%cum(r(1)))
!!$     end if
!!$   end function
!!$
!!$   ! Adjust the range r so that it only contains y0s within nsig of y.
!!$   ! r(1) may be larger than n, and r(2) may be 0, to ensure that
!!$   ! the empty range still works at the edges.
!!$   subroutine find_closest(y, y0, sig, r)
!!$     implicit none
!!$     real(dp)     :: y, y0(:), sig
!!$     integer(i4b) :: r(:)
!!$     if(r(1) > size(y0)) r(1) = size(y0)
!!$     if(r(2) < 1)        r(2) = 1
!!$     do while(r(1) > 1)
!!$        if(y0(r(1)-1) <= y-sig*nsig) exit
!!$        r(1) = r(1)-1
!!$     end do
!!$     do while(r(1) <= size(y0))
!!$        if(y0(r(1)) >=  y-sig*nsig) exit
!!$        r(1) = r(1)+1
!!$     end do
!!$     do while(r(2) < size(y0))
!!$        if(y0(r(2)+1) >=  y+sig*nsig) exit
!!$        r(2) = r(2)+1
!!$     end do
!!$     do while(r(2) >= 1)
!!$        if(y0(r(2)) <= y+sig*nsig) exit
!!$        r(2) = r(2)-1
!!$     end do
!!$   end subroutine
!!$
!!$   function f(x, sigma, x0, A, cum) result(y)
!!$     implicit none
!!$     real(dp)     :: x, sigma, x0(:), A(:), y, cum
!!$     y = x + sum(A*(erf((x-x0)/(sqrt(2d0)*sigma))+1)/2) + cum
!!$   end function
!!$
!!$   function df(x, sigma, x0, A) result(dy)
!!$     implicit none
!!$     real(dp) :: x, x0(:), sigma, A(:), dy
!!$     dy = 1 + sum(A*exp(-0.5*((x-x0)/sigma)**2))/(sqrt(2*pi)*sigma)
!!$   end function
!!$
!!$   function g(y, sigma, x0, A, cum) result(x)
!!$     implicit none
!!$     real(dp)     :: y, x0(:), sigma, A(:), x, cum
!!$     integer(i4b) :: i
!!$     x = y
!!$     do i = 1, 4
!!$        x = x - (f(x,sigma,x0,A,cum)-y)/df(x,sigma,x0,A)
!!$     end do
!!$   end function
!!$
!!$   subroutine read_binfos(fname, bs)
!!$     implicit none
!!$     type(binfo),       dimension(:), allocatable :: bs
!!$     character(len=*)   :: fname
!!$     character(len=256) :: line, word
!!$     character(len=64)  :: words(16)
!!$     integer(i4b)       :: nmod, ndi, d, i, j, k, m, n, unit, ntok
!!$     real(dp)           :: gtmp(1024,4), tmp(1024), fact
!!$     call free_binfos(bs)
!!$     nmod = size(quiet_horns)
!!$     ndi  = quiet_horns(1)%ndi
!!$     n    = nmod*ndi
!!$     fact = 2
!!$     allocate(bs(n))
!!$     unit = getlun()
!!$     open(unit,file=fname,status="old",action="read")
!!$     do
!!$        read(unit,'(a)',end=1) line
!!$        m = len_trim(line)
!!$        if(m == 0)           cycle
!!$        if(line(1:1) == "#") cycle
!!$        call get_tokens(line(1:m), " 	", words, ntok, maxnum=16)
!!$        select case(words(1))
!!$           case("__HEIGHT_FACTOR__")
!!$              read(words(2),*) fact
!!$           case("__PREAMP_RMS_FACTOR__")
!!$              do
!!$                 read(unit,'(a)') line
!!$                 if(line == "") exit
!!$                 if(line(1:1) == '#') cycle
!!$                 read(line,*) m, bs(m*ndi+1:(m+1)*ndi)%preamp
!!$              end do
!!$           case("__CENTER_HEIGHT_PATTERN__")
!!$              do
!!$                 read(unit,'(a)') line
!!$                 if(line == "") exit
!!$                 if(line(1:1) == '#') cycle
!!$                 read(line,*) m, bs(m*ndi+1:(m+1)*ndi)%pattern
!!$              end do
!!$           case("__CENTERS_HEIGHTS__")
!!$              read(words(2),*) i
!!$              j = 0
!!$              do
!!$                 read(unit,'(a)') line
!!$                 if(line == "") exit
!!$                 if(line(1:1) == '#') cycle
!!$                 j = j+1
!!$                 read(line,*) gtmp(j,:)
!!$              end do
!!$              gtmp(1:j,1) = gtmp(1:j,1)*bit2volt
!!$              gtmp(1:j,3) = gtmp(1:j,3)*bit2volt*fact
!!$              tmp(1) = 0
!!$              do k = 1, j
!!$                 tmp(k+1) = tmp(k) + gtmp(k,3)
!!$              end do
!!$              do k = 1, n
!!$                 if(bs(k)%pattern /= i) cycle
!!$                 allocate(bs(k)%x0(j), bs(k)%y0(j), bs(k)%A(j), bs(k)%cum(j+1))
!!$                 bs(k)%y0  = gtmp(1:j,1)
!!$                 bs(k)%x0  = gtmp(1:j,1) - tmp(1:j)
!!$                 bs(k)%A   = gtmp(1:j,3)
!!$                 bs(k)%cum = tmp
!!$              end do
!!$        end select
!!$     end do
!!$1    close(unit)
!!$     do k = 1, n
!!$        if(.not. allocated(bs(k)%x0))  allocate(bs(k)%x0(0))
!!$        if(.not. allocated(bs(k)%y0))  allocate(bs(k)%y0(0))
!!$        if(.not. allocated(bs(k)%A))   allocate(bs(k)%A(0))
!!$        if(.not. allocated(bs(k)%cum)) allocate(bs(k)%cum(1))
!!$     end do
!!$
!!$   end subroutine
!!$
!!$   subroutine free_binfos(bs)
!!$     implicit none
!!$     type(binfo), dimension(:), allocatable :: bs
!!$     integer(i4b)                           :: i
!!$     if(.not. allocated(bs)) return
!!$     do i = 1, size(bs)
!!$        if(allocated(bs(i)%x0))  deallocate(bs(i)%x0)
!!$        if(allocated(bs(i)%y0))  deallocate(bs(i)%y0)
!!$        if(allocated(bs(i)%A))   deallocate(bs(i)%A)
!!$        if(allocated(bs(i)%cum)) deallocate(bs(i)%cum)
!!$     end do
!!$     deallocate(bs)
!!$   end subroutine

end module
