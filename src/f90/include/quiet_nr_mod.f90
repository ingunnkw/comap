!****************************************************************************
!*                                                                          *
!*                Collection of Numerical Recipes routines                  *
!*                                                                          *
!*           Written by Sigurd K. Næss and Hans Kristian Eriksen, Oslo      *
!*                                                                          *
!****************************************************************************
module quiet_nr_mod
  use healpix_types
  use math_tools
  implicit none

  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

contains

  subroutine rtnewt_erfc(F_in, x, info)
    implicit none

    real(dp),     intent(in)  :: F_in
    real(dp),     intent(out) :: x
    integer(i4b), intent(out) :: info

    integer(i4b) :: j, JMAX = 10000
    real(dp) :: xacc, df, dx, f, one_over_sqrt2, sqrt_two_over_pi

    info             = 0
    xacc             = 1.e-6
    one_over_sqrt2   = 1.0/sqrt(2.0)
    sqrt_two_over_pi = sqrt(2.0/pi)

    x = 0.0
    do j = 1, JMAX
       f  = calerf(-x*one_over_sqrt2, 1) - 2. * F_in
       df = sqrt_two_over_pi * exp(-0.5*x**2)
       dx = f / df
       x  = x - dx
       if (abs(dx) < xacc) then
          return
       end if
    end do

    info = 1
!    write(*,*) 'Too many iterations in rtnewt_erfc'
!    stop

  end subroutine rtnewt_erfc



  subroutine MOV3(a1, b1, c1, a2, b2, c2) 
    implicit none

    real(dp), intent(inout) :: a1, b1, c1, a2, b2, c2

    real(dp) :: temp

    temp = a2
    a2   = a1
    a1   = temp

    temp = b2
    b2   = b1
    b1   = temp

    temp = c2
    c2   = c1
    c1   = temp

  end subroutine MOV3



  ! Computes the error function
  real(dp) function erf(x)
    implicit none

    real(dp), intent(in) :: x

    real(dp) :: ans, z, t, a

    erf = calerf(x, 0)
    return

!    if (abs(x) < 7.d0) then
       z = abs(x)
       t=1.d0/(1.d0+0.5d0*z)
       ans = t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(0.37409196d0+t*(0.09678418d0+&
            & t*(-0.18628806d0+t*(0.27886807d0+t*(-1.13520398d0+t*(1.48851587d0+&
            & t*(-0.82215223d0+t*0.17087277d0)))))))))
       if (x >= 0.d0) then
          erf = 1.d0-ans
       else
          erf = ans - 1.d0
       end if
 !   else
!       a = -8.d0 * (pi-3.d0) / (3.d0*pi *(pi-4.d0))
!       erf = sqrt(1.d0 - exp(-x**2 * (4.d0/pi + a*x**2) / (1.d0 + a*x**2)))
!       if (x < 0.d0) erf = -erf

!       z = abs(z)
!       t = 1

!    end if

  end function erf



  function calerf(arg, jint )

!*********************************************************************
!
!! CALERF !omputes various forms of the error fun!tion.
!
!  Dis!ussion:
!
!    This routine evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
!    for a real argument x.  
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    William Cody
!
!  Referen!e:
!
!    William Cody,
!    Rational Chebyshev Approximations for the Error Fun!tion,
!    Mathemati!s of !omputation,
!    Volume 23, Number 107, July 1969, pages 631-638.
!
!  Parameters:
!
!    Input, double precision ARG, the argument.  If JINT is 1, the
!    argument must be less than XBIG.  If JINT is 2, the argument
!    must lie between XNEG and XMAX.
!
!    Output, double precision RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = erf(x);
!    1, RESULT = erfc(x) = 1 - erf(x);
!    2, RESULT = exp(x*x)*erfc(x) = exp(x*x) - erf(x*x)*erf(x).
!
!    Input, integer JINT, chooses the function to be !omputed.
!    0, erf(x);
!    1, erfc(x);
!    2, exp(x*x)*erfc(x).

    implicit none

    real(dp),     intent(in) :: arg
    integer(i4b), intent(in) :: jint
    real(dp)                 :: calerf

    real(dp)     :: del
    real(dp)     :: x, xden
    real(dp)     :: xnum, small, y, ysq
    integer(i4b) :: i

!
!  Mathematical constants
!
    real(dp), parameter :: four   = 4.0d0
    real(dp), parameter :: one    = 1.0d0
    real(dp), parameter :: half   = 0.5d0
    real(dp), parameter :: two    = 2.0d0
    real(dp), parameter :: zero   = 0.0d0
    real(dp), parameter :: sqrpi  = 5.6418958354775628695d-1
    real(dp), parameter :: thresh = 0.46875d0
    real(dp), parameter :: sixten = 16.0d0
!
!  Machine-dependent constants
!
    real(dp), parameter :: xinf   = 1.79d308
    real(dp), parameter :: xneg   = -26.628d0
    real(dp), parameter :: xsmall = 1.11d-16
    real(dp), parameter :: xbig   = 26.543d0
    real(dp), parameter :: xhuge  = 6.71d7
    real(dp), parameter :: xmax = 2.53d307
!
!  Coefficients for approximation to  erf  in first interval
!
    real(dp), dimension(5), parameter ::  a = (/3.16112374387056560d00,1.13864154151050156d02, &
         &       3.77485237685302021d02,3.20937758913846947d03, &
         &       1.85777706184603153d-1/)
    real(dp), dimension(4), parameter ::  b = (/2.36012909523441209d01,2.44024637934444173d02, &
         &       1.28261652607737228d03,2.84423683343917062d03/)
!
!  Coefficients for approximation to  erfc  in second interval
!
    real(dp), dimension(9), parameter :: c = (/5.64188496988670089d-1,8.88314979438837594d0, &
         &       6.61191906371416295d01,2.98635138197400131d02, &
         &       8.81952221241769090d02,1.71204761263407058d03, &
         &       2.05107837782607147d03,1.23033935479799725d03, &
         &       2.15311535474403846d-8/)
    real(dp), dimension(8), parameter :: d = (/1.57449261107098347d01,1.17693950891312499d02, &
         &       5.37181101862009858d02,1.62138957456669019d03, &
         &       3.29079923573345963d03,4.36261909014324716d03, &
         &       3.43936767414372164d03,1.23033935480374942d03/)
!
!  !oefficients for approximation to  erfc  in third interval
!
    real(dp), dimension(6), parameter :: p = (/3.05326634961232344d-1,3.60344899949804439d-1, &
         &       1.25781726111229246d-1,1.60837851487422766d-2, &
         &       6.58749161529837803d-4,1.63153871373020978d-2/)
    real(dp), dimension(5), parameter :: q = (/2.56852019228982242d00,1.87295284992346047d00, &
         &       5.27905102951428412d-1,6.05183413124413191d-2, &
         &       2.33520497626869185d-3/)
    
    x = real(arg,dp)
    y = abs(x)
!
!  Evaluate erf for |X| <= 0.46875.
!
    if ( y .le. thresh ) then
       
       ysq = zero
       if ( xsmall .lt. y ) then
          ysq = y * y
       end if
       
       xnum = a(5) * ysq
       xden = ysq
       
       do i = 1, 3
          xnum = ( xnum + a(i) ) * ysq
          xden = ( xden + b(i) ) * ysq
       end do
       
       calerf = x * ( xnum + a(4) ) / ( xden + b(4) )
       
       if ( jint .ne. 0 ) then
          calerf = one - calerf
       end if
       
       if ( jint .eq. 2 ) then
          calerf = exp ( ysq ) * calerf
       end if
       
       return
!
!  Evaluate erfc for 0.46875 <= |X| <= 4.0.
!
    else if ( y .le. four ) then
       
       xnum = c(9) * y
       xden = y
       
       do i = 1, 7
          xnum = ( xnum + c(i) ) * y
          xden = ( xden + d(i) ) * y
       end do
       
       calerf = ( xnum + c(8) ) / ( xden + d(8) )
       
       if ( jint .ne. 2 ) then
          ysq = aint ( y * sixten ) / sixten
          del = ( y - ysq ) * ( y + ysq )
          calerf = exp ( -ysq * ysq ) * exp ( -del ) * calerf
       end if
!
!  Evaluate erfc for 4.0 < |X|.
!
    else
       
       calerf = zero
       
       if ( y .ge. xbig ) then
          
          if ( jint .ne. 2 .or. y .ge. xmax ) then
             go to 300
          end if
          
          if ( y .ge. xhuge ) then
             calerf = sqrpi / y
             go to 300
          end if
          
       end if
       
       ysq = one / ( y * y )
       xnum = p(6) * ysq
       xden = ysq
       do i = 1, 4
          xnum = ( xnum + p(i) ) * ysq
          xden = ( xden + q(i) ) * ysq
       end do
       
       calerf = ysq * ( xnum + p(5) ) / ( xden + q(5) )
       calerf = ( sqrpi -  calerf ) / y
       
       if ( jint .ne. 2 ) then
          ysq = aint ( y * sixten ) / sixten
          del = ( y - ysq ) * ( y + ysq )
          calerf = exp ( -ysq * ysq ) * exp ( -del ) * calerf
       end if
       
    end if
!
!  Fix up for negative argument, erf, etc.
!
300 continue
      
    if ( jint .eq. 0 ) then
       
       calerf = ( half - calerf ) + half
       if ( x .lt. zero ) then
          calerf = -calerf
       end if
       
    else if ( jint .eq. 1 ) then
       
       if ( x .lt. zero ) then
          calerf = two - calerf
       end if
       
    else
       
       if ( x .lt. zero ) then
          
          if ( x .lt. xneg ) then
             calerf = xinf
          else
             ysq = aint ( x * sixten ) / sixten
             del = ( x - ysq ) * ( x + ysq )
             y = exp ( ysq * ysq ) * exp ( del )
             calerf = ( y + y ) - calerf
          end if
          
       end if
       
    end if
    
    return
  end function calerf


  function gammp(a, x)
    implicit none

    real(dp), intent(in) :: a, x
    real(dp)             :: gammp

    if (x < a+1.d0) then
       gammp = gser(a,x)
    else
       gammp = 1.d0 - gcf(a,x)
    end if

  end function gammp
  
  function gser(a, x, gln)
    implicit none

    real(dp), intent(in)            :: a, x
    real(dp), intent(out), optional :: gln
    real(dp)                        :: gser

    integer(i4b), parameter :: ITMAX = 100
    real(dp), parameter :: eps = epsilon(x)
    integer(i4b) :: n
    real(dp)     :: ap, del, summ
    
    if (x == 0.d0) then
       gser = 0.d0
       return
    end if
    ap = a
    summ = 1.d0 / a
    del  = summ
    do n = 1, ITMAX
       ap   = ap + 1.d0
       del  = del * x/ap
       summ = summ + del
       if (abs(del) < abs(summ)*eps) exit
    end do
    if (n > ITMAX) then
       write(*,*) 'a too large, ITMAX too small in gser.'
       stop
    end if
    if (present(gln)) then
       gln = gammln_quiet(a)
       gser = summ*exp(-x+a*log(x)-gln)
    else
       gser = summ*exp(-x+a*log(x)-gammln_quiet(a))
    end if

  end function gser

!!$FUNCTION gser_s(a,x,gln)
!!$	USE nrtype; USE nrutil, ONLY : nrerror
!!$	USE nr, ONLY : gammln
!!$	IMPLICIT NONE
!!$	REAL(SP), INTENT(IN) :: a,x
!!$	REAL(SP), OPTIONAL, INTENT(OUT) :: gln
!!$	REAL(SP) :: gser_s
!!$	INTEGER(I4B), PARAMETER :: ITMAX=100
!!$	REAL(SP), PARAMETER :: EPS=epsilon(x)
!!$	INTEGER(I4B) :: n
!!$	REAL(SP) :: ap,del,summ
!!$	if (x == 0.0) then
!!$		gser_s=0.0
!!$		RETURN
!!$	end if
!!$	ap=a
!!$	summ=1.0_sp/a
!!$	del=summ
!!$	do n=1,ITMAX
!!$		ap=ap+1.0_sp
!!$		del=del*x/ap
!!$		summ=summ+del
!!$		if (abs(del) < abs(summ)*EPS) exit
!!$	end do
!!$	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
!!$	if (present(gln)) then
!!$		gln=gammln(a)
!!$		gser_s=summ*exp(-x+a*log(x)-gln)
!!$	else
!!$		gser_s=summ*exp(-x+a*log(x)-gammln(a))
!!$	end if
!!$	END FUNCTION gser_s



FUNCTION gcf(a,x,gln)
!	USE nrtype; USE nrutil, ONLY : nrerror
!	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: a,x
	REAL(dp), OPTIONAL, INTENT(OUT) :: gln
	REAL(dp) :: gcf
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(dp), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(dp) :: an,b,c,d,del,h
	if (x == 0.d0) then
		gcf=1.d0
		RETURN
	end if
	b=x+1.0_dp-a
	c=1.0_dp/FPMIN
	d=1.0_dp/b
	h=d
	do i=1,ITMAX
		an=-i*(i-a)
		b=b+2.0_dp
		d=an*d+b
		if (abs(d) < FPMIN) d=FPMIN
		c=b+an/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_dp) <= EPS) exit
	end do
!	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf')
	if (present(gln)) then
		gln=gammln_quiet(a)
		gcf=exp(-x+a*log(x)-gln)*h
	else
		gcf=exp(-x+a*log(x)-gammln_quiet(a))*h
	end if
END FUNCTION gcf

!!$  function gcf(a, x, gln)
!!$    implicit none
!!$
!!$    real(dp), intent(in)            :: a, x
!!$    real(dp), intent(out), optional :: gln
!!$    real(dp)                        :: gcf
!!$
!!$    integer(i4b), parameter :: ITMAX = 100
!!$    real(dp), parameter :: eps = epsilon(x), FPMIN=tiny(x)/EPS
!!$    integer(i4b) :: i
!!$    real(dp)     :: an, b, c, d, del, h
!!$    if (x == 0.d0) then
!!$       gcf = 0.d0
!!$       return
!!$    end if
!!$    b = x + 1.d0 - a
!!$    c = 1.d0 / FPMIN
!!$    d = 1.d0 / b
!!$    h = d
!!$    do i = 1, ITMAX
!!$       an = -i*(i-a)
!!$       b = b + 2.d0
!!$       d = an*d+b
!!$       if (abs(d) < FPMIN) d = FPMIN
!!$       c = b + an/c
!!$       if (abs(c) < FPMIN) c = FPMIN
!!$       d = 1.d0 / d
!!$       del = d*c
!!$       h = h*del
!!$       if (abs(del-1.d0) <= eps) exit
!!$    end do
!!$    if (i > ITMAX) then
!!$       write(*,*) 'a too large, ITMAX too small in gcf'
!!$       stop
!!$    end if
!!$    if (present(gln)) then
!!$       gln = gammln_quiet(a)
!!$       gcf = exp(-x+a*log(x)-gln)*h
!!$    else
!!$       gcf = exp(-x+a*log(x)-gammln_quiet(a))*h
!!$    end if
!!$
!!$  end function gcf

  function erf2(x)
    real(dp), intent(in) :: x
    real(dp)             :: erf2

    erf2 = gammp(0.5d0, x**2)
    if (x < 0.d0) erf2 = -erf2
  end function erf2

  function gammln_quiet(xx)
    implicit none

    real(dp), intent(in) :: xx
    real(dp)             :: gammln_quiet

    integer(i4b) :: i
    real(dp) :: tmp, x, tot
    real(dp) :: stp = 2.5066282746310005_dp
    real(dp), dimension(6) :: coef = (/76.18009172947146_dp, &
         -86.50532032941677_dp,24.01409824083091_dp, &
         -1.231739572450155_dp,0.1208650973866179e-2_dp, &
         -0.5395239384953e-5_dp/)
    x=xx
    tmp=x+5.5_dp
    tmp=(x+0.5_dp)*log(tmp)-tmp
    tot = 0.d0
    do i = 1, size(coef)
       tot = tot + coef(i) / (x+i)
    end do
    gammln_quiet=tmp+log(stp*(1.000000000190015_dp+tot)/x)

  end function gammln_quiet

!!$FUNCTION gammln_s(xx)
!!$	USE nrtype; USE nrutil, ONLY : arth,assert
!!$	IMPLICIT NONE
!!$	REAL(SP), INTENT(IN) :: xx
!!$	REAL(SP) :: gammln_s
!!$	REAL(DP) :: tmp,x
!!$	REAL(DP) :: stp = 2.5066282746310005_dp
!!$	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
!!$		-86.50532032941677_dp,24.01409824083091_dp,&
!!$		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
!!$		-0.5395239384953e-5_dp/)
!!$	call assert(xx > 0.0, 'gammln_s arg')
!!$	x=xx
!!$	tmp=x+5.5_dp
!!$	tmp=(x+0.5_dp)*log(tmp)-tmp
!!$	gammln_s=tmp+log(stp*(1.000000000190015_dp+&
!!$		sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
!!$	END FUNCTION gammln_s


  FUNCTION chebft(a,b,n,func)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a,b
    INTEGER(I4B), INTENT(IN) :: n
    REAL(dp), DIMENSION(n) :: chebft
    INTERFACE
       FUNCTION func(x)
         use healpix_types
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: x
         REAL(dp), DIMENSION(size(x)) :: func
       END FUNCTION func
    END INTERFACE
    REAL(DP) :: bma,bpa
    REAL(DP), DIMENSION(n) :: theta
    bma=0.5_dp*(b-a)
    bpa=0.5_dp*(b+a)
    theta(:)=PI*arth(0.5_dp,1.0_dp,n)/n
    chebft(:)=matmul(cos(outerprod(arth(0.0_dp,1.0_dp,n),theta)), &
         func(real(cos(theta)*bma+bpa,dp)))*2.0_dp/n
  END FUNCTION chebft

  FUNCTION chebev(a,b,c,x)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: a,b,x
    REAL(DP), DIMENSION(:), INTENT(IN) :: c
    REAL(DP) :: chebev
    INTEGER(I4B) :: j,m
    REAL(DP) :: d,dd,sv,y,y2
    if ((x-a)*(x-b) > 0.0) then
       write(*,*) 'x not in range in chebev', x-a, x-b
       stop
    end if
    m=size(c)
    d=0.0
    dd=0.0
    y=(2.0_dp*x-a-b)/(b-a)
    y2=2.0_dp*y
    do j=m,2,-1
       sv=d
       d=y2*d-dd+c(j)
       dd=sv
    end do
    chebev=y*d-dd+0.5d0*c(1)
  END FUNCTION chebev
  

  FUNCTION arth(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth(k)=arth(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth(k)=arth(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth

  FUNCTION outerprod(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod




end module quiet_nr_mod

 
