! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001
module ziggurat
   implicit none
   private

   type zig_table
      integer  :: kn(0:127), ke(0:255)
      real*8   :: wn(0:127), fn(0:127), we(0:255), fe(0:255)
      logical  :: initialized = .false.
   end type

   type zig_rng
      integer  :: jsr
   end type

   real*8, parameter  :: m1=2147483648.d0,   m2=2147483648.d0

   type(zig_table), save :: table

   interface zig_int
      module procedure shr3
   end interface

   public  :: zig_rng, zig_init, zig_int, zig_uni, zig_gauss, zig_exp, zig_init_tables

contains

  subroutine zig_init(zig, seed)
    implicit none
    type(zig_rng) :: zig
    integer       :: seed
     if(.not. table%initialized) call zig_init_tables
    zig%jsr = seed
  end subroutine

  subroutine zig_init_tables
     implicit none
     integer            :: i
     real*8             :: q, dn, tn, vn, de, te, ve

     if(table%initialized) return
     !$OMP CRITICAL

     dn=3.442619855899d0;      tn=3.442619855899d0
     vn=0.00991256303526217d0; de=7.697117470131487d0
     te=7.697117470131487d0;   ve=0.003949659822581572d0

     !  Tables for RNOR
     q = vn*exp(0.5d0*dn**2)
     table%kn(0) = (dn/q)*m1
     table%kn(1) = 0
     table%wn(0) = q/m1
     table%wn(127) = dn/m1
     table%fn(0) = 1.d0
     table%fn(127) = exp(-0.5d0*dn**2)
     do i = 126, 1, -1
        dn = sqrt(-2.d0 * log(vn/dn+exp(-0.5d0*dn**2)))
        table%kn(i+1) = (dn/tn)*m1
        tn = dn
        table%fn(i) = exp(-0.5d0*dn**2)
        table%wn(i) = dn/m1
     end do

     !  tables for rexp
     q = ve*exp(de)
     table%ke(0) = (de/q)*m2
     table%ke(1) = 0
     table%we(0) = q/m2
     table%we(255) = de/m2
     table%fe(0) = 1.d0
     table%fe(255) = exp(-de)
     do i = 254, 1, -1
        de = -log(ve/de + exp(-de))
        table%ke(i+1) = m2 * (de/te)
        te = de
        table%fe(i) = exp(-de)
        table%we(i) = de/m2
     end do
     table%initialized = .true.
     !$OMP END CRITICAL
  end subroutine

  !  Generate random 32-bit integers
  function shr3(zig) result(ival)
     implicit none
     type(zig_rng) :: zig
     integer       ::  ival, jz
     jz = zig%jsr
     zig%jsr = ieor(zig%jsr, ishft(zig%jsr, 13))
     zig%jsr = ieor(zig%jsr, ishft(zig%jsr,-17))
     zig%jsr = ieor(zig%jsr, ishft(zig%jsr,  5))
     ival = jz + zig%jsr
  end function

  !  generate uniformly distributed random numbers
  function zig_uni(zig) result(fn_val)
     implicit none
     type(zig_rng) :: zig
     real*8           :: fn_val
     fn_val = 0.5d0 + 0.2328306d-9 * shr3(zig)
  end function

  !  generate random normals
  function zig_gauss(zig) result(fn_val)
     implicit none
     type(zig_rng)        :: zig
     real*8             :: fn_val, x, y
     integer              :: hz, iz
     real*8, parameter  :: r = 3.44262d0

     if(.not. table%initialized) call zig_init_tables
     hz = shr3(zig)
     iz = iand(hz, 127)
     if(abs(hz) < table%kn(iz)) then
        fn_val = hz * table%wn(iz)
     else
        do
           if(iz == 0) then
              do
                 x = -0.2904764d0* log(zig_uni(zig))
                 y = -log(zig_uni(zig))
                 if(y+y >= x*x) exit
              end do
              fn_val = r+x
              if(hz <= 0) fn_val = -fn_val
              return
           end if
           x = hz * table%wn(iz)
           if(table%fn(iz) + zig_uni(zig)*(table%fn(iz-1)-table%fn(iz)) < exp(-0.5d0*x**2)) then
              fn_val = x
              return
           end if
           hz = shr3(zig)
           iz = iand(hz, 127)
           if(abs(hz) < table%kn(iz)) then
              fn_val = hz * table%wn(iz)
              return
           end if
        end do
     end if
  end function

  !  generate random exponentials
  function zig_exp(zig) result(fn_val)
     implicit none
     type(zig_rng) :: zig
     real*8        :: fn_val, x
     integer       :: jz, iz

     if(.not. table%initialized) call zig_init_tables
     jz = shr3(zig)
     iz = iand(jz, 255)
     if(abs(jz) < table%ke(iz)) then
        fn_val = abs(jz) * table%we(iz)
        return
     end if
     do
        if(iz == 0) then
           fn_val = 7.69711 - log(zig_uni(zig))
           return
        end if
        x = abs(jz) * table%we(iz)
        if(table%fe(iz) + zig_uni(zig)*(table%fe(iz-1) - table%fe(iz)) < exp(-x)) then
           fn_val = x
           return
        end if
        jz = shr3(zig)
        iz = iand(jz, 255)
        if(abs(jz) < table%ke(iz)) then
           fn_val = abs(jz) * table%we(iz)
           return
        end if
     end do
     return
  end function

end module
