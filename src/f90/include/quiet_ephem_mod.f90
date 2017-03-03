module quiet_ephem_mod
 use healpix_types
 use pix_tools
 implicit none
 integer(i4b), parameter :: &
  & EPHEM_MERCURY = 1, EPHEM_VENUS  = 2,  EPHEM_EARTH  = 3, EPHEM_MARS    = 4, &
  & EPHEM_JUPITER = 5, EPHEM_SATURN = 6,  EPHEM_URANUS = 7, EPHEM_NEPTUNE = 8, &
  & EPHEM_PLUTO   = 9, EPHEM_MOON   = 10, EPHEM_SUN    = 11

contains

 ! Provides Equ-coordinates of object as (phi,theta) in healpix angles.
 ! Use coord_convert later if you want another system.
 function ephem(object, mjd) result(pos)
   implicit none
   integer(i4b)           :: object
   real(dp)               :: mjd, pos(2), x(3)
   call ephem_equ_rect_c(object, mjd, x)
   call vec2ang(x, pos(2), pos(1))
 end function

 function ephem_dist(object, mjd) result(dist)
   implicit none
   integer(i4b)           :: object
   real(dp)               :: mjd, x(3), dist
   call ephem_equ_rect_c(object, mjd, x)
   dist = sqrt(sum(x**2))
 end function

 function ephem_ss(object, mjd) result(pos)
   implicit none
   integer(i4b)           :: object
   real(dp)               :: mjd, pos(2), x(3)
   call ephem_equ_rect_abs_c(object, mjd, x)
   call vec2ang(x, pos(2), pos(1))
 end function

 function name2eph(name) result(id)
   implicit none
   character(len=*) :: name
   integer(i4b)     :: id
   select case(name)
      case("sun");     id = EPHEM_SUN
      case("moon");    id = EPHEM_MOON
      case("mercury"); id = EPHEM_MERCURY
      case("venus");   id = EPHEM_VENUS
      case("earth");   id = EPHEM_EARTH
      case("mars");    id = EPHEM_MARS
      case("jupiter"); id = EPHEM_JUPITER
      case("saturn");  id = EPHEM_SATURN
      case("uranus");  id = EPHEM_URANUS
      case("neptune"); id = EPHEM_NEPTUNE
      case("pluto");   id = EPHEM_PLUTO
      case default;    id = 0
   end select
 end function
end module
