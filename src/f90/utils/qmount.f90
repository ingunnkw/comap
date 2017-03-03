! Read in a pointing scatter file, apply a mount model and spit out
! the residuals.
program qmount
  use quiet_pointing_mod
  implicit none

  type point_data
     real(dp)          :: mjd, p0(2), p(3)
     character(len=32) :: obj
     integer(i4b)      :: mod, ces
  end type

  character(len=512) :: parfile
  character(len=512) :: arg, toks(2), toks2(8)
  type(quiet_mount)  :: mount
  integer(i4b)       :: ai, n, i
  real(dp)           :: val, vals(8), psi, p(3)
  type(point_data)   :: point

  call getarg(1, parfile)
  call initialize_pmac_mod(parfile, "")
  call initialize_quiet_pointing_mod(parfile, .false.)

  ! Set up mount model
  mount = mount_none
  do ai = 2, iargc()
     call getarg(ai, arg)
     call get_tokens(arg, "=", toks, num=n, maxnum=size(toks))
     call assert(n == 2, "Invalid parameter: " // trim(arg))
     call get_tokens(toks(2), ",", toks2, num=n, maxnum=size(toks2))
     do i = 1, n; read(toks2(i),*) vals(i); end do
     vals = vals*DEG2RAD ! default is degrees
     select case(toks(1))
        case("col")
           if(n >= 1)   mount%theta_c = vals(1)
           if(n >= 2)   mount%psi_c   = vals(2)
        case("denc");   mount%enc_offset(:n) = vals(:n)
        case("daz");    mount%enc_offset(1)  = vals(1)
        case("del");    mount%enc_offset(2)  = vals(1)
        case("ddk");    mount%enc_offset(3)  = vals(1)
        case("flex");   mount%kf             = vals(1)
        case("atilt")
           if(n >= 1)   mount%theta          = vals(1)
           if(n >= 2)   mount%omega          = vals(2)
        case("etilt");  mount%theta_E        = vals(1)
        case("fp");     mount%fp_flex(:n)    = vals(:n)
        !case("test");   mount%test(:n)       = vals(:n)
        case("azcurr"); mount%azcurr       = vals(1)
        case default
           call assert(.false., "Unknown parameter: " // trim(toks(1)))
     end select
  end do
  call set_mount_override(.true., mount_model=mount)

  ! Transform points
  do
     read(*,*,end=1) point
     point%p  = point%p  * DEG2RAD
     point%p0 = point%p0 * DEG2RAD
     call swap_coordinate_convention(point%p(1),  point%p(2),  point%p(3), coord_tele)
     call swap_coordinate_convention(point%p0(1), point%p0(2), psi,        coord_hor)
     call coord_convert(coord_tele, point%p(1), point%p(2), point%p(3), &
      & coord_hor, p(1), p(2), p(3), mod=point%mod, mjd=point%mjd)
     write(*,'(f13.7,5f9.5,3e15.7,a10,2i8)') point%mjd, p(1:2), &
      & point%p0, point%p(3), (p(1:2)-point%p0(1:2))*RAD2DEG*60, &
      & (p(1)-point%p0(1))*RAD2DEG*60*sin(p(2)), &
      & " "//point%obj, point%ces, point%mod
  end do
  1 continue

end program
