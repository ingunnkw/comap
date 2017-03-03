! Uses the quiet pointing model to convert the (mjd,phi,theta,psi) coordinates
! passed as input from one system to another.
!
! Usage: [options] from to
program coord_convert_prog
  use quiet_utils
  use quiet_pointing_mod
  use powell_mod
  implicit none
  character(len=512) :: str, acc
  character(len=2048) :: parfile
  integer(i4b)       :: i, j, sys1, sys2, imod, diode, unit
  real(dp)           :: mjd, phi1, theta1, psi1, phi2, theta2, psi2
  logical(lgt)       :: ideg, odeg, iheal, oheal, posa, print_lst, model, invdk
  real(dp)           :: powell_dk_real(4)
  integer(i4b)       :: powell_dk_int(4)

  j = 0
  ideg = .false.
  odeg = .false.
  iheal = .false.
  oheal = .false.
  posa = .false.
  print_lst = .false.
  model = .true.
  invdk = .false.
  diode = -1
  imod = -1
  acc = "13.7"
  parfile = '/projects/quiet/auxilliary/std_point_param.txt'
  unit = getlun()
  i = 1
  do while(i <= iargc())
     call getarg(i, str)
     if(str == "-d") then; ideg = .true.; odeg = .true.
     elseif(str == "-id") then; ideg = .true.
     elseif(str == "-od") then; odeg = .true.
     elseif(str == "-r") then; ideg = .false.; odeg = .false.
     elseif(str == "-ir") then; ideg = .false.
     elseif(str == "-or") then; odeg = .false.
     elseif(str == "-m") then; iheal = .true.; oheal = .true.
     elseif(str == "-im") then; iheal = .true.
     elseif(str == "-om") then; oheal = .true.
     elseif(str == "-s") then; iheal = .false.; oheal = .false.
     elseif(str == "-is") then; iheal = .false.
     elseif(str == "-os") then; oheal = .false.
     elseif(str == "-invdk") then; invdk = .true.
     elseif(str == "-p") then; posa = .true.
     elseif(str == "-lst") then; print_lst = .true.
     elseif(str == "-nomount") then; model = .false.
     elseif(str == "-mod") then
        i = i+1
        call getarg(i, str)
        read(str,*) imod
     elseif(str == "-di") then
        i = i+1
        call getarg(i, str)
        read(str,*) diode
     elseif(str == "-params") then
        i = i+1
        call getarg(i, parfile)
     elseif(str == "-fmt") then
        i = i+1
        call getarg(i, acc)
     elseif(str == "-h") then
        call help
     else
        if(j == 0) then;     sys1 = parse_coord_name(str)
        elseif(j == 1) then; sys2 = parse_coord_name(str)
        end if
        j = j+1
     end if
     i = i+1
  end do
  if(j < 2) call help
  call initialize_quiet_pointing_mod(parfile, model)

  ! Initialization done. Now process the actual data on standard input
  do
     read(*,*,end=1) mjd, phi1, theta1, psi1
     if(ideg) then; phi1=phi1*DTOR; theta1=theta1*DTOR; psi1=psi1*DTOR; end if
     if(.not. iheal) call swap_coordinate_convention(phi1, theta1, psi1, sys1)
     if(.not. invdk) then
        call coord_convert(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, imod, diode)
     else
        psi2 = psi1
        call cc_invdk(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, imod, diode)
     end if
     if(.not. oheal) call swap_coordinate_convention(phi2, theta2, psi2, sys2)
     if(posa .and. phi2 < 0) phi2 = phi2+2*pi
     if(odeg) then; phi2=phi2*RTOD; theta2=theta2*RTOD; psi2=psi2*RTOD; end if
     write(*,fmt="(f20.12,3f" // trim(acc) // ")",advance="no") mjd, phi2, theta2, psi2
     if(print_lst) write(*,fmt="(f"//trim(acc)//")",advance="no") mjd2lst(mjd, QUIET_GEODETIC_LONGITUDE)
     write(*,*)
  end do
1 continue

contains
  subroutine help
     write(*,*) "Usage: [options] sys1 sys2. Pass coorinate sets on standard input."
     write(*,*) "Options are:"
     write(*,*) " -id, -od, -d: Coordinates in degrees in input, output, both."
     write(*,*) " -im, -om, -m: Mathematical (or healpix) convention in input, output, both."
     write(*,*) " -mod: Module in question. Only necessary for some conversions."
     write(*,*) " -di: Diode in question. Only necessary for some conversions."
     write(*,*) " -params: Location of parameter file for quiet pointing mod."
     write(*,*) " -h: Display this help."
     stop
  end subroutine
  ! Like coord convert, but the sys2 psi angle is known, while the sys1 psi angle isn't
  subroutine cc_invdk(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, mod, di)
    implicit none
    integer(i4b) :: sys1, sys2, mod, err, di
    real(dp)     :: phi1, phi2, theta1, theta2, psi1, psi2, mjd, p(1), foo
    p = 0
    powell_dk_real = [ phi1, theta1, psi2, mjd ]
    powell_dk_int  = [ sys1, sys2, mod, di ]
    call powell(p, powell_calc_dk, err)
    psi1 = p(1)
    call coord_convert(sys1, phi1, theta1, psi1, sys2, phi2, theta2, foo, mjd, mod, di)
  end subroutine

  function powell_calc_dk(p) result(res)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: foo(2), psi2, res
    call coord_convert(powell_dk_int(1), powell_dk_real(1), powell_dk_real(2), &
     & p(1), powell_dk_int(2), foo(1), foo(2), psi2, powell_dk_real(4), powell_dk_int(3), powell_dk_int(4))
    res = ang_diff(psi2, powell_dk_real(3))**2
  end function
end program
