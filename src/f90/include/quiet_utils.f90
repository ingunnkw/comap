module quiet_utils
  use healpix_types
  use pix_tools
  use sort_utils
  use spline_1D_mod
  use quiet_system_mod
  use quiet_defs
  use udgrade_nr
  use math_tools
  implicit none


  type rotmatrix
     real(dp), dimension(3,3) :: m
  end type rotmatrix

  interface set_ordering
     module procedure set_ordering_dp, set_ordering_sp, set_ordering_i4b, set_ordering_pix
  end interface

  interface healok
     module procedure healok, healok_sp
  end interface

  interface healnot
     module procedure healnot, healnot_sp
  end interface

  interface ifel
     module procedure ifeli, ifelc
  end interface

  interface irange
     module procedure irange_one, irange_two, irange_three
  end interface

  interface quantile
     module procedure quantile_single, quantile_multi
  end interface

  interface dump_matrix
     module procedure dump_matrix_mat, dump_matrix_vec
  end interface

  integer(i4b), parameter, private :: benchmax = 1024, benchomp = 32
  type benchmarker
     character(len=32) :: name(benchmax)
     real(dp)          :: t(2,benchomp,benchmax), dt(benchomp,benchmax)
     integer(i4b)      :: n = 0
     integer(i8b)      :: nsamp(benchomp,benchmax)
  end type

  type dmem_data
     integer(i4b) :: id = 0, level = 1
  end type

  integer(i4b), parameter :: stderr   = 0, stdout = 6
  integer(i4b), parameter :: RING     = 1, NEST = 2

  type(dmem_data), private, save :: ddata

contains

  subroutine get_map2mask_from_mask(mask, map2mask, mask2map)
    implicit none

    real(dp),              dimension(0:),  intent(in) :: mask
    integer(i4b), pointer, dimension(:)               :: map2mask, mask2map

    integer(i4b) :: i, j, n, npix

    n    = count(mask == 1.d0)
    npix = size(mask)

    allocate(map2mask(0:npix-1))
    allocate(mask2map(n))
    map2mask = -1

    j = 1
    do i = 0, npix-1
       if (mask(i) == 1.d0) then
          mask2map(j) = i
          map2mask(i) = j
          j           = j+1
       end if
    end do

  end subroutine get_map2mask_from_mask

!!$  ! unit argument is here for compatibility only. It is ignored.
!!$  ! This is an improved version of the old get_parameter. The main
!!$  ! part is the same, but it now uses a stack of unit numbers
!!$  ! and filenames to handle recursion into upto maxdepth levels of
!!$  ! include files. The filename stack is only needed for error
!!$  ! messages, and takes up 16kB of memory by default, but I don't
!!$  ! think that is a problem.
!!$  ! I am not fond of the way the open and read functions force
!!$  ! a goto-like structure on this subroutine.
!!$  subroutine get_parameter(unit, parfile, parname, par_int, par_char, &
!!$       & par_string, par_sp, par_dp, par_lgt, par_present)
!!$    implicit none
!!$    integer(I4b)               :: unit
!!$    character(len=*)           :: parfile, parname
!!$    integer(i4b),     optional :: par_int
!!$    character(len=1), optional :: par_char
!!$    character(len=*), optional :: par_string
!!$    real(sp),         optional :: par_sp
!!$    real(dp),         optional :: par_dp
!!$    logical(lgt),     optional :: par_lgt
!!$    logical(lgt),     optional :: par_present
!!$
!!$    integer(i4b), parameter    :: maxdepth = 256
!!$    integer(i4b)               :: depth, units(maxdepth), i
!!$    character(len=256)         :: key, value, filenames(maxdepth)
!!$    character(len=1)           :: equals
!!$
!!$    depth = 1
!!$    units(depth) = getlun()
!!$    !write(*,*) "Entering file " // trim(parfile)
!!$    filenames(depth) = parfile
!!$    open(units(depth),file=parfile,status="old",err=4)
!!$    do while(depth >= 1)
!!$       read(units(depth),*,end=1) key
!!$       if (key(1:1)=='#') cycle
!!$       backspace(units(depth))
!!$
!!$       if (key(1:1)=='@') then
!!$          if(key == '@INCLUDE') then
!!$             ! Recurse into the new file
!!$             read(units(depth),*,end=1) key, value
!!$             !write(*,*) "Entering file " // trim(value)
!!$             depth=depth+1
!!$             units(depth) = getlun()
!!$             filenames(depth) = value
!!$             open(units(depth),file=value,status="old",err=2)
!!$          else
!!$             goto 3
!!$          end if
!!$       else
!!$          read(units(depth),*) key, equals, value
!!$         if (trim(key) == trim(parname)) then
!!$             if (present(par_int)) then
!!$                read(value,*) par_int
!!$             else if (present(par_char)) then
!!$                read(value,*) par_char
!!$             else if (present(par_string)) then
!!$                read(value,'(a)') par_string
!!$             else if (present(par_sp)) then
!!$                read(value,*) par_sp
!!$             else if (present(par_dp)) then
!!$                read(value,*) par_dp
!!$             else if (present(par_lgt)) then
!!$                read(value,*) par_lgt
!!$             else
!!$                write(*,*) "get_parameter: Reached unreachable point!"
!!$             end if
!!$
!!$             ! Match found, so clean up and return.
!!$             do i = depth, 1, -1; close(units(i)); end do
!!$             if(present(par_present)) par_present = .true.
!!$             return
!!$          end if
!!$       end if
!!$       cycle
!!$       ! We get here if we reached the end of a file. Close it and
!!$       ! return to the file above.
!!$1      close(units(depth))
!!$       !write(*,*) "Exiting file " // filenames(depth)
!!$       depth = depth-1
!!$    end do
!!$
!!$    ! ===== Error handling section ======
!!$    ! Case 1: Failed to find matching parameter in file
!!$    if(present(par_present)) then
!!$       par_present = .false.
!!$       return
!!$    else
!!$       write(*,*) "get_parameter: Fatal error: Cannot find parameter " // trim(parname) // " in parameter file " // trim(parfile) // " or included files."
!!$       stop
!!$    end if
!!$
!!$    ! Case 2: Include file error
!!$2   write(*,*) "get_parameter: Fatal error: Cannot open include file '" // trim(value) // "'"
!!$    write(*,*) " in file " // trim(filenames(depth))
!!$    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
!!$    do i = depth-1, 1, -1; close(units(i)); end do
!!$    stop
!!$    ! Case 3: Directive error
!!$3   write(*,*) "get_parameter: Fatal error: Unrecognized directive '" // trim(key) //"'"
!!$    write(*,*) " in file " // trim(filenames(depth))
!!$    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
!!$    do i = depth, 1, -1; close(units(i)); end do
!!$    stop
!!$
!!$    ! Case 4: Top level parameter file unreadable
!!$4   write(*,*) "get_parameter: Fatal error: Cannot open parameter file '" // trim(parfile) // "'"
!!$    stop
!!$  end subroutine get_parameter






  ! Small utility for converting an integer to a string
  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string

  subroutine get_parameter(unit, parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    integer(I4b)               :: unit
    character(len=*)           :: parfile, parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found

    found = .false.
    call get_parameter_arg(parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, found, desc)
    if(found) then
       if(present(par_present)) par_present = .true.
    else
       call get_parameter_parfile(unit, parfile, parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt, par_present, desc)
    end if
  end subroutine

  subroutine get_parameter_parfile(unit, parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    integer(i4b)               :: unit
    character(len=*)           :: parfile, parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    logical(lgt)               :: found
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i
    character(len=512)         :: key, value, filenames(maxdepth), line
    character(len=1)           :: equals

    depth = 1
    units(depth) = getlun()
    !write(*,*) "Entering file " // trim(parfile)
    filenames(depth) = parfile
    open(units(depth),file=trim(parfile),status="old",err=4)
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       if (key(1:1)=='#') cycle
       backspace(units(depth))

       if (key(1:1)=='@') then
          if(key == '@INCLUDE') then
             ! Recurse into the new file
             read(units(depth),*,end=1) key, value
             !write(*,*) "Entering file " // trim(value)
             depth=depth+1
             units(depth) = getlun()
             filenames(depth) = value
             open(units(depth),file=value,status="old",err=2)
          else
             goto 3
          end if
       else
          read(units(depth),fmt="(a)") line
          call parse_parameter(line, parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
          if(found) then
             ! Match found, so clean up and return.
             do i = depth, 1, -1; close(units(i)); end do
             if(present(par_present)) par_present = .true.
             return
          end if
       end if
       cycle
       ! We get here if we reached the end of a file. Close it and
       ! return to the file above.
1      close(units(depth))
       !write(*,*) "Exiting file " // filenames(depth)
       depth = depth-1
    end do

    ! ===== Error handling section ======
    ! Case 1: Failed to find matching parameter in file
    if(present(par_present)) then
       par_present = .false.
       return
    else
       write(*,*) "get_parameter: Fatal error: Cannot find parameter " // trim(parname) &
            & // " in parameter file " // trim(parfile) // " or included files."
       if(present(desc)) write(*,*) trim(desc)
       stop
    end if

    ! Case 2: Include file error
2   write(*,*) "get_parameter: Fatal error: Cannot open include file '" // trim(value) // "'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth-1, 1, -1; close(units(i)); end do
    stop

    ! Case 3: Directive error
3   write(*,*) "get_parameter: Fatal error: Unrecognized directive '" // trim(key) //"'"
    write(*,*) " in file " // trim(filenames(depth))
    do i = depth-1, 1, -1; write(*,*) " included from " // trim(filenames(i)); end do
    do i = depth, 1, -1; close(units(i)); end do
    stop

    ! Case 4: Top level parameter file unreadable
4   write(*,*) "get_parameter: Fatal error: Cannot open parameter file '" // trim(parfile) // "'"
    stop
  end subroutine

  subroutine get_parameter_arg(parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present, desc)
    implicit none
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present
    character(len=*), optional :: desc

    character(len=512) :: line
    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, iargc()
       call getarg(i, line)
       if(line(1:2) /= "--") cycle
       call parse_parameter(line(3:), parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
       if(found) then
          if(present(par_present)) par_present = .true.
          return
       end if
    end do
    if(present(par_present)) then
       par_present = .false.
    else
       write(*,*) "get_parameter_arg: Fatal error: Cannot find parameter " // trim(parname) // " in argument list!"
       if(present(desc)) write(*,*) trim(desc)
       stop
    end if
  end subroutine

  subroutine get_parameter_arr(arr, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt, par_present)
    implicit none
    character(len=*)           :: arr(:)
    character(len=*)           :: parname
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt
    logical(lgt),     optional :: par_present

    integer(i4b)       :: i
    logical(lgt)       :: found
    do i = 1, size(arr)
       call parse_parameter(arr(i), parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
       if(found) then
          if(present(par_present)) par_present = .true.
          return
       end if
    end do
    if(present(par_present)) then
       par_present = .false.
    else
       write(*,*) "get_parameter_arr: Fatal error: Cannot find parameter " // trim(parname) // " in argument list!"
       stop
    end if
  end subroutine

  subroutine parse_parameter(line, parname, found, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
    implicit none
    character(len=*)           :: line, parname
    character(len=256)         :: toks(2), key, value, par
    logical(lgt)               :: found
    integer(i4b),     optional :: par_int
    character(len=*), optional :: par_char
    character(len=*), optional :: par_string
    real(sp),         optional :: par_sp
    real(dp),         optional :: par_dp
    logical(lgt),     optional :: par_lgt

    integer(i4b) :: i, n

    call get_tokens(trim(line), "=", group="''" // '""', maxnum=2, toks=toks, num=n)
    if(n < 2) then
       found = .false.
       return
    end if
    key = get_token(toks(1), "	 ", 1, group="''" // '""')
    value = get_token(toks(2), "	 ", 1, group="''" // '""')
    par = parname
    call tolower(key)
    call tolower(par)
    if (trim(key) == trim(par)) then
       if (present(par_int)) then
          read(value,*) par_int
       elseif (present(par_char)) then
          read(value,*) par_char
       elseif (present(par_string)) then
          read(value,*) par_string
       elseif (present(par_sp)) then
          read(value,*) par_sp
       elseif (present(par_dp)) then
          read(value,*) par_dp
       elseif (present(par_lgt)) then
          read(value,*) par_lgt
       else
          write(*,*) "get_parameter: Reached unreachable point!"
       end if
       found = .true.
    else
       found = .false.
    end if
  end subroutine

  ! Loops through the parameter files and children, counting lines.
  ! No error reporting.
  subroutine dump_expanded_paramfile(parfile, outfile)
    implicit none
    character(len=*)           :: parfile, outfile
    integer(i4b), parameter    :: maxdepth = 256
    integer(i4b)               :: depth, units(maxdepth), i, num, ounit
    character(len=1024)        :: key, value, arg
    character(len=1)           :: equals

    num = 0
    depth = 1
    ounit = getlun()
    open(ounit,file=outfile,action="write")
    write(ounit,fmt="(a)",advance="no") '# Arguments:'
    do i = 1, iargc()
       call getarg(i, arg)
       write(ounit,fmt="(a)",advance="no") " '" // trim(arg) // "'"
    end do
    write(ounit,*)

    units(depth) = getlun()
    open(units(depth),file=parfile,status="old",action="read")
    do while(depth >= 1)
       read(units(depth),*,end=1) key
       backspace(units(depth))
       if (key=='@INCLUDE') then
          ! Recurse into the new file
          read(units(depth),*,end=1) key, value
          write(ounit,fmt='(a)') pad("",depth-1," ") // "# File: " // trim(value)
          depth=depth+1
          units(depth) = getlun()
          open(units(depth),file=value,status="old")
       else
          read(units(depth),fmt="(a)") value
          write(ounit,fmt='(a)') pad("",depth-1," ") // trim(value)
       end if
       cycle
1      close(units(depth))
       depth = depth-1
    end do
    close(ounit)
  end subroutine

  ! The special filenames
  subroutine read_module_list(unit, filename, module_list)
    implicit none

    integer(i4b),                            intent(in) :: unit
    character(len=*),                        intent(in) :: filename
    integer(i4b),     pointer, dimension(:)             :: module_list

    integer(i4b) :: n, mod
    logical(lgt) :: exist

    inquire(file=trim(filename), exist=exist)
    if (.not. exist) then
       write(*,*) 'ERROR: Module list does not exist = ', trim(filename)
       stop
    end if

    n = 0
    open(unit, file=trim(filename))
    do while (.true.)
       read(unit,*,end=1)
       n = n+1
    end do
1   close(unit)

    allocate(module_list(n))
    open(unit, file=trim(filename))
    n = 1
    do while (.true.)
       read(unit,*,end=2) mod
       module_list(n) = mod
       n = n+1
    end do
2   close(unit)

  end subroutine read_module_list


  subroutine convert_EM2ang_zyz(EM, phi, theta, psi)
    implicit none

    real(dp), dimension(3,3), intent(in)   :: EM
    real(dp),                 intent(out)  :: phi, theta, psi

    theta = acos( EM(3,3))
    psi   = atan2(EM(3,2), -EM(3,1))
    phi   = atan2(EM(2,3),  EM(1,3))

  end subroutine convert_EM2ang_zyz

  subroutine convert_EM2ang(EM, phi, theta, psi)
    implicit none

    real(dp), dimension(3,3), intent(in)   :: EM
    real(dp),                 intent(out)  :: phi, theta, psi

    theta = acos( EM(3,3))
    phi   = atan2(EM(3,1), -EM(3,2))
    psi   = atan2(EM(1,3),  EM(2,3))

  end subroutine convert_EM2ang

  ! Convention: First phi around z, then theta around x, then psi around z
  subroutine compute_euler_matrix(phi, theta, psi, euler_matrix)
    implicit none

    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix

    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi

    sphi = sin(phi)
    cphi = cos(phi)

    sth  = sin(theta)
    cth  = cos(theta)

    spsi = sin(psi)
    cpsi = cos(psi)

    euler_matrix(1,1) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(1,2) =  sphi * cpsi + cth * cphi * spsi
    euler_matrix(1,3) =                sth * spsi
    euler_matrix(2,1) = -cphi * spsi - cth * sphi * cpsi
    euler_matrix(2,2) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(2,3) =                sth * cpsi
    euler_matrix(3,1) =                sth * sphi
    euler_matrix(3,2) =              - sth * cphi
    euler_matrix(3,3) =                cth

  end subroutine compute_euler_matrix

  ! Convention: First psi around z, then theta around y, then phi around z
  subroutine compute_euler_matrix_zyz(phi, theta, psi, euler_matrix)
    implicit none

    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix

    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi

    sphi = sin(phi)
    cphi = cos(phi)

    sth  = sin(theta)
    cth  = cos(theta)

    spsi = sin(psi)
    cpsi = cos(psi)

    euler_matrix(1,1) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(1,2) = -sphi * cpsi - cth * cphi * spsi
    euler_matrix(1,3) =                sth * cphi
    euler_matrix(2,1) =  cphi * spsi + cth * sphi * cpsi
    euler_matrix(2,2) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(2,3) =                sth * sphi
    euler_matrix(3,1) =              - sth * cpsi
    euler_matrix(3,2) =                sth * spsi
    euler_matrix(3,3) =                cth

  end subroutine compute_euler_matrix_zyz

  subroutine convert_euler_matrix_to_angles_zyz(M, phi, theta, psi)
    implicit none

    real(dp), dimension(3,3), intent(in)  :: M
    real(dp),                 intent(out) :: phi, theta, psi

    real(dp) :: sth, cth, spsi, sphi

    ! Get theta
    cth = M(3,3)
    if (cth > 1.d0) then
       theta = 0.d0
    else if (cth < -1.d0) then
       theta = pi
    else
       theta = acos(cth)
    end if

    if (abs(cth) < 0.9999998d0) then
       phi = atan2(M(2,3),M(1,3))
       psi = -atan2(-M(3,2),-M(3,1))
    else
       ! Psi and phi are degenerate; define phi = 0
       phi = 0.d0
       psi = atan2(-M(1,2),M(2,2))
    end if

  end subroutine convert_euler_matrix_to_angles_zyz






  ! NB! Only one of the arguments may be non-zero!!
  subroutine compute_rotation_matrix(axis, angle, matrix)
    implicit none

    character(len=1),         intent(in)   :: axis
    real(dp),                 intent(in)   :: angle
    real(dp), dimension(3,3), intent(out)  :: matrix

    matrix = 0.d0

    if (axis == 'x') then

       matrix(1,1) = 1.d0

       matrix(2,2) =  cos(angle)
       matrix(2,3) = -sin(angle)
       matrix(3,2) =  sin(angle)
       matrix(3,3) =  cos(angle)

    else if (axis == 'y') then

       matrix(2,2) = 1.d0

       matrix(1,1) =  cos(angle)
       matrix(1,3) =  sin(angle)
       matrix(3,1) = -sin(angle)
       matrix(3,3) =  cos(angle)

    else if (axis == 'z') then

       matrix(3,3) = 1.d0

       matrix(1,1) =  cos(angle)
       matrix(1,2) = -sin(angle)
       matrix(2,1) =  sin(angle)
       matrix(2,2) =  cos(angle)

    else

!       write(UNIT_CRITICAL,*) 'SCAN_GEN:   Invalid argument to compute_rotation_matrix routine'
       write(*,*)             'SCAN_GEN:   Invalid argument to compute_rotation_matrix routine'

    end if

  end subroutine compute_rotation_matrix


  ! Routine for computing the ecliptic to galactic rotation matrix
  subroutine compute_ecl2gal_matrix(matrix)
    implicit none

    real(dp), dimension(3,3) :: matrix

    matrix(1,1) =  -0.054882486d0
    matrix(1,2) =  -0.993821033d0
    matrix(1,3) =  -0.096476249d0
    matrix(2,1) =   0.494116468d0
    matrix(2,2) =  -0.110993846d0
    matrix(2,3) =   0.862281440d0
    matrix(3,1) =  -0.867661702d0
    matrix(3,2) =  -0.000346354d0
    matrix(3,3) =   0.497154957d0

  end subroutine compute_ecl2gal_matrix


  ! This routine is more or less stolen from IDL. Have no idea about copyrights issues,
  ! but I hope it is OK since we all have an IDL license...
  ! Conversion types:
  !   1 = Celestial to Galactic
  !   2 = Galactic  to Celestial
  !   3 = Celestial to Ecliptic
  !   4 = Ecliptic  to Celestial
  !   5 = Ecliptic  to Galactic
  !   6 = Galactic  to Ecliptic
  subroutine convert_coord_idl(convtype, vector)
    implicit none

    integer(i4b),           intent(in)     :: convtype
    real(dp), dimension(3), intent(inout)  :: vector

    integer(i4b) :: i
    real(dp)     :: deg2rad, ai, bi, sb, cb, cbsa, bo, ao
    real(dp),              dimension(3,3)       :: tempmat
    real(dp), allocatable, dimension(:,:), save :: rotmat
    real(dp),              dimension(6)         :: psi, stheta, ctheta, phi

    if (.not. allocated(rotmat)) then
       allocate(rotmat(3,3))

       psi = (/ 0.57477043300d0, 4.9368292465d0, &
             & 0.00000000000d0, 0.0000000000d0, &
             & 0.11142137093d0, 4.71279419371d0 /)
       stheta =(/ 0.88998808748d0,-0.88998808748d0, &
               & 0.39777715593d0,-0.39777715593d0, &
               & 0.86766622025d0,-0.86766622025d0 /)
       ctheta = (/ 0.45598377618d0, 0.45598377618d0, &
               & 0.91748206207d0, 0.91748206207d0, &
               & 0.49714719172d0, 0.49714719172d0 /)
       phi  = (/ 4.9368292465d0,  0.57477043300d0, &
              & 0.0000000000d0, 0.00000000000d0, &
              & 4.71279419371d0, 0.11142137093d0 /)

       deg2rad = pi / 180.d0

       ! Find each column of the rotation matrix
       do i = 1, 3

          if (i == 1) then
             ai = 0.d0
             bi = 0.d0
          else if (i == 2) then
             ai = 0.5d0 * pi
             bi = 0.d0
          else
             ai = 0.d0
             bi = 0.5d0 * pi
          end if

          ao   = ai - phi(convtype)
          bo   = bi
          sb   = sin(bo)
          cb   = cos(bo)
          cbsa = cb * sin(ao)
          bo   = -stheta(convtype) * cbsa + ctheta(convtype) * sb
          if (bo < 1.d0) bo  = asin(bo)

          ao =  atan2( ctheta(convtype) * cbsa + stheta(convtype) * sb, cb * cos(ao) )
          ao = mod(ao+psi(convtype)+4.d0*pi, 2.d0*pi)

          rotmat(1,i) = cos(ao) * sin(0.5d0*pi-bo)
          rotmat(2,i) = sin(ao) * sin(0.5d0*pi-bo)
          rotmat(3,i) = cos(0.5d0*pi-bo)

       end do

    end if

!    do i = 1, 3
!       write(*,*) real(rotmat(i,:),sp)
!    end do

    vector = matmul(rotmat, vector)


  end subroutine convert_coord_idl


  ! This is very stupid -- somebody has to rewrite this...
  subroutine convert_cel2gal_coord(vector)
    implicit none

    real(dp), dimension(3), intent(inout)  :: vector

    integer(i4b) :: i
    real(dp)  :: azimuth
    real(dp)  :: theta_gal_pole = 1.09955 ! theta = 63 degrees, inaccurate
    real(dp)  :: phi_gal_pole   = 3.36849 ! phi   = 193 degrees, inaccurate
    real(dp)  :: theta_gal_cent = 2.07694 ! theta = 119 degrees, inaccurate
    real(dp)  :: phi_gal_cent   = 4.64257 ! phi   = 266 degrees, inaccurate
    real(dp),              dimension(3)         :: center
    real(dp),              dimension(3,3)       :: tempmat
    real(dp), allocatable, dimension(:,:), save :: rotmat

    if (.not. allocated(rotmat)) then
       allocate(rotmat(3,3))

       center(1) = cos(phi_gal_cent) * sin(theta_gal_cent)
       center(2) = sin(phi_gal_cent) * sin(theta_gal_cent)
       center(3) = cos(theta_gal_cent)

       call compute_rotation_matrix('z', phi_gal_pole,   rotmat)
       call compute_rotation_matrix('y', theta_gal_pole, tempmat)
       rotmat = matmul(tempmat, rotmat)

       center = matmul(rotmat, center)

       azimuth = atan2(center(2), center(1))

       call compute_rotation_matrix('z', azimuth, tempmat)

       rotmat = matmul(tempmat, rotmat)

    end if

!    vector(1) = cos(phi_gal_pole) * sin(theta_gal_pole)
!    vector(2) = sin(phi_gal_pole) * sin(theta_gal_pole)
!    vector(3) = cos(theta_gal_pole)

!    vector(1) = cos(phi_gal_cent) * sin(theta_gal_cent)
!    vector(2) = sin(phi_gal_cent) * sin(theta_gal_cent)
!    vector(3) = cos(theta_gal_cent)

    vector = matmul(rotmat, vector)


  end subroutine convert_cel2gal_coord


  subroutine project_to_2D(vector, vector_cent, x, y)
    implicit none

    real(dp), dimension(3), intent(in)   :: vector, vector_cent
    real(dp),               intent(out)  :: x, y

    real(dp) :: phi, theta, cos_alpha, sgn
    real(dp) :: angular_dist, cos_dist, M(3,3), v(3)

    phi = atan2(vector_cent(2), vector_cent(1))
    theta = acos(vector_cent(3))
    call angdist2(vector, vector_cent, angular_dist)
    
    call compute_euler_matrix_zyz(0.d0, -theta, -phi, M)
    v = matmul(M, vector)
    v(1:2) = v(1:2) / sqrt(sum(v(1:2)**2))
    x = v(1) * angular_dist * 180.d0/pi
    y = v(2) * angular_dist * 180.d0/pi

  end subroutine project_to_2D

  function crossprod(u, v) RESULT(w)
    real(dp), dimension(3), intent(in) :: u, v
    real(dp), dimension(3)             :: w

    w(1) = u(2) * v(3) - u(3) * v(2)
    w(2) = u(3) * v(1) - u(1) * v(3)
    w(3) = u(1) * v(2) - u(2) * v(1)
  end function crossprod


  subroutine Make_Axis_Rotation_Transform(axis, angle, M)
    implicit none

    real(dp), dimension(3),   intent(in)  :: axis
    real(dp),                 intent(in)  :: angle
    real(dp), dimension(3,3), intent(out) :: M

    integer(i4b) :: i
    real(dp)     :: sa, ca, ica, t1, t2

    sa  = sin(angle)
    ca  = cos(angle)
    ica = 1.d0 - ca

    M(1,1) = axis(1)**2 * ica + ca
    M(2,2) = axis(2)**2 * ica + ca
    M(3,3) = axis(3)**2 * ica + ca

    t1 = axis(1)*axis(2)*ica
    t2 = axis(3)*sa
    M(2,1) = t1 + t2
    M(1,2) = t1 - t2

    t1 = axis(1)*axis(3)*ica
    t2 = axis(2)*sa
    M(3,1) = t1 - t2
    M(1,3) = t1 + t2

    t1 = axis(2)*axis(3)*ica
    t2 = axis(1)*sa
    M(2,3) = t1 - t2
    M(3,2) = t1 + t2

  end subroutine Make_Axis_Rotation_Transform

  ! Like angdist, but takes [ theta, phi ] coordinates instead of
  ! rectangular ones.
  function polangdist(p1, p2) result(dist)
    implicit none
    real(dp) :: p1(2), p2(2), v1(3), v2(3), dist
    call ang2vec(p1(1),p1(2),v1)
    call ang2vec(p2(1),p2(2),v2)
    call angdist(v1,v2,dist)
  end function

  !=======================================================================
  subroutine angdist2(v1, v2, dist)
    !=======================================================================
    ! call angdist(v1, v2, dist)
    ! computes the angular distance dist (in rad) between 2 vectors v1 and v2
    ! in general dist = acos ( v1 . v2 )
    ! except if the 2 vectors are almost aligned.
    !=======================================================================
    real(kind=DP), intent(IN), dimension(1:) :: v1, v2
    real(kind=DP), intent(OUT) :: dist

    real(kind=DP), dimension(1:3) :: r1, r2, vdiff
    real(kind=DP) :: diff, sprod
    !=======================================================================

    ! normalize both vectors
    r1(1:3) = v1(1:3) / sqrt(dot_product(v1,v1))
    r2(1:3) = v2(1:3) / sqrt(dot_product(v2,v2))

    sprod = DOT_PRODUCT(r1, r2)

    if (sprod > 0.999_dp) then
       ! almost colinear vectors
       vdiff(1:3) = r1(1:3) - r2(1:3)
       diff = sqrt(dot_product(vdiff,vdiff)) ! norm of difference
       dist = 2.0_dp * asin(diff * 0.5_dp)

    else if (sprod < -0.999_dp) then
       ! almost anti-colinear vectors
       vdiff(1:3) = r1(1:3) + r2(1:3)
       diff = sqrt(dot_product(vdiff,vdiff)) ! norm of sum
       dist = PI - 2.0_dp * asin(diff * 0.5_dp)

    else
       ! other cases
       dist = acos( sprod )
    endif


    return
  end subroutine angdist2

  function is_nan(a) result(res)
    implicit none
    real(dp) :: a
    logical(lgt) :: res
    res = (a .ne. a) .or. (a .eq. NaN)
  end function

  function getlun()
    implicit none
    integer(i4b) :: getlun
    logical(lgt) :: exists, isopen
    getlun = 9
    do
       getlun = getlun+1
       inquire(unit=getlun,exist=exists)
       if(exists) then
          inquire(unit=getlun,opened=isopen)
          if(.not. isopen) return
       end if
    end do
  end function getlun

  ! This will be slightly slow if the compiler doesn't optimize the
  ! allocation and deallocation here. Could destroy input array if this
  ! is a problem.
  function median(array) result(res)
    implicit none
    real(dp) :: array(:), res
    real(dp), dimension(:), allocatable :: tmp

    allocate(tmp(size(array)))
    tmp = array
    call QuickSort_real(tmp)
    res = tmp(size(tmp)/2+1)
    deallocate(tmp)
  end function median

  function quantile_single(array, q) result(res)
    implicit none
    real(dp) :: array(:), res, q
    real(dp), dimension(:), allocatable :: tmp
    allocate(tmp(size(array)))
    tmp = array
    call QuickSort_real(tmp)
    res = tmp(min(floor(size(tmp)*q)+1,size(tmp)))
    deallocate(tmp)
  end function

  ! Getting more at once is more efficient
  function quantile_multi(array, q) result(res)
    implicit none
    real(dp) :: array(:), q(:), res(size(q))
    real(dp), dimension(:), allocatable :: tmp
    allocate(tmp(size(array)))
    tmp = array
    call QuickSort_real(tmp)
    res = tmp(min(floor(size(tmp)*q)+1,size(tmp)))
    deallocate(tmp)
  end function

  function mean(array) result(res)
    implicit none
    real(dp) :: array(:), res
    res = sum(array)/size(array)
  end function mean

  function variance(array) result(res)
    implicit none
    real(dp) :: array(:), res
    res = sum((array-mean(array))**2)/(size(array)-1)
  end function

  subroutine median_filter(iarr, oarr, width, thinning)
    implicit none
    real(dp) :: iarr(0:), oarr(1:)
    integer(i4b) :: width, row, j, n, s, a, b, step
    integer(i4b), dimension(:), optional :: thinning

    n = size(iarr)
    do j = 1, size(oarr)
       if(present(thinning)) then
          row = thinning(j)-1
       else
          row = j-1
       end if
       s = width
       if(s/2 > row) s = 2*row+1
       if(row+s-s/2 >= n) s = 2*(n-1-row)+1
       a = s/2
       b = s-a
       oarr(j) = median(iarr(row-a:row+b-1))
    end do
  end subroutine median_filter

  subroutine average_filter(iarr, oarr, width, thinning)
    implicit none
    real(dp) :: iarr(0:), oarr(1:)
    integer(i4b) :: width, row, j, n, s, a, b, step
    integer(i4b), dimension(:), optional :: thinning

    n = size(iarr)
    do j = 1, size(oarr)
       if(present(thinning)) then
          row = thinning(j)-1
       else
          row = j-1
       end if
       s = width
       if(s/2 > row) s = 2*row+1
       if(row+s-s/2 >= n) s = 2*(n-1-row)+1
       a = s/2
       b = s-a
       oarr(j) = sum(iarr(row-a:row+b-1))/size(iarr(row-a:row+b-1))
    end do
  end subroutine average_filter

  ! Given something like (/ 1, 2, 5, 9 /) and (/ 1, 0, 0, 1 /)
  ! returns (/ 1, 9 /). Remember to deallocate the result!
  subroutine listmask(list,mask,res)
    implicit none
    integer(i4b), intent(in), dimension(:) :: list
    logical(lgt), intent(in), dimension(:) :: mask
    integer(i4b), pointer,    dimension(:) :: res
    integer(i4b) :: i, j, n

    n = 0
    do i = 1, size(mask)
       if(mask(i)) n = n+1
    end do
    allocate(res(n))
    j = 1
    do i = 1, size(list)
       if(mask(i)) then
          res(j) = list(i)
          j = j+1
       end if
    end do
  end subroutine

  ! First a first order polynomial through data points given by x and y
  subroutine linsqr(x, y, a, b)
    implicit none
    real(dp) :: x(:), y(:), a, b, t1, t2, t3, t4, t5, n
    n  = size(x)
    t1 = sum(x)
    t2 = sum(x**2)
    t3 = sum(y)
    t4 = sum(x*y)
    t5 = t1**2-n*t2

    a = (t3*t1-n*t4)/t5
    b = (t1*t4-t3*t2)/t5
  end subroutine linsqr

  function word_index(num, string) result(res)
    implicit none
    integer(i4b) :: num, res, i, j
    character(len=*) :: string

    res = 1
    if(num <= 0) return
    j = 0
    do i = 2, len(string)
       if(string(i:i) == " " .and. string(i-1:i-1) .ne. " ") then
          j = j+1
          if(j == num) then
             res = i
             return
          end if
       end if
    end do
    res = len(string)
  end function

  function skipw(num, string) result(res)
    implicit none
    integer(i4b) :: num, i
    character(len=*) :: string
    character(len=len(string)) :: res
    res = " " ! Make res blank
    i = word_index(num, string)
    res(1:len(string)-(i-1)) = string(i:)
  end function

  subroutine tokenize(string, sep, ext, group, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional   :: group
    character(len=256)  :: op, cl
    integer(i4b), save           :: level(256), nl
    integer(i4b), intent(inout)  :: ext(2)
    logical(lgt), optional       :: allow_empty

    integer(i4b) :: i, j, k, n, o, c, ng
    logical(lgt) :: intok, hit, empty

    empty = .false.; if(present(allow_empty)) empty = allow_empty

    if(ext(2) >= len(string)) then
       ext = (/ 0, -1 /)
       return
    end if
    ng = 0
    if(present(group)) then
       ng = len_trim(group)/2
       do i = 1, ng
          op(i:i) = group(2*i-1:2*i-1)
          cl(i:i) = group(2*i:2*i)
       end do
    end if
    if(ext(2) <= 0) then
       level = 0
       nl = 0
    end if
    intok = .false.
    do i = ext(2)+2, len(string)
       hit = .false.
       c = index(cl(1:ng), string(i:i))
       if(c /= 0) then; if(level(c) > 0) then
          level(c) = level(c) - 1
          if(level(c) == 0) nl = nl - 1
          hit = .true.
       end if; end if
       if(nl == 0) then
          ! Are we in a separator or not?
          if(index(sep, string(i:i)) == 0) then
             ! Nope, so we must be in a token. Register start of token.
             if(.not. intok) then
                j = i
                intok = .true.
             end if
          else
             ! Yes. This either means that a token is done, and we should
             ! return it, or that we are waiting for a new token, in
             ! which case do nothing.
             if(intok) then
                ext = (/ j, i-1 /)
                return
             elseif(empty) then
                ext = (/ i, i-1 /)
                return
             end if
          end if
       end if
       o = index(op(1:ng), string(i:i))
       if(o /= 0 .and. .not. hit) then
          if(level(o) == 0) nl = nl + 1
          level(o) = level(o) + 1
       end if
    end do
    ! Handle last token
    if(intok) then
       ext = (/ j, i-1 /)
    elseif(empty) then
       ext = (/ i, i-1 /)
    else
       ext = (/ 0, -1 /)
    end if
  end subroutine

  function get_token(string, sep, num, group, allow_empty) result(res)
    implicit none
    character(len=*)           :: string, sep
    character(len=len(string)) :: res
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)               :: i, num, ext(2)
    ext = -1
    do i = 1, num; call tokenize(string, sep, ext, group, allow_empty); end do
    res = string(ext(1):ext(2))
  end function

  ! Fill all tokens into toks, and the num filled into num
  subroutine get_tokens(string, sep, toks, num, group, maxnum, allow_empty)
    implicit none
    character(len=*) :: string, sep
    character(len=*) :: toks(:)
    character(len=*), optional :: group
    integer(i4b),     optional :: num, maxnum
    logical(lgt),     optional :: allow_empty
    integer(i4b) :: n, ext(2), nmax
    ext = -1
    n = 0
    nmax = size(toks); if(present(maxnum)) nmax = maxnum
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0 .and. n < nmax)
       n = n+1
       toks(n) = string(ext(1):ext(2))
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    if(present(num)) num = n
  end subroutine

  function has_token(token, string, sep, group, allow_empty) result(res)
    implicit none
    character(len=*) :: token, string, sep
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    logical(lgt)     :: res
    integer(i4b)     :: ext(2)
    res = .true.
    ext = -1
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0)
       if(string(ext(1):ext(2)) == trim(token)) return
       call tokenize(string, sep, ext, group, allow_empty)
    end do
    res = .false.
  end function

  function num_tokens(string, sep, group, allow_empty) result(res)
    implicit none
    character(len=*) :: string, sep
    character(len=*), optional :: group
    logical(lgt),     optional :: allow_empty
    integer(i4b)     :: res, ext(2)
    ext = -1
    res = 0
    call tokenize(string, sep, ext, group, allow_empty)
    do while(ext(1) > 0)
       res = res+1
       call tokenize(string, sep, ext, group, allow_empty)
    end do
  end function

  ! Shifts the array arr steps steps to the right. Steps can be
  ! A non-integer, as this subroutine interpolates. Extrapolation
  ! is used on the end that is shifted into existance.
  subroutine shift_array_interpol(arr, steps)
    implicit none
    real(dp) :: arr(:), steps
    real(dp), dimension(:), allocatable :: x, y2, tmp
    integer(i4b) :: i, n
    n = size(arr)
    if(n == 0) return
    allocate(x(n))
    allocate(y2(n))
    allocate(tmp(n))
    do i = 1, n
       x(i) = i
    end do
    call spline(x, arr, 1d30, 1d30, y2)
    do i = 1, n
       tmp(i) = splint(x,arr,y2,x(i)+steps)
    end do
    arr = tmp
    deallocate(x,y2,tmp)
  end subroutine

  function index_of(a,b) result(i)
    implicit none
    integer(i4b) :: a(:), b, i
    do i = 1, size(a); if(a(i) == b) return; end do
    i = 0
  end function

  ! Intersection of integer lists a and b into c, which must
  ! be deallocated afterwards
  subroutine intersection(a,b,c)
    implicit none
    integer(i4b) :: a(:), b(:)
    integer(i4b), dimension(:), allocatable :: c
    integer(i4b) :: i, j, k, n
    n = 0
    do i = 1, size(a)
       if(index_of(b,a(i)) > 0) n = n+1
    end do
    allocate(c(n))
    n = 0
    do i = 1, size(a)
       if(index_of(b,a(i)) > 0) then
          n = n+1
          c(n) = a(i)
       end if
    end do
  end subroutine

  function atoi(str, ok) result(res)
    implicit none
    character(len=*)       :: str
    logical(lgt), optional :: ok
    logical(lgt)           :: ok_
    integer(i4b)           :: res
    ok_ = .true.
    read(str,fmt=*, err=2) res
    return
2   ok_ = .false.
    if(present(ok)) ok = ok_
  end function

  ! Convert int to ascii string
  function itoa(i, npad) result(a)
    implicit none
    integer(i4b) :: i
    integer(i4b), optional :: npad
    character(len=32) :: a, tmp, tmp2
    write(tmp,*) i
    tmp2 = trim(adjustl(tmp))
    if(present(npad)) then
       if(i < 0) then
          a = "-" // pad(tmp2(2:), npad, '0')
       else
          a = pad(tmp2, npad, '0')
       end if
    else
       a = tmp2
    end if
  end function

  function atof(str, ok) result(res)
    implicit none
    character(len=*) :: str
    logical(lgt), optional :: ok
    logical(lgt)           :: ok_
    real(dp)               :: res
    ok_ = .true.
    read(str,fmt=*, err=2) res
    return
2   ok_ = .false.
    if(present(ok)) ok = ok_
  end function

  ! Convert float to ascii string
  function ftoa(f, ndig, npad) result(a)
    implicit none
    real(dp) :: f
    integer(i4b), optional :: npad, ndig
    integer(i4b) :: dig
    character(len=64) :: a, tmp, tmp2
    dig = 5; if(present(ndig)) dig = ndig
    write(tmp,fmt="(f64." // trim(itoa(dig)) // ")") f
    tmp2 = trim(adjustl(tmp))
    if(present(npad)) then
       if(f < 0) then
          a = "-" // pad(tmp2(2:), npad, '0')
       else
          a = pad(tmp2, npad, '0')
       end if
    else
       a = tmp2
    end if
  end function

  function pad(a, n, c) result(b)
    implicit none
    character(len=*) :: a
    character(len=1), optional :: c
    character(len=1):: k
    integer(i4b)    :: n, i, m
    character(len=max(len_trim(a),n)) :: b
    if(len_trim(a) > n) then
       b = a
    else
       k = '0'; if(present(c)) k = c
       m = len_trim(a)
       do i = 1, n-m; b(i:i) = k; end do
       b(n-m+1:n) = a(1:m)
    end if
  end function

  !---------------------------------------------------------------------------
  ! converts from anthenna temperature to thermodynamic, for given frequency
  !--------------------------------------------------------------------------
  function ant2thermo(frequency)
    implicit none

    real(dp), intent(in)  :: frequency
    real(dp) :: ant2thermo, x
    real(dp) :: k_B     = 1.3806503d-23
    real(dp) :: h       = 6.626068d-34
    real(dp) :: c       = 2.99792458d8
    real(dp) :: T_0     = 2.725d0

    x = frequency * 10**9 ! Input frequency is in GHz
    x = h*x/(k_B*T_0)
    ant2thermo = (exp(x)-1.d0)**2 / (x**2 * exp(x))

  end function ant2thermo

  subroutine sparse_degrade_nest(imaps, ipix, inside, omaps, opix, onside, counts)
    implicit none
    real(dp),     dimension(:,:)          :: imaps
    real(dp),     dimension(:,:), pointer :: omaps
    integer(i4b), dimension(:)            :: ipix
    integer(i4b), dimension(:),   pointer :: opix
    integer(i4b), dimension(:), allocatable :: map2mask
    integer(i4b), dimension(:), pointer, optional :: counts
    integer(i4b)             :: steps, inside, onside, nmaps
    integer(i4b)             :: i, j, k, l, mapn, n, m, inpix, onpix, step

    if(associated(omaps)) deallocate(omaps)
    if(associated(opix))  deallocate(opix)
    if(present(counts)) then; if(associated(counts)) deallocate(counts); end if

    nmaps = size(imaps(1,:))

    n = size(ipix)
    inpix = 12*inside**2
    onpix = 12*onside**2
    allocate(map2mask(0:inpix-1))
    map2mask = -1
    do i = 1, size(ipix)
       map2mask(ipix(i)) = i
    end do

    ! Determine the number of reduced pixels
    step = inpix/onpix
    m = 0; do i = 1, n; if(modulo(ipix(i),step) == 0) m = m+1; end do
    allocate(opix(m))
    allocate(omaps(m,nmaps))
    if(present(counts)) allocate(counts(m))
    omaps = 0
    j = 1
    do i = 1, n
       if(modulo(ipix(i),step) /= 0) cycle
       opix(j) = ipix(i)/step
       do mapn = 1, nmaps
          m = 0
          do k = 1, step
             l = map2mask(ipix(i)+k-1)
             if(l == -1) cycle
             if(imaps(l,mapn) /= hpx_dbadval) then
                m = m+1
                omaps(j,mapn) = omaps(j,mapn) + imaps(l,mapn)
             end if
          end do
          if(present(counts)) counts(j) = m
          if(m > 0) then
             omaps(j,mapn) = omaps(j,mapn) / m
          else
             omaps(j,mapn) = hpx_dbadval
          end if
       end do
       j = j+1
    end do
    deallocate(map2mask)
  end subroutine

  subroutine assert(condition, error_message)
    implicit none
    logical(lgt) :: condition
    character(len=*) :: error_message
    if(condition) return
    write(*,fmt="(a)") error_message
    stop
  end subroutine

  function concat_strs(arr, glue) result(res)
    implicit none
    character(len=*) :: arr(:), glue
    character(len=(len(arr)*size(arr)+len_trim(glue)*(size(arr)-1))) :: res
    integer(i4b)     :: i, j, n, m
    j = 1
    res = ""
    m = len_trim(glue)
    do i = 1, size(arr)
       n = len_trim(arr(i))
       res(j:j+n-1) = trim(arr(i))
       j = j+n
       if(i /= size(arr)) then
          res(j:j+m-1) = trim(glue)
          j = j+m
       end if
    end do
  end function

  function bor(arr) result(res)
    implicit none
    integer(i4b) :: arr(:), res, i
    res = 0
    do i = 1, size(arr); res = ior(res, arr(i)); end do
  end function
  !----------------------------------------------------------------------------------------------
  ! checks for healpix-not-a-number
  ! --------------------------------------------------------------------------------------------

  elemental function healnot(number) result(res)
    implicit none
    real(dp),     intent(in) :: number
    logical(lgt)             :: res
    res = abs(-1.6375d30 - number) < 1d25
  end function
  elemental function healnot_sp(number) result(res)
    implicit none
    real(sp),     intent(in) :: number
    logical(lgt)             :: res
    res = abs(-1.6375d30 - number) < 1d25
  end function
  elemental function healok(number) result(res)
    implicit none
    real(dp),     intent(in) :: number
    logical(lgt)             :: res
    res = .not. healnot(number)
  end function
  elemental function healok_sp(number) result(res)
    implicit none
    real(sp),     intent(in) :: number
    logical(lgt)             :: res
    res = .not. healnot(number)
  end function

  ! Shortest angular distance from b to a
  function ang_diff(a,b) result(d)
    implicit none
    real(dp) :: a, b, d, an, bn, dn
    dn = mod(a-b,2*pi)
    if(abs(dn) < pi) then
       d = dn
    else
       if(dn > 0) then
          d = dn-2*pi
       else
          d = dn+2*pi
       end if
    end if
  end function

  ! Yes, fortran actually needs these three do be different functions
  subroutine set_ordering_dp(new, old, map)
    implicit none
    integer(i4b) :: new, old, nside
    real(dp), dimension(0:,:) :: map
    nside = npix2nside(size(map,1))
    if(new /= old) then
       if(new == RING) then
          call convert_nest2ring(nside, map)
       else
          call convert_ring2nest(nside, map)
       end if
    end if
  end subroutine
  subroutine set_ordering_sp(new, old, map)
    implicit none
    integer(i4b) :: new, old, nside
    real(sp), dimension(0:,:) :: map
    nside = npix2nside(size(map,1))
    if(new /= old) then
       if(new == RING) then
          call convert_nest2ring(nside, map)
       else
          call convert_ring2nest(nside, map)
       end if
    end if
  end subroutine
  subroutine set_ordering_i4b(new, old, map)
    implicit none
    integer(i4b) :: new, old, nside
    integer(i4b), dimension(0:,:) :: map
    nside = npix2nside(size(map,1))
    if(new /= old) then
       if(new == RING) then
          call convert_nest2ring(nside, map)
       else
          call convert_ring2nest(nside, map)
       end if
    end if
  end subroutine
  subroutine set_ordering_pix(nside, new, old, from, to)
    implicit none
    integer(i4b) :: new, old, from, to, nside
    if(new /= old) then
       if(new == RING) then
          call nest2ring(nside, from, to)
       else
          call ring2nest(nside, from, to)
       end if
    else
       to = from
    end if
  end subroutine

  subroutine udgrade(imap, order, omap, onside)
    implicit none
    integer(i4b) :: inside, order, onside
    real(dp)     :: imap(:,:), omap(:,:)
    inside = npix2nside(size(imap,1))
    if(inside /= onside) then
       if(order == RING) then
          call udgrade_ring(imap, inside, omap, onside)
       else
          call udgrade_nest(imap, inside, omap, onside)
       end if
    else
       omap = imap
    end if
  end subroutine

  function asymmetry(mat, loc) result(res)
    implicit none
    real(dp)     :: mat(:,:), res, tmp
    real(dp), dimension(:), allocatable :: a
    integer(i4b), optional :: loc(2)
    integer(i4b) :: i, j, n, loc_(2)
    n   = size(mat,1)
    allocate(a(n))
    do i = 1, n; a(i) = 1/sqrt(abs(mat(i,i))); end do
    res = 0
    do i = 1, n
      do j = i, n
         tmp = abs(mat(i,j) - mat(j,i))*a(i)*a(j)
         if(tmp > res) then
            loc_ = (/i,j/)
            res = max(res, tmp)
         end if
      end do
    end do
    deallocate(a)
    if(present(loc)) loc = loc_
  end function

  function corrmed2mean(medcorr) result(res)
    implicit none
    real(dp) :: medcorr, res, x
    real(dp), parameter, dimension(0:8) :: coeff = [ -0.000977954415310288, &
     & 1.96168939384321, -0.252449424732305, -15.3520128943704, &
     & 101.915897258133, -326.802608791593, 567.461983026053, &
     & -509.678772635584, 185.574602788006 ]
    integer(i4b) :: i, s
    if(medcorr < 0) then
       s = -1; x = -medcorr
    else
       s =  1; x =  medcorr
    end if
    res = 0
    do i = 0, size(coeff)-1
       res = res + coeff(i) * x**i
    end do
    if(res > 1) then; res = 1; elseif(res < 0) then; res = 0; end if
    res = s*res
  end function

  function fmod(a,b) result(c)
    implicit none
    real(dp) :: a, b, c
    c = a - floor(a/b)*b
  end function

  subroutine normalize_angles(array, maxang)
    implicit none
    real(dp) :: array(:), maxang
    integer(i4b) :: i
    do i = 1, size(array)
      array(i) = fmod(array(i), maxang)
    end do
  end subroutine

  subroutine make_angles_safe(array, maxang)
    implicit none
    real(dp) :: array(:), offset, maxang, tol, tmp
    integer(i4b) :: i
    call normalize_angles(array, maxang)
    tol = maxang/4
    offset = 0
    do i = 2, size(array)
       tmp = array(i)+offset-array(i-1)
       if(abs(tmp) > tol) then
          if(tmp > 0) then
             offset = offset - maxang
          else
             offset = offset + maxang
          end if
       end if
       array(i) = array(i) + offset
    end do
  end subroutine

  ! Find the equivalent angle of a closest to b in value
  function nearest_ang(a,b) result(c)
    implicit none
    real(dp) :: a, b, c
    c = modulo(b-a,2*pi)
    if(c > pi) c = c-2*pi
    c = a+c
  end function

  subroutine shift_angle_array(array, shift, maxang)
    implicit none
    real(dp) :: array(:), shift, maxang
    call make_angles_safe(array, maxang)
    call shift_array_interpol(array, shift)
    call normalize_angles(array, maxang)
  end subroutine

  function average_ang(angs) result(res)
    implicit none
    real(dp) :: angs(:), work(size(angs)), res
    work = angs
    call make_angles_safe(work, 2*pi)
    res = mean(work)
  end function

  function average_pointing(point) result(res)
    implicit none
    real(dp)     :: point(:,:), res(2), va(3), v(3)
    integer(i4b) :: i
    va = 0
    do i = 1, size(point,2)
       call ang2vec(point(2,i), point(1,i), v)
       va = va + v
    end do
    call vec2ang(va, res(2), res(1))
  end function

  ! Returns cos(2*psi) and sin(2*psi) of the rotation induced by
  ! parallel transport from vec1 to vec2.
  subroutine qu_transport_rot(vec1, vec2, cos2psi, sin2psi)
    implicit none

    real(dp), dimension(3), intent(in)  :: vec1, vec2
    real(dp),               intent(out) :: cos2psi, sin2psi

    integer(i4b) :: i, j, nfield
    real(dp) :: len_u, len_v, z, cos_theta, sgn, c1, s1, c2, s2
    real(dp), dimension(3) :: u, v

    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       cos2psi = 1; sin2psi = 0
       return
    end if

    sgn    = 1.d0
    if (vec1(1)*vec2(2)-vec1(2)*vec2(1) < 0.d0) sgn = -1.d0

    ! Rotation from vec1 to vec 2
    u         = vec1(3) * vec1 
    u(3)      = u(3) - 1.d0
    v         = vec2 - z * vec1
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    ! Local angle from vertical
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    ! Double it and calculate cos and sin
    c1 = 2*cos_theta**2-1
    s1 = 2*sgn*sqrt(1-cos_theta**2)*cos_theta

    ! Rotation from vec2 to vec 1; sgn is opposite from 1->2
    u          = vec2(3) * vec2
    u(3)       = u(3) - 1.d0
    v          = vec1 - z * vec2
    len_u      = sqrt(sum(u*u))
    len_v      = sqrt(sum(v*v))
    cos_theta  = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    c2 =  2*cos_theta**2-1
    s2 = -2*sgn*sqrt(1-cos_theta**2)*cos_theta

    cos2psi = c1*c2+s1*s2
    sin2psi = c1*s2-s1*c2
  end subroutine

  subroutine convert_lin2mat(x, y)
    implicit none
    real(dp), dimension(6),   intent(in)  :: x
    real(dp), dimension(:,:), intent(out) :: y

    if (size(y(:,1)) == 1) then
       y(1,1) = x(1)
    else if (size(y(:,1)) == 2) then
       y(1,1) = x(4)
       y(1,2) = x(5)
       y(2,1) = x(5)
       y(2,2) = x(6)
    else
       y(1,1) = x(1)
       y(2,1) = x(2)
       y(1,2) = x(2)
       y(3,1) = x(3)
       y(1,3) = x(3)
       y(2,2) = x(4)
       y(2,3) = x(5)
       y(3,2) = x(5)
       y(3,3) = x(6)
    end if
  end subroutine convert_lin2mat

  ! This function implements our power spectrum ordering
  ! convention. It expands a linear power spectrum
  ! (TT), (TT,TE,EE) or (TT,TE,TB,EE,EB,BB) into
  ! a (T,E,B)x(T,E,B) matrix.
  subroutine convert_lin2mat2(x, y)
    implicit none
    real(dp), dimension(:),   intent(in)  :: x
    real(dp), dimension(3,3), intent(out) :: y

    y = 0
    if (size(x) == 1) then
       y(1,1) = x(1)
    else if (size(x) == 3) then
       y(1,1) = x(1)
       y(1,2) = x(2)
       y(2,1) = x(2)
       y(2,2) = x(3)
    else
       y(1,1) = x(1)
       y(2,1) = x(2)
       y(1,2) = x(2)
       y(3,1) = x(3)
       y(1,3) = x(3)
       y(2,2) = x(4)
       y(2,3) = x(5)
       y(3,2) = x(5)
       y(3,3) = x(6)
    end if
  end subroutine convert_lin2mat2

  subroutine convert_mat2lin2(y, x)
    implicit none
    real(dp), dimension(3,3), intent(in)  :: y
    real(dp), dimension(:),   intent(out) :: x

    if (size(x) == 1) then
       x(1) = y(1,1)
    else if (size(x) == 3) then
       x(1) = y(1,1)
       x(2) = y(1,2)
       x(3) = y(2,2)
    else
       x(1) = y(1,1)
       x(2) = y(2,1)
       x(3) = y(3,1)
       x(4) = y(2,2)
       x(5) = y(2,3)
       x(6) = y(3,3)
    end if
  end subroutine

  subroutine convert_ring2nest_sparse(nside, pixels)
    integer(i4b)                :: nside
    integer(i4b), dimension(1:) :: pixels
    integer(i4b)                :: i, ringpix
    do i = 1, size(pixels)
       ringpix=pixels(i)
       call ring2nest(nside, ringpix, pixels(i))
    end do
  end subroutine convert_ring2nest_sparse

  subroutine decimate_timestream_2d(n, trans, d_in, d_out, angle)
    implicit none

    integer(i4b),                 intent(in)             :: n
    logical(lgt),                 intent(in)             :: trans
    real(dp),     dimension(:,:), intent(inout)          :: d_in
    real(dp),     dimension(:,:), intent(out)            :: d_out
    logical(lgt),                 intent(in),   optional :: angle
    
    integer(i4b) :: i, j, p, q
    real(dp), allocatable, dimension(:) :: in, out

    if (trans) then
       p = size(d_in,2)
       q = size(d_out,2)
       allocate(in(p), out(q))
       do i = 1, size(d_in,1)
          in = d_in(i,:)
          call decimate_timestream_1d(n, in, out, angle)          
          d_out(i,:) = out
       end do
       deallocate(in, out)
    else
       do i = 1, size(d_out(1,:))
          call decimate_timestream_1d(n, d_in(:,i), d_out(:,i), angle)
       end do
    end if

  end subroutine decimate_timestream_2d

  subroutine decimate_timestream_1d(n, d_in, d_out, angle)
    implicit none
    
    integer(i4b),               intent(in)             :: n
    real(dp),     dimension(:), intent(inout)          :: d_in
    real(dp),     dimension(:), intent(out)            :: d_out
    logical(lgt),               intent(in),   optional :: angle
    
    integer(i4b) :: i, j

    if (n == 1) then
       d_out = d_in
       return
    end if

    if (present(angle)) then
       if (angle) call make_angles_safe(d_in, 2.d0*pi)
    end if
    do j = 1, size(d_out)
       d_out(j) = sum(d_in((j-1)*n+1:j*n)) / real(n,dp)
    end do

  end subroutine decimate_timestream_1d

  subroutine azbin_tod(tod, az, binned, counts_out)
    implicit none
    real(sp)     :: tod(:), az(:), binned(:)
    integer(i4b), optional :: counts_out(:)
    integer(i4b) :: nbin, i, j, k, l, m, n
    integer(i4b) :: counts(size(binned))
    binned = 0
    counts = 0
    do i = 1, size(tod)
       j = min(int(size(binned)*fmod(az(i)/(2*pi),1d0))+1, size(binned))
       binned(j) = binned(j) + tod(i)
       counts(j) = counts(j) + 1
    end do
    where(counts /= 0) binned = binned / counts
!    where(counts /= 0) binned = (binned - mean(binned)) / counts
    where(counts == 0) binned = NaN
    if (present(counts_out)) counts_out = counts
  end subroutine

  pure function ifeli(cond, a, b) result(r)
    implicit none
    logical(lgt), intent(in) :: cond
    integer(i4b), intent(in) :: a, b
    integer(i4b) :: r
    if(cond) then
       r = a
    else
       r = b
    end if
  end function

  function ifelc(cond, a, b) result(r)
    implicit none
    logical(lgt) :: cond
    character(len=*) :: a, b
    character(len=ifeli(cond,len(a),len(b))) :: r
    if(cond) then
       r = a
    else
       r = b
    end if
  end function

  subroutine dset(level, id)
    implicit none
    integer(i4b), optional :: level, id
    if(present(level)) ddata%level = level
    if(present(id))    ddata%id    = id
  end subroutine

  subroutine set_debug_level(level)
    implicit none
    integer(i4b) :: level
    ddata%level = level
  end subroutine

  function dtest(level)
    implicit none
    logical(lgt) :: dtest
    integer(i4b) :: level
    dtest = ddata%level >= level
  end function

  ! Debug routine
  subroutine dmem(name, level, quiet)
    implicit none
    character(len=*)       :: name
    integer(i4b), optional :: level
    logical(lgt), optional :: quiet
    logical(lgt)           :: q
    integer(i4b)           :: lev
    real(dp),     save     :: t1 = 0
    real(dp)               :: t2
    lev = 1; if(present(level)) lev = level
    q   = .false.; if(present(quiet)) q = quiet
    call wall_time(t2)
    if(t1 == 0) t1 = t2
    if(lev <= ddata%level .and. .not. q) write(*,fmt="(i4,2f11.3,f14.5,a)") ddata%id, &
     & get_mem_use()/1024d0**2, get_max_mem_use()/1024d0**2, t2-t1, " " // trim(name)
    t1 = t2
  end subroutine

  subroutine compute_conf_limit(vals, threshold, type, bounds)
    implicit none

    real(dp),         dimension(1:), intent(in)  :: vals
    real(dp),                        intent(in)  :: threshold
    character(len=*),                intent(in)  :: type
    real(dp),         dimension(2),  intent(out) :: bounds

    integer(i4b) :: i, j, n
    real(dp)     :: frac, median, upper, lower
    real(dp), allocatable, dimension(:) :: tmp

    ! ignore NaN
    n = count(vals == vals)
    allocate(tmp(n))
    j = 0
    do i = 1, size(vals)
       if(vals(i) /= vals(i)) cycle
       j = j+1
       tmp(j) = vals(i)
    end do

    ! And sort them
    call QuickSort_real(tmp)
    median = tmp(max(1,min(n,int(0.50000d0*n))))
    lower  = median - tmp(max(1,min(n,int(0.15865d0*n))))
    upper  = tmp(max(1,min(n,int(0.84134d0*n)))) - median

    if (trim(type) == 'two_sided') then
       bounds(1) = median - threshold*lower
       bounds(2) = median + threshold*upper
    else if (trim(type) == 'lower') then
       bounds(1) = median - threshold*lower
       bounds(2) = tmp(n)
    else if (trim(type) == 'upper') then
       bounds(1) = tmp(1)
       bounds(2) = median + threshold*upper
    else
       write(*,*) 'quiet_utils: Unknown split type!'
       stop
    end if
    deallocate(tmp)
  end subroutine compute_conf_limit

  subroutine abort(msg)
    implicit none
    character(len=*) :: msg
    write(stderr,'(a)') msg
    stop
  end subroutine

  ! Given a mask, produces a list of the matching indices
  subroutine wherei(mask, inds)
    implicit none
    logical(lgt)                            :: mask(:)
    integer(i4b), allocatable, dimension(:) :: inds
    integer(i4b)                            :: i, n
    if(allocated(inds)) deallocate(inds)
    n = count(mask)
    allocate(inds(n))
    n = 0
    do i = 1, size(mask)
       if(.not. mask(i)) cycle
       n = n+1
       inds(n) = i
    end do
  end subroutine

  ! Given a list of integers, produces a list of the unique values
  subroutine uniqi(vals, uval)
    implicit none
    integer(i4b) :: vals(:)
    integer(i4b), allocatable, dimension(:) :: s, uval
    integer(i4b) :: i, j, n
    if(allocated(uval)) deallocate(uval)
    if(size(vals) == 0) then
       allocate(uval(0))
       return
    end if
    allocate(s(size(vals)))
    s = vals
    call quicksort_int(s)
    n = 1
    do i = 2, size(s)
       if(s(i-1) /= s(i)) n = n+1
    end do
    allocate(uval(n))
    n = 1
    uval(1) = s(1)
    do i = 2, size(s)
       if(s(i-1)==s(i)) cycle
       n = n+1
       uval(n) = s(i)
    end do
    deallocate(s)
  end subroutine

  ! Why does healpix make this so hard?
  function pix2vec(nside, order, pix) result(vec)
    implicit none
    integer(i4b) :: nside, order, pix
    real(dp)     :: vec(3)
    if(order == ring) then
       call pix2vec_ring(nside, pix, vec)
    else
       call pix2vec_nest(nside, pix, vec)
    end if
  end function

  subroutine pix2ang(nside, order, pix, theta, phi)
    implicit none
    integer(i4b) :: nside, order, pix
    real(dp)     :: theta, phi
    if(order == ring) then
       call pix2ang_ring(nside, pix, theta, phi)
    else
       call pix2ang_nest(nside, pix, theta, phi)
    end if
  end subroutine

  function vec2pix(nside, order, vec) result(pix)
    implicit none
    integer(i4b) :: nside, order, pix
    real(dp)     :: vec(3)
    if(order == ring) then
       call vec2pix_ring(nside, vec, pix)
    else
       call vec2pix_nest(nside, vec, pix)
    end if
  end function

  function ang2pix(nside, order, theta, phi) result(pix)
    implicit none
    integer(i4b) :: nside, order, pix
    real(dp)     :: theta, phi
    if(order == ring) then
       call ang2pix_ring(nside, theta, phi, pix)
    else
       call ang2pix_nest(nside, theta, phi, pix)
    end if
  end function

  function irange_one(n) result(arr)
    implicit none
    integer(i4b) :: n, arr(n), i
    do i = 1, n; arr(i) = i; end do
  end function

  function irange_two(from, to) result(arr)
    implicit none
    integer(i4b) :: from, to, arr(to-from+1), i
    do i = 1, size(arr); arr(i) = from+i-1; end do
  end function

  function irange_three(from, to, step) result(arr)
    implicit none
    integer(i4b) :: from, to, step, arr((to-from+step)/step), i, j
    j = 0
    do i = from, to, step
       j = j+1
       arr(j) = i
    end do
  end function

  function lesssim(a,b) result(r)
    implicit none
    real(dp), intent(in) :: a, b
    logical(lgt)         :: r
    r = (a < b + maxval(abs([a,b]))*1d-12)
  end function

  function moresim(a,b) result(r)
    implicit none
    real(dp), intent(in) :: a, b
    logical(lgt)         :: r
    r = (a > b - maxval(abs([a,b]))*1d-12)
  end function

  subroutine dump_matrix_mat(mat, fname, unit, fmt, idx)
    implicit none
    real(dp),         intent(in)           :: mat(:,:)
    character(len=*), intent(in), optional :: fname, fmt
    integer(i4b),     intent(in), optional :: unit
    logical(lgt),     intent(in), optional :: idx
    character(len=256)                     :: fmt_
    integer(i4b)                           :: i, j, unit_
    logical(lgt)                           :: doclose, idx_
    fmt_ = '(e12.4)'; if(present(fmt)) fmt_ = fmt
    idx_ = .false.;   if(present(idx)) idx_ = idx
    doclose = .false.
    unit_ = stdout
    if(present(unit)) then
       unit_ = unit
    elseif(present(fname)) then
       unit_ = getlun()
       open(unit_,file=fname)
       doclose = .true.
    end if
    do i = 1, size(mat,1)
       if(idx_) write(unit_,fmt="(i9)",advance="no") i
       do j = 1, size(mat,2)
          write(unit_,fmt=fmt_,advance="no") mat(i,j)
       end do
       write(unit_,*)
    end do
    if(doclose) then
       close(unit_)
    end if
  end subroutine

  subroutine dump_matrix_vec(vec, fname, unit, fmt, idx)
    real(dp),         intent(in)           :: vec(:)
    character(len=*), intent(in), optional :: fname, fmt
    integer(i4b),     intent(in), optional :: unit
    logical(lgt),     intent(in), optional :: idx
    call dump_matrix_mat(reshape(vec,[size(vec),1]), fname, unit, fmt, idx)
  end subroutine

  function nside2order(nside) result(order)
    implicit none
    integer(i4b) :: nside, order
    order = nint(log(real(nside,dp))/log(2d0))
  end function

  function bench_init(bench, name) result(id)
    implicit none
    type(benchmarker) :: bench
    character(len=*)  :: name
    integer(i4b)      :: id
    id              = bench%n+1
    bench%nsamp(:,id) = 0
    bench%dt(:,id)    = 0
    bench%name(id)  = name
    bench%n         = id
  end function

  subroutine bench_start(bench, id)
    implicit none
    type(benchmarker) :: bench
    integer(i4b)      :: id, me, omp_get_thread_num
    me = omp_get_thread_num()+1
    call wall_time(bench%t(1,me,id))
  end subroutine

  subroutine bench_stop(bench, id)
    implicit none
    type(benchmarker) :: bench
    integer(i4b)      :: id, me, omp_get_thread_num
    me = omp_get_thread_num()+1
    call wall_time(bench%t(2,me,id))
    bench%dt(me,id)    = bench%dt(me,id) + bench%t(2,me,id) - bench%t(1,me,id)
    bench%nsamp(me,id) = bench%nsamp(me,id) + 1
  end subroutine

  ! Make bins with exponential steps. Each bin contains at least one element.
  subroutine make_exp_bins(in, bins)
    implicit none
    integer(i4b), intent(in)  :: in
    integer(i4b), intent(out) :: bins(:,:)
    integer(i4b)              :: i, j, n
    real(dp)                  :: x(0:size(bins,2))
    n = size(bins,2)
    call assert(n <= in, "Cannot bin into more bins than data points!")
    ! x contains the ideal ending points for each bin
    x = exp((irange(n+1)-1)*log(real(in,dp))/n)
    bins(1,1) = 1
    do i = 1, n-1
       bins(2,i)   = max(bins(1,i),floor(x(i)))
       bins(1,i+1) = bins(2,i)+1
    end do
    bins(2,n) = in
  end subroutine


  ! ============================================================================
  ! "WMAP_Read_NInv" reads a WMAP N-Inverse FITS file.
  !
  ! If the output array is  unassociated then this routine will allocate space
  ! for it.  This is why it is declared as a pointer.  If a destination array is
  ! supplied then be sure that it is large enough!
  !
  ! This routine requires the FITSIO library, available from
  ! http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html
  !
  ! Arguments:
  !	File      - The name of the file.
  !	Status    - A status code:  0=success.  Most of these are FITS status
  !	            codes, though memory allocation errors are also passed back
  !	            unchanged.
  !	NInv      - The square N-inverse array.
  !	NElements - The number of elements on each side of NInv.  Optional.
  !	NPixels   - The number of pixels in the maps that N-inverse applies to.
  !	            Optional.
  !
  ! Written by Michael R. Greason, SSAI, 13 March 2006.
  ! ============================================================================
  Subroutine WMAP_Read_NInv (File, Status, NInv, NElements, NPixels)
    !
    Implicit None
    !
    Character (*),                 Intent(In)	     :: File
    Integer (Kind=4),              Intent(Out)	     :: Status
    Real (Kind=4), Dimension(:,:), Pointer               :: NInv
    Integer (Kind=4),              Intent(Out), Optional :: NElements
    Integer (Kind=4),              Intent(Out), Optional :: NPixels
    !
    Character (80)                 :: comm
    Integer (Kind=4), Dimension(2) :: naxes
    Integer (Kind=4)               :: unit, lst, rwmode, tmp, npix
    Logical                        :: anyf
    ! ----------------------------------------------------------------------------
    If (Present(NPixels)  ) NPixels   = 0
    If (Present(NElements)) NElements = 0
    Status  = 0
    !
    !			Open the FITS file.  Leave it positioned at the
    !			primary header/data unit (HDU).
    !
    rwmode = 0
    Call FTGIOU (unit, Status)				! Get a free unit #.
    Call FTOPEN (unit, File, rwmode, tmp, Status)		! Open the file.
    If (Status .NE. 0) Return
    !
    !			How big is the array?
    !
    Call FTGISZ (unit, 2, naxes, Status)
    If ((naxes(1) .LE. 0) .OR. (naxes(1) .NE. naxes(2)) .OR. (Status .NE. 0)) GoTo 99

    !
    !			How many pixels are in the base map?  Start by looking
    !			at the NPIX keyword; if that isn't found, use LASTPIX.
    !			If neither is found, give up!
    !
    Call FTGKYJ (unit, 'NPIX', npix, comm, Status)
    If (Status .NE. 0) Then
       !
       Status = 0
       Call FTGKYJ (unit, 'LASTPIX', npix, comm, Status)
       If (Status .NE. 0) Go To 99
       npix = npix + 1
       !
    End If
    !
    !			Extract data from this first extension table.
    !
    If (.NOT. Associated(NInv)) Then
       Allocate(NInv(naxes(1), naxes(2)), Stat=Status)
       If (Status .NE. 0) Go To 99
    End If
    !
    tmp = naxes(1) * naxes(2)
    Call FTGPVE (unit, 0, 1, tmp, 0.0E0, NInv, anyf, Status)
    !
    !			Done!  Set the number of pixels, close the FITS file
    !			and return.
    !
    If (Present(NPixels)  ) NPixels   = npix
    If (Present(NElements)) NElements = naxes(1)
    !
99  Continue
    !
    lst = 0
    Call FTCLOS (unit, lst)					! Close the file.
    Call FTFIOU (unit, lst)					! Release the unit #.
    If ((lst .NE. 0) .AND. (Status .EQ. 0)) Status = lst
    !
    Return
    ! ----------------------------------------------------------------------------
  End Subroutine WMAP_Read_NInv

  function dzeroes(n) result(res)
    implicit none
    integer(i4b) :: n
    real(dp)     :: res(n)
    res = 0
  end function

end module quiet_utils
