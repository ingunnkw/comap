module quiet_module_mod
  use quiet_utils
  use quiet_detector_mod
  implicit none

  interface is_alive
     module procedure is_alive_abs, is_alive_rel
  end interface

  interface get_diode_angle
     module procedure get_diode_angle_abs, get_diode_angle_rel
  end interface

  interface get_diode_freq
     module procedure get_diode_freq_abs, get_diode_freq_rel
  end interface

contains

  ! Public interface. These all work with an internal modlist, but
  ! coudl easily be changed to take an optional modlist argument,
  ! so that one could have several sets of modules.

  subroutine initialize_module_mod(parfile)
    implicit none
    character(len=*)   :: parfile
    call init_detector_mod(parfile)
  end subroutine

  function get_num_diodes() result(res)
    implicit none
    integer(i4b) :: res
    res = quiet_horns(1)%ndi
  end function

  function get_num_modules() result(res)
    implicit none
    integer(i4b) :: res
    res = size(quiet_horns)
  end function

  function is_polarization_module(mod) result(res)
    implicit none
    integer(i4b) :: mod
    logical(lgt) :: res
    res = sum(quiet_diodes(quiet_horns(mod)%diodes)%stokes(2)) >= sum(quiet_diodes(quiet_horns(mod)%diodes)%stokes(1))
  end function

  function is_temperature_module(mod) result(res)
    implicit none
    integer(i4b) :: mod
    logical(lgt) :: res
    res = sum(quiet_diodes(quiet_horns(mod)%diodes)%stokes(2)) < sum(quiet_diodes(quiet_horns(mod)%diodes)%stokes(1))
  end function

  function get_partner(mod) result(res)
    implicit none
    integer(i4b) :: mod, res
    if(quiet_horns(mod)%n < 2) then
       res = -1
    else
       res = quiet_horns(mod)%groups(2)
    end if
  end function

  function module_valid(mod) result(res)
    implicit none
    integer(i4b) :: mod
    logical(lgt) :: res
    if(mod < 0 .or. mod >= size(quiet_horns)) then
       res = .false.
    else
       res = .not. any(quiet_diodes(quiet_horns(mod)%diodes)%ok)
    end if
  end function

  subroutine get_module_pos(mod, theta, phi)
    implicit none
    integer(i4b) :: mod
    real(dp)     :: theta, phi
    theta = quiet_horns(mod)%theta
    phi   = quiet_horns(mod)%phi
  end subroutine

  function get_module_fwhm(mod) result(fwhm)
    implicit none
    integer(i4b) :: mod
    real(dp)     :: fwhm
    fwhm = quiet_horns(mod)%fwhm
  end function

  function get_max_fwhm() result(fwhm)
    implicit none
    integer(i4b) :: mod
    real(dp)     :: fwhm
    fwhm = maxval(quiet_horns%fwhm)
  end function

  function get_diode_angle_rel(mod, di) result(res)
    implicit none
    integer(i4b) :: mod, di
    real(dp)     :: res
    res = quiet_diodes(quiet_horns(mod)%diodes(di))%psi
  end function

  function get_diode_angle_abs(diode) result(res)
    implicit none
    integer(i4b) :: mod, di, diode
    real(dp)     :: res
    res = quiet_diodes(diode)%psi
  end function

  function get_diode_freq_rel(mod, di) result(res)
    implicit none
    integer(i4b) :: mod, di
    real(dp)     :: res
    res = quiet_diodes(quiet_horns(mod)%diodes(di))%freq
  end function

  function get_diode_freq_abs(diode) result(res)
    implicit none
    integer(i4b) :: mod, di, diode
    real(dp)     :: res
    res = quiet_diodes(diode)%freq
  end function

  function get_module_freq(mod) result(res)
    implicit none
    integer(i4b) :: mod
    real(dp)     :: res
    res = mean(quiet_diodes(quiet_horns(mod)%diodes)%freq)
  end function

  function is_alive_rel(mod, di) result(res)
    implicit none
    integer(i4b) :: mod, di
    logical(lgt) :: res
    res = quiet_diodes(quiet_horns(mod)%diodes(di))%ok
  end function

  function is_alive_abs(diode) result(res)
    implicit none
    integer(i4b) :: diode, mod, di
    logical(lgt) :: res
    res = quiet_diodes(diode)%ok
  end function

  function get_focalplane_rad() result(res)
    implicit none
    real(dp) :: res
    res = maxval(quiet_horns%theta)
  end function

  function get_center_module() result(res)
    implicit none
    integer(i4b) :: res
    res = minloc(quiet_horns%theta,1)-1
  end function

  subroutine get_modmab(mod, board, slot)
    implicit none
    integer(i4b) :: mod, board, slot
    board = quiet_horns(mod)%board
    slot  = quiet_horns(mod)%slot
  end subroutine

  subroutine diode_abs2rel(diode, mod, di)
    implicit none
    integer(i4b) :: diode, mod, di
    mod = quiet_diodes(diode)%horn
    di  = quiet_diodes(diode)%sub
  end subroutine

  function diode_rel2abs(mod, di) result(diode)
    implicit none
    integer(i4b) :: diode, mod, di
    diode = quiet_horns(mod)%diodes(di)
  end function

  !! Helper functions

  !subroutine read_modules_from_focalplane(file, list)
  !  implicit none
  !  character(len=*)        :: file
  !  character(len=512)      :: line
  !  type(quiet_module_list) :: list
  !  type(quiet_module)      :: m
  !  integer(i4b)            :: nmod, ndi, n, unit, l, mod, diode
  !  character               :: p
  !  call count_mod_di(file, nmod, ndi)
  !  call init_module_list(list, nmod, ndi)
  !  call init_quiet_module(m, ndi)
  !  unit = getlun()
  !  open(unit,file=file,action="read")
  !  do
  !     read(unit,fmt='(a)',end=1) line
  !     line=trim(adjustl(line))
  !     if (line(1:1) == '#') cycle
  !     read(line,*) m%module_number, m%board, m%slot, m%theta, m%phi, p, m%partner, m%fwhm, m%psi, m%freq, m%alive
  !     m%theta = m%theta * pi/180
  !     m%phi   = m%phi   * pi/180
  !     m%psi   = m%psi   * pi/180
  !     m%fwhm  = m%fwhm  * pi/180/60
  !     m%polarization = p == 'p' .or. p == 'P'
  !     m%invalid = .false.

  !     list%modules(m%module_number) = m
  !     ! Temporary sanity check
  !     call sanity_check(m)
  !  end do
  !  1 close(unit)
  !  call free_quiet_module(m)
  !end subroutine

  !subroutine sanity_check(m)
  !  implicit none
  !  type(quiet_module)      :: m
  !  if(.not. m%polarization .and. all(m%psi == m%psi(0))) then
  !     write(*,*) "All temperature diodes have the same sign! Check focalplane file."
  !     stop
  !  end if
  !end subroutine

  !subroutine count_mod_di(file, nmod, ndi)
  !  implicit none
  !  character(len=*)   :: file
  !  character(len=512) :: line
  !  integer(i4b)       :: nmod, ndi, unit, n, mod, ndi2, nline
  !  logical(lgt)       :: first
  !  unit = getlun()
  !  first = .true.
  !  nline = 0
  !  nmod  = 0
  !  open(unit,file=file,action="read")
  !  do while (.true.)
  !     read(unit,'(a)', end=1) line
  !     line = trim(adjustl(line))
  !     if (line(1:1) == '#') cycle
  !     nline = nline + 1
  !     ndi2 = (num_tokens(line, " ") - 8)/3
  !     if(first) then; ndi = ndi2
  !     elseif(ndi2 /= ndi) then
  !        write(*,*) "Inconsistent number of diodes on line " // trim(itoa(nline)) // " in file '" // trim(file) // "'!"
  !        stop
  !     end if
  !     read(line,*) mod
  !     nmod = max(nmod, mod+1)
  !  end do
  !  1 close(unit)
  !end subroutine

  !subroutine init_quiet_module(m, ndi)
  !  implicit none
  !  type(quiet_module) :: m
  !  integer(i4b)       :: ndi
  !  allocate(m%psi(0:ndi-1), m%alive(0:ndi-1), m%freq(0:ndi-1))
  !  m%invalid = .true.
  !end subroutine

  !subroutine free_quiet_module(m)
  !  implicit none
  !  type(quiet_module) :: m
  !  if(allocated(m%psi))   deallocate(m%psi)
  !  if(allocated(m%freq))  deallocate(m%freq)
  !  if(allocated(m%alive)) deallocate(m%alive)
  !end subroutine

  !subroutine init_module_list(list, nmod, ndi)
  !  implicit none
  !  type(quiet_module_list) :: list
  !  integer(i4b)            :: nmod, ndi, i
  !  allocate(list%modules(0:nmod-1))
  !  do i = 0, nmod-1
  !     call init_quiet_module(list%modules(i), ndi)
  !  end do
  !  list%nmod = nmod
  !  list%ndi  = ndi
  !end subroutine

  !subroutine free_module_list(list)
  !  implicit none
  !  type(quiet_module_list) :: list
  !  integer(i4b)            :: i
  !  if(allocated(list%modules)) then
  !     do i = 0, size(list%modules)-1
  !        call free_quiet_module(list%modules(i))
  !     end do
  !     deallocate(list%modules)
  !  end if
  !end subroutine

end module
