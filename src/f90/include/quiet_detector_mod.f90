! This module replaces quiet_module_mod. It is more flexible, and should
! be able to represent both QUIET and planck detectors, for example. To do
! this, it separates the concept of horns from diodes. The horns describe
! the incoming signal on the sky, while the diodes describe components
! that are read off.
!
! Horn format:
! id board slot theta phi fwhm [partner amp [partner amp [...]]]
!
! Diode format
! id horn sub psi freq T P V ok
!
! The horn and diode arrays are made public instead of having a horde
! of accessor functions, and this makes array operations simpler.
! For example, to get the board corresponding to each diode, perhaps for
! jackknife purposes, you would simply do quiet_horns(quiet_diodes%horn)%board.
! Diode indices are absolute, but the relative diode number is available as "sub".
!
! To translate from (mod,sub) to di, use: quiet_horns(mod)%diodes(sub)

module quiet_detector_mod
  use quiet_utils
  implicit none

  type quiet_horn
     integer(i4b) :: id, board, slot, n, ndi
     real(dp)     :: phi, theta, fwhm
     integer(i4b), allocatable :: groups(:), diodes(:)
     real(dp),     allocatable :: amps(:)
  end type

  type quiet_diode
     integer(i4b) :: id, horn, sub
     logical(lgt) :: ok
     real(dp)     :: psi, freq, stokes(3) ! stokes: T, P, V
  end type

  type(quiet_horn),  dimension(:), allocatable, public :: quiet_horns  ! (0:nhorn-1)
  type(quiet_diode), dimension(:), allocatable, public :: quiet_diodes ! (1:ndi)

contains

  subroutine init_detector_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    character(len=512)           :: hfile, dfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    call get_parameter(0, parfile, "HORN_FILE",  par_string=hfile)
    call get_parameter(0, parfile, "DIODE_FILE", par_string=dfile)
    call read_horns (hfile, quiet_horns)
    call read_diodes(dfile, quiet_diodes)
    call setup_diodes(quiet_horns, quiet_diodes)
    initialized = .true.
  end subroutine

  ! Helper functions below

  subroutine read_horns(fname, horns)
    implicit none
    character(len=*),                            intent(in)    :: fname
    type(quiet_horn), dimension(:), allocatable, intent(inout) :: horns
    character(len=64),dimension(:), allocatable                :: tokens
    character(len=512) :: line
    integer(i4b)       :: i, j, k, m, n, id, unit, off
    call free_horns(horns)
    unit = getlun()
    open(unit,file=fname,status="old",action="read")
    ! Get the maximum horn number
    n = -1
    do
       read(unit,'(a)',end=1) line
       if(line(1:1) == "#" .or. line == "") cycle
       read(line,*) k
       n = max(n,k)
    end do
    1 rewind(unit)
    allocate(horns(0:n))
    do
       read(unit,'(a)',end=2) line
       if(line(1:1) == "#" .or. line == "") cycle
       n = num_tokens(line, "	 ")
       allocate(tokens(n))
       call get_tokens(line, "	 ", tokens)
       read(line,*) id
       call free_horn(horns(id))
       off = 6; horns(id)%n = (n-off)/2+1
       allocate(horns(id)%groups(horns(id)%n), horns(id)%amps(horns(id)%n))
       read(line,*) horns(id)%id, horns(id)%board, horns(id)%slot, &
        & horns(id)%theta, horns(id)%phi, horns(id)%fwhm
       horns(id)%groups(1) = id
       horns(id)%amps(1)   = 1d0
       do i = 2, horns(id)%n
          read(tokens((i-2)*2+off+1),*) horns(id)%groups(i)
          read(tokens((i-2)*2+off+2),*) horns(id)%amps(i)
       end do
       deallocate(tokens)
       horns(id)%theta = horns(id)%theta * DEG2RAD
       horns(id)%phi   = horns(id)%phi   * DEG2RAD
       horns(id)%fwhm  = horns(id)%fwhm  * DEG2RAD / 60
    end do
    2 close(unit)
    ! Make sure that missing horns are in a well-defined state
    do i = 0, size(horns)-1
       call assert(allocated(horns(i)%groups), "Error! Incomplete horn definition file!")
    end do
  end subroutine

  subroutine read_diodes(fname, diodes)
    implicit none
    character(len=*),                               intent(in) :: fname
    type(quiet_diode), dimension(:), allocatable, intent(inout) :: diodes
    character(len=512) :: line
    integer(i4b)       :: i, j, k, m, n, id, unit
    real(dp)           :: T, P, V
    call free_diodes(diodes)
    unit = getlun()
    open(unit,file=fname,status="old",action="read")
    ! Get the maximum diode number
    n = 0
    do
       read(unit,'(a)',end=1) line
       if(line(1:1) == "#" .or. line == "") cycle
       read(line,*) k
       n = max(n,k)
    end do
    1 rewind(unit)
    allocate(diodes(n))
    diodes%id = 0
    do
       read(unit,'(a)',end=2) line
       if(line(1:1) == "#" .or. line == "") cycle
       read(line,*) id
       read(line,*) diodes(id)%id, diodes(id)%horn, diodes(id)%sub, &
        & diodes(id)%psi, diodes(id)%freq, T, P, V, diodes(id)%ok
       diodes(id)%stokes = [ T, P, V ]
       diodes(id)%psi = diodes(id)%psi * DEG2RAD
    end do
    2 close(unit)
    do i = 1, size(diodes)
       call assert(diodes(i)%id > 0, "Error! Incomplete diode definition file!")
    end do
  end subroutine

  subroutine setup_diodes(horns, diodes)
    implicit none
    type(quiet_horn), intent(inout) :: horns(0:)
    type(quiet_diode),intent(in)    :: diodes(:)
    integer(i4b)                    :: i, j, k, n
    do i = 0, size(horns)-1
       if(allocated(horns(i)%diodes)) deallocate(horns(i)%diodes)
       n = count(diodes%horn == i)
       allocate(horns(i)%diodes(0:n-1))
       horns(i)%ndi    = n
       horns(i)%diodes = 0
       do j = 1, size(diodes)
          if(diodes(j)%horn /= i) cycle
          k = diodes(j)%sub
          call assert(k >= 0 .and. k < n, &
           & "Diode sub out of range for diode " // trim(itoa(j)))
          call assert(horns(i)%diodes(k) == 0, &
           & "Duplicate sub for diode " // trim(itoa(j)))
          horns(i)%diodes(k) = j
       end do
    end do
  end subroutine

  subroutine free_horn(horn)
    implicit none
    type(quiet_horn) :: horn
    if(allocated(horn%groups)) deallocate(horn%groups)
    if(allocated(horn%amps))   deallocate(horn%amps)
    if(allocated(horn%diodes)) deallocate(horn%diodes)
  end subroutine

  subroutine free_horns(horns)
    implicit none
    type(quiet_horn), dimension(:), allocatable :: horns
    integer(i4b)                                :: i
    if(.not. allocated(horns)) return
    do i = 1, size(horns)
       call free_horn(horns(i))
    end do
    deallocate(horns)
  end subroutine

  subroutine free_diodes(diodes)
    implicit none
    type(quiet_diode), dimension(:), allocatable :: diodes
    if(allocated(diodes)) deallocate(diodes)
  end subroutine

end module
