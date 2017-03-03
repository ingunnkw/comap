module quiet_gain_mod
  use healpix_types
  use quiet_utils
  use quiet_detector_mod
  use math_tools
  implicit none

  type gain_type
     character(len=1)      :: type
     ! For harmonic
     real(dp)              :: mjd0, mean, slope, amp, freq, phase
     ! For piecewise
     real(dp), allocatable :: pieces(:,:)
     ! For old Q-band format
     real(dp)              :: factor, alpha, offset
  end type

  type(gain_type), dimension(:),   allocatable, private :: gains
  real(dp),        dimension(:,:), allocatable, private :: temps ! Needed for old model

  ! Gain file format: We have two different gain behaviors:
  !  1. Piecewise linear
  !  2. Harmolinear
  !
  ! We will use this format:
  !  mod di type params
  ! Which for these can be:
  !  mod di L mjd  val [mjd val [mjd val [...]]]
  !  mod di H mjd0 mean slope amp freq phase

contains

  subroutine initialize_gain_mod(parfile)
    implicit none
    character(len=*)   :: parfile
    character(len=512) :: gainfile, tempfile
    integer(i4b)       :: i
    logical(lgt), save :: initialized = .false.
    if(initialized) return
    initialized = .true.
    call init_detector_mod(parfile)
    call get_parameter(0, parfile, "GAIN_FILE", par_string=gainfile)
    call read_gains(gainfile, gains)
  !print *, "BOut of here is "
    if(any(gains%type == 'Q')) then
       call get_parameter(0, parfile, "GAIN_P3T9_DATAFILE", par_string=tempfile)
       call read_temps(tempfile, temps)
    end if
  end subroutine

  ! Read gains from a file with format mjd, mod, di, gain, dev
  ! Each diode must have the same number of times, and the file
  ! is assumed to be sorted by module, diode, time, in that order.
  subroutine read_gains(gainfile, gains)
    implicit none
    character(len=*)     :: gainfile
    character(len=10000) :: line
    character(len=1)     :: t
    integer(i4b)         :: i, n, unit, m, d, ndi, di
    type(gain_type), allocatable, dimension(:) :: gains
    call free_gains(gains)
    allocate(gains(size(quiet_diodes)))
    gains%type = '0'
    unit = getlun()
    open(unit,file=gainfile,action="read",status="old")
    do
       read(unit,'(a)',end=1) line
       if(line(1:1)=="#") cycle
       n = num_tokens(line, "	 ")
       read(line,*) m, d, t
       di = quiet_horns(m)%diodes(d)
       gains(di)%type = t
       select case(t)
          case('L')
             allocate(gains(di)%pieces(2,(n-3)/2))
             read(line,*) m, d, t, gains(di)%pieces
          case('H')
             read(line,*) m, d, t, gains(di)%mjd0, gains(di)%mean, gains(di)%slope, &
              & gains(di)%amp, gains(di)%freq, gains(di)%phase
          case('Q') ! Old Q-band format
             read(line,*) m, d, t, gains(di)%factor, gains(di)%alpha, gains(di)%offset
          case default; call assert(.false., "Unrecognized gain type: '"//t//"'")
       end select
    end do
1   close(unit)
  end subroutine

  subroutine read_temps(tfile, temps)
    implicit none
    character(len=*),                              intent(in)    :: tfile
    real(dp),         dimension(:,:), allocatable, intent(inout) :: temps
    real(dp)     :: mjd, T
    integer(i4b) :: i, n, unit
    if(allocated(temps)) deallocate(temps)
    unit = getlun()
    open(unit,file=tfile,status="old",action="read")
    n = 0
    do
       read(unit,*,end=1) mjd, T
       n = n+1
    end do
    1 rewind(unit)
    allocate(temps(n,2))
    do i = 1, n
       read(unit,*) temps(i,:)
    end do
    close(unit)
  end subroutine

  subroutine get_gains(mjd, di, res)
    implicit none
    real(dp)     :: mjd(:), res(:)
    integer(i4b) :: di
    res = 0
    if(di < 1 .or. di > size(quiet_diodes)) return
    select case(gains(di)%type)
       case('L'); call lin_interpol(gains(di)%pieces(1,:), gains(di)%pieces(2,:), mjd, res)
       case('H'); res = gains(di)%mean + gains(di)%slope*(mjd-gains(di)%mjd0) + &
         & gains(di)%amp*sin(gains(di)%freq*(mjd-gains(di)%mjd0) + gains(di)%phase)
       case('Q')
          call lin_interpol(temps(:,1), temps(:,2), mjd, res)
          res = gains(di)%factor * (1 + gains(di)%alpha * (res - gains(di)%offset))
       case('0'); return
       !case('0'); call assert(.false., "Uninitialized gain: " // trim(itoa(di)))
       case default; call assert(.false., "Unrecognized gain type: '"//gains(di)%type//"'")
    end select
  end subroutine

  subroutine free_gains(gains)
    implicit none
    type(gain_type), dimension(:), allocatable :: gains
    integer(i4b)                               :: i
    if(.not. allocated(gains)) return
    do i = 1, size(gains)
       call free_gain(gains(i))
    end do
    deallocate(gains)
  end subroutine

  subroutine free_gain(gain)
    implicit none
    type(gain_type) :: gain
    if(allocated(gain%pieces)) deallocate(gain%pieces)
  end subroutine

  function get_gain(mjd, di) result(res)
    implicit none
    real(dp)     :: mjd, res, foo(1)
    integer(i4b) :: di
    call get_gains([mjd], di, foo)
    res = foo(1)
  end function

end module quiet_gain_mod
