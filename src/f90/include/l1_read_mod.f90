module l1_read_mod
  use healpix_types
  use spline_1D_mod
  use quiet_utils
  use quiet_l1_defs
  use quiet_apex_mod
  implicit none

  ! This is a direct F90 port of the QUIET IDL routine called l1_read:

!!$ ;; ========================================================================= ;;
!!$;; L1_read.pro                                                               ;;
!!$;; ========================================================================= ;;
!!$;; Purpose: Read QUIET Level-1 data                                          ;;
!!$;; ------------------------------------------------------------------------- ;;
!!$;; Syntax: L1_read, filename, [data=data, housekeeping=housekeeping, $       ;;
!!$;;                             pointing=pointing,                    $       ;;
!!$;;                             module=module,                        $       ;;
!!$;;                             time_convert=time_convert,            $       ;;
!!$;;                             adc_convert=adc_convert,              $       ;;
!!$;;                             hk_convert=hk_convert,                $       ;;
!!$;;                             convert=convert]                              ;;
!!$;; ------------------------------------------------------------------------- ;;
!!$;; Input Arguments :                                                         ;;
!!$;;   *filename - path to the Level-1 fits file                               ;;
!!$;; Output Arguments :                                                        ;;
!!$;;   *data - radiometer data read from the Level-1 data file                 ;;
!!$;;   *housekeeping - housekeeping data read from the Level-1 data file       ;;
!!$;;   *pointing - pointing data read from the Level-1 data file               ;;
!!$;; Optional Arguments :                                                      ;;
!!$;;   *module - set this keyword to a list of the module numbers you          ;;
!!$;;             wish to get radiometric/housekeeping data from                ;;
!!$;;             Note that the module numbers specify module physical sites in ;;
!!$;;             the array, not electronics sites.                             ;;
!!$;; Optional Keywords :                                                       ;;
!!$;;   *time_convert - When this keyword is set, all time values are           ;;
!!$;;                   converted to seconds since the start of the             ;;
!!$;;                   file.                                                   ;;
!!$;;   *adc_convert  - When this keyword is set, all radiometer data           ;;
!!$;;                   (demod, average, and quadrature) is converted           ;;
!!$;;                   into ADC voltage units.                                 ;;
!!$;;   *hk_convert   - When this keyword is set, housekeeping voltages are     ;;
!!$;;                   corrected for the gain and offset of the housekeeping   ;;
!!$;;                   ADC. Then, thermometers are converted to physical       ;;
!!$;;                   units, but bias monitoring data is left in units of     ;;
!!$;;                   housekeeping voltage.                                   ;;
!!$;;   *convert      - Setting this keyword is equivalent to setting           ;;
!!$;;                   /time_convert, /adc_convert, and /hk_convert.           ;;
!!$;; ------------------------------------------------------------------------- ;;
!!$;; Errors :                                                                  ;;
!!$;;                                                                           ;;
!!$;;   If you see the following error                                          ;;
!!$;;                                                                           ;;
!!$;;     FITS_READ ERROR: No extension for given extname, extver,              ;;
!!$;;     and/or extlevel found                                                 ;;
!!$;;                                                                           ;;
!!$;;   it means that the Level-1 fits file is missing some HDU that the        ;;
!!$;;   program was looking for. Some Level-1 files don't include housekeeping. ;;
!!$;;   Others don't include pointing information. Try rerunning L1_read        ;;
!!$;;   without requesting housekeeping or pointing information.                ;;
!!$;;                                                                           ;;
!!$;;   Another cause of this error is some old Level-1 files that use a        ;;
!!$;;   different naming convention for the housekeeping HDUs. If you know that ;;
!!$;;   housekeeping data is present, but that error shows up everytime you try ;;
!!$;;   to get the housekeeping, then use the /old keyword.                     ;;
!!$;;                                                                           ;;
!!$;; ------------------------------------------------------------------------- ;;
!!$;; To-do :                                                                   ;;
!!$;;   *Other HDUs to read include RQ_STAT, RQ_BIAS, RQ_DB,                    ;;
!!$;;                               RQ_HKRAW, RQ_MISC, RQ_INDEX, RQ_CLOCK       ;;
!!$;;   *Convert bias monitor values to real units (requires opto-isolator      ;;
!!$;;                                               conversion)                 ;;
!!$;;   *Allow user to specify housekeeping channels to read?                   ;;
!!$;; ========================================================================= ;;
!!$;; Last modified : 10/9/08 by colin                                          ;;
!!$;; ========================================================================= ;;
!!$;-

  integer(i4b), private :: num_det = 19

  type RQ_struct
     real(dp),     allocatable, dimension(:,:) :: demod, quad, avg
     integer(i2b), allocatable, dimension(:)   :: scale, quad_scale
     integer(i4b), allocatable, dimension(:,:) :: demod_raw, quad_raw, avg_raw
  end type RQ_struct

  type data_struct
     real(dp)                                     :: start_mjd
     real(dp),        allocatable, dimension(:)   :: time
     integer(i2b),    allocatable, dimension(:)   :: phase
     type(RQ_struct), allocatable, dimension(:)   :: RQ
  end type data_struct

  type point_struct
     real(dp)                                :: start_mjd
     integer(i2b), allocatable, dimension(:) :: mode
     real(dp),     allocatable, dimension(:) :: time
     real(dp),     allocatable, dimension(:) :: encoder_azimuth,   command_azimuth
     real(dp),     allocatable, dimension(:) :: encoder_elevation, command_elevation
     real(dp),     allocatable, dimension(:) :: encoder_deck,      command_deck
  end type point_struct

  integer(i4b), parameter :: default_scale = 6880

  ! This would look better with an enumeration. Do those exist in fortran?
  integer(i4b), parameter :: point_mode = 1
  integer(i4b), parameter :: point_time = 2
  integer(i4b), parameter :: point_encoder_azimuth = 3
  integer(i4b), parameter :: point_encoder_elevation = 4
  integer(i4b), parameter :: point_encoder_deck = 5
  integer(i4b), parameter :: point_command_azimuth = 6
  integer(i4b), parameter :: point_command_elevation = 7
  integer(i4b), parameter :: point_command_deck = 8
  integer(i4b), parameter :: data_time = 9
  integer(i4b), parameter :: data_phase = 10
  integer(i4b), parameter :: data_demod = 11
  integer(i4b), parameter :: data_quad = 12
  integer(i4b), parameter :: data_avg = 13
  integer(i4b), parameter :: data_scale = 14
  integer(i4b), parameter :: data_quad_scale = 15
  integer(i4b), parameter :: data_demod_raw = 16
  integer(i4b), parameter :: data_quad_raw = 17
  integer(i4b), parameter :: data_avg_raw = 18
  integer(i4b), parameter :: hk = 19
  ! Groups of these, for easy manipulation
  integer(i4b), parameter, dimension(8) :: point_all = [ point_mode, point_time, &
   & point_encoder_azimuth, point_encoder_elevation, point_encoder_deck, &
   & point_command_azimuth, point_command_elevation, point_command_deck ]
  integer(i4b), parameter, dimension(10):: data_all  = [ data_time, data_phase, &
   & data_demod, data_quad, data_avg, data_scale, data_quad_scale, data_demod_raw, &
   & data_quad_raw, data_avg_raw ]
  integer(i4b), parameter, dimension(5) :: point_std = [ point_mode, point_time, &
   & point_encoder_azimuth, point_encoder_elevation, point_encoder_deck ]
  integer(i4b), parameter, dimension(4) :: data_std  = [ data_time, data_phase, data_demod, data_avg ]
  integer(i4b), parameter, dimension(9) :: l1_std = [ point_mode, point_time, &
   & point_encoder_azimuth, point_encoder_elevation, point_encoder_deck, data_time, &
   & data_phase, data_demod, data_avg ]
  integer(i4b), parameter, dimension(1) :: hk_all = [ hk ]

  type l1_selector
     ! This has a fixed 'long enough' size instead of being
     ! allocatable, in order to make it more usable as a temporary
     ! type.
     logical(lgt), dimension(10000)   :: on, ron
  end type l1_selector

  integer(i4b), parameter :: l1ok = 0, l1error = 1, l1norq = 2, l1nohdu = 3, l1file = 4

  interface demod_arr_avg
     module procedure demod_arr_avg_single, demod_arr_avg_multi
  end interface

  interface demod_arr_first
     module procedure demod_arr_first_short, demod_arr_first_dp
  end interface
contains

function makesel(all, only, encoder, point_base, mods) result(sel)
  implicit none
  type(l1_selector) :: sel
  logical(lgt), intent(in), optional :: all, encoder, point_base
  integer(i4b), intent(in), optional, dimension(:) :: only, mods

  integer(i4b) :: i

  sel%on = .false.
  sel%ron = .false.

  if(present(all)) then
     sel%on = all
     sel%ron = all
  end if
  if(present(only))       sel%on(only) = .true.
  if(present(point_base)) sel%on((/ point_time, point_mode /)) = .true.
  if(present(encoder)) &
    & sel%on((/ point_encoder_azimuth, point_encoder_elevation, &
    & point_encoder_deck /)) = encoder
  if(present(mods)) sel%ron(mods) = .true.
end function makesel

function l1sel(a, b, c, d, e, f, g, h, i, mods) result(sel)
  implicit none
  integer(i4b), dimension(:), optional :: a, b, c, d, e, f, g, h, i, mods
  type(l1_selector)                    :: sel
  sel%on  = .false.
  sel%ron = .true.
  if(present(a)) sel%on(a) = .true.
  if(present(b)) sel%on(b) = .true.
  if(present(c)) sel%on(c) = .true.
  if(present(d)) sel%on(d) = .true.
  if(present(e)) sel%on(e) = .true.
  if(present(f)) sel%on(f) = .true.
  if(present(g)) sel%on(g) = .true.
  if(present(h)) sel%on(h) = .true.
  if(present(i)) sel%on(i) = .true.
  if(present(mods)) then
     sel%ron = .false.
     sel%ron(mods) = .true.
  end if
end function l1sel

  subroutine allocate_data_struct(num_RQ, num_samples, data, selector)
    implicit none

    integer(i4b),      intent(in)  :: num_RQ, num_samples
    type(data_struct), intent(out) :: data
    type(l1_selector), intent(in), optional :: selector

    integer(i4b) :: i
    type(l1_selector) :: sel
    if(present(selector)) then
       sel = selector
    else
       sel = makesel(all=.true.)
    end if

    if(sel%on(data_time))  allocate(data%time(num_samples))
    if(sel%on(data_phase)) allocate(data%phase(num_samples))
    allocate(data%RQ(num_RQ))

    do i = 1, num_RQ
       if(.not. sel%ron(i)) cycle
       if(sel%on(data_demod))      allocate(data%RQ(i)%demod(0:3,num_samples))
       if(sel%on(data_quad))       allocate(data%RQ(i)%quad(0:3,num_samples))
       if(sel%on(data_avg))        allocate(data%RQ(i)%avg(0:3,num_samples))
       if(sel%on(data_scale))      allocate(data%RQ(i)%scale(num_samples))
       if(sel%on(data_quad_scale)) allocate(data%RQ(i)%quad_scale(num_samples))
       if(sel%on(data_demod_raw))  allocate(data%RQ(i)%demod_raw(0:3,num_samples))
       if(sel%on(data_quad_raw))   allocate(data%RQ(i)%quad_raw(0:3,num_samples))
       if(sel%on(data_avg_raw))    allocate(data%RQ(i)%avg_raw(0:3,num_samples))
    end do

  end subroutine allocate_data_struct

  subroutine deallocate_data_struct(data)
    implicit none

    type(data_struct), intent(inout) :: data

    integer(i4b) :: i

    data%start_mjd = 0.d0
    if (allocated(data%time)) deallocate(data%time)
    if (allocated(data%phase)) deallocate(data%phase)
    if (allocated(data%RQ)) then
       do i = 1, size(data%RQ)
          if (allocated(data%RQ(i)%demod))      deallocate(data%RQ(i)%demod)
          if (allocated(data%RQ(i)%avg))        deallocate(data%RQ(i)%avg)
          if (allocated(data%RQ(i)%quad))       deallocate(data%RQ(i)%quad)
          if (allocated(data%RQ(i)%scale))      deallocate(data%RQ(i)%scale)
          if (allocated(data%RQ(i)%quad_scale)) deallocate(data%RQ(i)%quad_scale)
          if (allocated(data%RQ(i)%demod_raw))  deallocate(data%RQ(i)%demod_raw)
          if (allocated(data%RQ(i)%avg_raw))    deallocate(data%RQ(i)%avg_raw)
          if (allocated(data%RQ(i)%quad_raw))   deallocate(data%RQ(i)%quad_raw)
       end do
       deallocate(data%RQ)
    end if

  end subroutine deallocate_data_struct

  subroutine allocate_point_struct(num_samples, pointing, selector)
    implicit none

    integer(i4b),       intent(in)  :: num_samples
    type(point_struct), intent(out) :: pointing
    type(l1_selector), intent(in), optional :: selector

    type(l1_selector) :: sel
    if(present(selector)) then
       sel = selector
    else
       sel = makesel(all=.true.)
    end if

    if(sel%on(point_mode))              allocate(pointing%mode(num_samples))
    if(sel%on(point_time))              allocate(pointing%time(num_samples))
    if(sel%on(point_encoder_azimuth))   allocate(pointing%encoder_azimuth(num_samples))
    if(sel%on(point_encoder_elevation)) allocate(pointing%encoder_elevation(num_samples))
    if(sel%on(point_encoder_deck))      allocate(pointing%encoder_deck(num_samples))
    if(sel%on(point_command_azimuth))   allocate(pointing%command_azimuth(num_samples))
    if(sel%on(point_command_elevation)) allocate(pointing%command_elevation(num_samples))
    if(sel%on(point_command_deck))      allocate(pointing%command_deck(num_samples))

  end subroutine allocate_point_struct

  subroutine deallocate_point_struct(pointing) 
    implicit none

    type(point_struct), intent(inout) :: pointing

    pointing%start_mjd = 0.d0
    if (allocated(pointing%mode))              deallocate(pointing%mode)
    if (allocated(pointing%time))              deallocate(pointing%time)
    if (allocated(pointing%encoder_azimuth))   deallocate(pointing%encoder_azimuth)
    if (allocated(pointing%encoder_elevation)) deallocate(pointing%encoder_elevation)
    if (allocated(pointing%encoder_deck))      deallocate(pointing%encoder_deck)
    if (allocated(pointing%command_azimuth))   deallocate(pointing%command_azimuth)
    if (allocated(pointing%command_elevation)) deallocate(pointing%command_elevation)
    if (allocated(pointing%command_deck))      deallocate(pointing%command_deck)

  end subroutine deallocate_point_struct

  subroutine allocate_hk_struct(numsamp_bias, numsamp_cryo, numsamp_encl, numsamp_peri, hk)
    implicit none

    integer(i4b),   intent(in) :: numsamp_bias, numsamp_cryo, numsamp_encl, numsamp_peri
    type(hk_struct)            :: hk

    integer(i4b) :: i

    hk%bias%n   = num_bias_W
    hk%bias%n_t = numsamp_bias
    allocate(hk%bias%value(hk%bias%n_t, hk%bias%n), hk%bias%time(hk%bias%n_t), hk%bias%name(hk%bias%n))
    hk%bias%name = bias_W

    hk%cryo%n   = num_cryo_W
    hk%cryo%n_t = numsamp_cryo
    allocate(hk%cryo%value(hk%cryo%n_t, hk%cryo%n), hk%cryo%time(hk%cryo%n_t), hk%cryo%name(hk%cryo%n))
    hk%cryo%name = cryo_W

    hk%encl%n   = num_encl_W
    hk%encl%n_t = numsamp_encl
    allocate(hk%encl%value(hk%encl%n_t, hk%encl%n), hk%encl%time(hk%encl%n_t), hk%encl%name(hk%encl%n))
    hk%encl%name = encl_W

    hk%peri%n   = num_peri_W
    hk%peri%n_t = numsamp_peri
    allocate(hk%peri%value(hk%peri%n_t, hk%peri%n), hk%peri%time(hk%peri%n_t), hk%peri%name(hk%peri%n))
    hk%peri%name = periph_W

    hk%apex%n   = APEX_NUM_TYPE
    hk%apex%n_t = numsamp_bias
    allocate(hk%apex%value(hk%apex%n_t, hk%apex%n), hk%apex%time(hk%apex%n_t), hk%apex%name(hk%apex%n))
    hk%apex%name = apex_names
  end subroutine allocate_hk_struct

  subroutine deallocate_hk_struct(hk)
    implicit none
    type(hk_struct)          :: hk
    integer(i4b) :: i
    call deallocate_hk_type(hk%bias)
    call deallocate_hk_type(hk%cryo)
    call deallocate_hk_type(hk%encl)
    call deallocate_hk_type(hk%peri)
    call deallocate_hk_type(hk%apex)
  end subroutine deallocate_hk_struct

  subroutine deallocate_hk_type(hk)
    implicit none
    type(hk_type) :: hk
    if(allocated(hk%name))  deallocate(hk%name)
    if(allocated(hk%time))  deallocate(hk%time)
    if(allocated(hk%value)) deallocate(hk%value)
  end subroutine

  subroutine cat_hk_struct(in, out)
    implicit none
    type(hk_struct), intent(in)    :: in(:)
    type(hk_struct), intent(inout) :: out
    call cat_hk_type(in%bias, out%bias)
    call cat_hk_type(in%cryo, out%cryo)
    call cat_hk_type(in%encl, out%encl)
    call cat_hk_type(in%peri, out%peri)
    call cat_hk_type(in%apex, out%apex)
  end subroutine

  subroutine cat_hk_type(in, out)
    implicit none
    type(hk_type), intent(in)    :: in(:)
    type(hk_type), intent(inout) :: out
    integer(i4b)                 :: i, j
    if(.not. allocated(in(1)%time)) return
    call deallocate_hk_type(out)
    out%n_t = sum(in%n_t)
    out%n   = in(1)%n
    allocate(out%name(out%n), out%time(out%n_t), out%value(out%n_t,out%n))
    out%name = in(1)%name
    j = 0
    do i = 1, size(in)
       out%time (1+j:in(i)%n_t+j  ) = in(i)%time
       out%value(1+j:in(i)%n_t+j,:) = in(i)%value
       j = j + in(i)%n_t
    end do
  end subroutine

  ! Auxilliary functions

  ! convert electronics site numbers ('AXX') into physical site numbers 
  ! ('RQ_00YY')
  function map_site_electronic_to_physical(input)
    implicit none

    integer(i4b), intent(in)  :: input
    integer(i4b)              :: map_site_electronic_to_physical

    integer(i4b) :: i, n
    integer(i4b), dimension(num_det) :: electronics_site

    ! the electronics site number for the physical site corresponding to
    ! the index, if that makes sense
    electronics_site = [4,5,6,0,1,2,3,16,17,18,11,12,14,15,9,10,13,7,8]
    
    n = 0
    do i = 0, num_det-1
       if (electronics_site(i+1) == input) then
          map_site_electronic_to_physical = i
          return
       end if
    end do

    ! returns -1 if the input value is out of range
    map_site_electronic_to_physical = -1
    return

  end function map_site_electronic_to_physical


  ! convert physical site numbers ('RQ_00YY') into electronics site
  ! numbers ('AXX')
  function map_site_physical_to_electronic(input)
    implicit none

    integer(i4b), intent(in)  :: input
    integer(i4b)              :: map_site_physical_to_electronic

    integer(i4b) :: i, n
    integer(i4b), dimension(num_det) :: physical_site

    ! the physical site number for the electronics site corresponding to
    ! the index, if that makes sense
    physical_site = [3,4,5,6,0,1,2,17,18,14,15,10,11,16,12,13,7,8,9]
    
    n = 0
    do i = 0, num_det-1
       if (physical_site(i+1) == input) then
          map_site_physical_to_electronic = i
          return
       end if
    end do

    ! returns -1 if the input value is out of range
    map_site_physical_to_electronic = -1
    return

  end function map_site_physical_to_electronic


  ! correct housekeeping voltages,
  ! based on Ross' calibration of the housekeeping board ADC
  function hk_board_conversion(input)
    implicit none

    real(dp), intent(in) :: input
    real(dp)             :: hk_board_conversion

    ! conversion constants for the 19Q housekeeping board, from Ross
    real(dp) :: hk_gain   = -1.0145
    real(dp) :: hk_offset = -0.006

    hk_board_conversion = hk_offset + (hk_gain * input)
    return

  end function hk_board_conversion


  ! calibration curve for Lakeshore DT-670 Si diode thermometer
  function hk_dt670_conversion(input)
    implicit none

    real(dp), intent(in) :: input
    real(dp)             :: hk_dt670_conversion

    logical(lgt), save :: first_call = .true.
    real(dp), dimension(144), save :: t, v, y2

    if (first_call) then

       ! temperature values, in K
       t=[1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, &
            & 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, &
            & 3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, &
            & 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, &
            & 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, &
            & 18.0, 18.5, 19.0, 19.5, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, &
            & 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, &
            & 40.0, 42.0, 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 65.0, &
            & 70.0, 75.0, 77.35, 80.0, 85.0, 90.0, 100.0, 110.0, 120.0, 130.0, &
            & 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, &
            & 240.0, 250.0, 260.0, 270.0, 273.0, 280.0, 290.0, 300.0, 310.0, 320.0, &
            & 330.0, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, &
            & 430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0, 500.0] 
       
       ! voltage values, in V
       v=[1.644290, 1.642990, 1.641570, 1.640030, 1.638370, 1.636600, &
            & 1.634720, 1.632740, 1.630670, 1.628520, 1.626290, 1.624000, &
            & 1.621660, 1.619280, 1.616870, 1.614450, 1.612000, 1.609510, &
            & 1.606970, 1.604380, 1.601730, 1.599020, 1.596260, 1.59344, &
            & 1.59057, 1.58764, 1.58465, 1.57848, 1.57202, 1.56533, 1.55845, &
            & 1.55145, 1.54436, 1.53721, 1.53000, 1.52273, 1.51541, 1.49698, &
            & 1.47868, 1.46086, 1.44374, 1.42747, 1.41207, 1.39751, 1.38373, &
            & 1.37065, 1.35820, 1.34632, 1.33499, 1.32416, 1.31381, 1.30390, &
            & 1.29439, 1.28526, 1.27645, 1.26794, 1.25967, 1.25161, 1.24372, &
            & 1.23596, 1.22830, 1.22070, 1.21311, 1.20548, 1.197748, 1.181548, &
            & 1.162797, 1.140817, 1.125923, 1.119448, 1.115658, 1.112810, 1.110421, &
            & 1.108261, 1.106244, 1.104324, 1.102476, 1.100681, 1.098930, 1.097216, &
            & 1.095534, 1.093878, 1.092244, 1.090627, 1.089024, 1.085842, 1.082669, &
            & 1.079492, 1.076303, 1.073099, 1.069881, 1.066650, 1.063403, 1.060141, &
            & 1.056862, 1.048584, 1.040183, 1.031651, 1.027594, 1.022984, 1.014181, &
            & 1.005244, 0.986974, 0.968209, 0.949000, 0.929390, 0.909416, 0.889114, &
            & 0.868518, 0.847659, 0.826560, 0.805242, 0.783720, 0.762007, 0.740115, &
            & 0.718054, 0.695834, 0.673462, 0.650949, 0.628302, 0.621141, 0.605528, &
            & 0.582637, 0.559639, 0.536542, 0.513361, 0.490106, 0.466760, 0.443371, &
            & 0.419960, 0.396503, 0.373002, 0.349453, 0.325839, 0.302161, 0.278416, &
            & 0.254592, 0.230697, 0.206758, 0.182832, 0.159010, 0.135480, 0.112553, &
            & 0.090681]
       
       call spline(t, v, 1.d30, 1.d30, y2)

       first_call = .false.
    end if

    hk_dt670_conversion = splint(t, v, y2, input)
    
  end function hk_dt670_conversion


  function hk_murata_conversion(input)
    implicit none

    real(dp), intent(in) :: input
    real(dp)             :: hk_murata_conversion

    real(dp),     save :: bias_current
    logical(lgt), save :: first_call = .true.
    real(dp), dimension(34), save :: t, r, y2

    if (first_call) then

       bias_current = 10.d-6

       ! temperature values, in degrees C
       t = [-40.d0, -35.d0, -30.d0, -25.d0, -20.d0, -15.d0, -10.d0, -5.d0, 0.d0, 5.d0, 10.d0, &
            & 15.d0, 20.d0, 25.d0, 30.d0, 35.d0, 40.d0, 45.d0, 50.d0, 55.d0, 60.d0, 65.d0, &
            & 70.d0, 75.d0, 80.d0, 85.d0, 90.d0, 95.d0, 100.d0, 105.d0, 110.d0, 115.d0, 120.d0, &
            & 125.d0]

       ! resistance values, in kOhm
       r = [4397.119d0, 3088.599d0, 2197.225d0, 1581.881d0, 1151.037d0, 846.579d0, 628.988d0, &
            & 471.632d0, 357.012d0, 272.500d0, 209.710d0, 162.651d0, 127.080d0, 100.000d0, 79.222d0, &
            & 63.167d0, 50.677d0, 40.904d0, 33.195d0, 27.091d0, 22.224d0, 18.323d0, 15.184d0, 12.635d0, &
            & 10.566d0, 8.873d0, 7.481d0, 6.337d0, 5.384d0, 4.594d0, 3.934d0, 3.380d0, 2.916d0, 2.522d0]

       call spline(t, r, 1.d30, 1.d30, y2)

       first_call = .false.       

    end if

    hk_murata_conversion = splint(t, r, y2, input/bias_current/1.d3)
    
  end function hk_murata_conversion


  ! ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ! ;; MAIN l1_read PROCEDURE;;
  ! ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  subroutine l1_read(unit, filename, data, housekeeping, pointing, modules, selector, status_code)
    implicit none

    integer(i4b),                       intent(in)              :: unit
    character(len=*),                   intent(in)              :: filename
    integer(i4b),        dimension(1:), intent(in)              :: modules
    type(l1_selector),                  intent(in),    optional :: selector
    type(data_struct),                  intent(inout), optional :: data
    type(hk_struct),                    intent(inout), optional :: housekeeping
    type(point_struct),                 intent(inout), optional :: pointing
    integer(i4b),                       intent(out),   optional :: status_code

    logical(lgt) :: exist, casesens, raw, cooked
    integer(i4b) :: i, j, k, l, n_modules, mod, modn, nsamp
    integer(i4b) :: n_data
    integer(i4b), allocatable, dimension(:) :: mods
    character(len=2)   :: mod_text2
    character(len=4)   :: mod_text
    character(len=128) :: label
    character(len=128), dimension(3) :: names_data

    logical(lgt)                  :: anyf, global_scale_tested, global_scale_found
    integer(i4b)                  :: status, blocksize
    integer(i4b)                  :: naxis, bitpix, readwrite, hdutype
    integer(i4b)                  :: nrow, ncol, colnum, datacode, repeat, width
    integer(i4b), dimension(2)    :: naxes
    character(len=80)             :: comment
    real(sp)                      :: nullval

    integer(i2b), allocatable, dimension(:)   :: buffer_i2b
    character*1,  allocatable, dimension(:)   :: buffer_chr
    integer(i4b), allocatable, dimension(:)   :: buffer, scale_tmp, quad_scale_tmp
    integer(i4b), allocatable, dimension(:,:) :: buf_reform
    real(dp),     allocatable, dimension(:)   :: buffer_dp

    real(dp)                      :: starttime, scale
    character(5),  dimension(24)  :: hk_biasmon
    character(10), dimension(3)   :: hk_gnd
    character(8),  dimension(20)  :: hk_cryo
    character(8),  dimension(12)  :: hk_warm

    type(l1_selector) :: sel
    if(present(selector)) then
       sel = selector
    else
       sel = makesel(all=.true.)
    end if

    readwrite = 0 ! Read only
    casesens = .false.
    if(present(status_code)) status_code = l1ok
    status   = 0
    nsamp    = -1 ! nsamp unknown so far

    if(present(data))         call deallocate_data_struct(data)
    if(present(pointing))     call deallocate_point_struct(pointing)
    if(present(housekeeping)) call deallocate_hk_struct(housekeeping)

!    if (present(modules)) then
    n_modules = size(modules)
    allocate(mods(n_modules))
    mods      = modules

    global_scale_tested = .false.
    global_scale_found  = .false.

    ! Check that the file exists
    inquire(file=trim(filename), exist=exist)
    if (.not. exist) then
       write(stderr,*) 'ERROR: File = ', trim(filename), ' does not exist'
       status = l1file; goto 1
    end if
    
!    write(stderr,*) 'Reading L1 data file =', trim(filename)
    
    ! Open the existing file
    call ftopen(unit,trim(filename),readwrite,blocksize,status)
    if(status /= 0) then
        write(stderr,*) "Failed to open fits file " // trim(filename), status
        status = l1file; goto 1
    end if
!    write(stderr,*) 'status = ', status

    starttime = 0.d0

    ! read radiometer data
    if (present(data)) then
       if(sel%on(data_time)) then
          ! Get 100 Hz time stamps
          ! RQ_TIME is not in W-band. Fall back to TELESCOPE if missing.
          call ftmnhd(unit, 2, 'RQ_TIME', 0, status)
          if (status /= 0) then
             status = 0
             call ftmnhd(unit, 2, 'TELESCOPE', 0, status)
          end if
          if (status /= 0) then
             write(stderr,*) 'Both RQ_TIME and HDU TELESCOPE are missing!'
             status = l1file; goto 1
             stop
          end if
          call ftgnrw(unit, nrow, status)
          call ftgncl(unit, ncol, status)
          call ftgcno(unit,casesens,"time",colnum,status)
          call ftgtcl(unit,colnum, datacode,repeat,width,status)
          allocate(buffer(0:repeat*nrow-1))
          call ftgcvj(unit,colnum,1,1,repeat*nrow,0,buffer,anyf,status)
          allocate(data%time(nrow))
          do i = 0, nrow-1
             ! Milliseconds should be multiple of 10, so forcibly correct.
             data%time(i+1) = buffer(2*i) + 1.d-3/86400 * nint(0.1d0*buffer(2*i+1))*10
             data%start_mjd = data%time(1)
          end do
          deallocate(buffer)
       end if

       ! Read data
       allocate(data%RQ(n_modules))
       do modn = 1, n_modules
          status = 0

          if(.not. sel%ron(modn)) cycle

          ! In W-band, scale is stored globally (presumably it is the
          ! same for all RQ in Q-band too). We therefore first look for
          ! a global version, and fall back on the local one otherwise.
          if (.not. global_scale_tested) then
             call ftmnhd(unit, 2, "RQ_CLOCK", 0, status)
             if (status == 0) then
                call ftgnrw(unit, nrow, status)
                call ftgncl(unit, ncol, status)
                allocate(scale_tmp(nrow))
                call ftgcno(unit,casesens,"demodscale",colnum,status)
                call ftgtcl(unit,colnum, datacode,repeat,width,status)
                allocate(buffer_i2b(0:repeat*nrow-1))
                call ftgcvi(unit,colnum,1,1,repeat*nrow,0,buffer_i2b,anyf,status)
                scale_tmp(1:)      = buffer_i2b(0::repeat)
                deallocate(buffer_i2b)
                allocate(quad_scale_tmp(nrow))
                call ftgcno(unit,casesens,"quadscale",colnum,status)
                call ftgtcl(unit,colnum, datacode,repeat,width,status)
                allocate(buffer_i2b(0:repeat*nrow-1))
                call ftgcvi(unit,colnum,1,1,repeat*nrow,0,buffer_i2b,anyf,status)
                quad_scale_tmp(1:)      = buffer_i2b(0::repeat)
                deallocate(buffer_i2b)
                if(scale_tmp(1) == 0 .or. quad_scale_tmp(1) == 0 .or. status /= 0) then
                   deallocate(scale_tmp, quad_scale_tmp)
                   global_scale_found = .false.
                else
                   global_scale_found = .true.
                end if
             end if
             global_scale_tested = .true.
          end if

          mod = mods(modn)
          call int2string(mod, mod_text)
          label = 'RQ_' // trim(mod_text)
          call ftmnhd(unit, 2, label, 0, status)

          call ftgnrw(unit, nrow, status)
          call ftgncl(unit, ncol, status)

          if(sel%on(data_avg))        allocate(data%RQ(modn)%avg(0:3,nrow))
          if(sel%on(data_quad))       allocate(data%RQ(modn)%quad(0:3,nrow))
          if(sel%on(data_demod))      allocate(data%RQ(modn)%demod(0:3,nrow))
          if(sel%on(data_scale))      allocate(data%RQ(modn)%scale(nrow))
          if(sel%on(data_quad_scale)) allocate(data%RQ(modn)%quad_scale(nrow))
          if(sel%on(data_avg_raw))    allocate(data%RQ(modn)%avg_raw(0:3,nrow))
          if(sel%on(data_quad_raw))   allocate(data%RQ(modn)%quad_raw(0:3,nrow))
          if(sel%on(data_demod_raw))  allocate(data%RQ(modn)%demod_raw(0:3,nrow))

          if (status /= 0) then
             !write(stderr,*) 'HDU ' // trim(label) //  ' is missing'
             if(allocated(data%RQ(modn)%avg))        data%RQ(modn)%avg        = NaN
             if(allocated(data%RQ(modn)%quad))       data%RQ(modn)%quad       = NaN
             if(allocated(data%RQ(modn)%demod))      data%RQ(modn)%demod      = NaN
             if(allocated(data%RQ(modn)%scale))      data%RQ(modn)%scale      = 0
             if(allocated(data%RQ(modn)%quad_scale)) data%RQ(modn)%quad_scale = 0
             if(allocated(data%RQ(modn)%avg_raw))    data%RQ(modn)%avg_raw    = 0
             if(allocated(data%RQ(modn)%quad_raw))   data%RQ(modn)%quad_raw   = 0
             if(allocated(data%RQ(modn)%demod_raw))  data%RQ(modn)%demod_raw  = 0
             if(.not. global_scale_found) then
                if(allocated(scale_tmp))      deallocate(scale_tmp)
                if(allocated(quad_scale_tmp)) deallocate(quad_scale_tmp)
             end if
             cycle
             !status = l1norq; goto 1
          end if

          ! Check for inconsistent lengths
          if(nsamp < 0) then
             nsamp = nrow
          else if(nrow /= nsamp) then
             status = 1000
             write(stderr,'(a,2i8)') "L1read: Inconsistent lengths data:", nrow, nsamp
             goto 1
          end if

          ! First get the scale, as it is needed for converting the others
          if(.not. global_scale_found) then
             allocate(scale_tmp(nrow))
             allocate(quad_scale_tmp(nrow))
             call ftgcno(unit,casesens,"scale",colnum,status)
             call ftgtcl(unit,colnum, datacode,repeat,width,status)
             allocate(buffer_i2b(0:repeat*nrow-1))
             call ftgcvi(unit,colnum,1,1,repeat*nrow,0,buffer_i2b,anyf,status)
             scale_tmp(1:)      = buffer_i2b(0::2)
             quad_scale_tmp(1:) = buffer_i2b(1::2)
             deallocate(buffer_i2b)
             if(scale_tmp(1) == 0 .or. quad_scale_tmp(1) == 0 .or. status /= 0) then
                if(status == 0) status = 1001
                write(stderr,*) "Error: No scale data!"
                goto 1
             end if
          end if

          if(sel%on(data_scale))      data%RQ(modn)%scale      = scale_tmp
          if(sel%on(data_quad_scale)) data%RQ(modn)%quad_scale = quad_scale_tmp

          ! Read data
          scale         = scale_tmp(1) ! Scales are actually all the same
          names_data(1) = "demod"
          names_data(2) = "average"
          names_data(3) = "quad"
          do i = 1, 3
             select case(i)
                case(1); raw = sel%on(data_demod_raw);cooked = sel%on(data_demod)
                case(2); raw = sel%on(data_avg_raw);  cooked = sel%on(data_avg)
                case(3); raw = sel%on(data_quad_raw); cooked = sel%on(data_quad)
             end select
             if(.not. raw .and. .not. cooked) cycle
             call ftgcno(unit,casesens,names_data(i),colnum,status)
             call ftgtcl(unit,colnum, datacode,repeat,width,status)
             allocate(buffer(repeat*nrow))
             allocate(buf_reform(0:repeat-1,nrow))
             call ftgcvj(unit,colnum,1,1,repeat*nrow,0,buffer,anyf,status)
             do k = 0, repeat-1
                do j = 1, nrow
                   buf_reform(k,j) = buffer(repeat*(j-1)+k+1)
                end do
             end do
             deallocate(buffer)

             if(raw) then
                select case(i)
                   case(1); data%RQ(modn)%demod_raw = buf_reform
                   case(2); data%RQ(modn)%avg_raw   = buf_reform
                   case(3); data%RQ(modn)%quad_raw  = buf_reform
                end select
             end if
             if(cooked) then
                select case(i)
                   case(1); data%RQ(modn)%demod = buf_reform / 2d0**16 / scale
                   case(2); data%RQ(modn)%avg   = buf_reform / 2d0**16 / scale - 2.d0
                   case(3); data%RQ(modn)%quad  = buf_reform / 2d0**16 / scale - 2.d0
                end select
             end if
             deallocate(buf_reform)
          end do
          if(.not. global_scale_found) deallocate(scale_tmp, quad_scale_tmp)
       end do

       if(sel%on(data_phase)) then
          ! We also need the phase of the 50 Hz clock
          call ftmnhd(unit, 2, 'RQ_CLOCK', 0, status)
          if (status /= 0) then
             write(stderr,*) 'HDU RQ_CLOCK is missing'
             status = l1nohdu; goto 1
          end if
          call ftgnrw(unit, nrow, status)
          call ftgncl(unit, ncol, status)
          if(nsamp < 0) then
             nsamp = nrow
          else if(nrow /= nsamp) then ! Test for nrow consistency
             status = 1000
             write(stderr,'(a,2i8)') "L1read: Inconsistent lengths clock:", nrow, nsamp
             goto 1
          end if
          call ftgcno(unit,casesens,"FiftyHz",colnum,status)
          call ftgtcl(unit,colnum, datacode,repeat,width,status)
          allocate(buffer_chr(0:repeat*nrow-1))
          call ftgcvb(unit,colnum,1,1,repeat*nrow,0,buffer_chr,anyf,status)
          allocate(data%phase(nrow))
          do i = 1, nrow ! Loop necessary because ichar is stupid
             data%phase(i) = (-1)**ichar(buffer_chr(repeat*(i-1)))
          end do
          deallocate(buffer_chr)
       end if
    end if

    ! Read pointing data
    if (present(pointing)) then

!       write(stderr,*) 'Reading telescope pointing data'

       call ftmnhd(unit, 2, 'TELESCOPE', 0, status)
       if (status /= 0) then
          write(stderr,*) 'HDU TELESCOPE is missing'
          status = l1nohdu; goto 1
       end if
       call ftgnrw(unit, nrow, status)
       call ftgncl(unit, ncol, status)

       if(nsamp < 0) then
          nsamp = nrow
       else if(nrow /= nsamp) then
          status = 1000
          write(stderr,'(a,2i8)') "L1read: Inconsistent lengths pointing:", nrow, nsamp
          goto 1
       end if

       ! Read time stamps
       if(sel%on(point_time)) then
          call ftgcno(unit,casesens,"time",colnum,status)
          call ftgtcl(unit,colnum, datacode,repeat,width,status)
          allocate(buffer(0:repeat*nrow-1))
          call ftgcvj(unit,colnum,1,1,repeat*nrow,0,buffer,anyf,status)
          allocate(pointing%time(nrow))
          do i = 0, nrow-1
             ! Milliseconds should be multiple of 10, so forcibly correct.
             pointing%time(i+1) = buffer(2*i) + 1.d-3/86400 * nint(0.1d0*buffer(2*i+1))*10
             pointing%start_mjd = pointing%time(1)
          end do
          deallocate(buffer)
       end if

       ! Read observation mode
       if(sel%on(point_mode)) then
          call ftgcno(unit,casesens,"mode",colnum,status)
          call ftgtcl(unit,colnum, datacode,repeat,width,status)
          allocate(pointing%mode(nrow))
          call ftgcvi(unit,colnum,1,1,repeat*nrow,0,pointing%mode,anyf,status)
       end if

       ! Read encoder positions
       if(sel%on(point_encoder_azimuth))   allocate(pointing%encoder_azimuth(nrow))
       if(sel%on(point_encoder_elevation)) allocate(pointing%encoder_elevation(nrow))
       if(sel%on(point_encoder_deck))      allocate(pointing%encoder_deck(nrow))

       if(any(sel%on((/point_encoder_azimuth,point_encoder_elevation,point_encoder_deck/)))) then
          call ftgcno(unit,casesens,"encoder",colnum,status)
          call ftgtcl(unit,colnum, datacode,repeat,width,status)
          allocate(buffer_dp(repeat*nrow))
          call ftgcvd(unit,colnum,1,1,repeat*nrow,0,buffer_dp,anyf,status)
          if(sel%on(point_encoder_azimuth))   pointing%encoder_azimuth   = buffer_dp(1::3) * pi / 180.d0       
          if(sel%on(point_encoder_elevation)) pointing%encoder_elevation = buffer_dp(2::3) * pi / 180.d0
          if(sel%on(point_encoder_deck))      pointing%encoder_deck      = buffer_dp(3::3) * pi / 180.d0
          deallocate(buffer_dp)
       end if

       ! Read commanded positions
       if(sel%on(point_command_azimuth))   allocate(pointing%command_azimuth(nrow))
       if(sel%on(point_command_elevation)) allocate(pointing%command_elevation(nrow))
       if(sel%on(point_command_deck))      allocate(pointing%command_deck(nrow))

       if(any(sel%on((/point_command_azimuth,point_command_elevation,point_command_deck/)))) then
          call ftgcno(unit,casesens,"commanded",colnum,status)
          call ftgtcl(unit,colnum, datacode,repeat,width,status)
          allocate(buffer_dp(repeat*nrow))
          call ftgcvd(unit,colnum,1,1,repeat*nrow,0,buffer_dp,anyf,status)
          if(sel%on(point_command_azimuth))   pointing%command_azimuth   = buffer_dp(1::3) * pi / 180.d0
          if(sel%on(point_command_elevation)) pointing%command_elevation = buffer_dp(2::3) * pi / 180.d0
          if(sel%on(point_command_deck))      pointing%command_deck      = buffer_dp(3::3) * pi / 180.d0
          deallocate(buffer_dp)
       end if

    end if


    ! Read house keeping data
    if (present(housekeeping) .and. sel%on(hk)) then
       i = get_hk_numsamp(unit, 'HQ_' // bias_W)
       j = get_hk_numsamp(unit, 'HQ_' // cryo_W)
       k = get_hk_numsamp(unit, 'HQ_' // encl_W)
       l = get_hk_numsamp(unit, ['RQ_PERIPH'])
       if (i > 0 .and. j > 0 .and. k > 0 .and. l > 0) then

          call allocate_hk_struct(i, j, k, l, housekeeping)

          ! Read biases
          call read_hk_time(unit, 'HQ_' // housekeeping%bias%name(1), housekeeping%bias%time)
          call read_hk_values(unit, 'HQ_' // housekeeping%bias%name, housekeeping%bias%value)
          
          ! Read cryo tmperatures
          call read_hk_time(unit, 'HQ_' // housekeeping%cryo%name(1), housekeeping%cryo%time)
          call read_hk_values(unit, 'HQ_' // housekeeping%cryo%name, housekeeping%cryo%value)
          
          ! Read enclosure temperatures
          call read_hk_time(unit, 'HQ_' // housekeeping%encl%name(1), housekeeping%encl%time)
          call read_hk_values(unit, 'HQ_' // housekeeping%encl%name, housekeeping%encl%value)
          
          ! Read peripherals
          call read_hk_time(unit, 'RQ_PERIPH', housekeeping%peri%time)
          call read_peripherals(unit, 'RQ_PERIPH', housekeeping%peri%name, housekeeping%peri%value)

          ! Initialize APEX struct with the same time samples as bias, but don't read values
          call read_hk_time(unit, 'HQ_' // housekeeping%bias%name(1), housekeeping%apex%time)
          housekeeping%apex%value = 0.d0

       end if

    end if

    ! Clean up section. To avoid having to repeat this everywhere, I
    ! will go here with a goto when an error occurs
1   if(allocated(scale_tmp))      deallocate(scale_tmp)
    if(allocated(quad_scale_tmp)) deallocate(quad_scale_tmp)
    call ftclos(unit,status)
    if (allocated(mods)) deallocate(mods)
    if(status /= 0) then
       if(present(pointing)) call deallocate_point_struct(pointing)
       if(present(data))     call deallocate_data_struct(data)
       if(.not. present(status_code)) stop
       status_code = status
    end if
  end subroutine L1_read

  ! quite a bit less elegant now, as it has to care about which parts are allocated
  subroutine double_demodulate(data_in, data_out, pointing_in, pointing_out, offset, nstep)
    implicit none

    type(data_struct),  intent(in),  optional :: data_in
    type(data_struct),  intent(out), optional :: data_out
    type(point_struct), intent(in),  optional :: pointing_in
    type(point_struct), intent(out), optional :: pointing_out
    integer(i4b),       intent(in),  optional :: offset, nstep

    integer(i4b) :: i, j, k, n, n_demod, n_modules, sign, off, step

    step = 2; if(present(nstep)) step = nstep

    if (present(data_in)) then
       n_modules = size(data_in%RQ)

       if(allocated(data_in%time)) then
          data_out%start_mjd = data_in%start_mjd+offset/8640000d0
          call demod_arr_first(data_in%time, data_out%time, offset, step)
       end if

       if(allocated(data_in%phase)) then
          call demod_arr_first(data_in%phase, data_out%phase, offset, step)
       end if

       ! Do the double demodulation
       allocate(data_out%RQ(n_modules))
       do j = 1, n_modules
          if(allocated(data_in%RQ(j)%demod)) then
             n = (size(data_in%RQ(j)%demod(0,:))-offset)/step
             allocate(data_out%RQ(j)%demod(0:3,n))
             if(.not. allocated(data_in%phase)) then
                write(stderr,*) "Phase needed for demod double demodulation!"
                stop
             end if
             do k = 0, 3
                do i = 1, n
                   data_out%RQ(j)%demod(k,i) = mean(data_in%RQ(j)%demod(k, &
                    & 1+offset+step*(i-1):0+offset+step*i)*data_in%phase( &
                    & 1+offset+step*(i-1):0+offset+step*i))
                end do
             end do
          end if
          if(allocated(data_in%RQ(j)%quad)) then
             call demod_arr_avg(data_in%RQ(j)%quad, data_out%RQ(j)%quad, offset, step)
          end if
          if(allocated(data_in%RQ(j)%avg)) then
             call demod_arr_avg(data_in%RQ(j)%avg, data_out%RQ(j)%avg, offset, step)
          end if
          if(allocated(data_in%RQ(j)%scale)) then
             call demod_arr_first(data_in%RQ(j)%scale, data_out%RQ(j)%scale, offset, step)
          end if
          if(allocated(data_in%RQ(j)%quad_scale)) then
             call demod_arr_first(data_in%RQ(j)%quad_scale, data_out%RQ(j)%quad_scale, offset, step)
          end if
       end do
    end if

    if (present(pointing_in)) then
       if(allocated(pointing_in%mode)) then
          call demod_arr_first(pointing_in%mode, pointing_out%mode, offset, step)
       end if
       if(allocated(pointing_in%time)) then
          pointing_out%start_mjd = pointing_in%start_mjd+offset/8640000d0
          call demod_arr_first(pointing_in%time, pointing_out%time, offset, step)
       end if
       if(allocated(pointing_in%encoder_azimuth)) then
          call demod_arr_avg(pointing_in%encoder_azimuth, pointing_out%encoder_azimuth, offset, step)
       end if
       if(allocated(pointing_in%encoder_elevation)) then
          call demod_arr_avg(pointing_in%encoder_elevation, pointing_out%encoder_elevation, offset, step)
       end if
       if(allocated(pointing_in%encoder_deck)) then
          call demod_arr_avg(pointing_in%encoder_deck, pointing_out%encoder_deck, offset, step)
       end if
       if(allocated(pointing_in%command_azimuth)) then
          call demod_arr_avg(pointing_in%command_azimuth, pointing_out%command_azimuth, offset, step)
       end if
       if(allocated(pointing_in%command_elevation)) then
          call demod_arr_avg(pointing_in%command_elevation, pointing_out%command_elevation, offset, step)
       end if
       if(allocated(pointing_in%command_deck)) then
          call demod_arr_avg(pointing_in%command_deck, pointing_out%command_deck, offset, step)
       end if
    end if
  end subroutine double_demodulate

  subroutine demod_arr_avg_single(input, output, offset, navg)
    implicit none
    real(dp), dimension(:), intent(in) :: input
    real(dp), dimension(:), allocatable, intent(inout) :: output
    integer(i4b) :: i, n, offset, navg
    n = (size(input)-offset)/navg
    allocate(output(n))
    do i = 1, n
       output(i) = mean(input(1+offset+navg*(i-1):offset+navg*i))
    end do
  end subroutine

  subroutine demod_arr_avg_multi(input, output, offset, navg)
    implicit none
    real(dp), dimension(0:,:), intent(in) :: input
    real(dp), dimension(:,:),  allocatable, intent(inout) :: output
    integer(i4b) :: i, j, n, m, offset, navg
    m = size(input,1)
    n = (size(input,2)-offset)/navg
    allocate(output(0:m-1,n))
    do j = 0, m-1
       do i = 1, n
          output(j,i) = mean(input(j,1+offset+navg*(i-1):offset+navg*i))
       end do
    end do
  end subroutine

  subroutine demod_arr_first_short(input, output, offset, nstep)
    implicit none
    integer(i2b), dimension(:),              intent(in)    :: input
    integer(i2b), dimension(:), allocatable, intent(inout) :: output
    integer(i4b) :: i, n, offset, nstep
    n = (size(input)-offset)/nstep
    allocate(output(n))
    output = input(1+offset:nstep*n+offset:nstep)
  end subroutine

  subroutine demod_arr_first_dp(input, output, offset, nstep)
    implicit none
    real(dp), dimension(:),              intent(in)    :: input
    real(dp), dimension(:), allocatable, intent(inout) :: output
    integer(i4b) :: i, n, offset, nstep
    n = (size(input)-offset)/nstep
    allocate(output(n))
    output = input(1+offset:nstep*n+offset:nstep)
  end subroutine

  subroutine L1_info(filename, module_list, status)
    implicit none
    character(len=*),                        intent(in)  :: filename
    integer(i4b), dimension(:), allocatable, intent(out) :: module_list
    integer(i4b),                            intent(out) :: status

    integer(i4b) :: i,j,k,l,m,n, nhdu, unit, blocksize, readonly, what
    character(len=256) :: label
    character(len=4)   :: mod_text
    character(len=80)  :: key, record, comment
    integer(i4b), dimension(:), allocatable :: table
    logical(lgt) :: ok
    unit = getlun()
    call ftopen(unit,trim(filename),0,blocksize,status)
    if(status /= 0) then
       write(stderr,*) "Could not open fits file '" // trim(filename) // "'!"
       write(stderr,*) "Status is", status
       status = l1file; return
    end if
    ! Set up the module list by scanning through all the hdus
    call ftthdu(unit, nhdu, status)
    allocate(table(0:nhdu-1))
    table = -1
    do i = 0, nhdu-1
       status = 0
       call ftmahd(unit, i+1, what, status)
       call ftgkey(unit,"EXTNAME", record, comment, status)
       key = record(2:len_trim(record)-1)
       if(key(1:3) == "RQ_") then
          k = atoi(key(4:), ok)
          if(ok) table(i) = k
       end if
    end do
    allocate(module_list(count(table >= 0)))
    j = 1
    do i = 0, nhdu-1
       if(table(i) >= 0) then
          module_list(j) = table(i)
          j = j+1
       end if
    end do
    deallocate(table)
    call ftclos(unit, status)
  end subroutine L1_info

  subroutine L1_modlist_full(filename, module_list, status)
    implicit none
    character(len=*),                        intent(in)  :: filename
    integer(i4b), dimension(:), allocatable, intent(out) :: module_list
    integer(i4b),                            intent(out) :: status

    integer(i4b) :: i, nmod

    call L1_info(filename, module_list, status)
    nmod = maxval(module_list)
    deallocate(module_list)
    allocate(module_list(nmod+1))
    do i = 0, nmod
       module_list(i+1) = i
    end do

  end subroutine L1_modlist_full

  function get_hk_numsamp(unit, labels) result(n)
    implicit none

    integer(i4b),                   intent(in) :: unit
    character(len=*), dimension(:), intent(in) :: labels

    logical(lgt)                  :: exists, intel_bug
    integer(i4b)                  :: status, i, hdutype
    integer(i4b)                  :: nrow, n, nfound
    character(len=80)             :: name, comment

    n      = 0
    nfound = 0

    status = 0       
    call ftmahd(unit, 2, hdutype, status)
    do while (status == 0 .and. nfound < size(labels))
       call ftmrhd(unit, 1, hdutype, status)
       call ftgkys(unit, 'EXTNAME', name, comment, status)
       do i = 1, size(labels)
          intel_bug = trim(labels(i)) == trim(name)
          if (intel_bug) then
             call ftgnrw(unit, nrow, status)
             if(n <= 0 .or. nrow < n) n = nrow
             nfound = nfound+1
             exit
          end if
       end do
    end do
  end function get_hk_numsamp

  subroutine read_hk_time(unit, label, time)
    implicit none

    integer(i4b),                   intent(in)  :: unit
    character(len=*),               intent(in)  :: label
    real(dp),         dimension(:), intent(out) :: time

    logical(lgt) :: exist, casesens, anyf
    integer(i4b)                  :: blocksize
    integer(i4b)                  :: nrow, ncol, colnum, repeat, width, n, status, datacode, i
    integer(i4b), dimension(2)    :: naxes
    real(sp)                      :: nullval

    integer(i2b), allocatable, dimension(:)   :: buffer_i2b
    integer(i4b), allocatable, dimension(:)   :: buffer

    casesens = .false.
    status = 0
    n      = size(time)
    call ftmnhd(unit, 2, label, 0, status)
    call ftgnrw(unit, nrow, status)
    call ftgncl(unit, ncol, status)
    if (status /= 0 .or. nrow < n) then
       write(*,*) 'l1_read_mod: Too few samples in housekeeping data for time column'
       write(*,*) '             status = ', status
       write(*,*) '             n      = ', n
       write(*,*) '             nrow   = ', nrow
       write(*,*) '             label  = ', trim(label)
       time = 0
       return
    end if

    ! Read time stamps
    call ftgcno(unit,casesens,"time",colnum,status)
    call ftgtcl(unit,colnum, datacode,repeat,width,status)
    allocate(buffer(0:repeat*n-1))
    call ftgcvj(unit,colnum,1,1,repeat*n,0,buffer,anyf,status)
    do i = 0, n-1
       time(i+1) = buffer(2*i) + buffer(2*i+1) / (24.d0*3600.d0*1000.d0)
    end do
    deallocate(buffer)

  end subroutine read_hk_time

  subroutine read_hk_values(unit, labels, values)
    implicit none

    integer(i4b),                     intent(in)  :: unit
    character(len=*), dimension(:),   intent(in)  :: labels
    real(dp),         dimension(:,:), intent(out) :: values

    logical(lgt) :: exist, casesens, anyf, intel_bug
    integer(i4b)                  :: status, blocksize, i, j
    integer(i4b)                  :: naxis, bitpix, readwrite, hdutype
    integer(i4b)                  :: nrow, ncol, colnum, datacode, repeat, width, n
    integer(i4b), dimension(2)    :: naxes
    real(sp)                      :: nullval

    real(dp), allocatable, dimension(:)   :: buffer
    character(len=80)             :: name, comment

    nrow   = size(values,1)
    status = 0
    casesens = .false.
    call ftmahd(unit, 1, hdutype, status)
    do while (status == 0) 
       call ftmrhd(unit, 1, hdutype, status)
       call ftgkys(unit, 'EXTNAME', name, comment, status)
       do j = 1, size(labels)
          intel_bug = trim(labels(j)) == trim(name)
          if (intel_bug) then
             ! Read vales
             call ftgcno(unit,casesens,"value",colnum,status)
             call ftgtcl(unit,colnum, datacode,repeat,width,status)
             allocate(buffer(0:repeat*nrow-1))
             call ftgcvd(unit,colnum,1,1,repeat*nrow,0,buffer,anyf,status)
             do i = 0, nrow-1
                values(i+1,j) = buffer(i)
             end do
             deallocate(buffer)
             exit
          end if
       end do
    end do
  end subroutine read_hk_values

  subroutine read_peripherals(unit, ext, fields, values)
    implicit none

    integer(i4b),                     intent(in)  :: unit
    character(len=*),                 intent(in)  :: ext
    character(len=*), dimension(:),   intent(in)  :: fields
    real(dp),         dimension(:,:), intent(out) :: values

    logical(lgt) :: exist, casesens, anyf, intel_bug
    integer(i4b)                  :: status, blocksize, i, j
    integer(i4b)                  :: naxis, bitpix, readwrite, hdutype
    integer(i4b)                  :: nrow, ncol, colnum, datacode, repeat, width, n
    integer(i4b), dimension(2)    :: naxes
    real(sp)                      :: nullval

    integer(i4b), allocatable, dimension(:)   :: buffer
    character(len=80)             :: name, comment

    nrow   = size(values,1)
    ncol   = size(values,2)

    status = 0       
    casesens = .false.
    call ftmnhd(unit, 2, ext, 0, status)

    ! Read vales
    do j = 1, ncol
       call ftgcno(unit,casesens,trim(fields(j)),colnum,status)
       if (status /= 0) then
          write(*,*) 'Error -- could not read field = ', trim(fields(j))
       end if
       call ftgtcl(unit,colnum, datacode,repeat,width,status)
       allocate(buffer(0:repeat*nrow-1))
       call ftgcvj(unit,colnum,1,1,repeat*nrow,0,buffer,anyf,status)
       intel_bug = trim(fields(j)) == 'Peripheral_Time' .or. trim(fields(j)) == 'WPID_Heat_Power_UTC' .or. &
            & trim(fields(j)) == 'WPID_Fan_Power_UTC' .or. trim(fields(j)) == 'PSU_UTC'
       if (intel_bug) then
          ! Time field
          do i = 0, nrow-1
             values(i+1,j) = buffer(2*i) + buffer(2*i+1) / (24.d0*3600.d0*1000.d0)
          end do
       else 
          ! Data field; if repeat > 1, we keep only the first column
          do i = 0, nrow-1
             values(i+1,j) = buffer(repeat*i) 
          end do
       end if
       deallocate(buffer)
    end do

  end subroutine read_peripherals


  function hk_binary_to_raw(bits)
    integer(i2b) :: bits
    real(dp)     :: hk_binary_to_raw
    hk_binary_to_raw = -1.0145d0 * (bits/65536.d0 * 6.6d0 - 3.3d0) - 0.006d0
  end function hk_binary_to_raw

end module l1_read_mod
