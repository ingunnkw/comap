module tod2map_utils
  use healpix_types
  use pix_tools
  use math_tools
  use quiet_assembly_mod
  use quiet_pointing_mod
  use quiet_mpi_mod
  use quiet_patch_mod
  use quiet_system_mod
  use quiet_status_mod
  use quiet_fileutils
  use quiet_lx_mod
  use quiet_patch_detect_mod
  implicit none

  integer(i4b), parameter :: rhs_type       = sp, mpi_rhs_type   = mpi_real
  integer(i4b), parameter :: tod_type       = sp, mpi_tod_type   = mpi_real
  integer(i4b), parameter :: tod_cmplx_type = spc
  integer(i4b), parameter :: cov_type       = sp, mpi_cov_type   = mpi_real

  !integer(i4b), parameter :: rhs_type       = dp, mpi_rhs_type   = mpi_double_precision
  !integer(i4b), parameter :: tod_type       = dp, mpi_tod_type   = mpi_double_precision
  !integer(i4b), parameter :: tod_cmplx_type = dpc
  !integer(i4b), parameter :: cov_type       = dp, mpi_cov_type   = mpi_double_precision

  integer(i4b), parameter :: level_tot  = 0, level_ces = 1, level_ass = 2, level_num = 3
  integer(i4b), parameter :: itype_cov  = 0, itype_rhs = 1, itype_div = 2, &
    & itype_dir = 3, itype_ang = 4, itype_gfit = 5, itype_num = 6
  integer(i4b), parameter :: otype_eqn  = 0, otype_cov = 1, otype_rhs = 2, &
    & otype_bin = 3, otype_rms = 4, otype_cross = 5, otype_gfit = 6, otype_rawrhs = 7, &
    & otype_sigma = 8, otype_rhsdiv = 9, otype_num = 10
  character(len=16), parameter :: levelstrs(0:level_num-1) = (/ "tot", "ces", "ass" /)
  character(len=16), parameter :: itypestrs(0:itype_num-1) = (/ "COV ", "RHS ", "DIV ", "DIR ", "ANG ", "GFIT" /)
  character(len=16), parameter :: otypestrs(0:otype_num-1) = (/ "eqn   ", "cov   ", &
    & "rhs   ", "bin   ", "rms   ", "cross ", "gfit  ", "rawrhs", "sigma ", "rhsdiv" /)
  character(len=16)            :: ofiletypes(0:otype_num-1)

  integer(i4b),     parameter :: T = 1, Q = 2, U = 3, V = 4
  character(len=1), parameter, dimension(4) :: comp_names = [ "T", "Q", "U", "V" ]

  real(dp)                 :: lowpass_max_cap = 1000d0
  real(dp)                 :: highpass_min_cap = 0d0
  integer(i4b)             :: azorder_min_cap = 0

! Hack: Override noise parameters based on premeasured ones
real(dp), parameter :: hack_sig(364) = [ &
1.317908e-05, 1.469674e-05, 1.654284e-05, 1.389365e-05, 8.202692e-06, 1.089238e-05, & 
9.713501e-06, 1.022837e-05, 2.055504e-05, 2.946455e-05, 2.459417e-05, 1.859937e-05, &
1.638689e-05, 6.068068e-06, 7.817006e-06, 8.391395e-06, 3.380211e-05, 4.448666e-05, &
0.000000e+00, 4.099501e-05, 1.367600e-05, 8.632024e-06, 2.126063e-05, 2.201586e-05, &
1.240910e-05, 1.114500e-05, 7.918288e-06, 1.282297e-05, 0.000000e+00, 0.000000e+00, &
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, &
1.847843e-05, 2.399656e-05, 2.306458e-05, 1.580461e-05, 3.997898e-05, 3.595155e-05, &
4.559174e-05, 3.175894e-05, 6.649299e-06, 1.035415e-05, 8.492875e-06, 7.088927e-06, &
1.740201e-05, 2.011996e-05, 2.017162e-05, 1.960141e-05, 1.648340e-05, 1.184718e-05, &
1.591474e-05, 1.045570e-05, 3.246140e-05, 3.130401e-05, 3.574264e-05, 2.304617e-05, &
1.198408e-05, 1.547838e-05, 1.331151e-05, 8.898580e-06, 4.118595e-05, 4.535109e-05, &
4.488675e-05, 3.366597e-05, 3.061696e-05, 3.847982e-05, 4.019694e-05, 2.324471e-05, &
2.942727e-05, 3.866572e-05, 3.791519e-05, 3.302421e-05, 1.947536e-05, 3.123206e-05, &
2.587007e-05, 2.684453e-05, 3.227430e-05, 4.295871e-05, 3.863807e-05, 2.840108e-05, &
2.410769e-05, 2.671082e-05, 2.130747e-05, 1.743619e-05, 7.341083e-06, 6.578493e-06, &
9.520295e-06, 7.637437e-06, 3.511794e-05, 4.093507e-05, 3.362289e-05, 2.350032e-05, &
0.000000e+00, 1.703753e-05, 0.000000e+00, 1.334378e-05, 4.435525e-05, 5.268130e-05, &
4.847549e-05, 3.687810e-05, 4.493112e-06, 4.269168e-06, 5.117619e-06, 3.884560e-06, &
4.777285e-05, 5.061702e-05, 5.771336e-05, 4.038694e-05, 0.000000e+00, 0.000000e+00, &
0.000000e+00, 0.000000e+00, 2.025979e-05, 2.190597e-05, 2.363376e-05, 1.654174e-05, &
1.828566e-05, 1.944179e-05, 1.715122e-05, 1.234305e-05, 1.869658e-05, 2.249574e-05, &
2.862237e-05, 1.769535e-05, 3.943764e-05, 4.634981e-05, 6.745372e-05, 3.419732e-05, &
3.290565e-05, 3.411817e-05, 2.699083e-05, 2.747793e-05, 5.951269e-06, 6.063515e-06, &
5.768541e-06, 7.171299e-06, 7.332064e-06, 7.562677e-06, 7.134662e-06, 5.615656e-06, &
3.550378e-05, 4.623004e-05, 4.681000e-05, 2.614781e-05, 7.778778e-06, 8.196624e-06, &
1.175505e-05, 6.127354e-06, 0.000000e+00, 6.899604e-06, 7.369077e-06, 4.604605e-06, &
0.000000e+00, 2.037906e-05, 1.243301e-05, 0.000000e+00, 0.000000e+00, 1.924666e-05, &
2.577151e-05, 1.892996e-05, 3.644058e-05, 4.976639e-05, 5.328118e-05, 4.074794e-05, &
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.071492e-05, 1.050225e-05, &
1.070885e-05, 8.328035e-06, 6.211085e-06, 7.237996e-06, 6.392518e-06, 4.716800e-06, &
1.831761e-05, 2.110213e-05, 2.995644e-05, 1.750174e-05, 1.057104e-05, 1.202384e-05, &
1.346820e-05, 8.231896e-06, 2.737784e-05, 3.466611e-05, 2.985199e-05, 2.203867e-05, &
1.582412e-05, 1.349832e-05, 1.124682e-05, 1.066446e-05, 1.023342e-05, 1.284164e-05, &
1.228969e-05, 1.071107e-05, 1.795985e-05, 2.058875e-05, 2.061957e-05, 1.743927e-05, &
1.428268e-05, 2.079493e-05, 1.817553e-05, 1.078450e-05, 2.649010e-05, 2.604140e-05, &
3.486670e-05, 2.756030e-05, 2.018352e-05, 1.133781e-05, 2.356460e-05, 1.320916e-05, &
2.009725e-05, 2.442037e-05, 2.164438e-05, 1.569262e-05, 2.987847e-05, 3.168444e-05, &
3.387629e-05, 3.076669e-05, 1.417790e-05, 1.175154e-05, 1.475985e-05, 9.678346e-06, &
1.452784e-05, 1.500877e-05, 1.381703e-05, 9.697027e-06, 9.988320e-06, 1.038404e-05, &
1.234665e-05, 8.958848e-06, 4.773877e-06, 4.112324e-06, 9.776627e-06, 8.674899e-06, &
2.231653e-05, 1.786833e-05, 1.803863e-05, 1.818045e-05, 9.096671e-06, 1.071475e-05, &
9.940564e-06, 7.055880e-06, 3.432780e-06, 3.844260e-06, 3.929116e-06, 3.752448e-06, &
1.139183e-05, 1.359491e-05, 1.428448e-05, 1.368237e-05, 1.013066e-05, 1.247270e-05, &
1.258489e-05, 1.240198e-05, 7.790799e-06, 1.174374e-05, 1.176131e-05, 1.025224e-05, &
1.837579e-05, 1.852452e-05, 2.135598e-05, 1.758956e-05, 9.541828e-06, 1.015217e-05, &
1.077950e-05, 1.912507e-06, 1.850364e-05, 2.832304e-05, 2.833382e-05, 2.308018e-05, &
1.651000e-05, 8.521985e-06, 1.940733e-05, 1.194678e-05, 2.040567e-05, 2.752866e-05, &
2.412517e-05, 2.305934e-05, 2.676013e-06, 3.347170e-06, 2.890651e-06, 3.401678e-06, &
2.066343e-05, 2.352597e-05, 2.306417e-05, 1.799235e-05, 2.647169e-05, 2.998428e-05, &
3.046657e-05, 2.025698e-05, 3.101791e-05, 3.637777e-05, 4.062692e-05, 2.317002e-05, &
3.994508e-05, 3.654483e-05, 4.533967e-05, 3.766179e-05, 6.537603e-06, 8.334330e-06, &
9.441803e-06, 6.205050e-06, 7.681502e-06, 6.276767e-06, 8.401445e-06, 6.111410e-06, &
3.120531e-05, 3.155417e-05, 3.086934e-05, 2.480939e-05, 3.584718e-05, 4.420871e-05, &
4.388720e-05, 3.579160e-05, 4.067310e-05, 4.353573e-05, 3.769130e-05, 3.736240e-05, &
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 6.237876e-06, 6.769932e-06, &
8.499995e-06, 5.542190e-06, 2.043629e-05, 2.716968e-05, 2.679320e-05, 2.438654e-05, &
2.720379e-05, 3.483711e-05, 3.324316e-05, 3.150108e-05, 1.059818e-05, 0.000000e+00, &
0.000000e+00, 8.253163e-06, 3.949930e-05, 0.000000e+00, 0.000000e+00, 3.152268e-05, &
6.308765e-05, 0.000000e+00, 0.000000e+00, 6.961662e-05, 5.867689e-05, 0.000000e+00, &
0.000000e+00, 5.304593e-05, 2.268293e-05, 0.000000e+00, 0.000000e+00, 2.066214e-05, &
0.000000e+00, 0.000000e+00, 0.000000e+00, 6.552727e-06 ]

  type common_info
     ! Mpi
     integer(i4b)              ::             myid,       nproc
     integer(i4b)              :: comm_cov,   myid_cov,   nproc_cov
     integer(i4b)              :: comm_group, myid_group, nproc_group
     ! Analysis parameters
     character(len=64)         :: target
     integer(i4b)              :: cid, nside, coord
     integer(i4b)              :: ncomp
     integer(i4b), allocatable :: comps(:)
     ! Status monitor stuff
     type(status_file)         :: status
     type(benchmarker)  :: bench
  end type common_info

  type mapdata
     ! The index ordering here looks inefficient. Somebody should
     ! investigate this.
     real(rhs_type), pointer, dimension(:,:,:)   :: rhs => null() ! (npix,ncomp,nmap)
     real(rhs_type), pointer, dimension(:,:,:)   :: div => null() ! (npix,ncomp,ncomp)
     real(cov_type), pointer, dimension(:,:,:,:) :: cov => null() ! (npix,ncomp,npix,ncomp)
     real(rhs_type), pointer, dimension(:,:)     :: ang => null() ! (npix,3)
     real(rhs_type), pointer, dimension(:,:)     :: dir => null() ! (npix,3)
     integer(i8b) :: n, ncomp, nmap, col_start, col_stop
     logical(lgt) :: enabled = .false.
     logical(lgt) :: active(0:itype_num-1) = .false., delegated(0:itype_num-1) = .false.
     ! These constants map to component indices.
     integer(i4b) :: T, Q, U, S
  end type mapdata

  type jk_struct
     character(len=512) :: object, name
     integer(i4b), allocatable, dimension(:) :: nbin
     real(dp), allocatable, dimension(:,:) :: pte, chisq, chi, chisq_max
  end type jk_struct
  
  type report_struct
     real(dp)     :: tot_obs_time, tot_diode_time, tot_cpu_time
     integer(i4b) :: num_diode_ces, num_ces, orig_num_diode_ces
  end type report_struct

  interface add_mapdata
     module procedure add_mapdata_single, add_mapdata_multi
  end interface

  interface clear_mapdata
     module procedure clear_mapdata_single, clear_mapdata_multi
  end interface

  interface free_mapdata
     module procedure free_mapdata_single, free_mapdata_multi
  end interface

contains

  subroutine setup_filetypes(mapext)
    implicit none
    character(len=*) :: mapext
    ofiletypes = [ "unf ", "unf ", mapext, mapext, mapext, mapext, &
     & mapext, "unf ", mapext, "hdf" ]
  end subroutine

  function parse_output_options(string) result(res)
    implicit none
    character(len=*)     :: string
    logical(lgt) :: res(0:otype_num-1)
    integer(i4b) :: otype
    do otype = 0, otype_num-1
       res(otype) = has_token(otypestrs(otype), string, ", ")
    end do
  end function

  subroutine initialize_report_struct(summary)
    implicit none
    type(report_struct) :: summary
    summary%tot_obs_time   = 0.d0
    summary%tot_cpu_time   = 0.d0
    summary%tot_diode_time = 0.d0
    summary%num_diode_ces  = 0
    summary%num_ces        = 0
    summary%orig_num_diode_ces = 0
  end subroutine initialize_report_struct

  subroutine update_report_struct(summary, cpu_time, obs_time, num_diode_ces, num_ces, diode_time,orig_num_diode_ces)
    implicit none

    type(report_struct), intent(inout)           :: summary
    real(dp),            intent(in),    optional :: cpu_time, obs_time, diode_time
    integer(i4b),        intent(in),    optional :: num_diode_ces, num_ces,orig_num_diode_ces

    if (present(cpu_time))      summary%tot_cpu_time   = summary%tot_cpu_time   + cpu_time
    if (present(obs_time))      summary%tot_obs_time   = summary%tot_obs_time   + obs_time
    if (present(diode_time))    summary%tot_diode_time = summary%tot_diode_time + diode_time
    if (present(num_diode_ces)) summary%num_diode_ces  = summary%num_diode_ces  + num_diode_ces
    if (present(num_ces))       summary%num_ces        = summary%num_ces        + num_ces
    if (present(orig_num_diode_ces)) summary%orig_num_diode_ces  = summary%orig_num_diode_ces  + orig_num_diode_ces

  end subroutine update_report_struct

  ! Set up a mapdata structure. If data2 is passed, anything already
  ! allocated in data2 will be delgated to, and local allocation will
  ! only happen for the those missing in data2
  subroutine setup_mapdata(data, n, nmap, info, itypes, data2, delegate)
    implicit none

    type(common_info)                  :: info
    type(mapdata)                      :: data
    integer(i4b),           intent(in) :: n, nmap
    logical(lgt)                       :: itypes(0:)
    type(mapdata),optional, target     :: data2
    logical(lgt), optional             :: delegate(0:)
    type(mapdata),target               :: data_null
    type(mapdata),pointer              :: data_foo
    integer(i4b) :: ierr, i

    call free_mapdata(data)
    data_foo => data_null; if(present(data2)) data_foo => data2

    data%active    = itypes
    data%enabled   = any(itypes)
    data%n         = n
    data%nmap      = nmap
    data%ncomp     = info%ncomp

    data%col_start = info%myid_cov*n/info%nproc_cov+1
    data%col_stop  = (info%myid_cov+1)*n/info%nproc_cov
    if (info%myid_cov == info%nproc_cov-1) data%col_stop = n

    data%delegated = data_foo%active
    if(present(delegate)) data%delegated = data%delegated .and. delegate

    do i = 0, itype_num-1
       if(.not. data%active(i)) cycle
          if(data%delegated(i)) then
          select case(i)
             case(itype_rhs); data%rhs => data_foo%rhs
             case(itype_cov); data%cov => data_foo%cov
             case(itype_div); data%div => data_foo%div
             case(itype_ang); data%ang => data_foo%ang
             case(itype_dir); data%dir => data_foo%dir
          end select
       else
          select case(i)
             case(itype_rhs); allocate(data%rhs(n,info%ncomp,data%nmap))
             case(itype_cov); allocate(data%cov(n+1,info%ncomp,data%col_start:data%col_stop,info%ncomp))
             case(itype_div); allocate(data%div(n,info%ncomp,info%ncomp))
             case(itype_ang); allocate(data%ang(n,3))
             case(itype_dir); allocate(data%dir(n,3))
          end select
       end if
    end do
    call clear_mapdata(data)

  end subroutine

  ! Add data1 to data2, but only if data1 actually has contents itself
  subroutine add_mapdata_single(data1, data2)
    implicit none
    type(mapdata)  :: data1, data2
    integer(i4b)   :: i, j
    do i = 0, itype_num-1
       if(.not. data1%active(i) .or. .not. data2%active(i) .or. data1%delegated(i)) cycle
       select case(i)
          case(itype_rhs); data2%rhs = data2%rhs + data1%rhs
          case(itype_cov)
             do j = data1%col_start, data1%col_stop
                data2%cov(:,:,j,:) = data2%cov(:,:,j,:) + data1%cov(:,:,j,:)
             end do
          case(itype_div); data2%div = data2%div + data1%div
          case(itype_ang); data2%ang = data2%ang + data1%ang
          case(itype_dir); data2%dir = data2%dir + data1%dir
       end select
    end do
  end subroutine

  subroutine add_mapdata_multi(data1, data2)
    implicit none
    type(mapdata) :: data1(:), data2(:)
    integer(i4b)  :: set
    do set = 1, size(data1)
       call add_mapdata_single(data1(set), data2(set))
    end do
  end subroutine

  subroutine clear_mapdata_single(data)
    implicit none
    type(mapdata)  :: data
    integer(i4b)   :: i
    do i = 0, itype_num-1
       if(.not. data%active(i) .or. data%delegated(i)) cycle
       select case(i)
          case(itype_rhs); data%rhs = 0
          case(itype_cov); data%cov = 0
          case(itype_div); data%div = 0
          case(itype_ang); data%ang = 0
          case(itype_dir); data%dir = 0
       end select
    end do
  end subroutine

  subroutine clear_mapdata_multi(data)
    implicit none
    type(mapdata) :: data(:)
    integer(i4b)  :: set
    do set = 1, size(data)
       call clear_mapdata_single(data(set))
    end do
  end subroutine

  subroutine free_mapdata_single(data)
    implicit none
    type(mapdata)        :: data
    integer(i4b)   :: i
    do i = 0, itype_num-1
       if(.not. data%active(i) .or. data%delegated(i)) cycle
       select case(i)
          case(itype_rhs); deallocate(data%rhs)
          case(itype_cov); deallocate(data%cov)
          case(itype_div); deallocate(data%div)
          case(itype_ang); deallocate(data%ang)
          case(itype_dir); deallocate(data%dir)
       end select
    end do
    data%active    = .false.
    data%delegated = .false.
    nullify(data%rhs)
    nullify(data%cov)
    nullify(data%div)
    nullify(data%ang)
    nullify(data%dir)
  end subroutine

  subroutine free_mapdata_multi(data)
    implicit none
    type(mapdata) :: data(:)
    integer(i4b)  :: set
    do set = 1, size(data)
       call free_mapdata_single(data(set))
    end do
  end subroutine

  subroutine init_mon(info, filename)
    implicit none
    character(len=*)  :: filename
    type(common_info) :: info
    call init_status(info%status, filename)
  end subroutine

  subroutine free_mon(info)
    implicit none
    type(common_info):: info
    call free_status(info%status)
  end subroutine

  subroutine update_mon(info, desc)
    implicit none
    type(common_info):: info
    character(len=*) :: desc
    call update_status(info%status, desc)
  end subroutine

  subroutine get_fft3_magic_number(k, filename)
    implicit none

    integer(i4b),     intent(inout)          :: k
    character(len=*), intent(in),   optional :: filename

    integer(i4b)       :: i, unit, pos
    integer(i4b), save :: n
    integer(i4b), allocatable, dimension(:), save :: list

    if (present(filename)) then
       unit = getlun()
       open(unit, file=trim(filename))
       read(unit,*) n
       if (allocated(list)) deallocate(list)
       allocate(list(n))
       do i = 1, n
          read(unit,*) list(i)
       end do
    end if

    if (.not. allocated(list)) return

    pos = locate(list, k)
    if (pos > 1 .and. pos < n) then
       if (k - list(pos) < k/100 .and. k-list(pos) >= 0) then
          k = list(pos)
       end if
    end if
  end subroutine get_fft3_magic_number

  subroutine get_assembly_data(assembly, info, l3data)
    implicit none
    type(quiet_assembly)        :: assembly
    type(common_info)           :: info
    type(lx_struct)             :: l3data

    integer(i4b) :: i, j, h, k, m, n, nrow, ndi, di, partner, mod, diode, nhorn, d, p
    integer(i4b) :: ngroup, nbin, nfft
    real(dp)     :: mjd, phi, theta, psi, center(2)
    type(planck_rng) :: rng_handle
    type(patch_info) :: pinfo

    ! Define local parameters
    nrow            = size(l3data%time)
    call get_fft3_magic_number(nrow)  ! Truncate time stream to closest FFT3 magic number
    nfft            = nrow/2+1
    ndi             = assembly%num_diodes
    ! Allocate data structures
    allocate(assembly%time(nrow))
    allocate(assembly%tod(nrow,ndi))
    allocate(assembly%orig_point(3,nrow))
    allocate(assembly%diode_info(ndi))

    ! Fill in the correct pointing and TOD information
    assembly%samprate   = l3data%samprate
    assembly%scanfreq   = l3data%scanfreq
    assembly%time       = l3data%time
    assembly%tod        = l3data%tod(:,assembly%diodes_abs)
    assembly%decimation = l3data%decimation
    assembly%numsamp    = nrow
    assembly%nside      = info%nside
    assembly%coordinate_system = info%coord
    assembly%orig_point(1,:) =      - l3data%orig_point(1,:)
    assembly%orig_point(2,:) = pi/2 - l3data%orig_point(2,:)
    assembly%orig_point(3,:) =        l3data%orig_point(3,:)

    center = 0
    if(get_patch_info(info%target, pinfo)) then
       center = get_patch_pos_single(pinfo, l3data%time(1), coord_gal)
    end if

    ! The pointing must be converted to the right system
    do i = 1, ndi
       m     = assembly%diodes(i,1)
       d     = assembly%diodes(i,2)
       di    = assembly%diodes_abs(i)
       nhorn = size(quiet_horns(m)%groups)
       allocate(assembly%diode_info(i)%point(3,nrow,nhorn))
       allocate(assembly%diode_info(i)%pix(nhorn,nrow))
       allocate(assembly%diode_info(i)%gain(nrow))
       allocate(assembly%diode_info(i)%gamp(nhorn))
       assembly%diode_info(i)%gamp = quiet_horns(m)%amps

       do k = 1, nrow
          mjd   = l3data%time(k)
          do h = 1, nhorn
             p = quiet_horns(m)%groups(h)
             ! Must manually add diode angle, as level2-files do not contain them.
             call coord_convert(l3data%coord_sys, real(l3data%point(1,k,p+1),dp), &
               & real(l3data%point(2,k,p+1),dp), real(l3data%point(3,k,p+1),dp) + &
               & get_diode_angle(di), info%coord, phi, theta,&
               & psi, mjd, m, d, phic=center(1), thetac=center(2))
             assembly%diode_info(i)%point(:,k,h) = [ phi, theta, psi ]
             call ang2pix_nest(info%nside, theta, phi, assembly%diode_info(i)%pix(h,k))
          end do
       end do
    end do

    ! Extract the gain from the file, and interpolate it to full resolution
    do i = 1, ndi
       call lin_interpol(l3data%time_gain, real(l3data%gain(:,assembly%diodes_abs(i)),dp), &
            & l3data%time(1:nrow), assembly%diode_info(i)%gain)
    end do

    ! And fetch the precomputed noise stuff
    allocate(assembly%sigma0(ndi), assembly%alpha(ndi), assembly%fknee(ndi))
    allocate(assembly%corr(nfft,ndi,ndi))
    !write(*,*) 'Warning: Multiplying sigma0 by 1.02!!!'
! TMR: test stuff
!!$assembly%sigma0(:) = 1.d0
!!$assembly%fknee(:) = 1.d-5
!!$assembly%alpha(:) = -5.d0


    assembly%sigma0 = l3data%sigma0(assembly%diodes_abs) !* 1.02d0 !HACK
    assembly%alpha  = l3data%alpha (assembly%diodes_abs)
    assembly%fknee  = l3data%fknee (assembly%diodes_abs)
    call interpolate_corr(l3data%corr_freqs, l3data%samprate, &
     & l3data%corr(:,assembly%diodes_abs,assembly%diodes_abs), assembly%corr)

    ! Fetch filter information
    allocate(assembly%fft_low_freq(ndi), assembly%fft_low_alpha(ndi))
    allocate(assembly%fft_high_freq(ndi), assembly%fft_high_alpha(ndi))
    allocate(assembly%az_order(ndi))
    assembly%fft_low_freq   = minval(l3data%filter_par(assembly%diodes_abs,FILTER_LOW_NU))
    assembly%fft_low_alpha  = minval(l3data%filter_par(assembly%diodes_abs,FILTER_LOW_ALPHA))
    assembly%fft_high_freq  = maxval(l3data%filter_par(assembly%diodes_abs,FILTER_HIGH_NU_SCAN) * assembly%scanfreq)
    assembly%fft_high_alpha = minval(l3data%filter_par(assembly%diodes_abs,FILTER_HIGH_ALPHA))
    assembly%az_order       = maxval(l3data%filter_par(assembly%diodes_abs,FILTER_AZ_ORDER))

    ! Override filters
    assembly%fft_high_freq =max(assembly%fft_high_freq,assembly%scanfreq*highpass_min_cap)
    assembly%fft_low_freq  =min(assembly%fft_low_freq,lowpass_max_cap)
    assembly%az_order      =max(assembly%az_order,azorder_min_cap)

! Hack override noise parameters
if(info%target == "tau_a") then
assembly%sigma0 = hack_sig(assembly%diodes_abs)
assembly%alpha  = -5
assembly%fknee  = 0.005
if(info%myid == 0) write(*,*) "Taua noise override"
end if

    ! For tod2map to work, sigma0 and gain must be nonzero (especially sigma0),
    ! so check this here.
    if (any(assembly%sigma0 == 0)) then
       write(*,*) 'sigma0', assembly%sigma0
       write(*,*) 'sigma0 is zero for myid =', info%myid
       stop
    end if
    !call assert(.not. any(assembly%sigma0 == 0), "Zero sigma0 in get_assembly_data!")
  end subroutine

  subroutine interpolate_corr(freqs, srate, corr_in, corr_out)
    implicit none
    real(dp),       intent(in)  :: freqs(:), srate, corr_in(:,:,:)
    real(dp),       intent(out) :: corr_out(:,:,:)
    real(dp),       allocatable :: x(:), ox(:), variance(:), y(:), y2(:)
    integer(i4b)                :: d1, d2, ndi, n, i, m
    ndi = size(corr_in, 2)
    n   = size(corr_out,1)
    m   = size(freqs)-1
    allocate(x(m),ox(n),variance(m),y(m),y2(m))
    ! Assign measured corr to middle of interaval
    x = (freqs(2:m+1)+freqs(1:m))/2
    do i = 1, n
       ox(i) = ind2freq(i, srate, n)
    end do
    variance = srate**2/(freqs(2:m+1)-freqs(1:m))**2
    ! Our bins are exponentially spaced, to it makes sense to use a
    ! logarithmic x
    x = log(x); ox = log(ox)
    ox(1) = ox(2)-(ox(3)-ox(2))
    do d1 = 1, ndi
       corr_out(:,d1,d1) = 1
       do d2 = d1+1, ndi
          y = corr_in(2:m,d1,d2)
          !call smooth_spline("inv_var", 1d10, x, y, 1d30, 1d30, y2, variance)
          call spline_plain(x, y, 1d30, 1d30, y2)
          do i = 1, n
             corr_out(i,d1,d2) = splint_plain(x, y, y2, ox(i))
          end do
          corr_out(:,d2,d1) = corr_out(:,d1,d2)
       end do
    end do
    deallocate(x,ox,variance,y,y2)
  end subroutine

  subroutine extract_pixels(map2mask, pixels, mask)
    implicit none
    integer(i4b)           :: map2mask(0:)
    logical(lgt), optional :: mask(:)
    integer(i4b), dimension(:), allocatable :: pixels

    integer(i4b) :: i, n, m

    if(.not. present(mask)) then
       n = maxval(map2mask)
       allocate(pixels(n))
       do i = 0, size(map2mask)-1
          if(map2mask(i) <= 0) cycle
          pixels(map2mask(i)) = i
       end do
    else
       n = maxval(map2mask, mask=mask)
       allocate(pixels(n))
       do i = 0, size(map2mask)-1
          if(.not. mask(i)) cycle
          if(map2mask(i) <= 0) cycle
          pixels(map2mask(i)) = i
       end do
    end if
  end subroutine

  subroutine allocate_jk_struct(jk, object, name, nmap)
    implicit none

    type(jk_struct)  :: jk
    integer(i4b)     :: nmap
    character(len=*) :: object, name

    jk%object = object
    jk%name   = name
    allocate(jk%pte(nmap,3),jk%chisq(nmap,3),jk%chi(nmap,3),jk%chisq_max(nmap,3),jk%nbin(3))
    jk%nbin   = 0

  end subroutine allocate_jk_struct

  subroutine deallocate_jk_struct(jk)
    implicit none

    type(jk_struct)  :: jk

    jk%object = ''
    jk%name   = ''
    jk%nbin   = 0
    if (allocated(jk%pte))       deallocate(jk%pte)    
    if (allocated(jk%chisq))     deallocate(jk%chisq)    
    if (allocated(jk%chi))       deallocate(jk%chi)    
    if (allocated(jk%chisq_max)) deallocate(jk%chisq_max)    
    if (allocated(jk%nbin))      deallocate(jk%nbin)    
  end subroutine deallocate_jk_struct

  subroutine output_map_general(map, outname, nside, ordering, map2mask, sparse)
    implicit none
    real(dp),     dimension(:,:,:) :: map
    character(len=*)               :: outname
    integer(i4b)                   :: ordering, nside
    integer(i4b), dimension(0:)    :: map2mask
    logical(lgt),               optional   :: sparse

    integer(i4b), dimension(:),   allocatable :: pixels
    real(sp),     dimension(:,:,:), allocatable :: map_full
    integer(i4b)             :: i, j, k, m, n, npix
    logical(lgt)             :: sparse_

    sparse_ = .true.; if(present(sparse)) sparse_ = sparse

    if(sparse_) then
       ! Expand multiresolution
       n = count(map2mask > 0)
       allocate(map_full(n,3,size(map,3)), pixels(n))
       j = 0
       do i = 0, size(map2mask)-1
          if(map2mask(i) <= 0) cycle
          j = j+1
          map_full(j,:,:) = map(map2mask(i),:,:)
          pixels(j)       = i
       end do
       if(ordering /= NEST) then
          do i = 1, n; call nest2ring(nside, pixels(i), pixels(i)); end do
       end if
       call mkdirs(outname, .true.)
       call write_map(map_full, pixels, nside, ordering, outname)
       deallocate(map_full, pixels)
    else
       ! Not sparse, so fill the full pixel sets:
       npix = size(map2mask)
       allocate(map_full(0:npix-1,3,size(map,3)))
       map_full = hpx_dbadval
       do i = 0, npix-1
          if(map2mask(i) <= 0) cycle
          map_full(i,:,:) = map(map2mask(i),:,:)
       end do
       if(ordering /= 2) then
          do i = 1, size(map_full,3)
             call convert_nest2ring(nside, map_full(:,:,i))
          end do
       end if
       call write_map(map_full, ordering, outname)
       deallocate(map_full)
    end if
  end subroutine

  subroutine toks2inds(str,names,inds)
    implicit none
    character(len=*)                                     :: str, names(:)
    integer(i4b),              allocatable, dimension(:) :: inds
    character(len=len(names)), allocatable, dimension(:) :: tokens
    integer(i4b) :: i, j, n
    if(allocated(inds)) deallocate(inds)
    n = num_tokens(str,",")
    allocate(tokens(n),inds(n))
    call get_tokens(str,",",tokens)
    inds = 0
    do i = 1, n
       do j = 1, size(names)
          if(names(j) == tokens(i)) exit
       end do
       if(j <= size(names)) inds(i) = j
    end do
    deallocate(tokens)
  end subroutine

!  ! Generate a test assembly from scratch. We only generate what we
!  ! need for our tests.
!  subroutine get_assembly_test(assembly)
!    implicit none
!    type(quiet_assembly)        :: assembly
!    integer(i4b), parameter     :: n = 16, nf = n/2+1, ndi = 2, ncorr = 1
!    integer(i4b)                :: i
!    call free_assembly_data(assembly)
!    allocate(assembly%diodes_abs(2), assembly%diodes(2,2))
!    assembly%diodes_abs        = [1,2]
!    assembly%diodes            = reshape([1,2,1,1],[2,2])
!    assembly%num_diodes        = size(assembly%diodes_abs)
!    assembly%numsamp           = n
!    assembly%samprate          = 25
!    assembly%scanfreq          = 0.1
!    assembly%noise_corr_length = ncorr
!    allocate(assembly%N_corr_real(-ncorr:ncorr,ndi,ndi))
!    allocate(assembly%N_corr_F_sq_real(-ncorr:ncorr,ndi,ndi))
!    assembly%N_corr_real = reshape([&
!     & -2, 4, -2 &
!     &  1,-2,  1,&
!     &  1,-2,  1,&
!     & -2, 4, -2], [2*ncorr+1,ndi,ndi])
!    assembly%N_corr_F_sq_real = 2*assembly%N_corr_real
!    allocate(assembly%diode_info(ndi))
!    do i = 1, ndi
!       allocate(assembly%diode_info(i)%pix(n))
!       assembly(diode_info(i)%pix = [1,1,1,2,1,1,2,2,1,2,2,2,2,2,2,2]
!       allocate(assembly%diode_info
!
!    end do
!
!
!     real(dp),    allocatable, dimension(:,:,:) :: point ! ((phi,theta,psi),nsamp,ngroup)
!     integer(i4b),allocatable, dimension(:,:)   :: pix   ! (ngroup,nsamp)
!     real(dp),    allocatable, dimension(:)     :: gamp  ! (ngroup)
!     real(dp),    allocatable, dimension(:)     :: gain  ! (nsamp)
!
!
!
!             p1 = data%map2mask(assembly%diode_info(d1)%pix(h1,i1))
!             call bench_start(info%bench, bench_cov1)
!                   N_phase   = data%phase(d1)%val(:,h1,i1) * corr(d1,d2,k)
!                   N_phase_sq= data%phase(d1)%val(:,h1,i1) * corr_sq(d1,d2,k)
!                      p2 = data%map2mask(assembly%diode_info(d2)%pix(h2,i2))
!                      if(.not. (p2 >= out%col_start .and. p2 <= out%col_stop)) cycle
!                      do c1 = 1, info%ncomp
!                         do c2 = 1, info%ncomp
!                             & N_phase(c1) * data%phase(d2)%val(c2,h2,i2)
!                             & N_phase_sq(c1) * data%phase(d2)%val(c2,h2,i2)
!             call bench_stop(info%bench, bench_cov1)
!             call bench_start(info%bench, bench_cov2)
!             do col = out%col_start, out%col_stop
!
!
!
!
!    ! The pointing must be converted to the right system
!    do i = 1, ndi
!       m     = assembly%diodes(i,1)
!       d     = assembly%diodes(i,2)
!       di    = assembly%diodes_abs(i)
!       nhorn = size(quiet_horns(m)%groups)
!       allocate(assembly%diode_info(i)%point(3,nrow,nhorn))
!       allocate(assembly%diode_info(i)%pix(nhorn,nrow))
!       allocate(assembly%diode_info(i)%gain(nrow))
!       allocate(assembly%diode_info(i)%gamp(nhorn))
!       assembly%diode_info(i)%gamp = quiet_horns(m)%amps
!
!       do k = 1, nrow
!          mjd   = l3data%time(k)
!          do h = 1, nhorn
!             p = quiet_horns(m)%groups(h)
!             ! Must manually add diode angle, as level2-files do not contain them.
!             call coord_convert(l3data%coord_sys, real(l3data%point(1,k,p+1),dp), &
!               & real(l3data%point(2,k,p+1),dp), real(l3data%point(3,k,p+1),dp) + &
!               & get_diode_angle(di), info%coord, phi, theta,&
!               & psi, mjd, m, d)
!             assembly%diode_info(i)%point(:,k,h) = [ phi, theta, psi ]
!             call ang2pix_nest(info%nside, theta, phi, assembly%diode_info(i)%pix(h,k))
!          end do
!       end do
!    end do
!
!    ! Extract the gain from the file, and interpolate it to full resolution
!    do i = 1, ndi
!       call lin_interpol(l3data%time_gain, real(l3data%gain(:,assembly%diodes_abs(i)),dp), l3data%time(1:nrow), assembly%diode_info(i)%gain)
!    end do
!
!    ! And fetch the precomputed noise stuff
!    allocate(assembly%sigma0(ndi), assembly%alpha(ndi), assembly%fknee(ndi))
!    allocate(assembly%corr(nfft,ndi,ndi))
!    assembly%sigma0 = l3data%sigma0(assembly%diodes_abs)
!    assembly%alpha  = l3data%alpha (assembly%diodes_abs)
!    assembly%fknee  = l3data%fknee (assembly%diodes_abs)
!    call interpolate_corr(l3data%corr_freqs, l3data%samprate, &
!     & l3data%corr(:,assembly%diodes_abs,assembly%diodes_abs), assembly%corr)
!
!    ! Fetch filter information
!    allocate(assembly%fft_low_freq(ndi), assembly%fft_low_alpha(ndi))
!    allocate(assembly%fft_high_freq(ndi), assembly%fft_high_alpha(ndi))
!    allocate(assembly%az_order(ndi))
!    assembly%fft_low_freq   = minval(l3data%filter_par(assembly%diodes_abs,FILTER_LOW_NU))
!    assembly%fft_low_alpha  = minval(l3data%filter_par(assembly%diodes_abs,FILTER_LOW_ALPHA))
!    assembly%fft_high_freq  = maxval(l3data%filter_par(assembly%diodes_abs,FILTER_HIGH_NU_SCAN) * assembly%scanfreq)
!    assembly%fft_high_alpha = minval(l3data%filter_par(assembly%diodes_abs,FILTER_HIGH_ALPHA))
!    assembly%az_order       = maxval(l3data%filter_par(assembly%diodes_abs,FILTER_AZ_ORDER))
!
!    ! For tod2map to work, sigma0 and gain must be nonzero (especially sigma0),
!    ! so check this here.
!    call assert(.not. any(assembly%sigma0 == 0), "Zero sigma0 in get_assembly_data!")
!  end subroutine



end module
