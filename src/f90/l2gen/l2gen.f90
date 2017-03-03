! Given a runlist, extract CESes from level1 files, producing level2-files
! in hdf format.
program l2gen
  use fake_mpi_mod
  use l1_read_mod
  use quiet_utils
  use quiet_fileutils
  use quiet_detector_mod
  use quiet_ces_mod
  use quiet_typeb_mod
  implicit none

  character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile
  integer(i4b)         :: i, j, k, l, m, n, cnum, nces, unit, myid, nproc, ierr, nmod
  integer(i4b)         :: mstep, i2, phase_offset, decimation, mod, di, nsamp, status
  integer(i4b)         :: debug
  logical(lgt)         :: exist, reprocess, check_existing, gonext
  real(dp)             :: timing_offset, mjd_tol, mjd(2), dt_error, samprate, scanfreq
  integer(i4b), dimension(:),   allocatable :: mods
  real(sp),     dimension(:,:), allocatable :: buf
  type(quiet_ces_info) :: ces
  type(point_struct)   :: point_in, point
  type(data_struct)    :: data_in, data
  type(hk_struct)      :: hkeep
  type(hdf_file)       :: file

  call getarg(1, parfile)
  call get_parameter(unit, parfile, 'TIMING_OFFSET',       par_dp=timing_offset)
  call get_parameter(unit, parfile, 'L1_DECIMATION',       par_int=decimation)
  call get_parameter(unit, parfile, 'REPROCESS_ALL_FILES', par_lgt=reprocess)
  call get_parameter(unit, parfile, 'DEBUG',               par_int=debug)

  check_existing = .true.
  mjd_tol  = 10d0/60/60/24
  samprate = 100d0

  call initialize_ces_mod(parfile)
  call init_detector_mod(parfile)
  call initialize_typeb_mod(parfile)
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call dset(id=myid,level=debug)

  nmod = size(quiet_horns)
  allocate(mods(nmod))
  do i = 1, nmod; mods(i) = i-1; end do
  mstep = 10

  nces = get_num_ces()
  do cnum = 1+myid, nces, nproc
     call get_ces_info(cnum, ces)
     tmpfile = trim(ces%l2file) // ".part"
     inquire(file=tmpfile,exist=exist)
     if(exist) then
        write(*,fmt="(i3,a,2i5)") myid, " found incomplete run:", cnum, ces%cid
        call rm(tmpfile)
     end if
     inquire(file=ces%l2file,exist=exist)
     if(exist .and. .not. reprocess) then
        if(check_existing) then
           ! Check if the file is as long as it should be
           call get_l2_time_stats(ces%l2file, mjd, dt_error)
           gonext = dt_error < 0.25 .and. (mjd(2) - ces%mjd(2) < mjd_tol .or. ces%mjd(1) - mjd(1) < mjd_tol)
        else
           gonext = .true.
        end if
        if(gonext) then
           write(*,fmt="(i3,a,2i5,a)") myid, " skipping already finished ces:", cnum, ces%cid
           cycle
        end if
     end if
     write(*,fmt="(i3,a,i4,a)") myid, " processing ces ", ces%cid, " (" // trim(itoa(cnum)) // "/" // trim(itoa(nces)) // ")"
     call dmem("ces start")

     ! Read in the level1 data and massage it
     call l1_read_range(ces%l1files, ces%mjd, mods, point=point_in, hk=hkeep, &
          & data=data_in, sel=l1sel(point_std,data_std,hk_all,[data_scale]), status=status)
     call dmem("l1 read range")
     if(status /= 0) then
        write(*,fmt="(i3,a,i4,a)") myid, " error reading ", ces%cid, " skipping"
        cycle
     end if

     call typeb_correct(data_in); call dmem("typeb correct")
     if(timing_offset /= 0) call offset_pointing_interpol(point_in, timing_offset)
     call dmem("offset pointing interpol")

     phase_offset = get_phase_offset(data_in%phase)
     call double_demodulate(pointing_in=point_in, pointing_out=point, offset=phase_offset, nstep=decimation*2)
     call double_demodulate(data_in    =data_in,  data_out     =data, offset=phase_offset, nstep=decimation*2)
     call dmem("double demodulate")

     call deallocate_point_struct(point_in)
     call deallocate_data_struct(data_in)

     call pad_housekeeping(hkeep%bias)
     call pad_housekeeping(hkeep%cryo)
     call pad_housekeeping(hkeep%peri)
     call pad_housekeeping(hkeep%encl)
     call dmem("pad housekeeping")

     ! Write it to hdf files
     call mkdirs(tmpfile, .true.)

     call open_hdf_file(tmpfile, file, "w")
     call write_hdf     (file, "decimation", decimation)
     call write_hdf     (file, "samprate",   samprate/2/decimation)
     call write_hdf     (file, "time",       point%time)

     nsamp = size(point%time)
     allocate(buf(3,nsamp))
     buf(1,:) = point%encoder_azimuth
     buf(2,:) = point%encoder_elevation
     buf(3,:) = point%encoder_deck
     call write_hdf     (file, "orig_point", buf)
     deallocate(buf)

     call create_hdf_set(file, 'tod',        [nsamp,4*nmod], H5T_IEEE_F32LE)
     do mod = 1, nmod
        do di = 0, 3
           call write_hdf(file, 'tod', slice([1,-1],[(mod-1)*4+di+1]), data%RQ(mod)%demod(di,:))
        end do
     end do
     call create_hdf_set(file, 'tp',         [nsamp,4*nmod], H5T_IEEE_F32LE)
     do mod = 1, nmod
        do di = 0, 3
           call write_hdf(file, 'tp',  slice([1,-1],[(mod-1)*4+di+1]), data%RQ(mod)%avg(di,:))
        end do
     end do

     ! And the housekeeping
     call create_hdf_group(file, "bias")
     call write_hdf(file, "bias/time", hkeep%bias%time)
     call write_hdf(file, "bias/value", hkeep%bias%value)
     ! call write_hdf(file, "bias/name", hkeep%bias%name)

     call create_hdf_group(file, "cryo")
     call write_hdf(file, "cryo/time", hkeep%cryo%time)
     call write_hdf(file, "cryo/value", hkeep%cryo%value)
     ! call write_hdf(file, "cryo/name", hkeep%cryo%name)

     call create_hdf_group(file, "encl")
     call write_hdf(file, "encl/time", hkeep%encl%time)
     call write_hdf(file, "encl/value", hkeep%encl%value)
     ! call write_hdf(file, "encl/name", hkeep%encl%name)

     call create_hdf_group(file, "peri")
     call write_hdf(file, "peri/time", hkeep%peri%time)
     call write_hdf(file, "peri/value", hkeep%peri%value)
     ! call write_hdf(file, "peri/name", hkeep%peri%name)

     call close_hdf_file(file)
     call dmem("write")

     call mkdirs(trim(ces%l2file), .true.)
     call mv(tmpfile, ces%l2file)

     call deallocate_point_struct(point)
     call deallocate_data_struct(data)
     call deallocate_hk_struct(hkeep)
     call free_ces_info(ces)
  end do
  call mpi_finalize(ierr)

contains

  function get_phase_offset(phase) result(offset)
    implicit none
    integer(i2b) :: phase(:)
    integer(i4b) :: offset
    offset = 0
    if(size(phase) < 1) return
    if(phase(1) == 1)   return
    offset = 1
  end function

  subroutine offset_pointing_interpol(p, offset)
    implicit none
    type(point_struct) :: p
    real(dp) :: offset, shift, dt
    real(dp), parameter :: sample_rate = 100d0, maxang = 2*pi
    real(dp), dimension(:), allocatable :: mode_tmp
    shift = -offset*sample_rate

    ! Time must also be shifted, but we don't have to physically shift
    ! the array in this case, as we know what will happen to it.
    dt = -offset/60/60/24
    p%start_mjd = p%start_mjd + dt
    if(allocated(p%time)) p%time = p%time + dt

    ! Shift the encoder (and command, I guess?). Do not shift time
    if(allocated(p%encoder_azimuth)) &
      & call shift_angle_array(p%encoder_azimuth, shift, maxang)
    if(allocated(p%encoder_elevation)) &
      & call shift_angle_array(p%encoder_elevation, shift, maxang)
    if(allocated(p%encoder_deck)) &
      & call shift_angle_array(p%encoder_deck, shift, maxang)
    if(allocated(p%command_azimuth)) &
      & call shift_angle_array(p%command_azimuth, shift, maxang)
    if(allocated(p%command_elevation)) &
      & call shift_angle_array(p%command_elevation, shift, maxang)
    if(allocated(p%command_deck)) &
      & call shift_angle_array(p%command_deck, shift, maxang)
    ! Mode is a bit tricky. Turn it into an array of reals, interpoalte
    ! And then truncate back into integers
    if(allocated(p%mode)) then
       allocate(mode_tmp(size(p%mode)))
       mode_tmp = real(p%mode, dp)
       call shift_array_interpol(mode_tmp, shift)
       p%mode = nint(mode_tmp)
       deallocate(mode_tmp)
    end if
  end subroutine

  subroutine pad_housekeeping(kt)
    implicit none
    type(hk_type) :: kt
    if(kt%n_t == 0) then
       kt%n_t = 1
       deallocate(kt%time, kt%value)
       allocate(kt%time(1), kt%value(1, kt%n))
       kt%time = 0
       kt%value = 0
    end if
  end subroutine

  subroutine get_l2_time_stats(fname, mjd, dt_error)
    implicit none
    character(len=*) :: fname
    real(dp)         :: mjd(2), dt_ideal, dt_error, srate
    integer(i4b)     :: n(7)
    type(hdf_file)   :: file
    real(dp),     dimension(:), allocatable :: time
    call open_hdf_file(fname, file, "r")
    call get_size_hdf(file, "time", n)
    allocate(time(n(1)))
    call read_hdf(file, "time", time)
    call read_hdf(file, "samprate", srate)
    dt_ideal = 1000/srate
    dt_error = maxval(abs((time(2:)-time(1:size(time)-1))*24*60*60*1000 - dt_ideal))/dt_ideal
    mjd = [ time(1), time(size(time)) ]
    call close_hdf_file(file)
    deallocate(time)
  end subroutine

end program
