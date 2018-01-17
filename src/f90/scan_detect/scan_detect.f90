! Scan_detect: Scan through level1-files in chronological order, and look for
! continuous regions with no motion in deck and el, and continuous motion in
! az, as well as scan mode = 3 and uniform 100 Hz time-steps. Output these
! in the form "$id $ces $mjd_from $mjd_to"
!
! How to parallelize: Run for a subrange for each. Use no overlap, but always
! discard a ces that was already in progress at the start of the range, and
! at the end, read in enough data to finish the last ces, even if this would
! exceed the end time. This will ensure that all ceses are detected.
!
! In the end, sort ceses by starting time, and assign a ces number to them.
! This will be done by an external program.
!
! Inputs:
!   1. Time-sorted file list
!   2. Directory files are relative to
!
! FOR PEOPLE WHO WANT TO MODIFY THIS:
!
!  1. The iterators are somewhat complicated, I admit, but you
!     should not need to worry about them at all.
!  2. All the useful work is done in the function 'validate', so start
!     looking there.
!  3. If you need more data, change get_data, and possibly chunk_type and
!     its copying routines.
!
!  As an example, to find skydips instead of CESes, you would simply swap the
!  roles of az and el in validate.

program scan_detect
  use healpix_types
  use quiet_fileutils
  use comap_pointing_mod
  use comap_detector_mod
  use quiet_shared_output_mod
  use comap_patch_detect_mod
  use comap_lx_mod
  use comap_defs
  implicit none

  ! A chunk will always have a constant sample rate
  ! icont = chunk number in section.
  ! isect = section number
  ! ifile = file section started in
  type chunk_type
     integer(i4b)  :: n, ifile
     integer(i4b), dimension(:),   allocatable :: mode
     real(dp),     dimension(:),   allocatable :: az, el, dk, time
  end type chunk_type

  type chunk_iterator
     integer(i4b)       :: file_index, pos, target_dt, n, iglob, icont, isect, ifile
     logical(lgt)       :: last_cont
     real(dp)           :: mjd
     type(lx_struct)    :: data
     !type(data_struct)  :: data_extra ! Only used to check readability
     type(chunk_type)   :: work
     integer(i4b),       dimension(:),   allocatable :: mods
     character(len=512), dimension(:),   allocatable :: files
     real(dp),           dimension(:,:), allocatable :: ranges
  end type chunk_iterator

  type subchunk_iterator
     type(chunk_iterator) :: iterator
     type(chunk_type)     :: chunk1, chunk2, chunk
     integer(i4b)         :: icont, chunk_size, i, state
     logical(lgt)         :: ok, tmp_ok
  end type

  type angle_helper
     real(dp)     :: angsum, sqrsum, last, angmin, angmax, lastmax, lastmin
     real(dp)     :: possum, negsum
     integer(i4b) :: n, nsust
  end type

  type validator
     type(angle_helper) :: az_helper, el_helper, dk_helper
     integer(i4b)       :: reason, samples, filebeg, filefrom, fileend
     integer(i4b)       :: object
     real(dp)           :: avg_hvec(3), avg_gvec(3), dev, ces_delay
     real(dp)           :: ces_range(2), mindur, avgdk, az_changed
     logical(lgt)       :: in_ces, done
     character(len=5)   :: type
     character(len=512) :: ces_ofilename, range_ofilename
    type(shared_ofile)  :: ces_ofile, range_ofile
  end type

  real(dp)             :: max_az_dev, max_el_dev, max_dk_dev
  real(dp)             :: rst_el_max_dev
  real(dp)             :: az_stall_lim, el_stall_lim
  real(dp)             :: ces_delay, ces_mindur, cas_mindur, rst_mindur
  real(dp)             :: rst_az_tol, rst_az_stall, rst_az_amp_lim
  type(subchunk_iterator) :: iterator
  type(chunk_type)     :: chunk
  type(validator)      :: vals(16)
  integer(i4b)         :: target_dt, chunk_size, glitchtol, i, j, k, samples, n, isun
  integer(i4b)         :: subchunk_size, az_stall_timeout, el_stall_timeout, ierr, nproc, myid
  integer(i4b)         :: foo
  logical(lgt)         :: object_detect, ok
  character(len=512)   :: listfile, parfile, l1prefix, output_dir, range_ofilename
  character(len=512)   :: pswitch_file
  real(dp),           dimension(:), allocatable :: pswitch
  character(len=512), dimension(:), allocatable :: filelist

  call getarg(1, parfile)

  ! These are implementation details that probably don't need to be
  ! configurable.
  chunk_size  = 5000 ! 10000 frames = 100 s
  target_dt   = 100    ! 10 ms
  glitchtol   = 300   ! ignore mode /= 3 glitches shorter than this

  call get_parameter(0, parfile, "CES_RESOLUTION", par_int=subchunk_size, &
   & desc="Block size for pass/fail test in ces detection, in samples. 100 is sensible.")
  call get_parameter(0, parfile, "LEVEL1_DIR", par_string=l1prefix, &
   & desc="Directory where the level1 file hierarchy can be found.")
  call get_parameter(0, parfile, "LEVEL1_FILELIST", par_string=listfile, &
   & desc="List of level1-files to consider, relative to LEVEL1_DIR. Fills' // &
   & 'the same role as L1_DATABASE, but does not rely on run/seg/obj classification from Chicago.")
  call get_parameter(0, parfile, "CES_MIN_DUR", par_dp=ces_mindur, &
   & desc="Minimum duration of a ces, in seconds. 100 is sensible.")
  call get_parameter(0, parfile, "CAS_MIN_DUR", par_dp=cas_mindur, &
   & desc="Minimum duration of a cas, in seconds. 10 is sensible.")
  call get_parameter(0, parfile, "RASTER_MIN_DUR", par_dp=rst_mindur, &
   & desc="Minimum duration of a raster, in seconds. 10 is sensible.")
  call get_parameter(0, parfile, "AZ_STALL_TIMEOUT", par_int=az_stall_timeout, &
   & desc="Number of consecutive frames below az_stall_lim before a stall is registred (30).")
  call get_parameter(0, parfile, "AZ_RASTER_TIMEOUT", par_dp=rst_az_stall, &
   & desc="Number of consecutive frames below AZ_RASTER_LIM before a stall is registred (300).")
  call get_parameter(0, parfile, "EL_STALL_TIMEOUT", par_int=el_stall_timeout, &
   & desc="Number of consecutive frames below el_stall_lim before a stall is registred (100).")
  call get_parameter(0, parfile, "AZ_CHANGE_LIM", par_dp=max_az_dev, &
  & desc="Rate of change of az in rads per sample needed for az to be considered changing. 1e-5 is sensible.")
  call get_parameter(0, parfile, "EL_CHANGE_LIM", par_dp=max_el_dev, &
  & desc="Max deviation in elevation in radians allowed. 1 arc min is sensible.")
  call get_parameter(0, parfile, "DK_CHANGE_LIM", par_dp=max_dk_dev, &
  & desc="Max deviation in deck in radians allowed. 10 arc min is sensible.")
  call get_parameter(0, parfile, "AZ_STALL_LIM", par_dp=az_stall_lim, &
  & desc="Rate of change of az in rads per sample needed for az to be considered stalled.")
  call get_parameter(0, parfile, "AZ_RASTER_LIM", par_dp=rst_az_tol, &
  & desc="Rate of change of az in rads per sample needed for az to be considered stalled for raster scans (1e-5).")
  call get_parameter(0, parfile, "EL_STALL_LIM", par_dp=el_stall_lim, &
  & desc="Rate of change of el in rads per sample needed for el to be considered stalled.")
  call get_parameter(0, parfile, "OUTPUT_DIR",    par_string=output_dir, &
  & desc="Directory where the output files will be placed.")
  call get_parameter(0, parfile, "CES_DELAY", par_dp=ces_delay, &
  & desc="Number of seconds to cut out at the start of the CES.")
  call get_parameter(0, parfile, "CES_OBJECT", par_lgt=object_detect, &
  & desc="Whether to run the object detection or not. This is costly.")
  call get_parameter(0, parfile, "PHASE_SWITCH_FILE", par_string=pswitch_file, &
  & desc="A file with sorted mjds for phase switch turn-on events. Used with CES_DELAY.")
  call get_parameter(0, parfile, "AZ_RASTER_AMP_LIM", par_dp=rst_az_amp_lim, &
  & desc="Max allowed raster scan amplitude (0.035).")
  call get_parameter(0, parfile, "EL_RASTER_LIM", par_dp=rst_el_max_dev, &
  & desc="Max deviation in elevation in radians allowed for raster scans (0.035).")

  ces_delay = ces_delay/24/60/60
  call initialize_detector_mod(parfile)

  call initialize_comap_pointing_mod(parfile)
  call initialize_patch_detect_mod(parfile)
!  isun = lookup_patch("sun", patches)

  call read_alloc_filelist(listfile, filelist)
  !call read_pswitch(pswitch_file, pswitch)

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid,  ierr)
  call dset(id=myid)

  ! I will have responsibility for file [n*myid/nproc+1,n*(myid+1)/nproc]
  ! Start one file earlier, to allow the mode fixing to work properly
  ! even at the start. My responsibility starts at the first section
  ! which started inside my region, and ends at the end of the last chunk
  ! that started inside my region.

  ! Raster scans... Short steps in el every few seconds. El step duration is
  ! about 50-100 samples long, reaching amplitudes of 0.1 arc minutes per sample.
  ! During this, regular az scanning is going on, with a ~100 sample long pause
  ! when the elevation changes. Mode is 3.
  !
  ! Between these scans, there are 10-30 second long pauses in az-scanning, with
  ! large, long-lasting (3s+) elevation changes, and mode briefly switches from 3.

!open(40,file="foo.txt")
  n = 0
  n=n+1; call init_validator(vals(n),"ces",output_dir,size(filelist),ces_mindur,ces_delay)
  n=n+1; call init_validator(vals(n),"cas",output_dir,size(filelist),cas_mindur,ces_delay)
  !n=n+1; call init_validator(vals(n),"sun",output_dir,size(filelist),0d0,       0d0)
  n=n+1; call init_validator(vals(n),"rst",output_dir,size(filelist),rst_mindur,ces_delay)
  n=n+1; call init_validator(vals(n),"liss",output_dir,size(filelist),ces_mindur,ces_delay)
  call init_subchunk_iterator(iterator, filelist(vals(1)%filefrom:), chunk_size, subchunk_size, target_dt)
  write(*,*) "so far so good"
  ok = .true.
  do while(ok .and. any(.not. vals(1:n)%done))
     ok = next_subchunk(iterator, chunk)
     write(*,*) "done?: ", vals(4)%done, "ok?:", ok, "chunk%n", chunk%n
     do i = 1, n
        call process_chunk(chunk, vals(i))
        write(*,*) "reason", vals(i)%reason
        write(*,*) "ok?", ok
        write(*,*) "done?", vals(4)%done
     end do
  end do
  do i = 1, n
     call free_validator(vals(i))
  end do
  write(*,*) "Better"
  deallocate(filelist)
  call free_chunk(chunk)
  call free_subchunk_iterator(iterator)
  call mpi_finalize(ierr)
contains

  subroutine process_chunk(chunk, val)
    implicit none
    type(chunk_type)   :: chunk
    type(validator)    :: val
    integer(i4b)       :: ifile, i, j, k, m, n
    logical(lgt)       :: valid, split
    real(dp)           :: orange(2), avgel, avgaz, avgtheta, avgphi, avgpsi
    real(dp)           :: hvec(3), gvec(3), az, el, dk, phi, theta, psi, dist
    character(len=512) :: line, objname
    integer(i4b), dimension(:), allocatable :: matches
    
    if(val%done) return
    ifile = chunk%ifile+val%filefrom  ! There may be something wrong with the file indices!?!
    write(*,*) "ifile", ifile, "chunkifile", chunk%ifile
         
    ! Loop until we enter our area of responsibility
    if(ifile < val%filebeg) return
    write(*,*) "inside pr.ch", chunk%n
         
    if(chunk%n == 0) then
       valid = .false.
       split = .false.
    else
       select case(val%type)
          case("ces");  call validate_ces(chunk, val, valid, split) 
          case("cas"); call validate_cas(chunk, val, valid, split)
!          case("sun");  call validate_sun(chunk, val, valid, split)
          case("rst");  call validate_rst(chunk, val, valid, split)
          case("liss"); call validate_liss(chunk, val, valid, split)
          case default
             write(stderr,*) "Invalid type to detect: " // trim(val%type)
             stop
       end select
    end if
!if(val%type=="rst") then
!do i = 1, chunk%n
!write(40,'(f13.7,3f15.10,i4,i7)') chunk%time(i), chunk%az(i), chunk%el(i), chunk%dk(i), val%reason, iterator%icont
!end do
!end if
    write(*,*) "valid", valid
    write(*,*) "split", split
    write(*,*) "done?:", val%done
    write(*,*) "in ces?:", val%in_ces
    ! Has a ces ended?
    if(val%in_ces .and. (split .or. .not. valid)) then
       ! Should we restrict the range?
       orange = val%ces_range
       write(*,*) "ces_range", orange
       write(*,*) "ces_delay", val%ces_delay
     !  if(val%type == "ces") call adjust_range(orange, val%ces_delay)  ! removed because we don't have pswitch (whatever that is
       if(orange(2)-orange(1) >= val%mindur/24/60/60) then
          call vec2ang(val%avg_hvec, avgel,    avgaz)
          call vec2ang(val%avg_gvec, avgtheta, avgphi)
          call swap_coordinate_convention(avgaz, avgel,    val%avgdk, COORD_HOR)
          call swap_coordinate_convention(avgphi,avgtheta, avgpsi,    COORD_GAL)
          objname = "unknown"
          if(val%object > 0) objname = patches(val%object)%name
          if(ifile >= val%filebeg) then
             write(line,fmt="(2f14.7,5f9.3,a)") orange, avgaz*RAD2DEG, &
              & avgel*RAD2DEG, val%avgdk*RAD2DEG, avgphi*RAD2DEG, &
              & avgtheta*RAD2DEG, " " // trim(objname)
             write(*,fmt="(i4,a4,a)") myid, " " // val%type, trim(line)
             call write_shared_ofile(val%ces_ofile, trim(line))
             call flush_shared_ofile(val%ces_ofile)
             call get_overlaps(iterator%iterator%ranges, orange, matches)
             do i = 1, size(matches)
                write(line,fmt="(2f14.7,a)") iterator%iterator%ranges(matches(i),:), &
                 & " " // trim(iterator%iterator%files(matches(i)))
                call write_shared_ofile(val%range_ofile, trim(line))
                call flush_shared_ofile(val%range_ofile)
             end do
             deallocate(matches)
          end if
       end if
       val%avg_hvec = 0
       val%avg_gvec = 0
       val%dev      = 0
       val%in_ces   = .false.
       val%samples  = 0
       val%object   = 0
    end if
    ! After making sure we output our last ces, exit if we are
    ! outside our area.
    if(.not. valid .and. ifile > val%fileend) then !was originally "ifile >= val%fileend) then", but we ended after first ces always. This seems to work'
       val%done = .true.
       write(*,*) "fileend", val%fileend, "ifile", ifile
       write(*,*) "now I'm done"
       return
    end if

    ! If we are inside a ces, take some stats
    if(valid) then
       if(.not. val%in_ces) val%ces_range(1) = chunk%time(1)
       do i = 1, chunk%n, 7
          az = chunk%az(i); el = chunk%el(i); dk = chunk%dk(i)
          call swap_coordinate_convention(az, el, dk, COORD_HOR)
          call coord_convert(COORD_HOR, az, el, dk, COORD_GAL, phi, theta, psi, chunk%time(i))
          call ang2vec(el,    az,  hvec)
          call ang2vec(theta, phi, gvec)
          val%avg_hvec = val%avg_hvec + hvec
          val%avg_gvec = val%avg_gvec + gvec
          val%avgdk    = dk
          call angdist2(val%avg_gvec,gvec,dist)
          val%dev  = val%dev + dist**2
          val%samples  = val%samples+1
       end do
       val%ces_range(2) = chunk%time(chunk%n)
!       if(object_detect) call find_object(chunk, patches, val%object)
    end if
    val%in_ces   = valid
  end subroutine process_chunk

  subroutine init_validator(val, type, odir, filelen, mindur, delay)
    implicit none
    type(validator)  :: val
    character(len=*) :: type, odir
    integer(i4b)     :: filelen
    real(dp)         :: delay, mindur
    call reset_angle_helper(val%az_helper)
    call reset_angle_helper(val%el_helper)
    call reset_angle_helper(val%dk_helper)
    val%mindur   = mindur
    val%avg_hvec = 0
    val%avg_gvec = 0
    val%dev      = 0
    val%in_ces   = .false.
    val%samples  = 0
    val%object   = 0
    val%az_changed = 0
    val%ces_delay= delay
    val%type     = type
    val%filebeg  = filelen*myid/nproc+1
    val%filefrom = max(1, val%filebeg-1)
    val%fileend  = filelen*(myid+1)/nproc+1 ! The first index not part of what we want
!    write(*,*) "filebeg, fileend", val%filebeg, val%fileend
    val%done     = val%filebeg == val%fileend
!    write(*,*) val%done
    val%ces_ofilename   = trim(odir) // "/" // trim(type) // "_list.txt"
    val%range_ofilename = trim(odir) // "/" // trim(type) // "_fileranges.txt"
    call mkdirs(val%ces_ofilename,   .true.)
    call mkdirs(val%range_ofilename, .true.)
    call open_shared_ofile(val%ces_ofile,   val%ces_ofilename,   MPI_COMM_WORLD)
    call open_shared_ofile(val%range_ofile, val%range_ofilename, MPI_COMM_WORLD)
  end subroutine

  subroutine free_validator(val)
    implicit none
    type(validator) :: val
    call close_shared_ofile(val%ces_ofile)
    call close_shared_ofile(val%range_ofile)
  end subroutine

  ! The chunk_iterator iterates over the level1-data as a continuous stream,
  ! skipping data with the wrong sampling rate. The chunks will have size
  ! chunk_size unless a sample rate change happens, which will make them stop
  ! early. All data within a chunk is continuous, and the icont counter increases
  ! by one for every chunk which continuously follows the previous one, and
  ! is reset when a jump happens.
  subroutine init_chunk_iterator(iterator, filelist, chunk_size, target_dt)
    implicit none
    type(chunk_iterator) :: iterator
    character(len=*)     :: filelist(:)
    integer(i4b)         :: chunk_size, target_dt, i, unit, status
    call free_chunk_iterator(iterator)
    allocate(iterator%files(size(filelist)))
    allocate(iterator%ranges(size(filelist),2))
    iterator%files = filelist
    call init_chunk(iterator%work, chunk_size)
    iterator%target_dt = target_dt
    iterator%iglob     = 0
    iterator%icont     = 0
    iterator%isect     = 0
    iterator%ifile     = 0
!    iterator%ifile     = 1
    iterator%ranges    = 0
    iterator%last_cont = .true.
    allocate(iterator%mods(get_num_dets()))
    do i = 1, size(iterator%mods); iterator%mods(i) = i-1; end do
    iterator%file_index = 1
    ! Read in the data, and scan until we reach the correct time
    call get_data(l1prefix, iterator%files, iterator%file_index, iterator%mods, iterator%data)
    iterator%n = size(iterator%data%scanmode_l1)
    write(*,*) "filesize (n_t):", size(iterator%data%time_point)
    i = 1
!    write(*,*) [1,iterator%n]
!    write(*,*) iterator%data%time([1,iterator%n])
    iterator%ranges(iterator%file_index,:) = iterator%data%time_point([1,iterator%n])
    iterator%pos = i-1
    ! Assume that the previous step had the correct length
    iterator%mjd = iterator%data%time_point(i) - iterator%target_dt/24d0/60/60/1000
  end subroutine

  function next_chunk(iterator, chunk) result(ok)
    implicit none
    type(chunk_iterator) :: iterator
    type(chunk_type)     :: chunk
    logical(lgt)         :: ok
    integer(i4b)         :: i, j, k, m, n, unit, status
    real(dp)             :: dt
    ok   = iterator%file_index <= size(iterator%files)
    !write(*,*) "ok in next chunk:", ok, iterator%file_index, size(iterator%files)
    if(.not. ok) return
    n    = iterator%work%n
    unit = getlun()
    write(*,*) "icont in next_chunk", iterator%icont, "last_cont", iterator%last_cont
    if(iterator%last_cont) then
       iterator%icont = 0
       iterator%isect = iterator%isect + 1
    end if
    iterator%iglob = iterator%iglob + 1
    iterator%icont = iterator%icont + 1
    iterator%last_cont = .false.
    ! Copy into work as long as the frame rate is correct, and
    ! get more data if necessary
    do i = 1, n
       iterator%pos = iterator%pos + 1
       if(iterator%pos > iterator%n) then
          ! are we outside end of file
          iterator%file_index = iterator%file_index + 1
          call get_data(l1prefix, iterator%files, iterator%file_index, iterator%mods, iterator%data)
          if(iterator%file_index > size(iterator%files)) exit
          iterator%n   = size(iterator%data%scanmode_l1)
          iterator%ranges(iterator%file_index,:) = iterator%data%time_point([1,iterator%n])
          iterator%pos = 1
       end if
       dt = (iterator%data%time_point(iterator%pos) - iterator%mjd)*24*60*60*1000
       iterator%mjd = iterator%data%time_point(iterator%pos)
       !write(*,*) "dt - comparison", nint(dt), iterator%target_dt
       if(abs(nint(dt) - iterator%target_dt) > 3) then
          iterator%last_cont = .true.
          exit
       end if
       if(i == 1 .and. iterator%icont == 1) iterator%ifile = iterator%file_index
       iterator%work%ifile   = iterator%ifile
       iterator%work%mode(i) = iterator%data%scanmode_l1(iterator%pos)
       iterator%work%time(i) = iterator%data%time_point(iterator%pos)
       iterator%work%az(i)   = iterator%data%point_tel(1, iterator%pos) * pi / 180.0
       iterator%work%el(i)   = iterator%data%point_tel(2, iterator%pos) * pi / 180.0
       iterator%work%dk(i)   = iterator%data%point_tel(3, iterator%pos) * pi / 180.0
    end do
    ! Now copy over what we got.
    call copy_chunk(iterator%work, chunk, i-1)
  end function

  subroutine free_chunk_iterator(iterator)
    implicit none
    type(chunk_iterator) :: iterator
    call free_lx_struct(iterator%data)
    !call deallocate_data_struct(iterator%data_extra)
    call free_chunk(iterator%work)
    if(allocated(iterator%files)) deallocate(iterator%files)
    if(allocated(iterator%mods))  deallocate(iterator%mods)
  end subroutine

  ! Subchunk iterators work basically as a wrapper around chunk iterators. They
  ! are needed to fix problems with the data stream that cannot be handled by
  ! simply looking at a single chunk at a time. Right now, all this code
  ! is there basically to perform fix_mode. If fix_mode had not been needed,
  ! of of this could be skipped.
  subroutine init_subchunk_iterator(iterator, filelist, chunk_size, subchunk_size, target_dt)
    implicit none
    type(subchunk_iterator) :: iterator
    character(len=*)        :: filelist(:)
    integer(i4b)            :: chunk_size, subchunk_size, target_dt
    real(dp)                :: mjd_start
    logical(lgt)            :: foo
    call free_subchunk_iterator(iterator)
    call init_chunk_iterator(iterator%iterator, filelist, chunk_size, target_dt)
    iterator%chunk_size = subchunk_size
    write(*,*) "icont in init_subchunk_iterator", iterator%iterator%icont, "last_cont", iterator%iterator%last_cont
    iterator%ok         = next_chunk(iterator%iterator, iterator%chunk1)
    write(*,*) "icont after first next_chunk", iterator%iterator%icont
    iterator%state      = 0
  end subroutine

  subroutine free_subchunk_iterator(iterator)
    implicit none
    type(subchunk_iterator) :: iterator
    call free_chunk_iterator(iterator%iterator)
    call free_chunk(iterator%chunk1)
    call free_chunk(iterator%chunk2)
    call free_chunk(iterator%chunk)
  end subroutine

  function next_subchunk(iterator, chunk) result(ok)
    implicit none
    type(subchunk_iterator) :: iterator
    type(chunk_type)        :: chunk
    logical(lgt)            :: ok
    integer(i4b)            :: i, j, k, m, n

    chunk%n = 0
    ok = iterator%ok
    if(.not. ok) return
    write(*,*) "iterator%state", iterator%state
    if(iterator%state == 1) goto 1
    if(iterator%state == 2) goto 2
    write(*,*) "Got this far"
    do while(iterator%ok)
       iterator%i     = 1
       iterator%icont = 0
       do ! moving window
          iterator%tmp_ok = next_chunk(iterator%iterator, iterator%chunk2)
          !write(*,*) "icont", iterator%iterator%icont
          if(.not. iterator%tmp_ok .or. iterator%iterator%icont == 1) exit ! done with this section
          write(*,*) "But never this far"
          call append_chunk(iterator%chunk1, iterator%chunk2, iterator%chunk)
          !write(*,*) "mode", iterator%chunk%mode
          call fix_mode(iterator%chunk%mode, glitchtol)
1         iterator%state = 1
          m = iterator%chunk1%n+iterator%chunk2%n/2
          do while(iterator%i < m)
             n = min(iterator%chunk_size, iterator%chunk%n-iterator%i+1)
             call copy_chunk(iterator%chunk, chunk, from=iterator%i, len=n)
             iterator%i     = iterator%i     + n
             iterator%icont = iterator%icont + 1
             return
          end do
          iterator%i = iterator%i - iterator%chunk1%n
          call copy_chunk(iterator%chunk2, iterator%chunk1)
       end do
       ! we got here either because we ran out of chunks, or because
       ! we reached the end of the section. in both cases, there might
       ! be something left to process at the end of chunk1
       call copy_chunk(iterator%chunk1, iterator%chunk)
       call fix_mode(iterator%chunk%mode, glitchtol)
2      iterator%state = 2
       do while(iterator%i < iterator%chunk%n)
          n = min(iterator%chunk_size, iterator%chunk%n-iterator%i+1)
          call copy_chunk(iterator%chunk, chunk, from=iterator%i, len=n)
          iterator%i = iterator%i         + n
          iterator%icont = iterator%icont + 1
          return
       end do
       call copy_chunk(iterator%chunk2, iterator%chunk1)
       iterator%ok    = iterator%tmp_ok
       iterator%state = 0
    end do
    ok = .false.
  end function

  subroutine free_chunk(chunk)
    implicit none
    type(chunk_type) :: chunk
    if(allocated(chunk%mode)) deallocate(chunk%mode)
    if(allocated(chunk%time)) deallocate(chunk%time)
    if(allocated(chunk%az))   deallocate(chunk%az)
    if(allocated(chunk%el))   deallocate(chunk%el)
    if(allocated(chunk%dk))   deallocate(chunk%dk)
  end subroutine

  subroutine init_chunk(chunk, n)
    implicit none
    type(chunk_type) :: chunk
    integer(i4b)     :: n
    call free_chunk(chunk)
    allocate(chunk%mode(n))
    allocate(chunk%time(n))
    allocate(chunk%az(n))
    allocate(chunk%el(n))
    allocate(chunk%dk(n))
    chunk%n = n
  end subroutine

  subroutine copy_chunk(chunk1, chunk2, len, from)
    implicit none
    type(chunk_type)       :: chunk1, chunk2
    integer(i4b), optional :: len, from
    integer(i4b)           :: n, f
    n = chunk1%n; if(present(len))  n = len
    f = 1;        if(present(from)) f = from
    call init_chunk(chunk2, n)
    chunk2%ifile= chunk1%ifile
    chunk2%mode = chunk1%mode(f:f+n-1)
    chunk2%time = chunk1%time(f:f+n-1)
    chunk2%az   = chunk1%az  (f:f+n-1)
    chunk2%el   = chunk1%el  (f:f+n-1)
    chunk2%dk   = chunk1%dk  (f:f+n-1)
  end subroutine

  subroutine append_chunk(chunk1, chunk2, res)
    implicit none
    type(chunk_type)       :: chunk1, chunk2, res
    call init_chunk(res, chunk1%n + chunk2%n)
    res%mode(1:chunk1%n)  = chunk1%mode
    res%mode(chunk1%n+1:) = chunk2%mode
    res%time(1:chunk1%n)  = chunk1%time
    res%time(chunk1%n+1:) = chunk2%time
    res%az(1:chunk1%n)  = chunk1%az
    res%az(chunk1%n+1:) = chunk2%az
    res%el(1:chunk1%n)  = chunk1%el
    res%el(chunk1%n+1:) = chunk2%el
    res%dk(1:chunk1%n)  = chunk1%dk
    res%dk(chunk1%n+1:) = chunk2%dk
    res%ifile           = chunk1%ifile
  end subroutine

  subroutine read_alloc_filelist(listfile, list)
    implicit none
    character(len=*)   :: listfile
    character(len=512) :: line
    character(len=*),  allocatable, dimension(:)  :: list
    integer(i4b)     :: i, j, n, unit
    if(allocated(list)) deallocate(list)
    unit = getlun()
    open(unit,file=listfile,status="old",action="read")
    n = 0
    do
       read(unit,*,end=1)
       n = n+1
    end do
1   continue
    rewind(unit)
    allocate(list(n))
    do i = 1, n
       read(unit,fmt="(a)") list(i)
    end do
    close(unit)
  end subroutine

  subroutine fix_mode(mode, mindur)
    implicit none
    integer(i4b) :: mode(:), mindur, i, j, k, n
    n = size(mode)
    i = 1
    do while(i < n)
       ! Find first deviant sample
       do
          if(i >= n) exit
          if(mode(i+1) /= mode(i)) exit
          i = i+1
       end do
       ! Check if deviation is shorter than mindur
       j = 1
       do
          if(i+j > n) exit
          if(j > mindur) exit
          if(mode(i+j) == mode(i)) exit
          j = j+1
       end do
       if(i+j > n) exit
       ! If to short, overwrite deviant samples
       if(j < mindur) mode(i+1:i+j-1) = mode(i)
       i = i+j
    end do
  end subroutine

  subroutine validate_liss(chunk, val, ok, split)
    implicit none
    type(chunk_type) :: chunk
    type(validator)  :: val
    logical(lgt)     :: ok, split
    real(dp), dimension(:), allocatable :: work, subtracted_data
    integer(i4b)     :: nstall, nstall_tmp, iii
    split = iterator%icont == 1
    ok = .false.
    write(*,*) "inside liss, n: ", chunk%n, size(chunk%el), pi
    !val%reason = 1
    !if(any(chunk%mode /= 3)) return
    val%reason = 2
    write(*,*) "sizechunk",size(chunk%el)
    !do iii= 1,size(chunk%el)
       !write(*,*) chunk%el(iii)
    !   chunk%el(iii) = chunk%el(iii) * pi / 180.0
    !end do
    if(any(chunk%el < -pi/2) .or. any(chunk%el > pi/2))    return
    val%reason = 3
    if(angle_const(chunk%el, val%el_helper, max_el_dev)) return
    !if(stalled(chunk%az, az_stall_lim, az_stall_timeout))             return
    !val%reason = 4
    !if(stalled(chunk%el, el_stall_lim, el_stall_timeout))         return
    val%reason = 5
    allocate(subtracted_data(chunk%n))
    call remove_sinusoid(chunk%el, chunk%time, chunk%n, subtracted_data)
    if(.not. angle_const(subtracted_data, val%el_helper, max_el_dev)) return
    call remove_sinusoid(chunk%az, chunk%time, chunk%n, subtracted_data)
    if(.not. angle_const(subtracted_data, val%az_helper, max_az_dev)) return
    val%reason = 0
    ok = .true. !(should be true)
  end subroutine

  subroutine validate_ces(chunk, val, ok, split)
    implicit none
    type(chunk_type) :: chunk
    type(validator)  :: val
    logical(lgt)     :: ok, split
    real(dp), dimension(:), allocatable :: work
    integer(i4b)     :: nstall, nstall_tmp
    split = iterator%icont == 1
    ok = .false.
    !val%reason = 1
    !if(any(chunk%mode /= 3)) return
    val%reason = 2
    if(any(chunk%el < -pi/2) .or. any(chunk%el > pi/2))    return
    val%reason = 3
    if(stalled(chunk%az, az_stall_lim, az_stall_timeout))             return
    val%reason = 4
    if(.not. angle_const(chunk%el, val%el_helper, max_el_dev)) return
    val%reason = 5
    if(.not. angle_const(chunk%dk, val%dk_helper, max_dk_dev)) return
    val%reason = 0
    ok = .true.
  end subroutine

  subroutine validate_cas(chunk, val, ok, split)
    implicit none
    type(chunk_type) :: chunk
    type(validator)  :: val
    logical(lgt)     :: ok, split
    real(dp), dimension(:), allocatable :: work
    integer(i4b)     :: nstall, nstall_tmp
    split = iterator%icont == 1
    ok = .false.
    val%reason = 1
    if(any(chunk%mode /= 3)) return
    val%reason = 2
    if(any(chunk%el < -pi/2) .or. any(chunk%el > pi/2))    return
    val%reason = 3
    if(stalled(chunk%el, el_stall_lim, el_stall_timeout))             return
    val%reason = 4
    if(.not. angle_const(chunk%az, val%az_helper, max_az_dev)) return
    val%reason = 5
    if(.not. angle_const(chunk%dk, val%dk_helper, max_dk_dev)) return
    val%reason = 0
    ok = .true.
  end subroutine

!!$  subroutine validate_sun(chunk, val, ok, split)
!!$    implicit none
!!$    type(chunk_type) :: chunk
!!$    type(validator)  :: val
!!$    logical(lgt)     :: ok, split
!!$    integer(i4b)     :: obj
!!$    split = .false.
!!$    obj = 0
!!$    !call find_object(chunk, patches(isun:isun), obj)
!!$    ok = obj /= 0
!!$  end subroutine

  subroutine validate_rst(chunk, val, ok, split)
    implicit none
    type(chunk_type) :: chunk
    type(validator)  :: val
    logical(lgt)     :: ok, split
    real(dp), dimension(:), allocatable :: work
    integer(i4b)     :: nstall, nstall_tmp, n, pslope, nslope
    split = iterator%icont == 1
    val%reason = 1
    if(any(chunk%mode /= 3))                                          goto 1
    val%reason = 2
    if(any(chunk%el < -pi/2) .or. any(chunk%el > pi/2))               goto 1
    val%reason = 3
    call update_angle_helper(val%az_helper, chunk%az)
    call update_angle_helper(val%el_helper, chunk%el)
    call update_angle_helper(val%dk_helper, chunk%dk)
    if(maxangdev(val%el_helper) > rst_el_max_dev)                     goto 1
    val%reason = 4
    if(maxangdev(val%dk_helper) > max_dk_dev)                         goto 1
    val%reason = 5
    ! Ok, this was very similar to ces so far. But azimuth is a bit different, as
    ! we have relatively long periods with no az change, which are nevertheless ok.
    ! This is a weakness of the chunk method. A look-back method would be better.
    call rst_az_changed(chunk, val%az_changed)
    if((chunk%time(size(chunk%time)) - val%az_changed)*24*3600*100 > rst_az_stall) return
    val%reason = 6
    if(angdev(val%az_helper) > rst_az_amp_lim)                        goto 1
    val%reason = 0
    ok = .true.
    return
    1 ok = .false.
    call reset_angle_helper(val%az_helper)
    call reset_angle_helper(val%el_helper)
    call reset_angle_helper(val%dk_helper)
  end subroutine

  function angdev(helper)
    implicit none
    type(angle_helper) :: helper
    real(dp) :: angdev
    angdev = sqrt(helper%sqrsum/helper%n-(helper%angsum/helper%n)**2)
  end function

  function maxangdev(helper)
    implicit none
    type(angle_helper) :: helper
    real(dp) :: maxangdev
    maxangdev = maxval(abs(helper%angsum/helper%n-[helper%lastmax,helper%lastmin]))
  end function

  function stalled(arr, lim, timeout)
    implicit none
    real(dp)         :: arr(:), lim
    logical(lgt)     :: stalled
    real(dp), dimension(:), allocatable :: work, safe
    integer(i4b)     :: nstall, nstall_tmp, n, timeout, i
    allocate(safe(size(arr)),work(size(arr)-1))
    safe = arr
    call make_angles_safe(safe, 2*pi)
    work   = abs(safe(2:size(arr)) - safe(1:size(arr)-1))
    nstall = 0
    nstall_tmp = 0
    do i = 1, size(work)
       if(work(i) < lim) then
          nstall_tmp = nstall_tmp +1
       else
          nstall = max(nstall,nstall_tmp)
          nstall_tmp = 0
       end if
    end do
    nstall = max(nstall,nstall_tmp)
    deallocate(work, safe)
    stalled = nstall > timeout
  end function

  subroutine rst_az_changed(chunk, last_changed)
    implicit none
    type(chunk_type), intent(in)    :: chunk
    real(dp),         intent(inout) :: last_changed
    real(dp)                        :: safe(size(chunk%az))
    integer(i4b)                    :: i
    safe = chunk%az
    call make_angles_safe(safe, 2*pi)
    do i = size(safe), 2, -1
       if(abs(safe(i)-safe(i-1)) > rst_az_tol) exit
    end do
    if(i > 1) last_changed = chunk%time(i)
  end subroutine

  ! Try to read file ind, but possibly advance a bit if ind cannot be found.
  subroutine get_data(prefix, filelist, ind, mods, point)
    implicit none
    character(len=*)   :: prefix, filelist(:)
    character(len=512) :: file
    integer(i4b)       :: ind, mods(:), status, unit
    type(lx_struct)    :: point

    unit = getlun()
    do
       if(ind > size(filelist)) exit
       file = trim(prefix) // "/" // trim(filelist(ind))
       call free_lx_struct(point)
       status = 0
       call read_l1_file(file, data=point, only_point=.true.)
       if(status /= 0) write(stderr,*) "Error reading '" // trim(file) // "'!"
       if(status == 0) exit
       ind = ind + 1
    end do
  end subroutine

  subroutine get_overlaps(ranges, mrange, indices)
    implicit none
    real(dp)     :: ranges(:,:), mrange(:)
    integer(i4b), dimension(:), allocatable :: indices, work
    integer(i4b) :: i, n, m
    if(allocated(indices)) deallocate(indices)
    n = size(ranges,1)
    m = 0
    allocate(work(n))
    do i = 1, n
       if(ranges(i,2) > mrange(1) .and. ranges(i,1) <= mrange(2)) then
          m = m+1
          work(m) = i
       end if
    end do
    allocate(indices(m))
    indices = work(1:m)
    deallocate(work)
  end subroutine

  function angle_const(angs, helper, dev) result(ok)
    implicit none
    real(dp)           :: angs(:), dev
    type(angle_helper) :: helper
    logical(lgt)       :: ok
    call update_angle_helper(helper, angs)
    if(any(abs(helper%angsum/helper%n-[helper%lastmax,helper%lastmin]) > dev)) then
       call reset_angle_helper(helper)
       ok = .false.
    else
       ok = .true.
    end if
  end function

  ! This routine is heavy, because it needs to do a full coordinate
  ! conversion, but identifies the object being
  ! scanned with high accuracy. Since a chunk is pretty short,
  ! we only use 1 sample for ephemeris. In a way it would make
  ! more sense to do this object detection in l3gen, since
  ! pointing conversions depend on the mount model, which we might
  ! want to change often. But we only need approximate pointing here.
  subroutine find_object(chunk, patches, obj)
    implicit none
    type(chunk_type) :: chunk
    type(patch_info) :: patches(:)
    integer(i4b)     :: obj, i, j, k, m, n, mod, nmod, red, tobj
    real(dp)         :: patch_pos(3,1,size(patches)), fprad
    real(dp)         :: op(3), np(3), mat(3,3), red_blur, samprate, maxbeam, beam_rad
    logical(lgt)     :: mask(size(patches))
    real(dp),     dimension(:,:,:), allocatable :: mats
    real(dp),     dimension(:,:),   allocatable :: point
    integer(i4b), dimension(:),     allocatable :: mobjs, objs, inds
    call get_patch_pos_multi(patches, [chunk%time(1)], COORD_GAL, patch_pos)
    nmod     = get_num_dets()
    samprate = 100
    red      = 4
    n        = (chunk%n-1)/red+1
    red_blur = 2.0*red/samprate/2*pi/180 ! 2 degree per second scanning
    fprad    = get_focalplane_radius()
    maxbeam  = get_max_fwhm()*fwhm2sigma*5;
    allocate(mobjs(0:nmod-1), objs(n), point(3,n), mats(3,3,n))
    j = 0
    do i = 1, chunk%n, red
      j = j+1
      op = [ chunk%az(i), chunk%el(i), chunk%dk(i) ]
      call swap_coordinate_convention(op(1), op(2), op(3), COORD_TELE)
      call coord_convert(COORD_TELE, op(1), op(2), op(3), &
       & COORD_GAL, np(1), np(2), np(3), mjd=chunk%time(i), euler=mats(:,:,j))
      point(:,j) = np
    end do
    ! Find the distance to each object, so that we can exclude some of them
    call get_patch_hits(patches, patch_pos, point(1,:), point(2,:), &
     & fprad + red_blur + maxbeam, objs)
    ! Which patches are relevant?
    mask = .false.
    do i = 1, n
       if(objs(i) == 0) cycle
       mask(objs(i)) = .true.
    end do
    m = count(mask)
    allocate(inds(0:m))
    inds = 0
    j    = 0
    do i = 1, size(mask)
       if(.not. mask(i)) cycle
       j = j+1
       inds(j) = i
    end do
    ! Ok, from now on, only consider those in inds.
    tobj = 0
    if(m > 0) then
       do mod = 0, nmod-1
          do j = 1, n
             call rot2angles(matmul(mats(:,:,j), rot_detector2boresight(mod)), &
              & point(1,j),point(2,j),point(3,j))
          end do
          beam_rad = get_detector_fwhm(mod)*fwhm2sigma*5
          call get_patch_hits(patches(inds(1:)), patch_pos(:,:,inds(1:)), point(1,:), &
           & point(2,:), beam_rad+red_blur, objs)
          mobjs(mod) = inds(reduce_hits(patches, objs))
       end do
       tobj   = reduce_hits(patches, mobjs)
    end if
    obj = reduce_hits(patches, [obj,tobj])
    deallocate(point,mobjs,objs,mats,inds)
  end subroutine

!!$  subroutine read_pswitch(file, data)
!!$    implicit none
!!$    character(len=*), intent(in)                :: file
!!$    real(dp),         dimension(:), allocatable :: data
!!$    integer(i4b)                                :: i, j, k, n, unit
!!$    real(dp)                                    :: val
!!$    unit = getlun()
!!$    open(unit,file=file,status="old",action="read")
!!$    n = 0
!!$    do
!!$       read(unit,*,end=1) val
!!$       n = n+1
!!$    end do
!!$1   rewind(unit)
!!$    allocate(data(n))
!!$    do i = 1, n
!!$       read(unit,*) data(i)
!!$    end do
!!$    close(unit)
!!$    call quicksort_real(data)
!!$  end subroutine

  subroutine adjust_range(range, delay)
    implicit none
    real(dp), intent(inout) :: range(2)
    real(dp), intent(in)    :: delay
    integer(i4b)            :: i
    ! Find the last phase switch event before range(2)
    !write(*,*) "range",range(2)
    !write(*,*) "pswitch",pswitch
    i = locate(pswitch, range(2))
    if(i <= 0 .or. i > size(pswitch)) return
    range(1) = max(range(1), pswitch(i) + delay)
  end subroutine

  subroutine reset_angle_helper(helper)
    implicit none
    type(angle_helper) :: helper
    helper%n = 0
    helper%angsum = 0
    helper%sqrsum = 0
    helper%angmin = infinity
    helper%angmax = -infinity
    helper%last = 0
    helper%possum = 0
    helper%negsum = 0
  end subroutine

  subroutine update_angle_helper(helper, arr)
    implicit none
    type(angle_helper) :: helper
    integer(i4b)       :: i
    real(dp)           :: arr(:), safe(size(arr)), prev
    safe = arr
    call make_angles_safe(safe, 2*pi)
    safe = safe + nearest_ang(helper%last,safe(1))-safe(1)
    helper%n      = helper%n      + size(safe)
    helper%angsum = helper%angsum + sum(safe)
    helper%sqrsum = helper%sqrsum + sum(safe**2)
    helper%angmax = max(helper%angmax, maxval(safe))
    helper%angmin = min(helper%angmin, minval(safe))
    helper%lastmax= maxval(safe)
    helper%lastmin= minval(safe)
    prev = helper%last
    do i = 1, size(safe)
       helper%possum = helper%possum + max(0d0,safe(i)-prev)
       helper%negsum = helper%negsum + max(0d0,prev-safe(i))
       prev = safe(i)
    end do
    helper%last   = safe(size(safe))
  end subroutine

  subroutine remove_sinusoid(func, x, n, func2)
        implicit none
        REAL(dp), PARAMETER :: PI = 3.1415926535
        integer(i4b)                                         :: i, n_periods, n
        real(dp), intent(in)                                 :: func(:), x(:)
        real(dp)                                             :: mean, max_val, min_val, intersect, omega, omega_calc, first_intersect, dx, f2i
        logical                                              :: is_first
        real(dp), intent(out), allocatable, dimension(:)     :: func2
        if (allocated(func2)) deallocate(func2)
        allocate(func2(n))
        max_val = maxval(func)
        mean = max_val - 0.5 * (max_val - minval(func))
        dx = x(2) - x(1)
        intersect = 1.0 ! just to avoid float error, fix this later
        first_intersect = 0.0
        !write(*,*) "n: ", n
        n_periods = -1
        is_first  = .true.
        func2(1) = func(1) - mean
        do i=2,n
            func2(i) = func(i) - mean
            f2i = func2(i)
            if ((f2i * func2(i-1) < 0) .and. (f2i>func2(i-1))) then
                !write(*,*) "in if test"
                n_periods = n_periods + 1
                intersect = x(i) - f2i/(f2i - func2(i-1))*dx
                if (is_first) then
                    first_intersect = intersect
                    is_first = .false.
                end if
            end if
        end do
        if (abs(intersect - first_intersect) > 0.00000001 * (2 + n_periods)) then 
           omega_calc = 2.0 * PI / ((intersect - first_intersect)/n_periods)
        else 
           omega_calc = 1.0
        end if
        write(*,*) "Omega_calc", omega_calc
        write(*,*) "maxval", max_val
        write(*,*) "mean", mean
        write(*,*) "n_periods", n_periods
        
        

        do i=1,n
            func2(i) = func2(i) - (max_val - mean) * sin(omega_calc * (x(i) - first_intersect))
            write(*,*) func2(i), func(i), x(i)
!            if (mod(i,100) == 0) then
!                write(*,*) func2(i)
!            end if
        end do
    end subroutine
end program scan_detect
