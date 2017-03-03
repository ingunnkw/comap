program l1_dump
  use healpix_types
  use quiet_fileutils
  use quiet_module_mod
  use l1_read_mod
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
     type(point_struct) :: data
     type(data_struct)  :: data_extra ! Only used to check readability
     type(chunk_type)   :: work
     integer(i4b),       dimension(:),   allocatable :: mods
     character(len=512), dimension(:),   allocatable :: files
     real(dp),           dimension(:,:), allocatable :: ranges
  end type chunk_iterator

  type(chunk_type)     :: chunk
  type(chunk_iterator) :: iterator
  integer(i4b)         :: unit, i, j, k, chunk_size, target_dt, ierr
  real(dp)             :: start
  character(len=512)   :: listfile, parfile, l1prefix, ofile
  character(len=512), dimension(:), allocatable :: filelist

  call getarg(1, parfile)

  ! These are implementation details that probably don't need to be
  ! configurable.
  chunk_size  = 10000 ! 10000 frames = 100 s
  target_dt   = 10    ! 10 ms
  start       = nan

  call get_parameter(0, parfile, "LEVEL1_DIR", par_string=l1prefix, &
   & desc="Directory where the level1 file hierarchy can be found.")
  call get_parameter(0, parfile, "LEVEL1_FILELIST", par_string=listfile, &
   & desc="List of level1-files to consider, relative to LEVEL1_DIR. Fills '\\&
   & 'the same role as L1_DATABASE, but does not rely on run/seg/obj classification from Chicago.")

  call initialize_module_mod(parfile)
  call read_alloc_filelist(listfile, filelist)

  call init_chunk_iterator(iterator, filelist, chunk_size, target_dt)
  do while(next_chunk(iterator, chunk))
     if(chunk%n > 0 .and. start /= start) start = chunk%time(1)
     do j = 1, chunk%n
        write(*,'(f13.5,2i3,i5,2i8,3f12.8)') (chunk%time(j)-start)*24*3600, chunk%mode(j), &
         & iterator%file_index, iterator%isect, iterator%icont, iterator%iglob, &
         & chunk%az(j), chunk%el(j), chunk%dk(j)
     end do
  end do
  call free_chunk(chunk)
  call free_chunk_iterator(iterator)
  deallocate(filelist)
contains

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
    iterator%ranges    = 0
    iterator%last_cont = .true.
    allocate(iterator%mods(get_num_modules()))
    do i = 1, size(iterator%mods); iterator%mods(i) = i-1; end do
    iterator%file_index = 1
    ! Read in the data, and scan until we reach the correct time
    call get_data(l1prefix, iterator%files, iterator%file_index, iterator%mods, iterator%data, iterator%data_extra)
    iterator%n = size(iterator%data%mode)
    i = 1
    iterator%ranges(iterator%file_index,:) = iterator%data%time([1,iterator%n])
    iterator%pos = i-1
    ! Assume that the previous step had the correct length
    iterator%mjd = iterator%data%time(i) - iterator%target_dt/24d0/60/60/1000
  end subroutine

  function next_chunk(iterator, chunk) result(ok)
    implicit none
    type(chunk_iterator) :: iterator
    type(chunk_type)     :: chunk
    logical(lgt)         :: ok
    integer(i4b)         :: i, j, k, m, n, unit, status
    real(dp)             :: dt
    ok   = iterator%file_index <= size(iterator%files)
    if(.not. ok) return
    n    = iterator%work%n
    unit = getlun()
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
          iterator%file_index = iterator%file_index + 1
          call get_data(l1prefix, iterator%files, iterator%file_index, iterator%mods, iterator%data, iterator%data_extra)
          if(iterator%file_index > size(iterator%files)) exit
          iterator%n   = size(iterator%data%mode)
          iterator%ranges(iterator%file_index,:) = iterator%data%time([1,iterator%n])
          iterator%pos = 1
       end if
       dt = (iterator%data%time(iterator%pos) - iterator%mjd)*24*60*60*1000
       iterator%mjd = iterator%data%time(iterator%pos)
       if(abs(nint(dt) - iterator%target_dt) > 3) then
          iterator%last_cont = .true.
          exit
       end if
       if(i == 1 .and. iterator%icont == 1) iterator%ifile = iterator%file_index
       iterator%work%ifile   = iterator%ifile
       iterator%work%mode(i) = iterator%data%mode(iterator%pos)
       iterator%work%time(i) = iterator%data%time(iterator%pos)
       iterator%work%az(i)   = iterator%data%encoder_azimuth(iterator%pos)
       iterator%work%el(i)   = iterator%data%encoder_elevation(iterator%pos)
       iterator%work%dk(i)   = iterator%data%encoder_deck(iterator%pos)
    end do
    ! Now copy over what we got.
    call copy_chunk(iterator%work, chunk, i-1)
  end function

  subroutine free_chunk_iterator(iterator)
    implicit none
    type(chunk_iterator) :: iterator
    call deallocate_point_struct(iterator%data)
    call deallocate_data_struct(iterator%data_extra)
    call free_chunk(iterator%work)
    if(allocated(iterator%files)) deallocate(iterator%files)
    if(allocated(iterator%mods))  deallocate(iterator%mods)
  end subroutine

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

  ! Try to read file ind, but possibly advance a bit if ind cannot be found.
  subroutine get_data(prefix, filelist, ind, mods, point, data)
    implicit none
    character(len=*)   :: prefix, filelist(:)
    character(len=512) :: file
    integer(i4b)       :: ind, mods(:), status, unit
    type(point_struct) :: point
    type(data_struct)  :: data

    unit = getlun()
    do
       if(ind > size(filelist)) exit
       file = trim(prefix) // "/" // trim(filelist(ind))
       call deallocate_point_struct(point)
       call deallocate_data_struct(data)
       status = 0
       call l1_read(unit, file, &
        & pointing=point, data=data, modules=mods, selector=l1sel(&
        & [point_mode, point_time, point_encoder_azimuth, &
        & point_encoder_elevation, point_encoder_deck, data_phase, data_demod], &
        & mods=[1]), status_code=status)
       if(status /= 0) write(stderr,*) "Error reading '" // trim(file) // "'!"
       if(status == 0) exit
       ind = ind + 1
    end do
  end subroutine

end program l1_dump
