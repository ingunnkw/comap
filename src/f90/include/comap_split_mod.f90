! A new split module, which also includes some kind of accept lists.

module comap_split_mod
  use quiet_utils
  use quiet_hdf_mod
  !use comap_scan_mod
  !use comap_detector_mod
  !use comap_frequency_mod
  implicit none

  !interface is_accepted_scan_any
  !   module procedure is_accepted
  !end interface

  type split_type
     integer(i4b) :: nscans, nsplit, n_feed, n_coadd, nmultisplit
     character(len=4),   allocatable, dimension(:)       :: split_name                 ! (nsplit)
     integer(i4b),       allocatable, dimension(:)       :: scan_list, feedmap      ! (nscans) / (nsplit)
     integer(i4b),       allocatable, dimension(:,:,:,:) :: split                   ! (nsplit,nsb,nfeed,nscans)
     integer(i4b),       allocatable, dimension(:,:,:)   :: split_list, accept_list    ! (nsb,nfeed,nscans)
  end type split_type

  character(len=512) :: acceptlist, split_definition

contains

  ! Initialize an empty split_type
  subroutine initialize_empty_split(split)
    implicit none
    type(split_type), intent(inout) :: split

    integer(i4b) :: nfeed, nsb
 
    nfeed = 20
    nsb = 4

    call free_split_type(split)

    split%nsplit = 0
    split%nscans = 0

    allocate(split%split_list(nsb,nfeed,1), split%accept_list(nsb,nfeed,1))

    split%scan_list = 0000002
    split%split_list = 1
    split%accept_list = .true.

    allocate(split%split(split%nsplit, nsb, nfeed, split%nscans))
    split%split = 0
    
  end subroutine initialize_empty_split


  ! Read acceptlist
  subroutine read_acceptlist(split, acceptlist, split_definition)
    implicit none
    type(split_type),      intent(inout)  :: split
    character(len=512), intent(in)     :: acceptlist, split_definition

    type(hdf_file)     :: h5file
    character(len=512) :: line
    integer(i4b)       :: i, j, feed, sb, num, nfeed, nsb, temp, ext(3), unit, feedmap
    
    call free_split_type(split)
    
    call open_hdf_file(trim(acceptlist), h5file, "r")

    call get_size_hdf(h5file, "jk_list", ext)
    nsb = ext(1); nfeed = ext(2); split%nscans = ext(3)
    
    allocate(split%split_list(nsb, nfeed, split%nscans), &
         & split%accept_list(nsb, nfeed, split%nscans), &
         & split%scan_list(split%nscans))

    call read_hdf(h5file, "jk_list", split%split_list)
    call read_hdf(h5file, "scan_list", split%scan_list)
    call read_hdf(h5file, "accept_list", split%accept_list)

    call close_hdf_file(h5file)

    ! Read split definition file
    split%n_feed = 0; split%n_coadd = 0; split%nmultisplit = 0
    unit = getlun()
    open(unit, file=split_definition, action="read",status="old")
    read(unit,*) num
    split%nsplit = num-1
    allocate(split%split_name(split%nsplit), split%feedmap(split%nsplit))
    read(unit,*) line
    do i = 2, num
       read(unit,*) split%split_name(i-1), split%feedmap(i-1)!line
       print *, split%split_name(i-1), split%feedmap(i-1)
       !split%split_name(i-1) = line(1:4)
       if (split%feedmap(i-1) == 0) then
          split%n_coadd = split%n_coadd + 1
       else if (split%feedmap(i-1) == 2) then
          split%nmultisplit = split%nmultisplit + 1
       else
          split%n_feed = split%n_feed + 1
       end if
    end do
    close(unit)
    !write(*,*) split%n_coadd, split%n_feed
    !stop

    allocate(split%split(split%nsplit, nsb, nfeed, split%nscans))
    
    ! Define split
    do j = 1, split%nscans
        do feed = 1, nfeed
           do sb = 1, nsb
              temp = split%split_list(sb,feed,j)
              do i = split%nsplit, 1, -1
                 if (temp > 2**i) then 
                    split%split(i,sb,feed,j) = 1
                    temp = temp - 2**i
                 else
                    split%split(i,sb,feed,j) = 0
                 end if
              end do
           end do
        end do
     end do
     split%nsplit = split%nsplit - split%nmultisplit 

     !write(*,*) split%split
  end subroutine read_acceptlist


  ! Free split_type
  subroutine free_split_type(split)
    implicit none
    type(split_type) :: split

    if (allocated(split%split_list))     deallocate(split%split_list)
    if (allocated(split%scan_list))   deallocate(split%scan_list)
    if (allocated(split%accept_list)) deallocate(split%accept_list)
    if (allocated(split%split_name))     deallocate(split%split_name)
    if (allocated(split%feedmap))     deallocate(split%feedmap)
    if (allocated(split%split))       deallocate(split%split)

  end subroutine free_split_type


end module comap_split_mod
