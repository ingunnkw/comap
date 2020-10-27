! A new jackknife module, which also includes some kind of accept lists.

module comap_jackknife_mod
  use quiet_utils
  use quiet_hdf_mod
  !use comap_scan_mod
  !use comap_detector_mod
  !use comap_frequency_mod
  implicit none

  !interface is_accepted_scan_any
  !   module procedure is_accepted
  !end interface

  type jk_type
     integer(i4b) :: nscans, njk, n_feed, n_coadd, n_split
     character(len=4),   allocatable, dimension(:)       :: jk_name                 ! (njk)
     integer(i4b),       allocatable, dimension(:)       :: scan_list, feedmap      ! (nscans) / (njk)
     integer(i4b),       allocatable, dimension(:,:,:,:) :: split                   ! (njk,nsb,nfeed,nscans)
     integer(i4b),       allocatable, dimension(:,:,:)   :: jk_list, accept_list    ! (nsb,nfeed,nscans)
  end type jk_type

  character(len=512) :: acceptlist, jk_definition

contains

  ! Initialize an empty jk_type
  subroutine initialize_empty_jk(jk)
    implicit none
    type(jk_type), intent(inout) :: jk

    integer(i4b) :: nfeed, nsb

    nfeed = 20
    nsb = 4

    call free_jk_type(jk)

    jk%njk = 0
    jk%nscans = 0

    allocate(jk%jk_list(nsb,nfeed,1), jk%accept_list(nsb,nfeed,1))

    jk%scan_list = 0000002
    jk%jk_list = 1
    jk%accept_list = .true.

    allocate(jk%split(jk%njk, nsb, nfeed, jk%nscans))
    jk%split = 0
    
  end subroutine initialize_empty_jk


  ! Read acceptlist
  subroutine read_acceptlist(jk, acceptlist, jk_definition)
    implicit none
    type(jk_type),      intent(inout)  :: jk
    character(len=512), intent(in)     :: acceptlist, jk_definition

    type(hdf_file)     :: h5file
    character(len=512) :: line
    integer(i4b)       :: i, j, feed, sb, num, nfeed, nsb, temp, ext(3), unit, feedmap
    
    call free_jk_type(jk)
    
    call open_hdf_file(trim(acceptlist), h5file, "r")

    call get_size_hdf(h5file, "jk_list", ext)
    nsb = ext(1); nfeed = ext(2); jk%nscans = ext(3)
    
    allocate(jk%jk_list(nsb, nfeed, jk%nscans), &
         & jk%accept_list(nsb, nfeed, jk%nscans), &
         & jk%scan_list(jk%nscans))

    call read_hdf(h5file, "jk_list", jk%jk_list)
    call read_hdf(h5file, "scan_list", jk%scan_list)
    call read_hdf(h5file, "accept_list", jk%accept_list)

    call close_hdf_file(h5file)

    ! Read jackknife definition file
    jk%n_feed = 0; jk%n_coadd = 0; jk%n_split = 0
    unit = getlun()
    open(unit, file=jk_definition, action="read",status="old")
    read(unit,*) num
    jk%njk = num-1
    allocate(jk%jk_name(jk%njk), jk%feedmap(jk%njk))
    read(unit,*) line
    do i = 2, num
       read(unit,*) jk%jk_name(i-1), jk%feedmap(i-1)!line
       !jk%jk_name(i-1) = line(1:4)
       if (jk%feedmap(i-1) == 0) then
          jk%n_coadd = jk%n_coadd + 1
       else if (jk%feedmap(i-1) == 2) then
          jk%n_split = jk%n_split + 1
       else
          jk%n_feed = jk%n_feed + 1
       end if
    end do
    close(unit)
    !write(*,*) jk%n_coadd, jk%n_feed
    !stop

    allocate(jk%split(jk%njk, nsb, nfeed, jk%nscans))

    ! Define jackknife splitting
    do j = 1, jk%nscans
        do feed = 1, nfeed
           do sb = 1, nsb
              temp = jk%jk_list(sb,feed,j)
              do i = jk%njk, 1, -1
                 if (temp > 2**i) then 
                    jk%split(i,sb,feed,j) = 1
                    temp = temp - 2**i
                 else
                    jk%split(i,sb,feed,j) = 0
                 end if
              end do
           end do
        end do
     end do
     !write(*,*) jk%split
  end subroutine read_acceptlist


  ! Free jk_type
  subroutine free_jk_type(jk)
    implicit none
    type(jk_type) :: jk

    if (allocated(jk%jk_list))     deallocate(jk%jk_list)
    if (allocated(jk%scan_list))   deallocate(jk%scan_list)
    if (allocated(jk%accept_list)) deallocate(jk%accept_list)
    if (allocated(jk%jk_name))     deallocate(jk%jk_name)
    if (allocated(jk%feedmap))     deallocate(jk%feedmap)
    if (allocated(jk%split))       deallocate(jk%split)

  end subroutine free_jk_type


end module comap_jackknife_mod
