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
     integer(i4b) :: nscans, njk
     character(len=4),   allocatable, dimension(:)     :: jk_name
     integer(i4b),       allocatable, dimension(:)     :: scan_list
     integer(i4b),       allocatable, dimension(:,:)   :: split
     integer(i4b),       allocatable, dimension(:,:,:) :: jk_list, accept_list
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

    allocate(jk%split(jk%njk, jk%nscans))
    jk%split = 0
    
  end subroutine initialize_empty_jk


  ! Read acceptlist
  subroutine read_acceptlist(jk, acceptlist, jk_definition)
    implicit none
    type(jk_type),      intent(inout)  :: jk
    character(len=512), intent(in)     :: acceptlist, jk_definition

    type(hdf_file)     :: h5file
    character(len=512) :: line
    integer(i4b)       :: i, j, feed, sb, num, nfeed, nsb, temp, ext(3), unit 
    
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
    unit = getlun()
    open(unit, file=jk_definition, action="read",status="old")
    read(unit,*) num
    jk%njk = num-1; allocate(jk%jk_name(jk%njk))
    read(unit,*) line
    do i = 2, num
       read(unit,*) line
       jk%jk_name(i-1) = line(1:4)
    end do
    close(unit)

    allocate(jk%split(jk%njk, jk%nscans))

    ! Define jackknife splitting
    do j = 1, jk%nscans
        do feed = 1, nfeed
           do sb = 1, nsb
              temp = jk%jk_list(sb,feed,j)
              if (temp > 0) goto 1
           end do
        end do
        if (temp == 0) write(*,*) "Entire scan is rejected", jk%scan_list(j)
1       do i = jk%njk, -1, 1
           if (temp > 2**i) then 
              jk%split(i,j) = 1
              temp = temp - 2**i
           else
              jk%split(i,j) = 0
           end if
        end do
     end do
  
  end subroutine read_acceptlist


  ! Free jk_type
  subroutine free_jk_type(jk)
    implicit none
    type(jk_type) :: jk

    if (allocated(jk%jk_list))     deallocate(jk%jk_list)
    if (allocated(jk%scan_list))   deallocate(jk%scan_list)
    if (allocated(jk%accept_list)) deallocate(jk%accept_list)
    if (allocated(jk%jk_name))     deallocate(jk%jk_name)
    if (allocated(jk%split))       deallocate(jk%split)

  end subroutine free_jk_type


end module comap_jackknife_mod
