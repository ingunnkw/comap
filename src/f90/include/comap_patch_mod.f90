module comap_patch_mod
  use healpix_types
  use quiet_utils
!  use quiet_ephem_mod
  implicit none

  character(len=512), private :: patchfile

  type patch_info
     character(len=64) :: name
     real(dp)          :: pos(2), obj_rad, img_rad, resolution, priority
!     integer(i4b)      :: eph
     logical(lgt)      :: fixed, noise_dominated, ignore
  end type patch_info

  ! Publically visible store of patches
  type(patch_info), dimension(:), allocatable :: patches

  real(dp), parameter :: patch_strong_lim = 1d0

contains

  subroutine initialize_comap_patch_mod(parfile)
    implicit none
    integer(i4b)                 :: i, n
    character(len=*), intent(in) :: parfile
    logical(lgt)                 :: exist
    logical(lgt), save           :: initialized = .false.

    if(initialized) return
    initialized = .true.
    call get_parameter(0, parfile,  'PATCH_DEFINITION_FILE', par_string=patchfile)
    inquire(file=trim(patchfile),exist=exist)
    if(.not. exist) then
       write(*,*) 'Warning: Could not find patch definition file ' // trim(patchfile) // '. Patch checking disabled.'
       allocate(patches(0))
       return
    end if
    call read_patch_info(patchfile, patches)
  end subroutine

  function get_patch_info(name, pinfo) result(found)
    implicit none
    character(len=*)    :: name
    type(patch_info)    :: pinfo
    logical(lgt)        :: found
    integer(i4b)        :: i
    found = .false.
    i = lookup_patch(name, patches)
    if(i<1) return
    pinfo = patches(i)
    found = .true.
  end function

  subroutine cleanup_comap_patch_mod()
    implicit none
    deallocate(patches)
  end subroutine

  subroutine read_patch_info(file, pinfo)
    implicit none
    type(patch_info), dimension(:), allocatable :: pinfo
    integer(i4b)       :: i, n, unit
    character(len=512) :: file, name, theta, phi
    real(dp)           :: orad, irad, res, pri
    logical(lgt)       :: exist

    if(allocated(pinfo)) deallocate(pinfo)

    unit = getlun()
    open(unit,file=trim(file),action="read",status="old")
    i = 0
    do while(.true.)
       read(unit,*,end=1,err=1) name, phi, theta, orad, irad, res, pri
       if(name(1:1) == "#") cycle
       i = i+1
    end do
1   n = i
    rewind(unit)
    allocate(pinfo(n))
    i = 0
    do while(.true.)
       read(unit,*,end=2,err=2) name, phi, theta, orad, irad, res, pri
       if(name(1:1) == "#") cycle
       i = i+1
       pinfo(i)%name = name
       pinfo(i)%obj_rad = orad * DEG2RAD
       pinfo(i)%img_rad = irad * DEG2RAD
       pinfo(i)%priority = pri
       if(trim(theta) == "x" .or. trim(phi) == "x") then
          pinfo(i)%pos = 0
          pinfo(i)%fixed = .false.
          pinfo(i)%eph = name2eph(name)
       else
          read(theta,*) pinfo(i)%pos(2)
          read(phi,*)   pinfo(i)%pos(1)
          pinfo(i)%pos = pinfo(i)%pos * DEG2RAD
          pinfo(i)%pos(2) = PI/2 - pinfo(i)%pos(2)
          pinfo(i)%fixed = .true.
          pinfo(i)%eph = 0
       end if
       pinfo(i)%noise_dominated = pinfo(i)%priority <  1
       pinfo(i)%ignore          = pinfo(i)%priority <= 0
    end do
2   close(unit)
  end subroutine

  function lookup_patch(name, pinfos) result(i)
    implicit none
    character(len=*) :: name
    type(patch_info) :: pinfos(:)
    integer(i4b)     :: i
    do i = size(pinfos), 1, -1
       if(pinfos(i)%name == name) exit
    end do
  end function

end module
