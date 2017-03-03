! Makes standard masks for each object available, based on precomputed fits files
module quiet_mask_mod
  use quiet_utils
  use quiet_fileutils
  implicit none

  character(len=512), private :: mask_dir

contains

  subroutine initialize_mask_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    call get_parameter(0, parfile, "MASK_DIR", par_string=mask_dir, desc=&
     & "Directory where standard masks for tod2map reside.")
    initialized = .true.
  end subroutine

  subroutine get_mask(name, nside, order, mask, found)
    implicit none
    character(len=*),       intent(in) :: name
    integer(i4b),           intent(in) :: nside, order
    real(dp),               intent(out):: mask(:,:)
    logical(lgt), optional, intent(out) :: found
    character(len=512)                  :: fname
    integer(i4b)                        :: i, j, k, m, n, unit, inside, iorder
    logical(lgt)                        :: exist
    real(dp),               allocatable :: tmpmap(:,:), tmpmap3(:,:)
    fname = trim(mask_dir) // "/mask_" // trim(name) // '.fits'
    inquire(file=fname,exist=exist)
    if(exist) then
       call read_map(tmpmap, iorder, fname, nside=inside)
       call set_ordering(order, iorder, tmpmap)

       ! Mask has been defined with (0:npix,3). 
       ! Some maskfiles have only 1 component, giving tmpmap dimension 1 and segfault in udgrade.
       if (size(tmpmap,2) < 3) then
          allocate(tmpmap3(size(tmpmap,1),3))
          tmpmap3 = 0.d0
          tmpmap3(:,1) = tmpmap(:,1)
          call udgrade(tmpmap3,order,mask,nside)
          deallocate(tmpmap3)
       else
          call udgrade(tmpmap, order, mask, nside)
          deallocate(tmpmap)
       end if
    else
       mask = 1
    end if
    if(present(found)) found = exist
  end subroutine

end module
