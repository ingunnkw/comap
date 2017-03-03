! Module for handling sidelobes. These are defined as maps in
! focalplane-coordinates. Ideally, they would have a proper amplitude
! at each point, but for now they are just masks. They are read in
! from a specially formatted directory:
!
!  info.txt:
!    mjd_from1 mjd_to1 dirname1
!    ...
!  dirname1:
!    mask_000.fits
!    mask_001.fits
!    ...
! Where there is one mask file per module. Missing modules are allowed.
! These are assumed to be completely accepted. Only the first component
! of the masks are used.

module quiet_sidelobe_mod
  use quiet_utils
  use quiet_fileutils
  use quiet_module_mod
  use pix_tools
  implicit none

  type sidelobe_mask
     integer(i4b) :: nside, nmod, npix, order
     logical(dp), dimension(:,:,:), allocatable :: masks   ! (npix,nhorn,nrange) T if hit
     real(dp),    dimension(:,:,:), allocatable :: centers
     real(dp),    dimension(:,:),   allocatable :: ranges
  end type

  type(sidelobe_mask) :: sidelobes

contains

  subroutine initialize_sidelobe_mod(parfile)
    implicit none
    character(len=*)   :: parfile
    character(len=512) :: sidelobe_dir
    logical(lgt), save :: initialized = .false.
    if(initialized) return
    call initialize_module_mod(parfile)
    call get_parameter(0, parfile, "SIDELOBE_DIR", par_string=sidelobe_dir)
    call read_sidelobe(sidelobe_dir, sidelobes)
  end subroutine

  subroutine read_sidelobe(dirname, mask)
    implicit none
    character(len=*), intent(in)    :: dirname
    type(sidelobe_mask),    intent(inout) :: mask
    character(len=512)              :: name, maskfile
    integer(i4b)                    :: i, j, m, n, u1, nmod, nside, order, nmaps, npix
    logical(lgt)                    :: exist, anynull
    real(dp)                        :: range(2), nullval
    real(dp),         allocatable   :: map(:,:)
    call free_sidelobe(mask)
    nmod = get_num_modules()
    u1   = getlun()
    n    = 0
    if(dirname == '') goto 2
    open(u1, file=trim(dirname) // "/info.txt")
    ! Count number of ranges
    n = 0
    do
       read(u1,*,end=1) range, name
       n = n+1
    end do
1   rewind(u1)
    allocate(mask%ranges(2,n))
    do i = 1, n
       read(u1,*) mask%ranges(:,i), name
       do j = 1, nmod
          maskfile = trim(dirname) // "/" // trim(name) // "/mask_" // trim(itoa(j-1,3)) // ".fits"
          inquire(file=maskfile,exist=exist)
          if(.not. exist) cycle
          m = getsize_fits(maskfile, nside=nside, ordering=order, nmaps=nmaps)
          npix = 12*nside**2
          if(.not. allocated(mask%masks)) then
             allocate(mask%masks(0:npix-1,0:nmod-1,n))
             mask%masks = .false.
             mask%nside = nside
             mask%order = order
             mask%npix  = npix
          else
             call assert(nside == mask%nside, "Inconsistent nside in sunmaps!")
             call assert(order == mask%order, "Inconsistent order in sunmaps!")
          end if
          allocate(map(0:npix-1,1))
          call read_bintab(maskfile, map, npix, 1, nullval, anynull)
          mask%masks(:,j-1,i) = map(:,1) == 0
          deallocate(map)
       end do
    end do
    close(u1)
2   if(.not. allocated(mask%masks)) allocate(mask%masks(0:npix-1,0:nmod-1,n))
    allocate(mask%centers(3,0:nmod-1,n))
    call calc_sidelobe_centers(mask)
  end subroutine

  subroutine free_sidelobe(mask)
    implicit none
    type(sidelobe_mask) :: mask
    if(allocated(mask%masks))   deallocate(mask%masks)
    if(allocated(mask%ranges))  deallocate(mask%ranges)
  end subroutine

  subroutine calc_sidelobe_centers(mask)
    implicit none
    type(sidelobe_mask) :: mask
    integer(i4b)   :: i, j, k, n
    real(dp)       :: v(3)
    do i = 1, size(mask%masks,3)
       do j = 0, size(mask%masks,2)-1
          v = 0
          n = 0
          do k = 0, mask%npix-1
             if(.not. mask%masks(k,j,i)) cycle
             v = v + pix2vec(mask%nside, mask%order, k)
             n = n + 1
          end do
          if(n > 0) then
             mask%centers(:,j,i) = v/sqrt(sum(v**2))
          else
             mask%centers(:,j,i) = nan
          end if
       end do
    end do
  end subroutine

  function lookup_sidelobe_range(mask, mjd) result(ind)
    implicit none
    type(sidelobe_mask):: mask
    real(dp)     :: mjd
    integer(i4b) :: ind
    do ind = 1, size(mask%ranges,2)
       if(mask%ranges(1,ind) <= mjd .and. mask%ranges(2,ind) > mjd) exit
    end do
    if(ind > size(mask%ranges,2)) ind = 0
  end function

  ! Given focal-plane phi and theta, return T if we hit the sidelobe
  function lookup_sidelobe(mask, mjd, phi, theta, mod) result(res)
    implicit none
    type(sidelobe_mask):: mask
    real(dp)     :: mjd, phi, theta
    logical(lgt) :: res
    integer(i4b) :: pix, ind, mod
    res = .false.
    ind = lookup_sidelobe_range(mask, mjd)
    if(ind == 0) return
    pix = ang2pix(mask%nside, mask%order, theta, phi)
    res = mask%masks(pix, mod, ind)
  end function
end module
