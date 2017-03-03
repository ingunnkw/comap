! Memory-heavy (i.e. loads whole covar) but flexible covariance matrix
! plotter.

program covslice
  use quiet_pixspace_mod
  use quiet_fileutils
  implicit none
  character(len=512) :: command

  call getarg(1, command)
  select case(command)
     case("diagonal"); call command_diag
     case("vertical"); call command_vert
     case("scale");    call command_scale
     case("help");     call command_help
     case default
        write(*,*) "Unrecognized command : '" // trim(command) // "'"
        call command_help
  end select

contains

  ! Read a covariance matrix and turn its diagonal into a map
  ! with components in diag + col-order
  subroutine command_diag
    implicit none
    character(len=512) :: ifname, ofname
    integer(i4b)       :: nside, order, ncomp, n, nuniq, i, j
    type(pixinfo)      :: pinfo
    integer(i4b), dimension(:),       allocatable :: pixels
    integer(i4b), dimension(:,:),     allocatable :: remap
    real(dp),     dimension(:,:,:,:), allocatable :: icov
    real(dp),     dimension(:,:),     allocatable :: omap

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, ofname)

    call read_covmat(icov, pixels, nside, order, ifname)
    n = size(icov,1); ncomp = size(icov,2)
    nuniq = ncomp*(ncomp+1)/2
    ! slice and remap
    allocate(omap(n,nuniq), remap(nuniq,2))
    call vec2smat_map(remap,ncomp)
    do i = 1, n
       do j = 1, nuniq
          omap(i,j) = icov(i,remap(j,1),i,remap(j,2))
       end do
    end do
    call write_map(omap, pixels, nside, order, ofname)
    deallocate(icov,omap,remap)
  end subroutine

  subroutine command_vert
    implicit none
    character(len=512) :: ifname, ofname, arg
    integer(i4b)       :: nside, order, ncomp, n, nuniq, i, j, k, pix, ncoord, c
    logical(lgt)       :: corr
    type(pixinfo)      :: pinfo
    integer(i4b), dimension(:),       allocatable :: pixels
    integer(i4b), dimension(:,:),     allocatable :: remap
    real(dp),     dimension(:,:,:,:), allocatable :: icov
    real(dp),     dimension(:,:,:),   allocatable :: omap
    real(dp),     dimension(:,:),     allocatable :: coords

    if(iargc() < 5) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    corr   = .false.
    ncoord = (iargc()-3)/2
    allocate(coords(2,ncoord))
    call getarg(2, ifname)
    do i = 1, ncoord
       call getarg(1+i*2, arg); read(arg,*) coords(1,i)
       call getarg(2+i*2, arg); read(arg,*) coords(2,i)
    end do
    call getarg(1+i*2, ofname)
    if(iargc() >= 2+i*2) then
       call getarg(2+i*2, arg)
       corr = arg == "corr"
    end if

    call read_covmat(icov, pixels, nside, order, ifname)
    n = size(icov,1); ncomp = size(icov,2)
    nuniq = ncomp*(ncomp+1)/2

    allocate(omap(n,nuniq,ncoord), remap(nuniq,2))
    call vec2smat_map(remap,ncomp)
    do c = 1, ncoord
       ! Find the pixel corresponding to the points asked for
       if(order == RING) then
          call ang2pix_ring(nside, (90-coords(2,c))*DEG2RAD, coords(1,c)*DEG2RAD, i)
       else
          call ang2pix_nest(nside, (90-coords(2,c))*DEG2RAD, coords(1,c)*DEG2RAD, i)
       end if
       do pix = 1, size(pixels)
          if(pixels(pix) == i) exit
       end do
       call assert(pix <= size(pixels), "Indicated point #"//trim(itoa(c))//" is not in cov!")

       ! slice and remap
       do j = 1, nuniq
          do k = 1, n
             omap(k,j,c) = icov(k,remap(j,1),pix,remap(j,2))
             if(corr) omap(k,j,c) = omap(k,j,c)/sqrt(icov(k,remap(j,1),k,remap(j,1))*icov(pix,remap(j,2),pix,remap(j,2)))
          end do
       end do
    end do
    call write_map(omap, pixels, nside, order, ofname)
    deallocate(icov,omap,remap)
  end subroutine

  subroutine command_scale
    implicit none
    character(len=512) :: ifname, ofname, arg
    integer(i4b)       :: nside, order
    integer(i4b), dimension(:),       allocatable :: pixels
    real(dp),     dimension(:,:,:,:), allocatable :: cov
    real(dp)           :: mul

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, arg); read(arg,*) mul
    call getarg(4, ofname)

    call read_covmat (cov, pixels, nside, order, ifname)
    cov = cov * mul
    call write_covmat(cov, pixels, nside, order, ofname)
    deallocate(cov, pixels)
  end subroutine

  subroutine command_help
    implicit none
    write(*,*) "covslice: Slices a covariance matrix. The covariance"
    write(*,*) " matrix must contain pixel information. If it does not,"
    write(*,*) " use addpix_cov to add them. In the case of multiple"
    write(*,*) " components, the output map will have one component for"
    write(*,*) " each component, in diag-first ordering."
    write(*,*) "Commands:"
    write(*,*) " diagonal input_cov.unf output_map.fits"
    write(*,*) " vertical input_cov.unf lon lat output_map.fits"
  end subroutine

  !--------- Helper functions below this -----------

  subroutine error(str)
    implicit none
    character(len=*) str
    write(*,*) str
    stop
  end subroutine

  ! Constructs a list of the unique component combinations
  ! for n components.
  !  n=1: [1,1]
  !  n=2: [1,1],[2,2],[1,2]
  !  n=3: [1,1],[2,2],[3,3],[1,2],[1,3],[2,3]
  ! and so on.
  subroutine vec2smat_map(map, n)
    implicit none
    integer(i4b) :: map(:,:), i, j, k, n
    do i = 1, n; map(i,:) = [ i, i]; end do
    k = n
    do i = 1, n
       do j = i+1, n
          k = k+1
          map(k,:) = [ i, j ]
       end do
    end do
  end subroutine

end program
