program pixsdeg
  use quiet_pixspace_mod
  use quiet_fileutils
  use quiet_healpix_mod
  implicit none
  character(len=512) :: command

  call getarg(1, command)
  select case(command)
     case("map");  call command_map
     case("rms");  call command_rms
     case("block");call command_block
     case("cov");  call command_cov
     case("help"); call command_help
     case default
        write(*,*) "Unrecognized command : '" // trim(command) // "'"
        call command_help
  end select

contains

  subroutine command_map
    implicit none
    character(len=512) :: ifname, arg, ofname
    integer(i4b)       :: nside, order, onside, steps
    type(pixinfo)      :: pinfo
    type(degrade_info) :: deg
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: imap, omap

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, arg); read(arg,*) onside
    call getarg(4, ofname)

    call read_map(imap, pixels, nside, order, ifname)
    call assert(onside <= nside, "Upgrading is not supported!")
    call prepare_degrade(deg, pixels, order, nside, onside)
    allocate(omap(size(deg%opix), size(imap,2)))
    call degrade_map(imap, deg, omap)
    call write_map(omap, deg%opix, onside, order, ofname)
    deallocate(imap, omap, pixels)
    call free_degrade_info(deg)
  end subroutine

  ! Read an N-component map, with a possibility to override the number of components.
  ! 1, 2 and 3 correspond to T, QU and TQU, with larger numbers undefined. The output
  ! is 1: TT, 2: QQ, UU, QU, 3: TT, QQ, UU, TQ, TU, QU. For higher numbers of components,
  ! the result is undefined.
  subroutine command_rms
    implicit none
    character(len=512) :: ifname, arg, ofname
    integer(i4b)       :: nside, order, onside, steps, ncomp2, ncomp, a, b, off, i
    type(pixinfo)      :: pinfo
    type(degrade_info) :: deg
    integer(i4b), dimension(:),     allocatable :: pixels
    real(dp),     dimension(:,:),   allocatable :: irms, omap
    real(dp),     dimension(:,:,:), allocatable :: oblock

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, arg); read(arg,*) onside
    call getarg(4, ofname)
    if(iargc() >= 5) then
       call getarg(5, arg); read(arg,*) ncomp2
    else
       ncomp2 = 0
    end if

    call read_map(irms, pixels, nside, order, ifname)
    call assert(onside <= nside, "Upgrading is not supported!")
    call prepare_degrade(deg, pixels, order, nside, onside)

    ! Allow overriding number of components
    ncomp = count(healok(irms(1,:)))
    if(ncomp2 > 0) then
       if(ncomp < ncomp2) then
          call error("Asked for " // trim(itoa(ncomp2)) // " components, but only " // trim(itoa(ncomp)) // " present!")
       else
          ncomp = ncomp2
       end if
    end if
    ! Hack: Treat 1, 2, 3 components as T, QU, TQU
    off   = 0; if(ncomp == 2) off = 1
    allocate(oblock(deg%n, ncomp, ncomp))
    ! Degrade applies to variance, not deviation.
    irms = irms ** 2
    call degrade_noise(irms(:,1+off:), deg, oblock)

    ! Flatten to healpix format
    allocate(omap(deg%n, ncomp*(ncomp+1)/2))
    do i = 1, deg%n; call smat2vec(oblock(i,:,:), omap(i,:)); end do
    deallocate(oblock, irms)

    call write_map(omap, deg%opix, onside, order, ofname)
    deallocate(omap, pixels)
    call free_degrade_info(deg)
  end subroutine

  subroutine command_block
    implicit none
    character(len=512) :: ifname, arg, ofname
    integer(i4b)       :: nside, order, onside, steps, ncomp, a, b, off, i
    type(pixinfo)      :: pinfo
    type(degrade_info) :: deg
    integer(i4b), dimension(:),     allocatable :: pixels
    real(dp),     dimension(:,:),   allocatable :: map
    real(dp),     dimension(:,:,:), allocatable :: iblock, oblock

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, arg); read(arg,*) onside
    call getarg(4, ofname)

    call read_map(map, pixels, nside, order, ifname)
    call assert(onside <= nside, "Upgrading is not supported!")
    call prepare_degrade(deg, pixels, order, nside, onside)

    ncomp = (nint(sqrt(real(1+8*size(map,2),dp)))-1)/2
    allocate(iblock(size(pixels), ncomp, ncomp), oblock(deg%n, ncomp, ncomp))
    do i = 1, size(pixels); call vec2smat(map(i,:), iblock(i,:,:)); end do
    call degrade_noise(iblock, deg, oblock)
    do i = 1, deg%n; call smat2vec(oblock(i,:,:), map(i,:)); end do

    call write_map(map, deg%opix, onside, order, ofname)
    deallocate(map, pixels, iblock, oblock)
    call free_degrade_info(deg)
  end subroutine

  subroutine command_cov
    implicit none
    character(len=512) :: ifname, arg, ofname
    integer(i4b)       :: nside, order, onside, steps
    type(pixinfo)      :: pinfo
    type(degrade_info) :: deg
    integer(i4b), dimension(:),       allocatable :: pixels
    real(dp),     dimension(:,:,:,:), allocatable :: icov, ocov

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, arg); read(arg,*) onside
    call getarg(4, ofname)

    call dmem("start")
    call read_covmat(icov, pixels, nside, order, ifname)
    call dmem("read")
    call assert(onside <= nside, "Upgrading is not supported!")
    call prepare_degrade(deg, pixels, order, nside, onside)
    allocate(ocov(size(deg%opix), size(icov,2), size(deg%opix), size(icov,4)))
    call dmem("prepare " // trim(itoa(size(ocov,1))) // "," // trim(itoa(size(ocov,2))) // &
         & "," // trim(itoa(size(ocov,3))) // "," // trim(itoa(size(ocov,4))))
    call degrade_noise(icov, deg, ocov)
    call dmem("degrade")
    call write_covmat(ocov, deg%opix, onside, order, ofname)
    call dmem("write")
    deallocate(icov, ocov, pixels)
    call free_degrade_info(deg)
  end subroutine


  subroutine command_help
    implicit none
    write(*,*) "pixdeg: Degrades maps and noise."
    write(*,*) "Commands:"
    write(*,*) "  map inmap.fits nside outmap.fits"
    write(*,*) "    Degrade inmap.fits to nside into outmap.fits."
    write(*,*) "  rms irms.fits nside oblock.fits [ncomp]"
    write(*,*) "    Degrade irms.fits to nside. The result is in block-diagonal"
    write(*,*) "    format, with one small covmat per pixel. This is stored as"
    write(*,*) "    [TT] for 1 component, [QQ,UU,QU] for 2 and [TT,QQ,UU,TQ,TU,QU] for 3."
    write(*,*) "    Use ncomp to override the number of components considered."
    write(*,*) "  block iblock.fits nside oblock.fits"
    write(*,*) "    Degrade iblock.fits to nside into oblock.fits, where both"
    write(*,*) "    input and output is in block-diagonal format."
    write(*,*) "  cov incov.unf  nside outcov.unf"
    write(*,*) "    Degrade incov.unf to nside into outcov.unf."
  end subroutine


  !--------- Helper functions below this -----------

  subroutine error(str)
    implicit none
    character(len=*) str
    write(*,*) str
    stop
  end subroutine

end program
