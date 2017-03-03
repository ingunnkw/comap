! A pixel-space smoother of both signal and noise. Handles cut sky
! very well, becoming faster the smaller the patch and the smaller the beam,
! without any of the aliasing problems associated with harmonic
! space smoothing. On the other hand, it is horribly slow for the full-sky
! case unless the beam is tiny.

! We really should have a covmat format that includes pixel, nside, order
! and field information. So I will assume that we have that.

program pixsmooth
  use quiet_pixspace_mod
  use quiet_fileutils
  implicit none
  character(len=512)      :: command
  integer(i4b), parameter :: maxncomp = 256

  call getarg(1, command)
  select case(command)
     case("beam"); call command_beam
     case("map");  call command_map
     case("rms");  call command_rms
     case("block");call command_block
     case("cov");  call command_cov
     case("nside");call command_nside
     case("help"); call command_help
     case default
        write(*,*) "Unrecognized command : '" // trim(command) // "'"
        call command_help
  end select

contains

  subroutine command_beam
    implicit none
    character(len=512) :: beam_from, beam_to, ofname
    real(dp), dimension(:,:), allocatable :: rbeam

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, beam_from)
    call getarg(3, beam_to)
    call getarg(4, ofname)

    call setup_beam(beam_from, beam_to, rbeam)
    call write_beam_real(ofname, rbeam)
    deallocate(rbeam)
  end subroutine

  subroutine command_map
    implicit none
    character(len=512) :: ifname, bfile, ofname
    integer(i4b)       :: nside, order, i
    type(pixinfo)      :: pinfo
    integer(i4b), dimension(:),     allocatable :: pixels
    real(dp),     dimension(:,:),   allocatable :: beam, weights
    real(dp),     dimension(:,:,:), allocatable :: imap, omap

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(2, ifname)
    call getarg(3, bfile)
    call getarg(4, ofname)

    call read_map(imap, pixels, nside, order, ifname)
    call setup(bfile, nside, order, pixels, pinfo, weights)
    allocate(omap(size(imap,1),size(imap,2),size(imap,3)))
    do i = 1, size(imap,3)
       call smooth_map(imap(:,:,i), pinfo, weights, omap(:,:,i))
    end do
    call write_map(omap, pinfo%pix, pinfo%nside, pinfo%order, ofname)
    deallocate(imap,omap,weights)
    call free_pixinfo(pinfo)
  end subroutine

  subroutine command_rms
    implicit none
    character(len=512) :: ifname, bfile, ofname, arg, opt
    integer(i4b)       :: nside, order, ncomp, n, ai, comps(maxncomp), i, j
    type(pixinfo)      :: pinfo
    integer(i4b), dimension(:),       allocatable :: pixels
    real(dp),     dimension(:,:),     allocatable :: irms, beam, weights
    real(dp),     dimension(:,:,:,:), allocatable :: ocov

    ! Parse options. Why does this have to be so verbose?
    ai = 2; ncomp = 0
    do while(ai < iargc())
       call getarg(ai,opt)
       select case(opt)
          case("-sig")
             ncomp = ncomp+1
             ai    = ai+1
             call getarg(ai, arg); read(arg,*) comps(ncomp)
          case default
             call assert(opt(1:1) /= "-", "Unrecognized option " // trim(opt))
             exit
       end select
       ai = ai+1
    end do

    if(iargc()-ai+1 /= 3) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(ai, ifname)
    call getarg(ai+1, bfile)
    call getarg(ai+2, ofname)
    call read_map(irms, pixels, nside, order, ifname)
    call setup(bfile, nside, order, pixels, pinfo, weights)
    ! Find the actual number of components. Only T, QU and TQU supported.
    n     = size(pixels)
    if(ncomp == 0) then
       ncomp = size(irms,2)
       comps(1:ncomp) = irange(ncomp)
    end if
    call assert(all(comps(:ncomp) > 0 .and. comps(:ncomp) <= size(irms,2)), "Invalid component specified")
    allocate(ocov(n,ncomp,n,ncomp))
    irms = irms**2 ! we need the variance
    call smooth_noise(irms(:,comps(1:ncomp)), pinfo, weights, ocov)
    call write_covmat(ocov, pinfo%pix, pinfo%nside, pinfo%order, ofname)
    deallocate(irms,ocov,weights)
    call free_pixinfo(pinfo)
  end subroutine

  subroutine command_block
    implicit none
    character(len=512) :: ifname, bfile, ofname, arg, opt
    integer(i4b)       :: nside, order, ncomp, n, ai, comps(maxncomp),i,j
    type(pixinfo)      :: pinfo
    integer(i4b), dimension(:),       allocatable :: pixels
    real(dp),     dimension(:,:),     allocatable :: beam, weights
    real(dp),     dimension(:,:,:),   allocatable :: icov
    real(dp),     dimension(:,:,:,:), allocatable :: ocov

    ! Parse options. Why does this have to be so verbose?
    ai = 2; ncomp = 0
    do while(ai < iargc())
       call getarg(ai,opt)
       select case(opt)
          case("-sig")
             ncomp = ncomp+1
             ai    = ai+1
             call getarg(ai, arg); read(arg,*) comps(ncomp)
          case default
             call assert(opt(1:1) /= "-", "Unrecognized option " // trim(opt))
             exit
       end select
       ai = ai+1
    end do

    if(iargc()-ai+1 /= 3) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(ai,   ifname)
    call getarg(ai+1, bfile)
    call getarg(ai+2, ofname)
    call read_map(icov, pixels, nside, order, ifname)
    call setup(bfile, nside, order, pixels, pinfo, weights)
    ! Find the actual number of components. Only T, QU and TQU supported.
    n     = size(pixels)
    if(ncomp == 0) then
       ncomp = size(icov,2)
       comps(1:ncomp) = irange(ncomp)
    end if
    call assert(all(comps(:ncomp) > 0 .and. comps(:ncomp) <= size(icov,2)), "Invalid component specified")
    allocate(ocov(n,ncomp,n,ncomp))
    call smooth_noise(icov(:,comps(:ncomp),comps(:ncomp)), pinfo, weights, ocov)
    call write_covmat(ocov, pinfo%pix, pinfo%nside, pinfo%order, ofname)
    deallocate(icov,ocov,weights)
    call free_pixinfo(pinfo)
  end subroutine

  subroutine command_cov
    implicit none
    character(len=512) :: ifname, bfile, ofname, opt, arg
    integer(i4b)       :: nside, order, ai, comps(maxncomp), ncomp
    type(pixinfo)      :: pinfo
    integer(i4b), dimension(:),       allocatable :: pixels
    real(dp),     dimension(:,:,:,:), allocatable :: icov, ocov
    real(dp),     dimension(:,:),     allocatable :: beam, weights

    ai = 2; ncomp = 0
    do while(ai < iargc())
       call getarg(ai,opt)
       select case(opt)
          case("-sig")
             ncomp = ncomp+1
             ai    = ai+1
             call getarg(ai, arg); read(arg,*) comps(ncomp)
          case default
             call assert(opt(1:1) /= "-", "Unrecognized option " // trim(opt))
             exit
       end select
       ai = ai+1
    end do

    if(iargc()-ai+1 /= 3) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(ai,   ifname)
    call getarg(ai+1, bfile)
    call getarg(ai+2, ofname)

    call read_covmat(icov, pixels, nside, order, ifname)
    call setup(bfile, nside, order, pixels, pinfo, weights)
    allocate(ocov(size(icov,1),size(icov,ncomp),size(icov,3),size(icov,ncomp)))
    call smooth_noise(icov(:,comps(:ncomp),:,comps(:ncomp)), pinfo, weights, ocov)
    call write_covmat(ocov, pinfo%pix, pinfo%nside, pinfo%order, ofname)
    deallocate(icov,ocov,weights)
    call free_pixinfo(pinfo)
  end subroutine

  subroutine command_nside
    implicit none
    character(len=512) :: arg
    call getarg(2,arg)
    select case(arg)
       case("block"); call command_nside_block
       case default;  call error("Command nside " // trim(arg) // " not implemented")
    end select
  end subroutine

  subroutine command_nside_block
    implicit none
    character(len=512) :: ifname, bfile, ofname, arg, opt
    integer(i4b)       :: nside, order, ncomp, n, ai, comps(maxncomp), nside2
    type(degrade_info) :: deg
    integer(i4b), dimension(:),       allocatable :: pixels
    real(dp),     dimension(:,:,:),   allocatable :: icov, ocov

    ! Parse options. Why does this have to be so verbose?
    ai = 3; ncomp = 0
    do while(ai < iargc())
       call getarg(ai,opt)
       select case(opt)
          case("-sig")
             ncomp = ncomp+1
             ai    = ai+1
             call getarg(ai, arg); read(arg,*) comps(ncomp)
          case default
             call assert(opt(1:1) /= "-", "Unrecognized option " // trim(opt))
             exit
       end select
       ai = ai+1
    end do

    if(iargc()-ai+1 /= 3) then
       write(*,*) "Wrong number of arguments!"
       call command_help
       return
    end if
    call getarg(ai,   ifname)
    call getarg(ai+1, arg);   read(arg,*) nside2
    call getarg(ai+2, ofname)
    call read_map(icov, pixels, nside, order, ifname)
    call prepare_degrade(deg, pixels, order, nside, nside2)
    ! Find the actual number of components. Only T, QU and TQU supported.
    n     = size(pixels)
    if(ncomp == 0) then
       ncomp = size(icov,2)
       comps(1:ncomp) = irange(ncomp)
    end if
    call assert(all(comps(:ncomp) > 0 .and. comps(:ncomp) <= size(icov,2)), "Invalid component specified")
    allocate(ocov(deg%n,ncomp,ncomp))
    call degrade_noise(icov, deg, ocov)
    call write_map(ocov, deg%opix, nside2, order, ofname)
    deallocate(icov,ocov)
    call free_degrade_info(deg)
  end subroutine

  subroutine command_help
    implicit none
    write(*,*) "pixsmooth: Performs smoothing operations in pixel space."
    write(*,*) "Commands:"
    write(*,*) "  beam oldbeam newbeam outbeam.txt"
    write(*,*) "    Produces an effective real-space beam corresponding"
    write(*,*) "    to going from oldbeam to newbeam, which are expressed"
    write(*,*) "    in harmonic space, and are either fits files or fwhm'."
    write(*,*) "    Newbeam must be more than about 10% larger than oldbeam"
    write(*,*) "    to avoid ringing in the resulting beam."
    write(*,*) "  map inmap.fits beam.txt outmap.fits"
    write(*,*) "    Smooth inmap.fits with beam, producing outmap.fits."
    write(*,*) "  rms inrms.fits beam.txt outcov.unf [ncomp]"
    write(*,*) "    Smooth rms map inrms.fits. This produces correlations,"
    write(*,*) "    so the output is a covariance matrix."
    write(*,*) "  block inblock.hdf beam.txt outcov.unf [ncomp]"
    write(*,*) "    Smooth block-diagonal block.fits, which is (pix,comp,comp)."
    write(*,*) "  cov incov.unf  beam.txt outcov.unf"
    write(*,*) "    Smooths a covariance matrix, producing a new, even more"
    write(*,*) "    correlated covariance matrix."
  end subroutine


  !--------- Helper functions below this -----------


  subroutine setup(bfile, nside, order, pixels, pinfo, weights)
    implicit none
    character(len=*) :: bfile
    integer(i4b)     :: nside, order, pixels(:), i
    type(pixinfo)    :: pinfo
    real(dp), dimension(:,:), allocatable :: weights, beam
    real(dp)         :: maxrad
    call read_beam_real(bfile, beam)
    maxrad = radial_to_rmax(beam(:,1), beam(:,2))
    call init_pixinfo(pinfo, nside, order, maxrad, pixels)
    call alloc_weights(pinfo, weights)
    call calc_weights_radial(pinfo, beam(:,1), beam(:,2), weights, normalize=.true.)
    deallocate(beam)
  end subroutine

  subroutine setup_beam(ibeam, obeam, beam)
    implicit none
    character(len=*)                      :: ibeam, obeam
    real(dp), dimension(:,:), allocatable :: beam
    real(dp), dimension(:),   allocatable :: ilbeam, olbeam, lbeam
    integer(i4b)                          :: lmax, i
    real(dp)                              :: mval, rmax
    call beam_from_fwhm(ibeam, ilbeam)
    if(.not. allocated(ilbeam)) call beam_from_file(ibeam, ilbeam)
    if(.not. allocated(ilbeam)) call error("Error reading input beam!")

    call beam_from_fwhm(obeam, olbeam)
    if(.not. allocated(olbeam)) call beam_from_file(obeam, olbeam)
    if(.not. allocated(olbeam)) call error("Error reading output beam!")

    lmax = min(ubound(ilbeam,1),ubound(olbeam,1))
    allocate(lbeam(0:lmax))
    lbeam = olbeam(0:lmax)/ilbeam(0:lmax)

    ! Now transform to real space. But to do that, we need an idea
    ! about the size of the beam in real space. We get that by
    ! making a low-res realization first
    allocate(beam(1000,2))
    do i = 1, size(beam,1)
       beam(i,1) = pi*(i-1)/size(beam,1)
    end do
    call beam_to_radial(lbeam, beam(:,1), beam(:,2))
    mval = maxval(abs(beam))
    do i = size(beam,1),1,-1
       if(abs(beam(i,2)) > mval*1e-3) exit
    end do
    rmax = beam(i,1)
    ! Ok, now make the final beam. 1000 points is probably overkill,
    ! but we can adjust that later. If the splint is too slow, find
    ! a way to use the constant interval splint.
    do i = 1, size(beam,1)
       beam(i,1) = rmax*(i-1)/(size(beam,1)-1)
    end do
    call beam_to_radial(lbeam, beam(:,1), beam(:,2))
    deallocate(ilbeam, olbeam, lbeam)
  end subroutine

  ! Produces an alm-type beam from a string containing fwhm.
  ! Normally, one wouldn't want to send in a string, but rather
  ! a number, but this is just a helper function for the specific
  ! program above. The factor 0.5 ensures that this is an alm-type beam.
  subroutine beam_from_fwhm(fwhmstr, beam)
    implicit none
    character(len=*)                    :: fwhmstr
    real(dp), dimension(:), allocatable :: beam
    real(dp)                            :: fwhm, sigma
    integer(i4b)                        :: lmax, l
    read(fwhmstr,*,err=1) fwhm
    sigma = fwhm*pi/60/180/sqrt(8*log(2d0))
    lmax = int(10/(sigma+1e-5))
    allocate(beam(0:lmax))
    do l = 0, lmax
       beam(l) = exp(-0.5*l*(l+1)*sigma**2)
    end do
1   return
  end subroutine

  ! Our existing functions are inflexible, they only
  ! handle fixed-length beams in fits format.
  ! This could be improved to remedy that, but for now
  ! it just uses the existing one.
  ! Note 2: This framework only supports a single beam
  ! for all components.
  subroutine beam_from_file(fname, beam)
    implicit none
    character(len=*)                      :: fname
    real(dp), dimension(:),   allocatable :: beam
    real(dp), dimension(:,:), allocatable :: ibeam
    real(dp)                              :: mval
    integer(i4b)                          :: l
    allocate(ibeam(0:100000,1))
    ibeam = 0
    call read_beam(fname, ibeam)
    mval = maxval(abs(ibeam))
    do l = ubound(ibeam,1),lbound(ibeam,1),-1
       if(ibeam(l,1) /= 0) exit
    end do
    allocate(beam(0:l))
    beam = ibeam(0:l,1)
    deallocate(ibeam)
  end subroutine

  subroutine dump_mat(mat)
    implicit none
    real(dp) :: mat(:,:)
    integer(i4b) :: i, j
    do i = 1, size(mat,1)
      do j = 1, size(mat,2)
        write(*,fmt="(f8.4)",advance="no") mat(i,j)
      end do
      write(*,*)
    end do
  end subroutine

  subroutine error(str)
    implicit none
    character(len=*) str
    write(*,*) str
    stop
  end subroutine

end program
