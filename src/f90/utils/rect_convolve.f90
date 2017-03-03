! Given a map, convolute it with a rectangular top-hat function
! with top and bottom along lines of constant latitude, and with
! constant physical width (i.e. its width in terms of latitude
! increases towards the poles).
program rect_convolve
  use quiet_utils
  use quiet_fft_mod
  use quiet_mapfile_mod
  implicit none
  character(len=512)             :: ifile, ofile, arg
  real(dp),          allocatable :: map(:,:), hat(:), rdata(:), flat(:,:)
  complex(dpc),      allocatable :: rcdata(:), chat(:)
  integer(i4b),      allocatable :: pixlist(:)
  integer(i4b)                   :: i, j, k, m, n, np, nring, nside, order, npmax
  integer(i4b)                   :: npc, i2, j2, nrow, cnrow
  real(dp)                       :: dphi, dtheta, theta, phi
  call getarg(1, ifile)
  call getarg(2, arg);  read(arg, *) dphi;   dphi   = dphi  *pi/180
  call getarg(3, arg);  read(arg, *) dtheta; dtheta = dtheta*pi/180
  call getarg(4, ofile)

  call read_map(map, order, ifile)
  nside = npix2nside(size(map,1))
  call assert(order == ring, "Map must be in ring order!")

  map = 1

  where(map > 1e10 .or. map < -1e10) map = 0

  ! We will do this by first transforming into a flat rectangle,
  ! then convolve, and finally transform back.
  nring = 4*nside-1
  npmax = 4*nside
  allocate(flat(npmax,(nring-1)*2))
  allocate(pixlist(npmax))
  allocate(rcdata(npmax/2+1), chat(npmax/2+1))
  do i = 1, nring
     call dmem("Expand " // trim(itoa(i)))
     call in_ring(nside, i, 0d0, pi, pixlist, np)
     allocate(hat(np), rdata(np))
     ! Get data for this ring
     rdata = map(pixlist(:np),1)
     ! Build smoothing profile
     call pix2ang_ring(nside, pixlist(1), theta, phi)
     call build_hat(hat, dphi, 1/sin(theta))
     ! Convolute them
     npc = np/2+1
     call fft(rdata, rcdata(:npc), 1)
     call fft(hat,   chat(:npc),   1)
     rcdata = rcdata*chat
     rcdata(npc+1:) = 0
     ! And expand to full resolution flat grid
     call fft(flat(:,i), rcdata, -1)
     ! Compensate for fft-scaling. Would normally be
     ! division by sqrt(N), but a bit more complicated
     ! now due to resolution differences.
     flat(:,i) = flat(:,i) * np**0.5d0/npmax
     deallocate(rdata, hat)
  end do
  deallocate(chat, rcdata)
  ! Complete to make wrapping work
  do i = 2, nring-1
     i2 = size(flat,2)+1-i
     do j = 1, npmax
        j2 = modulo(j+npmax/2-1,npmax)+1
        flat(j2,i2) = flat(j,i)
     end do
  end do
  ! Now we are ready to smooth in the other direction
  nrow  = size(flat,2)
  cnrow = nrow/2+1
  allocate(chat(cnrow), rcdata(cnrow), hat(nrow))
  call build_hat(hat, dtheta, 1d0)
  call fft(hat, chat, 1)
  do i = 1, npmax
     call dmem("Vertical " // trim(itoa(i)))
     call fft(flat(i,:), rcdata,  1)
     rcdata = rcdata * chat
     call fft(flat(i,:), rcdata, -1)
  end do
  deallocate(chat, rcdata, hat)
  ! And finally convert back to the sphere
  allocate(rcdata(npmax/2+1))
  do i = 1, nring
     call dmem("Contract " // trim(itoa(i)))
     call in_ring(nside, i, 0d0, pi, pixlist, np)
     allocate(rdata(np))
     call fft(flat(:,i), rcdata,           1)
     call fft(rdata,     rcdata(:np/2+1), -1)
     map(pixlist(:np),1) = rdata
     deallocate(rdata)
  end do
  ! Output the resulting map
  call write_map(map, order, ofile)

contains

  subroutine build_hat(hat, dphi, scale)
    implicit none
    real(dp)      :: dphi, hat(:), scale
    integer(i4b)  :: hat_width, np
    np = size(hat)
    ! Build smoothing profile (hat) for this ring
    hat_width = nint(np*scale*dphi/(4*pi)) ! half the hat width
    hat = 0
    hat(:min(np,1+hat_width)) = 1
    hat(max(1,np-1-hat_width):) = 1
  end subroutine

end program
