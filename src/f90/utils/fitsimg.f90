! Read an N-dimensional fits primary image, and dump it as an ascii table.
! Syntax:
!  fitsimg file [slice]
!  The slice is a comma-separated list of fortran style range-strides,
!  for example :,:,1:100:4,2:,::3
program fitsimg
  use quiet_utils
  implicit none
  character(len=512) :: arg, file
  integer(i4b)       :: ndim, i, j, k, m, n, unit, status, bs, naxis, naxes(32)
  integer(i4b)       :: bitpix, simple, extended, gcount, pcount
  logical(lgt)       :: anyf
  integer(i4b), dimension(:,:), allocatable :: slice
  integer(i4b), dimension(:),   allocatable :: it, cum
  real(dp),     dimension(:),   allocatable :: data

  call getarg(1, file)
  unit   = getlun()
  status = 0
  call ftopen(unit, trim(file), 0, bs, status)
  call assert(status == 0, "Error opening file: " // trim(itoa(status)))
  call ftghpr(unit, size(naxes), simple, bitpix, naxis, naxes, pcount, &
   & gcount, extended, status)
  call assert(status == 0, "Error reading header: " // trim(itoa(status)))

  allocate(slice(3,naxis))
  arg = ""
  if(iargc() > 1) call getarg(2, arg)
  call parse_slice(arg, slice)
  call expand_slice(naxes, slice)

  ! Read everything into a one-dimensional array
  allocate(data(product(naxes(1:naxis))))
  call ftgpvd(unit, 1, 1, size(data), 0, data, anyf, status)
  call assert(status == 0, "Error reading data: " // trim(itoa(status)))
  ! Iterate through the slices and output
  allocate(it(naxis), cum(naxis))
  cum(1) = 1
  do i = 2, naxis; cum(i) = cum(i-1)*naxes(i-1); end do
  it = slice(1,:)
  do
     do j = 1, naxis
        write(*,fmt="(i6)",advance="no") it(j)
     end do
     write(*,fmt="(e17.8)") data(sum((it-1)*cum)+1)
     i = 1
     do while(i <= naxis)
        it(i) = it(i) + slice(3,i)
        if(it(i) <= slice(2,i)) exit
        it(i) = slice(1,i)
        i = i+1
     end do
     if(i > naxis) exit
  end do
  deallocate(it, cum, data)
  call ftclos(unit, status)

contains

  subroutine parse_slice(str, slice)
    implicit none
    character(len=*)  :: str
    integer(i4b)      :: slice(:,:), i, j, n, nctok, nktok
    character(len=64) :: ctoks(size(slice,2)), ktoks(3)
    n = size(slice, 2)
    ! Init to default state. 0=*
    do i = 1, n; slice(:,i) = (/ 0, 0, 1 /); end do
    call get_tokens(str, ",", ctoks, nctok)
    do i = 1, min(n, nctok)
       call get_tokens(ctoks(i), ":", ktoks, nktok, allow_empty=.true.)
       do j = 1, nktok
          if(ktoks(j) == "") then
             slice(j,i) = 0
          else
             read(ktoks(j),*) slice(j,i)
          end if
       end do
       slice(nktok+1:size(slice,1),i) = 0
       if(nktok == 1) slice(2,i) = slice(1,i)
    end do
  end subroutine

  subroutine expand_slice(lens, slice)
    implicit none
    integer(i4b) :: lens(:), slice(:,:)
    integer(i4b) :: i, n
    n = size(slice,2)
    do i = 1, n
       if(slice(3,i) == 0) slice(3,i) = 1
       if(slice(1,i) == 0) slice(1,i) = 1
       if(slice(2,i) == 0) slice(2,i) = lens(i)
    end do
  end subroutine

end program
