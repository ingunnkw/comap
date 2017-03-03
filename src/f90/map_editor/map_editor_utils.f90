module map_editor_utils
  use healpix_types
  use fitstools
  use udgrade_nr
  use head_fits
  use pix_tools
  implicit none

contains

  ! Routine for reading the first map
  subroutine initialize_single_map(filename, nside, ordering, nmaps, header, map)
    implicit none

    character(len=128),                         intent(in)  :: filename
    integer(i4b),                               intent(out) :: nside, ordering, nmaps
    real(dp),           pointer, dimension(:,:)             :: map
    character(len=80),           dimension(180)             :: header

    integer(i4b) :: temp, npix
    logical(lgt) :: anynull
    real(dp)     :: nullval

    temp = getsize_fits(filename, nside=nside, ordering=ordering, nmaps=nmaps)

    npix = nside2npix(nside)

    allocate(map(0:npix-1,nmaps))

    call read_bintab(filename, map, npix, nmaps, nullval, anynull, header=header)

  end subroutine initialize_single_map

  ! Routine for reading a second map, and convert its properties to those of the first map
  subroutine initialize_second_map(filename, nside, ordering, nmaps, map)
    implicit none

    character(len=128),                         intent(in)  :: filename
    integer(i4b),                               intent(in)  :: nside, ordering
    real(dp),           pointer, dimension(:,:)             :: map

    real(dp)     :: nullval
    integer(i4b) :: i, temp, nmaps, npix, npix_in, nside_in, ordering_in, nmaps_in
    logical(lgt) :: anynull

    real(dp), allocatable, dimension(:,:) :: map_in

    ! Get general information for the file
    temp = getsize_fits(filename, nside=nside_in, ordering=ordering_in, nmaps=nmaps_in)

    npix_in = nside2npix(nside_in)
    npix    = nside2npix(nside)

    if (nmaps_in /= nmaps) then
       write(*,*) 'Warning: Nmaps differ -- setting output Nmaps to ', nmaps_in
       nmaps = nmaps_in
    end if

    if (nside_in /= nside) then

       allocate(map_in(0:npix_in-1,nmaps))
       allocate(map(0:npix-1,nmaps))
       call read_bintab(filename, map_in, npix_in, nmaps, nullval, anynull)

       write(*,*) 'Warning: Nsides differ -- setting output Nside to ', nside
       if (ordering_in == 1) then
          do i = 1, nmaps
             call udgrade_ring(map_in(:,i), nside_in, map(:,i), nside)
          end do
       else
          do i = 1, nmaps
             call udgrade_nest(map_in(:,i), nside_in, map(:,i), nside)
          end do
       end if

       deallocate(map_in)

    else

       allocate(map(0:npix-1,nmaps))
       call read_bintab(filename, map, npix, nmaps, nullval, anynull)

    end if

    if (ordering_in /= ordering) then
       write(*,*) 'Warning: Orderings differ -- setting output ordering to ', ordering
       
       do i = 1, nmaps
          if (ordering_in == 1) then
             call convert_ring2nest(nside, map(0:npix-1,i))
          else
             call convert_nest2ring(nside, map(0:npix-1,i))             
          end if
       end do
    end if

  end subroutine initialize_second_map

  ! Routine for reading a mask, and convert properties the those of the first map
  subroutine initialize_mask(filename, nside, ordering, mask)
    implicit none

    character(len=128),                       intent(in)  :: filename
    integer(i4b),                             intent(in)  :: nside, ordering
    logical(lgt),       pointer, dimension(:)             :: mask

    real(dp)     :: nullval
    integer(i4b) :: i, j, unit, npix, npix_in, nside_in, ordering_in, fac, nmaps_in
    logical(lgt) :: anynull

    real(dp), allocatable, dimension(:,:) :: mask_in

    npix    = nside2npix(nside)

    npix_in = getsize_fits(filename, nside=nside_in, ordering=ordering_in, nmaps=nmaps_in)
    ! npix_in = nside2npix(nside_in)
    allocate(mask_in(0:npix_in-1,nmaps_in))
    call read_bintab(trim(filename), mask_in, npix_in, nmaps_in, nullval, anynull)

    if (ordering_in /= ordering) then
       write(*,*) 'Warning: Orderings differ -- setting output ordering to ', ordering
       
       do i = 1, nmaps_in
          if (ordering == 1) then
             call convert_nest2ring(nside, mask_in(:,i))
          else
             call convert_ring2nest(nside, mask_in(:,i))
          end if
       end do
    end if

    allocate(mask(0:npix-1))
    if (nside_in /= nside) then

       write(*,*) 'Warning: Nsides differ -- setting output Nside to ', nside
       if (nside_in > nside) then
          ! Degrade mask -- require that all pixels inside a large pixels are included
          fac = (nside_in/nside)**2
          mask = .true.
          do i = 0, npix-1
             do j = 0, fac-1
                if (.not. mask_in(i*fac+j,1) == 1.) then
                   mask(i) = .false.
                   exit
                end if
             end do
          end do
       else
          ! Upgrade mask -- set all pixels inside an included large pixel to .true.
          fac = (nside/nside_in)**2
          do i = 0, npix-1
             do j = 0, fac-1
                if (mask_in(i/fac,1) == 1.) then
                   mask(i) = .true.
                else
                   mask(i) = .false.
                end if
             end do
          end do          
       end if

    else

       mask = (mask_in(:,1) == 1.)

    end if
    deallocate(mask_in)

  end subroutine initialize_mask


  ! Routine for writing the output map
  subroutine write_result_map(filename, nside, ordering, header, map, double_precision)
    implicit none

    character(len=*),                     intent(in)    :: filename
    integer(i4b),                         intent(in)    :: nside, ordering
    character(len=80),  dimension(180),   intent(inout) :: header
    real(dp),           dimension(0:,1:), intent(in)    :: map
    logical(lgt),                         intent(in), optional :: double_precision

    integer(i4b) :: i, nlheader, nmaps, npix, j
    character(len=80)                   :: line
    character(len=80), dimension(1)     :: line2
    character(len=256) :: outfile

    npix  = size(map(:,1))
    nmaps = size(map(0,:))

    ! Modify necessery header keywords; leave all others as they were
    nlheader = size(header)
    do i=1,nlheader
       line  = header(i)
       line2 = ''
       if (line(1:8) == 'ORDERING') then
          if (ordering == 1) then
             call add_card(line2,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
          else
             call add_card(line2,"ORDERING","NESTED",  "Pixel ordering scheme, either RING or NESTED")
          end if
       end if

       if (line(1:5) == 'NSIDE') then
          call add_card(line2,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
       end if

       if (line(1:7) == 'LASTPIX') then
          call add_card(line2,"LASTPIX",npix-1,"Last pixel # (0 based)")
       end if

       if (line(1:7) == 'POLAR') then
          if (nmaps == 1) then
             call add_card(line2,"POLAR",".false."," Polarisation included (True/False)")
          else
             call add_card(line2,"POLAR",".true."," Polarisation included (True/False)")
          end if
       end if

       if (line(1:5) == 'TUNIT') then
          read(line(6:7),*) j
          if (j > nmaps) header(i) = ''
       end if

       if (trim(line2(1)) /= '') header(i) = line2(1)
    enddo

!!$    call add_card(header,"COMMENT","-----------------------------------------------")
!!$    call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
!!$    call add_card(header,"COMMENT","-----------------------------------------------")
!!$    call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
!!$    if (ordering == 1) then
!!$       call add_card(header,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
!!$    else
!!$       call add_card(header,"ORDERING","NESTED",  "Pixel ordering scheme, either RING or NESTED")
!!$    end if
!!$    call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
!!$    call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
!!$    call add_card(header,"LASTPIX",npix-1,"Last pixel # (0 based)")
!!$    call add_card(header) ! blank line

!!$    call add_card(header,"COMMENT","-----------------------------------------------")
!!$    call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
!!$    call add_card(header,"COMMENT","-----------------------------------------------")
!!$    call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
!!$    call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
!!$    call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
!!$    call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
!!$    call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
!!$    call add_card(header) ! blank line
!!$    if (nmaps == 1) then
!!$       call add_card(header,"POLAR",".false."," Polarisation included (True/False)")
!!$    else
!!$       call add_card(header,"POLAR",".true."," Polarisation included (True/False)")
!!$    end if
!!$    call add_card(header) ! blank line
!!$    call add_card(header,"TTYPE1", "TEMPERATURE","Temperature map")
!!$    call add_card(header,"TUNIT1", "n/a","map unit")
!!$    call add_card(header)

!!$    if (nmaps == 2) then
!!$
!!$       call add_card(header,"TTYPE2", "Q-POLARISATION","Q Polarisation map")
!!$       call add_card(header,"TUNIT2", "n/a","map unit")
!!$       call add_card(header)
!!$
!!$    else if (nmaps == 3) then
!!$       call add_card(header,"TTYPE2", "Q-POLARISATION","Q Polarisation map")
!!$       call add_card(header,"TUNIT2", "n/a","map unit")
!!$       call add_card(header)
!!$
!!$       call add_card(header,"TTYPE3", "U-POLARISATION","U Polarisation map")
!!$       call add_card(header,"TUNIT3", "n/a","map unit")
!!$       call add_card(header)
!!$    endif
!!$    call add_card(header,"COMMENT","*************************************")

    outfile = '!' // trim(filename)
    if (present(double_precision)) then
       call write_bintab(map, npix, nmaps, header, nlheader, outfile)
    else
       call write_bintab(real(map,sp), npix, nmaps, header, nlheader, outfile)
    end if
    
  end subroutine write_result_map


  ! Subroutine for reading Healpix ring quadrature weights 
  subroutine read_ringweights(nside, nmaps, weights)
    implicit none

    integer(i4b),                          intent(in)  :: nside, nmaps
    real(dp),     pointer, dimension(:,:)              :: weights
 
    character(len=128)  :: weight_file
    character(len=5)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    allocate(weights(1:2*nside,nmaps))
    
    call int2string_v1(nside, nside_text)
    weight_file = 'weight_ring_n' // nside_text // '.fits'
    inquire(file=weight_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(weight_file, weights, 2*nside, nmaps, nullval, anynull)
       weights = 1.d0 + weights
    else
       weights = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Weight file ', trim(weight_file), ' not found. '
       write(*,*) 'Using unity weights in the spherical harmonic transforms.'
       write(*,*) ''
    end if

  end subroutine read_ringweights

  subroutine read_pixwin(nside, nmaps, pixwin)
    implicit none

    integer(i4b),                          intent(in)  :: nside, nmaps
    real(dp),     pointer, dimension(:,:)              :: pixwin

    integer(i4b)        :: nc
    character(len=128)  :: pixwin_file
    character(len=4)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    allocate(pixwin(0:4*nside,nmaps))
    
    if (nmaps == 3) then
       nc = 2
    else
       nc = 1
    end if

    call int2string_v1(nside, nside_text)
    pixwin_file = 'pixel_window_n' // nside_text // '.fits'
    inquire(file=pixwin_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(pixwin_file, pixwin(0:4*nside,1:nc), 4*nside+1, nc, nullval, anynull)
       if (nmaps == 3) pixwin(:,3) = pixwin(:,2)
    else
       pixwin = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Pixel window file ', trim(pixwin_file), ' not found. '
       write(*,*) 'Using unity weights.'
       write(*,*) ''
    end if

  end subroutine read_pixwin


  ! Small utility for converting an integer to a string
  subroutine int2string_v1(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string_v1

end module map_editor_utils
