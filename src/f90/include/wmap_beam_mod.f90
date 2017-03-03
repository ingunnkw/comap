module wmap_beam_mod
  use healpix_types
  use pix_tools
  use quiet_utils
  use spline_2D_mod
  implicit none

  integer(i4b),                                      private    :: DA_ind

  real(dp),         allocatable, dimension(:),       private :: x_a, y_a, x_b, y_b
  real(dp),         allocatable, dimension(:,:,:,:), private :: coeff_a, coeff_b

  real(dp),                      dimension(10),      private  :: x_centers_a, y_centers_a, x_centers_b, y_centers_b, radius
  character(len=3),              dimension(10),      private  :: DAs

contains


!-------------------------------------------------------------------------------
!  subroutine initialize_wmap_beam_mod(unit, paramfile)
!------------------------------------------------------------------------------
 

  subroutine initialize_wmap_beam_mod(unit, paramfile)
    implicit none

    integer(i4b),     intent(in) :: unit
    character(len=*), intent(in) :: paramfile

    integer(i4b)                 :: i, j, m, n, nbeam, ind_start_x, ind_start_y, ind_stop_x, ind_stop_y
    character(len=128)           :: beamfile 
    character(len=3)             :: DA
    real(dp), allocatable, dimension(:)   :: x, y
    real(dp), allocatable, dimension(:,:) :: beam, beam_small
    
    DAs         = ['K1 ',       'Ka1',      'Q1 ',       'Q2 ',       'V1 ',         &
         & 'V2 ',        'W1 ',        'W2 ',        'W3 ',        'W4 '       ]
    radius       = [     7.d0,      5.5d0,       5.d0,       5.d0,       4.d0,         &
         & 4.d0,        3.5d0,       3.5d0,        3.5d0,       3.5d0]
    x_centers_a = [-2.29277d0,  2.19830d0,  1.80821d0, -1.82996d0,  1.90017d0,   -1.91269d0,   &
         & 0.526052d0,  0.543939d0, -0.562617d0, -0.562357d0]
    x_centers_b = [ 2.17751d0, -2.29527d0, -1.91339d0,  1.72696d0, -2.00839d0,    1.80151d0,  &
         & -0.657912d0, -0.663845d0,  0.440276d0,  0.430521d0]
    y_centers_a = [-2.79768d0, -2.66094d0,  1.80101d0,  1.80585d0, -0.0963288d0, -0.106879d0, &
         & -0.537606d0,  0.567084d0,  0.550837d0, -0.551170d0]
    y_centers_b = [-2.89328d0, -2.77340d0,  1.73234d0,  1.74449d0, -0.192072d0,  -0.177803d0, &
         & -0.634589d0,  0.480468d0,  0.484744d0, -0.629839d0]

    call get_parameter(unit, paramfile, 'BEAM_FILE',    par_string=beamfile)  
    call get_parameter(unit, paramfile, 'DA',           par_string=DA)  

    write(*,*) DA, 'DA'
    write(*,*) beamfile, 'beamfile'


    ! Find DA index
    do i = 1, 10
       if (trim(DA) == trim(DAs(i))) then
          DA_ind = i
          exit
       end if
    end do
    if (i > 10) then
       write(*,*) 'Unknown DA = ', trim(DA)
       stop
    end if


    ! Find number of beam pixels
    open(unit, file=trim(beamfile))
    i = 0
    do while (.true.)
       read(unit,*,end=1) 
       i = i + 1
    end do
1   close(unit)

    nbeam = nint(sqrt(real(i/2,dp)))
    write(*,*) 'nbeam = ', nbeam

    allocate(x(nbeam))
    allocate(y(nbeam))
    allocate(beam(nbeam, nbeam))

    open(unit, file=trim(beamfile))
    do i = 1, nbeam
       do j = 1, nbeam
          read(unit,*,end=1) x(i), y(j), beam(i,j)
       end do
    end do

    ind_start_x = 1
    ind_stop_x  = nbeam
    ind_start_y = 1
    ind_stop_y  = nbeam

    do while (sum(abs(beam(ind_start_x,:))) == 0.d0)
       ind_start_x = ind_start_x + 1
    end do

    do while (sum(abs(beam(ind_stop_x,:))) == 0.d0) 
       ind_stop_x = ind_stop_x - 1
    end do

    do while (sum(abs(beam(:,ind_start_y))) == 0.d0) 
       ind_start_y = ind_start_y + 1
    end do

    do while (sum(abs(beam(:,ind_stop_y))) == 0.d0)
       ind_stop_y = ind_stop_y - 1
    end do
    
    m = ind_stop_x - ind_start_x + 1
    n = ind_stop_y - ind_start_y + 1

    write(*,*) 'Number of pixels in A side beam = ', m, ' * ', n

    allocate(beam_small(m, n))
    allocate(x_a(m))
    allocate(y_a(n))
    allocate(coeff_a(4, 4, m, n))
    
    beam_small = beam(ind_start_x:ind_stop_x, ind_start_y:ind_stop_y)
    x_a      = x(ind_start_x:ind_stop_x)
    y_a      = y(ind_start_y:ind_stop_y)
    call splie2_full_precomp(x_a, y_a, beam_small, coeff_a)
    deallocate(beam_small)

    ! Read B-side beam
    do i = 1, nbeam
       do j = 1, nbeam
          read(unit,*,end=1) x(i), y(j), beam(i,j)
       end do
    end do

    ind_start_x = 1
    ind_stop_x  = nbeam
    ind_start_y = 1
    ind_stop_y  = nbeam

    do while (sum(abs(beam(ind_start_x,:))) == 0.d0) 
       ind_start_x = ind_start_x + 1
    end do

    do while (sum(abs(beam(ind_stop_x,:))) == 0.d0) 
       ind_stop_x = ind_stop_x - 1
    end do

    do while (sum(abs(beam(:,ind_start_y))) == 0.d0) 
       ind_start_y = ind_start_y + 1
    end do

    do while (sum(abs(beam(:,ind_stop_y))) == 0.d0) 
       ind_stop_y = ind_stop_y - 1
    end do

    m = ind_stop_x - ind_start_x + 1
    n = ind_stop_y - ind_start_y + 1

    write(*,*) 'Number of pixels in B side beam = ', m, ' * ', n

    allocate(beam_small(m, n))
    allocate(x_b(m))
    allocate(y_b(n))
    allocate(coeff_b(4, 4, m, n))
    
    beam_small = beam(ind_start_x:ind_stop_x, ind_start_y:ind_stop_y)
    x_b      = x(ind_start_x:ind_stop_x)
    y_b      = y(ind_start_y:ind_stop_y)

    call splie2_full_precomp(x_b, y_b, beam_small, coeff_b)

    close(unit)

    deallocate(beam_small)
    deallocate(x)
    deallocate(y)
    deallocate(beam)

  end subroutine initialize_wmap_beam_mod




! -----------------------------------------------------------    
! subroutine cleanup_wmap_beam_mod
!-------------------------------------------------------------

  subroutine cleanup_wmap_beam_mod
    implicit none

    if (allocated(coeff_a)) deallocate(coeff_a)
    if (allocated(coeff_b)) deallocate(coeff_b)
    if (allocated(x_a))     deallocate(x_a)
    if (allocated(y_a))     deallocate(y_a)
    if (allocated(x_b))     deallocate(x_b)
    if (allocated(y_b))     deallocate(y_b)

  end subroutine cleanup_wmap_beam_mod


!---------------------------------------------------------------------------
!  function get_beam(u, v, aorb)
!--------------------------------------------------------------------------
  
  function get_beam(u, v, aorb)
    implicit none
    
    real(dp),     intent(in) :: u, v   ! distance in degrees from beam center
    integer(i4b), intent(in) :: aorb
    real(dp)                 :: get_beam

    ! Do spline interpolation
    if (aorb == 1) then
       get_beam = splin2_full_precomp(x_a, y_a, coeff_a, u+x_centers_a(DA_ind), v+y_centers_a(DA_ind))       
    else if (aorb == 2) then
       get_beam = splin2_full_precomp(x_b, y_b, coeff_b, u+x_centers_b(DA_ind), v+y_centers_b(DA_ind))       
    end if

  end function get_beam



  function get_radius()
    implicit none

    real(dp) :: get_radius

    get_radius = radius(DA_ind)

  end function get_radius


end module wmap_beam_mod
