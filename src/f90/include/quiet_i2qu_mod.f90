module quiet_i2qu_mod
  use healpix_types
  use pix_tools
  use quiet_fileutils
  use quiet_gain_mod
  implicit none

  integer(i4b),                          private   :: nmaps, nside, ordering
  real(dp), dimension(:,:), allocatable, private   :: map
  real(dp), dimension(:,:), allocatable, private   :: leakage

contains

  subroutine initialize_i2qu(unit, parfile)
    implicit none
    integer(i4b),            intent(in) :: unit
    character(len=128),      intent(in) :: parfile

    character(len=256)                  :: tempmap_file 
    integer(i4b)                        :: all_polmod, i , j

    call assert(.false., "i2qu mod needs fixing")
    all_polmod = 17          !OBS: modify for W-band
 
    allocate(leakage(0:all_polmod-1, 0:3))

    write(*,*) 'hei'

    do i = 0, all_polmod-1
       do j = 0, 3
          !call get_calibration(mjd=0.d0, module=i, diode=j, i2qu_leak=leakage(i,j))
       end do
    end do
    leakage = leakage * 0.01      ! i2qu_leak is given in % 
    call get_parameter(unit, parfile, 'TEMPMAP_FILE', par_string=tempmap_file)
    call read_map(map, ordering, tempmap_file, nside=nside, nmap=nmaps)
    map=map/1000.d0
   end subroutine initialize_i2qu

!--------------------------------

  subroutine i2qu_correct(data, run, segment)
    implicit none

    type(module_struct), dimension(:)                      :: data
    integer(i4b),                              intent(in)  :: run, segment

    integer(i4b), dimension(:), pointer :: polmod, true_polmod
    integer(i4b)                        :: i, j , k, pix
    real(dp)                            :: temp, phi, theta
    real(dp), dimension(0:3)            :: gain
 

    call get_polmod(data, polmod, true_polmod)

!open(13, file='before_i2qu.txt')
!do i = 1, size(data(1)%time)
!   write(13,*) data(1)%tod(0,i)
!end do
!close(13)

    do i = 1, size(polmod)
       do j = 0, 3
          ! Get gain
!          call get_calibration(run=run, segment=segment, module=true_polmod(i), diode=j, gain=gain(j))
          gain(j) = get_gain(data(1)%time(1), quiet_horns(true_polmod(i))%diodes(j))
          ! Convert from mV/K to V/mK
          if (gain(j) /= 1.d0) gain(j) = gain(j) * 1.d-6
       end do
       do k = 1, size(data(1)%time)
          
          ! Get pointing
          phi   = data(polmod(i))%pointing(1,k)
          theta = data(polmod(i))%pointing(2,k)
          phi   = mod(phi + 2.d0*pi, 2.d0*pi )
              
          !Get corresponding pixel in map
          if (ordering == 1) then
             call ang2pix_ring(nside, theta, phi, pix)
          else if (ordering == 2) then
             call ang2pix_nest(nside, theta, phi, pix)
          else
             write(*,*) 'Input map is neither ringed nor nested. Quiting'
             stop
          end if
          
          ! Get temperature at that pixel/pointing
          temp = map(pix,1)

          do j = 0, 3
             ! Leakage correction
             data(i)%tod(j,k) = data(i)%tod(j,k) - leakage(true_polmod(i),j) * temp*gain(j)
          end do
       end do
    end do

!open(13, file='after_i2qu.txt')
!do i = 1, size(data(1)%time)
!   write(13,*) data(1)%tod(0,i)
!end do
!close(13)

    deallocate(polmod)
    deallocate(true_polmod)

  end subroutine i2qu_correct

!--------------------------------

  subroutine i2qu_cleanup
    implicit none
    
    if (allocated(map))     deallocate(map)
    if (allocated(leakage)) deallocate(leakage)

  end subroutine i2qu_cleanup

!---------------------------

  subroutine get_polmod(data, polmod, true_polmod)
    implicit none
    
    type(module_struct), dimension(:)                      :: data
    integer(i4b), dimension(:), pointer                    :: polmod, true_polmod

    integer(i4b)                                           :: num_temppairs, num_modules, num_polmod, i, j, k, m
    logical(lgt)                                           :: Qband, Wband, found_temp, temperature, polarization
    integer(i4b), dimension(:,:), allocatable              :: tempmod, true_tempmod
   
    temperature  = .false.
    polarization = .true.

    ! Initializing main parameters
    Qband = .true.
    Wband = .false.

    ! Initialize band-specific values
    if (Qband) then
       num_temppairs = 1
    else if (Wband) then
       num_temppairs = 3
    else
       write(*,*)'You have chosen neither Qband nor Wband. Quiting'
       stop
    end if

    allocate(tempmod(num_temppairs,2))
    allocate(true_tempmod(num_temppairs,2))
    
    if (Qband) then
       true_tempmod(1,1) = 17
       true_tempmod(1,2) = 18
    else if (Wband) then
       write(*,*)'Add Wband info here'            ! OBS Add Wband info here
    else
       write(*,*)'You have chosen neither Qband nor Wband. Quiting'
       stop
    end if
    
    ! Search for the temperature modules, putting them into tempmod
    ! Putting polarization modules into polmod

    num_modules = size(data)
    num_polmod = num_modules-2*num_temppairs
    allocate(polmod(num_polmod))
    allocate(true_polmod(num_polmod))
    m = 1
    tempmod = -1
    polmod  = -1
    
    do i = 1, num_modules
       found_temp = .false.
       do j = 1, num_temppairs
          do k = 1, 2
             if (data(i)%module_number == true_tempmod(j,k)) then
                tempmod(j,k) = i    ! Put tempmod to the module no who has the true temp module_number
                found_temp = .true.
             end if
          end do
       end do
       if (.not. found_temp) then
          true_polmod(m) = data(i)%module_number
          polmod(m) = i
          m = m + 1
       end if
    end do
    
    ! Checking to see if all temperature modes are there
    if (temperature) then
       write(*,*) 'Number of temperature modules  =', 2*num_temppairs
       do j = 1, num_temppairs
          do k = 1, 2
             if (tempmod(j,k)==-1) then
                write(*,*)'Data does not contain all temperature modules'
                stop
             end if
          end do
       end do
    end if
    
    ! Checking to see if all polarization modes are there
    if (polarization) then
!       write(*,*) 'Number of polarization modules =', num_polmod
       do j = 1, num_polmod
          if (polmod(j)==-1) then
             write(*,*)'Data does not contain all polarization modules'
             stop
          end if
       end do
    end if
    
  deallocate(tempmod)  
  deallocate(true_tempmod)  
    
    
  end subroutine get_polmod
  
end module quiet_i2qu_mod
