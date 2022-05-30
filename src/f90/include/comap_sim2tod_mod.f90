module comap_sim2tod_mod
    use healpix_types
    use comap_scan_mod
    use comap_defs
    use comap_detector_mod
    use quiet_mpi_mod
    use quiet_hdf_mod
    use quiet_fft_mod
    use quiet_utils
    implicit none 
    
    type simulation_struct
       integer(i4b)         :: nx, ny, nx_hr, ny_hr, nfreq, nsb        ! Number of pixels and number of frequency channel
       real(dp)             :: boost                                   ! Factor to boost signal
       real(dp),     allocatable, dimension(:,:,:,:) :: simcube        ! (nx, ny, nfreq, nsb) data cube
       real(dp),     allocatable, dimension(:,:,:,:) :: simcube_hr ! (nx_hr, ny_hr, nfreq, nsb) data cube
       real(dp),     allocatable, dimension(:) :: x, y, edgex, edgey   ! (nx), (ny) RA, Dec grid
       real(dp),     allocatable, dimension(:, :, :, :, :, :) :: allcoeff    ! (4, 4, nx, ny, nsb, nfreq) interpolation coefficients 
       integer(i4b), allocatable, dimension(:, :)             :: freqidx ! (nsb, nfreq) array to store flipped frequency indcies
    end type simulation_struct
  
  contains
  
    subroutine read_sim_file(filename, data)
      implicit none
      character(len=*), intent(in)           :: filename
      type(simulation_struct)                :: data
      type(hdf_file)                         :: file
      integer(i4b)                           :: i, j, k, l
      integer(i4b), allocatable, dimension(:)       :: buffer_int
      real(dp),     allocatable, dimension(:)       :: buffer_1d
      real(sp),     allocatable, dimension(:,:,:)   :: buffer_3d
      real(sp),     allocatable, dimension(:,:,:,:) :: buffer_4d
      integer(i4b)                           :: ext4(4)

      call free_simulation_struct(data)
      call open_hdf_file(filename, file, "r")
      call get_size_hdf(file, "simulation", ext4)
      data%nx = ext4(1); data%ny = ext4(2); data%nfreq = ext4(3); data%nsb = ext4(4)
    
      allocate(data%simcube(data%nx, data%ny, data%nfreq, data%nsb))
      allocate(buffer_4d(data%nx, data%ny, data%nfreq, data%nsb))
      allocate(data%edgex(data%nx + 1))
      allocate(data%edgey(data%ny + 1))

      ! Read simulation cube

      call read_hdf(file, "simulation",     data%simcube)
      data%simcube = data%simcube / 1.d6 ! muK to K
      call read_hdf(file, "x",              data%edgex)
      call read_hdf(file, "y",              data%edgey)
   
      deallocate(buffer_4d)

      if (.not. allocated(data%freqidx)) allocate(data%freqidx(data%nsb, data%nfreq))

      do i = 1, data%nfreq
        data%freqidx(1, i) = data%nfreq + 1 - i! flipped sideband
        data%freqidx(3, i) = data%nfreq + 1 - i! flipped sideband

        data%freqidx(2, i) = i
        data%freqidx(4, i) = i
      end do

    end subroutine read_sim_file
  
    subroutine free_simulation_struct(data)
      implicit none
      type(simulation_struct) :: data
      if(allocated(data%simcube))   deallocate(data%simcube)
      if(allocated(data%x))   deallocate(data%x)
      if(allocated(data%y))   deallocate(data%y)
      if(allocated(data%edgex))   deallocate(data%edgex)
      if(allocated(data%edgey))   deallocate(data%edgey)
      if(allocated(data%allcoeff))   deallocate(data%allcoeff)
      if(allocated(data%freqidx)) deallocate(data%freqidx)
    end subroutine
  
  
  
    subroutine copy_simulation_struct(sim_in, sim_out)
      type(simulation_struct), intent(in)    :: sim_in
      type(simulation_struct), intent(inout) :: sim_out
      call free_simulation_struct(sim_out)
      
  !    write(*,*) "start"
      sim_out%nfreq = sim_in%nfreq
      sim_out%nsb = sim_in%nsb
      sim_out%nx = sim_in%nx
      sim_out%ny = sim_in%ny
    
      if(allocated(sim_in%simcube))        then
         allocate(sim_out%simcube(sim_in%nx, sim_in%ny, sim_in%nfreq, sim_in%nsb))
         sim_out%simcube = sim_in%simcube
      end if

      if(allocated(sim_in%simcube_hr))        then
        allocate(sim_out%simcube_hr(sim_in%nx_hr, sim_in%ny_hr, sim_in%nfreq, sim_in%nsb))
        sim_out%simcube_hr = sim_in%simcube_hr
     end if
     
    end subroutine copy_simulation_struct
    
  end module comap_sim2tod_mod
  
  
