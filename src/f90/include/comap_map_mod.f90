module comap_map_mod
  use healpix_types
  use quiet_utils
  use quiet_hdf_mod
  implicit none

  real(dp), parameter :: MAP_BASE_PIXSIZE = 1.d0 ! Arcmin

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, nsb, ndet, ndet_tot, n_k, ntheta, nside, nsim, nsplit, nmultisplit ! 2^ntheta
     !real(dp)     :: x0, y0, f0, 
     real(dp)     :: dthetax, dthetay, df
     real(dp)     :: mean_az, mean_el, time(2), center(2)
     character(len=512) :: name

     character(len=4), allocatable, dimension(:) :: split_def !(nsplit), jack0, jack1

     integer(i4b), allocatable, dimension(:)           :: feeds, split_feed, split_split
     real(dp),     allocatable, dimension(:)           :: x, y, k                               ! (n_x or n_y or n_k)
     real(dp),     allocatable, dimension(:,:)         :: freq                                  ! (nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:,:)   :: m, rms, dsum, div                     ! (n_x, n_y, nfreq, nsb, ndet)
     real(sp),     allocatable, dimension(:,:,:,:)     :: m_co, rms_co, dsum_co, div_co         ! (n_x, n_y, nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:,:,:) :: m_split, rms_split, dsum_split, div_split         ! (n_x, n_y, nfreq, nsb, ndet, 2*nsplit)
     real(sp),     allocatable, dimension(:,:,:,:,:,:) :: m_multisplit, rms_multisplit, dsum_multisplit, div_multisplit    ! (n_x, n_y, nfreq, nsb, ndet, 2**nmultisplit)
     real(sp),     allocatable, dimension(:,:,:,:,:)   :: m_splitco, rms_splitco, dsum_splitco, div_splitco ! (n_x, n_y, nfreq, nsb, 2*nsplit)
     real(sp),     allocatable, dimension(:,:,:,:,:,:) :: m_sim, rms_sim, dsum_sim, div_sim     ! (n_x, n_y, nfreq, nsb, ndet, nsim)
     integer(i4b), allocatable, dimension(:,:,:,:,:)   :: nhit, nhit_splitco                       ! (n_x, n_y, nfreq, nsb, ndet/2*nsplit)
     integer(i4b), allocatable, dimension(:,:,:,:)     :: nhit_co                               ! (n_x, n_y, nfreq, nsb)
     integer(i4b), allocatable, dimension(:,:,:,:,:,:) :: nhit_split, nhit_multisplit                    ! (n_x, n_y, nfreq, nsb, ndet, 2*nsplit/2**nmultisplit)

  end type map_type


contains

  subroutine initialize_map_mod()
    implicit none

  end subroutine initialize_map_mod


  subroutine copy_map(map1, map2)
    implicit none
    type(map_type), intent(in)    :: map1
    type(map_type), intent(inout) :: map2

    call free_map_type(map2)
    map2%n_x     = map1%n_x
    map2%n_y     = map1%n_y
    map2%x       = map1%x
    map2%y       = map1%y
    map2%m       = map1%m
    map2%rms     = map1%rms
    map2%nhit    = map1%nhit
    map2%m_co    = map1%m_co
    map2%rms_co  = map1%rms_co
    map2%nhit_co = map1%nhit_co
    map2%nsplit     = map1%nsplit
    if (allocated(map1%m_split)) then
       map2%split_def    = map1%split_def
       map2%m_split      = map1%m_split
       map2%rms_split    = map1%rms_split
       map2%nhit_split   = map1%nhit_split
       map2%m_splitco    = map1%m_splitco
       map2%rms_splitco  = map1%rms_splitco
       map2%nhit_splitco = map1%nhit_splitco
    end if
    
    if (allocated(map1%m_multisplit)) then
       map2%m_multisplit    = map1%m_multisplit
       map2%rms_multisplit  = map1%rms_multisplit
       map2%nhit_multisplit = map1%nhit_multisplit
    end if


    map2%nsim   = map1%nsim
    if (allocated(map1%m_sim)) then
       map2%m_sim   = map1%m_sim
       map2%rms_sim = map1%rms_sim
    end if
    map2%freq    = map1%freq
    map2%time    = map1%time
    map2%mean_az = map1%mean_az
    map2%mean_el = map1%mean_el
    map2%feeds   = map1%feeds
    map2%center  = map1%center
    map2%nside   = map1%nside

  end subroutine copy_map


  ! Writes an h5 file with maps/rms/nhit
  subroutine output_map_h5(filename, map, det, sb)
    implicit none
    character(len=*), intent(in)    :: filename
    type(map_type),   intent(inout) :: map
    integer(i4b), optional :: det, sb
    integer(i4b)       :: i, nf, nc
    
    character(len=120) :: map_name, rms_name, hit_name  
    type(hdf_file)     :: file
    
    call open_hdf_file(trim(filename), file, "w")
    call write_hdf(file, "n_x",          map%n_x)
    call write_hdf(file, "n_y",          map%n_y)
    call write_hdf(file, "x",            map%x)
    call write_hdf(file, "y",            map%y)
    call write_hdf(file, "patch_center", map%center)
    call write_hdf(file, "nside",        map%nside)
    call write_hdf(file, "freq",         map%freq)
    call write_hdf(file, "mean_az",      map%mean_az)
    call write_hdf(file, "mean_el",      map%mean_el)
    call write_hdf(file, "time",         map%time)
    call write_hdf(file, "feeds",        map%feeds)
    call write_hdf(file, "nsplit",          map%nsplit)
    call write_hdf(file, "nsim",         map%nsim)
    if (present(det)) then
       call write_hdf(file, "map",  map%m(:,:,:,sb:sb,det:det))
       call write_hdf(file, "rms",  map%rms(:,:,:,sb:sb,det:det))
       call write_hdf(file, "nhit", map%nhit(:,:,:,sb:sb,det:det))
    else
       call write_hdf(file, "map",  map%m)
       call write_hdf(file, "rms",  map%rms)
       call write_hdf(file, "nhit", map%nhit)
       ! Co-added over feeds
       call write_hdf(file, "map_coadd",  map%m_co)
       call write_hdf(file, "rms_coadd",  map%rms_co)
       call write_hdf(file, "nhit_coadd", map%nhit_co)
    end if
    if (map%nsim > 0) then
       call write_hdf(file, "map_sim", map%m_sim)
       call write_hdf(file, "rms_sim", map%rms_sim)
    end if
    if (map%nsplit > 0) then
       call create_hdf_group(file, "jackknives")
       call write_hdf(file, "jackknives/split_def",  map%split_def)
       call write_hdf(file, "jackknives/split_feedmap",  map%split_feed)
       nf = 1; nc = 1
       do i = 1, map%nsplit
          map_name = "jackknives/map_"  // map%split_def(i)
          rms_name = "jackknives/rms_"  // map%split_def(i)
          hit_name = "jackknives/nhit_" // map%split_def(i)
          if (any(map%split_feed == i)) then
             call write_hdf(file, trim(map_name), map%m_split(:,:,:,:,:,2*nf-1:2*nf))
             call write_hdf(file, trim(rms_name), map%rms_split(:,:,:,:,:,2*nf-1:2*nf))
             call write_hdf(file, trim(hit_name), map%nhit_split(:,:,:,:,:,2*nf-1:2*nf))
             nf = nf + 1 
          else
             call write_hdf(file, trim(map_name), map%m_splitco(:,:,:,:,2*nc-1:2*nc))
             call write_hdf(file, trim(rms_name), map%rms_splitco(:,:,:,:,2*nc-1:2*nc))
             call write_hdf(file, trim(hit_name), map%nhit_splitco(:,:,:,:,2*nc-1:2*nc))
             nc = nc + 1
          end if
       end do

       if (map%nmultisplit > 0) then
          call write_hdf(file, "jackknives/map_split", map%m_multisplit)
          call write_hdf(file, "jackknives/rms_split", map%rms_multisplit)
          call write_hdf(file, "jackknives/nhit_split", map%nhit_multisplit)
       end if
    end if

    call close_hdf_file(file)

    !call free_map_type(map)

  end subroutine output_map_h5


  subroutine output_submap_sim_h5(filename, map, sim)
    implicit none
    character(len=*), intent(in) :: filename
    type(map_type),   intent(in) :: map
    integer(i4b),     intent(in) :: sim

    type(hdf_file)     :: file
 

    call open_hdf_file(trim(filename), file, "w")

    ! For simulated data 
    call write_hdf(file, "map_sim", map%m_sim) 
    call write_hdf(file, "rms_sim", map%rms_sim)
    call write_hdf(file, "sim", sim)  
    ! call write_hdf(file, "map_sim_beam", sum(map%m_sim, dim=5))
    ! call write_hdf(file, "rms_sim_beam", sum(map%rms_sim, dim=5))

    call close_hdf_file(file)

  end subroutine output_submap_sim_h5


  subroutine output_submap_h5(filename, map)
    implicit none
    character(len=*), intent(in) :: filename
    type(map_type),   intent(in) :: map

    type(hdf_file)     :: file
    
    call open_hdf_file(trim(filename), file, "w")
    call write_hdf(file, "n_x", map%n_x)
    call write_hdf(file, "n_y", map%n_y)
    call write_hdf(file, "x",   map%x)
    call write_hdf(file, "y",   map%y)
    call write_hdf(file, "map", map%m)!(:,:,1,1,1))
    call write_hdf(file, "rms", map%rms)!(:,:,1,1,1))
    call write_hdf(file, "nhit", map%nhit)!(:,:,1,1,1))
    call write_hdf(file, "freq", map%freq)
    call write_hdf(file, "mean_az", map%mean_az)
    call write_hdf(file, "mean_el", map%mean_el)
    call write_hdf(file, "time", map%time)
    call write_hdf(file, "feeds", map%feeds)
    call write_hdf(file, "patch_center", map%center)
    call write_hdf(file, "nside", map%nside)
    call close_hdf_file(file)

  end subroutine output_submap_h5


  ! Reads an h5 file
  subroutine read_map_h5(filename,map)
    implicit none
    character(len=*), intent(in)  :: filename
    type(map_type),   intent(out) :: map

    type(hdf_file) :: file
    integer(i4b)   :: i, nx, ny, nfreq, nsb, ndet, ext(7)
    integer(i4b)   :: nf, nc, n_feed, n_coadd
    character(len=120) :: map_name, rms_name, hit_name

    call free_map_type(map)

    call open_hdf_file(trim(filename), file, "r")

    call get_size_hdf(file, "map", ext)
    nx = ext(1); ny = ext(2); nfreq = ext(3); nsb = ext(4); ndet = ext(5)

    allocate(map%x(nx), map%y(ny))
    allocate(map%m(nx,ny,nfreq,nsb,ndet), map%rms(nx,ny,nfreq,nsb,ndet), &
         & map%nhit(nx,ny,nfreq, nsb,ndet))
    allocate(map%m_co(nx,ny,nfreq,nsb), map%rms_co(nx,ny,nfreq,nsb), &
         & map%nhit_co(nx,ny,nfreq,nsb))
    allocate(map%freq(nfreq,nsb), map%feeds(ndet))
    

    call read_hdf(file, "n_x",          map%n_x)
    call read_hdf(file, "n_y",          map%n_y)
    call read_hdf(file, "x",            map%x)
    call read_hdf(file, "y",            map%y)
    call read_hdf(file, "map",          map%m)
    call read_hdf(file, "rms",          map%rms)
    call read_hdf(file, "nhit",         map%nhit)
    call read_hdf(file, "freq",         map%freq)
    call read_hdf(file, "time",         map%time)
    call read_hdf(file, "mean_az",      map%mean_az)
    call read_hdf(file, "mean_el",      map%mean_el)
    call read_hdf(file, "feeds",        map%feeds)
    call read_hdf(file, "patch_center", map%center)
    call read_hdf(file, "nside",        map%nside)
    call read_hdf(file, "nsim",         map%nsim)
    call read_hdf(file, "njk",          map%nsplit)

    ! Read co-added over feeds
    call read_hdf(file, "map_coadd", map%m_co)
    call read_hdf(file, "rms_coadd", map%rms_co)
    call read_hdf(file, "nhit_coadd", map%nhit_co)
    
    ! Read jackknives
    if (map%nsplit > 0) then
       call get_size_hdf(file, "jackknives", ext)
       n_feed = ext(1); n_coadd = map%nsplit - n_feed
       allocate(map%split_def(map%nsplit), map%split_feed(n_feed))       
       call read_hdf(file, "jackknives/jk_def",  map%split_def)
       call read_hdf(file, "jackknives/jk_feedmap",  map%split_feed)
       allocate(map%m_splitco(nx,ny,nfreq,nsb,2*n_coadd), &
            & map%rms_splitco(nx,ny,nfreq,nsb,2*n_coadd), &
            & map%nhit_splitco(nx,ny,nfreq,nsb,2*n_coadd), &
            & map%m_split(nx,ny,nfreq,nsb,ndet,2*n_feed), &
            & map%rms_split(nx,ny,nfreq,nsb,ndet,2*n_feed), &
            & map%nhit_split(nx,ny,nfreq,nsb,ndet,2*n_feed))
       nf = 1; nc = 1
       do i = 1, map%nsplit
          map_name = "jackknives/map_"  // map%split_def(i)
          rms_name = "jackknives/rms_"  // map%split_def(i)
          hit_name = "jackknives/nhit_" // map%split_def(i)
          if (any(map%split_feed == i)) then
             call read_hdf(file, trim(map_name), map%m_split(:,:,:,:,:,2*nf-1:2*nf))
             call read_hdf(file, trim(rms_name), map%rms_split(:,:,:,:,:,2*nf-1:2*nf))
             call read_hdf(file, trim(hit_name), map%nhit_split(:,:,:,:,:,2*nf-1:2*nf))
             nf = nf + 1
          else
             call read_hdf(file, trim(map_name), map%m_splitco(:,:,:,:,2*nc-1:2*nc))
             call read_hdf(file, trim(rms_name), map%rms_splitco(:,:,:,:,2*nc-1:2*nc))
             call read_hdf(file, trim(hit_name), map%nhit_splitco(:,:,:,:,2*nc-1:2*nc))
             nc = nc + 1
          end if
       end do
    end if

    ! Read simulated data
    if (map%nsim > 0) then
       allocate(map%m_sim(nx,ny,nfreq,nsb,ndet, map%nsim), map%rms_sim(nx,ny,nfreq,nsb,ndet, map%nsim))
       call read_hdf(file, "map_sim", map%m_sim)
       call read_hdf(file, "rms_sim", map%rms_sim)
    end if


    call close_hdf_file(file)

  end subroutine read_map_h5


!   ! Creates a .dat file with the maps/rms/nhit for each frequency
!   subroutine output_maps(prefix, map)
!     implicit none
!     character(len=*), intent(in)    :: prefix
!     type(map_type),   intent(inout) :: map

!     integer(i4b)       :: i, j, k, unit
!     character(len=4)   :: itext
!     character(len=512) :: filename

!     unit = getlun()
!     do i = 6, 6
!     !do i = 1, map%nfreq
!        call int2string(i,itext)
!        filename = trim(prefix)//'_freq'//itext//'_map.dat'
!        open(unit, file=trim(filename), recl=100000)
!        write(unit,*) '# n_x = ', map%n_x
!        write(unit,*) '# n_y = ', map%n_y
!        write(unit,*) '# x   = ', real(map%x,sp)
!        write(unit,*) '# y   = ', real(map%y,sp)
!        do j = 1, map%n_x
!           do k = 1, map%n_y
!              write(unit,fmt='(e16.8)',advance='no') map%m(j,k,i,:)
!           end do
!           write(unit,*)
!        end do
!        close(unit)
!     end do

!     unit = getlun()
!     do i = 6, 6
! !    do i = 1, map%nfreq
!        call int2string(i,itext)
!        filename = trim(prefix)//'_freq'//itext//'_rms.dat'
!        open(unit, file=trim(filename), recl=100000)
!        write(unit,*) '# n_x = ', map%n_x
!        write(unit,*) '# n_y = ', map%n_y
!        write(unit,*) '# x   = ', real(map%x,sp)
!        write(unit,*) '# y   = ', real(map%y,sp)
!        do j = 1, map%n_x
!           do k = 1, map%n_y
!              write(unit,fmt='(e16.8)',advance='no') map%rms(j,k,i,:)
!           end do
!           write(unit,*)
!        end do
!        close(unit)
!     end do

!     unit = getlun()
!     do i = 6, 6
! !    do i = 1, map%nfreq
!        call int2string(i,itext)
!        filename = trim(prefix)//'_freq'//itext//'_nhit.dat'
!        open(unit, file=trim(filename), recl=100000)
!        write(unit,*) '# n_x = ', map%n_x
!        write(unit,*) '# n_y = ', map%n_y
!        write(unit,*) '# x   = ', real(map%x,sp)
!        write(unit,*) '# y   = ', real(map%y,sp)
!        do j = 1, map%n_x
!           do k = 1, map%n_y
!              write(unit,fmt='(e16.8)',advance='no') map%nhit(j,k,i,:)
!           end do
!           write(unit,*)
!        end do
!        close(unit)
!     end do

!     call free_map_type(map)

!   end subroutine output_maps

  subroutine nullify_map_type(map)
    implicit none
    type(map_type), intent(inout) :: map

    map%m       = 0.0
    map%rms     = 0.0
    map%dsum    = 0.0
    map%div     = 0.0
    map%nhit    = 0

    ! Co-added over feeds
    map%m_co    = 0.0
    map%rms_co  = 0.0
    map%dsum_co = 0.0
    map%div_co  = 0.0
    map%nhit_co = 0

    ! Jackknives
    map%m_split      = 0.0
    map%rms_split    = 0.0
    map%dsum_split   = 0.0
    map%div_split    = 0.0
    map%nhit_split   = 0
    map%m_splitco    = 0.0
    map%rms_splitco  = 0.0
    map%dsum_splitco = 0.0
    map%div_splitco  = 0.0
    map%nhit_splitco = 0

    ! Successive splits
    map%m_multisplit    = 0.0
    map%rms_multisplit  = 0.0
    map%dsum_multisplit = 0.0
    map%div_multisplit  = 0.0
    map%nhit_multisplit = 0

    ! Simulated data
    map%m_sim    = 0.0
    map%rms_sim  = 0.0
    map%dsum_sim = 0.0
    map%div_sim  = 0.0 

  end subroutine nullify_map_type


  subroutine free_map_type(map)
    implicit none
    type(map_type), intent(inout) :: map

    if (allocated(map%x))       deallocate(map%x) 
    if (allocated(map%y))       deallocate(map%y)
    if (allocated(map%freq))    deallocate(map%freq)
    if (allocated(map%k))       deallocate(map%k)
    if (allocated(map%m))       deallocate(map%m)
    if (allocated(map%rms))     deallocate(map%rms)
    if (allocated(map%dsum))    deallocate(map%dsum)
    if (allocated(map%nhit))    deallocate(map%nhit)
    if (allocated(map%div))     deallocate(map%div)
    if (allocated(map%feeds))   deallocate(map%feeds)

    ! Co-added over feeds
    if (allocated(map%m_co))    deallocate(map%m_co)
    if (allocated(map%rms_co))  deallocate(map%rms_co)
    if (allocated(map%nhit_co)) deallocate(map%nhit_co)
    if (allocated(map%div_co))  deallocate(map%div_co)
    if (allocated(map%dsum_co)) deallocate(map%dsum_co)

    ! Jackknives
    if (allocated(map%m_split))    deallocate(map%m_split)
    if (allocated(map%rms_split))  deallocate(map%rms_split)
    if (allocated(map%nhit_split)) deallocate(map%nhit_split)
    if (allocated(map%div_split))  deallocate(map%div_split)
    if (allocated(map%dsum_split)) deallocate(map%dsum_split)
    if (allocated(map%m_splitco))    deallocate(map%m_splitco)
    if (allocated(map%rms_splitco))  deallocate(map%rms_splitco)
    if (allocated(map%nhit_splitco)) deallocate(map%nhit_splitco)
    if (allocated(map%div_splitco))  deallocate(map%div_splitco)
    if (allocated(map%dsum_splitco)) deallocate(map%dsum_splitco)
    
    ! successive splits
    if (allocated(map%m_multisplit))      deallocate(map%m_multisplit)
    if (allocated(map%rms_multisplit))    deallocate(map%rms_multisplit)
    if (allocated(map%nhit_multisplit))   deallocate(map%nhit_multisplit)
    if (allocated(map%div_multisplit))    deallocate(map%div_multisplit)
    if (allocated(map%dsum_multisplit))   deallocate(map%dsum_multisplit)
    
    
    ! Simulated data 
    if (allocated(map%m_sim))    deallocate(map%m_sim)
    if (allocated(map%rms_sim))  deallocate(map%rms_sim) 
    if (allocated(map%dsum_sim)) deallocate(map%dsum_sim) 
    if (allocated(map%div_sim))  deallocate(map%div_sim) 
 
  end subroutine free_map_type

end module comap_map_mod
