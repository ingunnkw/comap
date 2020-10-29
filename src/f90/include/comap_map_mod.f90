module comap_map_mod
  use healpix_types
  use quiet_utils
  use quiet_hdf_mod
  implicit none

  real(dp), parameter :: MAP_BASE_PIXSIZE = 1.d0 ! Arcmin

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, nsb, ndet, ndet_tot, n_k, ntheta, nside, nsim, njk, nsplit ! 2^ntheta
     !real(dp)     :: x0, y0, f0, 
     real(dp)     :: dthetax, dthetay, df
     real(dp)     :: mean_az, mean_el, time(2), center(2)
     character(len=512) :: name

     character(len=4), allocatable, dimension(:) :: jk_def !(njk), jack0, jack1

     integer(i4b), allocatable, dimension(:)           :: feeds, jk_feed, succ_split
     real(dp),     allocatable, dimension(:)           :: x, y, k                               ! (n_x or n_y or n_k)
     real(dp),     allocatable, dimension(:,:)         :: freq                                  ! (nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:,:)   :: m, rms, dsum, div                     ! (n_x, n_y, nfreq, nsb, ndet)
     real(sp),     allocatable, dimension(:,:,:,:)     :: m_co, rms_co, dsum_co, div_co         ! (n_x, n_y, nfreq, nsb)
     real(sp),     allocatable, dimension(:,:,:,:,:,:) :: m_jk, rms_jk, dsum_jk, div_jk         ! (n_x, n_y, nfreq, nsb, ndet, 2*njk)
     real(sp),     allocatable, dimension(:,:,:,:,:,:) :: m_succ, rms_succ, dsum_succ, div_succ    ! (n_x, n_y, nfreq, nsb, ndet, 2**nsplit)
     real(sp),     allocatable, dimension(:,:,:,:,:)   :: m_jkco, rms_jkco, dsum_jkco, div_jkco ! (n_x, n_y, nfreq, nsb, 2*njk)
     real(sp),     allocatable, dimension(:,:,:,:,:,:) :: m_sim, rms_sim, dsum_sim, div_sim     ! (n_x, n_y, nfreq, nsb, ndet, nsim)
     integer(i4b), allocatable, dimension(:,:,:,:,:)   :: nhit, nhit_jkco                       ! (n_x, n_y, nfreq, nsb, ndet/2*njk)
     integer(i4b), allocatable, dimension(:,:,:,:)     :: nhit_co                               ! (n_x, n_y, nfreq, nsb)
     integer(i4b), allocatable, dimension(:,:,:,:,:,:) :: nhit_jk, nhit_succ                    ! (n_x, n_y, nfreq, nsb, ndet, 2*njk/2**nsplit)

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
    map2%njk     = map1%njk
    if (allocated(map1%m_jk)) then
       map2%jk_def    = map1%jk_def
       map2%m_jk      = map1%m_jk
       map2%rms_jk    = map1%rms_jk
       map2%nhit_jk   = map1%nhit_jk
       map2%m_jkco    = map1%m_jkco
       map2%rms_jkco  = map1%rms_jkco
       map2%nhit_jkco = map1%nhit_jkco
    end if
    
    if (allocated(map1%m_succ)) then
       map2%m_succ    = map1%m_succ
       map2%rms_succ  = map1%rms_succ
       map2%nhit_succ = map1%nhit_succ
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
    call write_hdf(file, "njk",          map%njk)
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
    if (map%njk > 0) then
       call create_hdf_group(file, "jackknives")
       call write_hdf(file, "jackknives/jk_def",  map%jk_def)
       call write_hdf(file, "jackknives/jk_feedmap",  map%jk_feed)
       nf = 1; nc = 1
       do i = 1, map%njk
          map_name = "jackknives/map_"  // map%jk_def(i)
          rms_name = "jackknives/rms_"  // map%jk_def(i)
          hit_name = "jackknives/nhit_" // map%jk_def(i)
          if (any(map%jk_feed == i)) then
             call write_hdf(file, trim(map_name), map%m_jk(:,:,:,:,:,2*nf-1:2*nf))
             call write_hdf(file, trim(rms_name), map%rms_jk(:,:,:,:,:,2*nf-1:2*nf))
             call write_hdf(file, trim(hit_name), map%nhit_jk(:,:,:,:,:,2*nf-1:2*nf))
             nf = nf + 1 
          else
             call write_hdf(file, trim(map_name), map%m_jkco(:,:,:,:,2*nc-1:2*nc))
             call write_hdf(file, trim(rms_name), map%rms_jkco(:,:,:,:,2*nc-1:2*nc))
             call write_hdf(file, trim(hit_name), map%nhit_jkco(:,:,:,:,2*nc-1:2*nc))
             nc = nc + 1
          end if
       end do

       if (map%nsplit > 0) then
          call write_hdf(file, "jackknives/map_succ", map%m_succ)
          call write_hdf(file, "jackknives/rms_succ", map%rms_succ)
          call write_hdf(file, "jackknives/nhit_succ", map%nhit_succ)
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
    call read_hdf(file, "njk",          map%njk)

    ! Read co-added over feeds
    call read_hdf(file, "map_coadd", map%m_co)
    call read_hdf(file, "rms_coadd", map%rms_co)
    call read_hdf(file, "nhit_coadd", map%nhit_co)
    
    ! Read jackknives
    if (map%njk > 0) then
       call get_size_hdf(file, "jackknives", ext)
       n_feed = ext(1); n_coadd = map%njk - n_feed
       allocate(map%jk_def(map%njk), map%jk_feed(n_feed))       
       call read_hdf(file, "jackknives/jk_def",  map%jk_def)
       call read_hdf(file, "jackknives/jk_feedmap",  map%jk_feed)
       allocate(map%m_jkco(nx,ny,nfreq,nsb,2*n_coadd), &
            & map%rms_jkco(nx,ny,nfreq,nsb,2*n_coadd), &
            & map%nhit_jkco(nx,ny,nfreq,nsb,2*n_coadd), &
            & map%m_jk(nx,ny,nfreq,nsb,ndet,2*n_feed), &
            & map%rms_jk(nx,ny,nfreq,nsb,ndet,2*n_feed), &
            & map%nhit_jk(nx,ny,nfreq,nsb,ndet,2*n_feed))
       nf = 1; nc = 1
       do i = 1, map%njk
          map_name = "jackknives/map_"  // map%jk_def(i)
          rms_name = "jackknives/rms_"  // map%jk_def(i)
          hit_name = "jackknives/nhit_" // map%jk_def(i)
          if (any(map%jk_feed == i)) then
             call read_hdf(file, trim(map_name), map%m_jk(:,:,:,:,:,2*nf-1:2*nf))
             call read_hdf(file, trim(rms_name), map%rms_jk(:,:,:,:,:,2*nf-1:2*nf))
             call read_hdf(file, trim(hit_name), map%nhit_jk(:,:,:,:,:,2*nf-1:2*nf))
             nf = nf + 1
          else
             call read_hdf(file, trim(map_name), map%m_jkco(:,:,:,:,2*nc-1:2*nc))
             call read_hdf(file, trim(rms_name), map%rms_jkco(:,:,:,:,2*nc-1:2*nc))
             call read_hdf(file, trim(hit_name), map%nhit_jkco(:,:,:,:,2*nc-1:2*nc))
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
    map%m_jk      = 0.0
    map%rms_jk    = 0.0
    map%dsum_jk   = 0.0
    map%div_jk    = 0.0
    map%nhit_jk   = 0
    map%m_jkco    = 0.0
    map%rms_jkco  = 0.0
    map%dsum_jkco = 0.0
    map%div_jkco  = 0.0
    map%nhit_jkco = 0

    ! Successive splits
    map%m_succ    = 0.0
    map%rms_succ  = 0.0
    map%dsum_succ = 0.0
    map%div_succ  = 0.0
    map%nhit_succ = 0

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
    if (allocated(map%m_jk))    deallocate(map%m_jk)
    if (allocated(map%rms_jk))  deallocate(map%rms_jk)
    if (allocated(map%nhit_jk)) deallocate(map%nhit_jk)
    if (allocated(map%div_jk))  deallocate(map%div_jk)
    if (allocated(map%dsum_jk)) deallocate(map%dsum_jk)
    if (allocated(map%m_jkco))    deallocate(map%m_jkco)
    if (allocated(map%rms_jkco))  deallocate(map%rms_jkco)
    if (allocated(map%nhit_jkco)) deallocate(map%nhit_jkco)
    if (allocated(map%div_jkco))  deallocate(map%div_jkco)
    if (allocated(map%dsum_jkco)) deallocate(map%dsum_jkco)
    
    ! successive splits
    if (allocated(map%m_succ))      deallocate(map%m_succ)
    if (allocated(map%rms_succ))    deallocate(map%rms_succ)
    if (allocated(map%nhit_succ))   deallocate(map%nhit_succ)
    if (allocated(map%div_succ))    deallocate(map%div_succ)
    if (allocated(map%dsum_succ))   deallocate(map%dsum_succ)
    
    
    ! Simulated data 
    if (allocated(map%m_sim))    deallocate(map%m_sim)
    if (allocated(map%rms_sim))  deallocate(map%rms_sim) 
    if (allocated(map%dsum_sim)) deallocate(map%dsum_sim) 
    if (allocated(map%div_sim))  deallocate(map%div_sim) 
 
  end subroutine free_map_type

end module comap_map_mod
