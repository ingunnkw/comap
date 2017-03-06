module comap_Lx_mod
  use healpix_types
  use quiet_defs
  use quiet_hdf_mod
  use l1_read_mod
  use quiet_fft_mod
  use quiet_utils
  implicit none

  type Lx_struct
     !! The level1 and level2 part, present in all files.
     integer(i4b)                                :: decimation
     real(dp)                                    :: samprate
     type(hk_struct)                             :: hk
     real(dp), allocatable, dimension(:)         :: time
     real(sp), allocatable, dimension(:,:,:)     :: tod         ! (time, freq, detector)
     real(sp), allocatable, dimension(:,:,:,:)   :: tod_l1      ! (time, freq, sideband, detector)
     real(sp), allocatable, dimension(:,:)       :: orig_point  ! Hor; (az/el/dk, time)

     !! The level3 part, which is only present in level3-files
     integer(i4b)                                 :: coord_sys
     real(dp)                                     :: scanfreq(2), pixsize, mapsize(2,2)
     real(sp),     allocatable, dimension(:,:,:)  :: point        ! Gal; (phi/theta/psi,time,mod)

     real(dp),     allocatable, dimension(:)        :: time_gain            ! (time)
     real(sp),     allocatable, dimension(:,:,:)    :: gain                 ! (time, freq, detector)
     real(dp),     allocatable, dimension(:,:)      :: sigma0, alpha, fknee ! (freq, detector)
     real(dp),     allocatable, dimension(:,:,:,:)  :: corr                 ! (freq, freq, detector, detector)
     real(dp)                                       :: stats(STAT_NUM)
     real(dp),     allocatable, dimension(:,:,:)    :: det_stats            ! (freq, detector, stat)
     real(dp),     allocatable, dimension(:,:,:)    :: filter_par           ! (freq, detector, param)

  end type Lx_struct

contains

  subroutine read_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: nsamp, nfreq, ndet, npoint, ext(7)
    call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "tod", ext)
    nsamp = ext(1); nfreq = ext(2); ndet = ext(3)
    call get_size_hdf(file, "orig_point", ext)
    npoint = ext(1)
    allocate(data%time(nsamp), data%tod(nsamp,nfreq,ndet))
    allocate(data%orig_point(npoint,nsamp))
    call read_hdf(file, "decimation", data%decimation)
    call read_hdf(file, "samprate",   data%samprate)
    call read_hdf(file, "time",       data%time)
    call read_hdf(file, "tod",        data%tod)
    call read_hdf(file, "orig_point", data%orig_point)
    call read_hk_hdf(file, data%hk)
    call close_hdf_file(file)
  end subroutine read_L2_file

! Version edited for Planck
  subroutine read_l2_planck(filename, data, parfile)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=512)           :: parfile
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: nsamp, ndi, npoint, ext(7), i
    real(dp), allocatable        :: time(:), tod(:,:), point(:,:,:)
    real(dp), allocatable        :: time_full(:), tod_full(:,:), point_full(:,:,:), tmp(:,:)
    call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "tod", ext) !tod : Dataset {2807349, 4} according to h5ls
    nsamp = ext(2); ndi = ext(1)
    call get_size_hdf(file, "point", ext)
    npoint = ext(1) 
    call get_parameter(0, parfile, 'DECIMATION_PLANCK', par_int=data%decimation)
    call get_parameter(0, parfile, 'SAMPRATE_ORIG_PLANCK', par_dp=data%samprate)

    ! Reading undecimated data
    allocate(time_full(nsamp), tod_full(nsamp,ndi), point_full(npoint,nsamp,0:ndi-1))
    call read_hdf(file, "time",   time_full)
    call read_hdf(file, "tod",     tod_full)
    call read_hdf(file, "point", point_full)
!    call read_hk_hdf(file, data%hk) ! no housekeeping as of yet
    call close_hdf_file(file)

    ! Scaling time to seconds rather than Planck OBT. Do we need it in MJD?
    time_full = time_full/2.0**16
    ! Ordering in l2files is theta, phi, psi, while for QUIET we have used phi,theta,psi. 
    ! There is probably a better way of doing this?
    allocate(tmp(nsamp,ndi))
    tmp = point_full(1,:,:)
    point_full(1,:,:) = point_full(2,:,:)
    point_full(2,:,:) = tmp
    deallocate(tmp)
!!$
!!$open(42,file="tod_undecimated.txt")
!!$do i=1,nsamp
!!$   write(42,*) time_full(i), tod_full(i,1)
!!$end do
!!$close(42)

    ! Decimating data
    if (data%decimation>0) then
       nsamp = nsamp/data%decimation !integer division - rounding down (nsamp is not in l3file, btw)
       data%samprate = data%samprate/data%decimation ! srate is real
    end if

    allocate(time(nsamp), tod(nsamp,ndi), point(npoint,nsamp,0:ndi-1))
    call decimate(time, time_full, tod, tod_full, point, point_full, data%decimation)
    deallocate(time_full, tod_full, point_full)
!!$write(*,*) "after decimation", nsamp
!!$open(42,file="decim_303.txt")
!!$do i=1,nsamp
!!$   write(42,'(3g15.7)') point(:,i,0)
!!$end do
!!$close(42)

    ! Shortening arrays to ensure periodicity (for pretty FFT's)
    ! Choosing the phi component for this
    call chop_pointing_data(point(2,:,0), nsamp)
!!$write(*,*) "after chop, should be shorter than former", nsamp
!!$open(42,file="chopped_303.txt")
!!$do i=1,nsamp
!!$   write(42,'(3g15.7)') point(:,i,0)
!!$end do
!!$close(42)

    ! Storing in data structure
    allocate(data%time(nsamp), data%tod(nsamp,ndi))! Now nsamp is updated
    allocate(data%orig_point(npoint,nsamp), data%point(npoint,nsamp,0:ndi-1))
    data%time = time(1:nsamp)
    data%tod = tod(1:nsamp,:)
    data%point = point(:,1:nsamp,:)
    data%orig_point = data%point(:,:,0) ! calc_scanfreq needs orig_point - just use the first diode
    deallocate(time, tod, point)

  end subroutine read_l2_planck

! where should this sub logically be?
  subroutine chop_pointing_data(array, i)
    implicit none
    real(dp), dimension(:), intent(in)    :: array
    integer(i4b),           intent(inout) :: i
    real(dp), allocatable                 :: diff(:)
    real(dp)                              :: delta
    allocate(diff(i))
    diff = abs(array - array(1))
    delta = diff(2)
    do while (diff(i) > delta)
       i = i-1
    end do
    deallocate(diff)
  end subroutine chop_pointing_data


  ! Where should this sub logically be?
  subroutine decimate(time, time_full, tod, tod_full, point, point_full, dec)
    implicit none
    integer(i4b)                                      :: n, ndi, npt, dec, i, j, k, ind
    real(dp), dimension(:),      intent(inout)        :: time, time_full
    real(dp), dimension(:,:),    intent(inout)        :: tod, tod_full
    real(dp), dimension(:,:,0:), intent(inout)        :: point, point_full

    if (dec>0) then
       n   = size(tod,1) 
       ndi = size(tod,2)
       npt = size(point,1)

       ! Averaging over every (dec) elements of the time dimension of time, tod and pointing arrays
       ind = 1
       do i=1,n
          time(i) = mean(time_full(ind:ind+dec-1))
          do j=1,ndi
             tod(i,j) = mean(tod_full(ind:ind+dec-1,j))
             do k=1,npt
                ! do I need to make the angles safe?
                point(k,i,j-1) = mean(point_full(k,ind:ind+dec-1,j-1))
             end do
          end do
          ind = ind + dec
       end do

    else
       time = time_full
       tod = tod_full
       point = point_full
    end if
  end subroutine decimate


  subroutine read_l3_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b) :: npix, nsamp, ndi, nmod, nobj, nbin, ext(7)
    call free_lx_struct(data)
    call read_l2_file(filename, data)
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "nside",     data%nside)
    call read_hdf(file, "scanfreq",  data%scanfreq)
    ! Read pointing
    call get_size_hdf(file, "point", ext)
    nsamp = ext(2); nmod = ext(3)
    allocate(data%point(3,nsamp,nmod))
    call read_hdf(file, "point",     data%point)
    call read_hdf(file, "coord_sys", data%coord_sys)
    ! Read object-relative pointing
    call get_size_hdf(file, "point_objrel", ext)
    nsamp = ext(2); nobj = ext(3)
    allocate(data%point_objrel(3,nsamp,nobj))
    call read_hdf(file, "point_objrel", data%point_objrel)
    ! Read pixel list
    call get_size_hdf(file, "pixels", ext)
    npix = ext(1)
    allocate(data%pixels(npix))
    call read_hdf(file, "pixels", data%pixels)
    ! Read gain
    call get_size_hdf(file, "time_gain", ext)
    nsamp = ext(1)
    allocate(data%time_gain(nsamp))
    call read_hdf(file, "time_gain", data%time_gain)
    call get_size_hdf(file, "gain", ext)
    nsamp = ext(1); ndi = ext(2)
    allocate(data%gain(nsamp,ndi))
    call read_hdf(file, "gain", data%gain)
    ! Read noise parameters
    call get_size_hdf(file, "corr", ext)
    nbin = ext(1)
    allocate(data%sigma0(ndi),data%alpha(ndi),data%fknee(ndi))
    allocate(data%corr(nbin,ndi,ndi),data%corr_freqs(nbin))
    allocate(data%diode_stats(ndi,NUM_DIODE_STATS), data%filter_par(ndi,NUM_FILTER_PAR))
    call read_hdf(file, "sigma0", data%sigma0)
    call read_hdf(file, "alpha",  data%alpha)
    call read_hdf(file, "fknee",  data%fknee)
    call read_hdf(file, "corr",   data%corr)
    call read_hdf(file, "corr_freqs", data%corr_freqs)
    ! Read APEX data
    call read_hk_elem(file, "apex", data%hk%apex)
    ! Read stats
    call read_hdf(file, "stats",       data%stats)
    call read_hdf(file, "diode_stats", data%diode_stats)
    ! Read filter parameters
    call read_hdf(file, "filter_par",  data%filter_par)
    call close_hdf_file(file)
  end subroutine read_L3_file

  subroutine read_hk_hdf(file, hk)
    implicit none
    type(hdf_file)  :: file
    type(hk_struct) :: hk
    call read_hk_elem(file, "bias", hk%bias)
    call read_hk_elem(file, "cryo", hk%cryo)
    call read_hk_elem(file, "encl", hk%encl)
    call read_hk_elem(file, "peri", hk%peri)
  end subroutine

  subroutine read_hk_elem(file, name, hkt)
    implicit none
    type(hdf_file)   :: file
    type(hk_type)    :: hkt
    character(len=*) :: name
    integer(i4b)     :: nsamp, ntype, ext(7)
    call get_size_hdf(file, name // "/value", ext)
    nsamp = ext(1); ntype = ext(2)
    allocate(hkt%name(ntype), hkt%time(nsamp), hkt%value(nsamp,ntype))
    hkt%n   = ntype
    hkt%n_t = nsamp
    call read_hdf(file, trim(name) // "/time",  hkt%time)
    call read_hdf(file, trim(name) // "/value", hkt%value)
  end subroutine

  subroutine free_lx_struct(data)
    implicit none
    type(lx_struct) :: data
    if(allocated(data%time))        deallocate(data%time)
    if(allocated(data%tod))         deallocate(data%tod)
    if(allocated(data%tp))          deallocate(data%tp)
    if(allocated(data%orig_point))  deallocate(data%orig_point)
    call deallocate_hk_struct(data%hk)

    if(allocated(data%point))       deallocate(data%point)
    if(allocated(data%pixels))      deallocate(data%pixels)
    if(allocated(data%time_gain))   deallocate(data%time_gain)
    if(allocated(data%gain))        deallocate(data%gain)
    if(allocated(data%sigma0))      deallocate(data%sigma0)
    if(allocated(data%alpha))       deallocate(data%alpha)
    if(allocated(data%fknee))       deallocate(data%fknee)
    if(allocated(data%corr))        deallocate(data%corr)
    if(allocated(data%corr_freqs))  deallocate(data%corr_freqs)
    if(allocated(data%diode_stats)) deallocate(data%diode_stats)
    if(allocated(data%filter_par))  deallocate(data%filter_par)
    !if(allocated(data%point_sidelobe))deallocate(data%point_sidelobe)
    if(allocated(data%point_objrel))deallocate(data%point_objrel)
  end subroutine

  subroutine write_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "time",         data%time)
    call write_hdf(file, "decimation",   data%decimation)
    call write_hdf(file, "samprate",     data%samprate)
    call write_hdf(file, "tod",          data%tod)
    call write_hdf(file, "tp",           data%tp)
    call write_hdf(file, "orig_point",   data%orig_point)
    call write_hk_hdf(file, data%hk)
    call close_hdf_file(file)
  end subroutine

  subroutine write_l3_file(filename, data, planck)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: i
    logical(lgt), optional       :: planck
    logical(lgt)                 :: planck_
    planck_  = .false.; if(present(planck)) planck_ = planck

    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "time",         data%time)
    call write_hdf(file, "decimation",   data%decimation)
    call write_hdf(file, "nside",        data%nside)
    call write_hdf(file, "scanfreq",     data%scanfreq)
    call write_hdf(file, "samprate",     data%samprate)
    call write_hdf(file, "coord_sys",    data%coord_sys)
    call write_hdf(file, "tod",          data%tod)
    call write_hdf(file, "orig_point",   data%orig_point)
    call write_hdf(file, "point",        data%point)
    call write_hdf(file, "pixels",       data%pixels)
    call write_hdf(file, "sigma0",       data%sigma0)
    call write_hdf(file, "alpha",        data%alpha)
    call write_hdf(file, "fknee",        data%fknee)
    call write_hdf(file, "corr",         data%corr)
    call write_hdf(file, "corr_freqs",   data%corr_freqs)
    call write_hdf(file, "time_gain",    data%time_gain)
    call write_hdf(file, "gain",         data%gain)
    call write_hdf(file, "stats",       data%stats)
    call write_hdf(file, "diode_stats", data%diode_stats)
    call write_hdf(file, "filter_par",  data%filter_par)
    if (.not. planck_) then
       call write_hdf(file, "tp",           data%tp)
       call write_hdf(file, "point_objrel", data%point_objrel)
       call write_hk_hdf(file, data%hk)
       call create_hdf_group(file, "apex")
       call write_hk_elem   (file, "apex", data%hk%apex)
    end if
    call close_hdf_file(file)
  end subroutine

  subroutine write_hk_hdf(file, hk)
    implicit none
    type(hdf_file)  :: file
    type(hk_struct) :: hk
    if(allocated(hk%bias%time)) then
      call create_hdf_group(file, "bias")
      call write_hk_elem   (file, "bias", hk%bias)
    end if
    if(allocated(hk%cryo%time)) then
      call create_hdf_group(file, "cryo")
      call write_hk_elem   (file, "cryo", hk%cryo)
    end if
    if(allocated(hk%encl%time)) then
      call create_hdf_group(file, "encl")
      call write_hk_elem   (file, "encl", hk%encl)
    end if
    if(allocated(hk%peri%time)) then
      call create_hdf_group(file, "peri")
      call write_hk_elem   (file, "peri", hk%peri)
    end if
  end subroutine

  subroutine write_hk_elem(file, name, hkt)
    implicit none
    type(hdf_file)   :: file
    type(hk_type)    :: hkt
    character(len=*) :: name
    call write_hdf(file, trim(name) // "/time",  hkt%time)
    call write_hdf(file, trim(name) // "/value", hkt%value)
  end subroutine

  ! Concatenate a set of lx_structs. We assume that:
  ! The structs in in are compatible, and that they are
  ! either level2 or level3 (not something in between).
  ! This means that the level2-only parts always will be allocated,
  ! and that if one level3 part exists, they all do
  subroutine cat_lx_structs(in, out)
    implicit none
    type(lx_struct), intent(in)    :: in(:)
    type(lx_struct), intent(inout) :: out
    integer(i4b)                   :: i, j, jg, n, m, mg, ndi, nmod, npix, ng
    logical(lgt)                   :: l3
    real(dp)                       :: f
    logical(lgt), allocatable      :: hit(:)

    call free_lx_struct(out)

    ndi  = size(in(1)%tod,2)
    l3   = allocated(in(1)%point)
    n    = 0
    ng   = 0
    do i = 1, size(in)
      n = n + size(in(i)%time)
      if(l3) ng = ng + size(in(i)%time_gain)
    end do
    allocate(out%time(n), out%tod(n,ndi), out%orig_point(3,n), out%tp(n,ndi))
    if(l3) then
       allocate(out%point(3,n,size(in(1)%point,3)))
       allocate(out%point_objrel(3,n,size(in(1)%point_objrel,3)))
       allocate(out%time_gain(ng), out%gain(ng,ndi))
    end if
    j  = 0
    jg = 0
    mg = 0
    do i = 1, size(in)
       m  = size(in(i)%time)
       out%time       (  1+j:m+j  ) = in(i)%time
       out%tod        (  1+j:m+j,:) = in(i)%tod
       out%tp         (  1+j:m+j,:) = in(i)%tp
       out%orig_point (:,1+j:m+j  ) = in(i)%orig_point
       if(l3) then
          mg = size(in(i)%time_gain)
          out%point(:,1+j:m+j,:)        = in(i)%point
          out%point_objrel(:,1+j:m+j,:) = in(i)%point_objrel
          out%time_gain(1+jg:mg+jg)     = in(i)%time_gain
          out%gain(1+jg:mg+jg,:)        = in(i)%gain
       end if
       j  = j+m
       jg = jg + mg
    end do

    ! These have fixed length
    out%decimation = in(1)%decimation
    out%samprate   = in(1)%samprate
    out%nside      = in(1)%nside
    out%coord_sys  = in(1)%coord_sys
    out%scanfreq   = in(1)%scanfreq

    ! These could actually differ between the files. We make
    ! a weighted average
    if(l3) then
       allocate(out%sigma0(ndi), out%alpha(ndi), out%fknee(ndi), out%corr(ndi,ndi,2))
       allocate(out%diode_stats(ndi,size(in(1)%diode_stats,2)))
       allocate(out%filter_par(ndi,size(in(1)%filter_par,2)))
       out%sigma0 = 0; out%alpha = 0; out%fknee = 0; out%corr = 0
       out%stats  = 0; out%diode_stats = 0; out%filter_par = 0
       do i = 1, size(in)
          f = real(size(in(i)%time),dp)/n
          out%sigma0 = out%sigma0 + in(i)%sigma0 * f
          out%alpha  = out%alpha  + in(i)%alpha  * f
          out%fknee  = out%fknee  + in(i)%fknee  * f
          out%corr   = out%corr   + in(i)%corr   * f
          out%stats  = out%stats  + in(i)%stats  * f
          out%diode_stats = out%diode_stats + in(i)%diode_stats*f
          out%filter_par  = out%filter_par  + in(i)%filter_par*f
       end do
    end if

    ! Pixels need to be completely rebuilt
    if(l3) then
       allocate(hit(0:12*out%nside**2))
       hit = .false.
       do i = 1, size(in)
          hit(in(i)%pixels) = .true.
       end do
       allocate(out%pixels(count(hit)))
       j = 0
       do i = 0, size(hit)-1
          if(.not. hit(i)) cycle
          j = j+1
          out%pixels(j) = i
       end do
       deallocate(hit)
    end if

    ! Ok, only housekeeping is left
    call cat_hk_struct(in%hk, out%hk)
  end subroutine

end module comap_Lx_mod

