module comap_lx_mod
  use healpix_types
  use comap_defs
  use quiet_mpi_mod
  use quiet_hdf_mod
  use quiet_fft_mod
  use quiet_utils
  implicit none

  type lx_struct
     ! Level 1 fields
     real(dp)                                        :: mjd_start
     real(dp)                                        :: samprate
     real(dp),     allocatable, dimension(:)         :: time
     real(dp),     allocatable, dimension(:)         :: sec         ! (time) in seconds since start
     real(dp),     allocatable, dimension(:,:,:)     :: nu          ! (freq, sideband, detector)
     real(sp),     allocatable, dimension(:,:,:,:)   :: tod         ! (time, freq, sideband, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: tod_mean    ! (freq, sideband, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: point_cel   ! Celestial; (RA/dec/psi, time, det)
     real(sp),     allocatable, dimension(:,:,:)     :: point_tel   ! Horizon; (az/el/dk, time, det)
     integer(i4b), allocatable, dimension(:)         :: flag        ! Status flag per time sample
     integer(i4b), allocatable, dimension(:)         :: pixels      ! Active pixels/detectors
     integer(i4b), allocatable, dimension(:,:,:)     :: n_nan       ! number of nan values for each frequency
     real(dp),     allocatable, dimension(:,:,:,:)   :: Tsys        ! (start/stop or middle, freq, sb,detector)
     real(dp),     allocatable, dimension(:)         :: t_hot       ! Ambient temperature in K
     integer(i4b), allocatable, dimension(:)         :: amb_state   ! Ambient load in/out
     real(dp),     allocatable, dimension(:)         :: amb_time    ! Ambient time in MJD

     ! Level 2 fields
     integer(i4b)                                    :: polyorder     ! Polynomial order for frequency filter
     integer(i4b)                                    :: n_pca_comp    ! Number of leading pca-components to subtract
     integer(i4b)                                    :: decimation_time, decimation_nu
     integer(i4b)                                    :: mask_outliers
     real(sp),     allocatable, dimension(:,:,:)     :: freqmask_full ! Full-resolution mask; (freq, sideband, detector)
     integer(i4b), allocatable, dimension(:,:,:)     :: freqmask_reason ! the (first) reason for masking a specific frequency
     real(sp),     allocatable, dimension(:,:,:)     :: freqmask      ! Reduced resolution mask; (freq, sideband, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: point         ! Sky coordinates; (phi/theta/psi,time,det)
     real(sp),     allocatable, dimension(:,:,:)     :: mean_tp
     real(sp),     allocatable, dimension(:,:,:,:)   :: tod_poly      ! Poly-filter TOD coefficients (time,0:poly,sb,det)
     real(sp),     allocatable, dimension(:,:,:)     :: var_fullres   ! Full-resolution variance (freq,sb,det)
     real(sp),     allocatable, dimension(:,:,:,:)   :: pca_ampl      ! amplitudes of pca-components (frec,sb,det,comp)
     real(sp),     allocatable, dimension(:,:)       :: pca_comp      ! actual pca component timestreams (time,comp)
     real(sp),     allocatable, dimension(:)         :: pca_eigv      ! eigenvalues of pca components (comp)
     real(sp),     allocatable, dimension(:,:)       :: acceptrate    ! fraction of freqs not masked (sb,det)
     real(sp),     allocatable, dimension(:,:,:,:)   :: diagnostics   ! various diagnostics used to make freqmask
     real(sp),     allocatable, dimension(:,:)       :: cut_params    ! means and stds used for the different diagnostics
     real(dp),     allocatable, dimension(:,:,:)     :: sigma0, alpha, fknee ! (freq, nsb, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: gain                 ! (freq_fullres, nsb, detector)
     real(dp),     allocatable, dimension(:,:,:)     :: Tsys_lowres   ! (freq, sb,detector)

     ! Level 3 fields
!!$     integer(i4b)                                    :: coord_sys
!!$     real(dp)                                        :: scanfreq(2), pixsize 
!!$     real(dp)                                        :: point_lim(4)         ! (RA_min, RA_max, dec_min, dec_max)
!!$     real(dp),     allocatable, dimension(:)         :: time_gain            ! (time)
!!$     real(dp)                                        :: stats(ST_NUM)
!!$     real(dp),     allocatable, dimension(:,:,:,:)   :: det_stats            ! (freq, detector, nsb, stat)
!!$     real(dp),     allocatable, dimension(:,:,:,:)   :: filter_par           ! (freq, detector, nsb, param)
!!$     real(sp),     allocatable, dimension(:,:,:)     :: weather_temp         ! (nsamp, nsb, detector)
!!$     real(sp),     allocatable, dimension(:,:,:)     :: rel_gain             ! (freq, nsb, detector)
!!$     real(dp),     allocatable, dimension(:,:,:)     :: sigma0_poly, alpha_poly, fknee_poly ! (poly, nsb, detector)

  end type lx_struct

contains

  subroutine read_l1_file(filename, data, id, only_point, freqmask, init)
    implicit none
    character(len=*), intent(in)           :: filename
    integer(i4b),     intent(in)           :: id
    logical(lgt),     intent(in), optional :: only_point, init
    real(sp), dimension(:,:,:), intent(inout), optional :: freqmask
    type(lx_struct)                        :: data
    type(hdf_file)                         :: file
    integer(i4b)                           :: i, j, k, l, nsamp, nsamp_tot, nfreq, ndet, npoint, nsb, ext4(4), ext1(1), numbad, temp_samp
    logical(lgt)                           :: all, ok, init_
    real(dp)                               :: t1, t2
    integer(i4b), allocatable, dimension(:)       :: buffer_int
    real(dp),     allocatable, dimension(:)       :: buffer_1d
    real(sp),     allocatable, dimension(:,:,:,:) :: buffer_4d
    all = .true.; if (present(only_point)) all = .not. only_point
    init_ = .true.; if (present(init)) init_ = init
    if (init_) call free_lx_struct(data) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "spectrometer/tod", ext4)
    nsamp_tot = ext4(1); nfreq = ext4(2) ; nsb = ext4(3); ndet = ext4(4)

    call get_size_hdf(file, "hk/antenna0/env/ambientLoadTemp", ext1)
    temp_samp = ext1(1)

    allocate(data%t_hot(temp_samp))
    allocate(data%amb_state(temp_samp))
    allocate(data%amb_time(temp_samp))
    allocate(data%point_tel(3,nsamp_tot,ndet))
    allocate(data%point_cel(3,nsamp_tot,ndet))
    allocate(data%pixels(ndet))
    if (all) allocate(data%nu(nfreq,nsb,ndet))

    ! Read telescope coordinates
    call read_hdf(file, "spectrometer/pixel_pointing/pixel_az",            data%point_tel(1,:,:))
    call read_hdf(file, "spectrometer/pixel_pointing/pixel_el",            data%point_tel(2,:,:))
    data%point_tel(3,:,:) = 0.d0

    ! Read ecliptic coordinates
    call read_hdf(file, "spectrometer/pixel_pointing/pixel_ra",            data%point_cel(1,:,:))
    call read_hdf(file, "spectrometer/pixel_pointing/pixel_dec",           data%point_cel(2,:,:))
    data%point_cel(3,:,:) = 0.d0

    ! Read ambient temp
    call read_hdf(file, "hk/antenna0/env/ambientLoadTemp", data%t_hot)
    call read_hdf(file, "hk/antenna0/vane/state",          data%amb_state)
    call read_hdf(file, "hk/antenna0/vane/utc",            data%amb_time)

    ! Read feed information
    call read_hdf(file, "spectrometer/feeds",               data%pixels)
    if (all) call read_hdf(file, "spectrometer/frequency",       data%nu(:,:,1))
    do i = 2, ndet
       data%nu(:,:,i) = data%nu(:,:,1)
    end do

    ! Do elements that may have NaNs
    allocate(buffer_1d(nsamp_tot))
    if (all) allocate(buffer_int(nsamp_tot))
    if (all) allocate(buffer_4d(nsamp_tot,nfreq,nsb,ndet))

    ! Find number of samples at end of file with NaNs
    if (all) call read_hdf(file, "spectrometer/tod",  buffer_4d)
    if (present(freqmask)) then
       allocate(data%n_nan(nfreq,nsb,ndet))
       ! Update frequency mask with channels that are all NaNs
       do i = 1, ndet
          do j = 1, nsb
             do k = 1, nfreq
                !write(*,*) "hei", freqmask(k,j,i), k, j, i
                if (freqmask(k,j,i) == 0.) cycle
                numbad = count(buffer_4d(:,k,j,i) .ne. buffer_4d(:,k,j,i))
                data%n_nan(k,j,i) = numbad
                if (numbad > 0.1*nsamp_tot) then
                   write(*,fmt='(a,i8,i6,i4,i8)') '  Removing frequency with >10% NaNs -- ', id, i, j, k
                   freqmask(k,j,i) = 0.
                end if
             end do
          end do
       end do

       ! Trim end for NaNs
       nsamp = nsamp_tot
       do while (nsamp > 0)
          ok = .true.
          do i = 1, ndet
             do j = 1, nsb
                do k = 1, nfreq
                   if (freqmask(k,j,i) == 0.) cycle
                   if (buffer_4d(nsamp,k,j,i) .ne. buffer_4d(nsamp,k,j,i)) then
                      ok = .false.
                      exit
                   end if
                end do
             end do
          end do
          if (ok) then
             exit
          else
             nsamp = nsamp-1
          end if
       end do
    else
       nsamp = nsamp_tot
    end if

    if (nsamp /= nsamp_tot) then
       write(*,*) '  Number of NaN elements in ', id, ' = ', nsamp_tot-nsamp, ' of ', nsamp_tot
    end if

    allocate(data%time(nsamp))
    if (all) allocate(data%tod(nsamp,nfreq,nsb,ndet))
    if (all) allocate(data%tod_mean(nfreq,nsb,ndet))
    if (all) allocate(data%flag(nsamp))

    call read_hdf(file, "spectrometer/MJD",           buffer_1d)
    data%time         = buffer_1d(1:nsamp)
    if (all) data%tod = buffer_4d(1:nsamp,:,:,:)
    data%mjd_start    = minval(data%time)
    data%samprate     = 1.d0 / (3600.d0*24.d0*(data%time(2)-data%time(1)))
    data%flag         = 0
    call close_hdf_file(file)
    deallocate(buffer_1d)
    if (all) deallocate(buffer_4d,buffer_int)
  end subroutine read_l1_file


  subroutine read_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: nsamp, nfreq, nfreq_full, nsb, ndet, npoint, ext(7), poly
    call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "tod", ext)
    nsamp = ext(1); nfreq = ext(2) ; nsb = ext(3); ndet = ext(4)
    allocate(data%time(nsamp), data%tod(nsamp,nfreq,nsb,ndet), data%pixels(ndet), data%tod_mean(1024, nsb, ndet))
    allocate(data%flag(nsamp))
    call get_size_hdf(file, "point_tel", ext)
    npoint = ext(1); nsamp = ext(2)
    allocate(data%point_tel(npoint,nsamp,ndet), data%point_cel(npoint,nsamp,ndet), data%nu(nfreq,nsb,ndet))
    call read_hdf(file, "decimation_time",  data%decimation_time)
    call read_hdf(file, "decimation_nu",    data%decimation_nu)
    call read_hdf(file, "samprate",         data%samprate)
    call read_hdf(file, "time",             data%time)
    call read_hdf(file, "nu",               data%nu)
    call read_hdf(file, "tod",              data%tod)
    call read_hdf(file, "point_tel",        data%point_tel)
    !call read_hdf(file, "tod_mean",         data%tod_mean )
    call read_hdf(file, "point_cel",        data%point_cel)
    !call read_hdf(file, "flag",             data%flag)
    nfreq_full = nfreq*data%decimation_nu
    allocate(data%freqmask_full(nfreq_full,nsb,ndet), data%freqmask(nfreq,nsb,ndet), data%mean_tp(nfreq_full,nsb,ndet))
    call read_hdf(file, "freqmask",         data%freqmask)    
    call read_hdf(file, "freqmask_full",    data%freqmask_full)
    call read_hdf(file, "freqmask_reason",  data%freqmask_reason)
    call read_hdf(file, "mean_tp",          data%mean_tp)
    allocate(data%tsys(2,nfreq_full,nsb,ndet), data%tsys_lowres(nfreq,nsb,ndet))
    call read_hdf(file, "Tsys",             data%tsys)
    call read_hdf(file, "Tsys_lowres",      data%tsys_lowres)
    allocate(data%n_nan(nfreq_full,nsb,ndet))
    call read_hdf(file, "n_nan",            data%n_nan)    

    call read_hdf(file, "polyorder",        data%polyorder)
    !if (data%polyorder >= 0) then
    !   allocate(data%tod_poly(nsamp,0:data%polyorder,nsb,ndet))
    !   call read_hdf(file, "tod_poly",         data%tod_poly)
    !end if
    call read_hdf(file, "pixels",           data%pixels)
    allocate(data%var_fullres(nfreq_full,nsb,ndet))
    call read_hdf(file, "var_fullres",      data%var_fullres)
    call read_hdf(file, "n_pca_comp",         data%n_pca_comp)
    if (data%n_pca_comp > 0) then
       allocate(data%pca_ampl(nfreq_full,nsb,ndet,data%n_pca_comp)) 
       allocate(data%pca_comp(nsamp,data%n_pca_comp))
       allocate(data%pca_eigv(data%n_pca_comp))
       call read_hdf(file, "pca_ampl",         data%pca_ampl)
       call read_hdf(file, "pca_comp",         data%pca_comp)
       call read_hdf(file, "pca_eigv",         data%pca_eigv)
    end if

    call read_hdf(file, "mask_outliers",    data%mask_outliers)
    if (data%mask_outliers == 1) then
       allocate(data%diagnostics(nfreq_full,nsb,ndet,5))
       allocate(data%cut_params(2,5))
       allocate(data%acceptrate(nsb,ndet))
       call read_hdf(file, "acceptrate",       data%acceptrate)
       call read_hdf(file, "diagnostics",       data%diagnostics)
       call read_hdf(file, "cut_params",        data%cut_params)
    end if
    allocate(data%sigma0(nfreq,nsb,ndet), data%alpha(nfreq,nsb,ndet), data%fknee(nfreq,nsb,ndet))
    call read_hdf(file, "sigma0", data%sigma0)
    call read_hdf(file, "alpha",  data%alpha)
    call read_hdf(file, "fknee",  data%fknee)    
    call close_hdf_file(file)
    !write(*,*) nfreq,nsb,ndet
  end subroutine read_l2_file


  ! ! Where should this sub logically be?
  ! subroutine decimate(time, time_full, tod, tod_full, point, point_full, dec)
  !   implicit none
  !   integer(i4b)                                      :: n, nfreq, ndet, npt, dec, i, j, k, l, ind
  !   real(dp), dimension(:),      intent(inout)        :: time, time_full
  !   real(dp), dimension(:,:,:),  intent(inout)        :: tod, tod_full
  !   real(dp), dimension(:,:,0:), intent(inout)        :: point, point_full

  !   if (dec>0) then
  !      n     = size(tod,1)
  !      nfreq = size(tod,2)
  !      ndet  = size(tod,3)
  !      npt   = size(point,1)

  !      ! Averaging over every (dec) elements of the time dimension of time, tod and pointing arrays
  !      ind = 1
  !      do i=1,n
  !         time(i) = mean(time_full(ind:ind+dec-1))
  !         do j=1,ndet
  !            do l=1, nfreq
  !               tod(i,l,j) = mean(tod_full(ind:ind+dec-1,l,j))
  !            end do
  !            do k=1,npt
  !               ! do I need to make the angles safe?
  !               point(k,i,j-1) = mean(point_full(k,ind:ind+dec-1,j-1))
  !            end do
  !         end do
  !         ind = ind + dec
  !      end do

  !   else
  !      time = time_full
  !      tod = tod_full
  !      point = point_full
  !   end if
  ! end subroutine decimate


!!$  subroutine read_l3_file(filename, data)
!!$    implicit none
!!$    character(len=*), intent(in) :: filename
!!$    type(lx_struct)              :: data
!!$    type(hdf_file)               :: file
!!$    integer(i4b) :: npoint, nsamp, nfreq, nfreq_full, nsb, ndet, nmod, ext(7)
!!$
!!$    call free_lx_struct(data)
!!$    call read_l2_file(filename, data)
!!$    call open_hdf_file(filename, file, "r")
!!$    call read_hdf(file, "scanfreq",  data%scanfreq)    
!!$    !call read_hdf(file, "pixsize", data%pixsize)
!!$    call read_hdf(file, "point_lim", data%point_lim)
!!$    ! Read pointing
!!$    call get_size_hdf(file, "point", ext)
!!$    npoint = ext(1); nsamp = ext(2); ndet = ext(3)
!!$    allocate(data%point(npoint,nsamp,ndet))
!!$    call read_hdf(file, "point",     data%point)
!!$
!!$    call read_hdf(file, "coord_sys", data%coord_sys)
!!$    ! Read gain
!!$    call get_size_hdf(file, "time_gain", ext)
!!$    nsamp = ext(1)
!!$    allocate(data%time_gain(nsamp))
!!$    call read_hdf(file, "time_gain", data%time_gain)
!!$    call get_size_hdf(file, "gain", ext)
!!$    nsamp = ext(1); nfreq = ext(2); nsb = ext(3); ndet = ext(4)
!!$    allocate(data%gain(nsamp,nfreq,nsb,ndet))
!!$    call read_hdf(file, "gain", data%gain)
!!$    !write(*,*) data%gain
!!$    ! Read noise parameters
!!$    !call get_size_hdf(file, "corr", ext)
!!$    call get_size_hdf(file, "sigma0", ext)
!!$    nfreq = ext(1); nsb = ext(2); ndet = ext(3)
!!$    allocate(data%sigma0(nfreq,nsb,ndet),data%alpha(nfreq,nsb,ndet),data%fknee(nfreq,nsb,ndet))
!!$    !allocate(data%corr(nfreq,nfreq,ndet,ndet))
!!$    !allocate(data%det_stats(nfreq,ndet,NUM_DET_STATS), data%filter_par(nfreq,ndet,NUM_FILTR_PAR))
!!$    !allocate(data%det_stats(ndet,nfreq,1))
!!$    !allocate(data%filter_par(4,nfreq,ndet))
!!$    call read_hdf(file, "sigma0", data%sigma0)
!!$    call read_hdf(file, "alpha",  data%alpha)
!!$    call read_hdf(file, "fknee",  data%fknee)
!!$    !call read_hdf(file, "corr",   data%corr)
!!$    ! Read stats
!!$    call read_hdf(file, "stats", data%stats)
!!$    nfreq_full = nfreq*data%decimation_nu
!!$    ! Read filter parameters
!!$    !call read_hdf(file, "filter_par",  data%filter_par)
!!$    call close_hdf_file(file)
!!$  end subroutine read_l3_file

  subroutine free_lx_struct(data)
    implicit none
    type(lx_struct) :: data
    if(allocated(data%time))        deallocate(data%time)
    if(allocated(data%nu))          deallocate(data%nu)
    if(allocated(data%tod))         deallocate(data%tod)
    if(allocated(data%tod_mean))    deallocate(data%tod_mean)
    if(allocated(data%point_tel))   deallocate(data%point_tel)
    if(allocated(data%point_cel))   deallocate(data%point_cel)
    if(allocated(data%flag))        deallocate(data%flag)

    if(allocated(data%point))         deallocate(data%point)
    if(allocated(data%gain))          deallocate(data%gain)
    if(allocated(data%sigma0))        deallocate(data%sigma0)
    if(allocated(data%alpha))         deallocate(data%alpha)
    if(allocated(data%Tsys))          deallocate(data%Tsys)
    if(allocated(data%Tsys_lowres))   deallocate(data%Tsys_lowres)
    if(allocated(data%fknee))         deallocate(data%fknee)
    if(allocated(data%freqmask))      deallocate(data%freqmask)
    if(allocated(data%freqmask_full)) deallocate(data%freqmask_full)
    if(allocated(data%freqmask_reason)) deallocate(data%freqmask_reason)
    if(allocated(data%mean_tp))       deallocate(data%mean_tp)
    if(allocated(data%n_nan))         deallocate(data%n_nan)
    if(allocated(data%tod_poly))      deallocate(data%tod_poly)
    if(allocated(data%tod_mean))      deallocate(data%tod_mean)
    if(allocated(data%pixels))        deallocate(data%pixels)
    if(allocated(data%var_fullres))   deallocate(data%var_fullres)
    if(allocated(data%pca_ampl))      deallocate(data%pca_ampl)
    if(allocated(data%pca_comp))      deallocate(data%pca_comp)
    if(allocated(data%pca_eigv))      deallocate(data%pca_eigv)
    if(allocated(data%acceptrate))    deallocate(data%acceptrate)
    if(allocated(data%diagnostics))   deallocate(data%diagnostics)
    if(allocated(data%cut_params))    deallocate(data%cut_params)
  end subroutine

  subroutine write_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "samprate",          data%samprate)
    call write_hdf(file, "mjd_start",         data%mjd_start)
    call write_hdf(file, "decimation_time",   data%decimation_time)
    call write_hdf(file, "decimation_nu",     data%decimation_nu)
    call write_hdf(file, "time",              data%time)
    call write_hdf(file, "nu",                data%nu)
    call write_hdf(file, "tod",               data%tod)
    call write_hdf(file, "point_cel",         data%point_cel)
    call write_hdf(file, "point_tel",         data%point_tel)
    !call write_hdf(file, "flag",              data%flag)
    !write(*,*) "right before", data%Tsys(1, 1, 1, 1)
    call write_hdf(file, "Tsys",              data%Tsys)
    call write_hdf(file, "Tsys_lowres",       data%Tsys_lowres)
    if (allocated(data%sigma0))    call write_hdf(file, "sigma0",            data%sigma0)
    if (allocated(data%alpha))     call write_hdf(file, "alpha",             data%alpha)
    if (allocated(data%fknee))     call write_hdf(file, "fknee",             data%fknee)   
    call write_hdf(file, "freqmask",          data%freqmask)
    call write_hdf(file, "freqmask_full",     data%freqmask_full)
    call write_hdf(file, "freqmask_reason",   data%freqmask_reason)
    call write_hdf(file, "n_nan",             data%n_nan)
    if (allocated(data%mean_tp)) call write_hdf(file, "mean_tp",           data%mean_tp)
    call write_hdf(file, "polyorder",         data%polyorder)
    if (data%polyorder >= 0) then
       call write_hdf(file, "tod_poly",         data%tod_poly)
    end if
    call write_hdf(file, "pixels",            data%pixels)
    call write_hdf(file, "var_fullres",       data%var_fullres)
    call write_hdf(file, "n_pca_comp",        data%n_pca_comp)
    if (data%n_pca_comp > 0) then
       call write_hdf(file, "pca_ampl",       data%pca_ampl)
       call write_hdf(file, "pca_comp",       data%pca_comp)
       call write_hdf(file, "pca_eigv",       data%pca_eigv)
    end if
    call write_hdf(file, "mask_outliers",     data%mask_outliers)
    if (data%mask_outliers == 1) then
       call write_hdf(file, "acceptrate",     data%acceptrate)
       call write_hdf(file, "diagnostics",    data%diagnostics)
       call write_hdf(file, "cut_params",     data%cut_params)
    end if
    call close_hdf_file(file)
  end subroutine

!!$  subroutine write_l3_file(filename, data)
!!$    implicit none
!!$    character(len=*), intent(in) :: filename
!!$    type(lx_struct)              :: data
!!$    type(hdf_file)               :: file
!!$    integer(i4b)                 :: i
!!$
!!$    call open_hdf_file(filename, file, "w")
!!$    call write_hdf(file, "time",              data%time)
!!$    call write_hdf(file, "nu",                data%nu)
!!$    call write_hdf(file, "decimation_time",   data%decimation_time)
!!$    call write_hdf(file, "decimation_nu",     data%decimation_nu)
!!$    call write_hdf(file, "mjd_start",         data%mjd_start)
!!$    call write_hdf(file, "sec",               data%sec)
!!$    call write_hdf(file, "scanfreq",          data%scanfreq)
!!$    call write_hdf(file, "samprate",          data%samprate)
!!$    call write_hdf(file, "coord_sys",         data%coord_sys)
!!$    call write_hdf(file, "pixsize",           data%pixsize)
!!$    call write_hdf(file, "point_lim",         data%point_lim)
!!$    call write_hdf(file, "tod",               data%tod)
!!$    call write_hdf(file, "point_tel",         data%point_tel)
!!$    call write_hdf(file, "point_cel",         data%point_cel)
!!$    call write_hdf(file, "point",             data%point)
!!$    if (allocated(data%sigma0))    call write_hdf(file, "sigma0",            data%sigma0)
!!$    if (allocated(data%alpha))     call write_hdf(file, "alpha",             data%alpha)
!!$    if (allocated(data%fknee))     call write_hdf(file, "fknee",             data%fknee)
!!$    if (allocated(data%time_gain)) call write_hdf(file, "time_gain",         data%time_gain)
!!$    if (allocated(data%gain))      call write_hdf(file, "gain",              data%gain)
!!$    call write_hdf(file, "stats",             data%stats)
!!$    !call write_hdf(file, "det_stats",         data%det_stats)
!!$    !call write_hdf(file, "filter_par",        data%filter_par)
!!$    !call write_hdf(file, "flag",              data%flag)
!!$    call write_hdf(file, "pixels",            data%pixels)
!!$    call write_hdf(file, "freqmask",          data%freqmask)
!!$    call write_hdf(file, "freqmask_full",     data%freqmask_full)
!!$    call write_hdf(file, "n_nan",             data%n_nan)
!!$    if (allocated(data%mean_tp)) call write_hdf(file, "mean_tp",           data%mean_tp)
!!$    call write_hdf(file, "polyorder",         data%polyorder)
!!$    call write_hdf(file, "var_fullres",       data%var_fullres)
!!$    call write_hdf(file, "n_pca_comp",        data%n_pca_comp)
!!$    if (data%n_pca_comp > 0) then
!!$       call write_hdf(file, "pca_ampl",          data%pca_ampl)
!!$       call write_hdf(file, "pca_comp",          data%pca_comp)
!!$       call write_hdf(file, "pca_eigv",          data%pca_eigv)
!!$    end if
!!$    call write_hdf(file, "acceptrate",        data%acceptrate)
!!$    call write_hdf(file, "diagnostics",       data%diagnostics)
!!$    call write_hdf(file, "cut_params",        data%cut_params)
!!$    
!!$    if (data%polyorder >= 0) then
!!$       call write_hdf(file, "tod_poly",       data%tod_poly)
!!$       call write_hdf(file, "sigma0_poly",    data%sigma0_poly)
!!$       call write_hdf(file, "alpha_poly",     data%alpha_poly)
!!$       call write_hdf(file, "fknee_poly",     data%fknee_poly)
!!$    end if
!!$    call close_hdf_file(file)
!!$  end subroutine


  ! Concatenate a set of lx_structs. We assume that:
  ! The structs in in are compatible, and that they are
  ! either level2 or level3 (not something in between).
  ! This means that the level2-only parts always will be allocated,
  ! and that if one level3 part exists, they all do
  subroutine cat_lx_structs(in, out)
    implicit none
    type(lx_struct), intent(in)    :: in(:)
    type(lx_struct), intent(inout) :: out
    integer(i4b)                   :: i, j, jg, n, m, mg, ndet, nmod, npix, ng, nf, nsb
    logical(lgt)                   :: l3
    real(dp)                       :: f

    call free_lx_struct(out)
    
    
!!$    nf   = size(in(1)%tod,2)
!!$    nsb  = size(in(1)%tod,3)
!!$    ndet = size(in(1)%tod,4)
!!$    l3   = allocated(in(1)%point)
!!$    n    = 0
!!$    ng   = 0
!!$    do i = 1, size(in)
!!$      n = n + size(in(i)%time)
!!$      if(l3) ng = ng + size(in(i)%time_gain)
!!$    end do
!!$    allocate(out%time(n), out%tod(n,nf,nsb,ndet), out%point_tel(3,n), out%point_cel(3,n))
!!$    if(l3) then
!!$       allocate(out%point(3,n))!,size(in(1)%point)))
!!$       allocate(out%time_gain(ng), out%gain(ng,nf,ndet))
!!$    end if
!!$    j  = 0
!!$    jg = 0
!!$    mg = 0
!!$    do i = 1, size(in)
!!$       m  = size(in(i)%time)
!!$       out%time       (  1+j:m+j  )   = in(i)%time
!!$       out%tod        (  1+j:m+j,:,:) = in(i)%tod
!!$       out%point_tel  (:,1+j:m+j  )   = in(i)%point_tel
!!$       out%point_cel  (:,1+j:m+j  )   = in(i)%point_cel
!!$       if(l3) then
!!$          mg = size(in(i)%time_gain)
!!$          out%point(:,1+j:m+j)          = in(i)%point
!!$          out%time_gain(1+jg:mg+jg)     = in(i)%time_gain
!!$          out%gain(1+jg:mg+jg,:,:)      = in(i)%gain
!!$       end if
!!$       j  = j+m
!!$       jg = jg + mg
!!$    end do
!!$
!!$    ! These have fixed length
!!$    out%decimation_time = in(1)%decimation_time
!!$    out%decimation_nu   = in(1)%decimation_nu
!!$    out%samprate        = in(1)%samprate
!!$    out%coord_sys       = in(1)%coord_sys
!!$    out%scanfreq        = in(1)%scanfreq
!!$    out%pixsize         = in(1)%pixsize
!!$    out%point_lim       = in(1)%point_lim
!!$
!!$    ! These could actually differ between the files. We make
!!$    ! a weighted average
!!$    if(l3) then
!!$       allocate(out%sigma0(nf,ndet), out%alpha(nf,ndet), out%fknee(nf,ndet))
!!$       allocate(out%det_stats(ndet,nf,size(in(1)%det_stats,2)))
!!$       allocate(out%filter_par(ndet,nf,size(in(1)%filter_par,2)))
!!$       out%sigma0 = 0; out%alpha = 0; out%fknee = 0
!!$       out%stats  = 0; out%det_stats = 0; out%filter_par = 0
!!$       do i = 1, size(in)
!!$          f = real(size(in(i)%time),dp)/n
!!$          out%sigma0 = out%sigma0 + in(i)%sigma0 * f
!!$          out%alpha  = out%alpha  + in(i)%alpha  * f
!!$          out%fknee  = out%fknee  + in(i)%fknee  * f
!!$          out%stats  = out%stats  + in(i)%stats  * f
!!$          out%det_stats  = out%det_stats  + in(i)%det_stats*f
!!$          out%filter_par = out%filter_par + in(i)%filter_par*f
!!$       end do
!!$    end if

  end subroutine


  subroutine copy_lx_struct(lx_in, lx_out)
    type(lx_struct), intent(in)    :: lx_in
    type(lx_struct), intent(inout) :: lx_out
    call free_lx_struct(lx_out)
    
!    write(*,*) "start"
    lx_out%mjd_start = lx_in%mjd_start
    lx_out%samprate = lx_in%samprate
    lx_out%polyorder = lx_in%polyorder
    lx_out%decimation_time = lx_in%decimation_time
    lx_out%decimation_nu = lx_in%decimation_nu
    lx_out%n_pca_comp = lx_in%n_pca_comp
    lx_out%mask_outliers = lx_in%mask_outliers

    if(allocated(lx_in%time))        then
       allocate(lx_out%time(size(lx_in%time,1)))
       lx_out%time = lx_in%time
    end if
    if(allocated(lx_in%nu))          then  
       allocate(lx_out%nu(size(lx_in%nu,1),size(lx_in%nu,2),size(lx_in%nu,3)))
       lx_out%nu = lx_in%nu
    end if
!    write(*,*) "Before tod"
    if(allocated(lx_in%tod))         then  
       allocate(lx_out%tod(size(lx_in%tod,1),size(lx_in%tod,2),size(lx_in%tod,3),size(lx_in%tod,4)))
       lx_out%tod = lx_in%tod
    end if
!    write(*,*) "after tod"
    if(allocated(lx_in%point_tel))   then 
       allocate(lx_out%point_tel(size(lx_in%point_tel,1),size(lx_in%point_tel,2),size(lx_in%point_tel,3)))
       lx_out%point_tel = lx_in%point_tel
    end if
    if(allocated(lx_in%point_cel))   then
       allocate(lx_out%point_cel(size(lx_in%point_cel,1),size(lx_in%point_cel,2),size(lx_in%point_cel,3)))
       lx_out%point_cel = lx_in%point_cel
    end if
!    write(*,*) "after point_tel/cel"
    if(allocated(lx_in%flag))        then
       allocate(lx_out%flag(size(lx_in%flag,1)))
       lx_out%flag = lx_in%flag
    end if
    if(allocated(lx_in%point))         then
       allocate(lx_out%point(size(lx_in%point,1),size(lx_in%point,2),size(lx_in%point,3)))
       lx_out%point = lx_in%point
    end if
    if(allocated(lx_in%gain))          then
       allocate(lx_out%gain(size(lx_in%gain,1),size(lx_in%gain,2),size(lx_in%gain,3)))
       lx_out%gain = lx_in%gain
    end if
!    write(*,*) "after gain"
    if(allocated(lx_in%sigma0))        then
       allocate(lx_out%sigma0(size(lx_in%sigma0,1),size(lx_in%sigma0,2),size(lx_in%sigma0,3)))
       lx_out%sigma0 = lx_in%sigma0
    end if
    if(allocated(lx_in%alpha))         then
       allocate(lx_out%alpha(size(lx_in%alpha,1),size(lx_in%alpha,2),size(lx_in%alpha,3)))
       lx_out%alpha = lx_in%alpha
    end if
    if(allocated(lx_in%fknee))         then
       allocate(lx_out%fknee(size(lx_in%fknee,1),size(lx_in%fknee,2),size(lx_in%fknee,3)))
       lx_out%fknee = lx_in%fknee
    end if
    if(allocated(lx_in%freqmask))      then
       allocate(lx_out%freqmask(size(lx_in%freqmask,1),size(lx_in%freqmask,2),size(lx_in%freqmask,3)))
       lx_out%freqmask = lx_in%freqmask
    end if
    if(allocated(lx_in%Tsys_lowres))      then
       allocate(lx_out%Tsys_lowres(size(lx_in%Tsys_lowres,1),size(lx_in%Tsys_lowres,2),size(lx_in%Tsys_lowres,3)))
       lx_out%Tsys_lowres = lx_in%Tsys_lowres
    end if
    if(allocated(lx_in%Tsys))      then
       allocate(lx_out%Tsys(size(lx_in%Tsys,1),size(lx_in%Tsys,2),size(lx_in%Tsys,3),size(lx_in%Tsys,4)))
       lx_out%Tsys = lx_in%Tsys
    end if
    if(allocated(lx_in%freqmask_full)) then
       allocate(lx_out%freqmask_full(size(lx_in%freqmask_full,1),size(lx_in%freqmask_full,2),size(lx_in%freqmask_full,3)))
       lx_out%freqmask_full = lx_in%freqmask_full
    end if
    if(allocated(lx_in%freqmask_reason)) then
       allocate(lx_out%freqmask_reason(size(lx_in%freqmask_reason,1),size(lx_in%freqmask_reason,2),size(lx_in%freqmask_reason,3)))
       lx_out%freqmask_reason = lx_in%freqmask_reason
    end if
    if(allocated(lx_in%n_nan)) then
       allocate(lx_out%n_nan(size(lx_in%n_nan,1),size(lx_in%n_nan,2),size(lx_in%n_nan,3)))
       lx_out%n_nan = lx_in%n_nan
    end if
    if(allocated(lx_in%mean_tp))       then
       allocate(lx_out%mean_tp(size(lx_in%mean_tp,1),size(lx_in%mean_tp,2),size(lx_in%mean_tp,3)))
       lx_out%mean_tp = lx_in%mean_tp
    end if
    if(allocated(lx_in%tod_poly))      then
       allocate(lx_out%tod_poly(size(lx_in%tod_poly,1),size(lx_in%tod_poly,2),size(lx_in%tod_poly,3),size(lx_in%tod_poly,4)))
       lx_out%tod_poly = lx_in%tod_poly
    end if
    if(allocated(lx_in%pixels))        then
       allocate(lx_out%pixels(size(lx_in%pixels,1)))
       lx_out%pixels = lx_in%pixels
    end if
    if(allocated(lx_in%var_fullres))   then
       allocate(lx_out%var_fullres(size(lx_in%var_fullres,1),size(lx_in%var_fullres,2),size(lx_in%var_fullres,3)))
       lx_out%var_fullres = lx_in%var_fullres
    end if
    
    if(allocated(lx_in%pca_ampl))   then
       allocate(lx_out%pca_ampl(size(lx_in%pca_ampl,1),size(lx_in%pca_ampl,2),size(lx_in%pca_ampl,3),size(lx_in%pca_ampl,4)))
       lx_out%pca_ampl = lx_in%pca_ampl
    end if
    if(allocated(lx_in%pca_comp))   then
       allocate(lx_out%pca_comp(size(lx_in%pca_comp,1),size(lx_in%pca_comp,2)))
       lx_out%pca_comp = lx_in%pca_comp
    end if
    if(allocated(lx_in%pca_eigv))   then
       allocate(lx_out%pca_eigv(size(lx_in%pca_eigv,1)))
       lx_out%pca_eigv = lx_in%pca_eigv
    end if
    if(allocated(lx_in%acceptrate))   then
       allocate(lx_out%acceptrate(size(lx_in%acceptrate,1),size(lx_in%acceptrate,2)))
       lx_out%acceptrate = lx_in%acceptrate
    end if
    if(allocated(lx_in%diagnostics))   then
       allocate(lx_out%diagnostics(size(lx_in%diagnostics,1),size(lx_in%diagnostics,2),size(lx_in%diagnostics,3),size(lx_in%diagnostics,4)))
       lx_out%diagnostics = lx_in%diagnostics
    end if
    if(allocated(lx_in%cut_params))   then
       allocate(lx_out%cut_params(size(lx_in%cut_params,1),size(lx_in%cut_params,2)))
       lx_out%cut_params = lx_in%cut_params
    end if
    if(allocated(lx_in%sec))   then
       allocate(lx_out%sec(size(lx_in%sec,1)))
       lx_out%sec = lx_in%sec
    end if
 !   write(*,*) "end"
  end subroutine copy_lx_struct
  

end module comap_lx_mod

