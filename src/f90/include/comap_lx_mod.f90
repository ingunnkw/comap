module comap_lx_mod
  use healpix_types
  use comap_scan_mod
  use comap_defs
  use comap_detector_mod
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
     real(sp),     allocatable, dimension(:,:,:)     :: sb_mean     ! (time, sideband, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: point_cel   ! Celestial; (RA/dec/psi, time, det)
     real(sp),     allocatable, dimension(:,:,:)     :: point_tel   ! Horizon; (az/el/dk, time, det)
     integer(i4b), allocatable, dimension(:)         :: flag        ! Status flag per time sample
     integer(i4b), allocatable, dimension(:)         :: pixels      ! Active pixels/detectors (i.e. ind2pix)
     integer(i4b), allocatable, dimension(:)         :: pix2ind     ! Which pixel does ind corresp. to
     integer(i4b), allocatable, dimension(:,:,:)     :: n_nan       ! number of nan values for each frequency
     integer(i4b), allocatable, dimension(:)         :: amb_state   ! Ambient load in/out
     real(dp),     allocatable, dimension(:)         :: amb_time    ! Ambient time in MJD

     ! Level 2 fields
     integer(i4b)                                    :: polyorder     ! Polynomial order for frequency filter
     integer(i4b)                                    :: n_pca_comp    ! Number of leading pca-components to subtract
     integer(i4b)                                    :: n_pca_comp_feed    ! Number of leading feed-pca-components to subtract
     integer(i4b)                                    :: decimation_time, decimation_nu
     integer(i4b)                                    :: mask_outliers
     integer(i4b)                                    :: irun          ! Run number identification
     integer(i4b)                                    :: n_cal         ! Number of successful ambient measurments
     integer(i4b)                                    :: n_freq_downsamp   ! Number of frequencies in temporary PCA downsampling
     integer(i4b)                                    :: cal_method  ! (1 = old, 2 = Jonas, 3 = from l1)
     logical(lgt)                                    :: use_freq_filter ! use the new frequency filter (T/F)
     real(dp),     allocatable, dimension(:,:)       :: Thot        ! (nmethods=3, start/end)
     real(dp),     allocatable, dimension(:,:)       :: time_hot    ! (nmethods=3, start/end)
     real(dp),     allocatable, dimension(:,:,:,:,:) :: Phot        ! (nmethods=3, start/end, freq, sb, detector)
!     real(dp),     allocatable, dimension(:,:,:)     :: Pcold       ! (nmethods=3, freq, sb, detector)
     real(dp),     allocatable, dimension(:,:,:)     :: Tsys        ! (start/stop or middle, freq, sb,detector)
     real(dp),     allocatable, dimension(:)         :: t_hot       ! Ambient temperature in K
     real(sp),     allocatable, dimension(:,:,:)     :: freqmask_full ! Full-resolution mask; (freq, sideband, detector)
     integer(i4b), allocatable, dimension(:,:,:)     :: freqmask_reason ! the (first) reason for masking a specific frequency
     real(sp),     allocatable, dimension(:,:,:)     :: freqmask      ! Reduced resolution mask; (freq, sideband, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: point         ! Sky coordinates; (phi/theta/psi,time,det)
     real(sp),     allocatable, dimension(:,:,:)     :: mean_tp
     real(sp),     allocatable, dimension(:,:,:,:)   :: tod_poly      ! Poly-filter TOD coefficients (time,0:poly,sb,det)
     real(sp),     allocatable, dimension(:,:,:,:)   :: T_cont        ! Continuum temperature comp from freq filter (0:ncomp,time,sb,det)
     real(sp),     allocatable, dimension(:,:,:)     :: dg            ! gain fluctuation comp from freq filter (time,sb,det)
     real(sp),     allocatable, dimension(:,:,:)     :: var_fullres   ! Full-resolution variance (freq,sb,det)
     real(sp),     allocatable, dimension(:,:,:,:)   :: pca_ampl      ! amplitudes of pca-components (frec,sb,det,comp)
     real(sp),     allocatable, dimension(:,:)       :: pca_comp      ! actual pca component timestreams (time,comp)
     real(sp),     allocatable, dimension(:)         :: pca_eigv      ! eigenvalues of pca components (comp)

     real(sp),     allocatable, dimension(:, :, :, :)   :: pca_ampl_feed   ! amplitudes of pca-components (frec, sb, det, comp)
     real(sp),     allocatable, dimension(:, :, :)       :: pca_comp_feed     ! actual pca component timestreams (det, time, comp)
     real(sp),     allocatable, dimension(:, :)         :: pca_eigv_feed      ! eigenvalues of pca components (det, comp)

     real(sp),     allocatable, dimension(:,:)       :: acceptrate    ! fraction of freqs not masked (sb,det)
     real(sp),     allocatable, dimension(:,:,:,:)   :: diagnostics   ! various diagnostics used to make freqmask
     real(dp),     allocatable, dimension(:,:,:,:,:) :: spike_data    ! spike and jump data (n_spikes,spike/jump,info) info = (amp,mjd,samp,sb,feed)
     real(sp),     allocatable, dimension(:,:)       :: cut_params    ! means and stds used for the different diagnostics
     real(dp),     allocatable, dimension(:,:,:)     :: AB_mask       ! aliasing between A and B (in db suppression) (freq, sb,detector)
     real(dp),     allocatable, dimension(:,:,:)     :: leak_mask     ! aliasing from ouside range (in db suppression) (freq, sb,detector)

     real(dp),     allocatable, dimension(:,:,:)     :: corr_templ_ampl ! correlation template amplitudes (n_ampl, band, detector)

     real(dp),     allocatable, dimension(:,:,:)     :: sigma0, alpha, fknee ! (freq, nsb, detector)
     real(dp),     allocatable, dimension(:,:,:)     :: chi2          ! (freq, nsb, detector)
     real(sp),     allocatable, dimension(:,:,:)     :: gain                 ! (freq_fullres, nsb, detector)
     real(dp),     allocatable, dimension(:,:,:)     :: Tsys_lowres   ! (freq, sb,detector)
     real(sp),     allocatable, dimension(:,:,:,:,:) :: el_az_stats ! (g/a, n_chunks, freq, sb, feed)
     integer(lgt)                                    :: import_freqmask, import_sigma ! bool whether to import freqmask from existing l2 file
     
     ! Baseline Template file
     real(sp),    allocatable, dimension(:, :, :, :) :: tod_baseline, amplitudes ! (time, freq, sideband, detector) Baseline template of tod
     integer(sp),    allocatable, dimension(:)       :: Nperbaseline ! (time, freq, sideband, detector) Baseline template of tod
      
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

  subroutine read_l1_file(filename, data, id, verb, only_point, freqmask, init)
    implicit none
    character(len=*), intent(in)           :: filename
    integer(i4b),     intent(in)           :: id
    logical(lgt),     intent(in), optional :: only_point, init
    logical(lgt),     intent(in)           :: verb
    real(sp), dimension(:,:,:), intent(inout), optional :: freqmask  ! should not be used
    type(lx_struct)                        :: data
    type(hdf_file)                         :: file
    integer(i4b)                           :: i, j, k, l, nsamp, nsamp_tot, nfreq, ndet, num_masked
    integer(i4b)                           :: npoint, nsb, ext4(4), ext1(1), numbad, temp_samp, maxdet
    logical(lgt)                           :: all, ok, init_
    real(dp)                               :: t1, t2
    integer(i4b), allocatable, dimension(:)       :: buffer_int
    real(dp),     allocatable, dimension(:)       :: buffer_1d
    real(sp),     allocatable, dimension(:,:,:)   :: buffer_3d
    real(sp),     allocatable, dimension(:,:,:,:) :: buffer_4d
    all = .true.; if (present(only_point)) all = .not. only_point
    init_ = .true.; if (present(init)) init_ = init
    if (init_) call free_lx_struct(data) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "spectrometer/tod", ext4)
    nsamp_tot = ext4(1); nfreq = ext4(2) ; nsb = ext4(3); ndet = ext4(4)
    
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

    
    ! if (verb) then
    !    write(*,*) 'Warning: Adding 3 sec delay to amb_time!'
    ! end if
    ! data%amb_time = data%amb_time + 3.d0 /(24.d0*3600.d0)

    ! Read feed information
    call read_hdf(file, "spectrometer/feeds",               data%pixels)
    if (all) call read_hdf(file, "spectrometer/frequency",       data%nu(:,:,1))
    do i = 2, ndet
       data%nu(:,:,i) = data%nu(:,:,1)
    end do
    maxdet = maxval(data%pixels(:))
    allocate(data%pix2ind(maxdet))
    data%pix2ind(:) = 0
    do i = 1, ndet
       data%pix2ind(data%pixels(i)) = i
    end do

    ! Do elements that may have NaNs
    allocate(buffer_1d(nsamp_tot))
    if (all) allocate(buffer_int(nsamp_tot))
    if (all) allocate(buffer_4d(nsamp_tot,nfreq,nsb,ndet))
    if (all) allocate(buffer_3d(nsamp_tot,nsb,ndet))
    ! Find number of samples at end of file with NaNs
    if (all) call read_hdf(file, "spectrometer/tod",  buffer_4d)
    ! if (present(freqmask)) then
    !    allocate(data%n_nan(nfreq,nsb,ndet))
    !    ! Update frequency mask with channels that are all NaNs
    !    do i = 1, ndet
    !       do j = 1, nsb
    !          num_masked = 0
    !          do k = 1, nfreq
    !             if (freqmask(k,j,i) == 0.) cycle
    !             numbad = count(buffer_4d(:,k,j,i) .ne. buffer_4d(:,k,j,i))
    !             data%n_nan(k,j,i) = numbad
    !             if (numbad > 0.1*nsamp_tot) then
    !                num_masked = num_masked + 1
    !                if (verb) then
    !                   if (num_masked < 5) then
    !                      write(*,fmt='(a,i8,i6,i4,i8)') '  Removing frequency with >10% NaNs -- ', id, data%pixels(i), j, k
    !                   end if
    !                   if (num_masked == 4) then
    !                      write(*,*) "Suppressing NaN warnings for the rest of this sideband"
    !                   end if
    !                end if
    !                freqmask(k,j,i) = 0.
    !             end if
    !          end do
    !       end do
    !    end do

    !    ! Trim end for NaNs
    !    nsamp = nsamp_tot
    !    do while (nsamp > 0)
    !       ok = .true.
    !       do i = 1, ndet
    !          do j = 1, nsb
    !             do k = 1, nfreq
    !                if (freqmask(k,j,i) == 0.) cycle
    !                if (buffer_4d(nsamp,k,j,i) .ne. buffer_4d(nsamp,k,j,i)) then
    !                   ok = .false.
    !                   exit
    !                end if
    !             end do
    !          end do
    !       end do
    !       if (ok) then
    !          exit
    !       else
    !          nsamp = nsamp-1
    !       end if
    !    end do
    ! else
    !    nsamp = nsamp_tot
    ! end if
    allocate(data%n_nan(nfreq,nsb,ndet))
    data%n_nan = 0
    ! Update frequency mask with channels that are all NaNs
    do i = 1, ndet
       do j = 1, nsb
          num_masked = 0
          do k = 1, nfreq
             numbad = count(buffer_4d(:,k,j,i) .ne. buffer_4d(:,k,j,i))
             data%n_nan(k,j,i) = numbad
             
             if (numbad > 0.1*nsamp_tot) then
                buffer_4d(:,k,j,i) = 0.d0
                num_masked = num_masked + 1
                if (verb) then
                   if (num_masked < 5) then
                      write(*,fmt='(a,i8,i6,i4,i8)') '  Removing frequency with >10% NaNs -- ', id, data%pixels(i), j, k
                   end if
                   if (num_masked == 4) then
                      write(*,*) "Suppressing NaN warnings for the rest of this sideband"
                   end if
                end if
             end if
          end do
       end do
    end do

    ! Trim end for NaNs
    nsamp = nsamp_tot
    do while (nsamp > nsamp_tot - 500)
       ok = .true.
       do i = 1, ndet
          if (.not. is_alive(data%pixels(i))) cycle
          do j = 1, nsb
             do k = 1, nfreq
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

    if (verb) then
       if (nsamp /= nsamp_tot) then
          write(*,*) '  Number of NaN elements in ', id, ' = ', nsamp_tot-nsamp, ' of ', nsamp_tot
       end if
    end if
    allocate(data%time(nsamp))
    if (all) allocate(data%tod(nsamp,nfreq,nsb,ndet))
    if (all) allocate(data%tod_mean(nfreq,nsb,ndet))
    if (all) allocate(data%sb_mean(nsamp,nsb,ndet))
    if (all) allocate(data%flag(nsamp))
    
    call read_hdf(file, "spectrometer/band_average",  buffer_3d)
    data%sb_mean      = buffer_3d(1:nsamp,:,:)
    call read_hdf(file, "spectrometer/MJD",           buffer_1d)
    data%time         = buffer_1d(1:nsamp)
    if (all) data%tod = buffer_4d(1:nsamp,:,:,:)
    data%mjd_start    = minval(data%time)
    data%samprate     = 1.d0 / (3600.d0*24.d0*(data%time(2)-data%time(1)))
    data%flag         = 0
    
    ! tsys_stuff,  registers changed at time: 58712.03706
    
    if (data%time(1) > 58712.03706) then
       call get_size_hdf(file, "hk/antenna0/vane/Tvane", ext1)
    else
       call get_size_hdf(file, "hk/antenna0/env/ambientLoadTemp", ext1)
    end if
    temp_samp = ext1(1)

    allocate(data%t_hot(temp_samp))
    allocate(data%amb_state(temp_samp))
    allocate(data%amb_time(temp_samp))
  
    ! Read ambient temp
    if (data%time(1) > 58712.0370603) then
       call read_hdf(file, "hk/antenna0/vane/Tvane", data%t_hot)
    else
       call read_hdf(file, "hk/antenna0/env/ambientLoadTemp", data%t_hot)
    end if
    call read_hdf(file, "hk/antenna0/vane/state",          data%amb_state)
    call read_hdf(file, "hk/antenna0/vane/utc",            data%amb_time)

    call close_hdf_file(file)
    deallocate(buffer_1d)
    if (all) deallocate(buffer_4d,buffer_int)
  end subroutine read_l1_file


  subroutine read_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: nsamp, nfreq, nfreq_full, nsb, ndet, npoint, nsim, ext(4), poly
    call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "tod", ext)
    nsamp = ext(1); nfreq = ext(2) ; nsb = ext(3); ndet = ext(4)
    allocate(data%time(nsamp), data%tod(nsamp,nfreq,nsb,ndet), data%tod_mean(1024, nsb, ndet))
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
    allocate(data%freqmask_full(nfreq_full,nsb,ndet), data%freqmask(nfreq,nsb,ndet), data%mean_tp(nfreq_full,nsb,ndet), data%freqmask_reason(nfreq_full,nsb,ndet))
    call read_hdf(file, "freqmask",         data%freqmask)    
    call read_hdf(file, "freqmask_full",    data%freqmask_full)
    call read_hdf(file, "freqmask_reason",  data%freqmask_reason)
    call read_hdf(file, "mean_tp",          data%mean_tp)
    allocate(data%Tsys(nfreq_full,nsb,ndet), data%tsys_lowres(nfreq,nsb,ndet))
    call read_hdf(file, "Tsys",             data%Tsys)
    call read_hdf(file, "Tsys_lowres",      data%tsys_lowres)
    allocate(data%n_nan(nfreq_full,nsb,ndet))
    call read_hdf(file, "n_nan",            data%n_nan)

    call read_hdf(file, "polyorder",        data%polyorder)
    !if (data%polyorder >= 0) then
    !   allocate(data%tod_poly(nsamp,0:data%polyorder,nsb,ndet))
    !   call read_hdf(file, "tod_poly",         data%tod_poly)
    !end if
    call get_size_hdf(file, "pix2ind", ext)
    allocate(data%pixels(ndet), data%pix2ind(ext(1)))
    call read_hdf(file, "pixels",           data%pixels)
    call read_hdf(file, "pix2ind",          data%pix2ind)
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

    call read_hdf(file, "n_pca_comp_feed",         data%n_pca_comp_feed)
    !data%n_pca_comp_feed  = 0
    if (data%n_pca_comp_feed > 0) then
      allocate(data%pca_ampl_feed(nfreq_full, nsb, ndet, data%n_pca_comp_feed)) 
      allocate(data%pca_comp_feed(ndet, nsamp, data%n_pca_comp_feed))
      allocate(data%pca_eigv_feed(ndet, data%n_pca_comp_feed))
      call read_hdf(file, "pca_ampl_feed",         data%pca_ampl_feed)
      call read_hdf(file, "pca_comp_feed",         data%pca_comp_feed)
      call read_hdf(file, "pca_eigv_feed",         data%pca_eigv_feed)
    end if

    call read_hdf(file, "mask_outliers",    data%mask_outliers)
    !if (data%mask_outliers == 1) then
    
    if ((data%mask_outliers == 1) .and. (data%n_cal > 0)) then  ! if n_cal is 0 we don't run diagnostics
       
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

  subroutine read_baselines(filename, data)
  implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: nbasis, nfreq, nsb, ndet, ext(4)

    ! Reading in baseline fit of tod from file.

    call free_lx_struct(data)

    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "amplitudes", ext)

    nbasis = ext(1); nfreq = ext(2) ; nsb = ext(3); ndet = ext(4)

    allocate(data%Nperbaseline(nbasis))
    allocate(data%amplitudes(nbasis,nfreq,nsb,ndet))

    call read_hdf(file, "Nperbaseline", data%Nperbaseline)
    call read_hdf(file, "amplitudes", data%amplitudes)

    call close_hdf_file(file)
  end subroutine read_baselines

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
    if(allocated(data%AB_mask))       deallocate(data%AB_mask)
    if(allocated(data%leak_mask))     deallocate(data%leak_mask)
    if(allocated(data%corr_templ_ampl)) deallocate(data%corr_templ_ampl)
    if(allocated(data%gain))          deallocate(data%gain)
    if(allocated(data%sigma0))        deallocate(data%sigma0)
    if(allocated(data%alpha))         deallocate(data%alpha)
    if(allocated(data%Tsys))          deallocate(data%Tsys)
    if(allocated(data%Phot))          deallocate(data%Phot)
!    if(allocated(data%Pcold))         deallocate(data%Pcold)
    if(allocated(data%Thot))          deallocate(data%Thot)
    if(allocated(data%time_hot))      deallocate(data%time_hot)
    if(allocated(data%Tsys_lowres))   deallocate(data%Tsys_lowres)
    if(allocated(data%fknee))         deallocate(data%fknee)
    if(allocated(data%freqmask))      deallocate(data%freqmask)
    if(allocated(data%chi2))          deallocate(data%chi2)
    if(allocated(data%freqmask_full)) deallocate(data%freqmask_full)
    if(allocated(data%freqmask_reason)) deallocate(data%freqmask_reason)
    if(allocated(data%mean_tp))       deallocate(data%mean_tp)
    if(allocated(data%n_nan))         deallocate(data%n_nan)
    if(allocated(data%tod_poly))      deallocate(data%tod_poly)
    if(allocated(data%T_cont))        deallocate(data%T_cont)
    if(allocated(data%dg))            deallocate(data%dg)
    if(allocated(data%tod_mean))      deallocate(data%tod_mean)
    if(allocated(data%sb_mean))       deallocate(data%sb_mean)
    if(allocated(data%pixels))        deallocate(data%pixels)
    if(allocated(data%pix2ind))       deallocate(data%pix2ind)
    if(allocated(data%var_fullres))   deallocate(data%var_fullres)
    if(allocated(data%pca_ampl))      deallocate(data%pca_ampl)
    if(allocated(data%pca_comp))      deallocate(data%pca_comp)
    if(allocated(data%pca_eigv))      deallocate(data%pca_eigv)
    
    if(allocated(data%pca_ampl_feed))      deallocate(data%pca_ampl_feed)
    if(allocated(data%pca_comp_feed))      deallocate(data%pca_comp_feed)
    if(allocated(data%pca_eigv_feed))      deallocate(data%pca_eigv_feed)

    if(allocated(data%acceptrate))    deallocate(data%acceptrate)
    if(allocated(data%diagnostics))   deallocate(data%diagnostics)
    if(allocated(data%spike_data))    deallocate(data%spike_data)
    if(allocated(data%cut_params))    deallocate(data%cut_params)
    if(allocated(data%t_hot))         deallocate(data%t_hot)
    if(allocated(data%amb_state))     deallocate(data%amb_state)
    if(allocated(data%amb_time))      deallocate(data%amb_time)
    if(allocated(data%el_az_stats))   deallocate(data%el_az_stats)
    if(allocated(data%tod_baseline))   deallocate(data%tod_baseline)
    if(allocated(data%Nperbaseline))   deallocate(data%Nperbaseline)
    if(allocated(data%amplitudes))   deallocate(data%amplitudes)
  end subroutine

  subroutine write_l2_file(scan, k, data, name_append)
    implicit none
    type(comap_scan_info), intent(in) :: scan
    integer(i4b),     intent(in) :: k
    character(len=*), intent(in), optional :: name_append
    type(lx_struct)              :: data
    type(hdf_file)               :: file, l1_file
    integer(i4b)                 :: n_hk(1), hk_start_ind, hk_end_ind
    real(dp), allocatable, dimension(:) :: hk_buffer
    !write(*,*) "before all", name_append
    
    if (present(name_append)) then 
       call open_hdf_file(scan%ss(k)%l2file(1:len(trim(scan%ss(k)%l2file))-3)//adjustl(trim(name_append))//adjustl(trim(".h5")), file, "w")
    else
       call open_hdf_file(scan%ss(k)%l2file, file, "w")
    end if
    call write_hdf(file, "runID",             data%irun)
    call write_hdf(file, "samprate",          data%samprate)
    call write_hdf(file, "mjd_start",         data%mjd_start)
    call write_hdf(file, "decimation_time",   data%decimation_time)
    call write_hdf(file, "decimation_nu",     data%decimation_nu)
    call write_hdf(file, "time",              data%time)
    call write_hdf(file, "nu",                data%nu)
    call write_hdf(file, "tod",               data%tod)
    call write_hdf(file, "point_cel",         data%point_cel)
    call write_hdf(file, "point_tel",         data%point_tel)
    call write_hdf(file, "sb_mean",           data%sb_mean)
    call write_hdf(file, "tod_mean",          data%tod_mean)
    !call write_hdf(file, "flag",              data%flag)
    !write(*,*) "right before", data%Tsys(1, 1, 1, 1)
    call write_hdf(file, "cal_method",        data%cal_method)
    call write_hdf(file, "n_cal",             data%n_cal)
    call write_hdf(file, "Tsys",              data%Tsys)
    call write_hdf(file, "Thot",              data%Thot)
    call write_hdf(file, "Phot",              data%Phot)
!    call write_hdf(file, "Pcold",             data%Pcold)
    call write_hdf(file, "time_hot",          data%time_hot)
    call write_hdf(file, "Tsys_lowres",       data%Tsys_lowres)
    
    if (allocated(data%sigma0))    call write_hdf(file, "sigma0",            data%sigma0)
    if (allocated(data%alpha))     call write_hdf(file, "alpha",             data%alpha)
    if (allocated(data%fknee))     call write_hdf(file, "fknee",             data%fknee)   
    
    call write_hdf(file, "freqmask",          data%freqmask)
    call write_hdf(file, "freqmask_full",     data%freqmask_full)
    call write_hdf(file, "freqmask_reason",   data%freqmask_reason)
    
    call write_hdf(file, "n_nan",             data%n_nan)
    if (allocated(data%mean_tp)) call write_hdf(file, "mean_tp",           data%mean_tp)
    if (allocated(data%chi2)) call write_hdf(file, "chi2",           data%chi2)
    if (allocated(data%el_az_stats)) call write_hdf(file, "el_az_stats", data%el_az_stats)
    call write_hdf(file, "polyorder",         data%polyorder)
    
    if (data%use_freq_filter) then
       call write_hdf(file, "corr_templ_ampl", data%corr_templ_ampl)
       call write_hdf(file, "T_cont",         data%T_cont)
       call write_hdf(file, "dg",             data%dg)
    else if (data%polyorder >= 0) then
       call write_hdf(file, "tod_poly",         data%tod_poly)
    end if
       

    call write_hdf(file, "pixels",            data%pixels)
    call write_hdf(file, "pix2ind",           data%pix2ind)
    call write_hdf(file, "var_fullres",       data%var_fullres)
    call write_hdf(file, "n_pca_comp",        data%n_pca_comp)
    if (data%n_pca_comp > 0) then
       call write_hdf(file, "pca_ampl",       data%pca_ampl)
       call write_hdf(file, "pca_comp",       data%pca_comp)
       call write_hdf(file, "pca_eigv",       data%pca_eigv)
    end if
    
    call write_hdf(file, "n_pca_comp_feed",   data%n_pca_comp_feed)
    if (data%n_pca_comp_feed > 0) then
       call write_hdf(file, "pca_ampl_feed",       data%pca_ampl_feed)
       call write_hdf(file, "pca_comp_feed",       data%pca_comp_feed)
       call write_hdf(file, "pca_eigv_feed",       data%pca_eigv_feed)
    end if
    
    call write_hdf(file, "spike_data",        data%spike_data)
    call write_hdf(file, "mask_outliers",     data%mask_outliers)
    !write(*,*) "middle", data%mask_outliers
    
    if ((data%mask_outliers == 1) .and. (data%n_cal > 0)) then  ! if n_cal is 0 we don't run diagnostics
       !write(*,*) data%mask_outliers
       call write_hdf(file, "AB_aliasing",    data%AB_mask)
       call write_hdf(file, "leak_aliasing",  data%leak_mask)
       call write_hdf(file, "diagnostics",    data%diagnostics)
       call write_hdf(file, "cut_params",     data%cut_params)
       call write_hdf(file, "acceptrate",     data%acceptrate)
    end if 
    !write(*,*) "right after", data%Tsys(1, 1, 1, 1)
    
    ! scan-data
    ! call write_hdf(file, "l1file", teststr)
    call write_hdf(file, "field", scan%objectnum)
    call write_hdf(file, "time_of_day", scan%ss(k)%time_of_day)
    call write_hdf(file, "scanid", scan%ss(k)%id)
    call write_hdf(file, "feature", scan%ss(k)%feature)
    
    ! hk-data
    call open_hdf_file(scan%l1file, l1_file, "r")
    call get_size_hdf(l1_file, "hk/array/weather/utc", n_hk)
    allocate(hk_buffer(n_hk(1)))
    call read_hdf(l1_file, "hk/array/weather/utc", hk_buffer)
    
    hk_start_ind = minloc(abs(data%time(1) - hk_buffer), 1) + 1
    hk_end_ind = minloc(abs(data%time(size(data%time)) - hk_buffer), 1) - 1
    
    call write_hdf(file, "hk_mjd", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/airTemperature", hk_buffer)
    call write_hdf(file, "hk_airtemp", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/dewPointTemp", hk_buffer)
    call write_hdf(file, "hk_dewtemp", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/pressure", hk_buffer)
    call write_hdf(file, "hk_pressure", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/rainToday", hk_buffer)
    call write_hdf(file, "hk_rain", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/relativeHumidity", hk_buffer)
    call write_hdf(file, "hk_humidity", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/windDirection", hk_buffer)
    call write_hdf(file, "hk_winddir", hk_buffer(hk_start_ind:hk_end_ind))
    
    call read_hdf(l1_file, "hk/array/weather/windSpeed", hk_buffer)
    call write_hdf(file, "hk_windspeed", hk_buffer(hk_start_ind:hk_end_ind))
    
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
    lx_out%n_pca_comp_feed = lx_in%n_pca_comp_feed
    lx_out%mask_outliers = lx_in%mask_outliers
    lx_out%cal_method = lx_in%cal_method
    lx_out%n_cal = lx_in%n_cal
    lx_out%use_freq_filter = lx_in%use_freq_filter


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

    if(allocated(lx_in%tod_mean))         then  
       allocate(lx_out%tod_mean(size(lx_in%tod_mean,1),size(lx_in%tod_mean,2),size(lx_in%tod_mean,3)))
       lx_out%tod_mean = lx_in%tod_mean
    end if
    
    if(allocated(lx_in%corr_templ_ampl))         then  
       allocate(lx_out%corr_templ_ampl(size(lx_in%corr_templ_ampl,1),size(lx_in%corr_templ_ampl,2),size(lx_in%corr_templ_ampl,3)))
       lx_out%corr_templ_ampl = lx_in%corr_templ_ampl
    end if

    if(allocated(lx_in%sb_mean))         then  
       allocate(lx_out%sb_mean(size(lx_in%sb_mean,1),size(lx_in%sb_mean,2),size(lx_in%sb_mean,3)))
       lx_out%sb_mean = lx_in%sb_mean
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
    if(allocated(lx_in%chi2))      then
       allocate(lx_out%chi2(size(lx_in%chi2,1),size(lx_in%chi2,2),size(lx_in%chi2,3)))
       lx_out%chi2 = lx_in%chi2
    end if
    if(allocated(lx_in%Tsys_lowres))      then
       allocate(lx_out%Tsys_lowres(size(lx_in%Tsys_lowres,1),size(lx_in%Tsys_lowres,2),size(lx_in%Tsys_lowres,3)))
       lx_out%Tsys_lowres = lx_in%Tsys_lowres
    end if
    if(allocated(lx_in%Tsys))      then
       allocate(lx_out%Tsys(size(lx_in%Tsys,1),size(lx_in%Tsys,2),size(lx_in%Tsys,3)))
       lx_out%Tsys = lx_in%Tsys
    end if
    if(allocated(lx_in%Phot))      then
       allocate(lx_out%Phot(size(lx_in%Phot,1),size(lx_in%Phot,2),size(lx_in%Phot,3),size(lx_in%Phot,4),size(lx_in%Phot,5)))
       lx_out%Phot = lx_in%Phot
    end if
    ! if(allocated(lx_in%Pcold))      then
    !    allocate(lx_out%Pcold(size(lx_in%Pcold,1),size(lx_in%Pcold,2),size(lx_in%Pcold,3)))
    !    lx_out%Pcold = lx_in%Pcold
    ! end if
    if(allocated(lx_in%Thot))      then
       allocate(lx_out%Thot(size(lx_in%Thot,1),size(lx_in%Thot,2)))
       lx_out%Thot = lx_in%Thot
    end if
    if(allocated(lx_in%time_hot))      then
       allocate(lx_out%time_hot(size(lx_in%time_hot,1),size(lx_in%time_hot,2)))
       lx_out%time_hot = lx_in%time_hot
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
    if(allocated(lx_in%T_cont))      then
       allocate(lx_out%T_cont(size(lx_in%T_cont,1),size(lx_in%T_cont,2),size(lx_in%T_cont,3),size(lx_in%T_cont,4)))
       lx_out%T_cont = lx_in%T_cont
    end if
    if(allocated(lx_in%dg))      then
       allocate(lx_out%dg(size(lx_in%dg,1),size(lx_in%dg,2),size(lx_in%dg,3)))
       lx_out%dg = lx_in%dg
    end if
    if(allocated(lx_in%pixels))        then
       allocate(lx_out%pixels(size(lx_in%pixels,1)))
       lx_out%pixels = lx_in%pixels
    end if
    if(allocated(lx_in%pix2ind))        then
       allocate(lx_out%pix2ind(size(lx_in%pix2ind,1)))
       lx_out%pix2ind = lx_in%pix2ind
    end if
    if(allocated(lx_in%var_fullres))   then
       allocate(lx_out%var_fullres(size(lx_in%var_fullres,1),size(lx_in%var_fullres,2),size(lx_in%var_fullres,3)))
       lx_out%var_fullres = lx_in%var_fullres
    end if
    if(allocated(lx_in%AB_mask))   then
       allocate(lx_out%AB_mask(size(lx_in%AB_mask,1),size(lx_in%AB_mask,2),size(lx_in%AB_mask,3)))
       lx_out%AB_mask = lx_in%AB_mask
    end if
    if(allocated(lx_in%leak_mask))   then
       allocate(lx_out%leak_mask(size(lx_in%leak_mask,1),size(lx_in%leak_mask,2),size(lx_in%leak_mask,3)))
       lx_out%leak_mask = lx_in%leak_mask
    end if

    if(allocated(lx_in%spike_data))   then
       allocate(lx_out%spike_data(size(lx_in%spike_data,1),size(lx_in%spike_data,2),size(lx_in%spike_data,3),size(lx_in%spike_data,4),size(lx_in%spike_data,5)))
       lx_out%spike_data = lx_in%spike_data
    end if

    if(allocated(lx_in%el_az_stats))   then
       allocate(lx_out%el_az_stats(size(lx_in%el_az_stats,1),size(lx_in%el_az_stats,2),size(lx_in%el_az_stats,3),size(lx_in%el_az_stats,4),size(lx_in%el_az_stats,5)))
       lx_out%el_az_stats = lx_in%el_az_stats
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


    if(allocated(lx_in%pca_ampl_feed))   then
      allocate(lx_out%pca_ampl_feed(size(lx_in%pca_ampl_feed, 1), size(lx_in%pca_ampl_feed, 2), size(lx_in%pca_ampl_feed, 3), size(lx_in%pca_ampl_feed, 4)))
      lx_out%pca_ampl_feed = lx_in%pca_ampl_feed
   end if
   if(allocated(lx_in%pca_comp_feed))   then
      allocate(lx_out%pca_comp_feed(size(lx_in%pca_comp_feed, 1),size(lx_in%pca_comp_feed, 2), size(lx_in%pca_comp_feed, 3)))
      lx_out%pca_comp_feed = lx_in%pca_comp_feed
   end if
   if(allocated(lx_in%pca_eigv_feed))   then
      allocate(lx_out%pca_eigv_feed(size(lx_in%pca_eigv_feed, 1), size(lx_in%pca_eigv_feed, 2)))
      lx_out%pca_eigv_feed = lx_in%pca_eigv_feed
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

