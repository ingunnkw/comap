! Given a runlist, extract CESes from level1 files, producing level2-files
! in hdf format.
program l2gen
  use quiet_utils
  use quiet_fileutils
  use comap_scan_mod
  use comap_detector_mod
  use comap_Lx_mod
  use rngmod
  use quiet_fft_mod
  use spline_1D_mod
  use quiet_status_mod
  use comap_gain_mod
  use comap_patch_mod
  implicit none

  character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile, freqmaskfile, monitor_file_name, tsysfile
  character(len=9)     :: id_old
  integer(i4b)         :: i, j, k, l, m, n, snum, nscan, unit, myid, nproc, ierr, ndet, npercore
  integer(i4b)         :: mstep, i2, decimation, nsamp, numfreq, n_nb
  integer(i4b)         :: debug, num_l1_files, seed, bp_filter, bp_filter0, n_pca_comp, pca_max_iter
  real(dp)             :: todsize, nb_factor, min_acceptrate
  real(dp)             :: pca_err_tol, corr_cut, mean_corr_cut, mean_abs_corr_cut, med_cut, var_cut
  logical(lgt)         :: exist, reprocess, check_existing, gonext, found, rm_outliers, mask_outliers
  logical(lgt)         :: process
  real(dp)             :: timing_offset, mjd(2), dt_error, samprate_in, samprate, scanfreq, nu_gain, alpha_gain, t1, t2
  type(comap_scan_info) :: scan
  type(Lx_struct)      :: data_l1, data_l2_fullres, data_l2_decimated, data_l2_filter
  type(planck_rng)     :: rng_handle
  type(status_file)    :: status
  type(patch_info)     :: pinfo
  !real(dp),     allocatable, dimension(:,:,:,:)   :: store_l2_tod

  call getarg(1, parfile)
  call get_parameter(unit, parfile, 'TSYS_LOC',                  par_string=tsysfile)
  call get_parameter(unit, parfile, 'L2_SAMPRATE',               par_dp=samprate)
  call get_parameter(unit, parfile, 'NUMFREQ',                   par_int=numfreq)
  call get_parameter(unit, parfile, 'REPROCESS_ALL_FILES',       par_lgt=reprocess)
  call get_parameter(unit, parfile, 'DEBUG',                     par_int=debug)
  call get_parameter(unit, parfile, 'SEED',                      par_int=seed)
  call get_parameter(unit, parfile, 'FREQUENCY_MASK',            par_string=freqmaskfile)
  call get_parameter(unit, parfile, 'GAIN_NORMALIZATION_NU',     par_dp=nu_gain)
  call get_parameter(unit, parfile, 'GAIN_NORMALIZATION_ALPHA',  par_dp=alpha_gain)
  call get_parameter(unit, parfile, 'BANDPASS_FILTER_ORDER',     par_int=bp_filter0)
  call get_parameter(unit, parfile, 'N_PCA_COMPONENTS',          par_int=n_pca_comp)
  call get_parameter(unit, parfile, 'PCA_ERROR_TOLERANCE',       par_dp=pca_err_tol)
  call get_parameter(unit, parfile, 'PCA_MAX_ITERATIONS',        par_int=pca_max_iter)
  call get_parameter(unit, parfile, 'MASK_OUTLIERS',             par_lgt=mask_outliers)
  call get_parameter(unit, parfile, 'CORRELATION_CUT',           par_dp=corr_cut)
  call get_parameter(unit, parfile, 'MEAN_CORRELATION_CUT',      par_dp=mean_corr_cut)
  call get_parameter(unit, parfile, 'VARIANCE_CUT',              par_dp=var_cut)
  call get_parameter(unit, parfile, 'MEAN_ABS_CORRELATION_CUT',  par_dp=mean_abs_corr_cut)
  call get_parameter(unit, parfile, 'MEDIAN_CUT',                par_dp=med_cut)
  call get_parameter(unit, parfile, 'N_NEIGHBOR',                par_int=n_nb)
  call get_parameter(unit, parfile, 'NEIGHBOR_FACTOR',           par_dp=nb_factor)
  call get_parameter(unit, parfile, 'MIN_ACCEPTRATE',            par_dp=min_acceptrate)
  call get_parameter(unit, parfile, 'REMOVE_OUTLIERS',           par_lgt=rm_outliers)
  call get_parameter(unit, parfile, 'LEVEL2_DIR',                par_string=l2dir)
    

  check_existing = .true.
  call mkdirs(trim(l2dir), .false.)
  call initialize_scan_mod(parfile)
  call initialize_detector_mod(parfile)
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  call init_status(status, trim(l2dir)//'/l2gen_mon.txt'); 
  call update_status(status, 'init')
  call initialize_comap_patch_mod(parfile)
  !call dset(id=myid,level=debug)


  call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)

  nscan    = get_num_scans()
  do snum = myid+1, nscan, nproc
     call get_scan_info(snum, scan)
     process = .false.
     do k = 1, scan%nsub
        inquire(file=scan%ss(k)%l2file,exist=exist)
        if (reprocess .and. exist) call rm(scan%ss(k)%l2file)           
        if (reprocess .or. .not. exist) process = .true.
     end do
     if (.not. process) then
        write(*,fmt="(i3,a,2i5,i8)") myid, " skipping already finished scan:", snum, scan%id
        cycle
     end if

     write(*,fmt="(i3,a,i10,a)") myid, " processing scan ", scan%id, " (" // trim(itoa(snum)) // "/" // trim(itoa(nscan)) // ")"
     call update_status(status, 'scan_start')

     ! Initialize frequency mask
     call initialize_fullres_frequency_mask(freqmaskfile, data_l1)
     call update_status(status, 'freq_mask1')

     ! Read in Level 1 file
     call wall_time(t1)
     call read_l1_file(scan%l1file, data_l1, scan%id, freqmask=data_l1%freqmask_full, init=.false.)
     call update_status(status, 'read_l1')
     if (size(data_l1%tod,1) <100) then
        write(*,*) 'Too few samples in ', scan%id
        cycle
     end if
     !call correct_missing_time_steps(data_l1%time)
     call wall_time(t2)
     todsize = real(size(data_l1%tod,1),dp)*real(size(data_l1%tod(1,:,:,:)),dp)*4.d0/1024.d0**2 ! in MB
     write(*,fmt='(i4,a,f10.2,a)') myid, ' -- disk read speed = ', real(todsize,dp)/(t2-t1), ' MB/sec'

     ! Interpolate over single NaNs
     call interpolate_nans(data_l1, scan%id)
     call update_status(status, 'nan_interp')

     ! Finalize frequency mask
     call postprocess_frequency_mask(numfreq, data_l1, scan%id)
     call update_status(status, 'freq_mask2')

     ! Compute absolute calibration
     call compute_Tsys(tsysfile, data_l1)

     ! Get patch info
     found = get_patch_info(scan%object, pinfo) 
     if (.not. found) then
        write(*,*) 'Error: Patch not found in patchfile = ', trim(scan%object)
        call mpi_finalize(ierr)
        stop
     end if

     do k = 1, scan%nsub

        ! Reformat L1 data into L2 format, and truncate
        call excise_subscan(scan%ss(k)%mjd, data_l1, data_l2_fullres)
        call update_status(status, 'excise')

        ! Normalize gain
        if (trim(pinfo%type) == 'gal' .or. trim(pinfo%type) == 'cosmo') then
           call normalize_gain(data_l2_fullres, nu_gain, alpha_gain, scan%id)
           call update_status(status, 'gain_norm')
        end if
     
        ! Apply absolute calibration
        call calibrate_tod(data_l1, data_l2_fullres)


        if (mask_outliers) then
           ! Copy tod, run filtering, make new mask, then do filtering again on original (unfiltered) data
           call copy_lx_struct(data_l2_fullres, data_l2_filter)
           
           ! Poly-filter copied data
           call polyfilter_TOD(data_l2_filter, bp_filter0)
           call update_status(status, 'polyfilter0')
           
           ! pca filter copied data
           call pca_filter_TOD(data_l2_filter, n_pca_comp, pca_max_iter, pca_err_tol)
           call update_status(status, 'pca_filter0')
           
           ! flag correlations and variance
           call flag_correlations(data_l2_filter, scan%id, corr_cut, mean_corr_cut, mean_abs_corr_cut, med_cut, var_cut, n_nb, nb_factor)
           call update_status(status, 'flag_corr')

           ! replace freqmask in original tod
           data_l2_fullres%freqmask_full = data_l2_filter%freqmask_full 
           
           allocate(data_l2_fullres%diagnostics(&
                & size(data_l2_filter%diagnostics,1), &
                & size(data_l2_filter%diagnostics,2), &
                & size(data_l2_filter%diagnostics,3), &
                & size(data_l2_filter%diagnostics,4)))
           data_l2_fullres%diagnostics = data_l2_filter%diagnostics

           allocate(data_l2_fullres%cut_params(&
                & size(data_l2_filter%cut_params,1), &
                & size(data_l2_filter%cut_params,2)))
           data_l2_fullres%cut_params = data_l2_filter%cut_params

           call update_freqmask(data_l2_fullres, min_acceptrate, scan%id)
           call update_status(status, 'made_freqmask')

           call free_lx_struct(data_l2_filter)
        end if
     
        ! Poly-filter if requested
        bp_filter = -1; if (trim(pinfo%type) == 'cosmo') bp_filter = bp_filter0
        call polyfilter_TOD(data_l2_fullres, bp_filter)
        call update_status(status, 'polyfilter')
        
        ! pca filter after polyfilter
        if (trim(pinfo%type) .ne. 'cosmo') n_pca_comp = 0
        call pca_filter_TOD(data_l2_fullres, n_pca_comp, pca_max_iter, pca_err_tol)
        call update_status(status, 'pca_filter')
        
        
        ! Fourier transform frequency direction
        !call convert_GHz_to_k(data_l2_fullres(i))
        
        ! If necessary, decimate L2 file in both time and frequency
        call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
        call update_status(status, 'decimate')
        !write(*,*) 'c'

        ! Fit noise
        call fit_noise(data_l2_decimated)

        ! Replace TOD with simulated data
        if (.false.) call simulate_gain_data(rng_handle, data_l2_decimated)

        ! Write L2 file to disk
        write(*,*) 'Writing ', scan%id, ' to disk', trim(scan%ss(k)%l2file)
        call mkdirs(trim(scan%ss(k)%l2file), .true.)
        call write_l2_file(scan%ss(k)%l2file, data_l2_decimated)
        call update_status(status, 'write_l2')

        ! Clean up data structures
        call free_lx_struct(data_l2_decimated)
        call free_lx_struct(data_l2_fullres)
     
     end do        

  end do
  call mpi_finalize(ierr)

  call free_status(status)

contains

  subroutine update_freqmask(data_l2, min_acceptrate, id)
    implicit none
    type(Lx_struct),                            intent(inout) :: data_l2
    integer(i4b),                               intent(in)    :: id
    real(dp),                                   intent(in)    :: min_acceptrate
    integer(i4b) :: i, j, k, l, m, n, nsamp, nfreq, nfreq_full, nsb, ndet, dfreq
    
    
    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%freqmask,1) ! nfreq in lowres data
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    nfreq_full  = size(data_l2%freqmask_full,1)
    dfreq       = nfreq_full / nfreq
    if(.not. allocated(data_l2%acceptrate)) allocate(data_l2%acceptrate(nsb,ndet)) 

    ! calculate acceptrate
    do i = 1, ndet
       do j = 1, nsb
          data_l2%acceptrate(j,i) = sum(data_l2%freqmask_full(:,j,i)) / nfreq_full
          if (data_l2%acceptrate(j,i) < min_acceptrate) then !Mask bad sidebands
             data_l2%freqmask_full(:,j,i) = 0.d0
             write(*,*) "Rejecting entire sideband (too much was masked) det, sb, acceptrate, scanid:"
             write(*,*) i, j, data_l2%acceptrate(j,i), id
             data_l2%acceptrate(j,i) = 0.d0
          end if
       end do
    end do
    

    ! update lowres freqmask if all frequencies within a lowres bin are masked
    do k = 1, ndet
       do j = 1, nsb
          do i = 1, nfreq
             if (all(data_l2%freqmask_full((i-1)*dfreq+1:i*dfreq,j,k) == 0.d0)) data_l2%freqmask(i,j,k) = 0.d0
          end do
       end do
    end do
    
    
    
  end subroutine update_freqmask

  subroutine flag_correlations(data_l2, id, corr_cut, mean_corr_cut, mean_abs_corr_cut, median_cut, var_cut, n_neighbor, neighbor_factor)
    implicit none
    type(Lx_struct),          intent(inout) :: data_l2
    integer(i4b),             intent(in)    :: id
    real(dp),                 intent(in)    :: corr_cut, mean_corr_cut, var_cut, median_cut
    real(dp),                 intent(in)    :: mean_abs_corr_cut, neighbor_factor
    integer(i4b),             intent(in)    :: n_neighbor
    integer(i4b) :: i, j, k, l, m, n, nsamp, nfreq, nsb, ndet
    real(dp)     :: corr, mean, var_0, std_median
    real(dp)     :: mean_meancorr, sigma_meancorr, mean_vars, sigma_vars, mean_maxcorr, sigma_maxcorr
    real(dp)     :: mean_meanabscorr, sigma_meanabscorr, mean_corrsum, sigma_corrsum
    real(dp),     allocatable, dimension(:, :, :)   :: means, vars, maxcorr, meancorr, meanabscorr
    real(dp),     allocatable, dimension(:, :, :)   :: corrsum
    real(sp),     allocatable, dimension(:, :, :)   :: corrsum_mask
    real(dp),     allocatable, dimension(:)         :: subt, median, smedian, minimum
    real(dp),     allocatable, dimension(:,:)       :: corrs
    real(dp),     allocatable, dimension(:,:)       :: outlier_mask
    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    
    allocate(means(nfreq, nsb, ndet), vars(nfreq, nsb, ndet))
    allocate(maxcorr(nfreq, nsb, ndet), meancorr(nfreq, nsb, ndet))
    allocate(meanabscorr(nfreq, nsb, ndet))
    if(.not. allocated(data_l2%diagnostics)) allocate(data_l2%diagnostics(nfreq,nsb,ndet,5)) 
    if(.not. allocated(data_l2%cut_params)) allocate(data_l2%cut_params(2,5)) 
    allocate(corrs(nfreq,nfreq))
    allocate(subt(nsamp-1), median(nfreq))
    
    vars = 0.d0
    maxcorr = 0.d0
    meancorr = 0.d0
    meanabscorr = 0.d0
    
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             means(k,j,i) = sum(data_l2%tod(:,k,j,i)) / nsamp
             vars(k,j,i) = sum(data_l2%tod(:,k,j,i) ** 2) / nsamp - means(k,j,i) ** 2
          end do
       end do
    end do
    
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          corrs = 0.d0
          !$OMP PARALLEL PRIVATE(k,l,m,n,corr)
          !$OMP DO SCHEDULE(guided)    
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             !do l = i, ndet
             !   if (.not. is_alive(l)) cycle
             !   do m = j, nsb
             l = i
             m = j
             do n = k + 1, nfreq 
                if (data_l2%freqmask_full(n,m,l) == 0.d0) cycle
                corr = sum((data_l2%tod(:,k,j,i) - means(k,j,i)) * (data_l2%tod(:,n,m,l) - means(n,m,l))) / nsamp
                corr = corr / sqrt(vars(k,j,i) * vars(n,m,l))
                corrs(k,n) = corr
                corrs(n,k) = corr
             end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             meancorr(k,j,i)    = sum(corrs(k,:)) / (sum(data_l2%freqmask_full(:,j,i))-1.d0)
             maxcorr(k,j,i)     = maxval(abs(corrs(k,:)))
             meanabscorr(k,j,i) = sum(abs(corrs(k,:))) / (sum(data_l2%freqmask_full(:,j,i))-1.d0)
          end do
       end do
       !write(*,*) "Done with correlation masking in detector", i
    end do
    
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             subt = (data_l2%tod(2:nsamp,k,j,i) - data_l2%tod(1:nsamp-1,k,j,i)) / sqrt(2.d0)
             mean = sum(subt) / (nsamp - 1)
             var_0 = sum(subt ** 2) / (nsamp - 1) - mean ** 2
             vars(k,j,i) = vars(k,j,i) / var_0
          end do
       end do
    end do
    call get_mean_and_sigma(maxcorr, data_l2%freqmask_full, mean_maxcorr, sigma_maxcorr, .true.)
    call get_mean_and_sigma(meancorr, data_l2%freqmask_full, mean_meancorr, sigma_meancorr, .true.)
    call get_mean_and_sigma(meanabscorr, data_l2%freqmask_full, mean_meanabscorr, sigma_meanabscorr, .true.)
    call get_mean_and_sigma(vars, data_l2%freqmask_full, mean_vars, sigma_vars, .true.)
    
    median = 0.0055d0
    std_median = 0.0001d0 
    deallocate(subt)
    allocate(subt(nfreq-1))
    
    
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             median(k) = median(k) + 0.00005 * sign(1.d0, meanabscorr(k,j,i) - median(k))
          end do
          subt = (meanabscorr(2:nfreq,j,i) - meanabscorr(1:nfreq-1,j,i)) / sqrt(2.d0)
          var_0 = sum(merge(subt ** 2, 0.d0, ((data_l2%freqmask_full(2:nfreq,j,i) == 1) .and. (data_l2%freqmask_full(1:nfreq-1,j,i) == 1))))
          if (sum(merge(1.d0, 0.d0, ((data_l2%freqmask_full(2:nfreq,j,i) == 1) .and. (data_l2%freqmask_full(1:nfreq-1,j,i) == 1)))) == 0.d0) then
             write(*,*) "OVERFLOW AVERTED", i, j, id
             var_0 = 0.d0
          else 
             var_0 = var_0 / sum(merge(1.d0, 0.d0, ((data_l2%freqmask_full(2:nfreq,j,i) == 1) .and. (data_l2%freqmask_full(1:nfreq-1,j,i) == 1))))
          end if
          std_median = std_median + 0.00002 * sign(1.d0, sqrt(var_0) - std_median)
       end do
    end do
    ! Mask outlier frequencies
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          do k = 1, nfreq
             ! cut on neighbours
             if (k < nfreq - n_neighbor) then
                if (all(var_cut * sigma_vars * neighbor_factor < abs(vars(k:k+n_neighbor,j,i) - mean_vars))) then
                   data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                else if (all(corr_cut * sigma_maxcorr * neighbor_factor < abs(maxcorr(k:k+n_neighbor,j,i) - mean_maxcorr))) then
                   data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                else if (all(mean_corr_cut * sigma_meancorr * neighbor_factor < abs(meancorr(k:k+n_neighbor,j,i) - mean_meancorr))) then
                   data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                else if (all(mean_abs_corr_cut * sigma_meanabscorr * neighbor_factor < abs(meanabscorr(k:k+n_neighbor,j,i) - mean_meanabscorr))) then
                   data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                else if (all(median_cut * std_median * neighbor_factor < meanabscorr(k:k+n_neighbor,j,i) - median(k))) then
                   data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                end if
             end if
             
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             ! cut on single outliers
             if (var_cut * sigma_vars < abs(vars(k,j,i) - mean_vars)) then
                data_l2%freqmask_full(k,j,i) = 0.d0
             else if (corr_cut * sigma_maxcorr < abs(maxcorr(k,j,i) - mean_maxcorr)) then
                data_l2%freqmask_full(k,j,i) = 0.d0
             else if (mean_corr_cut * sigma_meancorr < abs(meancorr(k,j,i) - mean_meancorr)) then
                data_l2%freqmask_full(k,j,i) = 0.d0
             else if (mean_abs_corr_cut * sigma_meanabscorr < abs(meanabscorr(k,j,i) - mean_meanabscorr)) then
                data_l2%freqmask_full(k,j,i) = 0.d0
             else if (median_cut * std_median < meanabscorr(k,j,i) - median(k)) then
                data_l2%freqmask_full(k,j,i) = 0.d0
             end if                       
          end do
       end do
    end do
    
    data_l2%diagnostics(:,:,:,1) = means
    data_l2%diagnostics(:,:,:,2) = vars
    data_l2%diagnostics(:,:,:,3) = maxcorr
    data_l2%diagnostics(:,:,:,4) = meancorr
    data_l2%diagnostics(:,:,:,5) = meanabscorr
    data_l2%cut_params(1,1) = mean_vars
    data_l2%cut_params(2,1) = sigma_vars
    data_l2%cut_params(1,2) = mean_maxcorr
    data_l2%cut_params(2,2) = sigma_maxcorr
    data_l2%cut_params(1,3) = mean_meancorr
    data_l2%cut_params(2,3) = sigma_meancorr
    data_l2%cut_params(1,4) = mean_meanabscorr
    data_l2%cut_params(2,4) = sigma_meanabscorr
    data_l2%cut_params(1,5) = sum(median) / nfreq
    data_l2%cut_params(2,5) = std_median
    deallocate(vars)
    deallocate(means)
    deallocate(median)
    deallocate(corrs)
    deallocate(subt)
    deallocate(maxcorr)
    deallocate(meancorr)
    deallocate(meanabscorr)  
  end subroutine flag_correlations

  subroutine get_mean_and_sigma(data, mask, mean, sigma, remove_outliers)
    real(dp), dimension(:,:,:), intent(in)  :: data
    real(sp), dimension(:,:,:), intent(in)  :: mask
    real(dp),                   intent(out) :: mean, sigma
    logical(lgt),               intent(in)  :: remove_outliers
    real(dp), allocatable, dimension(:,:,:) :: outlier_mask
    
    mean = sum(data) / sum(mask)
    sigma  = sqrt(sum(data ** 2) / sum(mask) - mean ** 2)
    
    if (remove_outliers) then
       allocate(outlier_mask(size(data,1),size(data,2),size(data,3)))
       outlier_mask = mask * merge(1.d0,0.d0,abs(data - mean) <= 3.0 * sigma)
       mean = sum(data * outlier_mask) / sum(outlier_mask)
       sigma = sqrt(sum((data * outlier_mask) ** 2) / sum(outlier_mask) - mean ** 2)
       deallocate(outlier_mask)
    end if
  end subroutine get_mean_and_sigma
  
  subroutine interpolate_nans(data, id)
    implicit none
    type(Lx_struct),  intent(inout) :: data
    integer(i4b),     intent(in)    :: id

    integer(i4b) :: i, j, k, l, nsamp, nfreq, nsb, ndet, n, ntot

    nsamp       = size(data%tod,1)
    nfreq       = size(data%tod,2)
    nsb         = size(data%tod,3)
    ndet        = size(data%tod,4)

    ntot = 0
    !$OMP PARALLEL PRIVATE(k,l,i,j,n)
    n = 0
    !$OMP DO SCHEDULE(guided)
    do k = 1, nfreq
       do l = 2, nsamp-1
          do i = 1, ndet
             do j = 1, nsb
                if (data%tod(l,k,j,i) .ne. data%tod(l,k,j,i)) then
                   ! Only interpolate over single missing NaNs, not multiple in a row; these must be dealt with
                   ! in other ways
                   if (data%tod(l-1,k,j,i) .eq. data%tod(l-1,k,j,i) .and. &
                        & data%tod(l+1,k,j,i) .eq. data%tod(l+1,k,j,i)) then
                      data%tod(l,k,j,i) = 0.5 * (data%tod(l-1,k,j,i)+data%tod(l+1,k,j,i))
                      n                 = n+1
                   end if
                end if
             end do
          end do
       end do
    end do
    !$OMP END DO
    !$OMP ATOMIC
    ntot = ntot + n
    !$OMP END PARALLEL
    if (ntot > 0) write(*,*) '  Interpolated over NaN samples in ', id, ', n_tot = ', ntot

  end subroutine interpolate_nans

  subroutine normalize_gain(data_l2, nu_gain, alpha_gain, id)
    implicit none
    type(Lx_struct),                            intent(inout) :: data_l2
    real(dp),                                   intent(in)    :: nu_gain, alpha_gain
    integer(i4b),                               intent(in)    :: id

    integer(i4b) :: i, j, k, l, nomp, nsamp, nfreq, nsb, ndet, err
    integer*8    :: plan_fwd, plan_back
    real(dp)     :: samprate, nu
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv

    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    samprate    = data_l2%samprate    
    n           = nsamp+1

    if (nsamp == 0) then
       allocate(data_l2%mean_tp(nfreq,nsb,ndet))
       data_l2%mean_tp = 0.d0
       return
    end if

    ! Set up OpenMP environment and FFTW plans
    nomp = 1
    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)

    allocate(dt(2*nsamp), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*nsamp, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, 2*nsamp, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)

    allocate(data_l2%mean_tp(nfreq,nsb,ndet))
    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       !write(*,*) '    Normalizing gains for det = ', i
       do j = 1, nsb
          !$OMP PARALLEL PRIVATE(k,l,dt,dv,nu)
          allocate(dt(2*nsamp), dv(0:n-1))
          !$OMP DO SCHEDULE(guided)
          do k = 1, nfreq
             if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
             dt(1:nsamp)            = data_l2%tod(:,k,j,i)
             dt(2*nsamp:nsamp+1:-1) = dt(1:nsamp)
             !call fft(dt, dv, 1)
             call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
             ! Apply lowpass filter
             do l = 0, n-1
                nu = ind2freq(l+1, samprate, n)
                dv(l) = dv(l) * 1.d0/(1.d0 + (nu/nu_gain)**alpha_gain)
             end do

             !call fft(dt, dv, -1)
             call sfftw_execute_dft_c2r(plan_back, dv, dt)
             dt                       = dt / (2*nsamp)
             data_l2%tod(:,k,j,i)     = data_l2%tod(:,k,j,i) / dt(1:nsamp) - 1.d0 
             data_l2%mean_tp(k,j,i)   = mean(dt(1:nsamp))
          end do
          !$OMP END DO
          deallocate(dt, dv)
          !$OMP END PARALLEL
       end do
    end do

    call sfftw_destroy_plan(plan_fwd)
    call sfftw_destroy_plan(plan_back)

  end subroutine normalize_gain
  
  subroutine pca_filter_TOD(data_l2, n_pca_comp, pca_max_iter, pca_err_tol)
    implicit none
    type(Lx_struct),           intent(inout) :: data_l2
    integer(i4b),                 intent(in)    :: n_pca_comp, pca_max_iter
    real(dp),                  intent(in)    :: pca_err_tol
    integer(i4b) :: i, j, k, l, nsamp, nfreq, nsb, ndet, stat, iters
    real(dp)     :: eigenv, dotsum, amp, err
    real(sp),     allocatable, dimension(:)   :: r, s
    CHARACTER(LEN=128) :: number
    
    data_l2%n_pca_comp = n_pca_comp
    if (n_pca_comp == 0) return

    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)

    allocate(r(nsamp), s(nsamp))
            
    if(.not. allocated(data_l2%pca_ampl)) allocate(data_l2%pca_ampl(nfreq,nsb,ndet,n_pca_comp)) 
    if(.not. allocated(data_l2%pca_comp)) allocate(data_l2%pca_comp(nsamp,n_pca_comp))
    if(.not. allocated(data_l2%pca_eigv)) allocate(data_l2%pca_eigv(n_pca_comp))
    
    do l = 1, n_pca_comp 
       err = 1.d0
       r(:) = data_l2%tod(:,5,2,5)
       iters = 0
       do while ((err > pca_err_tol) .and. (iters < pca_max_iter))
          s = 0.d0
          do i = 1, ndet
             if (.not. is_alive(i)) cycle
             do j = 1, nsb
                do k = 1, nfreq
                   if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle   
                   dotsum = sum(data_l2%tod(:,k,j,i) * r(:))
                   !dotsum = sdot(nfreq,data_l2%tod(:,k,j,i), 1, r(:), 1)
                   s(:) = s(:) + dotsum * data_l2%tod(:,k,j,i)
                end do
             end do
          end do
          
          eigenv = sum(s(:) * r(:))
          
          err = sqrt(sum((eigenv * r(:) - s(:)) ** 2))
          r(:) = s(:)/sqrt(sum(s(:) ** 2))
          iters = iters + 1
       end do
       data_l2%pca_eigv(l) = eigenv
       data_l2%pca_comp(:,l) = r(:)
       do i = 1, ndet
          if (.not. is_alive(i)) cycle
          do j = 1, nsb
             do k = 1, nfreq
                if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
                data_l2%pca_ampl(k,j,i,l) = sum(r(:) * data_l2%tod(:,k,j,i))
                data_l2%tod(:,k,j,i) = data_l2%tod(:,k,j,i) - data_l2%pca_ampl(k,j,i,l) * r(:)
             end do
          end do
       end do
    end do
    deallocate(r, s)

  end subroutine pca_filter_TOD

  subroutine polyfilter_TOD(data_l2, bp_filter)
    implicit none
    type(Lx_struct),                            intent(inout) :: data_l2
    integer(i4b),                               intent(in)    :: bp_filter

    integer(i4b) :: i, j, k, l, n, nsamp, nfreq, nsb, ndet, p, stat
    real(dp)     :: samprate, nu, mu
    real(dp),     allocatable, dimension(:)   :: x
    real(dp),     allocatable, dimension(:,:) :: T, A

    data_l2%polyorder = bp_filter
    if (bp_filter < 0) return

    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    p           = bp_filter

    if (nsamp == 0) then
       allocate(data_l2%tod_poly(nsamp,0:p,nsb,ndet))
       data_l2%tod_poly = 0.d0
       return
    end if

    allocate(T(nfreq,0:p), A(0:p,0:p), x(0:p))
    allocate(data_l2%tod_poly(nsamp,0:p,nsb,ndet))

    ! Precompute polynomial basis
    do k = 1, nfreq
       mu = max(min(2.d0*real(k-1,dp)/real(nfreq-1,dp)-1.d0,1.d0),-1.d0)
       call get_legendre_polynomials(mu,T(k,:))
    end do

    do i = 1, ndet
       if (.not. is_alive(i)) cycle
       do j = 1, nsb
          !write(*,*) 'Polyfiltering det, sb = ', i, j

          if (all(data_l2%freqmask_full(:,j,i) == 0.d0)) cycle

          ! Pre-compute Cholesky factor of coupling matrix for current sideband
          do m = 0, p
             do n = 0, m
                A(m,n) = sum(T(:,m)*T(:,n)*data_l2%freqmask_full(:,j,i))
             end do
          end do
          call dpotrf('L', p+1, A, p+1, stat )

          ! Solve for polynomial coefficients
          !$OMP PARALLEL PRIVATE(k,m,x,l,stat)
          !$OMP DO SCHEDULE(guided)
          do k = 1, nsamp
             do m = 0, p
                x(m) = 0.d0
                do l = 1, nfreq
                   if (data_l2%freqmask_full(l,j,i) == 0.d0) cycle
                   x(m) = x(m) + data_l2%tod(k,l,j,i)*T(l,m)
                end do
             end do
             call dpotrs( 'L', p+1, 1, A, p+1, x, p+1, stat)

             ! Store poly coeffs as separate TOD
             data_l2%tod_poly(k,:,j,i) = x

             ! Subtract fitted polynomial from main TOD
             do l = 1, nfreq
                if (data_l2%freqmask_full(l,j,i) == 0.d0) cycle
                data_l2%tod(k,l,j,i) = data_l2%tod(k,l,j,i) - sum(T(l,:)*x)
             end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
       end do
    end do
    deallocate(T, A, x)

  end subroutine polyfilter_TOD


  subroutine get_l2_time_stats(filename, mjd, dt_error)
    implicit none
    character(len=*) :: filename
    real(dp)         :: mjd(2), dt_ideal, dt_error, srate
    integer(i4b)     :: n(7)
    type(hdf_file)   :: file
    real(dp),     dimension(:), allocatable :: time
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "time", n)
    allocate(time(n(1)))
    call read_hdf(file, "time", time)
    call read_hdf(file, "samprate", srate)
    dt_ideal = 1000/srate
    dt_error = maxval(abs((time(2:)-time(1:size(time)-1))*24*60*60*1000 - dt_ideal))/dt_ideal
    mjd = [ time(1), time(size(time)) ]
    call close_hdf_file(file)
    deallocate(time)
  end subroutine

  subroutine compute_time_offset(time_el, el, time, tod, offset)
    implicit none
    real(dp), dimension(1:), intent(in)  :: time_el, el, time, tod
    real(dp),                intent(out) :: offset 

    integer(i4b) :: i, j, p, q, n_tod, n_el, start_tod, end_tod, start_el, end_el
    integer(i4b) :: off(1), left, center, right, n_val
    real(dp)     :: mjd_min, mjd_max, dmjd, threshold, max_slew_length, delta_tod, delta_el, delta_cross
    real(dp), allocatable, dimension(:) :: el_spline, tod_int, del, del2, numsteps, time_tod, tod_resamp
    integer(i4b), allocatable, dimension(:) :: slew_sample_id
    type(spline_type)  :: point_spline
    complex(dpc), dimension(:), allocatable :: fft_el, fft_tod

    mjd_min = max(minval(time_el), minval(time))
    mjd_max = min(maxval(time_el), maxval(time))
    dmjd    = mjd_max - mjd_min
    mjd_min = mjd_min + 0.01d0*dmjd
    mjd_max = mjd_max - 0.01d0*dmjd
    start_tod = locate(time,    mjd_min)
    end_tod   = locate(time,    mjd_max)
    start_el  = locate(time_el, mjd_min)
    end_el    = locate(time_el, mjd_max)
    n_el      = end_el  - start_el  + 1
    n_tod     = end_tod - start_tod + 1
    threshold = 0.3d0
    max_slew_length = 2000 ! Maximum number of samples between CES's

    ! Resample pointing on radiometer time grid
    allocate(el_spline(n_tod))
    call spline(point_spline, time_el(start_el:end_el), el(start_el:end_el))
    do i = start_tod, end_tod
       el_spline(i-start_tod+1) = splint(point_spline, time(i))
    end do

    ! Search for jumps, using the derivative of the elevation
    allocate(del(n_tod), del2(n_tod))
    del(1) = el_spline(2) - el_spline(1)
    do i = 2, n_tod-1
       del(i) = 0.5d0*(el_spline(i+1)-el_spline(i-1))
    end do
    del(n_tod) = el_spline(n_tod)-el_spline(n_tod-1)
    del        = abs(del)

    ! Threshold derivative
    where (del > threshold*maxval(del))
       del = 1.d0
    elsewhere
       del = 0.d0
    end where
    call median_filter(del, del2, 10)
    del = del2

    ! Identify valid slews
    allocate(slew_sample_id(10000))
    j = 0
    i = 1
    do while (i < n_tod)
       if (del(i) == 0.d0) then
          i = i+1
       else
          k = 0
          do while (del(i) == 1.d0)
             k = k+1
             i = i+1
             if (i == n_tod) exit
          end do
          if (k <= max_slew_length) then
             if (j == 0) then
                j                 = j+1
                slew_sample_id(j) = i - k/2
             else if (i-k/2-slew_sample_id(j) > 2*max_slew_length) then
                j                 = j+1
                slew_sample_id(j) = i - k/2
             end if
          end if
       end if
    end do

    ! Allocate shifted time and tod arrays 
    allocate(time_tod(n_tod), tod_resamp(n_tod))
    time_tod = time(start_tod:end_tod)
    tod_resamp = tod(start_tod:end_tod)

    ! Search for crude offset, giving equal weight to all points
    allocate(numsteps(-n_tod:n_tod))
    numsteps = 0.d0
    off          = -10000000
    do i = -n_tod, n_tod
       n_val    = 0
       do k = 1, j
          left   = slew_sample_id(k) + i - 0.5*max_slew_length
          center = slew_sample_id(k) + i 
          right  = slew_sample_id(k) + i + 0.5*max_slew_length
          if (left < 1 .or. right > n_tod)     cycle
          if (left-i < 1 .or. right-i > n_tod) cycle
          n_val = n_val+1
          delta_tod   = (tod_resamp(right)-tod_resamp(center))*(tod_resamp(center)-tod_resamp(left))
          delta_el    = (el_spline(right-i)-el_spline(center-i))*(el_spline(center-i)-el_spline(left-i))
          delta_cross = (tod_resamp(right)-tod_resamp(left))*(el_spline(right-i)-el_spline(left-i))
          if (delta_tod > 0.d0 .and. delta_el > 0.d0 .and. delta_cross < 0.d0) then
             numsteps(i) = numsteps(i) + 1.d0
          end if
       end do
       if (n_val > j/2) then 
          numsteps(i) = numsteps(i) / n_val
       else
          numsteps(i) = 0.d0
       end if
    end do

    where (numsteps == maxval(numsteps))
       numsteps = 1.d0
    elsewhere
       numsteps = 0.d0
    end where
    deallocate(del2)
    allocate(del2(-n_tod:n_tod))
    do i = -n_tod, n_tod
       p = max(i-100,-n_tod)
       q = min(i+100,n_tod)
       del2(i) = maxval(numsteps(p:q))
    end do
    numsteps = del2

    p = -n_tod
    do while (numsteps(p) == 0.d0)
       p = p+1
       if (p == n_tod) exit
    end do
       
    q = p
    do while (numsteps(q) == 1.d0)
       q = q+1
       if (q == n_tod) exit
    end do

    off = (q+p)/2

    ! Return offset
    offset = off(1) * (time(2)-time(1))

!!$    open(58,file='test.dat')
!!$    do i = 1, n_tod
!!$       write(58,*) time_tod(i), tod_resamp(i)
!!$    end do
!!$    write(58,*)
!!$    do i = 1, n_tod
!!$       if (i+off(1)>1 .and. i+off(1) < n_tod) then
!!$          write(58,*) time_tod(i), el_spline(i+off(1))
!!$       end if
!!$    end do
!!$    close(58)
!!$    call mpi_finalize(ierr)
!!$    stop

    deallocate(el_spline, del, del2, numsteps, slew_sample_id, time_tod, tod_resamp)
    call free_spline(point_spline)

  end subroutine compute_time_offset


  subroutine excise_subscan(mjd, data_l1, data_l2)
    implicit none
    real(dp),                                   intent(in)    :: mjd(2)
    type(Lx_struct),                            intent(in)    :: data_l1
    type(Lx_struct),                            intent(inout) :: data_l2

    
    integer(i4b) :: i, j, k, l, m, n, nsamp, nsamp_tot, num_l1_files, nfreq, nsb, ndet, nsamp_point, buffer, ndet0, err
    integer(i4b) :: ind(2), ind_point(2)
    real(dp)     :: samprate, mjd_min, mjd_max
    type(Lx_struct)   :: data_l2_fullres
    type(spline_type), allocatable, dimension(:) :: point_tel_spline, point_cel_spline
    
    ! Find basic information
    nsamp       = size(data_l1%tod,1)
    nfreq       = size(data_l1%tod,2)
    nsb         = size(data_l1%tod,3)
    ndet        = size(data_l1%tod,4)
    mjd_min     = max(mjd(1), minval(data_l1%time))
    mjd_max     = min(mjd(2), maxval(data_l1%time))
    
    ! Find start position
    ind(1) = 1
    do while (ind(1) <= nsamp)
       if (data_l1%time(ind(1)) > mjd_min) exit
       ind(1) = ind(1) + 1
    end do
    
    ! Find end position
    ind(2) = nsamp
    do while (ind(2) >= 1)
       if (data_l1%time(ind(2)) < mjd_max) exit
       ind(2) = ind(2) - 1
    end do
    
    nsamp_tot = ind(2)-ind(1)+1

    ! Allocate full-resolution L2 structure
    allocate(data_l2%time(nsamp_tot), stat=err)
    allocate(data_l2%nu(nfreq,nsb,ndet))
    allocate(data_l2%tod(nsamp_tot, nfreq, nsb, ndet))
    allocate(data_l2%point_tel(3,nsamp_tot,ndet))
    allocate(data_l2%point_cel(3,nsamp_tot,ndet))
    allocate(data_l2%pixels(ndet))
    allocate(data_l2%tod_mean(nfreq, nsb, ndet))
    allocate(data_l2%freqmask(nfreq,nsb,ndet))
    allocate(data_l2%freqmask_full(size(data_l1%nu,1,1),nsb,ndet))
    !allocate(data_l2%flag(nsamp_tot))

    ! Merge L1 data
    data_l2%decimation_time = 1
    data_l2%decimation_nu   = 1
    data_l2%samprate        = data_l1%samprate
    data_l2%nu              = data_l1%nu
    data_l2%pixels          = data_l1%pixels
    data_l2%point_tel       = data_l1%point_tel(:,ind(1):ind(2),:)
    data_l2%point_cel       = data_l1%point_cel(:,ind(1):ind(2),:)
    data_l2%time            = data_l1%time(ind(1):ind(2))
    data_l2%freqmask        = data_l1%freqmask
    data_l2%freqmask_full   = data_l1%freqmask_full

    do j = 1, ndet
       do m = 1, nsb
          do n = 1, nfreq
             data_l2%tod(:,n,m,j)    = data_l1%tod(ind(1):ind(2),n,m,j)
             data_l2%tod_mean(n,m,j) = mean(data_l1%tod(ind(1):ind(2),n,m,j))
          end do
       end do
    end do

  end subroutine excise_subscan



  subroutine decimate_L2_data(samprate_out, numfreq_out, data_in, data_out)
    implicit none
    real(dp),                          intent(in)  :: samprate_out
    integer(i4b),                      intent(in)  :: numfreq_out
    type(Lx_struct),                   intent(in)  :: data_in
    type(Lx_struct),                   intent(out) :: data_out

    integer(i4b) :: i, j, k, l, m, n, nsamp_in, nsamp_out, ndet, dt, dnu, nsb, numfreq_in
    real(dp)     :: w, weight
    real(dp), allocatable, dimension(:,:,:) :: sigmasq

    nsb                      = size(data_in%tod,3)
    ndet                     = size(data_in%tod,4)
    numfreq_in               = size(data_in%tod,2)
    data_out%samprate        = samprate_out
    dt                       = nint(data_in%samprate/samprate_out)
    data_out%decimation_time = dt
    call assert(data_out%decimation_time >= 1, 'Cannot ask for higher output sample rate than input')

    dnu                    = size(data_in%nu,1)/numfreq_out
    data_out%decimation_nu = dnu
    call assert(data_out%decimation_nu >= 1, 'Cannot ask for more frequencies than in input files')

    nsamp_out = int(size(data_in%time)/data_out%decimation_time)
 
    allocate(data_out%time(nsamp_out))
    allocate(data_out%nu(numfreq_out,nsb,ndet))
    allocate(data_out%tod(nsamp_out, numfreq_out, nsb, ndet))
    allocate(data_out%tod_mean(numfreq_in, nsb, ndet))
    allocate(data_out%point_tel(3,nsamp_out,ndet))
    allocate(data_out%point_cel(3,nsamp_out,ndet))
    allocate(data_out%flag(nsamp_out))
    allocate(data_out%freqmask(numfreq_out,nsb,ndet))
    allocate(data_out%freqmask_full(size(data_in%nu,1,1),nsb,ndet))
    allocate(data_out%mean_tp(size(data_in%nu,1),nsb,ndet))
    allocate(data_out%var_fullres(size(data_in%nu,1),nsb,ndet))
    allocate(data_out%n_nan(size(data_in%nu,1),nsb,ndet))
    
    allocate(data_out%acceptrate(nsb,ndet))
    allocate(data_out%diagnostics(size(data_in%nu,1),nsb,ndet,size(data_in%diagnostics,4)))
    allocate(data_out%cut_params(size(data_in%cut_params,1),size(data_in%cut_params,2)))
    allocate(data_out%pixels(ndet))
    if (allocated(data_in%mean_tp)) then
       data_out%mean_tp = data_in%mean_tp
    else
       data_out%mean_tp = 0.d0
    end if
    data_out%freqmask      = data_in%freqmask
    data_out%freqmask_full = data_in%freqmask_full
    data_out%pixels        = data_in%pixels
    data_out%acceptrate    = data_in%acceptrate
    
    data_out%n_pca_comp    = data_in%n_pca_comp
    if (data_in%n_pca_comp > 0) then
       allocate(data_out%pca_ampl(size(data_in%nu,1),nsb,ndet,size(data_in%pca_ampl,4)))
       allocate(data_out%pca_comp(size(data_in%pca_comp,1),size(data_in%pca_comp,2)))
       allocate(data_out%pca_eigv(size(data_in%pca_eigv,1)))

       data_out%pca_ampl      = data_in%pca_ampl
       data_out%pca_comp      = data_in%pca_comp
       data_out%pca_eigv      = data_in%pca_eigv
    end if
    ! Make angles safe for averaging
    do j = 1, ndet
       if (.not. is_alive(j)) cycle
       call make_angles_safe(data_in%point_tel(1,:,j), real(360.d0,sp)) ! Phi
       call make_angles_safe(data_in%point_tel(3,:,j), real(360.d0,sp)) ! Psi
       call make_angles_safe(data_in%point_cel(1,:,j), real(360.d0,sp)) ! Phi
       call make_angles_safe(data_in%point_cel(3,:,j), real(360.d0,sp)) ! Psi
    end do

    ! Compute variance per frequency channel
    !open(58,file='variance.dat')
!    open(58,file='freqmask_2036.dat')
    do k = 1, ndet
       if (nsamp_out == 0) cycle
       if (.not. is_alive(k)) cycle
       do j = 1, nsb
          data_out%tod_mean(:,j,k) = data_in%tod_mean(:,j,k)
          do i = 1, size(data_in%nu,1)
             if (data_out%freqmask_full(i,j,k) == 0) then
                data_out%var_fullres(i,j,k) = 2.8d-5
                cycle
             end if
             data_out%var_fullres(i,j,k) = variance(data_in%tod(:,i,j,k))
!             if (data_out%var_fullres(i,j,k) > 2.8d-5) write(58,*) '   ', i, j, k
             !write(58,*) k, j, i, data_out%var_fullres(i,j,k)
          end do
          !write(58,*)
       end do
    end do
!    close(58)
!    call mpi_finalize(ierr)
!    stop

    !$OMP PARALLEL PRIVATE(i,j,k,l,n,m,weight,w)
    !$OMP DO SCHEDULE(guided)    
    do i = 1, nsamp_out
       data_out%time(i) = mean(data_in%time((i-1)*dt+1:i*dt))  ! Time
       data_out%flag(i) = data_in%time((i-1)*dt+1)             ! Pick first flag in segment
       do j = 1, ndet
          data_out%point_tel(1,i,j) = mean(data_in%point_tel(1,(i-1)*dt+1:i*dt,j)) ! Phi
          data_out%point_tel(2,i,j) = mean(data_in%point_tel(2,(i-1)*dt+1:i*dt,j)) ! Theta
          data_out%point_tel(3,i,j) = mean(data_in%point_tel(3,(i-1)*dt+1:i*dt,j)) ! Psi
          data_out%point_cel(1,i,j) = mean(data_in%point_cel(1,(i-1)*dt+1:i*dt,j)) ! Phi
          data_out%point_cel(2,i,j) = mean(data_in%point_cel(2,(i-1)*dt+1:i*dt,j)) ! Theta
          data_out%point_cel(3,i,j) = mean(data_in%point_cel(3,(i-1)*dt+1:i*dt,j)) ! Psi
       end do

       do j = 1, nsb
          do k = 1, numfreq_out
             do l = 1, ndet           ! Time-ordered data
                if (.not. is_alive(l)) cycle
                data_out%nu(k,j,l)    = mean(data_in%nu((k-1)*dnu+1:k*dnu,j,l)) ! Frequency
                data_out%tod(i,k,j,l) = 0.d0
                weight                = 0.d0
                do n = (k-1)*dnu+1, k*dnu
                   if (data_out%freqmask_full(n,j,l) == 0) cycle
                   if (data_out%var_fullres(n,j,l) <= 0) then
                      w = 0.d0
                   else
                      w      = 1.d0 / data_out%var_fullres(n,j,l) * data_in%freqmask_full(n,j,l)
                   end if
                   weight = weight + w 
                   do m = (i-1)*dt+1, i*dt
                      data_out%tod(i,k,j,l) = data_out%tod(i,k,j,l) + w * data_in%tod(m,n,j,l) 
                   end do

                end do
                if (weight > 0.d0) then
                   data_out%tod(i,k,j,l) = data_out%tod(i,k,j,l) / weight
                else
                   data_out%tod(i,k,j,l) = 0.d0
                end if
             end do
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Polyfiltered TOD
    data_out%polyorder = data_in%polyorder
    if (data_out%polyorder >= 0) then
       allocate(data_out%tod_poly(nsamp_out, 0:data_out%polyorder, nsb, ndet))
       do l = 1, ndet           
          if (.not. is_alive(l)) then
             data_out%tod_poly(:,:,:,l) = 0.d0
             cycle
          end if
          do j = 1, nsb
             do k = 0, data_out%polyorder
                do i = 1, nsamp_out
                   data_out%tod_poly(i,k,j,l) = mean(data_in%tod_poly((i-1)*dt+1:i*dt,k,j,l))
                end do
             end do
          end do
       end do
    end if

    allocate(data_out%sec(nsamp_out))
    data_out%mjd_start = minval(data_out%time)  ! Starting MJD
    data_out%sec       = (data_out%time - data_out%mjd_start) * 24.d0 * 3600.d0

  end subroutine decimate_L2_data

  subroutine simulate_gain_data(rng_handle, data)
    implicit none
    type(planck_rng), intent(inout) :: rng_handle
    type(Lx_struct),  intent(inout) :: data

    real(dp)     :: T_0, Tsys, bw, gain, sigma, rms_drift, max_drift_freq
    integer(i4b) :: i, j, k, l, s, n, ind_cut
    real(sp),     dimension(:), allocatable :: drift
    complex(spc), dimension(:), allocatable :: ffts

    n         = size(data%time)
    T_0       = 8.d0   ! Atmosphere brightness temperature at zenith in K
    Tsys      = 35.d0  ! System temperature in K
    bw        = 8.d9 / 1024.d0   ! Assume 8 GHz bandwidth and 1024 channels
    gain      = 5d8    ! Gain in ADU/K
    sigma     = Tsys / sqrt(bw/data%samprate)
    rms_drift = 5d-6
    max_drift_freq = 1.d0/3.d0 ! High-frequency cutoff in Hz
    ind_cut   = freq2ind(max_drift_freq, data%samprate, n)

    ! Replace TOD with simulated data
    allocate(ffts(0:n/2), drift(1:n))
    do k = 1, size(data%tod,4)        ! Detector  
       do s = 1, size(data%tod,3)     ! Sideband
          do j = 1, size(data%tod,2)     ! Frequency
             ! Start with large-scale drifts
             !data%tod(:,j,k) = 36.d0 ! Monopole
             data%tod(:,j,s,k) = 0.d0 ! Monopole
             ffts = 0.d0
             do l = 1, ind_cut
                ffts(l) = sqrt(rms_drift * (real(l,dp)/real(ind_cut,dp))**(-3.d0)) * &
                     & cmplx(rand_gauss(rng_handle),rand_gauss(rng_handle))
             end do
             call fft(drift, ffts, -1)
             !data%tod(:,j,k) = data%tod(:,j,k) + drift
             
             do i = 1, size(data%tod,1)  ! Time
                data%tod(i,j,s,k) = data%tod(i,j,s,k) + T_0 / sin(data%point_tel(2,i,k)*DEG2RAD)  ! Co-secant model
                !data%tod(i,j,k) = data%tod(i,j,k) + sigma * rand_gauss(rng_handle) ! Radiometer equation noise
                data%tod(i,j,s,k) = data%tod(i,j,s,k) * gain  ! Apply gain
             end do
          end do
       end do
    end do
    deallocate(ffts, drift)

    if (myid == 0) then
       open(58,file='gain_sim.dat')
       do i = 1, size(data%tod,1)
          write(58,*) data%time(i), data%tod(i,1,1,1), data%point_tel(2,i,1)
          !write(58,*) i, data%tod(i,1,1), data%point_tel(2,i)
       end do
       close(58)
    end if

  end subroutine simulate_gain_data

  subroutine convert_GHz_to_k(data)
    implicit none
    type(Lx_struct), intent(inout)  :: data


  end subroutine convert_GHz_to_k

  subroutine correct_missing_time_steps(time)
    implicit none
    real(dp), dimension(:), intent(inout) :: time

    integer(i4b) :: i, n, step
    real(dp)     :: dt

    ! Check time array
    n  = size(time)
    dt = time(2)-time(1)
    if (dt == 0) then
       write(*,*) 'Error: First time sample is buggy'
       call mpi_finalize(ierr)
       stop
    end if

    ! Check each sample
    do i = 2, n
       if (abs(time(i)-(time(i-1)+dt))/dt > 0.01d0) time(i) = time(i-1) + dt
    end do

  end subroutine correct_missing_time_steps

  subroutine initialize_fullres_frequency_mask(freqmaskfile, data)
    implicit none
    character(len=*),                                intent(in)    :: freqmaskfile
    type(Lx_struct),                                 intent(inout) :: data

    integer(i4b) :: i, j, k, nfreq_full, nsb, ndet, unit, det, sb, freq, dfreq, ierr
    logical(lgt) :: first
    character(len=1024) :: line, val, equal

    ndet = get_num_dets()
    nsb  = get_num_sideband()
    unit = getlun()
    open(unit, file=trim(freqmaskfile), recl=1024)
    first = .true.
    do while (.true.)
       read(unit,'(a)', end=99) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (first) then
          read(line,*) nfreq_full
          allocate(data%freqmask_full(nfreq_full,nsb,ndet))
          data%freqmask_full = 1.d0
          first = .false.
       else
          read(line,*) freq, sb, det
          if (det == 0 .and. sb == 0 .and. freq == 0) then
             write(*,*) 'ERROR: All frequencies removed by freqmask!'
             stop
          else if (det == 0 .and. sb == 0) then
             data%freqmask_full(freq,:,:) = 0.d0
          else if (det == 0 .and. freq == 0) then
             data%freqmask_full(:,sb,:) = 0.d0
          else if (sb == 0 .and. freq == 0) then
             data%freqmask_full(:,:,det) = 0.d0
          else if (det == 0) then
             data%freqmask_full(freq,sb,:) = 0.d0
          else if (sb == 0) then
             data%freqmask_full(freq,:,det) = 0.d0
          else if (freq == 0) then
             data%freqmask_full(:,sb,det) = 0.d0
          else 
             data%freqmask_full(freq,sb,det) = 0.d0
          end if
          if (all(data%freqmask_full == 0)) then
             write(*,*) 'ERROR: All frequencies removed by freqmask!'
             stop
          end if
       end if
    end do

99  close(unit)

    do k = 1, ndet
       if (.not. is_alive(k)) data%freqmask_full(:,:,k) = 0.d0
    end do

  end subroutine initialize_fullres_frequency_mask

  subroutine postprocess_frequency_mask(nfreq, data, sid)
    implicit none
    integer(i4b),                                    intent(in)    :: sid, nfreq
    type(Lx_struct),                                 intent(inout) :: data
    

    integer(i4b) :: i, j, k, nfreq_full, nsb, ndet, unit, det, sb, freq, dfreq, ierr
    logical(lgt) :: first
    character(len=1024) :: line, val, equal

    ndet       = get_num_dets()
    nsb        = get_num_sideband()
    nfreq_full = size(data%freqmask_full,1)

    ! Exclude frequencies with any NaNs
    do k = 1, ndet
       do j = 1, nsb
          do i = 1, nfreq_full
             if (data%freqmask_full(i,j,k) == 0.d0) cycle
             if (any(data%tod(:,i,j,k) .ne. data%tod(:,i,j,k))) then
                write(*,fmt='(a,i8,3i6)') '   Rejecting NaNs, (sid,det,sb,freq) = ', sid, k,j,i
                data%freqmask_full(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
    
    ! Downgrade mask
    allocate(data%freqmask(nfreq,nsb,ndet))
    dfreq = nfreq_full/nfreq
    if (dfreq == 1) then
       data%freqmask = data%freqmask_full
    else
       data%freqmask = 1.
       do k = 1, ndet
          do j = 1, nsb
             do i = 1, nfreq
                if (all(data%freqmask_full((i-1)*dfreq+1:i*dfreq,j,k) == 0.d0)) data%freqmask(i,j,k) = 0.0
             end do
          end do
       end do
    end if

  end subroutine postprocess_frequency_mask

  subroutine remove_elevation_gain(data) 
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: n, ndet, nsb, nfreq, i, j, k, l, m, tmin, tmax, parfile_time
    real(dp)        :: g, a, sigma0, chisq, tsys, tau, dnu, const
    real(dp), dimension(:), allocatable :: el, dat

!    m=a*data%samprate/data%scanfreq(2)      ! Number of tod-samples per gain estimate
    m   = (size(data%time))
    write(*,*) 'Number of samples per gain estimate    =', m
!    n   = (size(data%time)+m-1)/m           ! Number of gain samples
    n   = (size(data%time))/m               ! Number of gain samples
    
    write(*,*) 'Number of gain estimates for this scan =', n
    write(*,*) 'n*m                                    =', n*m
    write(*,*) 'Total number of samples                =', size(data%time)
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)
    write(*,*) nfreq, '= nfreq', nsb, '= nsb', ndet, '= ndet'
    write(*,*) '---------------------------------------------------------'
    !allocate(data%time_gain(n), data%gain(n,nfreq,nsb,ndet))
    allocate(el(m), dat(m))
    !data%time_gain = data%time(::m) ! OBS may reallocate length of time_gain!!!
    !open(13,file='gain.dat')
    !open(14,file='chisq.dat')
    do k = 1, ndet
       do l = 1, nsb
          do j = 1, nfreq
          if (data%freqmask_full(j,l,k) == 0.d0) cycle
          do i = 1, n
             tmin = (i-1)*m+1
             tmax = i*m
             el  = data%point_tel(2,tmin:tmax,k)
             if (any(el == 0.d0)) write(*,*) k, l, j, i
             dat = data%tod(tmin:tmax,j,l,k)
             !write(*,*) tmin, tmax, 'min max'
             call estimate_gain(el,dat,g)
             data%tod(tmin:tmax,j,l,k) = data%tod(tmin:tmax,j,l,k) - g*1/(sin(el*pi/180.))
!             write(13,*) g
!             write(13,*) (g-5.6d10)/5.6d10*100
!             write(14,*) chisq
          end do
          end do
       end do
    end do
    !close(13)
    !close(14)
    deallocate(el, dat)
  end subroutine remove_elevation_gain


  subroutine compute_Tsys(tsys_file, data)
    implicit none
    character(len=*),            intent(in)       :: tsys_file
    type(Lx_struct),             intent(inout)    :: data
    type(hdf_file)                                :: file
    real(dp), dimension(:,:,:,:), allocatable     :: tsys_fullres
    real(dp), dimension(:,:,:), allocatable       :: tsys_1, tsys_2, mean_tod
    integer(i4b)                                  :: nfreq, nfreq_fullres, nsb, ndet, i, j, k, l, mjd_index1,mjd_index2, nsamp, dnu
    integer(i4b)                                  :: nsamp_gain(7), num_bin, n, mean_count
    integer(i4b), dimension(:), allocatable       :: scanID
    real(dp)                                      :: mjd_high,w, sum_w_t, sum_w, t1, t2, tsys, tod_mean_dec
    real(dp), dimension(:), allocatable           :: time

    allocate(data%Tsys(2, nfreq_fullres, nsb, ndet))                                                                            
    nsamp = size(data%tod,1)
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)

    ! 1) get tsys-values                                                                                                                        
    call open_hdf_file(tsys_file, file, "r")
    call read_hdf(file, "nfreq", nfreq_fullres)
    call get_size_hdf(file, "MJD", nsamp_gain)

    allocate(tsys_fullres(ndet, nsb, nfreq_fullres, nsamp_gain(1)))
    allocate(tsys_1(nfreq_fullres, nsb, ndet))
    allocate(tsys_2(nfreq_fullres, nsb, ndet))
    allocate(mean_tod(nfreq_fullres, nsb, ndet))
    allocate(time(nsamp_gain(1)))
    call read_hdf(file, "MJD", time)
    call read_hdf(file, "tsys_pr_P", tsys_fullres)
    call close_hdf_file(file)

    ! finding closest time-value                                                                                                                
    mjd_index1 = max(locate(time, data%time(1)),1)

    if (mjd_index1 < data%time(1)) then
        mjd_index2 = mjd_index1 + 1
    else
       mjd_index2   = mjd_index1 - 1
    end if


    tsys_1 = tsys_fullres(mjd_index1,:,:,:)
    tsys_2 = tsys_fullres(mjd_index1,:,:,:)
    
    do i=1, ndet
       write(*,*) "det ", i
       if (.not. is_alive(i)) cycle
       do j=1, nsb
          do k=1, nfreq
             mean_tod = 0.d0
             mean_count = 0.d0
             do l=1, nsamp_gain(1)
                if ( .not. isnan(data%tod(l,k,j,i))) then
                   ! Calculating the nanmean in time of each tod
                   mean_tod = mean_tod + data%tod(l,k,j,i)
                   mean_count = mean_count + 1
                end if
             end do
             !mean(x[!isnan(x)])
             mean_tod = mean_tod / mean_count
             data%Tsys(1,k,j,i) = tsys_1(k,j,i)*mean_tod(k,j,i)*data%freqmask_full(k,j,i)
             data%Tsys(2,k,j,i) = tsys_2(k,j,i)*mean_tod(k,j,i)*data%freqmask_full(k,j,i)
          end do
       end do
    end do
  deallocate(tsys_fullres)
  deallocate(tsys_1)
  deallocate(tsys_2)
  deallocate(mean_tod)
  deallocate(time)

  end subroutine compute_Tsys

  subroutine calibrate_tod(data_l1, data_l2_fullres)
    implicit none
    type(lx_struct), intent(in)    :: data_l1
    type(lx_struct), intent(inout) :: data_l2_fullres

    integer(i4b) :: i, j, k, ndet, nfreq, nsb

    nfreq = size(data_l1%tod,2)
    nsb   = size(data_l1%tod,3)
    ndet  = size(data_l1%tod,4)

    ! Interpolate Tsys to current time for each detector
    do i = 1, ndet
       do j = 1, nsb
          do k = 1, nfreq
             data_l2_fullres%tod(:,k,j,i)  = (data_l1%Tsys(1,k,j,i) + data_l1%Tsys(2,k,j,i))/2.d0 * data_l1%tod(:,k,j,i)
          end do
       end do
    end do
  end subroutine calibrate_tod

  subroutine fit_noise(data_l2)
    implicit none

    type(lx_struct), intent(inout) :: data_l2

    integer(i4b) :: i, j, k, ndet, nfreq, nsb

    nfreq = size(data_l2%tod,2)
    nsb   = size(data_l2%tod,3)
    ndet  = size(data_l2%tod,4)    

    ! Approximate noise as white for now
    allocate(data_l2%sigma0(nfreq,nsb,ndet))
    allocate(data_l2%alpha(nfreq,nsb,ndet))
    allocate(data_l2%fknee(nfreq,nsb,ndet))
    do i = 1, ndet
       do j = 1, nsb
          do k = 1, nfreq
             data_l2%sigma0(k,j,i) = sqrt(variance(data_l2%tod(:,k,j,i)))
             data_l2%alpha(k,j,i)  = -2.d0
             data_l2%fknee(k,j,i)  = 1.d-6   ! Set to a very low value
          end do
       end do
    end do

  end subroutine fit_noise


end program
