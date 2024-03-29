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
   use spline_2D_mod
   use quiet_status_mod 
   use comap_gain_mod
   use comap_patch_mod
   use comap_ephem_mod
   use comap_sim2tod_mod
   implicit none

   character(len=512)   :: parfile, runlist, l1dir, l2dir, tmpfile, freqmaskfile, monitor_file_name, tsysfile, freq_import_dir, freq_dir, sigma_import_dir, sigma_dir
   character(len=9)     :: id_old
   character(len = 1024)  :: freq_import_name, sigma_import_name      
   character(len=512)   :: param_dir, runlist_in, simulation_path
   character(len=10)    :: target_name
   character(len=1024)  :: param_name, param_name_raw, runlist_name, runlist_name_raw, cal_db
   integer(i4b)         :: i, j, k, l, m, n, snum, nscan, unit, myid, nproc, ierr, ndet, npercore, irun
   integer(i4b)         :: mstep, i2, decimation, nsamp, numfreq, n_nb, mask_outliers, n_tsys, polyorder_store, n_pca_store, n_pca_store_feed, numfreq_pca_downsamp
   integer(i4b)         :: debug, num_l1_files, seed, bp_filter, bp_filter0, n_pca_comp, n_pca_comp_feed, pca_max_iter, tsys_ind(2), cal_method
   real(dp)             :: todsize, nb_factor, min_acceptrate, pca_sig_rem, var_max, corr_max, tsys_mjd_max, tsys_mjd_min, tsys_time(2)
   real(dp)             :: pca_err_tol, corr_cut, mean_corr_cut, mean_abs_corr_cut, med_cut, var_cut, sim_tsys, max_tsys
   logical(lgt)         :: exist, reprocess, check_existing, gonext, found, rm_outliers
   logical(lgt)         :: import_freqmask, import_sigma
   type(comap_scan_info) :: scan
   type(Lx_struct)      :: data_l1, data_l2_fullres, data_l2_decimated, data_l2_filter, data_l2_import
   type(planck_rng)     :: rng_handle 
   type(status_file)    :: status
   type(patch_info)     :: pinfo
   !real(dp),     allocatable, dimension(:,:,:,:)   :: store_l2_tod
   logical(lgt)         :: process, is_sim, rem_el, verb, diag_l2, use_freq_filter, is_sim2tod
   type(simulation_struct) :: simdata
   real(dp)             :: timing_offset, mjd(2), dt_error, samprate_in, samprate, scanfreq, nu_gain, alpha_gain, t1, t2, boost_factor
   
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
   call get_parameter(unit, parfile, 'N_FEED_PCA_COMPONENTS',     par_int=n_pca_comp_feed)
   call get_parameter(unit, parfile, 'PCA_ERROR_TOLERANCE',       par_dp=pca_err_tol)
   call get_parameter(unit, parfile, 'PCA_MAX_ITERATIONS',        par_int=pca_max_iter)
   call get_parameter(unit, parfile, 'MASK_OUTLIERS',             par_int=mask_outliers)
   call get_parameter(unit, parfile, 'MIN_ACCEPTRATE',            par_dp=min_acceptrate)
   call get_parameter(unit, parfile, 'LEVEL2_DIR',                par_string=l2dir)
   call get_parameter(unit, parfile, 'PCA_NSIGMA_REMOVE',         par_dp=pca_sig_rem)
   call get_parameter(unit, parfile, 'IS_SIM',                    par_lgt=is_sim)
   call get_parameter(unit, parfile, 'SIM_TSYS',                  par_dp=sim_tsys)
   call get_parameter(unit, parfile, 'MAX_TSYS',                  par_dp=max_tsys)
   call get_parameter(unit, parfile, 'CALIBRATION_METHOD',        par_int=cal_method)
   call get_parameter(unit, parfile, 'CAL_DATABASE_FILE',         par_string=cal_db)
   call get_parameter(unit, parfile, 'REMOVE_ELEVATION_TEMP',     par_lgt=rem_el)
   call get_parameter(unit, parfile, 'VERBOSE_PRINT',             par_lgt=verb)
   call get_parameter(unit, parfile, 'RETURN_DIAG_L2_FILES',      par_lgt=diag_l2)
   call get_parameter(unit, parfile, 'IMPORT_FREQMASK',           par_lgt=import_freqmask)
   call get_parameter(unit, parfile, 'IMPORT_SIGMA0',             par_lgt=import_sigma)
   call get_parameter(unit, parfile, 'TARGET_NAME',               par_string=target_name)
   call get_parameter(unit, parfile, 'RUNLIST',                   par_string=runlist_in)
   call get_parameter(unit, parfile, 'USE_FREQ_FILTER',           par_lgt=use_freq_filter)
   call get_parameter(unit, parfile, 'NUMFREQ_PCA_DOWNSAMPLER',   par_int=numfreq_pca_downsamp)
   call get_parameter(unit, parfile, 'RUN_SIM2TOD',               par_lgt = is_sim2tod)
   call get_parameter(unit, parfile, 'DATACUBE',                  par_string = simulation_path)
   call get_parameter(unit, parfile, 'BOOST_FACTOR',              par_dp = boost_factor)
   
   if (import_freqmask) then 
      call get_parameter(unit, parfile, 'FREQMASK_L2_FOLDER',     par_string=freq_import_dir)
   else if (import_sigma) then
      call get_parameter(unit, parfile, 'SIGMA0_L2_FOLDER',       par_string=sigma_import_dir)
   end if 
 
   ! Require that outliers are masked if mask is to be imported
   if (import_freqmask .and. .not. mask_outliers) then 
      mask_outliers = .true. 
   end if 
   
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
   
   if (myid == 0) then 
      irun = 1
      ! Copy parameter file and runlists to l2-file output directory 
      param_dir = trim(l2dir)//"/"//trim(target_name)//"/param4level2"
      inquire(directory=trim(param_dir), exist=exist)
      if (.not. exist) then 
         call execute_command_line("mkdir "//trim(param_dir), wait=.true.)
      end if
 
      param_name_raw = trim(param_dir)//"/param_"
      runlist_name_raw = trim(param_dir)//"/runlist_"
      write(param_name, "(A512, I6.6, A4)") trim(param_name_raw), 0, ".txt"
      write(runlist_name, "(A512, I6.6, A4)") trim(runlist_name_raw), 0, ".txt"
      
      inquire(file=trim(param_name), exist=exist)
      
      if (.not. exist) then
         call execute_command_line("cp "//trim(parfile)//" "//param_name, wait=.true.)
         call execute_command_line("cp "//trim(runlist_in)//" "//runlist_name, wait=.true.)
         irun = irun + 1
      else     
         irun = 1
         exist = .true.
         do while (exist)
            write(param_name, "(A512, I6.6, A4)") trim(param_name_raw), irun, ".txt"
            write(runlist_name, "(A512, I6.6, A4)") trim(runlist_name_raw), irun, ".txt"
            inquire(file=trim(param_name), exist=exist)
            irun = irun + 1
         end do
         call execute_command_line("cp "//trim(parfile)//" "//param_name, wait=.true.)
         call execute_command_line("cp "//trim(runlist_in)//" "//runlist_name, wait=.true.)
      end if
      print *, "Run number: ", irun - 1
   end if
   
   call mpi_bcast(irun,  1, mpi_integer, 0, mpi_comm_world, ierr)
   
   if (is_sim2tod) then
      write(*, *) "Reading signal simulation realization."
      call read_sim_file(simulation_path, simdata)
      simdata%boost = boost_factor
    end if

   nscan    = get_num_scans()
   do snum = myid+1, nscan, nproc     
      call get_scan_info(snum, scan)
      process = .false.
      do k = 2, scan%nsub-1
         inquire(file=scan%ss(k)%l2file,exist=exist)
         if (reprocess .and. exist) call rm(scan%ss(k)%l2file)           
         if (reprocess .or. .not. exist) process = .true.
      end do
      if (.not. process) then
         write(*,fmt="(i3,a,i8,i8)") myid, " skipping already finished obsid:", snum, scan%id
         cycle
      end if
 
      write(*,fmt="(i3,a,i10,a)") myid, " processing obsid ", scan%id, " (" // trim(itoa(snum)) // "/" // trim(itoa(nscan)) // ")"
      call update_status(status, 'scan_start')
 
 
      ! Read in Level 1 file
      call wall_time(t1)
      call read_l1_file(scan%l1file, data_l1, scan%id, verb, init=.false.)
      call update_status(status, 'read_l1')
      
      ! Initialize frequency mask
      !if (.not. import_freqmask) then
      call initialize_fullres_frequency_mask(freqmaskfile, data_l1, verb)
      call update_status(status, 'freq_mask1')
      !end if
      
      if (size(data_l1%tod,1) <100) then
         write(*,*) 'Too few samples in ', scan%id
         cycle
      end if
      !call correct_missing_time_steps(data_l1%time)
      call wall_time(t2)
      todsize = real(size(data_l1%tod,1),dp)*real(size(data_l1%tod(1,:,:,:)),dp)*4.d0/1024.d0**2 ! in MB
      if (verb) then
         write(*,fmt='(i4,a,f10.2,a)') myid, ' -- disk read speed = ', real(todsize,dp)/(t2-t1), ' MB/sec'
      end if
      ! Interpolate over single NaNs
      call interpolate_nans(data_l1, scan%id, verb)
      call update_status(status, 'nan_interp')
 
      ! Finding Phot, Thot and time_hot
      call prepare_calibration(scan, data_l1, is_sim, cal_method, cal_db, verb)
 
 
      ! Finalize frequency mask
      !if (.not. import_freqmask) then
      call postprocess_frequency_mask(numfreq, data_l1, scan%id, verb)
      call update_status(status, 'freq_mask2')
      
 !     n_tsys = 0
 !     nsamp = size(data_l1%tod,1)
 
      !write(*,*) 4
 
      !write(*,*) data_l1%tod(1751:2650,1,1,1)
      ! do i=1, scan%nsub
      !    if (scan%ss(i)%scanmode == 'tsys') then
      !       n_tsys = n_tsys + 1
      !       if (n_tsys > 2) exit
      !       tsys_mjd_min     = max(scan%ss(i)%mjd(1), minval(data_l1%time))
      !       tsys_mjd_max     = min(scan%ss(i)%mjd(2), maxval(data_l1%time))
      !       !write(*,*) "A", tsys_mjd_min, tsys_mjd_max
      !       tsys_ind(1) = 1
      !       do while (tsys_ind(1) <= nsamp)
      !          if (data_l1%time(tsys_ind(1)) > tsys_mjd_min) exit
      !          tsys_ind(1) = tsys_ind(1) + 1
      !       end do
 
      !       ! Find end position
      !       tsys_ind(2) = nsamp
      !       do while (tsys_ind(2) >= 1)
      !          if (data_l1%time(tsys_ind(2)) < tsys_mjd_max) exit
      !          tsys_ind(2) = tsys_ind(2) - 1
      !       end do
      !       !write(*,*) "B", tsys_ind
      !       !stop
 
      !       ! Compute absolute calibration
 
      !       call compute_P_hot(tsysfile, data_l1, tsys_ind, n_tsys, is_sim, verb)
      !       tsys_time(n_tsys) = 0.5*(tsys_mjd_min + tsys_mjd_max)
      !    end if
      ! end do
      ! if (n_tsys == 1) tsys_time(2) = data_l1%time(nsamp)
      ! if (n_tsys == 0) then
      !    write(*,*) "No ambient subscans!"
      !    cycle
      ! end if
 
      ! Get patch info
      found = get_patch_info(scan%object, pinfo) 
      if (.not. found) then
         write(*,*) 'Error: Patch not found in patchfile = ', trim(scan%object)
         call mpi_finalize(ierr)
         stop
      end if
      
      if (is_sim2tod) then
         call prepare_sim_interpolation(simdata, pinfo)
      end if 

      do k = 2, scan%nsub-1
         ! Reformat L1 data into L2 format, and truncate
         call excise_subscan(scan%ss(k)%mjd, data_l1, data_l2_fullres)
         call update_status(status, 'excise')
 
         call find_tsys(data_l2_fullres, is_sim, sim_tsys, scan%ss(k), verb, myid)
         call update_status(status, 'find tsys')
         
         if (is_sim2tod) then
             write(*, *)  "Rank:", myid, "Adding simulations to TODs. "  
             call sim2tod(data_l2_fullres, simdata)
             call update_status(status, 'Add simulated signal')
         end if 

         if (verb) then
            write(*,*) "Rank:", myid, "Starting analysis of scan", scan%ss(k)%id
            write(*,'(A, F18.7, F9.4)') " Time and duration (in mins) of scan: ", data_l2_fullres%time(1), (data_l2_fullres%time(size(data_l2_fullres%time, 1)) - data_l2_fullres%time(1)) * 24 * 60
            write(*,'(A, I8, 2A)') " Number of samples, observed patch: ", size(data_l2_fullres%time, 1), "  ", trim(scan%object)
            write(*,*) "Scan type: ", trim(scan%ss(k)%scanmode)
         end if
 
         
         ! Write diagnostic l2_file to disc
         if (diag_l2) then
            polyorder_store = data_l2_fullres%polyorder
            n_pca_store     = data_l2_fullres%n_pca_comp
            n_pca_store_feed     = data_l2_fullres%n_pca_comp_feed
            data_l2_fullres%polyorder = -1
            data_l2_fullres%n_pca_comp = 0
            data_l2_fullres%n_pca_comp_feed = 0
            data_l2_fullres%mask_outliers = 0
         
            call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
            
            ! Fit noise
            call fit_noise(data_l2_decimated)
            call mkdirs(trim(scan%ss(k)%l2file), .true.)
            write(*,*) "Writing out l2-data for diagnostcs ", adjustl(trim("_1_before_norm"))
            call write_l2_file(scan, k, data_l2_decimated, adjustl(trim("_1_before_norm")))
            data_l2_fullres%polyorder = polyorder_store
            data_l2_fullres%n_pca_comp = n_pca_store
            data_l2_fullres%n_pca_comp_feed = n_pca_store_feed
         end if
         
         ! Normalize gain
         if (trim(pinfo%type) == 'gal' .or. trim(pinfo%type) == 'cosmo') then
            call normalize_gain(data_l2_fullres, nu_gain, alpha_gain, scan%id)
            call update_status(status, 'gain_norm')
         end if
         
         ! Write diagnostic l2_file to disc
         if (diag_l2) then
            polyorder_store = data_l2_fullres%polyorder
            n_pca_store     = data_l2_fullres%n_pca_comp
            n_pca_store_feed     = data_l2_fullres%n_pca_comp_feed
            data_l2_fullres%polyorder = -1
            data_l2_fullres%n_pca_comp = 0
            data_l2_fullres%n_pca_comp_feed = 0
 
            call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
            
            ! Fit noise
            call fit_noise(data_l2_decimated)
            write(*,*) "Writing out l2-data for diagnostcs ", adjustl(trim("_2_after_norm"))
            call write_l2_file(scan, k, data_l2_decimated, adjustl(trim("_2_after_norm")))
            data_l2_fullres%polyorder = polyorder_store
            data_l2_fullres%n_pca_comp = n_pca_store
            data_l2_fullres%n_pca_comp_feed = n_pca_store_feed
         end if
         
         if (rem_el) then
            ! Remove elevation gain
            ! if ((scan%ss(k)%scanmode == 'circ') .or. (scan%ss(k)%scanmode == 'raster') &
            !      .or. (scan%ss(k)%scanmode == 'liss')) then 
            !    if (verb) then
            !       write(*,*) 'Removing elevation gain, id: ', scan%ss(k)%id
            !    end if
            !    call remove_elevation_gain(data_l2_fullres)
            !    call update_status(status, 'remove_elevation')
            ! end if
            if (verb) then
               write(*,*) 'Removing pointing templates, id: ', scan%ss(k)%id
            end if
            call subtract_pointing_templates(data_l2_fullres, scan%ss(k), verb)
            call update_status(status, 'subtract_pointing_templates')
         end if
         
         ! Write diagnostic l2_file to disc
         if (diag_l2) then
            polyorder_store = data_l2_fullres%polyorder
            n_pca_store     = data_l2_fullres%n_pca_comp
            n_pca_store_feed     = data_l2_fullres%n_pca_comp_feed
            data_l2_fullres%polyorder = -1
            data_l2_fullres%n_pca_comp = 0
            data_l2_fullres%n_pca_comp_feed = 0
 
            call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
            
            ! Fit noise
            call fit_noise(data_l2_decimated)
            write(*,*) "Writing out l2-data for diagnostcs ", adjustl(trim("_3_after_template_subtr"))
            call write_l2_file(scan, k, data_l2_decimated, adjustl(trim("_3_after_template_subtr")))
            data_l2_fullres%polyorder = polyorder_store
            data_l2_fullres%n_pca_comp = n_pca_store
            data_l2_fullres%n_pca_comp_feed = n_pca_store_feed
         end if
 
         if (verb) then           
            if (.not. all(data_l2_fullres%tod == data_l2_fullres%tod)) then
               write(*,*) "NaN in tod before filtering!", scan%ss(k)%id
            end if
         end if
 
 
         data_l2_fullres%import_sigma = import_sigma
         data_l2_decimated%import_sigma = import_sigma
 
         if (import_sigma .and. .not. import_freqmask) then
               write (*,*) "Importing sigma0 and var_fullres from existing l2 file"
               
               sigma_dir = scan%ss(k)%l2file
               sigma_import_name = sigma_dir(len_trim(sigma_dir)-20:)
               sigma_import_name = trim(sigma_import_dir)//trim(sigma_import_name)
               
               call read_l2_file(sigma_import_name, data_l2_import)
 
               call transfer_imported_sigma(data_l2_import, data_l2_fullres)
 
               call update_status(status, 'imported_sigma0')
               
               call free_lx_struct(data_l2_import)
         end if 
         
         if (sum(data_l2_fullres%freqmask_full) == 0.d0) then
            if (verb) then
               write(*,*) "All channels masked before masking!!!!!!! Scan: ", scan%ss(k)%id
            end if
            data_l2_fullres%polyorder = -1
            data_l2_fullres%n_pca_comp = 0
            data_l2_fullres%n_pca_comp_feed = 0
            data_l2_fullres%use_freq_filter = .false.
            data_l2_fullres%mask_outliers = 0
            
 
            ! If necessary, decimate L2 file in both time and frequency
            
            call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
            call update_status(status, 'decimate')
            
            ! Fit noise
            if (.not. import_freqmask .and. .not. import_sigma) then 
               call fit_noise(data_l2_decimated)
            else
               call transfer_imported_sigma(data_l2_fullres, data_l2_decimated)
            end if
            
            ! Write L2 file to disk
            if (verb) then
               write(*,*) 'Writing empty file ', scan%ss(k)%id, ' to disk', trim(scan%ss(k)%l2file)
            end if
            
            call mkdirs(trim(scan%ss(k)%l2file), .true.)
            !call write_l2_file(scan%ss(k)%l2file, data_l2_decimated)
            call write_l2_file(scan, k, data_l2_decimated)
            call update_status(status, 'write_l2')
            
            cycle
 
            
         end if
 
         data_l2_fullres%use_freq_filter = use_freq_filter !.true.
         data_l2_fullres%mask_outliers = mask_outliers
         data_l2_fullres%import_freqmask = import_freqmask
         data_l2_decimated%import_freqmask = import_freqmask
         if ((mask_outliers == 1) .and. (.not. (sum(data_l2_fullres%freqmask_full) == 0.d0))) then
            if (import_freqmask) then
               data_l2_import%n_cal = data_l2_fullres%n_cal
 
               freq_dir = scan%ss(k)%l2file
               freq_import_name = freq_dir(len_trim(freq_dir)-20:)
               freq_import_name = trim(freq_import_dir)//trim(freq_import_name)
               
               call read_l2_file(freq_import_name, data_l2_import)
               call transfer_diagnostics(data_l2_import, data_l2_fullres)
               call update_status(status, 'imported_freqmask')
 
               call free_lx_struct(data_l2_import)     
            else
               if (verb) then
                  write(*,*) 'Making frequency mask', scan%ss(k)%id
               end if
               ! Copy tod, run filtering, make new mask, then do filtering again on original (unfiltered) data
               call copy_lx_struct(data_l2_fullres, data_l2_filter)
               call update_status(status, 'copy_data')
 
               ! Poly-filter copied data
               if (data_l2_fullres%use_freq_filter) then
                  !write(*,*) "start"
                  call frequency_filter_TOD(data_l2_filter, nu_gain, alpha_gain, scan%ss(k))
                  !write(*,*) "stop"
               else
                  call polyfilter_TOD(data_l2_filter, bp_filter0)
               end if
               call update_status(status, 'polyfilter0')
               !           write(*,*) sum(data_l2_filter%freqmask_full) / 19.d0 / 4.d0 / 1024.d0 
 
               call find_spikes(data_l2_filter, verb, scan%ss(k)%id)
               call update_status(status, 'find_spikes')
               
               
               ! pca filter copied data
               call pca_filter_TOD(data_l2_filter, n_pca_comp, pca_max_iter, pca_err_tol, pca_sig_rem, verb)
               call update_status(status, 'pca_filter0')
               
               ! feed pca filter copied data
               data_l2_filter%n_freq_downsamp = numfreq_pca_downsamp
 
               call pca_filter_feed_TOD(data_l2_filter, n_pca_comp_feed, pca_max_iter, pca_err_tol, pca_sig_rem, verb)
               call update_status(status, 'pca_feed_filter0')
             
               if (diag_l2) then
                  data_l2_filter%mask_outliers = 0
                  call decimate_L2_data(samprate, numfreq, data_l2_filter, data_l2_decimated)
                  ! write(*,*) "done with decimation"
                  ! Fit noise
                  call fit_noise(data_l2_decimated)
                  write(*,*) "Writing out l2-data for diagnostcs ", adjustl(trim("_4_before_mask"))
                  call write_l2_file(scan, k, data_l2_decimated, adjustl(trim("_4_before_mask")))
                  data_l2_filter%mask_outliers = mask_outliers
               end if
 
               ! flag correlations and variance
               call flag_correlations(data_l2_filter, scan%ss(k)%id, parfile, myid)!corr_cut, mean_corr_cut, mean_abs_corr_cut, med_cut, var_cut, n_nb, nb_factor, var_max, corr_max)
               call update_status(status, 'flag_corr')
 
               ! replace freqmask in original tod
               call transfer_diagnostics(data_l2_filter, data_l2_fullres)
               
 
               call update_freqmask(data_l2_fullres, min_acceptrate, scan%ss(k)%id, .false.)
               call update_status(status, 'made_freqmask')
                       
               call free_lx_struct(data_l2_filter)
 
            end if 
         end if
 
 !        write(*,*) data_l2_fullres%freqmask_reason(:, 1, 17)
  !       write(*,*) data_l2_fullres%freqmask_reason(:, 4, 10)
 
         !!! All channels masked!!
         if (sum(data_l2_fullres%freqmask_full) == 0.d0) then
            
            if (verb) then
               write(*,*) "All channels masked! Scan: ", scan%ss(k)%id
            end if
            data_l2_fullres%polyorder = -1
            data_l2_fullres%n_pca_comp = 0
            data_l2_fullres%n_pca_comp_feed = 0
            data_l2_fullres%use_freq_filter = .false.
 
            ! If necessary, decimate L2 file in both time and frequency
            
            call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
            call update_status(status, 'decimate')
            
            ! Fit noise
            if (.not. import_freqmask .and. .not. import_sigma) then 
               call fit_noise(data_l2_decimated)
            else
               call transfer_imported_sigma(data_l2_fullres, data_l2_decimated)
            end if
            
            ! Write L2 file to disk
            if (verb) then
               write(*,*) 'Writing empty file ', scan%ss(k)%id, ' to disk', trim(scan%ss(k)%l2file)
            end if
            
            call mkdirs(trim(scan%ss(k)%l2file), .true.)
            !call write_l2_file(scan%ss(k)%l2file, data_l2_decimated)
            call write_l2_file(scan, k, data_l2_decimated)
            call update_status(status, 'write_l2')
            
            cycle
 
         end if
         
         ! Poly-filter if requested
         bp_filter = -1; if (trim(pinfo%type) == 'cosmo') bp_filter = bp_filter0
         
         if (data_l2_fullres%use_freq_filter) then
            call frequency_filter_TOD(data_l2_fullres, nu_gain, alpha_gain, scan%ss(k))
         else
            call polyfilter_TOD(data_l2_fullres, bp_filter)
         end if
            
         call update_status(status, 'polyfilter')
         
         ! Write diagnostic l2_file to disc
         if (diag_l2) then
            n_pca_store     = data_l2_fullres%n_pca_comp
            n_pca_store_feed     = data_l2_fullres%n_pca_comp_feed
            data_l2_fullres%n_pca_comp = 0
            data_l2_fullres%n_pca_comp_feed = 0
 
            call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
            
            ! Fit noise
            call fit_noise(data_l2_decimated)
            write(*,*) "Writing out l2-data for diagnostcs ", adjustl(trim("_5_after_poly"))
            call write_l2_file(scan, k, data_l2_decimated, adjustl(trim("_5_after_poly")))
            data_l2_fullres%n_pca_comp = n_pca_store
            data_l2_fullres%n_pca_comp_feed = n_pca_store_feed
         end if
         
         if ((mask_outliers == 0) .and. (bp_filter > -1)) then
            call find_spikes(data_l2_fullres, verb, scan%ss(k)%id)
            call update_status(status, 'find_spikes')
         end if
         !call freq_filter_TOD(data_l2_fullres)
 
         ! pca filter after polyfilter
         if (trim(pinfo%type) == 'cosmo') then
            call pca_filter_TOD(data_l2_fullres, n_pca_comp, pca_max_iter, pca_err_tol, pca_sig_rem, verb)
            call update_status(status, 'pca_filter')
            
            data_l2_fullres%n_freq_downsamp = numfreq_pca_downsamp
 
            call pca_filter_feed_TOD(data_l2_fullres, n_pca_comp_feed, pca_max_iter, pca_err_tol, pca_sig_rem, verb)
            call update_status(status, 'pca_feed_filter')
         end if
         
         ! Write diagnostic l2_file to disc
         if (diag_l2) then
             call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
             
             ! Fit noise
             call fit_noise(data_l2_decimated)
             write(*,*) "Writing out l2-data for diagnostcs ", adjustl(trim("_6_after_PCA"))
             call write_l2_file(scan, k, data_l2_decimated, adjustl(trim("_6_after_PCA")))
          end if
         
         if (trim(pinfo%type) == 'gal') then
            call copy_lx_struct(data_l2_fullres, data_l2_filter)
            
            if (data_l2_fullres%use_freq_filter) then
               call frequency_filter_TOD(data_l2_fullres, nu_gain, alpha_gain, scan%ss(k))
            else
               call polyfilter_TOD(data_l2_fullres, bp_filter)
            end if
            !call polyfilter_TOD(data_l2_filter, bp_filter0)
            
            ! pca filter copied data
            call pca_filter_TOD(data_l2_filter, n_pca_comp, pca_max_iter, pca_err_tol, pca_sig_rem, verb)
         
            call remove_pca_components(data_l2_filter, data_l2_fullres, pca_sig_rem)
         end if
         
         ! Fourier transform frequency direction
         !call convert_GHz_to_k(data_l2_fullres(i))               
 
         ! Apply absolute calibration
         call calibrate_tod(data_l2_fullres, max_tsys)
         
         ! Update freqmask after max_tsys cut
         call update_freqmask(data_l2_fullres, min_acceptrate, scan%ss(k)%id, verb)
 
         if (verb) then
            write(*,*) "Average acceptrate for scan: ", scan%ss(k)%id, &
                 & sum(data_l2_fullres%freqmask_full) &
                 & / (size(data_l2_fullres%pixels, 1) - 1.d0) &
                 & / 4.d0 / 1024.d0 
         end if
         
         if (verb) then
            if (.not. all(data_l2_fullres%tod == data_l2_fullres%tod)) then
               write(*,*) "NaN in tod after filtering!"
            end if
         end if
         
         ! If necessary, decimate L2 file in both time and frequency
         
         call decimate_L2_data(samprate, numfreq, data_l2_fullres, data_l2_decimated)
 
         call update_status(status, 'decimate')
         !write(*,*) 'c'
 
         ! Fit noise
         if (.not. import_freqmask .and. .not. import_sigma) then 
            call fit_noise(data_l2_decimated)
         else
            call transfer_imported_sigma(data_l2_fullres, data_l2_decimated)
         end if
 
         ! Replace TOD with simulated data
         if (.false.) call simulate_gain_data(rng_handle, data_l2_decimated)
 
 
         ! Write L2 file to disk
         if (verb) then
            write(*,*) 'Writing ', scan%ss(k)%id, ' to disk2', trim(scan%ss(k)%l2file)
         end if
         data_l2_decimated%irun = irun - 1
         call mkdirs(trim(scan%ss(k)%l2file), .true.)
         call write_l2_file(scan, k, data_l2_decimated)
         call update_status(status, 'write_l2')
 
         ! Clean up data structures
         call free_lx_struct(data_l2_decimated)
         call free_lx_struct(data_l2_fullres)
         
      end do        
      call free_lx_struct(data_l1)
   end do
 
   call mpi_finalize(ierr)
 
   call free_status(status)
 
 contains

   subroutine prepare_sim_interpolation(simdata, pinfo)
      implicit none 
      type(simulation_struct),      intent(inout)    :: simdata
      type(patch_info),             intent(in)    :: pinfo
      real(dp), allocatable, dimension(:,:,:,:)   :: coeff
      integer(i4b)    :: i, sb, freq, feed, nfreq, nsb, ndet, nsamp, nx, ny

      nx        = size(simdata%simcube, 1)
      ny        = size(simdata%simcube, 2)
      nfreq     = size(simdata%simcube, 3)
      nsb       = size(simdata%simcube, 4)
      
      
      if (.not. allocated(coeff)) allocate(coeff(4, 4, nx, ny))
      if (.not. allocated(simdata%allcoeff)) allocate(simdata%allcoeff(4, 4, nx, ny, nsb, nfreq))
      if (.not. allocated(simdata%x)) allocate(simdata%x(nx))
      if (.not. allocated(simdata%y)) allocate(simdata%y(ny))

      simdata%allcoeff = 0.d0
      coeff = 0.d0
      ! NB! This way of adding signal only makes sense for co2 field (field 1) due to it being at Dec = 0 giving no
      ! projection effects. Add projection before using co6 or co7.

      do i = 1, nx
         simdata%x(i) = 0.5 * (simdata%edgex(i) + simdata%edgex(i + 1)) + pinfo%pos(1) 
      end do
      do i = 1, ny
         simdata%y(i) = 0.5 * (simdata%edgey(i) + simdata%edgey(i + 1)) + pinfo%pos(2)
      end do

      write(*,*) "Preparing simulation interpolation."

      !$OMP PARALLEL PRIVATE(coeff, sb, freq, i)
      !$OMP DO SCHEDULE(guided)    
      do sb = 1, nsb 
         do i = 1, nfreq
            freq = simdata%freqidx(sb, i)   ! Using flipped frequency index for sb = 1 and 3
            call splie2_full_precomp(simdata%y, simdata%x, simdata%boost * simdata%simcube(:, :, freq, sb), coeff)
            simdata%allcoeff(:, :, :, :, sb, i) = coeff
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      deallocate(coeff)

   end subroutine prepare_sim_interpolation

   subroutine sim2tod(data_l2, simdata)
      implicit none 
      type(Lx_struct),              intent(inout) :: data_l2
      type(simulation_struct),      intent(in)    :: simdata
      real(dp), allocatable, dimension(:) :: ra, dec

      integer(i4b)    :: i, sb, freq, feed, nfreq, nsb, ndet, nsamp, nx, ny 
      real(dp)        :: signal

      nsamp       = size(data_l2%tod,1)
      nfreq       = size(data_l2%freqmask_full,1)
      nsb         = size(data_l2%tod,3)
      ndet        = size(data_l2%tod,4)
      
      nx        = size(simdata%x, 1)
      ny        = size(simdata%y, 1)
      
      allocate(ra(nsamp))
      allocate(dec(nsamp))
      !$OMP PARALLEL PRIVATE(ra, dec, feed, i, sb, freq)
      !$OMP DO SCHEDULE(guided)    
      do feed = 1, ndet 
         if (.not. is_alive(data_l2%pixels(feed))) cycle
         !write(*,*) "Adding signal to feed:", feed
         ra = 0.d0
         dec = 0.d0
         do i = 1, nsamp

            ra(i)  = data_l2%point_cel(1, i, feed)
            dec(i) = data_l2%point_cel(2, i, feed)
         end do
         do sb = 1, nsb 

            if (sum(data_l2%freqmask_full(:, sb, feed)) == 0.d0) cycle
            do freq = 1, nfreq
               if (data_l2%freqmask_full(freq, sb, feed) == 0.d0) cycle
               if (data_l2%Tsys(freq, sb, feed) /= data_l2%Tsys(freq, sb, feed)) then 
                  write(*,*) "Tsys became NaN"
                  cycle
               end if

               if (data_l2%Tsys(freq, sb, feed) == 0) then
                  write(*,*) "Tsys is zero! No signal added!", freq, sb, feed 
                  !data_l2%tod(i, freq, sb, feed) = 0
                  cycle
               end if 

               do i = 1, nsamp
                  signal  = splin2_full_precomp(simdata%y, simdata%x, simdata%allcoeff(:, :, :, :, sb, freq), dec(i), ra(i))
                  data_l2%tod(i, freq, sb, feed) = data_l2%tod(i, freq, sb, feed) * (1 + signal / data_l2%Tsys(freq, sb, feed))
               end do

            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      deallocate(ra)
      deallocate(dec)   
   end subroutine sim2tod
   
   subroutine find_spikes(data_l2, verb, id)
     implicit none
     type(Lx_struct),     intent(inout) :: data_l2
     logical(lgt),        intent(in)    :: verb
     integer(i4b),        intent(in)    :: id
     real(dp)        :: gamma, cutoff, Tmean
     integer(i4b)    :: i, j, k, nsamp, nfreq, nsb, ndet, is_spike, n_spikes(4)
     real(dp), allocatable, dimension(:,:) :: ampsum
     
     nsamp       = size(data_l2%tod,1)
     nfreq       = size(data_l2%freqmask_full,1)
     nsb         = size(data_l2%tod,3)
     ndet        = size(data_l2%tod,4)
     
     if (.not. allocated(data_l2%spike_data)) allocate(data_l2%spike_data(1000,4,nsb,ndet,3))
     allocate(ampsum(nsb,ndet))
     data_l2%spike_data(:,:,:,:,:) = 0.d0
 
     n_spikes(:) = 0
     gamma = 0.7d0
     cutoff = 0.0015d0 * 8.d0
     ampsum(:,:) = 0.d0
     Tmean = 55.d0
 
     i = 0
     outer:do 
        i = i+1
        if (i > nsamp - 1) exit
        do k = 1, ndet
           if (.not. is_alive(data_l2%pixels(k))) cycle
           do j = 1, nsb
              if (sum(data_l2%freqmask_full(:,j,k)) == 0.d0) cycle
 
              if (data_l2%use_freq_filter) then                
                 if ((isnan(data_l2%T_cont(i+1,j,k,1))) .or. (isnan(data_l2%T_cont(i,j,k,1)))) cycle
                 ampsum(j,k) = ampsum(j,k) * gamma + data_l2%T_cont(i+1,j,k,1) / Tmean - data_l2%T_cont(i,j,k,1) / Tmean 
              else
                 if ((isnan(data_l2%tod_poly(i+1,0,j,k))) .or. (isnan(data_l2%tod_poly(i,0,j,k)))) cycle
                 ampsum(j,k) = ampsum(j,k) * gamma + data_l2%tod_poly(i+1,0,j,k) - data_l2%tod_poly(i,0,j,k) 
              end if
              if (abs(ampsum(j,k)) > cutoff) then
                 ! write(*,*) "Spike at:", ampsum(j,k), j, k, i + 1, data_l2%tod_poly(i+1,0,j,k), data_l2%tod_poly(i,0,j,k)
                 call get_spike_data(data_l2,k,j,i+1,n_spikes)
                 ampsum(:,:) = 0.d0
                 i = i + 10
                 cycle outer
              end if
           end do
        end do
     end do outer
     if (verb) then
        write(*,*) "Spike data for scan", id
        write(*,*) "Found ", n_spikes(1), " spikes"
        write(*,*) "Found ", n_spikes(2), " jumps"
        write(*,*) "Found ", n_spikes(3), " anomalies"
        write(*,*) "Found ", n_spikes(4), " edge spikes"

     end if
     deallocate(ampsum)
 
   end subroutine find_spikes
   
   subroutine get_spike_data(data, k, j, i, n_spikes)
     implicit none
     type(Lx_struct),       intent(inout)     :: data
     integer(i4b),          intent(in)        :: i, j, k
     integer(i4b),          intent(inout)     :: n_spikes(4)
     real(dp), allocatable, dimension(:,:,:)  :: fwd
     real(dp)        :: gamma, cutoff, Tmean
     integer(i4b)    :: m, n, l, nsamp, nfreq, nsb, ndet, spike_type, max_ind,indices(3)
     integer(i4b)    :: jump_mean_dur
 
     nsamp       = size(data%tod,1)
     nfreq       = size(data%freqmask_full,1)
     nsb         = size(data%tod,3)
     ndet        = size(data%tod,4)
     
     if ((i < 61) .or. (i > nsamp - 61)) then
        spike_type = 4
        n_spikes(spike_type) = n_spikes(spike_type) + 1
        if (n_spikes(spike_type) .le. 1000) then
 !          write(*,*) n_spikes(spike_type), spike_type, j, k, i
           data%spike_data(n_spikes(spike_type), spike_type, j, k, 1) = 1
           data%spike_data(n_spikes(spike_type), spike_type, j, k, 2) = i
           data%spike_data(n_spikes(spike_type), spike_type, j, k, 3) = data%time(i)
        end if
        return
     end if
     
     gamma = 0.7d0
     jump_mean_dur = 20
     allocate(fwd(41,nsb,ndet))
 
     Tmean = 55.d0
 
     
     fwd(:,:,:) = 0.d0
           
     do n = 1, ndet
        if (.not. is_alive(data%pixels(n))) cycle
        do l = 1, nsb
           do m = -20, 19
              if (data%use_freq_filter) then
                 if ((isnan(data%T_cont(i+m+1,l,n,1))) .or. (isnan(data%T_cont(i+m,l,n,1)))) cycle
                 fwd(21 + m + 1,l,n) = fwd(21 + m,l,n) * gamma + data%T_cont(i + m + 1,l,n,1) / Tmean - data%T_cont(i + m,l,n,1) / Tmean
              else
                 if ((isnan(data%tod_poly(i+m+1,0,l,n))) .or. (isnan(data%tod_poly(i+m,0,l,n)))) cycle
                 fwd(21 + m + 1,l,n) = fwd(21 + m,l,n) * gamma + data%tod_poly(i + m + 1,0,l,n) - data%tod_poly(i + m,0,l,n)
              end if
           end do
        end do
     end do
     indices = maxloc(abs(fwd(15:25,:,:)))
     max_ind = i + indices(1) - 7
     l = indices(2)
     n = indices(3)
 
 !    write(*,*) max_ind, l, n, fwd(max_ind - i + 21,l,n)
 
     if (any(fwd(max_ind - i + 21:max_ind - i + 7 + 21,l,n) * sign(1.d0,fwd(max_ind - i + 21,l,n)) < -0.0015d0 * 2)) then
        spike_type = 1
     else if (.not. any(data%sb_mean(max_ind - jump_mean_dur:max_ind - 10,l,n) &
          /= data%sb_mean(max_ind - jump_mean_dur:max_ind - 10,l,n))) then 
        if (abs(sum(data%sb_mean(max_ind - jump_mean_dur:max_ind - 10,l,n)) &
             - sum(data%sb_mean(max_ind + 10:max_ind + jump_mean_dur,l,n))) &
             / max(abs(sum(data%sb_mean(max_ind - jump_mean_dur:max_ind - 10,l,n))), 1.d-4) > 0.010d0) then
           spike_type = 2                 
        else
           spike_type = 3
        end if
        
     else    
        spike_type = 3               
     end if
     n_spikes(spike_type) = n_spikes(spike_type) + 1
           
     do n = 1, ndet
        if (.not. is_alive(data%pixels(n))) cycle
        do l = 1, nsb
           max_ind = i + maxloc(abs(fwd(15:25,l,n)), dim=1) - 7
           !write(*,*) max_ind
           if (n_spikes(spike_type) .le. 1000) then
              data%spike_data(n_spikes(spike_type), spike_type, l, n, 1) = fwd(max_ind - i + 21,l,n)
              data%spike_data(n_spikes(spike_type), spike_type, l, n, 2) = max_ind
              data%spike_data(n_spikes(spike_type), spike_type, l, n, 3) = data%time(max_ind)
           end if
        end do
     end do
     
 
     deallocate(fwd)
   end subroutine get_spike_data
 
   subroutine transfer_diagnostics(data_l2_in, data_l2_out)
     implicit none
     type(Lx_struct),                            intent(in) :: data_l2_in
     type(Lx_struct),                            intent(inout) :: data_l2_out
     integer(i4b) :: i, j, k, l, m, n, nsamp, nfreq, nsb, ndet
     
     nsamp       = size(data_l2_in%tod,1)
     nfreq       = size(data_l2_in%freqmask_full,1) ! nfreq in highres data
     nsb         = size(data_l2_in%tod,3)
     ndet        = size(data_l2_in%tod,4)
     
     if (.not. allocated(data_l2_out%freqmask_full)) allocate(data_l2_out%freqmask_full(nfreq,nsb,ndet))
     data_l2_out%freqmask_full = data_l2_in%freqmask_full
     if (.not. allocated(data_l2_out%freqmask_reason)) allocate(data_l2_out%freqmask_reason(nfreq,nsb,ndet))
     data_l2_out%freqmask_reason = data_l2_in%freqmask_reason
     if (.not. allocated(data_l2_out%diagnostics)) allocate(data_l2_out%diagnostics(ndet,nsb,nfreq,size(data_l2_in%diagnostics,4)))
     data_l2_out%diagnostics = data_l2_in%diagnostics
     if (.not. allocated(data_l2_out%cut_params)) allocate(data_l2_out%cut_params(size(data_l2_in%cut_params,1),size(data_l2_in%cut_params,2)))
     data_l2_out%cut_params = data_l2_in%cut_params
     if (.not. allocated(data_l2_out%spike_data)) allocate(data_l2_out%spike_data(size(data_l2_in%spike_data,1),size(data_l2_in%spike_data,2),&
             &size(data_l2_in%spike_data,3),size(data_l2_in%spike_data,4),size(data_l2_in%spike_data,5)))
     data_l2_out%spike_data = data_l2_in%spike_data
     if (.not. allocated(data_l2_out%AB_mask)) allocate(data_l2_out%AB_mask(nfreq,nsb,20))
     data_l2_out%AB_mask = data_l2_in%AB_mask
     if (.not. allocated(data_l2_out%leak_mask)) allocate(data_l2_out%leak_mask(nfreq,nsb,20))
     data_l2_out%leak_mask = data_l2_in%leak_mask
     
     if (data_l2_out%use_freq_filter) then
        if (.not. allocated(data_l2_out%corr_templ_ampl)) allocate(data_l2_out%corr_templ_ampl(4,nsb/2,ndet))
        data_l2_out%corr_templ_ampl = data_l2_in%corr_templ_ampl
     end if

         
     if (data_l2_out%import_freqmask) then
        if (.not. allocated(data_l2_out%freqmask)) allocate(data_l2_out%freqmask(ndet,nsb,nfreq))  !!! is this order of arguments correct?
        data_l2_out%freqmask = data_l2_in%freqmask
 
        if (.not. allocated(data_l2_out%acceptrate)) allocate(data_l2_out%acceptrate(nsb,ndet))
        data_l2_out%acceptrate = data_l2_in%acceptrate
 
        if (.not. allocated(data_l2_out%sigma0)) allocate(data_l2_out%sigma0(nfreq,nsb,ndet))
        data_l2_out%sigma0 = data_l2_in%sigma0
       
        if (.not. allocated(data_l2_out%alpha)) allocate(data_l2_out%alpha(nfreq,nsb,ndet))
        data_l2_out%alpha  = data_l2_in%alpha
 
        if (.not. allocated(data_l2_out%fknee)) allocate(data_l2_out%fknee(nfreq,nsb,ndet))
        data_l2_out%fknee  = data_l2_in%fknee
 
        if (.not. allocated(data_l2_out%var_fullres)) allocate(data_l2_out%var_fullres(nfreq,nsb,ndet))
        data_l2_out%var_fullres = data_l2_in%var_fullres
     end if 
 
     call free_lx_struct(data_l2_in)
 
   end subroutine transfer_diagnostics
   
   subroutine transfer_imported_sigma(data_l2_in, data_l2_out)
     implicit none
     type(Lx_struct),                            intent(in)    :: data_l2_in
     type(Lx_struct),                            intent(inout) :: data_l2_out
     integer(i4b) :: i, j, k, l, m, n, nsamp, nfreq, nsb, ndet
     nsamp       = size(data_l2_in%tod,1)
     nfreq       = size(data_l2_in%freqmask_full,1) ! nfreq in lowres data
     nsb         = size(data_l2_in%tod,3)
     ndet        = size(data_l2_in%tod,4)
     

     if (.not. allocated(data_l2_out%sigma0)) allocate(data_l2_out%sigma0(nfreq,nsb,ndet))
     data_l2_out%sigma0 = data_l2_in%sigma0
 
     if (.not. allocated(data_l2_out%alpha)) allocate(data_l2_out%alpha(nfreq,nsb,ndet))
     data_l2_out%alpha  = data_l2_in%alpha
 
     if (.not. allocated(data_l2_out%fknee)) allocate(data_l2_out%fknee(nfreq,nsb,ndet))
     data_l2_out%fknee  = data_l2_in%fknee
 
     if (.not. allocated(data_l2_out%var_fullres)) allocate(data_l2_out%var_fullres(nfreq,nsb,ndet))
     data_l2_out%var_fullres = data_l2_in%var_fullres
 
     call free_lx_struct(data_l2_in)
   end subroutine 
 
   subroutine update_freqmask(data_l2, min_acceptrate, id, verb)
     implicit none
     type(Lx_struct),                            intent(inout) :: data_l2
     integer(i4b),                               intent(in)    :: id
     real(dp),                                   intent(in)    :: min_acceptrate
     logical(lgt),                               intent(in)    :: verb
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
              data_l2%freqmask_reason(:,j,i) = data_l2%freqmask_reason(:,j,i) + nint(15.d0 * data_l2%freqmask_full(:,j,i))
              data_l2%freqmask_full(:,j,i) = 0.d0
              if (verb) then
                 write(*,*) "Rejecting entire sideband (too much was masked) det, sb, acceptrate, scanid:"
                 write(*,*) data_l2%pixels(i), j, data_l2%acceptrate(j,i), id
              end if
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
 
   subroutine mask_specific_corr(data_l2, vars, means, i, j, k, l, m, n, cut, id, reason)
     implicit none
     type(Lx_struct),                            intent(inout) :: data_l2
     integer(i4b),                               intent(in)    :: id, reason
     real(dp), allocatable, dimension(:, :, :),  intent(in)    :: means, vars
     integer(i4b),                               intent(in)    :: i, j, k, l, m, n           
     real(dp),                                   intent(in)    :: cut
     real(dp)                                                  :: corr
     integer(i4b)                                              :: nsamp
     if (vars(k,j,i) < 1.d-14) return
     if (vars(n,m,l) < 1.d-14) return
     nsamp       = size(data_l2%tod,1)
     
     corr = sum((data_l2%tod(:,k,j,i) - means(k,j,i)) * (data_l2%tod(:,n,m,l) - means(n,m,l))) / nsamp
     ! if (.not. (vars(k,j,i) > 0)) then
     !    write(*,*) vars(k,j,i), k,j,i
     ! end if
     ! if (.not. (vars(n,m,l) > 0)) then
     !    write(*,*) vars(n,m,l), n,m,l
     ! end if
     
     corr = corr / sqrt(vars(k,j,i) * vars(n,m,l))
     
     if (corr > cut) then
        data_l2%freqmask_full(k,j,i) = 0.d0
        data_l2%freqmask_full(n,m,l) = 0.d0
        data_l2%freqmask_reason(k,j,i) = reason
        data_l2%freqmask_reason(n,m,l) = reason
      end if

   end subroutine mask_specific_corr
 
   subroutine remove_corr_templates(corrs, templates, ampl)
     ! Fits and removes the templates from the correlation matrix by linear regression
     ! Model: y = X ampl + noise, solution: ampl = (X^TX)^(-1)X^Ty, where y is corr, and X are
     ! the templates. 
     implicit none
     real(sp), allocatable, dimension(:, :),     intent(inout) :: corrs
     real(sp), allocatable, dimension(:,:,:),    intent(in)    :: templates
     ! real(dp), allocatable, dimension(:,:),      intent(in)    :: mask 
 !    integer(i4b),                               intent(in)    :: feed, band
     real(dp), allocatable, dimension(:)                       :: xty
     real(dp), allocatable, dimension(:),        intent(out)   :: ampl
     real(dp), allocatable, dimension(:,:)                     :: xtx
     integer(i4b)                                              :: nsamp, i, j, k, l, m, n, ntemp
     
     ntemp = size(templates, 3)  ! check this
     if(.not. allocated(ampl)) allocate(ampl(ntemp)) 
     allocate(xty(ntemp))
     allocate(xtx(ntemp,ntemp))
 
     if (sum(merge(1.d0,0.d0,corrs /= 0.d0)) == 0.0) then
        return
     end if
     
     do i = 1, ntemp
        xty(i) = sum(corrs(:,:) * templates(:,:,i))
        xtx(i,i) = sum(templates(:,:,i) * templates(:,:,i) * merge(1.d0,0.d0,corrs /= 0.d0))  ! last part is basically a mask
        do j = i+1, ntemp
           xtx(i,j) = sum(templates(:,:,i) * templates(:,:,j) * merge(1.d0,0.d0,corrs /= 0.d0))
           xtx(j,i) = xtx(i,j)
        end do
     end do
     
     if (xtx(3,1) == 0.0) then
        write(*,*) corrs(59, 321)
        write(*,*) templates(123, 413, 3)
     end if
     
 !    write(*,*) xtx
     call invert_matrix(xtx)
     
        
     
     do i = 1, ntemp
        ampl(i) = sum(xtx(:,i) * xty(:))
        corrs(:,:) = corrs(:,:) - ampl(i) * templates(:,:,i) * merge(1.d0,0.d0,corrs /= 0.d0)
     end do
 !    write(*,*) ampl
     deallocate(xty)
     deallocate(xtx)
 
   end subroutine remove_corr_templates
 
   
   subroutine get_pca_corr_templates(ampl, Delta, freqmask)
      ! Fits and removes the correlation pattern induced by Feed-PCA filter.
      implicit none
      real(sp), dimension(:, :, :),               intent(in)    :: ampl      !(freq, sb, comp)
      real(sp), dimension(:, :),                  intent(in)    :: freqmask  !(sb, freq)
      real(sp),  dimension(:,:),                intent(inout)   :: Delta
      real(sp), allocatable, dimension(:, :)                    :: B          !(nsb * nfreq, comp) Basis matrix
      ! real(sp), allocatable, dimension(:)                       :: magnitudes !(comp) magnitudes of ampl
      real(sp)                                                  :: magnitudes ! magnitudes of ampl
      real(sp), allocatable, dimension(:,:)                     :: L
      
      integer(i4b)                                              :: i, j, k, n, m, o, p, idx1, idx2, nfreq, nsb, nchl, ncomp, nu, comp
      
      nchl  = size(ampl, 1) ! Number of channels per sideband 
      nsb   = size(ampl, 2) ! Number of sidebands
      ncomp = size(ampl, 3) ! Number of subtracted PCA components

      nfreq  = nsb * nchl  ! Total number of frequency channels
      if (.not. allocated(B)) allocate(B(nfreq, ncomp))  
      if (.not. allocated(L)) allocate(L(nfreq, nfreq))
      !if (.not. allocated(magnitudes)) allocate(magnitudes(ncomp))
      
      B   = 0.d0
      L = 0.d0
      magnitudes = 1.d0


      ! Normalizing and flattening frequency amplitudes.
      
      do i = 1, ncomp
         magnitudes =  sqrt(sum(sum(ampl(:, :, i) * ampl(:, :, i), 1), 1))
         if (magnitudes == 0.d0) cycle
         if (magnitudes .ne. magnitudes) then
            cycle
         end if
         do j = 1, nsb
            do k = 1, nchl 
               if (freqmask(k, j) == 0.d0) cycle
               B((j - 1) * nchl + k, i) = ampl(k, j, i) / magnitudes
            end do
         end do
      end do

      do i = 1, nfreq
         L(i, i) = 1.0  ! Filling diagonal with 1, i.e. the identity
      end do

      call sgemm("N", "T", nfreq, nfreq, ncomp, -1.0, B, nfreq, B, nfreq, 1.0, L, nfreq)  ! BLAS computes L := I - LL^T
   
      do i = 1, nfreq
         Delta(i, i) = 1.0 ! Filling diagonal with 1, i.e. the identity
      end do

      call sgemm("N", "T", nfreq, nfreq, nfreq, -1.0, L, nfreq, L, nfreq, 1.0, Delta, nfreq) ! BLAS computes: Delta := I - LL^T

      deallocate(B)
      deallocate(L)
      ! deallocate(magnitudes)

   end subroutine get_pca_corr_templates


   subroutine flag_correlations(data_l2, id, parfile, rank)!corr_cut, mean_corr_cut, mean_abs_corr_cut, median_cut, var_cut, n_neighbor, neighbor_factor, var_max, corr_max)
     implicit none
     type(Lx_struct),          intent(inout) :: data_l2
     integer(i4b),             intent(in)    :: id, rank
     type(hdf_file)                          :: file
     character(len=512),       intent(in)    :: parfile
     integer(i4b) :: i, j, k, l, m, n, o, p, p2, q, q2, pp, qq, kk, nn, reason
     integer(i4b) :: nsamp, nfreq, nsb, ndet, n_offset, dof, nprod, jump, unit
     real(dp)     :: corr, mean, var_0, std_median, dnu, radiometer, prod, box_sum, edge_corr_cut
     real(dp)     :: mean_meancorr, sigma_meancorr, mean_vars, sigma_vars, mean_maxcorr, sigma_maxcorr
     real(dp)     :: mean_meanabscorr, sigma_meanabscorr, mean_corrsum, sigma_corrsum
     real(dp)     :: corr_cut, mean_corr_cut, var_cut, median_cut, nsigma_edge_corrs
     real(dp)     :: mean_abs_corr_cut, neighbor_factor, var_max, corr_max
     real(dp),     allocatable, dimension(:, :, :)   :: means, vars, maxcorr, meancorr, meanabscorr
     real(dp),     allocatable, dimension(:, :, :)   :: temp_freqmask !corrsum
     real(sp),     allocatable, dimension(:, :, :)   :: corrsum_mask
     real(dp),     allocatable, dimension(:)         :: subt, median, smedian, minimum, ampl
     real(sp),     allocatable, dimension(:,:)       :: corrs, corr_prod, tod_subset
     real(sp),     allocatable, dimension(:,:,:)     :: corr_templates
     real(sp),     allocatable, dimension(:,:)       :: corr_template
     real(sp),     allocatable, dimension(:,:)       :: corr_template_pca
     real(dp),     allocatable, dimension(:,:)       :: outlier_mask
     logical(lgt) :: mask_edge_corrs, rm_outliers, verb
     character(len=512) :: box_offset_str, stripe_offset_str, nsigma_prod_stripe_str
     character(len=512) :: nsigma_prod_box_str, nsigma_mean_box_str, aliasing_filename, filename, template_name
     real(dp)     :: nsigma_chi2_box, nsigma_chi2_stripe, aliasing_db_cutoff
     integer(i4b) :: n_neighbor, prod_offset
     integer(i4b), dimension(3) :: box_offsets, stripe_offsets
     real(dp),     dimension(3) :: nsigma_prod_box, nsigma_prod_stripe, nsigma_mean_box
     
     nsamp       = size(data_l2%tod,1)
     nfreq       = size(data_l2%tod,2)
     nsb         = size(data_l2%tod,3)
     ndet        = size(data_l2%tod,4)
 
     call get_parameter(unit, parfile, 'CORRELATION_CUT',           par_dp=corr_cut)
     call get_parameter(unit, parfile, 'MEAN_CORRELATION_CUT',      par_dp=mean_corr_cut)
     call get_parameter(unit, parfile, 'VARIANCE_CUT',              par_dp=var_cut)
     call get_parameter(unit, parfile, 'MEAN_ABS_CORRELATION_CUT',  par_dp=mean_abs_corr_cut)
     call get_parameter(unit, parfile, 'MEDIAN_CUT',                par_dp=median_cut)
     call get_parameter(unit, parfile, 'N_NEIGHBOR',                par_int=n_neighbor)
     call get_parameter(unit, parfile, 'VARIANCE_MAX',              par_dp=var_max)
     call get_parameter(unit, parfile, 'CORRELATION_MAX',           par_dp=corr_max)
     call get_parameter(unit, parfile, 'NEIGHBOR_FACTOR',           par_dp=neighbor_factor)
 !    call get_parameter(unit, parfile, 'NSIGMA_EDGE_CORRS',         par_dp=nsigma_edge_corrs)
     call get_parameter(unit, parfile, 'ALIASING_FILE',             par_string=aliasing_filename)
     call get_parameter(unit, parfile, 'ALIASING_DB_CUT',           par_dp=aliasing_db_cutoff)
     call get_parameter(unit, parfile, 'MASK_EDGE_CORRS',           par_lgt=mask_edge_corrs)
     call get_parameter(unit, parfile, 'REMOVE_OUTLIERS',           par_lgt=rm_outliers)
     call get_parameter(unit, parfile, 'BOX_OFFSETS',               par_string=box_offset_str)
     call get_parameter(unit, parfile, 'STRIPE_OFFSETS',            par_string=stripe_offset_str)
     call get_parameter(unit, parfile, 'NSIGMA_PROD_BOX',           par_string=nsigma_prod_box_str)
     call get_parameter(unit, parfile, 'NSIGMA_PROD_STRIPE',        par_string=nsigma_prod_stripe_str)
     call get_parameter(unit, parfile, 'NSIGMA_MEAN_BOX',           par_string=nsigma_mean_box_str)
     call get_parameter(unit, parfile, 'PROD_OFFSET',               par_int=prod_offset)
     call get_parameter(unit, parfile, 'NSIGMA_CHI2_BOX',           par_dp=nsigma_chi2_box)
     call get_parameter(unit, parfile, 'NSIGMA_CHI2_STRIPE',        par_dp=nsigma_chi2_stripe)
     call get_parameter(unit, parfile, 'VERBOSE_PRINT',             par_lgt=verb)
 
     read(box_offset_str,*) box_offsets
     read(stripe_offset_str,*) stripe_offsets
     read(nsigma_prod_box_str,*) nsigma_prod_box
     read(nsigma_prod_stripe_str,*) nsigma_prod_stripe
     read(nsigma_mean_box_str,*) nsigma_mean_box
 
     allocate(means(nfreq, nsb, ndet), vars(nfreq, nsb, ndet))
     allocate(maxcorr(nfreq, nsb, ndet), meancorr(nfreq, nsb, ndet))
     allocate(meanabscorr(nfreq, nsb, ndet))
     if(.not. allocated(data_l2%diagnostics)) allocate(data_l2%diagnostics(nfreq,nsb,ndet,5)) 
     if(.not. allocated(data_l2%cut_params)) allocate(data_l2%cut_params(2,5)) 
 !    allocate(corrs(nfreq,nfreq), corr_template(nfreq,nfreq))
     allocate(corrs(2 * nfreq, 2 * nfreq))
 !    allocate(corr_prod(2 * nfreq-1, 2 * nfreq))
     allocate(subt(nsamp-1), median(nfreq))
     allocate(temp_freqmask(nfreq,nsb,ndet))
     allocate(tod_subset(nsamp,2*nfreq))
 
     data_l2%diagnostics = 0.d0
     
     means = 0.d0
     vars = 0.d0
     maxcorr = 0.d0
     meancorr = 0.d0
     meanabscorr = 0.d0
     temp_freqmask = 1.d0
 
     if (data_l2%use_freq_filter) then
        allocate(corr_templates(2 * nfreq, 2 * nfreq, 4))
        if(.not. allocated(data_l2%corr_templ_ampl)) allocate(data_l2%corr_templ_ampl(4, nsb / 2, ndet)) 
     else if (data_l2%polyorder == 1) then
        allocate(corr_template(nfreq,nfreq))
 
        call open_hdf_file("/mn/stornext/d22/cmbco/comap/protodir/auxiliary/corr_template.h5", file, "r")
        call read_hdf(file, "corr", corr_template)
        call close_hdf_file(file)
     end if
 
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, nfreq
              if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
              means(k,j,i) = sum(data_l2%tod(:,k,j,i)) / nsamp
              vars(k,j,i) = sum(data_l2%tod(:,k,j,i) ** 2) / nsamp - means(k,j,i) ** 2
           end do
        end do
     end do
 
 !     allocate(corrs(4 * nfreq, 4 * nfreq))
 !     corrs = 0.d0
 !     i = 6
 !     l = i
 !     !!$OMP PARALLEL PRIVATE(j,k,l,m,n,p,q,corr)
 !     !!$OMP DO SCHEDULE(guided)    
 !     do p = 1, nsb * nfreq
 !        k = mod((p-1), nfreq) + 1
 !        j = (p-1) / nfreq + 1
 ! !       write(*,*) k, j, p
 !        if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
 !        do q = p+1, nsb * nfreq
 !           n = mod((q-1), nfreq) + 1
 !           m = (q-1) / nfreq + 1 
 ! !          write(*,*) k, j, n, m, p, q
 !           if (data_l2%freqmask_full(n,m,l) == 0.d0) cycle
 !           !write(*,*) k, j, n, m, p, q
 !           corr = sum((data_l2%tod(:,k,j,i) - means(k,j,i)) * (data_l2%tod(:,n,m,l) - means(n,m,l))) / nsamp
 ! !          write(*,*) k, j, n, m, p, q, corr          
 !           corr = corr / sqrt(vars(k,j,i) * vars(n,m,l))
 !           corrs(p,q) = corr
 !           corrs(q,p) = corr
 !        end do
 !     end do
 !     !!$OMP END DO
 !     !!$OMP END PARALLEL
 !     open(22, file="corr_06.unf", form="unformatted") ! Adjusted open statement
 !     write(22) corrs
 !     close(22)
     
 !     deallocate(corrs)
 !     allocate(corrs(2 * nfreq, 2 * nfreq))
 
     if (mask_edge_corrs) then
        
        call open_hdf_file(aliasing_filename, file, "r")
        if(.not. allocated(data_l2%AB_mask)) allocate(data_l2%AB_mask(nfreq,nsb,20)) 
        if(.not. allocated(data_l2%leak_mask)) allocate(data_l2%leak_mask(nfreq,nsb,20)) 

        
        ! Read telescope coordinates
        call read_hdf(file, "AB_mask",            data_l2%AB_mask)
        call read_hdf(file, "leak_mask",          data_l2%leak_mask)
 
        
        call close_hdf_file(file)
 
        do i = 1, ndet
           if (.not. is_alive(data_l2%pixels(i))) cycle
           do j = 1, nsb
              do k = 1, nfreq
                 if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
                 
                 if (data_l2%AB_mask(k,j,data_l2%pixels(i)) < aliasing_db_cutoff) then
                    data_l2%freqmask_full(k,j,i) = 0.d0
                    data_l2%freqmask_reason(k,j,i) = 40
                 end if
 
                 if (data_l2%leak_mask(k,j,data_l2%pixels(i)) < aliasing_db_cutoff) then
                    data_l2%freqmask_full(k,j,i) = 0.d0
                    data_l2%freqmask_reason(k,j,i) = 41
                 end if
              end do
           end do
        end do
        vars = vars * data_l2%freqmask_full       
     end if
     !    allocate(box_offsets(3), prod_offsets(3))
     
     if (data_l2%n_pca_comp_feed > 0) write(*,*) "Cleaning PCA template from correlation matrix",  " Rank:", rank, "id:", id
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle

        if (data_l2%n_pca_comp_feed > 0) then
            if (.not. allocated(corr_template_pca)) allocate(corr_template_pca(nsb * nfreq, nsb * nfreq))
            corr_template_pca = 0.d0
            if (sum(data_l2%pca_ampl_feed(:, :, i, :)) /= 0) then
               call get_pca_corr_templates(data_l2%pca_ampl_feed(:, :, i, :), corr_template_pca, data_l2%freqmask_full(:,:,i))         
            end if
        end if
      
      do o = 1, nsb / 2
         corrs = 0.d0
         tod_subset = 0.d0
         do p = 1, 2 * nfreq
            k = mod((p-1), nfreq) + 1
            j = (o-1) * 2 + 1 + (p-1) / nfreq
            tod_subset(:,p) = data_l2%tod(:,k,j,i) - means(k,j,i)
         end do
         call sgemm("T", "N", 2*nfreq, 2*nfreq, nsamp, 1.0, tod_subset, nsamp, tod_subset, nsamp, 0.0, corrs, 2*nfreq)


         do p = 1, 2 * nfreq
            k = mod((p-1), nfreq) + 1
            j = (o-1) * 2 + 1 + (p-1) / nfreq
            if (data_l2%freqmask_full(k,j,i) == 0.d0) then
               corrs(p,:) = 0.d0
               corrs(:,p) = 0.d0
               cycle
            end if
            do q = p + 2, 2 * nfreq
               m = (o-1) * 2 + 1 + (q-1) / nfreq
               n = mod((q-1), nfreq) + 1
               if (data_l2%freqmask_full(n,m,i) == 0.d0) cycle
               corrs(p,q) = corrs(p,q) / sqrt(vars(k,j,i) * vars(n,m,i))
               corrs(p,q) = corrs(p,q) / nsamp
               if (.not. data_l2%use_freq_filter .and. ((j == m) .and. (data_l2%polyorder == 0))) then 
                  corrs(p,q) = corrs(p,q) + 1.d0 / (nfreq - 3)  ! -3 because we mask 3 freqs always before poly-filter
               end if
               corrs(q,p) = corrs(p,q)
            end do
         end do
         do p = 1, 2*nfreq
            corrs(p,p) = 0.d0
         end do
         do p = 1, 2*nfreq - 1
            corrs(p,p+1) = 0.d0
            corrs(p+1,p) = 0.d0
         end do


           if (data_l2%n_pca_comp_feed > 0) then
               if (.false.) then
                  write(filename, "(A, I0.3, A, I0.3, A, I0.3, A)") 'corr_before_', data_l2%pixels(i), '_', o, '_', id,'_sp.unf' 
                  open(22, file=filename, form="unformatted") ! Adjusted open statement
                  write(22) corrs
                  close(22)
               end if
               corrs = corrs + corr_template_pca((o - 1) * 2 * nfreq + 1 : o * 2 * nfreq, &
                                               & (o - 1) * 2 * nfreq + 1 : o * 2 * nfreq)
               if (.false.) then
                  write(filename, "(A, I0.3, A, I0.3, A, I0.3, A)") 'corr_after_', data_l2%pixels(i), '_', o, '_', id,'_sp.unf' 
                  open(22, file=filename, form="unformatted") ! Adjusted open statement
                  write(22) corrs
                  close(22)
               end if
           end if
           
           ! if (i == 12 .and. o == 1) then
           !    open(22, file="corr_12_12_before.unf", form="unformatted") ! Adjusted open statement
           !    write(22) corrs
           !    close(22)
           ! end if
           if (.false.) then
              write(filename, "(A, I0.3, A, I0.3, A, I0.3, A)") 'corr_', data_l2%pixels(i), '_', o, '_', id,'_sp.unf' 
              open(22, file=filename, form="unformatted") ! Adjusted open statement
              write(22) corrs
              close(22)
           end if
 
           if (data_l2%use_freq_filter) then
              
              call open_hdf_file("/mn/stornext/d22/cmbco/comap/protodir/auxiliary/corr_templates_band.h5", file, "r")
              if (o == 1) then
                 write(template_name, "(A, I0.2, A)") 'corr_templates_', data_l2%pixels(i), '_A' 
              else 
                 write(template_name, "(A, I0.2, A)") 'corr_templates_', data_l2%pixels(i), '_B' 
              end if
              !write(*,*) template_name, id
              call read_hdf(file, template_name, corr_templates)
              call close_hdf_file(file)
 
              
              call remove_corr_templates(corrs, corr_templates, ampl)
              
              data_l2%corr_templ_ampl(:,o,i) = ampl(:)
              !write(*,*) data_l2%corr_templ_ampl(:,o,i)
              !corrs(1:nfreq,1:nfreq) = corrs(1:nfreq,1:nfreq) - corr_template * merge(1.d0,0.d0,corrs(1:nfreq,1:nfreq) /= 0.d0)
              !corrs(nfreq+1:,nfreq+1:) = corrs(nfreq+1:,nfreq+1:) - corr_template * merge(1.d0,0.d0,corrs(nfreq+1:,nfreq+1:) /= 0.d0)
           else if (data_l2%polyorder == 1) then
              corrs(1:nfreq,1:nfreq) = corrs(1:nfreq,1:nfreq) - corr_template * merge(1.d0,0.d0,corrs(1:nfreq,1:nfreq) /= 0.d0)
              corrs(nfreq+1:,nfreq+1:) = corrs(nfreq+1:,nfreq+1:) - corr_template * merge(1.d0,0.d0,corrs(nfreq+1:,nfreq+1:) /= 0.d0)
           end if
 
           if (.false.) then
              write(filename, "(A, I0.3, A, I0.3, A, I0.3, A)") 'corr_after_', data_l2%pixels(i), '_', o, '_', id,'_sp.unf' 
              open(22, file=filename, form="unformatted") ! Adjusted open statement
              write(22) corrs
              close(22)
           end if
 
           
           ! if (i == 12 .and. o == 1) then
           !    open(22, file="corr_12_12.unf", form="unformatted") ! Adjusted open statement
           !    write(22) corrs
           !    close(22)
           ! end if
 
           do p = 1, 2 * nfreq
              j = (o-1) * 2 + 1 + (p-1) / nfreq
              k = mod((p-1), nfreq) + 1
              if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
              if (sum(merge(1.d0,0.d0,corrs(p,:) /= 0.d0)) == 0.d0) cycle
              meancorr(k,j,i)    = sum(corrs(p,:)) / sum(merge(1.d0,0.d0,corrs(p,:) /= 0.d0))
              maxcorr(k,j,i)     = maxval(abs(corrs(p,:)))
              meanabscorr(k,j,i) = sum(abs(corrs(p,:))) / sum(merge(1.d0,0.d0,corrs(p,:) /= 0.d0))
           end do
           
           !!! Boxes
 !          box_offsets = (/32, 128, 512/)  ! must divide 1024
 !          prod_offsets = 16 !box_offsets / 2 + 1
           do l = 1, 3
              n_offset = box_offsets(l)
              jump = prod_offset
              !write(*,*) l, i, o
              do pp = 1, (2 * nfreq)/n_offset
                 p = (pp-1) * n_offset + 1
                 j = (o-1) * 2 + 1 + (p-1) / nfreq
                 k = mod((p-1), nfreq) + 1
                 p2 = pp * n_offset
                 kk = mod((p2-1), nfreq) + 1
                 do qq = pp + 1, (2 * nfreq)/n_offset
                    q = (qq-1) * n_offset + 1
                    m = (o-1) * 2 + 1 + (q-1) / nfreq
                    n = mod((q-1), nfreq) + 1
                    q2 = qq * n_offset
                    nn = mod((q2-1), nfreq) + 1
                    box_sum = sum(corrs(p:p2, q:q2) * sqrt(1.d0 * nsamp))
                    corr = sum(corrs(p:p2, q:q2) ** 2 * nsamp) 
                    dof = sum(merge(1,0,corrs(p:p2, q:q2) /= 0.d0))
                    prod = sum(corrs(p:p2-jump:2, q:q2) * corrs(p+jump:p2:2, q:q2) * nsamp)
                    nprod = sum(merge(1,0,corrs(p:p2-jump:2, q:q2) * corrs(p+jump:p2:2, q:q2) /= 0.d0))
                    if (nprod == 0) nprod = 1
                    if (dof == 0) dof = 1
                    if (corr - dof > nsigma_chi2_box * sqrt(2.d0 * dof)) then
                    !   if (i == 10 .and. o == 2) then
                    !      write(*,*) (corr - n_offset ** 2) / sqrt(2.d0 * n_offset ** 2), l, pp, qq, k, kk, n, nn, j, m
                    !   end if
                    
                       temp_freqmask(k:kk,j,i) = 0.d0
                       data_l2%freqmask_reason(k:kk,j,i) = 15 + l
                       temp_freqmask(n:nn,m,i) = 0.d0
                       data_l2%freqmask_reason(n:nn,m,i) = 15 + l
                    else if (prod / nprod > nsigma_prod_box(l) * sqrt(1.d0 / nprod)) then
                       !write(*,*) prod * sqrt(1.d0 / nprod), l, "box", i, j, k
                       temp_freqmask(k:kk,j,i) = 0.d0
                       data_l2%freqmask_reason(k:kk,j,i) = 21 + l
                       !data_l2%freqmask_full(n * n_offset:(n+1) * n_offset,j,i) = 0.d0
                       temp_freqmask(n:nn,m,i) = 0.d0
                       data_l2%freqmask_reason(n:nn,m,i) = 21 + l
                    else if (abs(box_sum / dof) > nsigma_mean_box(l) * sqrt(1.d0 / dof)) then
                       !write(*,*) abs(box_sum / dof) / sqrt(1.d0 / dof), l, "box_mean", i, j, k
                       temp_freqmask(k:kk,j,i) = 0.d0
                       data_l2%freqmask_reason(k:kk,j,i) = 30 + l
                       temp_freqmask(n:nn,m,i) = 0.d0
                       data_l2%freqmask_reason(n:nn,m,i) = 30 + l
                    end if
                 end do
              end do
           end do
           
           !!!! Stripes
           do l = 1, 3
              n_offset = stripe_offsets(l)
              jump = prod_offset
              do pp = 1, (2 * nfreq)/n_offset
                 corr = 0.d0
                 dof = 0
                 prod = 0.d0
                 nprod = 0
                 p = (pp-1) * n_offset + 1
                 j = (o-1) * 2 + 1 + (p-1) / nfreq
                 k = mod((p-1), nfreq) + 1
                 p2 = pp * n_offset
                 kk = mod((p2-1), nfreq) + 1
                 if (pp > 1) then 
                    corr = corr + sum(corrs(p:p2, :p-1) ** 2 * nsamp)
                    dof = dof + sum(merge(1,0,corrs(p:p2, :p-1) /= 0.d0))
                    prod = sum(corrs(p:p2, :p-1-jump:2) * corrs(p:p2, 1+jump:p-1:2) * nsamp)
                    nprod = sum(merge(1,0,corrs(p:p2, :p-1-jump:2) * corrs(p:p2, 1+jump:p-1:2) /= 0.d0))
                 end if
                 ! if (i == 12 .and. o == 1) then
                 !    write(*,*) "stripes, first"
                 !    if (nprod > 0) then
                 !       write(*,*) prod * sqrt(1.d0 / nprod), nprod, l, pp, k, kk, j
                 !       write(*,*) p, p2, p-2
                 !       write(*,*) sum(merge(1,0,corrs(p:p2, :p-2:jump) * corrs(p:p2, 2:p-1:jump) /= 0.d0))
                 !    end if
                 ! end if
                 
                 corr = corr + sum(corrs(p:p2,p:p2) ** 2 * nsamp) / 2.d0 
                 dof = dof + sum(merge(1,0,corrs(p:p2,p:p2) /= 0.d0)) / 2
                 prod = prod + sum(corrs(p:p2-jump:2, p:p2-jump:2) * corrs(p+jump:p2:2, p+jump:p2:2) * nsamp) / 2.d0
                 nprod = nprod + sum(merge(1,0,corrs(p:p2-jump:2, p:p2-jump:2) * corrs(p+jump:p2:2, p+jump:p2:2) /= 0.d0)) / 2
 
                 if (pp < (2 * nfreq)/n_offset) then
                    corr = corr + sum(corrs(p:p2, p2+1:) ** 2 * nsamp)
                    dof = dof + sum(merge(1,0,corrs(p:p2, p2+1:) /= 0.d0))
                    prod = prod + sum(corrs(p:p2, p2+1 + jump::2) * corrs(p:p2, p2+1:2 * nfreq - jump:2) * nsamp)
                    nprod = nprod + sum(merge(1,0,corrs(p:p2, p2+1 + jump::2) * corrs(p:p2, p2+1:2 * nfreq - jump:2) /= 0.d0))
                 end if
                 if (nprod == 0) nprod = 1
                 if (corr - dof > nsigma_chi2_stripe * sqrt(2.d0 * dof)) then
                    !   if (i == 10 .and. o == 2) then
                    !      write(*,*) (corr - n_offset ** 2) / sqrt(2.d0 * n_offset ** 2), l, pp, qq, k, kk, n, nn, j, m
                    !   end if
 
                    temp_freqmask(k:kk,j,i) = 0.d0
                    data_l2%freqmask_reason(k:kk,j,i) = 18 + l
                 else if (prod / nprod > nsigma_prod_stripe(l) * sqrt(1.d0 / nprod)) then
                    !write(*,*) prod * sqrt(1.d0 / nprod), l, "stripe", i, j, k
                    temp_freqmask(k:kk,j,i) = 0.d0
                    data_l2%freqmask_reason(k:kk,j,i) = 24 + l
                 end if
              end do
           end do
        end do
        !write(*,*) "Done with correlation masking in detector", i
     end do
     
     ! Thresholds for masking frequencies (radiometer noise is null-hypothesis
     dnu = (data_l2%nu(2, 1, 1) - data_l2%nu(3, 1, 1)) * 1d9
     radiometer = 1.0d0 / sqrt(dnu * 1.d0 / data_l2%samprate)
     !write(*,*) radiometer, dnu, data_l2%samprate
 
     ! Hard cuts
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, nfreq
              if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
              if (vars(k,j,i) > (radiometer * var_max) ** 2) then
                 data_l2%freqmask_full(k,j,i) = 0.d0
                 data_l2%freqmask_reason(k,j,i) = 3
              else if (maxcorr(k,j,i) > corr_max) then
                 data_l2%freqmask_full(k,j,i) = 0.d0
                 data_l2%freqmask_reason(k,j,i) = 4
              end if
           end do
        end do
     end do
     ! if (sum(data_l2%freqmask_full) == 0.0) then
     !    if (verb) then
     !       write(*,*) "All frequencies masked after hard cut, id: ", id
     !    end if
     !    ! fix this to "stop working on this file"
     ! end if
 
     ! convert to relative variance
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
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
 
     call get_mean_and_sigma(maxcorr, data_l2%freqmask_full, mean_maxcorr, sigma_maxcorr, rm_outliers)
     call get_mean_and_sigma(meancorr, data_l2%freqmask_full, mean_meancorr, sigma_meancorr, rm_outliers)
     call get_mean_and_sigma(meanabscorr, data_l2%freqmask_full, mean_meanabscorr, sigma_meanabscorr, rm_outliers)
     call get_mean_and_sigma(vars, data_l2%freqmask_full, mean_vars, sigma_vars, rm_outliers)
     
     median = 0.001d0!0.0055d0
     std_median = 0.0001d0 
     deallocate(subt)
     allocate(subt(nfreq-1))
   
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, nfreq
              if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
              median(k) = median(k) + 0.0001 * sign(1.d0, meanabscorr(k,j,i) - median(k))
           end do
           subt = (meanabscorr(2:nfreq,j,i) - meanabscorr(1:nfreq-1,j,i)) / sqrt(2.d0)
           var_0 = sum(merge(subt ** 2, 0.d0, ((data_l2%freqmask_full(2:nfreq,j,i) == 1) .and. (data_l2%freqmask_full(1:nfreq-1,j,i) == 1))))
           if (sum(merge(1.d0, 0.d0, ((data_l2%freqmask_full(2:nfreq,j,i) == 1) .and. (data_l2%freqmask_full(1:nfreq-1,j,i) == 1)))) == 0.d0) then
              if (verb) then
                 write(*,*) "OVERFLOW AVERTED", i, j, id
              end if
              var_0 = 0.d0
           else 
              var_0 = var_0 / sum(merge(1.d0, 0.d0, ((data_l2%freqmask_full(2:nfreq,j,i) == 1) .and. (data_l2%freqmask_full(1:nfreq-1,j,i) == 1))))
              std_median = std_median + 0.00002 * sign(1.d0, sqrt(var_0) - std_median)
           end if
        end do
     end do
 !    write(*,*) std_median, median
     ! Mask outlier frequencies
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, nfreq
              ! cut on neighbours
              if (k < nfreq - n_neighbor) then
                 if (all(var_cut * sigma_vars * neighbor_factor < abs(vars(k:k+n_neighbor,j,i) - mean_vars))) then
                    data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                    data_l2%freqmask_reason(k,j,i) = 5
                 else if (all(corr_cut * sigma_maxcorr * neighbor_factor < abs(maxcorr(k:k+n_neighbor,j,i) - mean_maxcorr))) then
                    data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                    data_l2%freqmask_reason(k,j,i) = 6
                 else if (all(mean_corr_cut * sigma_meancorr * neighbor_factor < abs(meancorr(k:k+n_neighbor,j,i) - mean_meancorr))) then
                    data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                    data_l2%freqmask_reason(k,j,i) = 7
                 else if (all(mean_abs_corr_cut * sigma_meanabscorr * neighbor_factor < abs(meanabscorr(k:k+n_neighbor,j,i) - mean_meanabscorr))) then
                    data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                    data_l2%freqmask_reason(k,j,i) = 8
                 ! else if (all(median_cut * std_median * neighbor_factor < meanabscorr(k:k+n_neighbor,j,i) - median(k))) then
                 !    data_l2%freqmask_full(k:k+n_neighbor,j,i) = 0.d0
                 !    data_l2%freqmask_reason(k,j,i) = 9
                 end if
              end if
              
              if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
              ! cut on single outliers
              if (var_cut * sigma_vars < abs(vars(k,j,i) - mean_vars)) then
                 data_l2%freqmask_full(k,j,i) = 0.d0
                 data_l2%freqmask_reason(k,j,i) = 10
              else if (corr_cut * sigma_maxcorr < abs(maxcorr(k,j,i) - mean_maxcorr)) then
                 data_l2%freqmask_full(k,j,i) = 0.d0
                 data_l2%freqmask_reason(k,j,i) = 11
              else if (mean_corr_cut * sigma_meancorr < abs(meancorr(k,j,i) - mean_meancorr)) then
                 data_l2%freqmask_full(k,j,i) = 0.d0
                 data_l2%freqmask_reason(k,j,i) = 12
              else if (mean_abs_corr_cut * sigma_meanabscorr < abs(meanabscorr(k,j,i) - mean_meanabscorr)) then
                 data_l2%freqmask_full(k,j,i) = 0.d0
                 data_l2%freqmask_reason(k,j,i) = 13
 !             else if (median_cut * std_median < meanabscorr(k,j,i) - median(k)) then
 !                data_l2%freqmask_full(k,j,i) = 0.d0
 !                data_l2%freqmask_reason(k,j,i) = 14
              end if                       
           end do
        end do
     end do
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, nfreq
              if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
              if (temp_freqmask(k,j,i) == 0.d0) then
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
     
     if (data_l2%use_freq_filter) then
        deallocate(corr_templates)
     else if (data_l2%polyorder == 1) then
        deallocate(corr_template)
     end if
     
     if (data_l2%n_pca_comp_feed > 0) then
        deallocate(corr_template_pca)
     end if 

     deallocate(vars)
     deallocate(means)
     deallocate(median)
     deallocate(corrs)
     deallocate(tod_subset)
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
     if (sum(mask) == 0) then
        mean = 0
        sigma = 1d20
     else
        mean = sum(data * mask) / sum(mask)
 !       write(*,*) "in gmas:", sum(data ** 2) / sum(mask), mean ** 2
        sigma  = sqrt(sum((data * mask) ** 2) / sum(mask) - mean ** 2)
 
        if (remove_outliers) then
           allocate(outlier_mask(size(data,1),size(data,2),size(data,3)))
           outlier_mask = mask * merge(1.d0,0.d0,abs(data - mean) <= 3.0 * sigma)
           mean = sum(data * outlier_mask) / sum(outlier_mask)
           sigma = sqrt(sum((data * outlier_mask) ** 2) / sum(outlier_mask) - mean ** 2)
           deallocate(outlier_mask)
        end if
     end if
   end subroutine get_mean_and_sigma
   
   subroutine interpolate_nans(data, id, verb)
     implicit none
     type(Lx_struct),  intent(inout) :: data
     integer(i4b),     intent(in)    :: id
     logical(lgt),     intent(in)    :: verb
 
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
        do i = 1, ndet
           do j = 1, nsb
              do l = 2, nsamp-1
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
     if (verb) then
        if (ntot > 0) write(*,*) '  Interpolated over NaN samples in ', id, ', n_tot = ', ntot
     end if
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
     data_l2%mean_tp = 0.d0
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
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
   
   subroutine remove_pca_components(data_l2_in, data_l2_out, pca_sig_rem)
     implicit none
     type(Lx_struct),            intent(in)    :: data_l2_in
     type(Lx_struct),            intent(inout) :: data_l2_out
     real(dp),                  intent(in)    :: pca_sig_rem
     integer(i4b) :: i, j, k, l, nsamp, nfreq, nsb, ndet, stat, iters
     real(dp)     :: std_tol, comp_std, amp_lim, dnu, radiometer
     
     nsamp       = size(data_l2_in%tod,1)
     nfreq       = size(data_l2_in%tod,2)
     nsb         = size(data_l2_in%tod,3)
     ndet        = size(data_l2_in%tod,4)
 
     ! Thresholds for removing PCA-components
     std_tol = pca_sig_rem / sqrt(real(nsamp))
     dnu = (data_l2_in%nu(2, 1, 1) - data_l2_in%nu(3, 1, 1)) * 1d9
     radiometer = 1 / sqrt(dnu * 1.0 / data_l2_in%samprate)
     amp_lim = std_tol * radiometer
 
     data_l2_out%n_pca_comp = data_l2_in%n_pca_comp
 
     if (data_l2_out%n_pca_comp == 0) return
 
     if (.not. allocated(data_l2_out%pca_comp)) allocate(data_l2_out%pca_comp(nsamp,data_l2_out%n_pca_comp))
     data_l2_out%pca_comp = data_l2_in%pca_comp
     if (.not. allocated(data_l2_out%pca_ampl)) allocate(data_l2_out%pca_ampl(nfreq,nsb,ndet,size(data_l2_in%pca_ampl,4)))
     data_l2_out%pca_ampl = data_l2_in%pca_ampl
     if (.not. allocated(data_l2_out%pca_eigv)) allocate(data_l2_out%pca_eigv(size(data_l2_in%pca_eigv,1)))
     data_l2_out%pca_eigv = data_l2_in%pca_eigv
     
     call free_lx_struct(data_l2_in)
     
     do l = 1, data_l2_out%n_pca_comp
        !write(*,*) 'removing component'
        comp_std = sqrt(sum(data_l2_out%pca_comp(:,l) ** 2) / nsamp - (sum(data_l2_out%pca_comp(:,l)) / nsamp) ** 2)
        
        do i = 1, ndet
           if (.not. is_alive(data_l2_out%pixels(i))) cycle
           do j = 1, nsb
              do k = 1, nfreq
                 if (data_l2_out%freqmask_full(k,j,i) == 0.d0) cycle
                 data_l2_out%tod(:,k,j,i) = data_l2_out%tod(:,k,j,i) - data_l2_out%pca_ampl(k,j,i,l) * data_l2_out%pca_comp(:,l)
              end do
           end do
        end do
     end do
     
   end subroutine remove_pca_components
   
   subroutine pca_filter_TOD(data_l2, n_pca_comp, pca_max_iter, pca_err_tol, pca_sig_rem, verb)
     implicit none
     type(Lx_struct),           intent(inout) :: data_l2
     integer(i4b),              intent(in)    :: n_pca_comp, pca_max_iter
     real(dp),                  intent(in)    :: pca_err_tol, pca_sig_rem
     logical(lgt),              intent(in)    :: verb
     integer(i4b) :: i, j, k, l, nsamp, nfreq, nsb, ndet, stat, iters
     real(dp)     :: eigenv, dotsum, amp, err, ssum 
     real(dp)     :: std_tol, comp_std, amp_lim, dnu, radiometer
     real(dp),     allocatable, dimension(:)   :: r, s, mys
     CHARACTER(LEN=128) :: number
     
 
     data_l2%n_pca_comp = n_pca_comp
     if (n_pca_comp == 0) return
 
     nsamp       = size(data_l2%tod,1)
     nfreq       = size(data_l2%tod,2)
     nsb         = size(data_l2%tod,3)
     ndet        = size(data_l2%tod,4)
 
     ! Thresholds for removing PCA-components
     std_tol = pca_sig_rem / sqrt(real(nsamp))
     dnu = (data_l2%nu(2, 1, 1) - data_l2%nu(3, 1, 1)) * 1d9
     radiometer = 1 / sqrt(dnu * 1.0 / data_l2%samprate)
     amp_lim = std_tol * radiometer
 !    write(*,*) dnu
     allocate(r(nsamp), s(nsamp))
             
     
     if(.not. allocated(data_l2%pca_ampl)) allocate(data_l2%pca_ampl(nfreq,nsb,ndet,n_pca_comp)) 
     if(.not. allocated(data_l2%pca_comp)) allocate(data_l2%pca_comp(nsamp,n_pca_comp))
     if(.not. allocated(data_l2%pca_eigv)) allocate(data_l2%pca_eigv(n_pca_comp))
     data_l2%pca_ampl = 0.d0
     data_l2%pca_comp = 0.d0
     data_l2%pca_eigv = 0.d0
     
     do l = 1, n_pca_comp 
        err = 1.d0
        r(:) = sum(sum(sum(data_l2%tod, 2), 2), 2) !sum of all freqs
        if (verb) then
           if (sum(r) == 0.d0) then
              write(*,*) "PCA initialized with zero vector"
           end if
           if (sum(r) .ne. sum(r)) then
              write(*,*) "NaN in initial PCA vector"
           end if
        end if       
        iters = 0
        do while ((err > pca_err_tol) .and. (iters < pca_max_iter))
           s = 0.d0
           !$OMP PARALLEL PRIVATE(k,i,j,dotsum,mys)
           allocate(mys(nsamp))
           mys = 0.d0
           !$OMP DO SCHEDULE(guided)
           do k = 1, nfreq
              do i = 1, ndet
                 if (.not. is_alive(data_l2%pixels(i))) cycle
                 do j = 1, nsb
                    if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle   
                    dotsum = sum(data_l2%tod(:,k,j,i) * r(:))
                    !write(*,*) "yo", k, i, j
           !         write(*,*) sum(data_l2%tod(:,k,j,i))
                    !dotsum = sdot(nfreq,data_l2%tod(:,k,j,i), 1, r(:), 1)
                    mys(:) = mys(:) + dotsum * data_l2%tod(:,k,j,i)
                 end do
              end do
           end do
           !$OMP END DO
           !$OMP CRITICAL
           s(:) = s(:) + mys(:)
           !$OMP END CRITICAL
           deallocate(mys)
           !$OMP END PARALLEL
           eigenv = sum(s(:) * r(:))
           err = sqrt(sum((eigenv * r(:) - s(:)) ** 2))
           ! write(*,*) sum(s(:) ** 2)
           ! write(*,*) s(1), s(170)
           ! write(*,*) iters, l
           ssum = sqrt(sum(s(:) ** 2))
           if (ssum == 0.d0) then
              write(*,*) "Weird stuff happening in PCA-filter"
              r(:) = 1.d0 / sqrt(1.d0 * nsamp)
           else
              r(:) = s(:)/ssum
           end if
           iters = iters + 1
        end do
        data_l2%pca_eigv(l) = eigenv
        data_l2%pca_comp(:,l) = r(:)
        !means(k,j,i) = sum(data_l2%tod(:,k,j,i)) / nsamp
        !vars(k,j,i) = sum(data_l2%tod(:,k,j,i) ** 2) / nsamp - means(k,j,i) ** 2
 
        comp_std = sqrt(abs(sum(r ** 2) / nsamp - (sum(r) / nsamp) ** 2))
        !$OMP PARALLEL PRIVATE(k,i,j,dotsum,mys)
        !$OMP DO SCHEDULE(guided)
        do k = 1, nfreq
           do i = 1, ndet
              if (.not. is_alive(data_l2%pixels(i))) cycle
              do j = 1, nsb
                 if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
                 data_l2%pca_ampl(k,j,i,l) = sum(r(:) * data_l2%tod(:,k,j,i))
                 !write(*,*) data_l2%pca_ampl(k,j,i,l), amp_lim / comp_std, l
                 !if (abs(data_l2%pca_ampl(k,j,i,l)) > amp_lim / comp_std) then
                 data_l2%tod(:,k,j,i) = data_l2%tod(:,k,j,i) - data_l2%pca_ampl(k,j,i,l) * r(:)
                 !end if
              end do
           end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        if (.not. (any(sum(abs(data_l2%pca_ampl(:, :, :, l)), 1) / nfreq > amp_lim / comp_std))) EXIT
     end do
     deallocate(r, s)
     
   end subroutine pca_filter_TOD
 
 
   
   subroutine pca_filter_feed_TOD(data_l2, n_pca_comp, pca_max_iter, pca_err_tol, pca_sig_rem, verb)

    implicit none
    type(Lx_struct),           intent(inout) :: data_l2
    integer(i4b),              intent(in)    :: n_pca_comp, pca_max_iter
    real(dp),                  intent(in)    :: pca_err_tol, pca_sig_rem
    logical(lgt),              intent(in)    :: verb
    integer(i4b) :: i, j, k, l, nsamp, nfreq, nsb, ndet, stat, iters, numfreq_out, delta_nu
    real(dp)     :: eigenv, dotsum, amp, err, ssum, w, weight, Tsys_downsamp
    real(dp)     :: std_tol, comp_std, amp_lim, dnu, radiometer
    real(dp),     allocatable, dimension(:)   :: r, s, mys
    real(sp),     allocatable, dimension(:,:,:,:)   :: tod_downsamp   ! (time, freq, sideband, detector) used for temporary downsampled tod
    real(sp),     allocatable, dimension(:,:,:,:)   :: pca_ampl_feed_lr   ! (freq, sideband, detector, ncomp) used for temporary downsampled tod
    integer(i4b), allocatable, dimension(:,:,:)   :: freqmask_lr   ! (freq, sideband, detector, ncomp) lowres frequency mask
 
    CHARACTER(LEN=128) :: number
    
    !=======================================================================
    !====                     Feed PCA filter                          =====
    !=======================================================================
    
 
 
    data_l2%n_pca_comp_feed = n_pca_comp
    if (n_pca_comp == 0) return
      
    write(*,*) "Feed-PCA filtering TOD"

    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
 
 
    numfreq_out = data_l2%n_freq_downsamp
    
    if (numfreq_out /= 0) then
       allocate(tod_downsamp(nsamp, numfreq_out, nsb, ndet))
       allocate(freqmask_lr(numfreq_out, nsb, ndet))

       tod_downsamp = 0.d0
       freqmask_lr = 1

       delta_nu                    = size(data_l2%nu, 1) / numfreq_out
       call assert(delta_nu >= 1, 'Cannot ask for more frequencies than in input files')
       
       
       do i = 1, nsamp
          do j = 1, nsb
             do k = 1, numfreq_out
                do l = 1, ndet           ! Time-ordered data
                   if (.not. is_alive(data_l2%pixels(l))) cycle
                   tod_downsamp(i,k,j,l) = 0.d0
                   weight                = 0.d0
                   if (sum(data_l2%freqmask_full((k-1)*delta_nu+1:k*delta_nu, j, l)) == 0) then
                     freqmask_lr(k, j, l) = 0
                   end if

                   do n = (k-1)*delta_nu+1, k*delta_nu
                      if (data_l2%freqmask_full(n,j,l) == 0) cycle
                      !write(*,*) "TSYS = ", data_l2%Tsys(n,j,l), n, j, l 
                      if (data_l2%Tsys(n,j,l) .ne. data_l2%Tsys(n,j,l)) then
                         write(*,*) "Tsys: I have become NaN, destroyer of codes:", i, j, k, l, data_l2%Tsys(n,j,l)
                      end if
                     
                      if (data_l2%Tsys(n,j,l) <= 0) then 
                        write(*,*) "Tsys less or equal to zero. Tsys = " 
                         w = 0.d0
                      else
                         w  = 1.d0 / (data_l2%Tsys(n,j,l)) * data_l2%freqmask_full(n,j,l)
                      end if
                      weight = weight + w * w 
                      tod_downsamp(i,k,j,l) = tod_downsamp(i,k,j,l) + w * data_l2%tod(i,n,j,l) 
                   end do
                   
                   if (weight .ne. weight) then
                      write(*,*) "weight: I have become NaN, destroyer of codes:", i, j, k, l, weight
                   end if
 
                   if (weight > 0.d0) then
                     Tsys_downsamp = sqrt(delta_nu / weight)  ! (delta_nu) = delta_nu_lowres / delta_nu_highres
                     !write(*,*) "Updating tod_downsamp", Tsys_downsamp
                      tod_downsamp(i,k,j,l) = tod_downsamp(i,k,j,l) / (weight * Tsys_downsamp) 
                   else
                      !write(*,*) "Zero weights!!!"
                      tod_downsamp(i,k,j,l) = 0.d0
                   end if
 
                end do
             end do
          end do
       end do
    
      !  do i = 1, nsamp
      !     do j = 1, nsb
      !        do k = 1, numfreq_out
      !           do l = 1, ndet           ! Time-ordered data
      !              if (.not. is_alive(data_l2%pixels(l))) cycle
      !              tod_downsamp(i,k,j,l) = 0.d0
      !              weight = 0.d0
      !              do n = (k-1)*delta_nu+1, k*delta_nu
      !                 if (data_l2%freqmask_full(n,j,l) == 0) cycle
      !                 tod_downsamp(i,k,j,l) = tod_downsamp(i,k,j,l) + data_l2%tod(i,n,j,l) 
      !                 weight = weight + 1
      !              end do
      !              if (weight > 0) then
      !                tod_downsamp(i,k,j,l) = tod_downsamp(i,k,j,l) / weight 
      !              else
      !                tod_downsamp(i,k,j,l) = 0.d0 
      !              end if
                
      !           end do
      !        end do
      !     end do
      !  end do 
    
 
       if(.not. allocated(pca_ampl_feed_lr)) allocate(pca_ampl_feed_lr(numfreq_out, nsb, ndet, n_pca_comp)) 
       pca_ampl_feed_lr = 0.d0

    end if
 
    ! Thresholds for removing PCA-components
    std_tol = pca_sig_rem / sqrt(real(nsamp))
    dnu = (data_l2%nu(2, 1, 1) - data_l2%nu(3, 1, 1)) * 1d9
    radiometer = 1 / sqrt(dnu * 1.0 / data_l2%samprate)
    amp_lim = std_tol * radiometer
 !    write(*,*) dnu
    allocate(r(nsamp), s(nsamp))
    
    if(.not. allocated(data_l2%pca_ampl_feed)) allocate(data_l2%pca_ampl_feed(nfreq, nsb, ndet, n_pca_comp)) 
    if(.not. allocated(data_l2%pca_comp_feed)) allocate(data_l2%pca_comp_feed(ndet, nsamp, n_pca_comp))
    if(.not. allocated(data_l2%pca_eigv_feed)) allocate(data_l2%pca_eigv_feed(ndet, n_pca_comp))
    data_l2%pca_ampl_feed = 0.d0
    data_l2%pca_comp_feed = 0.d0
    data_l2%pca_eigv_feed = 0.d0
    
    do i = 1, ndet
       if (.not. is_alive(data_l2%pixels(i))) cycle
       !if (sum(sum(data_l2%freqmask_full(:, :, i), 1), 1) <= n_pca_comp) cycle
       if (sum(sum(data_l2%freqmask_full(:, :, i), 1), 1) <= 0.1 * nsb * nfreq) cycle
       do l = 1, n_pca_comp 
          err = 1.d0
          if (numfreq_out == 0) then
             r(:) = sum(sum(data_l2%tod(:, :, :, i), 2), 2) !sum of all freqs
          else
             r(:) = sum(sum(tod_downsamp(:, :, :, i), 2), 2) !sum of all freqs
          end if

          if (verb) then
             if (sum(r) == 0.d0) then
                write(*,*) "Feed-PCA initialized with zero vector"
                !stop
                !cycle
             end if
             if (sum(r) .ne. sum(r)) then
                write(*,*) "NaN in initial Feed-PCA vector"
                !cycle
             end if
          end if 
 
          iters = 0
          do while ((err > pca_err_tol) .and. (iters < pca_max_iter))
             s = 0.d0
             !$OMP PARALLEL PRIVATE(k, j, dotsum, mys)
             allocate(mys(nsamp))
             mys = 0.d0
 
             if (numfreq_out == 0) then
                !$OMP DO SCHEDULE(guided)
                do k = 1, nfreq
                !do k = 1, numfreq_out
                   do j = 1, nsb
                      if (data_l2%freqmask_full(k, j, i) == 0) cycle   
                 
                      dotsum = sum(data_l2%tod(:, k, j, i) * r(:))
                      !dotsum = sum(tod_downsamp(:, k, j, i) * r(:))
 
                      mys(:) = mys(:) + dotsum * data_l2%tod(:, k, j, i)
                      !mys(:) = mys(:) + dotsum * tod_downsamp(:,k,j,i)
                   end do
                end do
                !$OMP END DO
             else
                !$OMP DO SCHEDULE(guided)
                !do k = 1, nfreq
                do k = 1, numfreq_out
                   do j = 1, nsb
                      if (freqmask_lr(k, j, i) == 0) cycle
                      !if (data_l2%freqmask(k, j, i) == 0.d0) cycle   
                      
                      !dotsum = sum(data_l2%tod(:, k, j, i) * r(:))
                      dotsum = sum(tod_downsamp(:, k, j, i) * r(:))
 
                      !mys(:) = mys(:) + dotsum * data_l2%tod(:,k,j,i)
                      mys(:) = mys(:) + dotsum * tod_downsamp(:, k, j, i)
                   end do
                end do
                !$OMP END DO
             end if
             !$OMP CRITICAL
             s(:) = s(:) + mys(:)
             !$OMP END CRITICAL
             deallocate(mys)
             !$OMP END PARALLEL
             eigenv = sum(s(:) * r(:))
             err = sqrt(sum((eigenv * r(:) - s(:)) ** 2))
             ! write(*,*) sum(s(:) ** 2)
             ! write(*,*) s(1), s(170)
             ! write(*,*) iters, l
             ssum = sqrt(sum(s(:) ** 2))
             
             if (ssum .ne. ssum) then
 
                write(*,*) "NaN in ssum"
             end if
             do k = 1, nsamp
                if (r(k) .ne. r(k)) then
 
                   write(*,*) "NaN in r(k)"
                end if
             end do
 
             if (ssum == 0.d0) then
                write(*,*) "Weird stuff happening in Feed-PCA-filter", i
                r(:) = 1.d0 / sqrt(1.d0 * nsamp)
             else
                r(:) = s(:)/ssum
             end if
             iters = iters + 1
          end do
 
          data_l2%pca_eigv_feed(i, l) = eigenv
          data_l2%pca_comp_feed(i, :, l) = r(:)
          !means(k,j,i) = sum(data_l2%tod(:,k,j,i)) / nsamp
          !vars(k,j,i) = sum(data_l2%tod(:,k,j,i) ** 2) / nsamp - means(k,j,i) ** 2
          comp_std = sqrt(abs(sum(r ** 2) / nsamp - (sum(r) / nsamp) ** 2))
          !$OMP PARALLEL PRIVATE(k,j,dotsum,mys)
          !$OMP DO SCHEDULE(guided)
          do j = 1, nsb
             do k = 1, nfreq
                if (data_l2%freqmask_full(k, j, i) == 0.d0) cycle
                data_l2%pca_ampl_feed(k, j, i, l) = sum(r(:) * data_l2%tod(:, k, j, i))
                !if (abs(data_l2%pca_ampl_feed(k, j, i, l)) > amp_lim / comp_std) then
                data_l2%tod(:, k, j, i) = data_l2%tod(:, k, j, i) - data_l2%pca_ampl_feed(k, j, i, l) * r(:)
                !end if
             end do
 
             if (numfreq_out /= 0) then
                do k = 1, numfreq_out
                   if (freqmask_lr(k, j, i) == 0) cycle
                   pca_ampl_feed_lr(k, j, i, l) = sum(r(:) * tod_downsamp(:, k, j, i))
                   tod_downsamp(:, k, j, i) = tod_downsamp(:, k, j, i) - pca_ampl_feed_lr(k, j, i, l) * r(:)
                end do
             end if
          end do
          !$OMP END DO
          !$OMP END PARALLEL
          !if (.not. (any(sum(abs(data_l2%pca_ampl_feed(:, :, i, l)), 1) / nfreq > amp_lim / comp_std))) EXIT
       end do
    end do
 
    deallocate(r, s)
    if (allocated(tod_downsamp)) deallocate(tod_downsamp)
    if (allocated(pca_ampl_feed_lr)) deallocate(pca_ampl_feed_lr)
    
   end subroutine pca_filter_feed_TOD
 
 
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
 
     allocate(T(nfreq,0:p), A(0:p,0:p))
     allocate(data_l2%tod_poly(nsamp,0:p,nsb,ndet))
 
     ! Precompute polynomial basis
     do k = 1, nfreq
        mu = max(min(2.d0*real(k-1,dp)/real(nfreq-1,dp)-1.d0,1.d0),-1.d0)
        call get_legendre_polynomials(mu,T(k,:))
     end do
 
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
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
           allocate(x(0:p))
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
           deallocate(x)
           !$OMP END PARALLEL
        end do
     end do
     deallocate(T, A)
     
 !    open(22, file="tod_after_poly_2_3.unf", form="unformatted") ! Adjusted open statement
 !    write(22) data_l2%tod(:,:,3,1)
 !    close(22)
 !    open(22, file="tod_after_poly_15_2.unf", form="unformatted") ! Adjusted open statement
 !    write(22) data_l2%tod(:,:,1,14)
 !    close(22)
 
   end subroutine polyfilter_TOD
 
   subroutine frequency_filter_TOD(data_l2, fknee_W, alpha_W, scan)
    implicit none
    type(Lx_struct),                            intent(inout) :: data_l2
    type(comap_subscan),                        intent(in)    :: scan
    real(dp),                                   intent(in)    :: fknee_W, alpha_W
    integer(i4b) :: i, j, k, l, n, nsamp, nfreq, nsb, ndet, stat, err, nomp, feed, sb, n_order
    real(dp) :: detinv, samprate, nu, mu, t3, t4, t5, t6, sigma0, dnu
    real(sp), allocatable, dimension(:,:) :: F, FTZ, a, P, PTP, PTP_inv, m, I_mat, Z, z_scalar, FTZY
    real(sp), allocatable, dimension(:) :: Cf, sigma0_prior, fknee_prior, alpha_prior
    real(sp),     allocatable, dimension(:) :: dt
    complex(spc), allocatable, dimension(:) :: dv
    integer*8    :: plan_fwd, plan_back
    type(hdf_file)   :: prior_file
 
    call wall_time(t5)
    
 
 
    nsamp       = size(data_l2%tod,1)
    nfreq       = size(data_l2%tod,2)
    nsb         = size(data_l2%tod,3)
    ndet        = size(data_l2%tod,4)
    n = nsamp + 1
     
    samprate = data_l2%samprate
 
    n_order = 1  
  
    if(.not. allocated(data_l2%T_cont)) allocate(data_l2%T_cont(nsamp,nsb,ndet,n_order+1)) 
    if(.not. allocated(data_l2%dg)) allocate(data_l2%dg(nsamp,nsb,ndet)) 
 
    data_l2%T_cont = 0.d0
    data_l2%dg = 0.d0
 
    allocate(P(nfreq, 2))
    allocate(PTP(2, 2))
    allocate(PTP_inv(2, 2))
    allocate(m(2, nsamp))
    allocate(F(nfreq, 1))
    allocate(FTZ(1, nfreq))
    allocate(a(1, nsamp))
    allocate(FTZY(1, nsamp))
    allocate(Z(nfreq, nfreq))
    allocate(I_mat(nfreq, nfreq))
    allocate(z_scalar(1,1))
    allocate(Cf(0:n-1))
    allocate(sigma0_prior(20))
    allocate(fknee_prior(20))
    allocate(alpha_prior(20))
 
    ! Load in the feed-specific 1/f prior parameters on the gain fluctuations.
    call open_hdf_file("/mn/stornext/d22/cmbco/comap/protodir/auxiliary/Cf_prior_data.hdf5", prior_file, "r")  ! get from parameter file 
    call read_hdf(prior_file, "sigma0_prior", sigma0_prior)
    call read_hdf(prior_file, "fknee_prior", fknee_prior)
    call read_hdf(prior_file, "alpha_prior", alpha_prior)
    call close_hdf_file(prior_file)
 
    ! Wiener filter normalization parameters.
 !   fknee_W = 0.01
 !   alpha_W = -4.0
    dnu = (data_l2%nu(2, 1, 1) - data_l2%nu(3, 1, 1)) * 1d9
    !radiometer = 1.0d0 / sqrt(dnu * 1.d0 / data_l2%samprate)
    
    sigma0 = 1.0d0 / sqrt(dnu * 1.d0 / data_l2%samprate)  !0.005  ! From radiometer eq
 
    ! Create identity matrix
    do i = 1, nfreq
       do j = 1, nfreq
             I_mat(j,i) = 0
       end do
       I_mat(i,i) = 1
    end do
 
    ! open(1, file = 'ftr_Cf.dat', status = 'new')
    ! do j=1,n-1
    !    write(1,*) Cf(j)
    ! end do
 
    ! Doing some setup stuff for the Fourier module. Don't really know what this does.
    nomp = 1
    call sfftw_init_threads(err)
    call sfftw_plan_with_nthreads(nomp)
    allocate(dt(2*nsamp), dv(0:n-1))
    call sfftw_plan_dft_r2c_1d(plan_fwd,  2*nsamp, dt, dv, fftw_estimate + fftw_unaligned)
    call sfftw_plan_dft_c2r_1d(plan_back, 2*nsamp, dv, dt, fftw_estimate + fftw_unaligned)
    deallocate(dt, dv)
 
    allocate(dt(2*nsamp), dv(0:n-1))
    do feed = 1, ndet
       if (.not. is_alive(data_l2%pixels(feed))) cycle
       do sb = 1, nsb
          if (all(data_l2%freqmask_full(:,sb,feed) == 0.d0)) cycle
          ! Create the 1/f prior "Cf" for the gain flucuations. Divide it by the Wiener filter PS, to compensate for the normalization of the data.
          do i = 1, n-1 
             nu = ind2freq(i+1, samprate, n)
             Cf(i) = sigma0_prior(data_l2%pixels(feed))**2*(nu/fknee_prior(data_l2%pixels(feed)))**alpha_prior(data_l2%pixels(feed))
             Cf(i) = Cf(i)/(1.d0 + (nu/fknee_W)**alpha_W)
          end do
          Cf(0) = 1.d0      ! Should this be zero?
 
          ! Create P and F matrices.
          do i = 1, nfreq
             if (data_l2%freqmask_full(i,sb,feed) == 0.d0) then
                P(i,1) = 0
                P(i,2) = 0
                F(i,1) = 0
             else
                P(i,1) = 1.0/data_l2%Tsys(i,sb,feed)
                P(i,2) = ((2.0d0*i)/nfreq - 1.0d0)/data_l2%Tsys(i,sb,feed)
                F(i,1) = 1
             end if
          end do
 
          ! Solve (PT*P)^-1
          PTP = MATMUL(TRANSPOSE(P), P)
          detinv = 1/(PTP(1,1)*PTP(2,2) - PTP(1,2)*PTP(2,1))
          PTP_inv(1,1) = +detinv * PTP(2,2)
          PTP_inv(2,1) = -detinv * PTP(2,1)
          PTP_inv(1,2) = -detinv * PTP(1,2)
          PTP_inv(2,2) = +detinv * PTP(1,1)
 
          ! Solve for the Z-matrix and the scalar quantity z.
          !!!!!! This is how this section used to look, which makes no sense to me, since Z is not initialized
          ! FTZ = MATMUL(TRANSPOSE(F), Z)
          ! Z = I_mat - MATMUL(MATMUL(P, PTP_inv), TRANSPOSE(P))
          ! z_scalar = MATMUL(FTZ, F)
          !!!!!! This is the new version that makes sense to me
          Z = I_mat - MATMUL(MATMUL(P, PTP_inv), TRANSPOSE(P))
          FTZ = MATMUL(TRANSPOSE(F), Z)
          z_scalar = MATMUL(FTZ, F)
 
 
          ! Solve F^T*Z*y, which is the RHS of the linear equation for a.
          FTZY = 0.d0
          !$OMP PARALLEL PRIVATE(i, j)
          !$OMP DO SCHEDULE(guided)
          do i=1,nsamp
             do j=1,nfreq
                if (data_l2%freqmask_full(j,sb,feed) == 0.d0) cycle
                FTZY(1,i) = FTZY(1,i) + FTZ(1,j)*data_l2%tod(i,j,sb,feed)
             end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
 
          ! Solve the fourier transform of FTZY, which we name dt.
          dt(1:nsamp)            = FTZY(1,:)
          dt(2*nsamp:nsamp+1:-1) = dt(1:nsamp)
          call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
 
          ! Put a prior on a in the fourier domain.
          do i = 0, n-1
             dv(i) = dv(i)/(z_scalar(1,1) + (sigma0**2)/Cf(i))
          end do
          
          ! Fourier transform it back to get a (gain fluctuations).
          call sfftw_execute_dft_c2r(plan_back, dv, dt)
          dt = dt / (2*nsamp)
          a(1,:) = dt
 
          ! Solve for m (containing the temperature fluctuations).
          m = 0.d0
          !$OMP PARALLEL PRIVATE(i, j)
          !$OMP DO SCHEDULE(guided)
          do i=1,nsamp
             do j=1,nfreq
                m(1,i) = m(1,i) + (data_l2%tod(i,j,sb,feed) - F(j,1)*a(1,i))*P(j,1)
                m(2,i) = m(2,i) + (data_l2%tod(i,j,sb,feed) - F(j,1)*a(1,i))*P(j,2)
             end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
          m = MATMUL(PTP_inv, m)
 
          if ((sum(a(1,:)) .ne. sum(a(1,:))) .or. (sum(m(1,:)) .ne. sum(m(1,:))) .and. (sum(m(2,:)) .ne. sum(m(2,:)))) then
             write(*,*) "Problem 2 in freq filter: ", a(1,1), m(1,1), m(2,1), feed, sb, scan%id
             data_l2%freqmask_full(:,sb,feed) = 0.d0
             cycle
          end if
 
          ! Subtract the gain and temperature fluctuations from the TOD.
          !$OMP PARALLEL PRIVATE(i, j)
          !$OMP DO SCHEDULE(guided)
          do i=1,nsamp
             data_l2%tod(i,:,sb,feed) = data_l2%tod(i,:,sb,feed) - F(:,1)*a(1,i) - P(:,1)*m(1,i) - P(:,2)*m(2,i)
             data_l2%dg(i,sb,feed) = a(1,i)
             data_l2%T_cont(i,sb,feed,1) = m(1,i)
             data_l2%T_cont(i,sb,feed,2) = m(2,i)
             ! if ((a(1,i) == a(1,i)) .and. (m(1,i) == m(1,i)) .and. (m(2,i) == m(2,i))) then
             !    ! do j=1,nfreq
             !    !    data_l2%tod(i,j,sb,feed) = data_l2%tod(i,j,sb,feed) - F(j,1)*a(1,i) - P(j,1)*m(1,i) - P(j,2)*m(2,i)
             !    ! end do
             !    data_l2%tod(i,:,sb,feed) = data_l2%tod(i,:,sb,feed) - F(:,1)*a(1,i) - P(:,1)*m(1,i) - P(:,2)*m(2,i)
             !    data_l2%dg(i,sb,feed) = a(1,i)
             !    data_l2%T_cont(i,sb,feed,1) = m(1,i)
             !    data_l2%T_cont(i,sb,feed,2) = m(2,i)
             ! else
             !    write(*,*) "Problem in freq filter: ", a(1,i), m(1,i), m(2,i), scan%id
             !    data_l2%freqmask_full(:,sb,feed) = 0.d0
             !    exit
             ! end if
          end do
          !$OMP END DO
          !$OMP END PARALLEL
 
       end do
    end do
    deallocate(dt, dv)
    deallocate(F, FTZ, a, P, PTP, PTP_inv, m, I_mat, Z, z_scalar, FTZY)
    deallocate(sigma0_prior, fknee_prior, alpha_prior, Cf)
    call wall_time(t6)
    write(*,*) "Time spent on frequency filter: ", t6 - t5
 
 end subroutine frequency_filter_TOD
 
 
   subroutine freq_filter_TOD(data_l2)
     implicit none
     type(Lx_struct),                            intent(inout) :: data_l2
     
     integer*8    :: plan_fwd, plan_back
     
     integer(i4b) :: i, j, k, l, m, n, nsamp, nfreq, nsb, ndet, p, stat, nomp, err
     real(dp)     :: samprate, nu, mu, freq_step, nu_knee, alpha, mean_val
     real(sp),     allocatable, dimension(:)   :: dt
     complex(spc), allocatable, dimension(:)   :: dv
     real(dp),     allocatable, dimension(:,:) :: T, A
 
     nsamp       = size(data_l2%tod,1)
     nfreq       = size(data_l2%tod,2)
     nsb         = size(data_l2%tod,3)
     ndet        = size(data_l2%tod,4)
     p = 1
     
     allocate(data_l2%tod_poly(nsamp,0:p,nsb,ndet))
     data_l2%tod_poly = 0.d0
        
     
     n = nfreq + 1
 
     ! Set up OpenMP environment and FFTW plans
     nomp = 1
     call sfftw_init_threads(err)
     call sfftw_plan_with_nthreads(nomp)
 
     allocate(dt(2*nfreq), dv(0:n-1))
     call sfftw_plan_dft_r2c_1d(plan_fwd,  2*nfreq, dt, dv, fftw_estimate + fftw_unaligned)
     call sfftw_plan_dft_c2r_1d(plan_back, 2*nfreq, dv, dt, fftw_estimate + fftw_unaligned)
     deallocate(dt, dv)
     freq_step = 31.25d6 / 16.d0 ! Hz
     nu_knee   = 1 / 4.d9    ! 1 / Hz
     alpha     = -4.d0
     do i = 1, ndet
        if (.not. is_alive(data_l2%pixels(i))) cycle
        do j = 1, nsb
           !!$OMP PARALLEL PRIVATE(k,l,m,mean_val,dt,dv,nu)
           allocate(dt(2*nfreq), dv(0:n-1))
           !!$OMP DO SCHEDULE(guided)
           
           do l = 1, nsamp
              mean_val = 0.d0
              do k = 1, nfreq
                 if (data_l2%freqmask_full(k,j,i) == 0.d0) then 
                    dt(k) = mean_val
                    cycle
                 end if
                 dt(k) = data_l2%tod(l,k,j,i)
              end do
              dt(2*nfreq:nfreq+1:-1) = dt(1:nfreq)
              call sfftw_execute_dft_r2c(plan_fwd, dt, dv)
              ! Apply highpass filter
              do m = 1, n-1
                 nu = ind2freq(m+1, 1.d0 / freq_step, n)
                 !write(*,*) m, nu, nu_knee
                 dv(m) = dv(m) * 1.d0/(1.d0 + (nu/nu_knee)**alpha)
              end do
              dv(1) = 0.d0
              !call fft(dt, dv, -1)
              call sfftw_execute_dft_c2r(plan_back, dv, dt)
              dt                       = dt / (2*nfreq)
              do k = 1, nfreq
                 if (data_l2%freqmask_full(k,j,i) == 0.d0) cycle
                 data_l2%tod(l,k,j,i) = dt(k)
              end do
              
           end do
           !!$OMP END DO
           deallocate(dt, dv)
           !!$OMP END PARALLEL
        end do
     end do
     call sfftw_destroy_plan(plan_fwd)
     call sfftw_destroy_plan(plan_back)
 
   end subroutine freq_filter_TOD
 
 
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
     call free_lx_struct(data_l2)
     ! Allocate full-resolution L2 structure
     allocate(data_l2%time(nsamp_tot), stat=err)
     allocate(data_l2%nu(nfreq,nsb,ndet))
     allocate(data_l2%tod(nsamp_tot, nfreq, nsb, ndet))
     allocate(data_l2%point_tel(3,nsamp_tot,ndet))
     allocate(data_l2%point_cel(3,nsamp_tot,ndet))
     allocate(data_l2%pixels(ndet))
     allocate(data_l2%pix2ind(size(data_l1%pix2ind,1)))
     allocate(data_l2%tod_mean(nfreq, nsb, ndet))
     allocate(data_l2%sb_mean(nsamp_tot, nsb, ndet))
     allocate(data_l2%freqmask(nfreq,nsb,ndet))
     allocate(data_l2%freqmask_full(size(data_l1%nu,1,1),nsb,ndet))
     allocate(data_l2%freqmask_reason(size(data_l1%nu,1,1),nsb,ndet))
     !allocate(data_l2%flag(nsamp_tot))
     allocate(data_l2%Phot(3,2,nfreq,nsb,ndet))
     allocate(data_l2%Thot(3,2))
     allocate(data_l2%time_hot(3,2))
     allocate(data_l2%n_nan(nfreq,nsb,ndet))
     
     ! Merge L1 data
     data_l2%decimation_time = 1
     data_l2%decimation_nu   = 1
     data_l2%samprate        = data_l1%samprate
     data_l2%nu              = data_l1%nu
     data_l2%pixels          = data_l1%pixels
     data_l2%pix2ind         = data_l1%pix2ind
     data_l2%point_tel       = data_l1%point_tel(:,ind(1):ind(2),:)
     data_l2%point_cel       = data_l1%point_cel(:,ind(1):ind(2),:)
     data_l2%time            = data_l1%time(ind(1):ind(2))
     data_l2%freqmask        = data_l1%freqmask
     data_l2%freqmask_full   = data_l1%freqmask_full
     data_l2%freqmask_reason = data_l1%freqmask_reason
     data_l2%cal_method      = data_l1%cal_method
     data_l2%n_cal           = data_l1%n_cal    
     data_l2%Phot            = data_l1%Phot
     data_l2%Thot            = data_l1%Thot
     data_l2%time_hot        = data_l1%time_hot
     data_l2%n_nan           = data_l1%n_nan
 
     do j = 1, ndet
        if (.not. is_alive(data_l2%pixels(j))) cycle
        do m = 1, nsb
           data_l2%sb_mean(:,m,j) = data_l1%sb_mean(ind(1):ind(2),m,j)
           do n = 1, nfreq
              if (data_l2%freqmask_full(n,m,j) == 0) cycle
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
 
     integer(i4b) :: i, j, k, l, m, n, nsamp_in, nsamp_out, ndet, dt, dnu, nsb, numfreq_in, nw
     real(dp)     :: w, weight, var
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
 
     if (allocated(data_out%time)) call free_lx_struct(data_out)
 
  
     allocate(data_out%time(nsamp_out))
     allocate(data_out%nu(numfreq_out,nsb,ndet))
     allocate(data_out%tod(nsamp_out, numfreq_out, nsb, ndet))
     allocate(data_out%tod_mean(numfreq_in, nsb, ndet))
     allocate(data_out%sb_mean(size(data_in%time), nsb, ndet))
     allocate(data_out%point_tel(3,nsamp_out,ndet))
     allocate(data_out%point_cel(3,nsamp_out,ndet))
     allocate(data_out%flag(nsamp_out))
     allocate(data_out%freqmask(numfreq_out,nsb,ndet))
     allocate(data_out%chi2(numfreq_out,nsb,ndet))
     allocate(data_out%Tsys(numfreq_in,nsb,ndet))
     allocate(data_out%Phot(3,2,numfreq_in,nsb,ndet))
     allocate(data_out%Thot(3,2))
     allocate(data_out%Time_hot(3,2))
     allocate(data_out%Tsys_lowres(numfreq_out,nsb,ndet))
     allocate(data_out%freqmask_full(size(data_in%nu,1,1),nsb,ndet))
     allocate(data_out%freqmask_reason(size(data_in%nu,1,1),nsb,ndet))
     allocate(data_out%mean_tp(size(data_in%nu,1),nsb,ndet))
     allocate(data_out%var_fullres(size(data_in%nu,1),nsb,ndet))
     allocate(data_out%n_nan(size(data_in%nu,1),nsb,ndet))
     data_out%mask_outliers = data_in%mask_outliers
     if (data_in%mask_outliers == 1) then
        allocate(data_out%acceptrate(nsb,ndet))
        allocate(data_out%diagnostics(size(data_in%nu,1),nsb,ndet,size(data_in%diagnostics,4)))
        allocate(data_out%cut_params(size(data_in%cut_params,1),size(data_in%cut_params,2)))
        allocate(data_out%AB_mask(size(data_in%nu,1,1),nsb,20))
        allocate(data_out%leak_mask(size(data_in%nu,1,1),nsb,20))
        data_out%acceptrate    = data_in%acceptrate
        data_out%diagnostics   = data_in%diagnostics
        data_out%cut_params    = data_in%cut_params
        data_out%AB_mask       = data_in%AB_mask
        data_out%leak_mask     = data_in%leak_mask
     end if
     allocate(data_out%pixels(ndet))
     allocate(data_out%pix2ind(size(data_in%pix2ind,1)))
     if (allocated(data_in%mean_tp)) then
        data_out%mean_tp = data_in%mean_tp
     else
        data_out%mean_tp = 0.d0
     end if
     if(allocated(data_in%spike_data))   then
        allocate(data_out%spike_data(size(data_in%spike_data,1),&
             &size(data_in%spike_data,2),size(data_in%spike_data,3),&
             &size(data_in%spike_data,4),size(data_in%spike_data,5)))
        data_out%spike_data = data_in%spike_data
     end if
 
     if(allocated(data_in%el_az_stats))   then
        allocate(data_out%el_az_stats(size(data_in%el_az_stats,1),&
             &size(data_in%el_az_stats,2),size(data_in%el_az_stats,3),&
             &size(data_in%el_az_stats,4),size(data_in%el_az_stats,5)))
        data_out%el_az_stats = data_in%el_az_stats
     end if
     
     data_out%freqmask      = data_in%freqmask
     data_out%freqmask_full = data_in%freqmask_full
     data_out%freqmask_reason = data_in%freqmask_reason
     data_out%pixels        = data_in%pixels
     data_out%pix2ind       = data_in%pix2ind
     data_out%cal_method    = data_in%cal_method
     data_out%n_cal         = data_in%n_cal
     data_out%Tsys          = data_in%Tsys
     data_out%Phot          = data_in%Phot
     data_out%Thot          = data_in%Thot
     data_out%time_hot      = data_in%time_hot    
     data_out%sb_mean       = data_in%sb_mean
 !    data_out%spike_data    = data_in%spike_data
     data_out%n_nan         = data_in%n_nan
     data_out%n_pca_comp    = data_in%n_pca_comp
     data_out%n_pca_comp_feed    = data_in%n_pca_comp_feed
 
     data_out%use_freq_filter = data_in%use_freq_filter
 
     if (data_in%n_pca_comp > 0) then
        allocate(data_out%pca_ampl(size(data_in%nu,1),nsb,ndet,size(data_in%pca_ampl,4)))
        allocate(data_out%pca_comp(size(data_in%pca_comp,1),size(data_in%pca_comp,2)))
        allocate(data_out%pca_eigv(size(data_in%pca_eigv,1)))
 
        data_out%pca_ampl      = data_in%pca_ampl
        data_out%pca_comp      = data_in%pca_comp
        data_out%pca_eigv      = data_in%pca_eigv
     end if
 
 
     if (data_in%n_pca_comp_feed > 0) then
       allocate(data_out%pca_ampl_feed(size(data_in%nu, 1), nsb, ndet, size(data_in%pca_ampl_feed, 4)))
       allocate(data_out%pca_comp_feed(ndet, size(data_in%pca_comp_feed, 2), size(data_in%pca_comp_feed, 3)))
       allocate(data_out%pca_eigv_feed(ndet, size(data_in%pca_eigv_feed, 2)))
 
       data_out%pca_ampl_feed      = data_in%pca_ampl_feed
       data_out%pca_comp_feed      = data_in%pca_comp_feed
       data_out%pca_eigv_feed      = data_in%pca_eigv_feed
    end if
 
     if (data_in%use_freq_filter) then
        if (data_in%mask_outliers == 1) then
           allocate(data_out%corr_templ_ampl(4, nsb/2, ndet))
           data_out%corr_templ_ampl = data_in%corr_templ_ampl
        end if
        allocate(data_out%T_cont(nsamp_out,nsb,ndet,size(data_in%T_cont,4))) 
        allocate(data_out%dg(nsamp_out,nsb,ndet))
        data_out%T_cont        = data_in%T_cont
        data_out%dg            = data_in%dg
     end if
 
     
     ! Make angles safe for averaging
     do j = 1, ndet
        if (.not. is_alive(data_in%pixels(j))) cycle
        call make_angles_safe(data_in%point_tel(1,:,j), real(360.d0,sp)) ! Phi
        call make_angles_safe(data_in%point_tel(3,:,j), real(360.d0,sp)) ! Psi
        call make_angles_safe(data_in%point_cel(1,:,j), real(360.d0,sp)) ! Phi
        call make_angles_safe(data_in%point_cel(3,:,j), real(360.d0,sp)) ! Psi
     end do
 
     ! Compute variance per frequency channel
     !open(58,file='variance.dat')
 !    open(58,file='freqmask_2036.dat')
     data_out%tod_mean = 0.d0
     if (data_in%import_freqmask) then 
        write(*,*) "Using imported tod variance"
        data_out%tod_mean    = data_in%tod_mean
        data_out%var_fullres = data_in%var_fullres
     else
        do k = 1, ndet
           if (nsamp_out == 0) cycle
           if (.not. is_alive(data_out%pixels(k))) cycle
           do j = 1, nsb
              do i = 1, size(data_in%nu,1)
                 if (data_out%freqmask_full(i,j,k) == 0) then
                 !data_out%var_fullres(i,j,k) = 2.8d-5
                 data_out%var_fullres(i,j,k) = 0.d0
                    cycle
                 end if
                 data_out%tod_mean(i,j,k) = data_in%tod_mean(i,j,k)
 
                 !if (.not. all(data_in%tod(:,i,j,k) == data_in%tod(:,i,j,k))) then
                 !   write(*,*) "NaN in tod"
                 !end if
                 data_out%var_fullres(i,j,k) = variance(data_in%tod(:,i,j,k))
 !               if (data_out%var_fullres(i,j,k) > 2.8d-5) write(58,*) '   ', i, j, k
                 !write(58,*) k, j, i, data_out%var_fullres(i,j,k)
              end do
              !write(58,*)
            end do
        end do
     end if
 !    close(58)
 !    call mpi_finalize(ierr)
 !    stop
     
     data_out%nu = 0.d0
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
                 if (.not. is_alive(data_out%pixels(l))) cycle
                 data_out%nu(k,j,l)    = mean(data_in%nu((k-1)*dnu+1:k*dnu,j,l)) ! Frequency
                 data_out%tod(i,k,j,l) = 0.d0
                 weight                = 0.d0
                 do n = (k-1)*dnu+1, k*dnu
                    if (data_out%freqmask_full(n,j,l) == 0) cycle
                    if (data_out%var_fullres(n,j,l) <= 0) then
                       w = 0.d0
                    else
                       w  = 1.d0 / data_out%var_fullres(n,j,l) * data_in%freqmask_full(n,j,l)
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
     
     data_out%Tsys_lowres = 0.d0
     ! Calculate properly weighted lowres tsys
     do i = 1, ndet
        if (.not. is_alive(data_out%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, numfreq_out
              weight = 0.d0
              nw     = 0
              do n = (k-1)*dnu+1, k*dnu
                 if (data_out%freqmask_full(n,j,i) == 0) cycle
                 if (data_out%var_fullres(n,j,i) <= 0) then
                    w = 0.d0
                 else
                    w  = 1.d0 / data_out%var_fullres(n,j,i) * data_in%freqmask_full(n,j,i)
                    nw = nw + 1
                 end if
                 weight = weight + w 
                 ! must say that I do not understand this line (but I will trust my former self):
                 !data_out%Tsys_lowres(k,j,i) = data_out%Tsys_lowres(k,j,i) + (w * (data_in%Tsys(1,n,j,i) + data_in%Tsys(2,n,j,i))/2.d0) ** 2
                 data_out%Tsys_lowres(k,j,i) = data_out%Tsys_lowres(k,j,i) + (w * data_in%Tsys(n,j,i)) ** 2  
              end do
              data_out%Tsys_lowres(k,j,i) = sqrt(data_out%Tsys_lowres(k,j,i))
              if (weight > 0.d0) then
                 data_out%Tsys_lowres(k,j,i) = data_out%Tsys_lowres(k,j,i) * sqrt(1.d0 * nw) / weight
              else
                 data_out%Tsys_lowres(k,j,i) = 0.d0
              end if
           end do
        end do
     end do    
     
     if (data_in%import_freqmask) then 
        data_out%chi2 = data_in%chi2
     else
        data_out%chi2 = 0.d0
     ! calculate chisquared statistics on decimated data
        do i = 1, ndet
           if (.not. is_alive(data_out%pixels(i))) cycle
           do j = 1, nsb
              do k = 1, numfreq_out
                 if (data_out%freqmask(k,j,i) == 0) cycle
                 ! var = sum(data_out%tod(:,k,j,i) ** 2) 
                 ! var = sqrt(2.d0*nsamp_out)
                 ! var = 1.d0 / variance(data_out%tod(2:,k,j,i) - data_out%tod(:nsamp_out-1,k,j,i))
                 var = variance(data_out%tod(2:,k,j,i) - data_out%tod(:nsamp_out-1,k,j,i)) / 2.d0
                 if (var == 0.d0) then 
                    data_out%freqmask(k,j,i) = 0
                 else
                    data_out%chi2(k,j,i) = (sum(data_out%tod(:,k,j,i) ** 2) / var - nsamp_out) / sqrt(2.d0*nsamp_out)
                    if (data_out%chi2(k,j,i) > 5.d0) then   !!!! add to parameter file
                       data_out%freqmask(k,j,i) = 0.d0
                       data_out%freqmask_full((k-1)*dnu+1:k*dnu,j,i) = 0.d0
                       data_out%freqmask_reason((k-1)*dnu+1:k*dnu,j,i) = 50
                    end if
                 end if
              end do
           end do
        end do
     end if
     
     ! Polyfiltered TOD
     data_out%polyorder = data_in%polyorder
     if ((data_out%polyorder >= 0) .and. (.not. data_in%use_freq_filter)) then
        allocate(data_out%tod_poly(nsamp_out, 0:data_out%polyorder, nsb, ndet))
        do l = 1, ndet           
           if (.not. is_alive(data_out%pixels(l))) then
              data_out%tod_poly(:,:,:,l) = 0.d0
              cycle
           end if
           do j = 1, nsb
              do k = 0, data_out%polyorder
                 do i = 1, nsamp_out
                    data_out%tod_poly(i,k,j,l) = data_in%tod_poly(i,k,j,l) ! MUST BE CHANGED IF WE WILL DECIMATE IN TIME
                    !if (all(data_in%tod_poly((i-1)*dt+1:i*dt,k,j,l) == data_in%tod_poly((i-1)*dt+1:i*dt,k,j,l))) then
                    !   data_out%tod_poly(i,k,j,l) = mean(data_in%tod_poly((i-1)*dt+1:i*dt,k,j,l))
                    !else
                    !   data_out%tod_poly(i,k,j,l) = 0.d0
                    !end if
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
 
   subroutine initialize_fullres_frequency_mask(freqmaskfile, data, verb)
     implicit none
     character(len=*),                                intent(in)    :: freqmaskfile
     type(Lx_struct),                                 intent(inout) :: data
     logical(lgt),                                    intent(in)    :: verb
     integer(i4b) :: i, j, k, nfreq_full, nsb, ndet, unit, det, sb, freq, dfreq, ierr, max_det
     logical(lgt) :: first
     character(len=1024) :: line, val, equal
 
     max_det = get_num_dets()
     ndet = size(data%tod,4)
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
           if (.not. allocated(data%freqmask_full)) allocate(data%freqmask_full(nfreq_full,nsb,ndet))
           if (.not. allocated(data%freqmask_reason)) allocate(data%freqmask_reason(nfreq_full,nsb,ndet))
           data%freqmask_full = 1.d0
           data%freqmask_reason = 0
           first = .false.
        else
           read(line,*) freq, sb, det
           if (det == 0 .and. sb == 0 .and. freq == 0) then
              if (verb) then
                 write(*,*) 'ERROR: All frequencies removed by freqmask!'
              end if
           else if (det == 0 .and. sb == 0) then
              data%freqmask_full(freq,:,:) = 0.d0
           else if (det == 0 .and. freq == 0) then
              data%freqmask_full(:,sb,:) = 0.d0
           else if (sb == 0 .and. freq == 0) then
              data%freqmask_full(:,:,data%pix2ind(det)) = 0.d0
           else if (det == 0) then
              data%freqmask_full(freq,sb,:) = 0.d0
           else if (sb == 0) then
              data%freqmask_full(freq,:,data%pix2ind(det)) = 0.d0
           else if (freq == 0) then
              data%freqmask_full(:,sb,data%pix2ind(det)) = 0.d0
           else 
              data%freqmask_full(freq,sb,data%pix2ind(det)) = 0.d0
           end if
           if (all(data%freqmask_full == 0)) then
              if (verb) then
                 write(*,*) 'ERROR: All frequencies removed by freqmask!'
              end if
           end if
        end if
     end do
 
 99  close(unit)
 
     do k = 1, ndet
        if (.not. is_alive(data%pixels(k))) data%freqmask_full(:,:,k) = 0.d0
     end do
     data%freqmask_reason = nint(1.d0 - data%freqmask_full)  ! give all masked frequency reason = 1
   end subroutine initialize_fullres_frequency_mask
 
   subroutine postprocess_frequency_mask(nfreq, data, sid, verb)
     implicit none
     integer(i4b),                                    intent(in)    :: sid, nfreq
     type(Lx_struct),                                 intent(inout) :: data
     logical(lgt),                                    intent(in)    :: verb
     
     integer(i4b) :: i, j, k, nfreq_full, nsb, ndet, unit, det, sb, freq, dfreq, ierr, n_print
     logical(lgt) :: first
     character(len=1024) :: line, val, equal
 
     ndet       = size(data%tod,4) !get_num_dets()
     nsb        = get_num_sideband()
     nfreq_full = size(data%freqmask_full,1)
 
     ! Exclude frequencies with any NaNs                                                                                                                                                                                           
     do k = 1, ndet
        if (.not. is_alive(data%pixels(k))) cycle
        do j = 1, nsb
           n_print = 0
           do i = 1, nfreq_full
              if (data%freqmask_full(i,j,k) == 0.d0) cycle
              if (any(data%tod(:,i,j,k) .ne. data%tod(:,i,j,k))) then
                 if (verb) then
                    if (n_print < 5) then
                       write(*,fmt='(a,i8,3i6)') '   Rejecting NaNs, (sid,det,sb,freq) = ', sid, k,j,i
                       n_print = n_print + 1
                       if (n_print == 5) then
                          write(*,fmt='(a,i8,2i6)') "Suppressing further NaN print from this sideband", sid, k, j
                       end if
                    end if
                 end if
                 data%freqmask_full(i,j,k) = 0.d0
                 data%freqmask_reason(i,j,k) = 2
              end if
              if (all(data%tod(:,i,j,k) == 0.d0)) then
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
     
     do i = 1, ndet
        if (.not. is_alive(data%pixels(i))) then
           data%tod(:,:,:,i) = 0.d0
        end if
     end do
 
   end subroutine postprocess_frequency_mask
 
   subroutine subtract_pointing_templates(data, scan, verb) 
     implicit none
     type(lx_struct),  intent(inout) :: data
     type(comap_subscan), intent(in) :: scan
     logical(lgt),              intent(in)    :: verb
     integer(i4b)    :: n, ndet, nsb, nfreq, i, j, k, l, m, tmin, tmax, parfile_time
     real(dp)        :: a, sigma0, chisq, tsys, tau, dnu, const, g
     real(dp), dimension(:), allocatable :: el, az, dat
     
     if (scan%scanmode == 'stationary') then
        if (verb) then
           write(*,*) 'Scan mode is stationary, so no pointing templates removed', scan%id
        end if
        return
     else if (scan%scanmode == 'ces') then
        if (verb) then
           write(*,*) 'Scan mode is ces, so only az-template is removed', scan%id
        end if
     else 
        if (verb) then
           write(*,*) 'Scan mode is ', trim(scan%scanmode), ', so we remove both el- and az-templates', scan%id
        end if
     end if
 
     ! new linear function fit roughly every 5 min
     n = max(1, floor((data%time(size(data%time)) - data%time(1)) * 24 * 60 / 5.d0))  
     m = (size(data%time)) / n
     nfreq = size(data%tod,2)
     nsb   = size(data%tod,3)
     ndet  = size(data%tod,4)
     if (allocated(data%el_az_stats)) deallocate(data%el_az_stats)
     allocate(data%el_az_stats(2,n,nfreq,nsb,ndet))
 !    if (.not. allocated(data%el_az_stats)) allocate(data%el_az_stats(2,n,nfreq,nsb,ndet))
     data%el_az_stats = 0.d0
     
     !$OMP PARALLEL PRIVATE(j,k,l,i,el,az,dat,g,a,tmin,tmax)
     allocate(el(m), az(m), dat(m))
     !$OMP DO SCHEDULE(guided)
     do j = 1, nfreq
        do k = 1, ndet
           if (.not. is_alive(data%pixels(k))) cycle
           do l = 1, nsb
              if (data%freqmask_full(j,l,k) == 0.d0) cycle
              do i = 1, n    ! I think we could change the order of these loops to improve things 
                 tmin = (i-1)*m+1
                 tmax = i*m
                 if (i == n) then
                    tmax = size(data%time)
                 end if
                 az  = data%point_tel(1,tmin:tmax,k) 
                 call make_angles_safe(az, real(360.d0,dp))
                 
                 dat = data%tod(tmin:tmax,j,l,k)
                 !write(*,*) tmin, tmax, 'min max'
                 if (scan%scanmode == 'ces') then
                    call fit_az_template(az, dat, a)
                    data%tod(tmin:tmax,j,l,k) = &
                         data%tod(tmin:tmax,j,l,k) &
                         - a * az &
                         + sum(a * az) / size(az)
                    data%el_az_stats(2,i,j,l,k) = a
                 else 
                    el  = data%point_tel(2,tmin:tmax,k)
                    if (any(el == 0.d0)) write(*,*) "El is zero!!!", k, l, j, i
                    call fit_pointing_templates(el, az, dat, g, a)
                    data%tod(tmin:tmax,j,l,k) = &
                         data%tod(tmin:tmax,j,l,k) &
                         - g / sin(el*pi/180.) &
                         - a * az &
                         + sum(g / sin(el*pi/180.) + a * az) / size(el)
                    data%el_az_stats(1,i,j,l,k) = g
                    data%el_az_stats(2,i,j,l,k) = a
                 end if
              end do
           end do
        end do
     end do
     !$OMP END DO
     deallocate(el,az,dat)    
     !$OMP END PARALLEL
   end subroutine subtract_pointing_templates
 
 
  subroutine remove_elevation_gain(data) 
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: n, ndet, nsb, nfreq, i, j, k, l, m, tmin, tmax, parfile_time
    real(dp)        :: a, sigma0, chisq, tsys, tau, dnu, const, g
    real(dp), dimension(:), allocatable :: el, dat

!    m=a*data%samprate/data%scanfreq(2)      ! Number of tod-samples per gain estimate
    ! new linear function fit roughly every 5 min
    !write(*,*) "lasttime: ", data%time(size(data%time))
    !write(*,*) "firsttime: ", data%time(1)
    n = max(1, floor((data%time(size(data%time)) - data%time(1)) * 24 * 60 / 5.d0))  
    m = (size(data%time)) / n
!    write(*,*) 'Number of samples per gain estimate    =', m
!    n   = (size(data%time)+m-1)/m           ! Number of gain samples
!    n   = (size(data%time))/m               ! Number of gain samples
    
!    write(*,*) 'Number of gain estimates for this scan =', n
!    write(*,*) 'n*m                                    =', n*m
!    write(*,*) 'Total number of samples                =', size(data%time)
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)
!    write(*,*) nfreq, '= nfreq', nsb, '= nsb', ndet, '= ndet'
!    write(*,*) '---------------------------------------------------------'
    !allocate(data%time_gain(n), data%gain(n,nfreq,nsb,ndet))
    !data%time_gain = data%time(::m) ! OBS may reallocate length of time_gain!!!
!    open(13,file='gain.dat')
!    open(14,file='data1.dat')
!    open(15,file='el.dat')
!    open(16,file='data2.dat')
    !write(14,*) data%tod(:,30,1,1)

    !$OMP PARALLEL PRIVATE(j,k,l,i,el,dat,g,tmin,tmax)
    !allocate(dt(2*nsamp), dv(0:n-1))
    allocate(el(m), dat(m))
    !$OMP DO SCHEDULE(guided)    
    do j = 1, nfreq
       do k = 1, ndet
          if (.not. is_alive(data%pixels(k))) cycle
          do l = 1, nsb
             if (data%freqmask_full(j,l,k) == 0.d0) cycle
             do i = 1, n
                tmin = (i-1)*m+1
                tmax = i*m
                if (i == n) then
                   tmax = size(data%time)
                end if
                el  = data%point_tel(2,tmin:tmax,k)
                if (any(el == 0.d0)) write(*,*) "El is zero!!!", k, l, j, i
                dat = data%tod(tmin:tmax,j,l,k)
                !write(*,*) tmin, tmax, 'min max'
                call estimate_gain(el,dat,g)
!                if ((k == 1) .and. (l == 1) .and. (j == 30) .and. (i == 2)) then
!                   write(13,*) g
!                   write(*,*) 'Yo'
!                   !             write(13,*) (g-5.6d10)/5.6d10*100
!                   write(14,*) dat
!                   write(15,*) el
!                end if
                data%tod(tmin:tmax,j,l,k) = data%tod(tmin:tmax,j,l,k) - g * 1 / sin(el*pi/180.) + sum(g * 1 / sin(el*pi/180.)) / size(el)
             end do
          end do
       end do
    end do
    !$OMP END DO
!    write(16,*) data%tod(:,30,1,1)
!    close(13)
!    close(14)
!    close(15)
!    close(16)
    deallocate(el,dat)
    
    !$OMP END PARALLEL
  end subroutine remove_elevation_gain


  subroutine init_vanemask(data, vanemask, Thot, tsys_ind)
    implicit none
    type(Lx_struct),                                 intent(inout) :: data    
    integer(i4b), dimension(:), allocatable,         intent(inout) :: vanemask
    integer(i4b),                                    intent(in)    :: tsys_ind(2)
    real(dp),                                        intent(out)   :: Thot

    integer(i4b) :: i, j, nsamp_lowres, nsamp_highres, n_mean
    real(dp)     :: mjd_start, mjd_stop

    nsamp_lowres  = size(data%amb_state, 1)
    nsamp_highres = size(data%tod,1)

    if (.not. allocated(vanemask)) allocate(vanemask(nsamp_highres))
    
    vanemask = -1
    n_mean = 0
    Thot = 0.d0
    i = 1 ! Counter for highres grid
    j = 1 ! Couner for lowres grid
    do j = 1, nsamp_lowres-1
       if (data%amb_state(j) == 1) then ! .or. data%amb_state(j) == 3) then
          mjd_start = data%amb_time(j)
          mjd_stop  = data%amb_time(j+1)
          if (mjd_start > data%time(nsamp_highres)) exit
          do while (data%time(i) < mjd_start)
             i = i+1
          end do
          do while (data%time(i) < mjd_stop .and. i <= nsamp_highres)
             if (i > tsys_ind(1) .and. i < tsys_ind(2)) then
                vanemask(i) = data%amb_state(j)
                n_mean = n_mean + 1
                Thot = Thot + data%t_hot(j)
             end if
             i            = i+1
          end do
       end if
    end do
    if (n_mean > 0) then
       Thot = Thot / n_mean / 100.d0 + 273.15  ! in Kelvin
    end if
  end subroutine init_vanemask


  subroutine prepare_calibration(scan, data, is_sim, cal_method, cal_db, verb)
    implicit none
    type(Lx_struct),             intent(inout)    :: data
    logical(lgt),                intent(in)       :: is_sim, verb
    integer(i4b),                intent(in)       :: cal_method
    character(len=1024),         intent(in)       :: cal_db
    logical(lgt)                                  :: vane_in, no_nans
    type(hdf_file)                                :: file
    real(dp)                                      :: mean_tod, P_hot, P_cold
    real(dp), dimension(:,:,:,:), allocatable     :: tsys_fullres, buffer
    real(dp), dimension(:,:), allocatable         :: buff_T
    real(dp), dimension(:,:,:), allocatable       :: buff_time
    integer(i4b), dimension(:,:), allocatable     :: buff_succ
    integer(i4b)                                  :: nfreq_fullres, nsb, ndet, i, j, k, l, mjd_index1,mjd_index2, nsamp, dnu, n_hot, n_cold
    integer(i4b)                                  :: nsamp_gain(7), num_bin, n, mean_count, tsys_ind(2), n_tsys, cal_id(2), n_cal
    !integer(i4b), dimension(:), allocatable       :: scanID, vane_in_index, vane_out_index
    integer(i4b)                                  :: vane_in_index, vane_out_index
    real(dp)                                      :: mjd_high,w, sum_w_t, sum_w, t1, t2, tsys, tod_mean_dec, t_cold, Thot
    real(dp), dimension(:), allocatable           :: time, Y, tod_hot, tod_cold
    integer(i4b), dimension(:), allocatable       :: vanemask
    character(len=1024)                           :: db_string
    type(comap_scan_info) :: scan

    nsamp         = size(data%tod,1)
    nfreq_fullres = size(data%tod,2)
    nsb           = size(data%tod,3)
    ndet          = size(data%tod,4)


    

!    if (.not. allocated(data%Tsys)) allocate(data%Tsys(2, nfreq_fullres, nsb, ndet))
    if (.not. allocated(data%Phot)) allocate(data%Phot(3, 2, nfreq_fullres, nsb, ndet))
    if (.not. allocated(data%Thot)) allocate(data%Thot(3, 2))
    if (.not. allocated(data%time_hot)) allocate(data%time_hot(3, 2))
    
    data%Phot = 0.d0
    data%Thot = 0.d0
    data%time_hot = 0.d0 
    
    cal_id(1) = 1
    cal_id(2) = scan%nsub
    do m = 1, 2 
       tsys_ind(1) = minloc(abs(data%time - scan%ss(cal_id(m))%mjd(1)), 1) + 1
       tsys_ind(2) = minloc(abs(data%time - scan%ss(cal_id(m))%mjd(2)), 1) - 1
       ! write(*,*) "tsys_ind ", tsys_ind
       ! write(*,*) "cal_id ", cal_id
    
                                                                           
       call init_vanemask(data, vanemask, Thot, tsys_ind)
       n_hot = 0
!       write(*,*) "vanemask ", vanemask(tsys_ind(1):tsys_ind(2))
       
       vane_in = .false.   
       ! Find indices when vane is in
       do l = tsys_ind(1), tsys_ind(2)
          if (vanemask(l) == 1) then
             n_hot = n_hot + 1
!             write(*,*) "n_hot ", n_hot
             if (.not. vane_in) then 
                vane_in_index = l
                vane_in = .true.
                !write(*,*) "vane_in", vane_in
             end if
          else if (vane_in) then
             vane_out_index = l
!             write(*,*) "vane_out", l, vane_out_index
             exit
          end if
       end do
       if (n_hot == 0) then
          if (verb) then
             write(*,*) "no n_hot", cal_id(m)
          end if
       else
!          write(*,*) "in/out ", vane_in_index, vane_out_index
          data%time_hot(1, m) = sum(data%time(vane_in_index:vane_out_index - 1)) / (vane_out_index - vane_in_index)
!          write(*,*) "data%time_hot ", data%time_hot          
          data%Thot(1, m) = Thot
          do i=1, ndet
             if (.not. is_alive(data%pixels(i))) cycle
             do j=1, nsb
                do k=1, nfreq_fullres
                   if (data%freqmask_full(k,j,i) == 0.d0) cycle
                   data%Phot(1, m, k, j, i) = sum(data%tod(vane_in_index:vane_out_index - 1,k,j,i)) / (vane_out_index - vane_in_index)
                   !write(*,*) data%pixels(i), j, k, data%Phot(1, m, k, j, i)
                   !if (sum(data%tod(vane_in_index:vane_out_index - 1,k,j,i)) == sum(data%tod(vane_in_index:vane_out_index - 1,k,j,i))) then
                      !data%Phot(1, m, k, j, i) = median(data%tod(vane_in_index:vane_out_index - 1,k,j,i))
                   !end if
                   ! P_hot = 0.d0
                   ! no_nans = .true.
                   ! do l = vane_in_index, 
                   !    if (.not. data%tod(l,k,j,i) == data%tod(l,k,j,i)) then  !! check for NAN
                   !       no_nans = .false.
                   !    end if
                   ! end do
                   
                   ! P_hot  = median(tod_hot(1:n_hot))   !P_hot  / n_hot
                   ! if (is_sim) then 
                   !    data%Tsys(n_tsys,k,j,i) = 1.d0
                   ! else if (P_hot == 0.d0) then
                   !    data%Tsys(n_tsys,k,j,i) = 1.d0
                   !    if (verb) then
                   !       write(*,*) "P_hot = 0", P_hot, i, j, k
                   !    end if
                   ! else
                   !    data%Tsys(n_tsys,k,j,i) = P_hot!(t_hot-t_cold)/(Y(k)-1.d0)/P_cold
                   ! end if
                   !    !write(*,*) (t_hot-t_cold)/(Y(k)-1.d0)
                   ! end if
                end do
             end do
          end do
       end if
    end do
    deallocate(vanemask)
    
    if (.not. (trim(cal_db) == '') .and. (cal_method == 2)) then
       allocate(buffer(2, nfreq_fullres, nsb, 20))
       allocate(buff_T(2, 20))
       allocate(buff_time(2, 2, 20))
       allocate(buff_succ(2, 20))
       call open_hdf_file(cal_db, file, "r")
       write(db_string, '(A,I7.7)') 'obsid/', scan%id !15376 !scan%id
!       write(*,*) db_string, data%pixels
       call read_hdf(file, trim(db_string) // '/Phot', buffer)
       data%Phot(2, :, :, :, :) = buffer(:,:,:,data%pixels)   !!!!! is this right now?
       call read_hdf(file, trim(db_string) // '/Thot', buff_T) 
       call read_hdf(file, trim(db_string) // '/successful', buff_succ) 
       data%Thot(2, :) = buff_T(:, 1)
       !call read_hdf(file, trim(db_string) // '/calib_times', buff_time)
       if (buff_succ(1, 1) == 1) then
          data%time_hot(2, 1) = 0.5d0 * (scan%ss(cal_id(1))%mjd(1) + scan%ss(cal_id(1))%mjd(2))
       end if
       if (buff_succ(2, 1) == 1) then
          data%time_hot(2, 2) = 0.5d0 * (scan%ss(cal_id(2))%mjd(1) + scan%ss(cal_id(2))%mjd(2))
       end if
       call close_hdf_file(file)
       deallocate(buffer, buff_T, buff_time, buff_succ)
    end if

    data%cal_method = cal_method
    n_cal = 0 
    do j = 1, 2
       if (.not. (data%time_hot(data%cal_method,j) == 0.d0)) then
          n_cal = n_cal + 1
       end if
    end do
    data%n_cal = n_cal 
    if (n_cal < 1) then
       data%freqmask_full(:,:,:) = 0.d0
       if (verb) then
          write(*,*) "No valid ambient load measurments, Obsid: ", scan%id
       end if
    end if

  end subroutine prepare_calibration

  ! subroutine compute_P_hot(tsys_file, data, tsys_ind, n_tsys, is_sim, verb)
  !   implicit none
  !   character(len=*),            intent(in)       :: tsys_file
  !   type(Lx_struct),             intent(inout)    :: data
  !   logical(lgt),                intent(in)       :: is_sim, verb 
  !   type(hdf_file)                                :: file
  !   real(dp)                                      :: mean_tod, P_hot, P_cold
  !   real(dp), dimension(:,:,:,:), allocatable     :: tsys_fullres
  !   integer(i4b)                                  :: nfreq_fullres, nsb, ndet, i, j, k, l, mjd_index1,mjd_index2, nsamp, dnu, n_hot, n_cold
  !   integer(i4b)                                  :: nsamp_gain(7), num_bin, n, mean_count, tsys_ind(2), n_tsys
  !   integer(i4b), dimension(:), allocatable       :: scanID, vane_in_index, vane_out_index
  !   real(dp)                                      :: mjd_high,w, sum_w_t, sum_w, t1, t2, tsys, tod_mean_dec, t_cold, t_hot
  !   real(dp), dimension(:), allocatable           :: time, Y, tod_hot, tod_cold
  !   integer(i4b), dimension(:), allocatable       :: vanemask
    

  !   nsamp         = size(data%tod,1)
  !   nfreq_fullres = size(data%tod,2)
  !   nsb           = size(data%tod,3)
  !   ndet          = size(data%tod,4)


  !   if (.not. allocated(data%Tsys)) allocate(data%Tsys(2, nfreq_fullres, nsb, ndet))
  !   allocate(Y(nfreq_fullres), tod_hot(nsamp), tod_cold(nsamp))

  !   data%Tsys(n_tsys,:,:,:) = 1.d0
    
  !   !call locate_ambient_indecies(data, vane_in_index, vane_out_index, tsys_ind)                                                                                                                                                  
  !   call init_vanemask(data, vanemask, tsys_ind)

  !   do i=1, ndet
  !      if (.not. is_alive(data%pixels(i))) cycle
  !      do j=1, nsb
  !         do k=1, nfreq_fullres
  !            if (data%freqmask_full(k,j,i) == 0.d0) cycle

  !            P_hot = 0.d0
  !            n_hot = 0
  !            do l = 1, nsamp
  !               if (vanemask(l) == 1 .and. data%tod(l,k,j,i) == data%tod(l,k,j,i)) then
  !                  !P_hot = P_hot + data%tod(l,k,j,i)
  !                  n_hot = n_hot + 1
  !                  tod_hot(n_hot) = data%tod(l,k,j,i)
  !               else
  !                  cycle
  !               end if
  !            end do
  !            if (n_hot == 0) then
  !               data%Tsys(n_tsys,k,j,i) = 1.d0!(t_hot-t_cold)/(Y(k)-1.d0)/P_cold
  !               if (verb) then
  !                  write(*,*) "no n_hot", i, j, k
  !               end if
  !            else
  !               P_hot  = median(tod_hot(1:n_hot))   !P_hot  / n_hot
  !               if (is_sim) then 
  !                  data%Tsys(n_tsys,k,j,i) = 1.d0
  !               else if (P_hot == 0.d0) then
  !                  data%Tsys(n_tsys,k,j,i) = 1.d0
  !                  if (verb) then
  !                     write(*,*) "P_hot = 0", P_hot, i, j, k
  !                  end if
  !               else
  !                  data%Tsys(n_tsys,k,j,i) = P_hot!(t_hot-t_cold)/(Y(k)-1.d0)/P_cold
  !               end if
  !               !write(*,*) (t_hot-t_cold)/(Y(k)-1.d0)
  !            end if
  !         end do
  !      end do
  !   end do

  !   deallocate(vanemask, tod_hot, tod_cold)

  ! end subroutine compute_P_hot



  ! subroutine compute_Tsys_per_tp(tsys_file, data, tsys_ind, n_tsys, is_sim, verb)
  !   implicit none
  !   character(len=*),            intent(in)       :: tsys_file
  !   type(Lx_struct),             intent(inout)    :: data
  !   logical(lgt),                intent(in)       :: is_sim, verb 
  !   type(hdf_file)                                :: file
  !   real(dp)                                      :: mean_tod, P_hot, P_cold
  !   real(dp), dimension(:,:,:,:), allocatable     :: tsys_fullres
  !   integer(i4b)                                  :: nfreq_fullres, nsb, ndet, i, j, k, l, mjd_index1,mjd_index2, nsamp, dnu, n_hot, n_cold
  !   integer(i4b)                                  :: nsamp_gain(7), num_bin, n, mean_count, tsys_ind(2), n_tsys
  !   integer(i4b), dimension(:), allocatable       :: scanID, vane_in_index, vane_out_index
  !   real(dp)                                      :: mjd_high,w, sum_w_t, sum_w, t1, t2, tsys, tod_mean_dec, t_cold, t_hot
  !   real(dp), dimension(:), allocatable           :: time, Y, tod_hot, tod_cold
  !   integer(i4b), dimension(:), allocatable       :: vanemask
    

  !   nsamp         = size(data%tod,1)
  !   nfreq_fullres = size(data%tod,2)
  !   nsb           = size(data%tod,3)
  !   ndet          = size(data%tod,4)


  !   if (.not. allocated(data%Tsys)) allocate(data%Tsys(2, nfreq_fullres, nsb, ndet))
  !   allocate(Y(nfreq_fullres), tod_hot(nsamp), tod_cold(nsamp))

  !   data%Tsys(n_tsys,:,:,:) = 40.d0
  !   t_cold = 2.73
  !   t_hot = mean(data%t_hot)/100 + 273.15

  !   !call locate_ambient_indecies(data, vane_in_index, vane_out_index, tsys_ind)                                                                                                                                                  
  !   call init_vanemask(data, vanemask, tsys_ind)

  !   do i=1, ndet
  !      if (.not. is_alive(data%pixels(i))) cycle
  !      do j=1, nsb
  !         do k=1, nfreq_fullres
  !            if (data%freqmask_full(k,j,i) == 0.d0) cycle

  !            P_hot = 0.d0
  !            P_cold = 0.d0
  !            n_hot = 0
  !            n_cold = 0
  !            do l = 1, nsamp
  !               if (vanemask(l) == 2 .and. data%tod(l,k,j,i) == data%tod(l,k,j,i)) then
  !                  !P_hot = P_hot + data%tod(l,k,j,i)
  !                  n_hot = n_hot + 1
  !                  tod_hot(n_hot) = data%tod(l,k,j,i)
  !               else if (vanemask(l) == 3 .and. data%tod(l,k,j,i) == data%tod(l,k,j,i)) then
  !                  !P_cold = P_cold + data%tod(l,k,j,i)
  !                  n_cold = n_cold + 1
  !                  tod_cold(n_cold) = data%tod(l,k,j,i)
  !               else
  !                  cycle
  !               end if
  !            end do
  !            if ((n_hot == 0) .or. (n_cold == 0)) then
  !               data%Tsys(n_tsys,k,j,i) = 1.d0!(t_hot-t_cold)/(Y(k)-1.d0)/P_cold
  !               if (verb) then
  !                  write(*,*) "no n_cold", i, j, k
  !               end if
  !            else
  !               P_hot  = median(tod_hot(1:n_hot))   !P_hot  / n_hot
  !               P_cold = median(tod_cold(1:n_cold)) !P_cold / n_cold
  !               Y(k)   = P_hot/P_cold
  !               if (is_sim) then 
  !                  data%Tsys(n_tsys,k,j,i) = 1.d0
  !               else if ((Y(k) == 1.d0) .or. (P_cold == 0.d0)) then
  !                  data%Tsys(n_tsys,k,j,i) = 1.d0
  !                  if (verb) then
  !                     write(*,*) "Y = 1 or P_cold = 0", P_hot, P_cold, i, j, k
  !                  end if
  !               else
  !                  data%Tsys(n_tsys,k,j,i) = (t_hot-t_cold)/(Y(k)-1.d0)/P_cold
  !               end if
  !               !write(*,*) (t_hot-t_cold)/(Y(k)-1.d0)
  !            end if
  !         end do
  !      end do
  !   end do

  !   deallocate(vanemask, tod_hot, tod_cold)

  ! end subroutine compute_Tsys_per_tp

 
  subroutine find_tsys(data_l2_fullres, is_sim, sim_tsys, scan, verb, rank)
   implicit none
   type(lx_struct), intent(inout) :: data_l2_fullres
   type(comap_subscan), intent(in) :: scan
   real(dp),        intent(in)    :: sim_tsys
   integer(i4b),        intent(in)    :: rank
!    real(dp),    intent(in)        :: tsys_time(2)
   logical(lgt), intent(in)       :: is_sim, verb
   logical(lgt)                   :: tsys_nan
   real(dp)                       :: interp1d_P_hot, interp1d_T_hot, mean_tod, y0, y1, x0, x1, x, scan_time
   real(dp)                       :: T_cold, T_hot(2), tsys, time_hot(2)
   integer(i4b) :: i, j, k, ndet, nfreq, nsb, l2_nsamp, n_print, cal_method
   l2_nsamp = size(data_l2_fullres%tod, 1)
   nfreq    = size(data_l2_fullres%tod,2)
   nsb      = size(data_l2_fullres%tod,3)
   ndet     = size(data_l2_fullres%tod,4)
   T_cold = 2.725
   scan_time = data_l2_fullres%time(l2_nsamp/2)
   
   tsys_nan = .false.
  
   if (.not. allocated(data_l2_fullres%Tsys)) allocate(data_l2_fullres%Tsys(nfreq, nsb, ndet))
   data_l2_fullres%Tsys = 0
   
   cal_method = data_l2_fullres%cal_method
   T_hot = data_l2_fullres%Thot(cal_method, :)
   time_hot = data_l2_fullres%time_hot(cal_method, :)
   if (data_l2_fullres%n_cal == 2) then
      x0 = time_hot(1); x1 = time_hot(2)
      if (x0 == x1) then
         data_l2_fullres%freqmask_full(:,:,:) = 0.d0
         write(*,*) "Should not have been n_cal=2!! ", x0, x1, cal_method, scan%id
         return
      end if
      y0 = T_hot(1); y1 = T_hot(2)
      interp1d_T_hot = (y0*(x1-scan_time) + y1*(scan_time - x0))/(x1-x0)
   end if

   do i = 1, ndet
      if (.not. is_alive(data_l2_fullres%pixels(i))) cycle
      do j = 1, nsb
         do k = 1, nfreq
            if (data_l2_fullres%freqmask_full(k,j,i) == 0.d0) cycle
            mean_tod = data_l2_fullres%tod_mean(k,j,i)
               
            if (is_sim) then
               data_l2_fullres%tsys(k,j,i)  = sim_tsys
            else if (mean_tod == 0.d0) then
               data_l2_fullres%freqmask_full(k,j,i) = 0.d0
            else if (x0 == x1) then
               data_l2_fullres%freqmask_full(k,j,i) = 0.d0
            else if (data_l2_fullres%n_cal == 2) then  ! Two good ambient load measurements
               ! Interpolate Tsys to current time for each detector
               
               y0 = data_l2_fullres%Phot(cal_method,1,k,j,i); y1 = data_l2_fullres%Phot(cal_method,2,k,j,i)
               interp1d_P_hot = (y0*(x1-scan_time) + y1*(scan_time - x0))/(x1-x0)
               if (interp1d_P_hot == mean_tod) then
                  data_l2_fullres%freqmask_full(k,j,i) = 0.d0
                  write(*,*) "Phot = Pcold", i, j, k, scan%id
                  cycle
               end if
               tsys = (interp1d_T_hot - T_cold)/(interp1d_P_hot/mean_tod-1.d0)
               data_l2_fullres%Tsys(k,j,i) = tsys
            else
               if (time_hot(2) == 0.d0) then  ! Only first ambient measurement is good
                  tsys = (T_hot(1) - T_cold) / (data_l2_fullres%Phot(cal_method,1,k,j,i) / mean_tod - 1.d0)
                  data_l2_fullres%Tsys(k,j,i) = tsys
               else                           ! Only second ambient measurement is good
                  tsys = (T_hot(2) - T_cold) / (data_l2_fullres%Phot(cal_method,2,k,j,i) / mean_tod - 1.d0)
                  data_l2_fullres%Tsys(k,j,i) = tsys
               end if
            end if
           
            ! Masking failed Tsys measurements
            if (tsys .ne. tsys) then
               write(*,*) "Rank: ", rank, "Found NaN in Tsys!!!"
               tsys_nan = .true.
               data_l2_fullres%freqmask_full(k,j,i) = 0.d0
               data_l2_fullres%Tsys(k,j,i) = tsys
            end if
         end do
      end do
   end do
   end subroutine find_tsys
 
   subroutine calibrate_tod(data_l2_fullres, max_tsys)
     implicit none
     type(lx_struct), intent(inout) :: data_l2_fullres
     real(dp), intent(in)           :: max_tsys
     integer(i4b) :: i, j, k, ndet, nfreq, nsb
     nfreq    = size(data_l1%tod,2)
     nsb      = size(data_l1%tod,3)
     ndet     = size(data_l1%tod,4)
 
     do i = 1, ndet
        if (.not. is_alive(data_l2_fullres%pixels(i))) cycle
        do j = 1, nsb
           do k = 1, nfreq
              if (data_l2_fullres%freqmask_full(k,j,i) == 0.d0) cycle
              if (data_l2_fullres%Tsys(k,j,i) .ne. data_l2_fullres%Tsys(k,j,i)) then
                 write(*,*) "tsys is nan!!", i, j, k
                 data_l2_fullres%freqmask_full(k,j,i) = 0.d0
                 data_l2_fullres%freqmask_reason(k,j,i) = 70
              else if (data_l2_fullres%Tsys(k,j,i) > max_tsys) then
                 data_l2_fullres%freqmask_full(k,j,i) = 0.d0
                 data_l2_fullres%freqmask_reason(k,j,i) = 60
              else
                 data_l2_fullres%tod(:,k,j,i)  = data_l2_fullres%Tsys(k,j,i) * data_l2_fullres%tod(:,k,j,i)
              end if
           end do
        end do
     end do
   end subroutine calibrate_tod
 
   subroutine fit_noise(data_l2)
     implicit none
 
     type(lx_struct), intent(inout) :: data_l2
 
     integer(i4b) :: i, j, k, nsamp, ndet, nfreq, nsb
     real(dp)     :: var
     nsamp = size(data_l2%tod,1)
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
              var = variance(data_l2%tod(2:,k,j,i) - data_l2%tod(:nsamp-1,k,j,i)) / 2
              data_l2%sigma0(k,j,i) = sqrt(var)!sqrt(variance(data_l2%tod(:,k,j,i)))
              data_l2%alpha(k,j,i)  = -2.d0
              data_l2%fknee(k,j,i)  = 1.d-6   ! Set to a very low value
           end do
        end do
     end do
 
   end subroutine fit_noise
 
 
 end program
 