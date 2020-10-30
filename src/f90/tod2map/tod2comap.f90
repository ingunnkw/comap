program tod2comap
  use comap_lx_mod
  use comap_map_mod
  use comap_scan_mod
  use comap_acceptlist_mod
  use comap_jackknife_mod
  use comap_patch_mod
  use comap_ephem_mod
  use quiet_fft_mod
  use quiet_hdf_mod
  use tod2comap_mapmaker
  use rngmod
  !use tod2comap_utils
  implicit none

!  include "mpif.h"

  !type tod_type
  !   real(dp)     :: samprate, Tsys
  !   integer(i4b) :: nsamp, ndet, nfreq, nsb
  !   real(dp)     :: fmin, fmax, df
  !   real(dp), allocatable, dimension(:)       :: t                ! (time or freq)
  !   real(dp), allocatable, dimension(:,:,:,:) :: d, d_raw, g, rms ! (time, freq,  sb, det) 
  !   real(dp), allocatable, dimension(:,:)     :: point, f         ! (3, time) or (sb, freq)
  !end type tod_type


  type(tod_type)        :: tod
  type(map_type)        :: map_tot, map_scan, map_obs, buffer!, map_split1, map_split2, buffer1, buffer2
  type(map_type), allocatable, dimension(:) :: map_jk, buffer_jk
  type(comap_scan_info) :: scan
  !type(acceptlist)      :: alist
  type(patch_info)      :: pinfo
  type(jk_type)         :: jk_info

  integer(i4b), allocatable, dimension(:,:) :: pixels
  character(len=512)    :: filename, map_filename, parfile, acceptfile, prefix, pre, map_name, object, coord_system, l1file
  character(len=512)    :: sim_filename, prefix_sim, pre_sim, sim_name, jackknife, jk_def_file, acc_id, map_file1, map_file2
  character(len=6)      :: obsid
  character(len=8)      :: scanid
  character(len=5)      :: sim_string
  integer(i4b)          :: nscan, nsub, i, j, k, det, sb, freq, sim, nsim, n1, nn1, n2, nn2, split, scan_index

  integer(i4b)          :: myid, nproc, ierr, root
  logical               :: binning_split, found, obs_map, scan_map, split_mode, verbose, use_acc
  real(dp), allocatable, dimension(:,:) :: offsets
  real(dp)              :: my_x_max, my_y_max

  type(planck_rng)      :: rng_handle
  integer(i4b)          :: seed
 
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid,  ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)
  root = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! mpi per scan, all detectors
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1, parfile)
  !call initialize_comap_ephem_mod(parfile)
  !write(*,*) get_obj_info('jupiter', 58526.592246526d0)
  !stop
  call get_parameter(0, parfile, 'MAP_DIR', par_string=pre)
  call get_parameter(0, parfile, 'SIM_DIR', par_string=pre_sim)
  call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)
  call get_parameter(0, parfile, 'MAP_NAME', par_string=map_name)
  call get_parameter(0, parfile, 'SIM_NAME', par_string=sim_name)
  call get_parameter(0, parfile, 'N_NOISE_SIMULATIONS', par_int=nsim)
  call get_parameter(0, parfile, 'SEED', par_int=seed)
  !call get_parameter(0, parfile, 'JACKKNIVES', par_string=jackknife)
  !call get_parameter(0, parfile, 'NUMBER_JK', par_int=njk)
  call get_parameter(0, parfile, 'OBSID_MAPS', par_lgt=obs_map)
  call get_parameter(0, parfile, 'SCAN_MAPS', par_lgt=scan_map)
  call get_parameter(0, parfile, 'VERBOSE_PRINT', par_lgt=verbose)
  call get_parameter(0, parfile, 'USE_ACCEPT', par_lgt=use_acc)
  if (use_acc) then 
     call get_parameter(0, parfile, 'ACCEPT_DATA_FOLDER', par_string=acceptfile)
     call get_parameter(0, parfile, 'ACCEPT_DATA_ID_STRING', par_string=acc_id)
     call get_parameter(0, parfile, 'JK_DEF_FILE', par_string=jk_def_file)
     if (trim(acc_id) .ne. '') acc_id = '_' // acc_id
     acceptfile = trim(acceptfile) // 'jk_data' // trim(acc_id) // '_' // trim(object) // '.h5'
  end if


  call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)

  !call get_parameter(0, parfile, 'BIN_SPLIT', par_)
  !call get_parameter()
  call initialize_comap_patch_mod(parfile)
  call initialize_detector_mod(parfile)
  call initialize_scan_mod(parfile, object)
  if (use_acc) then
     call read_acceptlist(jk_info, acceptfile, jk_def_file)
  else
     call initialize_empty_jk(jk_info)
  end if
  !call initialize_accept_list(trim(acceptfile), alist)
  found = get_patch_info(object, pinfo)
  if (.not. found) then
     write(*,*) "Error: patch not found"
     call mpi_finalize(ierr)
     stop
  end if
  nscan = get_num_scans()
  if (myid == 0) write(*,*) nscan
  !stop
  !allocate(tod(nscan))
  !nscan = 1
  !call free_map_type(map)


  !allocate(alist%status(tod(i)%nfreq, tod(i)%ndet))
  !alist%status = 0

  !if (jackknife .ne. 'none') then
  !   split_mode = .true.
  !else
  !   split_mode = .false.
  !end if
  !n1 = 0
  !n2 = 0

  !write(*,*) jk_info%split(:,1,2,75)
  !stop


  if (myid==0) write(*,*) "Initialising mapmaker"
  call initialize_mapmaker(map_scan, parfile, pinfo, jk_info)
  call initialize_mapmaker(map_obs,  parfile, pinfo, jk_info)
  call initialize_mapmaker(map_tot,  parfile, pinfo, jk_info)
  call initialize_mapmaker(buffer,   parfile, pinfo, jk_info)
  call nullify_map_type(map_tot); call nullify_map_type(buffer) ! Should not be necessary

  !!if (jk_info%njk .gt. 0) then
  !!   allocate(map_jk(2*jk_info%njk), buffer_jk(2*jk_info%njk))
  !!   do i = 1, 2*jk_info%njk
  !!      call initialize_mapmaker(map_jk(i), parfile, pinfo, jk_info)
  !!      call initialize_mapmaker(buffer_jk(i), parfile, pinfo, jk_info)
  !!   end do
  !!end if
 

  ! This loop currently requiers that all obsIDs are of the same patch
  do i = 1+myid, nscan, nproc 
     !write(*,*) myid, i, 'of', nscan
     !if (allocated(alist%status)) deallocate(alist%status)
     call get_scan_info(i,scan)
     !write(*,*) scan%id

     write(*,*) myid, scan%id
     call nullify_map_type(map_obs)

     nsub = scan%nsub
     do j = 2, nsub-1
        call nullify_map_type(map_scan)
     
        ! Get TOD / read level 2 file
        !nsub = scan%nsub
        !do j = 2, nsub-1
        !write(*,*) myid, 'scan', i, 'of', nscan, 'subscan', j, 'of', nsub
        filename = scan%ss(j)%l2file !write(*,*) trim(filename)
        
        call get_tod(trim(filename), tod, parfile)
        !elevation cuts here
        !if (tod%mean_el .lt. 35.d0) cycle
        !if (tod%mean_el .gt. 65.d0) cycle
        
        !call int2string(scan%ss(j)%id, scanid)
        !write(*,*) 'scan', scanid
     
        if (use_acc) then
           scan_index = findloc(jk_info%scan_list, scan%ss(j)%id, dim=1)
        else
           scan_index = 1
        end if  

        !call nullify_map_type(map_scan)
        !call time2pix(tod, map_scan, parfile, pinfo, jk_info%jk_list(:,:,scan_index))
        !call time2pix(tod, map_tot, parfile, pinfo, jk_info%jk_list(:,:,scan_index))
        if (use_acc) then
           scan_index = findloc(jk_info%scan_list, scan%ss(j)%id, dim=1)
        !   do k = 1, jk_info%njk
        !      !write(*,*) 2*k, 2*k-1
        !      if (jk_info%split(k,scan_index) .eq. 0) then
        !         call time2pix(tod, map_jk(2*k -1), parfile, pinfo, jk_info%jk_list(:,:,scan_index))
        !      else
        !         !if (k==2) write(*,*) 'yes'
        !         call time2pix(tod, map_jk(2*k), parfile, pinfo, jk_info%jk_list(:,:,scan_index))
        !      end if
        !   end do
        end if
        if (verbose) write(*,*) myid, "making maps, obsID", i, 'scan', scan%ss(j)%id
        call binning(map_scan, tod, i, parfile, pinfo, jk_info%jk_list(:,:,scan_index), jk_info%split(:,:,:,scan_index))
        !call finalize_scan_binning(map_scan)
        !prefix = trim(pre)//trim(scan%object)//'_'//trim(scanid)
        !call output_submap_h5(trim(prefix), map_scan)
        !!call free_map_type(map_scan)
        !call free_tod_type(tod)
        !end do

        if (obs_map) then   
           ! Add to obsID maps
           map_obs%dsum      = map_obs%dsum + map_scan%dsum
           map_obs%div       = map_obs%div  + map_scan%div
           map_obs%nhit      = map_obs%nhit + map_scan%nhit
           map_obs%dsum_co   = map_obs%dsum_co + map_scan%dsum_co
           map_obs%div_co    = map_obs%div_co  + map_scan%div_co
           map_obs%nhit_co   = map_obs%nhit_co + map_scan%nhit_co
           map_obs%dsum_jk   = map_obs%dsum_jk + map_scan%dsum_jk 
           map_obs%div_jk    = map_obs%div_jk  + map_scan%div_jk 
           map_obs%nhit_jk   = map_obs%nhit_jk + map_scan%nhit_jk
           map_obs%dsum_jkco = map_obs%dsum_jkco + map_scan%dsum_jkco 
           map_obs%div_jkco  = map_obs%div_jkco  + map_scan%div_jkco
           map_obs%nhit_jkco = map_obs%nhit_jkco + map_scan%nhit_jkco
           map_obs%dsum_split = map_obs%dsum_split + map_scan%dsum_split 
           map_obs%div_split  = map_obs%div_split  + map_scan%div_split
           map_obs%nhit_split = map_obs%nhit_split + map_scan%nhit_split
        end if

        ! Add to total map
        map_tot%dsum      = map_tot%dsum + map_scan%dsum
        map_tot%div       = map_tot%div  + map_scan%div
        map_tot%nhit      = map_tot%nhit + map_scan%nhit
        map_tot%dsum_co   = map_tot%dsum_co + map_scan%dsum_co
        map_tot%div_co    = map_tot%div_co  + map_scan%div_co
        map_tot%nhit_co   = map_tot%nhit_co + map_scan%nhit_co
        map_tot%dsum_jk   = map_tot%dsum_jk + map_scan%dsum_jk
        map_tot%div_jk    = map_tot%div_jk  + map_scan%div_jk
        map_tot%nhit_jk   = map_tot%nhit_jk + map_scan%nhit_jk
        map_tot%dsum_jkco = map_tot%dsum_jkco + map_scan%dsum_jkco
        map_tot%div_jkco  = map_tot%div_jkco  + map_scan%div_jkco
        map_tot%nhit_jkco = map_tot%nhit_jkco + map_scan%nhit_jkco
        map_tot%dsum_split = map_tot%dsum_split + map_scan%dsum_split
        map_tot%div_split  = map_tot%div_split  + map_scan%div_split
        map_tot%nhit_split = map_tot%nhit_split + map_scan%nhit_split
        !if (use_acc) then
        !   ! Add to jackknives
        !   do k = 1, jk_info%njk
        !      if (jk_info%split(k,scan_index) .eq. 0) then
        !         map_jk(2*k-1)%dsum    = map_jk(2*k-1)%dsum    + map_scan%dsum
        !         map_jk(2*k-1)%div     = map_jk(2*k-1)%div     + map_scan%div
        !         map_jk(2*k-1)%dsum_co = map_jk(2*k-1)%dsum_co + map_scan%dsum_co
        !         map_jk(2*k-1)%div_co  = map_jk(2*k-1)%div_co  + map_scan%div_co
        !         !n1 = n1 + 1
        !      else
        !         map_jk(2*k)%dsum    = map_jk(2*k)%dsum    + map_scan%dsum
        !         map_jk(2*k)%div     = map_jk(2*k)%div     + map_scan%div
        !         map_jk(2*k)%dsum_co = map_jk(2*k)%dsum_co + map_scan%dsum_co
        !         map_jk(2*k)%div_co  = map_jk(2*k)%div_co  + map_scan%div_co
        !         !n2 = n2 + 1
        !      end if
        !   end do
        !end if

        ! Maps per scan
        if (scan_map) then
           call finalize_binning(map_scan)
           call int2string(scan%ss(j)%id, scanid)
           prefix = trim(pre)//trim(scan%object)//'_'//trim(scanid)
           map_filename = trim(prefix)//'_'//trim(map_name)//'.h5'
           call output_map_h5(map_filename, map_scan)
        end if

     end do ! end loop over scans

     call int2string(scan%id, obsid)

     ! Maps per obsID
     if (obs_map) then
        call finalize_binning(map_obs)
        prefix = trim(pre)//trim(scan%object)//'_'//trim(obsid)
        map_filename = trim(prefix)//'_'//trim(map_name)//'.h5'
        call output_map_h5(map_filename, map_obs)
     end if
     call free_tod_type(tod)

     ! THIS NEEDS TO BE UPDATED !!!
     !do sim=1, nsim
     !   call nullify_map_type(map_scan)
     !   write(*,*) myid, "making sim", sim, 'of', nsim
     !   do j = 2, nsub-1

     !      filename = scan%ss(j)%l2file

     !      call get_sim(trim(filename), tod, parfile, rng_handle)
     !      call time2pix_sim(tod, map_scan, parfile, pinfo)
     !      call time2pix_sim(tod, map_tot, parfile, pinfo)
     !      call binning_sim(map_tot, map_scan, tod, i, parfile, pinfo)

     !   end do
     !   if (obs_map) then 
     !      call finalize_binning(map_scan)

     !      call int2string(sim, sim_string)
     !      prefix_sim = trim(pre_sim)//trim(scan%object)//'_'//trim(obsid)
     !      sim_filename = trim(prefix_sim)//'_'//trim(sim_name)//'_'//trim(sim_string)//'.h5'
     !      call output_submap_sim_h5(sim_filename, map_scan, sim)
     !   end if
     !   call free_tod_type(tod)
     !end do
     !write(*,*) maxval(map_tot%dsum), maxval(map_tot%dsum_jk)
     !call mpi_finalize(ierr)
     !stop
     
  end do ! end loop over obsIDs

  call free_map_type(map_scan)
  call free_map_type(map_obs)
  call free_jk_type(jk_info)


  !call mpi_reduce(map_tot%div, buffer%div, size(map_tot%div), MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%div,     buffer%div,     size(map_tot%div),     MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%dsum,    buffer%dsum,    size(map_tot%dsum),    MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%nhit,    buffer%nhit,    size(map_tot%nhit),    MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)

  ! Co-added over feeds 
  call mpi_allreduce(map_tot%div_co,  buffer%div_co,  size(map_tot%div_co),  MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%dsum_co, buffer%dsum_co, size(map_tot%dsum_co), MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%nhit_co, buffer%nhit_co, size(map_tot%nhit_co), MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)

  ! Jackknives
  call mpi_allreduce(map_tot%div_jk,  buffer%div_jk,  size(map_tot%div_jk),  MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%dsum_jk, buffer%dsum_jk, size(map_tot%dsum_jk), MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%nhit_jk, buffer%nhit_jk, size(map_tot%nhit_jk), MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%div_jkco,  buffer%div_jkco,  size(map_tot%div_jkco),  MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%dsum_jkco, buffer%dsum_jkco, size(map_tot%dsum_jkco), MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%nhit_jkco, buffer%nhit_jkco, size(map_tot%nhit_jkco), MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%div_split,  buffer%div_split,  size(map_tot%div_split),  MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%dsum_split, buffer%dsum_split, size(map_tot%dsum_split), MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(map_tot%nhit_split, buffer%nhit_split, size(map_tot%nhit_split), MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)

  !!do i = 1, 2*jk_info%njk
  !!   call mpi_allreduce(map_jk(i)%div,     buffer_jk(i)%div,     size(map_tot%div),     MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  !!   call mpi_allreduce(map_jk(i)%dsum,    buffer_jk(i)%dsum,    size(map_tot%dsum),    MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  !!   call mpi_allreduce(map_jk(i)%nhit,    buffer_jk(i)%nhit,    size(map_tot%nhit),    MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
  !!   call mpi_allreduce(map_jk(i)%div_co,  buffer_jk(i)%div_co,  size(map_tot%div_co),  MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  !!   call mpi_allreduce(map_jk(i)%dsum_co, buffer_jk(i)%dsum_co, size(map_tot%dsum_co), MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  !!   call mpi_allreduce(map_jk(i)%nhit_co, buffer_jk(i)%nhit_co, size(map_tot%nhit_co), MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
  !!end do

  ! Simulations
  !call mpi_allreduce(map_tot%div_sim,  buffer%div_sim,  size(map_tot%div_sim),  MPI_REAL, MPI_SUM, mpi_comm_world, ierr)
  !call mpi_allreduce(map_tot%dsum_sim, buffer%dsum_sim, size(map_tot%dsum_sim), MPI_REAL, MPI_SUM, mpi_comm_world, ierr)


  if (myid == 0) then
     map_tot%div  = buffer%div
     map_tot%dsum = buffer%dsum
     map_tot%nhit = buffer%nhit

     ! Co-added over feeds
     map_tot%div_co  = buffer%div_co
     map_tot%dsum_co = buffer%dsum_co
     map_tot%nhit_co = buffer%nhit_co
     !write(*,*) 'sum', sum(abs(map_tot%dsum)), sum(abs(map_tot%div))

     ! Jackknives
     map_tot%div_jk  = buffer%div_jk
     map_tot%dsum_jk = buffer%dsum_jk
     map_tot%nhit_jk = buffer%nhit_jk
     map_tot%div_jkco  = buffer%div_jkco
     map_tot%dsum_jkco = buffer%dsum_jkco
     map_tot%nhit_jkco = buffer%nhit_jkco
     map_tot%div_split  = buffer%div_split
     map_tot%dsum_split = buffer%dsum_split
     map_tot%nhit_split = buffer%nhit_split

     ! Simulations
     !map_tot%div_sim  = buffer%div_sim
     !map_tot%dsum_sim = buffer%dsum_sim


     write(*,*) "Finalising"
     ! finalize_jk ???????????
     call finalize_binning(map_tot)

     !write(*,*) maxval(map_tot%m), maxval(map_tot%m_jk)
     prefix = trim(pre)//trim(scan%object)
     map_filename = trim(prefix)//'_'//trim(map_name)//'.h5'
     call output_map_h5(map_filename, map_tot)
     call free_map_type(map_tot)

     if (allocated(map_jk)) deallocate(map_jk, buffer_jk)


  end if

  !deallocate(tod)
  if (myid == 0) write(*,*) 'Done'
  call mpi_finalize(ierr)
  

  !if (myid == 0) write(*,*) trim(scan%object), trim(itoa(scan%sid))
  !prefix = trim(pre)//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID
  !if (myid == 0) write(*,*) prefix
  !call mpi_finalize(ierr)
  !stop

  !if (.not. binning_split) then
  !call initialize_mapmaker(map_tot, tod, parfile)
  !call time2pix(tod, map)

  ! Compute and co-add maps
  !call binning(map, tod(i), alist)
  ! if (myid == 0) write(*,*) "mapmaker ..."
  ! do det = 3, 3
  !    do sb = 1, tod(1)%nsb
  !       if (myid == 0) write(*,*) "sb", sb
  !       !do freq = tod(1)%nfreq/4, tod(1)%nfreq/4
  !       !do freq = 30,30
  !       do freq = 1, tod(1)%nfreq
  !          if (myid == 0 .and. modulo(freq, 10) == 0) write(*,*) 'freq', freq, 'of', tod(1)%nfreq
  !          !call pcg_mapmaker(tod, map, alist, det, sb, freq, parfile)
  !          !call binning(map_tot, map_scan, tod, alist)
  !       end do
  !    end do
  ! end do
  ! call finalize_binning(map)

  ! !end do
  ! if (myid == 0) write(*,*) "Writing to file ..."
  ! !write(*,*) trim(itoa(scan%sid))
  ! prefix = trim(pre)//trim(scan%object)//'_'//trim(itoa(scan%sid)) ! patchID_scanID

  ! if (myid == 0) call output_map_h5(trim(prefix), map)

  !end if

  !end do

  !if (myid == 0) write(*,*) 'Done'
  !do i = 1, nscan
  !call free_tod_type(tod(i))
  !end do
  !deallocate(tod)
  !call free_map_type(map)
  !call mpi_finalize(ierr)

contains



  subroutine output_tod(prefix, det, tod)!, alist)
    implicit none
    character(len=*), intent(in) :: prefix
    type(tod_type),   intent(in) :: tod
    !type(acceptlist), intent(in) :: alist
    integer(i4b),     intent(in) :: det

    character(len=512) :: filename
    character(len=4)   :: jtext
    integer(i4b) :: unit, i, j, l

    unit = getlun()
    j=1!do j = 1, tod%nfreq
       !if (alist%status(j,det) == 0) then
          call int2string(j,jtext)
          filename = trim(prefix) // 'tod.dat' !trim(prefix) // '_freq' // jtext // '_tod.dat'
          open(unit, file=trim(filename), recl=4096)
          do i = 1, tod%nsamp-100000
             write(unit, fmt='(f16.8)', advance='no') tod%t(i)
             write(unit, fmt='(2f24.8)', advance='no') tod%d(i,j,:,det)!, tod%d_raw(i,j,det)
             write(unit,*)
          end do
          close(unit)
       !end if
       !end do
       
  end subroutine output_tod

     
end program tod2comap
