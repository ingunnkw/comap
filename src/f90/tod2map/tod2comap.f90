program tod2comap
  use quiet_ces_mod
  use quiet_assembly_mod
  use comap_lx_mod
  use quiet_module_mod
  use quiet_utils
  use quiet_mpi_utils
  use quiet_fileutils
  use quiet_task_mod
  use quiet_acceptlist_mod
  use quiet_filter_mod
  use quiet_target_mod
  use quiet_mask_mod
  use tod2comap_utils
  use tod2comap_mapmaker
  use tod2comap_cl_mod
  implicit none
  character(len=512)   :: parfile, acceptfile, output_dir, target_name, arg
  character(len=512)   :: jackknives, string, fft3_magic_file, task_method
  character(len=512)   :: tot_dir, other_dir, jackinfofile, outparamfile
  character(len=512)   :: error_dir, task_file, monitor_file_name, ces_prefix, knife_dir
  character(len=512)   :: tot_prefix, scan_task_file, component_str, benchfile
  character(len=512)   :: checkfile, map_filetype
  character(len=8)     :: checkpointing
  integer(i4b)         :: root, nside, mask_degrade, i, j, k, l, m, n
  integer(i4b)         :: coord_out, ierr, nexp, nmap, debug, n_jk
  integer(i4b)         :: subind, subind_init, subfrom, subto, serial_steps
  real(dp)             :: comp_threshold
  logical(lgt)         :: fullsky, sparse_maps, status_output, reduced = .false.
  logical(lgt)         :: approximate_hor, output_pseudo_cls, output_nobs
  logical(lgt)         :: serial
  type(common_info)    :: info
  type(acceptlist)     :: alist
  type(task_list)      :: tasks
  type(quiet_target)   :: target
  type(quiet_ces_info) :: ces
  type(lx_struct)      :: l3data
  type(quiet_target),    dimension(:),     allocatable :: targets_all, targets
  type(swiss_knife),     dimension(:),     allocatable :: knife_defs
  integer(i4b),          dimension(:,:,:), allocatable :: knife_res
  integer(i4b),          dimension(:,:),   allocatable :: pixlist
  integer(i4b),          dimension(:),     allocatable :: cnums, map2mask
  real(sp),              dimension(:,:),   allocatable :: nobs
  real(dp),              dimension(:,:),   allocatable :: binmap
  real(sp),              dimension(:,:),   allocatable :: cstat
  real(sp),              dimension(:,:,:), allocatable :: dstat
  real(dp),              dimension(:,:,:), allocatable :: binmaps1, binmaps2
  real(dp),              dimension(:,:),   allocatable :: weight1, weight2
  logical(lgt),          dimension(:),     allocatable :: done
  type(shared_assembly), dimension(:),     allocatable :: groups
  type(mapdata),         dimension(:,:),   allocatable :: data
  type(mapdata)                                        :: data_ass
  logical(lgt), dimension(0:otype_num-1,0:level_num-1) :: needed_otypes
  logical(lgt), dimension(0:itype_num-1,0:level_num-1) :: needed_itypes
  logical(lgt), dimension(0:itype_num-1,0:level_num-1) :: raw_itypes
  type(jk_struct), dimension(100)                      :: jacks

  integer(i4b) :: bench_l3scan, bench_l3read, bench_ces, bench_ass, bench_aget, bench_aproc, bench_out, bench_task

  call setup_mpi
  call check_invocation
  call read_parameters
  call setup_mpi_groups
  call setup_paths
  call initialize_modules
  call setup_components
  call setup_bench

  call init_mon(info, monitor_file_name);         call update_mon(info, "init_mon")
  call init_target(target, alist);                call update_mon(info, "init_target")
  call filter_object(target, target_name);        call update_mon(info, "filter_object")
  call filter_comp(target, info, comp_threshold); call update_mon(info, "filter_comp")
  call collect_ces_info(target, cstat, &
    & dstat, map2mask);                           call update_mon(info, "collect_ces_info")
  call list_pixels(map2mask, pixlist)
  call jackknife(target, jackknives, cstat, dstat, targets_all, knife_defs, knife_res)
  call print_knives(knife_defs, knife_res, cid_list)
  deallocate(knife_res)
  call assert(size(targets_all) > 0, "No accepted ceses in any targets!")

  if(info%myid==0 .and. serial) then
     write(*,*) "Running in serial mode to reduce maximum memory usage when making split maps."
     write(*,*) "Note that this mode is not compatible with nulltest machinery."
  end if

  serial_steps = 1
  if(serial) serial_steps = size(targets_all)
  subind_init = 1
  if(index(checkpointing,"R")>0) call read_checkpoint_prepare(checkfile)
  do subind = subind_init, serial_steps
     ! Extract our current subset of targets
     subfrom = (subind-1)*size(targets_all)/serial_steps+1
     subto   = (subind-0)*size(targets_all)/serial_steps+0

     allocate(targets(subto-subfrom+1))
     do i = 1, subto-subfrom+1
        call copy_target(targets_all(subfrom+i-1),targets(i))
     end do
     if(info%myid == 0 .and. serial) write(*,'(a,i2,a1,i2)') "Processing serial group ", subind, "/", serial_steps

     call get_accepted_ceses_multi(targets%alist, cnums); call update_mon(info, "get_accepted_ceses")

     if (output_nobs) call collect_nobs(target, target_name, nobs)

     needed_itypes = setup_itypes(needed_otypes)
     raw_itypes    = setup_itypes(needed_otypes, raw=.true.)
     nexp          = maxval(map2mask)
     nmap          = 1+get_num_sim()
     info%nside    = nside
     info%coord    = coord_out
     info%target   = target_name
     tot_prefix    = trim(tot_dir) // "/totmaps/" // trim(target_name)

     ! Print out useful information
     if (info%myid == 0) call print_info(info, nexp, nmap, size(targets))

     call dmem("Before main memory alloc")
     allocate(data(0:level_num-1, 1:size(targets)))
     call setup_mapdata_multi(data(level_tot,:),nexp,nmap,info,needed_itypes(:,level_tot))
     call update_mon(info, "setup_mapdata_tot")
     call setup_mapdata_multi(data(level_ces,:),nexp,nmap,info,needed_itypes(:,level_ces))
     call update_mon(info, "setup_mapdata_ces")
     ! To save memory, delegate assembly data to total if possible
!     write(*,*) "A", info%ncomp, nmap, nexp
     if(size(data,2) /= 1) then
        call setup_mapdata(data_ass, nexp, nmap, info, needed_itypes(:,level_ass))
     else
        call setup_mapdata(data_ass, nexp, nmap, info, needed_itypes(:,level_ass), data2=data(level_tot,1), &
             & delegate=.not. any(raw_itypes(:,[level_ces,level_ass]),2))
     end if
     call update_mon(info, "setup_mapdata_assembly")
     call dmem("Main memory allocation")

     allocate(done(size(cnums)))
     done = .false.

     call mkdirs(task_file, .true.)
     call init_task_list_hack(tasks, task_file, size(cnums), info%comm_group, task_method)
     if(index(checkpointing,"R")>0) call read_checkpoint(checkfile, data(level_tot,:), done, reduced)

     ! Loop through each scan
     do while(get_next_task_hack(tasks, i, info%comm_cov))
        if(done(i)) cycle ! Could already be done due to checkpointing
        call dmem("ces start " // trim(itoa(i)))
        call bench_start(info%bench, bench_ces)
        call update_mon(info, "task_loop_start")

        call get_ces_info(cnums(i), ces)
        write(*,fmt="(i3,a12,a,i5,a)") info%myid, " processing", " ces ", ces%cid, " (" // &
         & trim(itoa(i)) // "/" // trim(itoa(size(cnums))) // ")"
        info%cid = ces%cid
        ces_prefix = trim(other_dir) // "/" // trim(ces%object)

        call bench_start(info%bench, bench_l3read)
        call read_l3_file(ces%l3file, l3data)
        call bench_stop (info%bench, bench_l3read)

        call update_mon(info, "task_loop_l3data")
        call clear_mapdata(data(level_ces,:))
        ! Loop over assemblies
        do j = 1, num_assemblies
           if(all(assemblies(j)%tquv_frac(info%comps) < comp_threshold)) cycle
           call bench_start(info%bench, bench_ass)
           call split_assembly(assemblies(j), targets%alist, ces%cid, groups)
           do k = 1, size(groups)
              call dmem("ass group " // trim(itoa(j)) // " " // trim(itoa(k)))
              call clear_mapdata(data_ass)
              call bench_start(info%bench, bench_aget)
              call get_assembly_data(groups(k)%assembly, info, l3data)
              call bench_stop(info%bench, bench_aget)
              call bench_start(info%bench, bench_aproc)

              call process_assembly_data(groups(k)%assembly, map2mask, info, data_ass)
              call bench_stop(info%bench, bench_aproc)
              do l = 1, size(groups(k)%members)
                 call add_mapdata(data_ass, data(level_tot,groups(k)%members(l)))
                 call add_mapdata(data_ass, data(level_ces,groups(k)%members(l)))
                 call output_as_needed_single(ces_prefix, data_ass, &
                   & needed_otypes(:,level_ass), info, set=groups(k)%members(l)+subfrom-1, &
                   & ces=ces%cid, ass=j, sparse=sparse_maps)
              end do
           end do
           call free_shared_assemblies(groups)
           call update_mon(info, "loop assemblies")
           call bench_stop(info%bench, bench_ass)
        end do
        ! SKN: I don't get this. Why don't we just loop over size(data,2) here?
        ! Like this:
        do j = 1, size(data,2)
           call output_as_needed_single(ces_prefix, data(level_ces,j), &
                & needed_otypes(:,level_ces), info, set=j+subfrom-1, ces=ces%cid, sparse=sparse_maps, silent=.true.)
        end do
        !call output_as_needed_single(ces_prefix, data(level_ces,1), &
        !     & needed_otypes(:,level_ces), info, set=1, ces=ces%cid, sparse=sparse_maps, silent=.true.)
        !if (size(data,2) > 1) then
        !   call output_as_needed_single(ces_prefix, data(level_ces,2), &
        !        & needed_otypes(:,level_ces), info, set=2, ces=ces%cid, sparse=sparse_maps, silent=.true.)
        !end if
        !call output_as_needed(ces_prefix, data(level_ces,:), &
        !     & needed_otypes(:,level_ces), info, sparse=sparse_maps)
        call free_lx_struct(l3data)
        call bench_stop(info%bench, bench_ces)
        done(i) = .true.
        if(index(checkpointing,"C")>0) then
           call mpi_barrier(info%comm_cov, ierr) ! Must have whole cov before writing
           call write_checkpoint(checkfile, data(level_tot,:), done, .false.)
        end if
     end do
     call free_mapdata(data_ass)
     call free_mapdata(data(level_ces,:))
     ! Final output stuff here
     if(.not. reduced) call merge_data(data(level_tot,:))
     if(index(checkpointing,"E")>0) then
        call mpi_barrier(info%comm_cov, ierr)
        call write_checkpoint(checkfile, data(level_tot,:), done, .true.)
     end if

     if(info%myid_group == root) call output_as_needed(tot_prefix, data(level_tot,:), &
          & needed_otypes(:,level_tot), info, sparse=sparse_maps, simstoo=.true., set_off=subfrom-1)
     call output_benchmarks(info, benchfile)
     ! Output jackknife power spectra and PTEs
     if(output_pseudo_cls) then
        i    = 1
        n_jk = 0
        do while (i < size(targets))
           if (targets(i)%n == 2) then
              n_jk = n_jk+1
              call allocate_jk_struct(jacks(n_jk), targets(i)%object, targets(i)%knife%full_name, &
                   & size(data(level_tot,i)%rhs,3))
              if (info%myid==0) then
                 write(*,*) 'Compute jackknife for ', trim(targets(i)%object), ', ', trim(targets(i)%knife%full_name)
              end if
              call mpi_bcast(data(level_tot,i)%rhs, size(data(level_tot,i)%rhs), &
               & mpi_rhs_type, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(data(level_tot,i)%div, size(data(level_tot,i)%div), &
               & mpi_rhs_type, 0, MPI_COMM_WORLD, ierr)
              call compute_binned_maps(info, data(level_tot,i),   binmaps1, weight1)

              call mpi_bcast(data(level_tot,i+1)%rhs, size(data(level_tot,i+1)%rhs), &
               & mpi_rhs_type, 0, MPI_COMM_WORLD, ierr)
              call mpi_bcast(data(level_tot,i+1)%div, size(data(level_tot,i+1)%div), &
               & mpi_rhs_type, 0, MPI_COMM_WORLD, ierr)
              call compute_binned_maps(info, data(level_tot,i+1), binmaps2, weight2)

              if (minval(weight1) /= 1.d30 .and. minval(weight2) /= 1.d30) then
                 weight1 = weight1 + weight2
                 weight1 = abs(minval(weight1) / weight1)
                 weight1 = weight1 / sum(weight1) *(12*nside**2/(4.d0*pi)) ! Normalize to unit integral on full sky
                 call output_pseudo_spectra_and_chisq(info, map2mask, output_dir, &
                  & targets(i)%object, targets(i)%knife%full_name, binmaps1, binmaps2, &
                  & weight1, jacks(n_jk)%pte, jacks(n_jk)%chisq, jacks(n_jk)%chi, jacks(n_jk)%chisq_max, &
                  & jacks(n_jk)%nbin)
              else
                 if (info%myid == 0) write(*,*) 'One of the two maps were unexposed. Rejecting jackknife.'
                 call deallocate_jk_struct(jacks(n_jk))
                 n_jk = n_jk-1
              end if
              deallocate(binmaps1, binmaps2, weight1, weight2)
           end if
           i = i + targets(i)%n
        end do
        if (info%myid==0) call output_jk_summary(output_dir, 1, 'full', jacks(1:n_jk))
        if (info%myid==0) call output_jk_summary(output_dir, 2, 'EE',   jacks(1:n_jk))
        if (info%myid==0) call output_jk_summary(output_dir, 3, 'BB',   jacks(1:n_jk))
     end if
     if(index(checkpointing,"D")>0) call rm(checkfile, noerr=.true.)
     call free_mapdata(data(level_tot,:))
     deallocate(data)
     deallocate(targets, cnums, done)
  end do
  call mpi_finalize(ierr)

  if (info%myid == root) then
     write(*,*)
     write(*,*) 'tod2map completed successfully'
  end if

contains

  subroutine check_invocation
    implicit none
    integer(i4b)         :: iargc, ierr
    if (iargc() < 1) then
       if (info%myid == root) then
          write(*,*) 'Usage: [mpirun [-n N]] tod2map paramfile'
       end if
       call mpi_finalize(ierr)
       stop
    else if (info%myid == root) then
       write(*,*) '***** QUIET map maker (svn ' // trim(itoa(svn_revision())) // ') ******'
       write(*,*)
    end if
  end subroutine

  subroutine read_parameters
    call getarg(1,parfile)
    call get_parameter(0, parfile, 'DEBUG',                 par_int=debug)
    call set_debug_level(debug)
    call get_parameter(0, parfile, 'ACCEPTLIST',            par_string=acceptfile)
    call get_parameter(0, parfile, 'OUTPUT_DIR',            par_string=output_dir)
    call get_parameter(0, parfile, 'COMPONENTS',            par_string=component_str, &
     & desc="Comma-separated string of which of T, Q, U and V to include.")
    call get_parameter(0, parfile, 'COMPONENT_THRESHOLD',   par_dp=comp_threshold, &
     & desc="Skip assemblies where the chosen components make up less than this fraction.")
    call get_parameter(0, parfile, 'OUTPUT_PSEUDO_CLS',     par_lgt=output_pseudo_cls)
    call get_parameter(0, parfile, 'NSIDE_OUT',             par_int=nside)
    call get_parameter(0, parfile, 'TARGET_NAME',           par_string=target_name)
    call get_parameter(0, parfile, 'JACKKNIVES',            par_string=jackknives)
    call get_parameter(0, parfile, 'OUTPUT_TOTAL',          par_string=string)
    needed_otypes(:,level_tot) = parse_output_options(string)
    call get_parameter(0, parfile, 'OUTPUT_CES',            par_string=string)
    needed_otypes(:,level_ces) = parse_output_options(string)
    call get_parameter(0, parfile, 'OUTPUT_ASSEMBLY',       par_string=string)
    needed_otypes(:,level_ass) = parse_output_options(string)
    call get_parameter(0, parfile, 'OUTPUT_NOBS',           par_lgt=output_nobs)
    call get_parameter(0, parfile, 'NUMPROCS_PER_COVAR',    par_int=info%nproc_cov)
    call get_parameter(0, parfile, 'COORDINATE_SYSTEM',     par_string=string)
    coord_out = parse_coord_name(string)
    call get_parameter(0, parfile, 'SPARSE_MAPS',           par_lgt=sparse_maps)
    call get_parameter(0, parfile, 'OUTPUT_STATUS',         par_lgt=status_output)
    call get_parameter(0, parfile, 'FFT3_MAGIC_NUMBER_FILE',par_string=fft3_magic_file)
    call get_parameter(0, parfile, 'APPROXIMATE_GROUND',    par_lgt=approximate_hor)
    call get_parameter(0, parfile, 'TASK_ITERATOR',         par_string=task_method, &
         & desc="Either 'smart' or 'dumb'. Use 'smart' for titan and 'dumb' for hercules. " // &
         & "Smart uses file locking, which hercules does not support.")
    call get_parameter(0, parfile, 'MASK_DEGRADE',          par_int=mask_degrade, &
     & desc="Number of steps to degrade map outside the mask. Negative to crop completely.")
    call get_parameter(0, parfile, 'FULLSKY',               par_lgt=fullsky, &
     & desc="Whether to use all pixels, or just the ones actually exposed.")
    call get_parameter(0, parfile, 'MAP_FILETYPE',          par_string=map_filetype, &
    & desc="File type for output maps. Can be fits or hdf.")
    call get_parameter(0, parfile, 'HIGHPASS_MIN_CAP',      par_dp=highpass_min_cap, &
    & desc="Highpass in scan freq units will be at least this.")
    call get_parameter(0, parfile, 'LOWPASS_MAX_CAP',      par_dp=lowpass_max_cap, &
    & desc="Lowpass frequency will be no more than this.")
    call get_parameter(0, parfile, 'AZORDER_MIN_CAP',       par_int=azorder_min_cap, &
    & desc="The azimuth filter order will be at least this.")
    call get_parameter(0, parfile, 'CHECKPOINTING',         par_string=checkpointing, &
    & desc="R to resume, C to save after each CES, E to save at end. Or any combination of these")
    call get_parameter(0, parfile, 'SERIAL',                par_lgt=serial, &
    & desc="If true, subsets in a jaccknife will be done one after another instead of in paralell. This saves memory.")
  end subroutine

  subroutine setup_paths
    ! Make all the output paths relative to the output directory:
    call mkdirs(output_dir, .false.)
    tot_dir        = trim(output_dir)
    other_dir      = trim(output_dir) // "/maps"
    jackinfofile   = trim(output_dir) // "/jackknives.txt"
    outparamfile   = trim(output_dir) // "/parameters.txt"
    error_dir      = trim(output_dir) // "/error"
    task_file      = trim(output_dir) // "/locks/tasks.dat"
    scan_task_file = trim(output_dir) // "/locks/scan_tasks.dat"
    knife_dir      = trim(output_dir) // "/knives"
    benchfile      = trim(output_dir) // "/bench.txt"
    checkfile      = trim(output_dir) // "/checkpoints/checkpoint" // trim(itoa(info%myid,3)) // ".hdf"
    if(status_output) then
       monitor_file_name = trim(output_dir) // "/status.txt"
    else
       monitor_file_name = ""
    end if
  end subroutine

  subroutine setup_mpi
    integer(i4b) :: ierr, myid, numprocs ! necessary line ??
    ! Initialize MPI environment
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, info%myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, info%nproc, ierr)
    call dset(id=info%myid)
    call print_host_mapping
    root = 0
  end subroutine

  subroutine setup_mpi_groups
    implicit none
    integer(i4b) :: icol, irow, ierr, numgroups
    ! Check that the number of processors match the number of chains and maps
    numgroups = info%nproc / info%nproc_cov
    if (numgroups*info%nproc_cov /= info%nproc) then
       if (info%myid == root) then
          write(*,*) ''
          write(*,*) 'ERROR: Number of processors is not a multiple of nproc_cov'
          write(*,*) '       numgroups         = ', numgroups
          write(*,*) '       num_proc_per_band = ', info%nproc_cov
          write(*,*) '       numprocs          = ', info%nproc, ' /= ', numgroups*info%nproc_cov
          write(*,*) ''
       end if
       call mpi_finalize(ierr)
       stop
    end if

    ! Split processors into work groups
    icol = info%myid/info%nproc_cov
    irow = mod(info%myid,info%nproc_cov)
    call mpi_comm_split(MPI_COMM_WORLD, icol, irow, info%comm_cov,   ierr) ! comm for each work group
    call mpi_comm_split(MPI_COMM_WORLD, irow, icol, info%comm_group, ierr) ! comm for different work groups
    call mpi_comm_rank(info%comm_cov,   info%myid_cov,    ierr)
    call mpi_comm_size(info%comm_cov,   info%nproc_cov,   ierr)
    call mpi_comm_rank(info%comm_group, info%myid_group,  ierr)
    call mpi_comm_size(info%comm_group, info%nproc_group, ierr)
  end subroutine

  subroutine initialize_modules
    integer(i4b)     :: base_seed, i, unit
    type(planck_rng) :: rng_handle
    real(dp)         :: mem, t1, t2

    unit = getlun()
    ! Set up auxilliary modules
    call cpu_time(t2)
    i = 14716
    call get_fft3_magic_number(i, fft3_magic_file)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid, ' time: ', t1-t2, ' mem: ', mem, ' -- assembly module'
    call initialize_quiet_assembly_mod(unit, parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid, ' time: ', t1-t2, ' mem: ', mem, ' -- pointing module'
    call initialize_quiet_pointing_mod(parfile, apparent_correction = .not. approximate_hor)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- mapmaker module'
    call initialize_map_maker(parfile, info)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- ces module'
    call initialize_ces_mod(parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- acceptlist'
    call initialize_accept_list(acceptfile, alist)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- target module'
    call initialize_target_mod(parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- jackknife module'
    call initialize_filter_mod(parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- filter module'
    call initialize_tod2map_cl_mod(parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- cl module'
    call initialize_patch_detect_mod(parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- patch module'
    call initialize_mask_mod(parfile)
    call cpu_time(t1)
    mem = get_mem_use()/1024d0**3
    if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- mask module'
    call setup_filetypes(map_filetype)

    !if (apply_i2qu) then
    !   if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- i2qu module'
    !   call initialize_i2qu(unit, parfile)
    !   call cpu_time(t1)
    !   mem = get_mem_use()/1024d0**3
    !end if
    !if (analyze_simulated_data) then 
    !   if (dtest(1)) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', info%myid,' time: ', t1-t2, ' mem: ', mem, ' -- tod sim module'
    !   call get_parameter(unit, parfile, 'BASE_SEED', par_int=base_seed)
    !   call initialize_random_seeds(MPI_COMM_WORLD, base_seed, rng_handle)
    !   base_seed = nint(rand_uni(rng_handle)*1.d7)
    !   call initialize_todsim_mod(unit, base_seed, parfile)
    !   call zigset(base_seed+1)
    !   call cpu_time(t1)
    !   mem = get_mem_use()/1024d0**3
    !end if
  end subroutine initialize_modules

  subroutine setup_components
    implicit none
    call toks2inds(component_str, comp_names, info%comps)
    info%ncomp = size(info%comps)
  end subroutine

  ! get_next_task was faulty in the cause of multiple procs per covar.
  ! These complications do not really belong in quiet task mod, so I am
  ! making a wrapper here instead.
  function get_next_task_coop(tasks, i, comm_coop) result(ok)
    implicit none
    type(task_list)  :: tasks
    integer(i4b)     :: i, coid, ierr, comm_coop
    logical(lgt)     :: ok
    call mpi_comm_rank(comm_coop, coid, ierr)
    if(coid == 0) ok = get_next_task(tasks, i)
    call mpi_bcast(i,  1, mpi_integer, 0, comm_coop, ierr)
    call mpi_bcast(ok, 1, mpi_logical, 0, comm_coop, ierr)
  end function

  ! Hercules does not support file locking properly, which breaks get_next_task.
  ! This wrapper makes it configurable whether to use this or a simple incremental
  ! method. This hack interprets negative file descriptors as an indication that
  ! simple stepping should be used
  function get_next_task_hack(tasks, i, comm_coop) result(ok)
    implicit none
    type(task_list)  :: tasks
    integer(i4b)     :: i, ierr, nproc
    logical(lgt)     :: ok
    integer(i4b), optional :: comm_coop
    call bench_start(info%bench, bench_task)
    if(tasks%fd >= 0) then
       if(present(comm_coop)) then
          ok = get_next_task_coop(tasks, i, comm_coop)
       else
          ok = get_next_task(tasks, i)
       end if
    else
       ! -fd is the current progress
       call mpi_comm_size(tasks%comm, nproc, ierr)
       i  = -tasks%fd
       ok = i <= tasks%n
       tasks%fd = -(i+nproc)
    end if
    call bench_stop(info%bench, bench_task)
  end function

  subroutine init_task_list_hack(tasks, filename, n, comm, method)
    implicit none
    type(task_list)  :: tasks
    character(len=*) :: filename, method
    integer(i4b)     :: n, id, ierr, comm
    select case(method)
       case("smart")
          call init_task_list(tasks, filename, n, comm)
       case("dumb")
          call mpi_comm_rank(comm, id, ierr)
          tasks%comm = comm
          tasks%fd = -id-1
          tasks%n  = n
       case default
          write(*,*) "Unrecognized task method '" // trim(method) // "'!"
          stop
    end select
  end subroutine

  subroutine claim_tasks_hack(tasks, mask)
    implicit none
    type(task_list)  :: tasks
    logical(lgt)     :: mask(:)
    if(tasks%fd >= 0 .and. info%myid_cov == 0) call claim_tasks(tasks, mask)
  end subroutine

  subroutine collect_ces_info(target, cstat, dstat, map2mask)
    implicit none
    type(quiet_target)                          :: target
    type(quiet_ces_info)                        :: ces
    type(hdf_file)                              :: file
    type(task_list)                             :: tasklist
    logical(lgt)                                :: found
    integer(i4b)                                :: i,j,k,m,n,npix, fnside, steps, ierr, id, maskcomp
    integer(i4b)                                :: ndi, nces, nmod, ext(7)
    real(dp)                                    :: t1, t2
    real(sp),     dimension(:,:),   allocatable :: cstat, mycstat
    real(sp),     dimension(:,:,:), allocatable :: dstat, mydstat
    real(dp),     dimension(:,:),   allocatable :: mask
    integer(i4b), dimension(:),     allocatable :: map2mask, cnums, hits, myhits, pixels, mask_helper

    npix = 12*nside**2
    nces = get_num_ces()
    nmod = get_num_modules()
    ndi  = get_num_diodes()*nmod
    allocate(mycstat(nces,STAT_NUM),cstat(nces,STAT_NUM))
    allocate(mydstat(ndi, nces,NUM_DIODE_STATS),dstat (ndi, nces,NUM_DIODE_STATS))
    allocate(myhits(0:npix-1),hits(0:npix-1))
    mycstat = 0; mydstat = 0
    myhits  = 0

    call get_accepted_ceses(target%alist, cnums)
    call mkdirs(scan_task_file, .true.)
    call init_task_list_hack(tasklist, scan_task_file, size(cnums), mpi_comm_world, task_method)
    do while(get_next_task_hack(tasklist, i)) ! missing comm_coop ??
       call get_ces_info(cnums(i), ces)
       write(*,fmt="(i3,a12,a,i5,a)") info%myid, " scanning", " ces ", ces%cid, " (" // &
        & trim(itoa(i)) // "/" // trim(itoa(size(cnums))) // ")"
       call bench_start(info%bench, bench_l3scan)
       call open_hdf_file(ces%l3file, file, "r")
       call read_hdf(file, "stats",       mycstat(  cnums(i),:))
       call read_hdf(file, "diode_stats", mydstat (:,cnums(i),:))
       call get_size_hdf(file, "pixels", ext)
       allocate(pixels(ext(1)))
       call read_hdf(file, "pixels", pixels)
       call read_hdf(file, "nside", fnside)
       steps = nint(log(real(fnside/nside,dp))/log(2d0))*2
       do j = 1, size(pixels)
          k = pixels(j)/2**steps
          myhits(k) = myhits(k)+1
       end do
       deallocate(pixels)
       call close_hdf_file(file)
       call bench_stop(info%bench, bench_l3scan)
    end do
    deallocate(cnums)

    call mpi_allreduce(mycstat, cstat, size(cstat), mpi_real,    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(mydstat, dstat, size(dstat), mpi_real,    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(myhits,  hits,  npix,        mpi_integer, mpi_sum, mpi_comm_world, ierr)
    deallocate(mycstat, mydstat, myhits)

    ! Get the standard mask
    allocate(mask(0:npix-1,3))
    call get_mask(target%object, nside, NEST, mask, found)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    if(.not. found .and. id == 0) write(stderr,'(a)') "Warning: Could not find mask for " // trim(target%object)

    ! Note: We have been using the T layer of the mask to build map2mask, using the assumption that all the layers are equal.
    ! What other things do we use the mask files for? Do we need all 3 layers? 
    maskcomp=1
    if(found .and. (.not. any(mask(:,1) > 0))) then
       maskcomp=2
       if(id==0) write(stderr,'(a)') "Note: T layer of mask for " // trim(target%object) // &
            & " is empty. Using Q layer to build map2mask in collect_ces_info."
    end if

    allocate(map2mask(0:npix-1), mask_helper(0:npix/4**max(0,mask_degrade)-1))
    map2mask  = 0
    mask_helper = 0
    j = 0
    if(fullsky) hits = 1
    do i = 0, npix-1
       if(hits(i) == 0) then
          cycle
       elseif(mask(i,maskcomp) > 0) then
          j = j+1
          map2mask(i) = j
       else
          ! We are in the masked region, so either degrade or skip entierly
          if(mask_degrade < 0) cycle
          k = i/4**mask_degrade
          if(mask_helper(k) == 0) then
             j = j+1
             mask_helper(k) = j
          end if
          map2mask(i) = mask_helper(k)
       end if
    end do
    deallocate(hits, mask, mask_helper)
  end subroutine collect_ces_info

  subroutine collect_nobs(target, name, nobs)
    implicit none
    character(len=*)     :: name
    type(quiet_target)   :: target
    type(quiet_ces_info) :: ces
    type(hdf_file)       :: file
    type(task_list)      :: tasklist
    character(len=512)   :: ofname
    integer(i4b)         :: i, j, k, m, n, ierr, nsamp, nmod, pix, npix, ext(7)
    real(sp),     dimension(:,:),   allocatable :: nobs, my_nobs
    real(dp),     dimension(:,:,:), allocatable :: point
    integer(i4b), dimension(:),     allocatable :: cnums

    npix = 12*nside**2
    allocate(nobs(0:npix-1,3), my_nobs(0:npix-1,3))
    my_nobs = 0

    call get_accepted_ceses(target%alist, cnums)
    call mkdirs(scan_task_file, .true.)
    call init_task_list_hack(tasklist, scan_task_file, size(cnums), mpi_comm_world, task_method)
    do while(get_next_task_hack(tasklist, i)) ! missing comm_coop ??
       call get_ces_info(cnums(i), ces)
       write(*,fmt="(i3,a12,a,i5,a)") info%myid, " scanning", " ces ", ces%cid, " (" // &
        & trim(itoa(i)) // "/" // trim(itoa(size(cnums))) // ")"
       call open_hdf_file(ces%l3file, file, "r")
       call get_size_hdf(file, "point", ext)
       nsamp = ext(2); nmod = ext(3)
       allocate(point(3,nsamp,nmod))
       call read_hdf(file, "point", point)
       call close_hdf_file(file)
       do i = 1, nmod
          do j = 1, nsamp
             call ang2pix_nest(nside, point(2,j,i), point(1,j,i), pix)
             my_nobs(pix,:) = my_nobs(pix,:) + 1.d0
          end do
       end do
       deallocate(point)
    end do
    deallocate(cnums)

    call mpi_allreduce(my_nobs, nobs, size(nobs), mpi_real, mpi_sum, mpi_comm_world, ierr)
    deallocate(my_nobs)

    where(nobs==0) nobs = hpx_dbadval
    ofname = trim(output_dir)//'/nobs/' // trim(name) // '.' // map_filetype
    call mkdirs(ofname, .true.)
    if (info%myid==0) call write_map(nobs, NEST, ofname)
    call mpi_finalize(ierr)
    stop
  end subroutine

  subroutine setup_mapdata_multi(data_sets, n, nmap, info, itypes)
    implicit none
    type(mapdata) :: data_sets(:)
    type(common_info) :: info
    integer(i4b)  :: n, set, nmap
    logical(lgt)  :: itypes(0:)
    call update_mon(info, "setup_mapdata_int1")
    do set = 1, size(data_sets)
       call setup_mapdata(data_sets(set), n, nmap, info, itypes)
    end do
    call update_mon(info, "setup_mapdata_int2")
  end subroutine

  function setup_itypes(otypes, raw) result(itypes)
    implicit none
    logical(lgt), dimension(0:otype_num-1,0:level_num-1) :: otypes
    logical(lgt), dimension(0:itype_num-1,0:level_num-1) :: itypes
    logical(lgt), dimension(0:itype_num-1,0:otype_num-1) :: deps
    logical(lgt), optional :: raw
    logical(lgt) :: T = .true., F = .false., raw_
    raw_ = .false.; if(present(raw)) raw_ = raw
    ! The itypes' dependency on the otypes
    deps = transpose(reshape([&
     & T, T, F, F, F, F, F, F, F, F, &
     & T, F, T, T, F, F, F, T, T, T, &
     & F, F, F, T, T, F, F, F, T, T, &
     & F, F, F, F, F, T, F, F, F, F, &
     & F, F, F, F, F, T, F, F, F, F, &
     & F, F, F, F, F, F, T, F, F, F  &
     &], [ otype_num, itype_num ]))
    itypes([itype_cov, itype_rhs, itype_div, itype_dir, itype_ang, &
      & itype_gfit],:) = matmul(deps, otypes([otype_eqn, otype_cov, &
      & otype_rhs, otype_bin, otype_rms, otype_cross, otype_gfit, &
      & otype_rawrhs, otype_sigma, otype_rhsdiv],:))
    ! The itypes' dependency on each other
    if (output_pseudo_cls) itypes([itype_rhs,itype_div],level_tot) = .true.
    if(.not. raw_) itypes(:,level_ass) = itypes(:,level_ass) .or. itypes(:,level_ces) &
      & .or. itypes(:,level_tot)
  end function

  subroutine output_as_needed_single(prefix, data, outopts, info, set, ces, ass, sparse, silent, simstoo)
    implicit none
    character(len=*)       :: prefix
    type(mapdata)          :: data
    type(common_info)      :: info
    logical(lgt), dimension(0:otype_num-1) :: outopts
    logical(lgt), optional :: sparse, silent, simstoo
    integer(i4b)           :: set, i, j, k
!real(dp) :: sumvar !TMR
    integer(i4b), optional :: ces, ass

    real(dp),     dimension(:,:,:), allocatable :: rhs, map
    character(len=512) :: outname
    integer(i4b)       :: unit, nside, npix, which_rhs, ordering
    logical(lgt)       :: chatty, simstoo_
    type(hdf_file)     :: hfile

    npix      = size(map2mask)
    nside     = nint(sqrt(real(npix,dp)/12.d0))
    ordering  = NEST
    ! Output simulated data if there is any, otherwise normal output
    chatty    = .true.; if (present(silent)) chatty = .not. silent
    simstoo_  = .false.; if(present(simstoo)) simstoo_ = simstoo
    which_rhs = 1; if(simstoo_) which_rhs = 0

    if (associated(data%rhs)) then
       if (all(data%rhs == 0.d0)) return
    end if
    !if(mapdata_unexposed(data)) return
    call bench_start(info%bench, bench_out)

    if(associated(data%rhs)) call rhs2map(data, rhs, which_rhs, info%comps)
    if(outopts(otype_rhs) .and. info%myid_cov == root) then
       call assert(associated(data%rhs), "Outputting rhs, but rhs not computed!")
       outname = make_outname(prefix, info%comps, otype_rhs, set, ces, ass)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
! TMR
!!$do i=1,size(rhs,2)
!!$   sumvar = 0.0
!!$   do j=1,size(rhs,1)
!!$      do k=1,size(rhs,3)
!!$         if (rhs(j,i,k) .ne. -1.6375d30) sumvar = sumvar + rhs(j,i,k)
!!$      end do
!!$   end do
!!$   write(32,*) sumvar
!!$end do
!!$write(32,*) "------------------------"
!!$do i=1,size(rhs,3)
!!$   sumvar = 0.0
!!$   do j=1,size(rhs,1)
!!$      do k=1,size(rhs,2)
!!$         if (rhs(j,k,i) .ne. -1.6375d30) sumvar = sumvar + rhs(j,k,i)
!!$      end do
!!$   end do
!!$   write(32,*) sumvar
!!$end do
!!$write(32,*) "------------------------"
       call output_map_general(rhs, outname, nside, ordering, map2mask, sparse)
    end if
    if(outopts(otype_bin) .and. info%myid_cov == root) then
       call assert(associated(data%rhs), "Outputting bin, but rhs not computed!")
       call assert(associated(data%div), "Outputting bin, but div not computed!")
       outname = make_outname(prefix, info%comps, otype_bin, set, ces, ass)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
       call solve_rhs_div(data, map, which_rhs, info%comps)
       call output_map_general(map, outname, nside, ordering, map2mask, sparse)
       deallocate(map)
    end if
    if(outopts(otype_rms) .and. info%myid_cov == root) then
       call assert(associated(data%div), "Outputting rms, but div not computed!")
       outname = make_outname(prefix, info%comps, otype_rms, set, ces, ass)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
       call div2rms(data, map, info%comps)
       call output_map_general(map, outname, nside, ordering, map2mask, sparse)
       deallocate(map)
    end if
    if(outopts(otype_sigma) .and. info%myid_cov == root) then
       call assert(associated(data%rhs), "Outputting sigma, but rhs not computed!")
       call assert(associated(data%div), "Outputting sigma, but div not computed!")
       outname = make_outname(prefix, info%comps, otype_sigma, set, ces, ass)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
       call calc_sigma(data, map, which_rhs)
       call output_map_general(map, outname, nside, ordering, map2mask, sparse)
       deallocate(map)
    end if
    if(outopts(otype_cov)) then ! Really the inverse cov
       if(info%myid_cov == root) then
          call assert(associated(data%cov), "Outputting cov, but cov not computed!")
          outname = make_outname(prefix, info%comps, otype_cov, set, ces, ass)
          call mkdirs(outname, .true.)
          if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
          unit = getlun()
          open(unit,file=outname,form="unformatted")
          write(unit) int(data%n*info%ncomp,i4b)
          write(unit) int(2,i4b)
          write(unit) int(info%ncomp,i4b)
       end if
       call output_cov(unit, data) ! root gathers col by col and writes
       if(info%myid_cov == root) then
          write(unit) .true.
          close(unit)
       end if
    end if
    if(outopts(otype_eqn)) then
       if(info%myid_cov == root) then
          call assert(associated(data%cov), "Outputting eqn, but cov not computed!")
          call assert(associated(data%rhs), "Outputting eqn, but rhs not computed!")
          outname = make_outname(prefix, info%comps, otype_eqn, set, ces, ass)
          call mkdirs(outname, .true.)
          if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
          unit = getlun()
          open(unit,file=outname,form="unformatted")
          write(unit) int(info%ncomp*data%n,i4b)
          write(unit) int(2,i4b)
          write(unit) int(info%ncomp,i4b)
       end if
!       write(*,*) 'rhs sum = ', sum(abs(data%rhs(:,:,1)))
!       write(*,*) 'cov sum = ', sum(abs(data%cov))
       call output_cov(unit, data)
       if(info%myid_cov == root) then
          write(unit) int(data%nmap,i4b)
          do i = 1, data%nmap
             write(unit) real(data%rhs(:,:,i),dp) ! FIXME: do we want sp output?
          end do
          write(unit) npix
          write(unit) map2mask
          close(unit)
       end if
    end if
    if(outopts(otype_cross) .and. info%myid_cov == root) then
       call assert(associated(data%dir), "Outputting cross, but dir not computed!")
       call assert(associated(data%ang), "Outputting cross, but ang not computed!")
       call dirang2cross(data, map)
       outname = make_outname(prefix, info%comps, otype_cross, set, ces, ass)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
       call output_map_general(map, outname, nside, ordering, map2mask, sparse)
       deallocate(map)
    end if
    if(outopts(otype_rawrhs) .and. info%myid_cov == root) then
       call assert(associated(data%rhs), "Outputting rawrhs, but rhs not computed!")
       outname = make_outname(prefix, info%comps, otype_rawrhs, set, ces, ass)
       call mkdirs(outname, .true.)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
       unit = getlun()
       open(unit,file=outname,form="unformatted")
       write(unit) int(data%n*info%ncomp,i4b)
       write(unit) int(2,i4b)
       write(unit) int(info%ncomp,i4b)
       write(unit) int(data%nmap,i4b)
       do i = 1, data%nmap
          write(unit) data%rhs(:,:,i)
       end do
       close(unit)
    end if
    if(outopts(otype_rhsdiv) .and. info%myid_cov == root) then
       call assert(associated(data%rhs), "Outputting rhsdiv, but rhs not computed!")
       call assert(associated(data%div), "Outputting rhsdiv, but div not computed!")
       outname = make_outname(prefix, info%comps, otype_rhsdiv, set, ces, ass)
       call mkdirs(outname, .true.)
       if(dtest(1) .and. chatty) write(*,"(a)") "Writing: " // trim(outname)
       call open_hdf_file(outname, hfile, "w")
       call write_hdf(hfile, "npix",    int(data%n,i4b))
       call write_hdf(hfile, "nside",   info%nside)
       call write_hdf(hfile, "order",   nest)
       call write_hdf(hfile, "ncomp",   info%ncomp)
       call write_hdf(hfile, "comps",   info%comps)
       call write_hdf(hfile, "rhs",     data%rhs)
       call write_hdf(hfile, "div",     data%div)
       call write_hdf(hfile, "pixlist", pixlist)
       call close_hdf_file(hfile)
    end if

    if(allocated(rhs))      deallocate(rhs)
    call bench_stop(info%bench, bench_out)
  end subroutine

  subroutine output_as_needed(prefix, data, outopts, info, ces, ass, sparse, silent, simstoo, set_off)
    implicit none
    character(len=*)       :: prefix
    type(mapdata)          :: data(:)
    type(common_info)      :: info
    logical(lgt)           :: outopts(0:)
    logical(lgt), optional :: sparse, silent, simstoo
    integer(i4b)           :: set, sub, set_off_
    integer(i4b), optional :: ces, ass, set_off
    set_off_ = 0; if(present(set_off)) set_off_ = set_off
    do set = 1, size(data)
       call output_as_needed_single(prefix, data(set), outopts, &
         & info, set+set_off_, ces, ass, sparse, silent, simstoo)
    end do
  end subroutine

  ! Reforms rhs into a sparse map (*not* a full map). As with other
  ! maps here, it has 3 indices, where the first is used for temperature
  ! and the others for polarization.
  subroutine rhs2map(data, rhsmap, which, comps)
    implicit none
    type(mapdata)                           :: data
    real(dp), dimension(:,:,:), allocatable :: rhsmap
    integer(i4b)                            :: which, comps(:), r(2)
    r = [1, int(data%nmap)]; if(which > 0) r = [which,which]
    allocate(rhsmap(data%n,3,r(2)-r(1)+1))
    rhsmap = hpx_dbadval
    rhsmap(:,comps,:) = real(data%rhs(:,:,r(1):r(2)),dp)
  end subroutine

  subroutine solve_rhs_div(data, binmap, which, comps)
    implicit none
    type(mapdata)                          :: data
    real(dp), dimension(:,:,:), allocatable  :: binmap
    integer(i4b)                           :: which, comps(:), i, status, r(2), j
    real(dp)                               :: v(size(comps))
    r = [1, int(data%nmap)]; if(which > 0) r = [which,which]
    allocate(binmap(data%n,3,r(2)-r(1)+1))
    binmap = hpx_dbadval
    do j = r(1), r(2)
       do i = 1, data%n
          if(data%div(i,1,1) == 0) cycle
          call solve_system_eiglim(real(data%div(i,:,:),dp), real(data%rhs(i,:,j),dp), v, 1d-12,status)
          if(status /= 0) cycle
          binmap(i,comps,j-r(1)+1) = v
       end do
    end do
  end subroutine

  subroutine calc_sigma(data, sigmamap, which)
    implicit none
    type(mapdata)                          :: data
    real(dp), dimension(:,:,:), allocatable  :: sigmamap
    integer(i4b)                           :: which, i, j, n, status, r(2)
    real(dp)                               :: v(size(data%rhs,2))
    real(dp)                               :: A(size(data%rhs,2),size(data%rhs,2))
    real(dp)                               :: intelbug(size(data%rhs,2))
    real(dp)                               :: x, mean, dev
    r = [1, int(data%nmap)]; if(which > 0) r = [which,which]
    allocate(sigmamap(data%n,3,r(2)-r(1)+1))
    sigmamap = hpx_dbadval
    n    = size(v)
    mean = 1-2/(9d0*n)
    dev  = sqrt(2/(9d0*n))
    do i = 1, data%n
       A = data%div(i,:,:)
       if(any(A /= A)) cycle ! NaN
       intelbug = get_diag(A)
       if(any(intelbug == 0)) cycle ! unexposed
       status = 0
       call invert_matrix(A); if(status /= 0) cycle
       do j = r(1), r(2)
          v = data%rhs(i,:,j)
          x = dot_product(v,matmul(A,v))
          sigmamap(i,1,j-r(1)+1) = ((x/n)**(1d0/3)-mean)/dev
       end do
    end do
  end subroutine

  subroutine div2rms(data, rmsmap, comps)
    implicit none
    type(mapdata)                         :: data
    real(dp), dimension(:,:,:), allocatable :: rmsmap
    integer(i4b)                          :: i, status, comps(:)
    real(dp)                              :: A(size(comps),size(comps))
    allocate(rmsmap(data%n,3,1))
    rmsmap = hpx_dbadval
    do i = 1, data%n
       A = data%div(i,:,:)
       if(any(A /= A)) cycle ! NaN
       status = 0
       call invert_sym(A, status); if(status /= 0) cycle
       call cholesky  (A, status); if(status /= 0) cycle
       rmsmap(i,comps,1) = get_diag(A)
    end do
  end subroutine

  subroutine dirang2cross(data, map)
    implicit none
    type(mapdata)                          :: data
    real(dp), dimension(:,:,:), allocatable  :: map
    integer(i4b)                           :: i, status
    allocate(map(data%n,3,1))
    map = hpx_dbadval
    do i = 1, data%n
       if(data%dir(i,3) /= 0) map(i,1,1) = sum((data%dir(i,1:2)/data%dir(i,3))**2)
       if(data%ang(i,3) /= 0) map(i,2,1) = sum((data%ang(i,1:2)/data%ang(i,3))**2)
       if(data%ang(i,3) /= 0) map(i,3,1) = data%ang(i,3)
    end do
  end subroutine

  ! Writes the distributed covariance matrix cov to unit, which must
  ! be an already opened file on root. cov must be distributed in the
  ! mapdata manner.
  ! This becamoe less efficient when we went from (:,:) to (:,:,:,:)
  ! due to extra waiting.
  subroutine output_cov(unit, data)
    implicit none
    integer(i4b)  :: unit, i, j, p, q, owner, stat(5), ierr, col1, col2, n, m, l
    integer(i4b)  :: s1, s2, p1, p2
    type(mapdata) :: data
    real(cov_type), dimension(:), allocatable :: gcol, lcol

    if(dtest(1)) write(*,*) 'Processor ', info%myid_cov, ' ready to output covariance matrix'
    n = data%n*size(data%rhs,2)
    allocate(gcol(n+1), lcol(n+1))
    do j = 1, n
       lcol = 0
       s1 = (j-1)/data%n+1
       p1 = modulo(j-1,data%n)+1
       ! 2. Build up our local contribution to this column
       do i = 1, n ! F part
          s2 = (i-1)/data%n+1
          p2 = modulo(i-1,data%n)+1
          if(p1 >= p2) then ! At F-owned part in local.
             if(i <= j .and. p1 >= data%col_start .and. p1 <= data%col_stop) &
              & lcol(i) = data%cov(p2,s2,p1,s1) ! F part
             if(i >= j .and. p2 >= data%col_start .and. p2 <= data%col_stop) &
              & lcol(i+1) = data%cov(p1+1,s1,p2,s2) ! FF part
          else ! At FF-owned part. Must mirror if we want F-part.
             if(i <= j .and. p2 >= data%col_start .and. p2 <= data%col_stop) &
              & lcol(i) = data%cov(p1,s1,p2,s2) ! F part
             if(i >= j .and. p1 >= data%col_start .and. p1 <= data%col_stop) &
              & lcol(i+1) = data%cov(p2+1,s2,p1,s1) ! FF part
          end if
       end do
       ! 3. Mpi reduce these to get the full column, and output.
       call mpi_reduce(lcol, gcol, n+1, mpi_cov_type, MPI_SUM, &
        & root, info%comm_cov, ierr)
       ! FIXME: Change gcol type if you want output in cov_type
       if(info%myid_cov == root) write(unit) real(gcol,dp)
    end do
    deallocate(lcol, gcol)
    if(dtest(1)) write(*,*) 'Processor ', info%myid_cov, ' finished outputting covariance matrix'
  end subroutine

  ! Merge data in intergroup-direction.
  subroutine merge_data_single(data)
    implicit none
    type(mapdata) :: data
    integer(i4b)  :: i, ncomp
    real(rhs_type), dimension(:,:),   allocatable :: ibuf1, obuf1
    real(rhs_type), dimension(:,:,:), allocatable :: ibuf2, obuf2
    real(cov_type), dimension(:,:),   allocatable :: ibuf3, obuf3
    ncomp = size(data%rhs,2)
    call update_mon(info, "merge_data_A")
    if(associated(data%rhs) .and. info%myid_cov == root) then
       allocate(ibuf1(data%n,ncomp), obuf1(data%n,ncomp))
       do i = 1, data%nmap
          ibuf1 = data%rhs(:,:,i)
          call mpi_reduce(ibuf1, obuf1, size(ibuf1), mpi_rhs_type, &
            & mpi_sum, root, info%comm_group, ierr)
          data%rhs(:,:,i) = obuf1
       end do
       call update_mon(info, "merge_data_B")
       deallocate(ibuf1,obuf1)
    end if
    if(associated(data%div) .and. info%myid_cov == root) then
       allocate(ibuf2(data%n,ncomp,ncomp))
       allocate(obuf2(data%n,ncomp,ncomp))
       ibuf2 = data%div
       call mpi_reduce(ibuf2, obuf2, size(ibuf2), mpi_rhs_type, &
         & mpi_sum, root, info%comm_group, ierr)
       data%div = obuf2
       call update_mon(info, "merge_data_C")
       deallocate(ibuf2,obuf2)
    end if
    if(associated(data%cov)) then
       allocate(ibuf3(data%n+1,ncomp),obuf3(data%n+1,ncomp))
       do j = 1, ncomp
          do i = lbound(data%cov,3), ubound(data%cov,3)
             ibuf3 = data%cov(:,:,i,j)
             call mpi_reduce(ibuf3, obuf3, size(ibuf3), mpi_cov_type, &
               & mpi_sum, root, info%comm_group, ierr)
             data%cov(:,:,i,j) = obuf3
          end do
       end do
       call update_mon(info, "merge_data_D")
       deallocate(ibuf3,obuf3)
    end if
    if(associated(data%ang) .and. info%myid_cov == root) then
       allocate(ibuf1(size(data%ang,1),size(data%ang,2)), &
        & obuf1(size(data%ang,1),size(data%ang,2)))
       ibuf1 = data%ang
       call mpi_reduce(ibuf1, obuf1, size(ibuf1), mpi_rhs_type, &
         & mpi_sum, root, info%comm_group, ierr)
       data%ang = obuf1
       call update_mon(info, "merge_data_E")
       deallocate(ibuf1,obuf1)
    end if
    if(associated(data%dir) .and. info%myid_cov == root) then
       allocate(ibuf1(size(data%dir,1),size(data%dir,2)), &
        & obuf1(size(data%dir,1),size(data%dir,2)))
       ibuf1 = data%dir
       call mpi_reduce(ibuf1, obuf1, size(ibuf1), mpi_rhs_type, &
         & mpi_sum, root, info%comm_group, ierr)
       data%dir = obuf1
       call update_mon(info, "merge_data_F")
       deallocate(ibuf1,obuf1)
    end if
  end subroutine

  subroutine merge_data(data)
    implicit none
    type(mapdata) :: data(:)
    integer(i4b)  :: set
    do set = 1, size(data)
       call merge_data_single(data(set))
    end do
  end subroutine

  function make_outname(prefix, comps, type, set, ces, ass) result(outfile)
    implicit none
    character(len=*),         intent(in) :: prefix
    integer(i4b),             intent(in) :: comps(:)
    integer(i4b), optional,   intent(in) :: set, ces, ass
    character(len=512)                   :: outfile
    integer(i4b)                         :: type, i

    character(len=32) :: cesstr, assstr, tstr, setstr

    setstr = ""
    cesstr = ""
    assstr = ""
    tstr   = ""
    if(present(set)) setstr = "_set" // trim(itoa(set,3))
    if(present(ces)) cesstr = "_ces" // trim(itoa(ces,4))
    if(present(ass)) assstr = "_ass" // trim(itoa(ass,3))
    do i = 1, size(comps)
       tstr = trim(tstr) // comp_names(comps(i))
    end do
    outfile = trim(prefix) // trim(setstr) // trim(cesstr) // &
       & trim(assstr) // "_" // trim(tstr) // "_" // trim(otypestrs(type)) // "." // &
       & trim(ofiletypes(type))
  end function

  subroutine compute_binned_maps(info, data, binmaps, weight)
    implicit none

    type(common_info) :: info
    type(mapdata)     :: data
    real(dp), allocatable, dimension(:,:,:) :: binmaps
    real(dp), allocatable, dimension(:,:) :: weight

    integer(i4b) :: i, j, pix, npix, nmap, ncomp
    real(dp), allocatable, dimension(:,:,:) :: binmap, wtmp

    nmap   = data%nmap
    ncomp  = info%ncomp

    allocate(binmaps(data%n,3,nmap))
    do i = 1+info%myid, nmap, info%nproc
       call solve_rhs_div(data, binmap, i, info%comps)
       binmaps(:,:,i) = binmap(:,:,1)
       deallocate(binmap)
    end do

    call div2rms(data, wtmp, info%comps)
    allocate(weight(size(wtmp,1),size(wtmp,2)))
    weight = wtmp(:,:,1)
    deallocate(wtmp)
    do i = 1, size(weight,1)
       if (sum(data%div(i,:,:)) == 0.d0) then
          weight(i,:) = 1.d30
       else if (any(info%comps == Q) .and. any(info%comps == U)) then
          weight(i,T)   = 1.d30
          weight(i,[Q,U]) = weight(i,Q)**2 + weight(i,U)**2
       else if (any(info%comps == T)) then
          weight(i,T)   = weight(i,1)**2
          weight(i,[Q,U]) = 1.d30
       end if
    end do

  end subroutine compute_binned_maps

  subroutine output_jk_summary(dir, comp, tag, jacks)
    implicit none

    integer(i4b),                   intent(in) :: comp
    character(len=*),               intent(in) :: dir, tag
    type(jk_struct), dimension(1:), intent(in) :: jacks

    integer(i4b) :: i, j, k, l, m, unit, n, nmap, nbin
    real(dp)     :: mu, sigma
    logical(lgt) :: first
    character(len=16)  :: split
    character(len=512) :: filename
    real(dp), allocatable, dimension(:) :: chisq, chi, chisq_max

    n    = size(jacks)
    if(n < 1) return
    nmap = size(jacks(1)%pte(:,1))
    unit = getlun()
    filename = trim(dir) // '/cls/jk_summary_' // tag // '.dat'
    call mkdirs(filename, .true.)
    open(unit, file=trim(filename), recl=1024)
    write(unit,*) '******** TOD2MAP JACKKNIFE SUMMARY ***********'

    allocate(chisq(nmap), chi(nmap), chisq_max(nmap))
    do i = 1, size(patches)
       first = .true.
       chisq     = 0.d0
       chi       = 0.d0
       chisq_max = 0.d0
       nbin      = 0
       do j = 1, n
          if (trim(jacks(j)%object) == trim(patches(i)%name)) then
             if (first) then
                first = .false.
                write(unit,*)
                write(unit,*) 'Object = ', trim(jacks(j)%object)
             end if
             chisq = chisq + jacks(j)%chisq(:,comp)
             chi   = chi   + jacks(j)%chi(:,comp)
             nbin  = nbin  + jacks(j)%nbin(comp)
             do k = 1, nmap
                chisq_max(k) = max(chisq_max(k), jacks(j)%chisq_max(k,comp))
             end do
             split = ''
             split = trim(adjustl(jacks(j)%name))
             mu = mean(jacks(j)%chisq(2:,comp)/jacks(j)%nbin(comp))
             sigma = sqrt(variance(jacks(j)%chisq(2:,comp)/jacks(j)%nbin(comp))) 
             write(unit,fmt="(a,a,a,f6.4,a,f6.4,a,f7.4,a,f7.4,a,f6.4)") &
                  & ' ', split, ' val ', jacks(j)%chisq(1,comp)/jacks(j)%nbin(comp), &
                  & ' PTE ', jacks(j)%pte(1,comp), ' sigma ', (jacks(j)%chisq(1,comp)/jacks(j)%nbin(comp)-mu)/sigma, &
                  & ' sim ', mean(jacks(j)%chisq(2:,comp)/jacks(j)%nbin(comp)), ' +- ', &
                  & sqrt(variance(jacks(j)%chisq(2:,comp)/jacks(j)%nbin(comp)))
          end if
       end do
       if (.not. first) then
          chisq = chisq / nbin
          chi   = chi   / nbin
          mu = mean(chisq(2:))
          sigma = sqrt(variance(chisq(2:))) 
          write(unit,fmt="(a,f6.4,a,f6.4,a,f7.4,a,f7.4,a,f6.4)") &
               & ' Total chisq      red ', chisq(1), &
               & ' PTE ', count(chisq(2:) > chisq(1)) / real(nmap-1,dp), &
               & ' sigma ', (chisq(1)-mu)/sigma, &
               & ' sim ', mu, ' ± ', sigma
          mu = mean(chisq_max(2:))
          sigma = sqrt(variance(chisq_max(2:))) 
          write(unit,fmt="(a,f6.2,a,f6.4,a,f7.4,a,f7.4,a,f6.4)") &
               & ' Max chisq        val ', real(chisq_max(1),sp), &
               & ' PTE ', count(chisq_max(2:) > chisq_max(1)) / real(nmap-1,dp), &
               & ' sigma ', (chisq_max(1)-mu)/sigma, &
               & ' sim ', mu, ' ± ', sigma
          mu = mean(chi(2:))
          sigma = sqrt(variance(chi(2:))) 
          write(unit,fmt="(a,f6.4,a,f6.4,a,f7.4,a,f7.4,a,f6.4)") &
               & ' Mean chi shift   val ', chi(1), &
               & ' PTE ', count(chi(2:) > chi(1)) / real(nmap-1,dp), &
               & ' sigma ', (chi(1)-mu)/sigma, &
               & ' sim ', mu, ' ± ', sigma
       end if
    end do
    deallocate(chisq)
    close(unit)

  end subroutine output_jk_summary

  subroutine print_info(info, nexp, nmap, nset)
    implicit none

    type(common_info)    :: info
    integer(i4b)         :: nexp, nmap, nset, neff, ncov
    neff = nexp*info%ncomp
    ncov = count(needed_itypes(itype_cov,[level_tot,level_ces]))
    if(size(targets) /= 1 .or. any(raw_itypes(itype_cov,[level_ces,level_ass]))) then
       ncov = ncov + count(needed_itypes(itype_cov,[level_ass]))
    end if
    1 format(a,i9)
    2 format(a,f9.3)
    write(*,*)
    write(*,1) '   Number of exposed   pixels = ', nexp
    write(*,1) '   Number of effective pixels = ', neff
    write(*,1) '   Number of maps             = ', nmap
    write(*,1) '   Number of sets             = ', nset
    write(*,1) '   Number of covs             = ', ncov
    write(*,1) '   Tasks per cov              = ', info%nproc_cov
    write(*,2) '   Covmat GB per task         = ', real(neff,dp)**2*ncov/info%nproc_cov*4/1024**3
    write(*,*)
  end subroutine print_info

  subroutine print_knives(defs, groups, cids)
    implicit none
    type(swiss_knife), intent(in) :: defs(:)
    integer(i4b),      intent(in) :: groups(:,:,:)
    integer(i4b),      intent(in) :: cids(:)
    integer(i4b)                  :: i, j, k, n, unit, cnums
    character(len=512)            :: fname
    if (info%myid /= 0) return
    n = size(defs)
    unit = getlun()
    do i = 1, n
       fname = trim(knife_dir) // "/knife" // trim(itoa(i,3)) // "_" // trim(defs(i)%full_name) // ".txt"
       call mkdirs(fname,.true.)
       open(unit,file=fname)
       do j = 1, size(groups,2)
          if(.not. any(groups(:,j,i) /= 0)) cycle
          write(unit,'(i5)',advance="no") cids(j)
          do k = 1, size(groups,1)
             write(unit,'(i4)',advance="no") groups(k,j,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do
  end subroutine

  subroutine setup_bench
    implicit none
    bench_l3scan = bench_init(info%bench, "l3scan")
    bench_l3read = bench_init(info%bench, "l3read")
    bench_ces    = bench_init(info%bench, "ces")
    bench_ass    = bench_init(info%bench, "ass_loop")
    bench_aget   = bench_init(info%bench, "ass_get")
    bench_aproc  = bench_init(info%bench, "ass_proc")
    bench_out    = bench_init(info%bench, "output")
    bench_task   = bench_init(info%bench, "get_task")
  end subroutine

  !subroutine dump_bench
  !  implicit none
  !  integer(i4b) :: i
  !  write(*,'(i4)') info%myid
  !  do i = 1, info%bench%n
  !     write(*,'(a,f12.5)',advance="no") ' ', info%bench%dt(i)/info%bench%nsamp(i)
  !  end do
  !  write(*,*)
  !end subroutine

  subroutine output_benchmarks(info, benchfile)
    implicit none
    type(common_info), intent(in) :: info
    character(len=*),  intent(in) :: benchfile
    integer(i4b)                  :: i, j, n, ierr, unit
    real(dp),     dimension(:,:), allocatable :: rbuf, dt
    integer(i8b), dimension(:,:), allocatable :: ibuf, nsamp

    n = info%bench%n
    allocate(ibuf (n,info%nproc),rbuf(n,info%nproc))
    allocate(nsamp(n,info%nproc),dt  (n,info%nproc))
    ibuf = 0; rbuf = 0
    ! Sum up the OMP thread contributions
    ibuf(:,info%myid+1) = sum(info%bench%nsamp(:,1:n),1)
    rbuf(:,info%myid+1) = sum(info%bench%dt   (:,1:n),1)
    call mpi_reduce(ibuf, nsamp, size(ibuf), mpi_integer8, MPI_SUM, &
     & root, mpi_comm_world, ierr)
    call mpi_reduce(rbuf, dt, size(ibuf), mpi_double_precision, MPI_SUM, &
     & root, mpi_comm_world, ierr)
    if(info%myid == root) then
       unit = getlun()
       open(unit,file=benchfile)
       do i = 1, n
          write(unit,'(a16)',advance="no") info%bench%name(i)
          ! Total average
          write(unit,'(e15.7)',advance="no") sum(dt(i,:))/sum(nsamp(i,:))
          do j = 1, info%nproc
             write(unit,'(e15.7)',advance="no") dt(i,j)/nsamp(i,j)
          end do
          ! Total nsamp
          write(unit,'(e15.7)',advance="no") real(sum(nsamp(i,:)),dp)
          do j = 1, info%nproc
             write(unit,'(e15.7)',advance="no") real(nsamp(i,j),dp)
          end do
          write(unit,*)
       end do
       close(unit)
    end if
    deallocate(ibuf, rbuf, dt, nsamp)
  end subroutine

  subroutine read_checkpoint_prepare(fname)
    implicit none
    character(len=*), intent(in)    :: fname
    type(hdf_file)                  :: file
    logical(lgt)                    :: exist
    inquire(file=fname, exist=exist)
    if(exist) then
       call open_hdf_file(fname, file, "r")
       call read_hdf(file, "subind", subind_init)
       call close_hdf_file(file)
    end if
  end subroutine

  subroutine read_checkpoint(fname, data, done, reduced)
    implicit none
    character(len=*), intent(in)    :: fname
    type(mapdata),    intent(inout) :: data(:)
    logical(lgt),     intent(out)   :: done(:)
    type(hdf_file)                  :: file
    integer(i4b)                    :: npix, ncomp, nmap, i, idone(size(done)), ired
    character(len=32)               :: set
    logical(lgt)                    :: exist, reduced
    inquire(file=fname, exist=exist)
    if(exist) then
       call open_hdf_file(fname, file, "r")
       call read_hdf(file, "ncomp", ncomp)
       call read_hdf(file, "npix",  npix)
       call read_hdf(file, "nmap",  nmap)
       call assert(ncomp == info%ncomp,    "Inconsistent ncomp in checkpoint")
       call assert(npix  == data(1)%n,     "Inconsistent npix  in checkpoint")
       call assert(nmap  == data(1)%nmap,  "Inconsistent nmap  in checkpoint")
       call read_hdf(file, "done", idone)
       done = idone > 0
       call read_hdf(file, "reduced", ired)
       reduced = ired > 0
       if(.not. reduced .or. info%myid == root) then
          do i = 1, size(data)
             set = "set" // trim(itoa(i))
             if(associated(data(i)%rhs)) call read_hdf(file,trim(set)//"/rhs",data(i)%rhs)
             if(associated(data(i)%div)) call read_hdf(file,trim(set)//"/div",data(i)%div)
             if(associated(data(i)%cov)) call read_hdf(file,trim(set)//"/cov",data(i)%cov)
             if(associated(data(i)%ang)) call read_hdf(file,trim(set)//"/ang",data(i)%ang)
             if(associated(data(i)%dir)) call read_hdf(file,trim(set)//"/dir",data(i)%dir)
          end do
       end if
       call close_hdf_file(file)
       call claim_tasks_hack(tasks, done) ! Update task iterator to match done list
    end if
    call mpi_barrier(mpi_comm_world, i)
  end subroutine

  subroutine write_checkpoint(fname, data, done, reduced)
    implicit none
    character(len=*), intent(in)    :: fname
    type(mapdata),    intent(in)    :: data(:)
    logical(lgt),     intent(in)    :: done(:)
    character(len=512)              :: tmpname
    type(hdf_file)                  :: file
    logical(lgt)                    :: reduced
    integer(i4b)                    :: i, idone(size(done)), err, ired
    character(len=32)               :: set
    call mkdirs(fname, .true.)
    tmpname = trim(fname) // ".tmp"
    call open_hdf_file(tmpname, file, "w")
    call write_hdf(file, "subind",int(subind,i4b))
    call write_hdf(file, "ncomp", int(info%ncomp,i4b))
    call write_hdf(file, "npix",  int(data(1)%n,i4b))
    call write_hdf(file, "nmap",  int(data(1)%nmap,i4b))
    idone = 0; where(done) idone = 1
    call write_hdf(file, "done", idone)
    ired = 0; if(reduced) ired = 1
    call write_hdf(file, "reduced", ired)
    if(.not. reduced .or. info%myid == root) then
       do i = 1, size(data)
          set = "set" // trim(itoa(i))
          call create_hdf_group(file, trim(set))
          if(associated(data(i)%rhs)) call write_hdf(file,trim(set)//"/rhs",data(i)%rhs)
          if(associated(data(i)%div)) call write_hdf(file,trim(set)//"/div",data(i)%div)
          if(associated(data(i)%cov)) call write_hdf(file,trim(set)//"/cov",data(i)%cov)
          if(associated(data(i)%ang)) call write_hdf(file,trim(set)//"/ang",data(i)%ang)
          if(associated(data(i)%dir)) call write_hdf(file,trim(set)//"/dir",data(i)%dir)
       end do
    end if
    call close_hdf_file(file)
    call mv(tmpname, fname)
  end subroutine

  subroutine filter_comp(target, info, lim)
    implicit none
    type(quiet_target) :: target
    type(common_info)  :: info
    real(dp)           :: lim, tquv(4)
    integer(i4b)       :: di
    do di = 1, size(quiet_diodes)
       tquv = abs(quiet_diodes(di)%stokes([1,2,2,3]))
       tquv = tquv/maxval(tquv)
       if(all(tquv(info%comps) < comp_threshold)) then
          target%alist%status(quiet_diodes(di)%sub, quiet_diodes(di)%horn,:) = REJECTED_ALIST
       end if
    end do
  end subroutine

  ! Turns full(lpix) into sparse(lpix,gpix)
  subroutine list_pixels(map2mask, pixlist)
    implicit none
    integer(i4b), intent(in)   :: map2mask(0:)
    integer(i4b), allocatable  :: pixlist(:,:)
    integer(i4b)               :: i, j, n
    ! Count the number of exposed high-res pixels in map2mask
    n = count(map2mask > 0)
    if(allocated(pixlist)) deallocate(pixlist)
    allocate(pixlist(2,n))
    j = 0
    do i = 0, ubound(map2mask,1)
       if(map2mask(i) <= 0) cycle
       j = j+1
       pixlist(:,j) = [ map2mask(i), i ]
    end do
  end subroutine

end program tod2comap
