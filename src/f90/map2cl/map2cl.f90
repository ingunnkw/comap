! *******************************************************************************
!
!   Program for computing ML spectra from QUIET maps by Newton-Raphson optimization
!
!    Written by Sigurd K. Næss, Ingunn K. Wehus and Hans Kristian Eriksen, 2010
!
! *******************************************************************************
program map2cl
  use quiet_fileutils
  use quiet_utils
  use map2cl_nr_mod
  use quiet_task_mod
  use map2cl_parfit_mod
  implicit none
 
  integer(i4b)       :: i, j, k, numsim, seed
  integer(i4b)       :: ierr, task, status
  integer(i4b)       :: num_groups, num_procs_per_group, comm, tot_num_real
  logical(lgt)       :: analyze_data, compute_ML_spectrum, fit_qn_to_EE, fit_r_to_BB, compute_null_chisq
  logical(lgt)       :: plot_slices, fit_tabulated_amplitude
  character(len=256) :: parfile, outfile, simfile, meanfile, statfile, ksfile
  character(len=256) :: conditional_spectrum
  real(dp)           :: t1, t2, q, n, sigma_q, sigma_n, r, sigma_r
  type(map_struct), allocatable, dimension(:,:)   :: maps, my_maps
  real(dp),         allocatable, dimension(:,:,:) :: cls_fix
  real(dp),         allocatable, dimension(:)     :: lnL, chisq, stat_buff
  type(powspec),    allocatable, dimension(:,:)   :: x, x_buff
  type(planck_rng)   :: rng_handle
  type(task_list)    :: tasks
  type(common_info)  :: info
  type(scinfo)       :: info_sc
  type(scalamat)     :: mat1, mat2, mat3
  real(dp), allocatable, dimension(:,:) :: glob

  call getarg(1, parfile) 
  call get_parameter(0, parfile, 'OUTFILE',                 par_string=outfile)
  call get_parameter(0, parfile, 'SIMFILE',                 par_string=simfile)
  call get_parameter(0, parfile, 'SIM_MEAN_FILE',           par_string=meanfile)
  call get_parameter(0, parfile, 'STAT_FILE',               par_string=statfile)
  call get_parameter(0, parfile, 'KS_FILE',                 par_string=ksfile)
  call get_parameter(0, parfile, 'NUMSIM',                  par_int=numsim)
  call get_parameter(0, parfile, 'ANALYZE_DATA',            par_lgt=analyze_data)
  call get_parameter(0, parfile, 'SEED',                    par_int=seed)
  call get_parameter(0, parfile, 'NUM_PROCS_PER_GROUP',     par_int=num_procs_per_group)
  call get_parameter(0, parfile, 'NUM_GROUPS',              par_int=num_groups)
  call get_parameter(0, parfile, 'COMPUTE_ML_SPECTRUM',     par_lgt=compute_ML_spectrum)
  call get_parameter(0, parfile, 'FIT_QN_TO_EE',            par_lgt=fit_qn_to_EE)
  call get_parameter(0, parfile, 'FIT_R_TO_BB',             par_lgt=fit_r_to_BB)
  call get_parameter(0, parfile, 'FIT_TABULATED_AMPLITUDE', par_lgt=fit_tabulated_amplitude)
  call get_parameter(0, parfile, 'PLOT_SLICES',             par_lgt=plot_slices)
  call get_parameter(0, parfile, 'COMPUTE_NULL_CHISQ',      par_lgt=compute_null_chisq)
  call get_parameter(0, parfile, 'CONDITIONAL_SPECTRUM',    par_string=conditional_spectrum)
  call get_parameter(0, parfile, 'NOISE_LEVEL_TRANSFUNC',   par_dp=noise_level_transfunc)

  ! Initialize MPI groups and random seeds
  call setup_mpi(num_procs_per_group, num_groups, info, info_sc)
  call initialize_random_seeds(info%comm_world, seed, rng_handle)

  ! Initialize data module
  call initialize_data_mod(parfile, rng_handle, info_sc)  
  tot_num_real = numsim + numreal

  ! Initialize Newton-Raphson module
  call initialize_nr_mod(parfile, info_sc, info%myid_inter==0, rng_handle)

  ! Prepare data and simulations
  allocate(maps(num_data_set,tot_num_real), my_maps(num_data_set,tot_num_real))
  do i = 1, num_data_set
     do j = 1, tot_num_real
        allocate(maps(i,j)%data(data_sets(i)%ntot,1))
        allocate(my_maps(i,j)%data(data_sets(i)%ntot,1))
        maps(i,j)%data    = 0.d0
        my_maps(i,j)%data = 0.d0
     end do
  end do
  if (info%myid == 0) then
     do i = 1, num_data_set
        do j = 1, numreal
           ! First numreal entries are the input maps
           my_maps(i,j)%data = data_sets(i)%maps(:,j:j) 
           if (noise_level_transfunc > 0.d0) call add_noise_to_data(info_sc, i, data_sets(i)%maps(:,j))
        end do
     end do
  end if
  do i = 1+info%myid_inter, numsim, info%numprocs_inter
     call get_simulation(info_sc, my_maps(:,numreal+i))
     if (info_sc%myid /= 0) then
        do j = 1, num_data_set
           my_maps(j,numreal+i)%data = 0.d0
        end do
     end if
  end do
  do i = 1, num_data_set
     do j = 1, tot_num_real
        call mpi_allreduce(my_maps(i,j)%data, maps(i,j)%data, size(maps(i,j)%data), &
             & MPI_DOUBLE_PRECISION, MPI_SUM, info%comm_world, ierr)
        deallocate(my_maps(i,j)%data)
     end do
  end do
  deallocate(my_maps)
  
  ! Clean up unnecessary data structures before starting computations
  call cleanup_data_mod

  ! Allocate data structures
  allocate(chisq(tot_num_real))
  allocate(lnL(tot_num_real))
  allocate(stat_buff(tot_num_real))
  allocate(cls_fix(0:lmax,6,num_data_set))
  allocate(x(tot_num_real,2), x_buff(tot_num_real,2))
  chisq = 0.d0; lnL = 0.d0
  do i = 1, tot_num_real
     do j = 1, 2
        call copy_powspec(x_template, x(i,j))
        call copy_powspec(x_template, x_buff(i,j))
     end do
  end do

  ! Analyze data and simulations
  call init_task_list(tasks, 'map2cl_tasks.dat', tot_num_real, info%comm_world)
  t1 = 0.d0; t2 = 0.d0
  if (info_sc%myid == 0) then

     do while(get_next_task(tasks, i))
        if (.not. analyze_data .and. i <= numreal) cycle
        write(*,fmt='(a,i4,a,i5,a,i5,a,f6.2)') 'Myid = ', info%myid, ' -- data set ', i, &
             & ' of ', tot_num_real, ', CPU = ', t2-t1
        
        ! Let other nodes in this group know which data set to analyze
        call mpi_bcast(i, 1, MPI_INTEGER, 0, info_sc%comm, ierr)

        ! Compute ML spectrum
        if (compute_ML_spectrum) then
           call wall_time(t1)
           call compute_max_lnL_spectrum(info_sc, i, maps(:,i), x(i,1), x(i,2), &
                & lnL(i), chisq(i), status)
           if (status /= 0) then
              write(*,*) 'WARNING: NR search did not converge for realization = ', i
           end if

           call wall_time(t2)
           write(*,fmt='(a,i4,a,i5,a,f10.2,a,f10.2)') 'Myid = ', info%myid, ' -- data set no. ', i, &
                & ' chisq = ', chisq(i), ', lnL = ', lnL(i)
        end if

        ! Set up conditional power spectrum
        if (fit_qn_to_EE .or. fit_r_to_BB .or. plot_slices .or. compute_null_chisq) then
           if (info_sc%myid == 0) then
              if (trim(conditional_spectrum) == 'fiducial') then
                 do k = 1, num_data_set
                    cls_fix(:,:,k) = cl_fid 
                 end do
              else if (trim(conditional_spectrum) == 'precomputed') then
                 call input_powspec(cl_fid, cls_fix, realization=i)
              else if (trim(conditional_spectrum) == 'internal') then
                 call convert_x2cls(x(i,1), cls_fix)
              else
                 write(*,*) 'Unknown conditional spectrum time = ', trim(conditional_spectrum)
                 write(*,*) 'Valid options are {fiducial, precomputed, internal}'
                 stop
              end if
           end if
           call mpi_bcast(cls_fix, size(cls_fix), MPI_DOUBLE_PRECISION, 0, info_sc%comm, ierr)
        end if

        ! Fit q-n model to EE spectrum
        if (fit_qn_to_EE) then
           call fit_qn_model(info_sc, parfile, cls_fix, maps(:,i), i, q, sigma_q, n, sigma_n)
           write(*,fmt='(a,i4,a,i5,a,f5.2,a,f5.2,a,f5.2,a,f5.2)') 'Myid = ', info%myid, ' -- data set no. ', i, &
                & ' q = ', q, '+/-', sigma_q, ', n = ', n, '+/-', sigma_n
        end if

        ! Fit r to BB spectrum
        if (fit_r_to_BB) then
           call fit_r2(info_sc, info%myid_inter==0, parfile, cls_fix, maps(:,i), i, r, sigma_r)
           write(*,fmt='(a,i4,a,i5,a,f5.2,a,f5.2)') 'Myid = ', info%myid, ' -- data set no. ', i, &
                & ' r = ', r, '+/-', sigma_r
        end if

        ! Fit amplitude versus tabulated spectra
        if (fit_tabulated_amplitude) then
           call fit_tab_amplitude(info_sc, parfile, maps(:,i), i)
        end if

        ! Plot slices through the likelihood
        if (plot_slices) then
           call plot_likelihood_slices(info_sc, parfile, cls_fix, maps(:,i), i)
        end if

        ! Compute chi-squares relative to null spectrum
        if (compute_null_chisq) then
           call compute_chisq_vs_null(info_sc, parfile, cls_fix, maps(:,i), i)
        end if

     end do

     ! Free up processors
     call mpi_bcast(-1, 1, MPI_INTEGER, 0, info_sc%comm, ierr)

  else

     call mpi_bcast(i, 1, MPI_INTEGER, 0, info_sc%comm, ierr)
     do while (i >= 0)
        if (compute_ML_spectrum) then
           call compute_max_lnL_spectrum(info_sc, i, maps(:,i), x(i,1), x(i,2), &
                & lnL(i), chisq(i), status)
        end if
        if (fit_qn_to_EE .or. fit_r_to_BB .or. plot_slices .or. compute_null_chisq) then
           call mpi_bcast(cls_fix, size(cls_fix), MPI_DOUBLE_PRECISION, 0, info_sc%comm, ierr)
        end if
        if (fit_qn_to_EE) then
           call fit_qn_model(info_sc, parfile, cls_fix, maps(:,i), i, q, sigma_q, n, sigma_n)
        end if
        if (fit_r_to_BB) then
           call fit_r2(info_sc, info%myid_inter==0, parfile, cls_fix, maps(:,i), i, r, sigma_r)
        end if
        if (plot_slices) then
           call plot_likelihood_slices(info_sc, parfile, cls_fix, maps(:,i), i)
        end if
        if (compute_null_chisq) then
           call compute_chisq_vs_null(info_sc, parfile, cls_fix, maps(:,i), i)
        end if
        if (fit_tabulated_amplitude) then
           call fit_tab_amplitude(info_sc, parfile, maps(:,i), i)
        end if
        call mpi_bcast(i, 1, MPI_INTEGER, 0, info_sc%comm, ierr)
     end do
     
  end if

  ! Collect and post-process results 
  if (info%myid_group == info%root) then
     do i = 1, tot_num_real
        do j = 1, 2
           call mpi_reduce(x(i,j)%coeff, x_buff(i,j)%coeff, size(x(i,j)%coeff), MPI_DOUBLE_PRECISION, MPI_SUM, &
                & info%root, info%comm_inter, ierr)
           x(i,j)%coeff = x_buff(i,j)%coeff
        end do
     end do
     call mpi_reduce(chisq, stat_buff, size(chisq), MPI_DOUBLE_PRECISION, MPI_SUM, &
          & info%root, info%comm_inter, ierr)
     chisq = stat_buff
     call mpi_reduce(lnL, stat_buff, size(lnL), MPI_DOUBLE_PRECISION, MPI_SUM, &
          & info%root, info%comm_inter, ierr)
     lnL = stat_buff

     if (info%myid == info%root) then
        ! Output results
        do i = 1, tot_num_real
           call output_powspec(simfile, i/=1, x(i,1), x(i,2))
        end do
        if (numsim > 0) call output_mean_powspec(meanfile, x(numreal+1:numreal+numsim,1))
        call output_statistics(statfile, lnL, chisq, sum(np))
!        if (analyze_data) call output_p_values(ksfile, x(num_lowl_junk_bins+1:npar-num_highl_junk_bins,:,1))
     end if

  end if
  

  ! Clean up
  deallocate(chisq, lnL, stat_buff, x, x_buff, maps, cls_fix)
  call free_task_list(tasks)
  call mpi_finalize(ierr)

contains

  subroutine setup_mpi(num_procs_per_group, num_groups, info, info_sc)
    implicit none
    
    integer(i4b),      intent(in)  :: num_procs_per_group, num_groups
    type(common_info), intent(out) :: info
    type(scinfo),      intent(out) :: info_sc
    
    integer(i4b) :: icol, irow, ierr

    info%comm_world = MPI_COMM_WORLD    
    call mpi_init(ierr)
    call mpi_comm_rank(info%comm_world, info%myid, ierr)
    call mpi_comm_size(info%comm_world, info%numprocs, ierr)
    info%root = 0
    
    ! Check that the number of processors match the number of chains and maps
    if (num_groups*num_procs_per_group /= info%numprocs) then
       if (info%myid == info%root) then
          write(*,*) ''
          write(*,*) 'ERROR: Number of processors is not a multiple of numprocs_covar'
          write(*,*) '       numgroups           = ', num_groups
          write(*,*) '       num_proc_per_groups = ', num_procs_per_group
          write(*,*) '       numprocs            = ', info%numprocs, ' /= ', num_groups*num_procs_per_group
          write(*,*) ''
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
    ! Split processors into work groups
    icol = info%myid/num_procs_per_group
    irow = mod(info%myid,num_procs_per_group)
    call mpi_comm_split(info%comm_world, icol, irow, info%comm_group, ierr) ! comm for each work group
    call mpi_comm_split(info%comm_world, irow, icol, info%comm_inter, ierr) ! comm for different work groups
    call mpi_comm_rank(info%comm_group, info%myid_group, ierr)
    call mpi_comm_rank(info%comm_inter, info%myid_inter, ierr)
    call mpi_comm_size(info%comm_group, info%numprocs_group, ierr)
    call mpi_comm_size(info%comm_inter, info%numprocs_inter, ierr)

    ! Set up Scalapack communicator
    call sc_init(info_sc, startcomm=info%comm_group, active=.true.)

  end subroutine setup_mpi

end program map2cl
