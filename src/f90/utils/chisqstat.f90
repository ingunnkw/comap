program chisqstat
  use quiet_fileutils
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_filter_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_acceptlist_mod
  use quiet_defs
  use quiet_lx_mod
  implicit none
  character(len=512)             :: parfile, acceptfile, odir, statfile, lockfile
  integer(i4b)                   :: myid, nproc, ierr, nbin, cnum, ext(5), debug
  integer(i4b)                   :: i, j, k, m, n, b, d, ndi, nces, w
  integer(i4b)                   :: nmodel, nmoment, nweight
  logical(lgt)                   :: exist, use_templates
  real(dp)                       :: srate, sfreq
  type(quiet_ces_info)           :: ces
  type(acceptlist)               :: alist
  type(task_list)                :: tasks
  type(hdf_file)                 :: file
  complex(dpc),      allocatable :: ft(:,:)
  real(dp),          allocatable :: chisq(:,:,:), model(:,:,:), sigmas(:,:,:,:), pars(:,:)
  real(dp),          allocatable :: tod(:,:), mchisqs(:,:,:,:,:), chisqs(:,:,:,:,:)
  real(dp),          allocatable :: counts(:,:,:), mcounts(:,:,:), fpar(:,:)
  real(dp),          allocatable :: means(:,:,:,:), devs(:,:,:,:), weights(:,:,:)
  real(sp),          allocatable :: temps(:,:) ! TMR adding support for weather templates
  integer(i4b),      allocatable :: cnums(:)

  call getarg(1, parfile)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call print_host_mapping

  call get_parameter(0, parfile, 'OUTPUT_DIR',            par_string=odir)
  call get_parameter(0, parfile, 'ACCEPTLIST',            par_string=acceptfile)
  call get_parameter(0, parfile, 'DEBUG',                 par_int=debug)

  ! TMR adding templates
  call get_parameter(0, parfile, 'USE_TEMPLATES',       par_lgt=use_templates)

  statfile = trim(odir) // '/stats.hdf'
  lockfile = trim(odir) // "/lock.dat"
  call dset(id=myid,level=debug)

  call dmem("init")
  call mkdirs(trim(lockfile), .true.)

  call initialize_ces_mod(parfile);           call dmem("ces mod")
  call initialize_module_mod(parfile);        call dmem("module mod")
  call initialize_filter_mod(parfile);        call dmem("filter mod")
  call initialize_accept_list(acceptfile, alist)

  ndi     = get_num_diodes() * get_num_modules()
  nbin    = 100000
  nmodel  = 2 ! Models are 1/f-profile and white noise only - might consider adding a layer for with-template model?
  nmoment = 2 ! Moments of the distribution to compute
  nweight = 2 ! Weighting schemes: Unweighted or filter weighted
  call dmem("init done")

  call get_accepted_ceses(alist, cnums)
  nces = size(cnums)
  ! We will be doing several operations in parallel on chisquares
  ! and their variances for various models. These become much
  ! more compact when handled using one array than with multiple
  ! separate arrays. The cost is the somewhat confusiong array
  ! notation. I try to make this as understandable as possible
  ! by always using explicit slicing.
  allocate(mchisqs(nbin,ndi,nmodel,nmoment,nweight))
  allocate(mcounts(nbin,ndi,nweight))
  allocate(pars(ndi,3))
  mcounts = 0; mchisqs = 0
  call dmem("alloc done")
  call init_task_list(tasks, lockfile, size(cnums), MPI_COMM_WORLD)
  do while(get_next_task(tasks, j))
     cnum = cnums(j)
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " scanning ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(get_num_ces())) // ")"
     inquire(file=trim(ces%l3file), exist=exist)
     if (.not. exist) cycle
     call open_hdf_file(ces%l3file, file, 'r')
     call read_hdf(file, 'sigma0',      pars(:,1))
     call read_hdf(file, 'fknee',       pars(:,2))
     call read_hdf(file, 'alpha',       pars(:,3))
     call read_hdf(file, 'samprate',    srate)
     call read_hdf(file, 'scanfreq',    sfreq)
     call get_size_hdf(file, 'tod', ext)
     n = ext(1); m = n/2+1
     allocate(tod(n,ndi), ft(m,ndi), chisq(m,ndi,2), model(m,ndi,2), weights(m,ndi,2))
     call read_hdf(file, 'tod', tod)
     call read_alloc_hdf(file, 'filter_par', fpar)
     call close_hdf_file(file)
     call fft_multi(tod, ft, 1)
     ! Set up the models
     do i = 1, ndi
        call get_N_filter(srate, pars(i,1), pars(i,2), pars(i,3), 1d0, model(:,i,1), .true.)
        model(:,i,2) = pars(i,1)**2
     end do
     ! Set up the weights. :,:,1 = 1, :,:,2 = filtered
     weights = 1
     do i = 1, ndi
        ! We're fetching the inverted filter here (filtered freqs-> 1e100)
        call apodize_filter_fft(.false., n, srate, fpar(i,FILTER_LOW_NU), &
         & fpar(i,FILTER_LOW_ALPHA), .false., weights(:,i,2))
        call apodize_filter_fft(.false., n, srate, fpar(i,FILTER_HIGH_NU_SCAN)*sfreq, &
         & fpar(i,FILTER_HIGH_ALPHA), .true., weights(:,i,2))
     end do
     
     weights = weights**(-2)
     ! This produces the chisquares for each model. The spread duplicates
     ! the array into the model direction, so that it is ready to be divided
     ! by the model array. The factor of two at the end is there to ensure
     ! that each entry of the array is a properly normalized chisq distribution
     ! two degrees of freedom. It has two degrees of freedom because of the
     ! two components of the complex numbers. By dividing by the exptectd model,
     ! we will end up with a number with a mean of 1 instead of 2, which needs
     ! to be compensated for.

     if (use_templates) then
        call get_templates(temps, m, srate)
        chisq = spread(real(ft * conjg(ft),dp)/temps,3,nmodel)/model*2
     else
        chisq = spread(real(ft * conjg(ft),dp),3,nmodel)/model*2
     end if

     do d = 1, ndi
        if(.not. is_accepted(alist, ces%cid, quiet_diodes(d)%horn, quiet_diodes(d)%sub)) cycle
        ! Loop over all frequency bins
        do i = 2, m
           b = int(i-1,i8b)*nbin/m+1
           do w = 1, nweight
              mcounts(b,d,w) = mcounts(b,d,w) + 2*weights(i,d,w)
              do j = 1, nmoment
                 do k = 1, nmodel
                    mchisqs(b,d,k,j,w) = mchisqs(b,d,k,j,w) + weights(i,d,w)*chisq(i,d,k)**j
                 end do
              end do
           end do
        end do
     end do
     deallocate(tod, ft, chisq, model, fpar, weights)
     if(allocated(temps)) deallocate(temps)
  end do
  allocate(chisqs(nbin,ndi,nmodel,nmoment,nweight))
  allocate(counts(nbin,ndi,nweight))
  call mpi_allreduce(mchisqs, chisqs, size(chisqs), mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(mcounts, counts, size(counts), mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  deallocate(mchisqs, mcounts)
  call free_task_list(tasks)

  ! Ok, all information has been gathered. chisqs holds the various chisqs,
  ! while counts holds the corresponding effective number of degrees of
  ! freedom. For the weighted version, this number can be fractional, as
  ! some modes have been given lower weight. In this case the number can't
  ! really be interpreted as a number of degrees of freedom.

  ! For each model, we wish to find: The chisquares and counts themselves
  ! (already available); the mean and stddev of the chisquares; and the
  ! number of standard deviations away from the model we are
  allocate(sigmas(nbin,ndi,nmodel,nweight))
  allocate(means (nbin,ndi,nmodel,nweight))
  allocate(devs  (nbin,ndi,nmodel,nweight))
  do k = 1, nmodel
     sigmas(:,:,k,:) = (chisqs(:,:,k,1,:)-counts)/(2*counts)**0.5d0
     means (:,:,k,:) = chisqs(:,:,k,1,:)/counts
     devs  (:,:,k,:) = sqrt(chisqs(:,:,k,2,:)-means(:,:,k,:)**2)
  end do

  ! Ok, all the info is available. Now we just have to output it
  if(myid == 0) then
     call dmem("output")
     call open_hdf_file(statfile, file, "w")
     call write_hdf(file, "chisqs", chisqs(:,:,:,1,:))
     call write_hdf(file, "counts", counts)
     call write_hdf(file, "sigmas", sigmas)
     call write_hdf(file, "means",  means)
     call write_hdf(file, "devs",   devs)
     call close_hdf_file(file)
  end if
  call mpi_finalize(ierr)
end program
