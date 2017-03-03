program corrstat
  use quiet_fileutils
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_acceptlist_mod
  use quiet_gain_mod
  use quiet_lx_mod
  implicit none
  character(len=512)             :: parfile, acceptfile, odir, statfile, lockfile
  integer(i4b)                   :: myid, nproc, ierr, nbin, cnum, ext(5), debug
  integer(i4b)                   :: i, j, k, m, n, b, d, d1, d2, ndi, nces, nmod, mdi
  logical(lgt)                   :: exist
  type(quiet_ces_info)           :: ces
  type(acceptlist)               :: alist
  type(task_list)                :: tasks
  type(hdf_file)                 :: file
  integer(i4b),      allocatable :: cnums(:)
  real(dp),          allocatable :: cov(:,:,:,:,:), corr(:,:,:), sigma0(:), gtime(:)
  real(dp),          allocatable :: gain(:), noise(:), ocov(:,:,:,:), icov(:,:,:)
  real(dp),          allocatable :: ocorr(:,:,:,:), modcov(:,:,:,:,:,:)

  call getarg(1, parfile)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call print_host_mapping

  call get_parameter(0, parfile, 'OUTPUT_DIR',            par_string=odir)
  call get_parameter(0, parfile, 'ACCEPTLIST',            par_string=acceptfile)
  call get_parameter(0, parfile, 'DEBUG',                 par_int=debug)
  statfile = trim(odir) // '/stats.hdf'
  lockfile = trim(odir) // "/lock.dat"
  call dset(id=myid,level=debug)

  call dmem("init")
  call mkdirs(trim(lockfile), .true.)

  call initialize_ces_mod(parfile);           call dmem("ces mod")
  call initialize_gain_mod(parfile);          call dmem("gain mod")
  call initialize_accept_list(acceptfile, alist)

  ndi  = size(quiet_diodes)
  nmod = size(quiet_horns)
  mdi  = get_num_diodes()
  nbin = 25
  allocate(cov(ndi,ndi,nbin,5,2), noise(ndi), gain(ndi), sigma0(ndi))
  cov = 0

  call get_accepted_ceses(alist, cnums)
  nces = size(cnums)
  call init_task_list(tasks, lockfile, size(cnums), MPI_COMM_WORLD)
  do while(get_next_task(tasks, j))
     cnum = cnums(j)
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " scanning ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(get_num_ces())) // ")"
     inquire(file=trim(ces%l3file), exist=exist)
     if (.not. exist) cycle
     call open_hdf_file(ces%l3file, file, "r")
     call read_alloc_hdf(file, "time_gain", gtime)
     call read_alloc_hdf(file, "corr", corr)
     call read_hdf(file, "sigma0", sigma0)
     do d = 1, ndi
        gain(d) = get_gain(gtime(1), d) / ant2thermo(quiet_diodes(d)%freq) / 1e3
     end do
     noise = sigma0/gain
     do d1 = 1, ndi
        if(.not. is_accepted(alist, ces%cid, quiet_diodes(d1)%horn,&
          quiet_diodes(d1)%sub)) cycle
        do d2 = 1, ndi
           if(.not. is_accepted(alist, ces%cid, quiet_diodes(d2)%horn,&
             quiet_diodes(d2)%sub)) cycle
           cov(d1,d2,:,1,1) = cov(d1,d2,:,1,1) + corr(:,d1,d2)*noise(d1)*noise(d2)
           cov(d1,d2,:,2,1) = cov(d1,d2,:,2,1) +(corr(:,d1,d2)*noise(d1)*noise(d2))**2
           cov(d1,d2,:,3,1) = cov(d1,d2,:,3,1) + corr(:,d1,d2)
           cov(d1,d2,:,4,1) = cov(d1,d2,:,4,1) + corr(:,d1,d2)**2
           cov(d1,d2,:,5,1) = cov(d1,d2,:,5,1) + 1
        end do
     end do
     call close_hdf_file(file)
     deallocate(gtime)
  end do
  call dmem("A")
  call mpi_allreduce(cov(:,:,:,:,1), cov(:,:,:,:,2), size(cov(:,:,:,:,1)), mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  call dmem("B")
  if(myid == 0) then
     allocate(ocov(ndi,ndi,nbin,4),icov(ndi,ndi,nbin),ocorr(ndi,ndi,nbin,2))
     allocate(modcov(mdi,mdi,nmod,nbin,2,2))
     call dmem("C")
     ! Reduce to mean and dev of cov entries, with 0 as placeholder
     ocov = cov(:,:,:,1:4,2)/cov(:,:,:,[5,5,5,5],2)
     call dmem("D")
     where(cov(:,:,:,5,2)==0)
        ocov(:,:,:,1) = 0
        ocov(:,:,:,2) = 0
        ocov(:,:,:,3) = 0
        ocov(:,:,:,4) = 0
     elsewhere
        ocov(:,:,:,2) = sqrt(ocov(:,:,:,2)-ocov(:,:,:,1)**2)
        ocov(:,:,:,4) = sqrt(ocov(:,:,:,4)-ocov(:,:,:,3)**2)
     end where
     call dmem("E")
     ! Compute corresponding correlation matrix
     ocorr = 0
     call dmem("F")
     do d1 = 1, ndi
        do d2 = 1, ndi
           if(all(ocov(d1,d1,:,1:2) > 0 .and. ocov(d2,d2,:,1:2) > 0)) then
              ocorr(d1,d2,:,:) = ocov(d1,d2,:,1:2)/sqrt(ocov(d1,d1,:,[1,1])*ocov(d2,d2,:,[1,1]))
           end if
        end do
     end do
     call dmem("G")
   
     ! Inverse correlation matrix based on full matrix
     icov = ocov(:,:,:,1)
     do i = 1, nbin
        call invert_matrix_with_mask(icov(:,:,i))
     end do
     call dmem("H")
   
     ! Per module
     do i = 1, nmod
        modcov(:,:,i,:,:,1) = ocov((i-1)*mdi+1:i*mdi,(i-1)*mdi+1:i*mdi,:,1:2)
        modcov(:,:,i,:,:,2) = modcov(:,:,i,:,:,1)
        do j = 1, nbin
           call invert_matrix_with_mask(modcov(:,:,i,j,1,2))
        end do
     end do
     call dmem("I")

     call open_hdf_file(statfile, file, "w")
     call write_hdf(file,  "cov", ocov)
     call write_hdf(file, "icov", icov)
     call write_hdf(file, "corr", ocorr)
     call write_hdf(file, "mcov", modcov(:,:,:,:,:,1))
     call write_hdf(file, "imcov", modcov(:,:,:,:,1,2))
  end if
  call mpi_finalize(ierr)
end program
