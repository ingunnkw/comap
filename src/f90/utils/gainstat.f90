program gainstat
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
  integer(i4b)                   :: i, j, k, m, n, b, d, ndi, nces
  logical(lgt)                   :: exist
  type(quiet_ces_info)           :: ces
  type(acceptlist)               :: alist
  type(task_list)                :: tasks
  type(hdf_file)                 :: file
  real(dp),       allocatable    :: gain(:,:), mystat(:,:), stat(:,:), a2t(:), tmp(:,:)
  real(dp),       allocatable    :: gtime(:), sigma0(:), w(:), corr(:,:,:), Smat(:,:)
  real(dp),       allocatable    :: mySmat(:,:), gtmp(:), corr2(:,:)
  integer(i4b),   allocatable    :: cnums(:), mynsmat(:,:), nsmat(:,:)

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

  call get_accepted_ceses(alist, cnums)
  nces = size(cnums)
  call init_task_list(tasks, lockfile, size(cnums), MPI_COMM_WORLD)
  allocate(mystat(ndi,5), stat(ndi,5), sigma0(ndi), Smat(ndi,ndi), mySmat(ndi,ndi))
  allocate(gtmp(ndi),corr2(ndi,ndi),nsmat(ndi,ndi),mynsmat(ndi,ndi))
  mystat = 0
  mySmat = 0
  mynsmat = 0
  do while(get_next_task(tasks, j))
     cnum = cnums(j)
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " scanning ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(get_num_ces())) // ")"
     inquire(file=trim(ces%l3file), exist=exist)
     if (.not. exist) cycle
     call open_hdf_file(ces%l3file, file, "r")
     !call read_alloc_hdf(file, "gain", gain)
     call read_alloc_hdf(file, "time_gain", gtime)
     call read_alloc_hdf(file, "corr", corr)
     call read_hdf(file, "sigma0", sigma0)
     call close_hdf_file(file)
     where(corr/=corr) corr = 0
     m = size(corr,1)
     corr2 = 0
     do d = 1, ndi, 4
        corr2(d:d+3,d:d+3) = corr(m,d:d+3,d:d+3)
        call invert_matrix_with_mask(corr2(d:d+3,d:d+3))
     end do
     do d = 1, ndi
        gtmp(d) = get_gain(gtime(1), d) / ant2thermo(quiet_diodes(d)%freq) / 1e3
     end do
     do i = 1, ndi; do j = 1, ndi
        if(sigma0(i) > 0 .and. sigma0(j) > 0 .and. gtmp(i) > 0 .and. gtmp(j) > 0) then
           corr2(i,j) = corr2(i,j)/sigma0(i)/sigma0(j)*gtmp(i)*gtmp(j)
           mynsmat(i,j) = mynsmat(i,j) + 1
        else
           corr2(i,j) = 0
        end if
     end do; end do
     mySmat = mySmat + corr2

     allocate(gain(size(gtime),ndi),w(size(gtime)))
     do d = 1, ndi
        if(.not. is_accepted(alist, ces%cid, quiet_diodes(d)%horn, quiet_diodes(d)%sub)) cycle
        call get_gains(gtime, d, gain(:,d))
        w = gain(:,d)/sigma0(d)**2
        mystat(d,1) = mystat(d,1) + size(w)
        mystat(d,2) = mystat(d,2) + sum(gain(:,d)*w)
        mystat(d,3) = mystat(d,3) + sum(gain(:,d)**2*w)
        mystat(d,4) = mystat(d,4) + sum(w)
        mystat(d,5) = mystat(d,5) + sum(sigma0(d)/(gain(:,d)*1d-3)*w)
     end do
     deallocate(gain, gtime, w, corr)
  end do
  call mpi_allreduce(mystat, stat, size(stat), mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(mySmat, Smat, size(Smat), mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(mynsmat, nsmat, size(nsmat), mpi_integer, MPI_SUM, mpi_comm_world, ierr)
  deallocate(mystat,mySmat,mynsmat)
  Smat = Smat/nsmat
  call free_task_list(tasks)
  allocate(a2t(ndi),tmp(ndi,2))
  do d = 1, ndi
     a2t(d) = ant2thermo(get_diode_freq(d))
  end do
  if(myid == 0) then
     call dmem("output")
     call open_hdf_file(statfile, file, "w")
     call write_hdf(file, "counts", int(stat(:,1)))
     stat(:,2:3) = stat(:,2:3)/stat(:,[4,4])
     tmp(:,1) = stat(:,2); tmp(:,2) = (stat(:,3)-stat(:,2)**2)**0.5
     call write_hdf(file, "gain",   tmp)
     call write_hdf(file, "thermo", tmp/spread(a2t,2,2))
     call write_hdf(file, "SN", stat(:,4)/stat(:,1))
     call write_hdf(file, "sens", stat(:,5)/stat(:,4)*a2t/sqrt(25.0)*1e6)
     call write_hdf(file, "smat", Smat)
     call close_hdf_file(file)
  end if
  deallocate(a2t, stat, tmp)
  call mpi_finalize(ierr)
end program
