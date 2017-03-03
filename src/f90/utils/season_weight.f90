! Output a file with information about the CES-diodes and their weights
! in the data processing. This can then be used to quickly compute
! various season averages.
program season_weight
  use quiet_fileutils
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_acceptlist_mod
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
  real(dp),          allocatable :: time(:), dur(:), weight(:,:), myweight(:,:)
  real(dp),          allocatable :: sigma0(:)
  integer(i4b),      allocatable :: cid(:), cnums(:)

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
  call init_detector_mod(parfile);           call dmem("detector mod")
  call initialize_accept_list(acceptfile, alist)

  ndi  = size(quiet_diodes)

  call get_accepted_ceses(alist, cnums)
  nces = size(cnums)
  call init_task_list(tasks, lockfile, size(cnums), MPI_COMM_WORLD)
  allocate(time(nces), dur(nces), weight(nces,ndi), myweight(nces,ndi), cid(nces))
  allocate(sigma0(ndi))
  time = (ces_db%ceses(cnums)%mjd(2)+ces_db%ceses(cnums)%mjd(1))/2
  dur  = ces_db%ceses(cnums)%mjd(2)-ces_db%ceses(cnums)%mjd(1)
  cid  = ces_db%ceses(cnums)%cid
  weight   = 0
  myweight = 0
  do while(get_next_task(tasks, j))
     cnum = cnums(j)
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " scanning ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(get_num_ces())) // ")"
     inquire(file=trim(ces%l3file), exist=exist)
     if (.not. exist) cycle
     call open_hdf_file(ces%l3file, file, "r")
     call read_hdf(file, "sigma0", sigma0)
     call close_hdf_file(file)
     do d = 1, ndi
        if(.not. is_accepted(alist, ces%cid, quiet_diodes(d)%horn, quiet_diodes(d)%sub)) cycle
        myweight(j,d) = dur(j)/sigma0(d)**2
     end do
  end do
  call mpi_allreduce(myweight, weight, size(weight), mpi_double_precision, MPI_SUM, &
   & mpi_comm_world, ierr)
  deallocate(myweight)
  call free_task_list(tasks)
  if(myid == 0) then
     call dmem("output")
     call open_hdf_file(statfile, file, "w")
     call write_hdf(file, "ces", cid)
     call write_hdf(file, "dur", dur)
     call write_hdf(file, "time", time)
     call write_hdf(file, "weight", weight)
     call close_hdf_file(file)
  end if
  deallocate(cid, dur, time, weight, sigma0)
  call mpi_finalize(ierr)
end program
