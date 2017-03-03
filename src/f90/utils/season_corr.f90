! Scan through level3-files in a runlist, and accumulate
! correlation info
program season_corr
  use quiet_ces_mod
  use quiet_detector_mod
  use quiet_utils
  use quiet_hdf_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_fileutils
  use quiet_fft_mod
  implicit none
  character(len=512)   :: parfile, arg, line, lockfile, outfile, odir
  type(task_list)      :: tasks
  type(quiet_ces_info) :: ces
  type(hdf_file)       :: file
  integer(i4b)         :: cnum, err, myid, nproc, debug, nbin, ndi, n, i, j, a, b
  integer(i4b)         :: d1, d2, nactloc, nactive, ext(7)
  real(dp)             :: srate, sloc, nloc, ntot
  integer(i4b), allocatable :: bins(:,:), hits(:,:,:), tothits(:,:,:)
  real(sp),     allocatable :: tod(:,:), corr(:,:,:), totcorr(:,:,:), norm(:,:)
  real(dp),     allocatable :: freqs(:)
  complex(spc), allocatable :: ffts(:,:)

  call getarg(1, parfile)
  nbin = 1000

  call get_parameter(0, parfile, 'OUTPUT_DIR', par_string=odir)
  call get_parameter(0, parfile, 'DEBUG',      par_int=debug)

  call mpi_init(err)
  call mpi_comm_rank(MPI_COMM_WORLD, myid,  err)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, err)
  call dset(id=myid,level=debug)
  call print_host_mapping

  call initialize_ces_mod(parfile);            call dmem("ces mod")
  call init_detector_mod(parfile)

  ndi      = size(quiet_diodes)
  lockfile = trim(odir) // "/lock.dat"
  outfile  = trim(odir) // "/corr.dat"
  call mkdirs(trim(lockfile), .true.)

  allocate(bins(2,nbin),corr(nbin,ndi,ndi),hits(nbin,ndi,ndi),freqs(nbin))

  corr = 0
  hits = 0
  nloc = 0
  sloc = 0
  nactloc = 0

  call dmem("init")

  ! Process all CES's
  call init_task_list(tasks, lockfile, get_num_ces(), MPI_COMM_WORLD)
  do while(get_next_task(tasks, cnum))
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " processing ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(get_num_ces())) // ")"

     ! Read tod and extract fft
     call open_hdf_file(ces%l2file, file, "r")
     call get_size_hdf(file, "tod", ext)
     allocate(tod(ext(1),ext(2)),ffts(ext(1)/2+1,ext(2)))
     call read_hdf(file, "tod",      tod)
     call read_hdf(file, "samprate", srate)
     call close_hdf_file(file)
     call dmem("read")
     call fft_multi(tod, ffts, 1)
     deallocate(tod)
     call dmem("fourier")

     ! Calculate the correlations
     n    = size(ffts,1)
     nloc = nloc + n
     sloc = sloc + srate
     nactloc = nactloc + 1
     if(n < nbin) goto 1
     call make_exp_bins(n, bins)
     allocate(norm(n,ndi))
     norm = sqrt(real(ffts*conjg(ffts)))

     do d1 = 1, ndi
        corr(:,d1,d1) = corr(:,d1,d1) + bins(2,:) - bins(1,:) + 1
        hits(:,d1,d1) = hits(:,d1,d1) + bins(2,:) - bins(1,:) + 1
        do d2 = d1+1, ndi
           do i = 1, nbin
              a = bins(1,i); b = bins(2,i)
              corr(i,d1,d2) = corr(i,d1,d2) + sum(real(ffts(a:b,d1)*conjg(ffts(a:b,d2)))/(norm(a:b,d1)*norm(a:b,d2)))
              hits(i,d1,d2) = hits(i,d1,d2) + b-a+1
           end do
        end do
     end do
     call dmem("corr")

1    deallocate(ffts, norm)
  end do
  call dmem("reducing")

  allocate(totcorr(nbin,ndi,ndi),tothits(nbin,ndi,ndi))
  call mpi_reduce(corr, totcorr, size(corr), mpi_real,    mpi_sum, 0, mpi_comm_world, err)
  call mpi_reduce(hits, tothits, size(hits), mpi_integer, mpi_sum, 0, mpi_comm_world, err)
  call mpi_reduce(nloc, ntot,    1,          mpi_double_precision, mpi_sum, 0, mpi_comm_world, err)
  call mpi_reduce(sloc, srate,   1,          mpi_double_precision, mpi_sum, 0, mpi_comm_world, err)
  call mpi_reduce(nactloc, nactive, 1,       mpi_integer, mpi_sum, 0, mpi_comm_world, err)
  if(myid == 0) then
     do d1 = 1, ndi
        do d2 = d1+1, ndi
           totcorr(:,d2,d1) = totcorr(:,d1,d2)
           tothits(:,d2,d1) = tothits(:,d1,d2)
        end do
     end do
     totcorr = totcorr / tothits
     srate   = srate   / nactive
     n       = nint(ntot / nactive)
     call make_exp_bins(n, bins)
     do i = 1, nbin
        freqs(i) = (ind2freq(bins(1,i), srate, n)+ind2freq(bins(2,i), srate, n))/2
     end do
     call open_hdf_file(outfile, file, "w")
     call write_hdf(file, "corr", totcorr)
     call write_hdf(file, "hits", tothits)
     call write_hdf(file, "freq", freqs)
     call close_hdf_file(file)
  end if
  call dmem("done")

  deallocate(corr, totcorr, hits, tothits, bins, freqs)
  call mpi_barrier(mpi_comm_world, err)
  if(myid == 0) write(*,*) "Finished"
  call mpi_finalize(err)

end program
