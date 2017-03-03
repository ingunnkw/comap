! Building a weather template because the 1/f QUIET noise model is seemingly unable to account for what looks like weather excess in the area around 0.5-1 Hz (in Q-band)
! Including templates, powerspectrum is modelled as P(nu,diode) = s0**2 (1 + (f/fknee)**alpha) * T(nu,diode)

program wtemplate
  use quiet_utils
  use quiet_fileutils
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_fft_mod
  use quiet_acceptlist_mod
  implicit none

  character(len=512)                           :: parfile, baseline_alist, full_alist, odir, lockfile, path, dumpfile
  integer(i4b)                                 :: myid, nproc, ierr, debug, ndi, nces, cnum, ext(2), nsamp, nf, numbins, diode, mod
  integer(i4b)                                 :: i, di, bin, unit, j, a, b, iter, bin_lowlimit, bin_highlimit(1)
  integer(i4b), dimension(:),      allocatable :: acc_ceses, n
  integer(i4b), dimension(:,:),    allocatable :: mask
  real(dp)                                     :: srate, fact, template_freq_limit_high
  real(dp)                                     :: template_freq_limit_low, ymin, ymax, dy, totmean, totsd
  real(dp),     dimension(:),      allocatable :: s0, alpha, fknee, f
  real(dp),     dimension(:,:),    allocatable :: tods, pows, ratio, mratios, msquares, mcounts
  real(dp),     dimension(:,:),    allocatable :: ratio_tot, sd_tot, counts_tot, runav
  real(dp),     dimension(:,:,:),  allocatable :: smooth
  complex(dpc), dimension(:,:),    allocatable :: ffts
  logical(lgt)                                 :: exists
  type(task_list)                              :: tasks
  type(quiet_ces_info)                         :: ces
  type(hdf_file)                               :: file
  type(acceptlist)                             :: alist

  call getarg(1, parfile)
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call print_host_mapping

  call get_parameter(0, parfile, 'OUTPUT_DIR',            par_string=odir)
  call get_parameter(0, parfile, 'DEBUG',                 par_int=debug)
  !  call get_parameter(0, parfile, 'ACCEPT_LIST_INPUT',     par_string=baseline_alist)
  call get_parameter(0, parfile, 'ACCEPTLIST',            par_string=full_alist)
  lockfile = trim(odir) // "/lock.dat"
  call dset(id=myid,level=debug)

  call dmem("init")
  call mkdirs(trim(lockfile), .true.)
  ! Extra folder for raw templates
  call mkdirs(trim(trim(odir)//'/unsmoothed'), .false.)

  ! Initializing necessary modules
  call initialize_ces_mod(parfile);           call dmem("ces mod")
  call initialize_module_mod(parfile);        call dmem("module mod")

  ! Adding support for input acceptlist to avoid the crazy stuff
  ! on second thought, use full acceptlist - why would we want to use all the bad scans to make templates?
  ! Oh, this is great. In acceptlists, 0 means reject, 1 means accept. In acceptlist_mod, alist%status=0 means accepted, 1 (or larger) means rejected.
  !  call initialize_accept_list(baseline_alist, alist, REJECTED_NONE)
  call initialize_accept_list(full_alist, alist)
  call dmem("after init")

  ! dont want to bother with the ones we don't accept
  call get_accepted_ceses(alist,acc_ceses)
  nces = size(acc_ceses)
  ndi  = get_num_diodes() * get_num_modules()
  numbins = 1000 ! or whatever - in effect we have numbins+1 bins
  ymin = log10(1.d-4)
  ymax = log10(12.5d0)
  dy = ymax - ymin

  template_freq_limit_high = 4.5 ! Q-band specific. Might want to un-hardcode this
  template_freq_limit_low = 5.d-3 

  
  allocate(s0(ndi),alpha(ndi),fknee(ndi))
!  allocate(mratios(numbins,ndi),msquares(numbins,ndi),mcounts(numbins,ndi))
  allocate(mratios(0:numbins,ndi),msquares(0:numbins,ndi),mcounts(0:numbins,ndi))
  mratios=0.d0
  msquares=0.d0
  mcounts=0.d0
  exists = .false.

  call dmem("after alloc")
  call init_task_list(tasks, lockfile, nces, MPI_COMM_WORLD)
  do while(get_next_task(tasks,j))     
     cnum = acc_ceses(j)
     call get_ces_info(cnum, ces)                   
     write(*,fmt="(i3,a,i4,a)") myid, " scanning ces ", ces%cid, " (" // trim(itoa(j)) // "/" //  trim(itoa(nces)) // ")"
     inquire(file=trim(ces%l3file), exist=exists)
     if (.not. exists) cycle
     call open_hdf_file(ces%l3file, file, 'r')
     call read_hdf(file, 'sigma0', s0)
     call read_hdf(file, 'alpha', alpha)
     call read_hdf(file, 'fknee', fknee)
     call read_hdf(file, 'samprate', srate)
     call get_size_hdf(file, 'tod', ext)
     nsamp = ext(1)
     nf = nsamp/2 + 1
     allocate(tods(nsamp,ndi), ffts(nf,ndi), pows(nf,ndi))
     call read_hdf(file, 'tod', tods)
     call close_hdf_file(file)

     ! Computing powerspectra
     call fft_multi(tods, ffts, 1)
     do di = 1, ndi
        call extract_powspec(ffts(:,di), pows(:,di))
     end do
     deallocate(ffts)

     allocate(ratio(nf,ndi),f(nf))
     call dmem("Building templates")

     ! Load frequencies
     do i=1,nf
        f(i) = ind2freq(i,srate,nf)
     end do

     ! Compute ratio of powerspectrum and 1/f noise model
     do di = 1,ndi
        if (fknee(di)>0.d0) then
           ratio(:,di) = pows(:,di)/(s0(di)**2*(1.d0+(f/fknee(di))**alpha(di)))
        else
           ! This check still hits. can still be bad scans in an accepted ces! 
           ratio(:,di) = 0.d0
        end if
     end do

     ! Building binned ratio and counts for each diode
     do di = 1,ndi
        mod = (di-1)/get_num_diodes()
        diode  = modulo(di-1,get_num_diodes())
        if(.not. is_accepted(alist, ces%cid, mod, diode)) cycle
        do i = 1,nf 
!!$           if (log10(f(i))<ymin) then
!!$              bin=1
!!$           else
!!$              ! This needs to run from 2 to numbins
!!$              bin = int((log10(f(i))-ymin)*(numbins-2)/dy) + 2 ! stuffing everything below min freq into lowest 2 bins
!!$           end if
           if (f(i)> 10.d0**ymin) then ! Skipping whatevers below ymin (not planning on using the lowest freqs anyway)
              bin = int((log10(f(i))-ymin)*(numbins-1)/dy) + 1
           else
              cycle
           end if

           !           bin = int(i-1,i8b)*numbins/nf +1 
           mratios(bin,di) = mratios(bin,di) + ratio(i,di)
           msquares(bin,di)= msquares(bin,di) + ratio(i,di)**2
           mcounts(bin,di) = mcounts(bin,di) + 1.d0
        end do
     end do
     deallocate(tods,pows,ratio,f)
  end do

!  allocate(ratio_tot(numbins,ndi),sd_tot(numbins,ndi),counts_tot(numbins,ndi))
  allocate(ratio_tot(0:numbins,ndi),sd_tot(0:numbins,ndi),counts_tot(0:numbins,ndi), mask(0:numbins,ndi))

! Reduce ces-binned-ratios to total ratios, sum of squares, and counts
  call mpi_allreduce(mratios, ratio_tot,  size(ratio_tot),  mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(msquares, sd_tot,    size(sd_tot),     mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  call mpi_allreduce(mcounts, counts_tot, size(counts_tot), mpi_double_precision, MPI_SUM, mpi_comm_world, ierr)
  deallocate(mratios,msquares,mcounts)
  call free_task_list(tasks)
  deallocate(s0,alpha,fknee)
!!$
!!$  do i=0,numbins
!!$     do j=1,ndi
!!$        if((myid==0).and. (sd_tot(i,j)<1.d-20)) then
!!$           if((j.ne.17) .and. (j.ne.35) .and. (j.ne.65) .and. (j.ne.66) .and. (j.ne.67) .and. (j.ne.68)) write(*,*) i,j,sd_tot(i,j)
!!$        end if
!!$     end do
!!$  end do

  ! fixing the lowest freqs (I want all my lowfreq bins intact, so do this before masking)
  bin_lowlimit = int((log10(template_freq_limit_low)-ymin)*(numbins-1)/dy) + 1
  if(myid==0) write(*,*) 'lowbin:', bin_lowlimit
  counts_tot(0:bin_lowlimit,:) = 1.d0

  ! Computing mean and stddev for each bin
  ratio_tot = ratio_tot/counts_tot
  sd_tot = sqrt(sd_tot/counts_tot - ratio_tot**2)

  mask = 1
  where(counts_tot == 0.d0) mask = 0
  ! this works - I checked

  ! Iterative smoothing of the templates
  allocate(smooth(0:numbins,ndi,3)) ! freqs,means,sds
!  allocate(smooth(numbins,ndi,3)) ! freqs,means,sds
  allocate(n(ndi)) ! Individual template lengths
  
!!$  do i=1,numbins
!!$     smooth(i,:,1) = ind2freq(i,srate,numbins)
!!$  end do

! Setting the frequencies
  smooth(0,:,1) = 0.d0
!  smooth(1,:,1) = 1.d-3
  do i=1,numbins
     smooth(i,:,1) = (10**(i*dy/numbins + ymin) + 10**((i+1)*dy/numbins + ymin))/2.d0
  end do
  smooth(:,:,2) = ratio_tot
  smooth(:,:,3) = sd_tot

  ! Setting lowest and highest freqs to 1
  smooth(0:bin_lowlimit,:,2:3) = 1.d0

  bin_highlimit = minloc(abs(smooth(:,1,1)-template_freq_limit_high))
  smooth(bin_highlimit(1):numbins,:,2) = 1.d0
  smooth(bin_highlimit(1):numbins,:,3) = 1.d0

  n = numbins  

  do iter = 1,4
     if (myid==0) write(*,*) 'Smoothing templates, iteration ', trim(itoa(iter))
     do di=1,ndi
        a = 1 
        i = 1

        ! This used to count from 2 (when I used linbins)
!        do i=2,n(di)-1,2
        do while (i<n(di))
           if (((iter==1) .and. (mask(i,di)==1)) .or. (iter>1)) then ! In 1st it we skip all the empty bins, so they are no longer present
              j=i+1
              do while ((iter==1) .and. (mask(j,di)==0))
                 if(myid==0) write(*,*) iter,di, i,j
                 j = j+1
              end do
              fact = abs(smooth(j,di,2)-smooth(i,di,2)) / sqrt(smooth(j,di,3)*smooth(i,di,3))
              if (fact>1.0) then ! copy directly
!                 if((myid==0) .and. (fact>50)) write(*,*) 'large factor:', iter, fact, di, smooth(i,di,1), smooth(i+1,di,2), smooth(i,di,2), smooth(i,di,3), smooth(i+1,di,3), sd_tot(i,di), sd_tot(i+1,di)
 
                 smooth(a:a+1,di,1) = (/smooth(i,di,1), smooth(j,di,1)/)
                 smooth(a:a+1,di,2) = (/smooth(i,di,2), smooth(j,di,2)/)
                 smooth(a:a+1,di,3) = (/smooth(i,di,3), smooth(j,di,3)/)
                 a=a+2
              else
                 smooth(a,di,1) = smooth(i,di,1) + (smooth(j,di,1)-smooth(i,di,1))/2.d0
                 smooth(a,di,2) = ( smooth(i,di,2)/smooth(i,di,3)**2 + &
                      & smooth(j,di,2)/smooth(j,di,3)**2 ) / (1.d0/smooth(i,di,3)**2 + &
                      & 1.d0/smooth(j,di,3)**2) 
                 smooth(a,di,3) = sqrt(1.d0 / (1.d0/smooth(i,di,3)**2 + 1.d0/smooth(j,di,3)**2) )
                 a = a+1
              end if

              i=j+1
           else
              i=i+1
           end if
        end do
        n(di) = a-1
     end do
     if (myid==0) write(*,*) n             
  end do

  ! Temp diode templates are currently bunk - overwrite them just so they don't make any trouble
  ! According to chisqstats the temp diodes are fine anyway. Do we even use the pol signal from the temp horns?
  smooth(:,69:ndi,2:3) = 1.d0

  ! Remove the worst outliers: for each template, calculate total mean and stddev. If value is more than 0.1sigma off, exchange it for the average of the neighbours
  ! (Visually validated. Limit is ad hoc, but the resulting templates look ok)
  do di= 1,ndi
     totmean = mean(smooth(:,di,2))
     totsd = sqrt(mean(smooth(:,di,1)**2) - totmean**2)
     if (myid==0) write(*,*) di, totmean, totsd
     do i = 2,n(di)-1
        if (abs(smooth(i,di,2)-totmean) > 0.1*totsd) then
           if (myid==0) write(*,*) di, i, abs(smooth(i,di,2)-totmean), 0.1*totsd
           a = 1; b = 1
           do while (abs(smooth(i-a,di,2)-totmean) > 0.1*totsd)
              a = a-1
           end do
           do while (abs(smooth(i+b,di,2)-totmean) > 0.1*totsd)
              b = b+1
           end do
           if (myid==0) write(*,*) a, b
           if (myid==0) write(*,*)  ' '
           smooth(i-a+1:i+b-1,di,2) = (smooth(i-a,di,2) + smooth(i+b,di,2))/2. 

        end if
     end do
  end do


!!$ This was unstable, now using the above instead.
!!$! Removing the worst of the single-bin spikes - note that threshold is ad-hoc
!!$  do di = 1,ndi
!!$     do i=1,n(di)-2
!!$        if ((abs(smooth(i,di,2)-smooth(i+1,di,2))>0.4) .and. (abs(smooth(i+1,di,2)-smooth(i+2,di,2))>0.4)) then
!!$           smooth(i+1,di,2) = (smooth(i,di,2) + smooth(i+2,di,2))/2.d0
!!$           if (myid==0) write(*,*) 'spike!', di,i
!!$        end if
!!$     end do
!!$  end do


! Median filtering
!!$  allocate(medfiltered(numbins,ndi))
!!$  do di=1,ndi
!!$     call median_filter(ratio_tot(:,di),medfiltered(:,di),numbins/100)
!!$  end do

! running average to smooth some more
  allocate(runav(0:numbins,ndi))
  do di = 1,ndi
     call average_filter(smooth(0:n(di),di,2),runav(:,di),2)
  end do


  if (myid==0) then
     ! Writing info file that filter mod can read - since the temps have varying lengths
     unit=getlun()
     dumpfile = trim(odir) // '/info.dat'
     open(unit,file=dumpfile,action='write')
     write(unit,'(2i6)') ndi, maxval(n)
     close(unit)

     ! Dumping the templates
     do di=1,ndi
        unit=getlun()
        dumpfile = trim(odir) // '/diode_' // trim(itoa(di,3)) // '.dat'
        open(unit,file=dumpfile,action='write')
        write(unit,'(i6)') n(di)
        do i=0,n(di)
           !write(unit,'(f10.6, 2g16.7)') smooth(i,di,:)
           write(unit,'(f10.6, 2g16.7)') smooth(i,di,1), runav(i,di), smooth(i,di,3)
        end do
        close(unit)
     end do

     ! Dumping the original ratios too, just in case
     do di=1,ndi
        unit=getlun()
        dumpfile = trim(odir) // '/unsmoothed/diode_' // trim(itoa(di,3)) // '.dat'
        open(unit,file=dumpfile,action='write')
        write(unit,'(i6)') numbins
        do i=0,numbins
           write(unit,'(f10.6, 2g16.7, f14.1)') (10**(i*dy/numbins + ymin) + &
                & 10**((i+1)*dy/numbins + ymin))/2.d0, ratio_tot(i,di), sd_tot(i,di), counts_tot(i,di)
        end do
        close(unit)
     end do
  end if
   
  deallocate(ratio_tot,counts_tot,sd_tot)
  deallocate(n,smooth)
!!$  deallocate(medfiltered)
  call mpi_finalize(ierr)
end program wtemplate
