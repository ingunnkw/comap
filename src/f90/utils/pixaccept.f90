! Given an acceptlist, a runlist and a target name, produce a list
! of possible pixel ganis / residual signal losses for removal of various
! files from the acceptlist.
program pixaccept
  use quiet_utils
  use quiet_Lx_mod
  use quiet_acceptlist_mod
  use quiet_mpi_mod
  implicit none

  ! Describes a single reduced pixel
  type pix_info
     integer(i4b)              :: n
     integer(i4b), allocatable :: pix(:), cnums(:)
  end type

  type scan_info
     real(dp)                  :: weight
     integer(i4b)              :: cid
     integer(i4b), allocatable :: pix(:), inds(:)
  end type

  type cut_info
     integer(i4b)              :: npix, nces
     integer(i4b), allocatable :: cnums(:)
  end type

  character(len=512)   :: parfile, afile, target_name, ofname
  type(acceptlist)     :: alist
  integer(i4b)         :: myid, nproc, err, nside, nmax, cnum, i, n
  type(quiet_ces_info) :: ces
  type(scan_info), allocatable :: scans(:)
  type(pix_info),  allocatable :: regions(:)
  type(cut_info),  allocatable :: cuts(:)
  integer(i4b),    allocatable :: pixels(:), cnums(:), pnums(:)

  if(iargc() < 2) then
     write(*,*) "syntax: pixaccept parfile ofile"
     stop
  end if

  call mpi_init(err)
  call mpi_comm_rank(MPI_COMM_WORLD, myid,  err)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, err)

  call getarg(1, parfile)
  call get_parameter(0, parfile, 'ACCEPTLIST',  par_string=afile)
  call get_parameter(0, parfile, 'TARGET_NAME', par_string=target_name)
  call get_parameter(0, parfile, 'NSIDE_OUT',   par_int   =nside)
  call get_parameter(0, parfile, 'MAXCUT_NCES', par_int   =nmax)
  call getarg(2, ofname)

  if(myid==0) call dmem("start")
  call init_detector_mod(parfile)
  call initialize_ces_mod(parfile)
  call initialize_accept_list(afile, alist)

  ! Prune to the set of interesting scans
  do cnum = 1, get_num_ces()
     call get_ces_info(cnum, ces)
     if(target_name == "any" .or. target_name == ces%object) cycle
     alist%status(:,:,cnum) = REJECTED_TARGET
     call free_ces_info(ces)
  end do
  call get_accepted_ceses(alist, cnums);            if(myid==0) call dmem("accept")
  call read_scans(cnums, nside, scans);             if(myid==0) call dmem("read")
  call select_prune(scans, nmax, pnums, pixels, n); if(myid==0) call dmem("prune")
  call regionize(scans(pnums), pixels, regions);    if(myid==0) call dmem("regionize")
  call optimize(scans(pnums), regions, cuts);       if(myid==0) call dmem("optimize")
  call output(scans(pnums), n, cuts, ofname);       if(myid==0) call dmem("output")

  call mpi_finalize(err)

contains

  ! Get all accepted scans while reducing pixels to the given nside
  subroutine read_scans(cnums, nside, scans)
     integer(i4b),   intent(in)   :: nside, cnums(:)
     type(scan_info),allocatable  :: scans(:)
     type(quiet_ces_info)         :: ces
     type(hdf_file)               :: hfile
     integer(i4b)                 :: cnum, i, j, k, n, owner, nside_file, err, nscan
     integer(i4b)                 :: asize(3), nsamp
     integer(i4b),   allocatable  :: hits(:), pixels(:)
     nscan = size(cnums)
     allocate(scans(nscan))
     ! Now fetch the exposure information from each of these
     allocate(hits(0:12*nside**2-1))
     do i = 1+myid, size(cnums), nproc
        call get_ces_info(cnums(i), ces)
        write(*,fmt="(i3,a12,a,i5,a)") myid, " scanning", " ces ", ces%cid, " (" // &
         & trim(itoa(i)) // "/" // trim(itoa(size(cnums))) // ")"
        call open_hdf_file(ces%l3file, hfile, "r")
        call get_size_hdf(hfile, "time", asize)
        nsamp = asize(1)
        scans(i)%weight = nsamp*real(count(alist%status(:,:,cnums(i))==REJECTED_NONE),dp)/size(alist%status(:,:,cnums(i)))
        call read_alloc_hdf(hfile, "pixels", pixels)
        call read_hdf(hfile, "nside",  nside_file)
        ! Reduce to requested nside
        hits = 0
        do j = 1, size(pixels)
           k = pixels(j)/(nside_file/nside)**2
           hits(k) = hits(k) + 1
        end do
        call wherei(hits>0, scans(i)%pix)
        scans(i)%cid = ces%cid
        scans(i)%pix = scans(i)%pix-1 ! 0 based
        scans(i)%weight = scans(i)%weight / size(scans(i)%pix)
        call free_ces_info(ces)
     end do
     ! Broadcast the info to everybody
     do i = 1, nscan
        owner = modulo(i-1,nproc)
        if(owner == myid) n = size(scans(i)%pix)
        call mpi_bcast(n, 1, mpi_integer, owner, mpi_comm_world, err)
        if(owner /= myid) allocate(scans(i)%pix(n))
        call mpi_bcast(scans(i)%pix, n, mpi_integer, owner, mpi_comm_world, err)
        call mpi_bcast(scans(i)%weight, 1, mpi_double_precision, owner, mpi_comm_world, err)
        call mpi_bcast(scans(i)%cid, 1, mpi_integer, owner, mpi_comm_world, err)
     end do
  end subroutine

  ! Given a set of pixels per ces, build an exposure map and find an exposure
  ! threshold so that at most nmax levels of scans are below it. Allocate the pixels
  ! array and fill it with pixels *below* the threshold.
  subroutine select_pixels(scans, npixtot, nmax, pixels)
    implicit none
    type(scan_info),intent(in) :: scans(:)
    integer(i4b),   intent(in) :: nmax
    integer(i4b),   intent(out):: npixtot
    integer(i4b),   allocatable:: pixels(:)
    integer(i4b),   allocatable:: counts(:)
    integer(i4b)               :: i, j, k, m, n, p
    n= 0
    do i = 1, size(scans)
       n= max(n,maxval(scans(i)%pix))
    end do
    allocate(counts(0:n))
    counts = 0
    do i = 1, size(scans)
       do j = 1, size(scans(i)%pix)
          p = scans(i)%pix(j)
          counts(p) = counts(p) + 1
       end do
    end do
    npixtot = count(counts>0)
    ! Now find the areas where counts is between 1 and nmax
    call wherei(counts > 0 .and. counts <= nmax, pixels)
    pixels = pixels-1 ! 0 based
    deallocate(counts)
  end subroutine


  ! Given a set of pixels per ces and a pixel set, find the indices of each of their
  ! pixels into the pixel set, removing those that fall outside. Note %pix is
  ! not changed.
  subroutine prune(scans, cnums, pixels)
    implicit none
    type(scan_info) :: scans(:)
    integer(i4b)    :: pixels(:)
    integer(i4b)    :: i, j, k, n, m, nuseful
    integer(i4b),   allocatable :: mask(:), cnums(:)
    n = maxval(pixels)
    allocate(mask(0:n))
    mask = 0
    do i = 1, size(pixels)
       mask(pixels(i)) = i
    end do
    nuseful = 0
    do i = 1, size(scans)
       k = 0
       do j = 1, size(scans(i)%pix)
          if(mask(scans(i)%pix(j)) <= 0) cycle
          k = k+1
       end do
       if(allocated(scans(i)%inds)) deallocate(scans(i)%inds)
       allocate(scans(i)%inds(k))
       k = 0
       do j = 1, size(scans(i)%pix)
          if(mask(scans(i)%pix(j)) <= 0) cycle
          k = k+1
          scans(i)%inds(k) = mask(scans(i)%pix(j))
       end do
       if(size(scans(i)%inds) > 0) nuseful = nuseful+1
    end do
    allocate(cnums(nuseful))
    j = 0
    do i = 1, size(scans)
       if(size(scans(i)%inds) == 0) cycle
       j = j+1
       cnums(j) = i
    end do
  end subroutine

  ! Like above, but make sure there are no more than nmax ceses left.
  ! Does this by estimating the number of pixels available for
  ! cutting for each ces as npot = sum_pix(1/nexp), and picks the
  ! nmax best of these.
  subroutine select_prune(scans, nmax, cnums, pixels, npixtot)
    implicit none
    type(scan_info),intent(inout) :: scans(:)
    integer(i4b),   intent(out)   :: npixtot
    integer(i4b),   intent(in)    :: nmax
    integer(i4b),   allocatable   :: pixels(:), cnums(:), cnums2(:), inds(:)
    integer(i4b),   allocatable   :: pinds(:)
    real(dp),       allocatable   :: lens(:), hits(:), scores(:)
    logical(lgt),   allocatable   :: ign(:)
    integer(i4b)                  :: nces, i, j, k, n, m, nleft
    call select_pixels(scans, npixtot, nmax, pixels)
    call prune(scans, cnums2, pixels)
    n = size(cnums2)
    ! Build an exposure map
    allocate(hits(size(pixels)))
    hits = 0
    do i = 1, n
       k = cnums2(i)
       do j = 1, size(scans(k)%inds)
          hits(scans(k)%inds(j)) = hits(scans(k)%inds(j))+scans(k)%weight
       end do
    end do
    ! And give each scan a score (negative to get the best ones first)
    allocate(scores(n), inds(n))
    scores = 0
    do i = 1, n
       k = cnums2(i)
       do j = 1, size(scans(k)%inds)
          scores(i) = scores(i) - scans(k)%weight/hits(scans(k)%inds(j))
       end do
    end do
    inds = irange(n)
    ! And pick out the best ones
    call quicksort(inds, scores)
    m = min(n, nmax)
    ! We must now update the pixel information for the remaining scans,
    ! as pixels with hits from the ignored ceses will not be removable.
    allocate(ign(size(pixels)))
    ign = .false.
    do i = m+1, n
       k = cnums2(inds(i))
       do j = 1, size(scans(k)%inds)
          ign(scans(k)%inds(j)) = .true.
       end do
    end do
    do i = 1, size(scans)
       allocate(pinds(size(scans(i)%inds)))
       pinds = scans(i)%inds
       deallocate(scans(i)%inds)
       nleft = count(.not. ign(pinds))
       allocate(scans(i)%inds(nleft))
       j = 0
       do k = 1, size(pinds)
          if(ign(pinds(k))) cycle
          j = j+1
          scans(i)%inds(j) = pinds(k)
       end do
       deallocate(pinds)
    end do
    ! And create our final cnums, the set of ceses to consider
    allocate(cnums(m))
    cnums = cnums2(inds(:m))
    deallocate(hits, inds, cnums2, ign, scores)
  end subroutine

  ! Given a set of scans with active inds, reduce to a minimal set of
  ! homogeneous regions, where each region is hit by the same ceses.
  ! This is a pretty heavy operation, scaling as O(npix**2). With
  ! hash_maps, an O(npix) version would be possible. Another brute
  ! force implementation is O(npix*nreg), which should be better,
  ! since the number of regions must be smaller than the number of pixels.
  ! I'll try that one now.
  ! Only send in the scans you want to include in the analysis,
  ! i.e. those indicated by the cnums array from prune. The region
  ! cnum array will correspond to indices in scans.
  subroutine regionize(scans, pixels, regions)
    implicit none
    type(scan_info),     intent(in)   :: scans(:)
    integer(i4b),        intent(in)   :: pixels(:)
    type(pix_info),      allocatable  :: regions(:)
    integer(i4b),        allocatable  :: pixces(:,:), pixn(:)
    type(pix_info),      allocatable  :: regs(:)
    integer(i4b)                      :: i, j, k, m, n, nces, npix, p, regn
    ! Build up a list of the (sorted) ceses for each pixel
    nces = size(scans)
    npix = size(pixels)
    allocate(pixces(nces,npix),pixn(npix))
    pixn = 0
    do i = 1, nces
       do j = 1, size(scans(i)%inds)
          k = scans(i)%inds(j)
          pixn(k) = pixn(k)+1
          pixces(pixn(k),k) = i
       end do
    end do
    ! For each pixel, find which region it belongs to O(nreg) (possibly
    ! creating a new one) and add it to it.
    allocate(regs(npix))
    regn = 0
    do p = 1, npix
       if(pixn(p) == 0) cycle ! Pixel not hit at all?
       do i = 1, regn
          if(size(regs(i)%cnums) /= pixn(p)) cycle
          if(any(regs(i)%cnums /= pixces(:pixn(p),p))) cycle
          exit
       end do
       if(i > regn) then
          regn = i
          allocate(regs(i)%cnums(pixn(p)))
          regs(i)%cnums = pixces(:pixn(p),p)
          allocate(regs(i)%pix(npix))
          regs(i)%n = 0
       end if
       regs(i)%n = regs(i)%n + 1
       regs(i)%pix(regs(i)%n) = pixels(p)
    end do
    ! Ok, we have found all regions. Copy over the necessary parts to the
    ! output array
    allocate(regions(regn))
    do i = 1, regn
       regions(i)%n = regs(i)%n
       allocate(regions(i)%pix(regs(i)%n))
       regions(i)%pix = regs(i)%pix(:regs(i)%n)
       allocate(regions(i)%cnums(size(regs(i)%cnums)))
       regions(i)%cnums = regs(i)%cnums
    end do
    ! And clean up
    do i = 1, regn
       deallocate(regs(i)%pix, regs(i)%cnums)
    end do
    deallocate(regs, pixces, pixn)
  end subroutine

  ! For each number of sacrificed ceses, find the maximal number of
  ! freed pixels and which ceses this corresponds to.
  ! This one is extremely heavy, scaling as O(2^nces). However, for
  ! realistic cases nces is pretty small.
  subroutine optimize(scans, regions, cuts)
    implicit none
    type(scan_info),   intent(in) :: scans(:)
    type(pix_info),    intent(in) :: regions(:)
    type(cut_info),    allocatable:: cuts(:)
    integer(i4b),      allocatable:: regs(:), best(:,:,:)
    integer(i4b)                  :: i, j, k, m, n, nreg, nscan, err, which(32), pos(1)
    integer(i4b)                  :: mask
    ! Encode the region information in binary for efficiency
    nreg = size(regions)
    nscan= size(scans)
    allocate(regs(nreg), best(2,0:nscan,0:nproc-1))
    do i = 1, nreg
       regs(i) = 0
       do j = 1, size(regions(i)%cnums)
          regs(i) = ior(regs(i),2**(regions(i)%cnums(j)-1))
       end do
    end do
    ! Now loop through every possibility, and count the number of pixels for each.
    ! This is a very long loop, so its content must be fast!
    best  = 0
    do i = myid, 2**nscan-1, nproc
       ! How many scans do we remove here?
       m = count_set_bits(i)
       ! Count number of pixels
       n = 0
       do j = 1, nreg
          if(iand(i,regs(j))==regs(j)) n = n + size(regions(j)%pix)
       end do
       if(n > best(1,m,myid)) then
          best(1,m,myid) = n
          best(2,m,myid) = i
       end if
    end do
    ! And gather the results
    do i = 0, nproc-1
       call mpi_bcast(best(:,:,i), size(best(:,:,i)), mpi_integer, i, mpi_comm_world, err)
    end do
    ! And output
    allocate(cuts(nscan+1))
    do i = 0, nscan
       n = 0
       pos = maxloc(best(1,i,:))
       k = pos(1)-1
       ! Expand binary representation
       mask = best(2,i,k)
       do j = 1, 32
          if(iand(mask,1)==1) then
             n=n+1
             which(n) = j
          end if
          mask=mask/2
       end do
       allocate(cuts(i+1)%cnums(n))
       cuts(i+1)%cnums = which(:n)
       cuts(i+1)%nces  = n
       cuts(i+1)%npix  = best(1,i,k)
    end do
    deallocate(regs, best)
  end subroutine

  ! Output the pixels saved per cut
  subroutine output(scans, npixtot, cuts, ofname)
    implicit none
    type(scan_info),   intent(in)    :: scans(:)
    integer(i4b),      intent(in)    :: npixtot
    type(cut_info),    intent(in)    :: cuts(:)
    character(len=*),  intent(in)    :: ofname
    integer(i4b)                     :: i, j, k, m, n, unit
    if(myid /= 0) return
    unit = getlun()
    open(unit,file=ofname)
    write(unit,'(a)') 'n: ' // trim(itoa(npixtot))
    do i = 1, size(cuts)
       write(unit,'(i3,i8)',advance="no") cuts(i)%nces, cuts(i)%npix
       do j = 1, cuts(i)%nces
          write(unit,'(i5)',advance="no") scans(cuts(i)%cnums(j))%cid
       end do
       write(unit,*)
    end do
    close(unit)
  end subroutine

end program
