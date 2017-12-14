! Scan validation tool.
!Run through each scan, and determine whether a frequency-detector-scan 
!is acceptable for map making

program scan_validate
  use quiet_utils
  use quiet_fileutils
  !use quiet_system_mod
  use quiet_fft_mod
  use comap_scan_mod
  use comap_detector_mod
  use comap_frequency_mod
  !use quiet_filter_mod
  use quiet_mpi_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_shared_output_mod
  use comap_acceptlist_mod
  use comap_lx_mod
  use scan_validate_mod
  use comap_patch_mod
  implicit none

  type info_struct
     type(shared_ofile) :: statfile, acceptfile
     integer(i4b)       :: id, nproc
     logical(lgt)       :: sigma0, alpha, f_knee, corr, validate
  end type info_struct

  character(len=512)    :: parfile, odir, lockfile, acceptlist_in, acceptlist_out
  character(len=512)    :: statfile, runlist_file, l3dir, summary_file, outfile, donefile
  integer(i4b)          :: ierr, sid, ndet, i, j, nscan, snum, det, nfreq, unit, debug
  integer(i4b)          :: firstscan, lastscan
  logical(lgt)          :: accepted, apply_scan_det_cuts, bright_object, use_templates
  logical(lgt)          :: use_filter_highpass, use_filter_lowpass, baseline
  type(task_list)       :: tasks
  type(info_struct)     :: info
  type(acceptlist)      :: alist
  type(comap_runlist)   :: runlist 
  type(comap_scan_info) :: scan
  type(lx_struct)       :: data
  type(hdf_file)        :: file
  type(validation_struct), dimension(:,:),   allocatable :: stats
  integer(i4b), dimension(:,:),   allocatable :: my_status, status_tot
  real(sp),     dimension(:,:),   allocatable :: powspecs, tods, templates
  real(dp),     dimension(:),     allocatable :: F_filter, N_fft, scanfreqs, sftmp 
  real(dp),     dimension(:,:),   allocatable :: ftmp
  real(sp),     dimension(:),     allocatable :: stmp
  real(sp),     dimension(:,:),   allocatable :: scan_stats, scan_buffer, dtmp
  real(sp),     dimension(:,:,:), allocatable :: det_stats, det_buffer
  complex(spc), dimension(:,:),   allocatable :: ffts
  logical(lgt), dimension(:),     allocatable :: exists, existbuf

  call getarg(1, parfile)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, info%id, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, info%nproc, ierr)
  call dset(id=info%id)
  call print_host_mapping

  call get_parameter(0, parfile, 'OUTPUT_DIR',            par_string=odir)
  call get_parameter(0, parfile, 'ACCEPT_LIST_INPUT',     par_string=acceptlist_in)
  call get_parameter(0, parfile, 'RUNLIST',               par_string=runlist_file)
  call get_parameter(0, parfile, 'LEVEL3_DIR',            par_string=l3dir)
  call get_parameter(0, parfile, 'APPLY_SCAN_DET_CUTS', par_lgt=apply_scan_det_cuts)
  call get_parameter(0, parfile, 'DEBUG',                 par_int=debug)
  call get_parameter(0, parfile, 'FIRST_SCAN',            par_int=firstscan)
  call get_parameter(0, parfile, 'LAST_SCAN',             par_int=lastscan)
  call get_parameter(0, parfile, 'USE_TEMPLATES',         par_lgt=use_templates)
  call get_parameter(0, parfile, 'APPLY_HIGHPASS_FILTER', par_lgt=use_filter_highpass)
  call get_parameter(0, parfile, 'APPLY_LOWPASS_FILTER',  par_lgt=use_filter_lowpass)
  call get_parameter(0, parfile, 'MAKE_BASELINE_ACCLIST', par_lgt=baseline)

  if (baseline .and. (info%id == 0)) then
     write(*,*) 'Building baseline acceptlist. All but the initial obvious cuts disabled.'
     write(*,*) '***Note: summary.txt and scan_statistics.txt will not be valid.'
  end if

  acceptlist_out = trim(odir) // '/accept.txt'
  summary_file   = trim(odir) // '/summary.txt'
  lockfile       = trim(odir) // "/lock.dat"
  statfile       = trim(odir) // "/scan_statistics.txt"
  !call dset(level=debug)

  !call dmem("init")

  call mkdirs(trim(lockfile), .true.)
  call open_shared_ofile(info%statfile, trim(statfile), MPI_COMM_WORLD)
  call open_shared_ofile(info%acceptfile, trim(acceptlist_out), MPI_COMM_WORLD)

  call initialize_comap_patch_mod(parfile)
  call dmem("patch mod")
  call initialize_scan_mod(parfile)
  call dmem("scan mod")
  call initialize_detector_mod(parfile)
  call dmem("detector mod")
  !call initialize_filter_mod(parfile)
  !call dmem("filter mod")

  ndet  = get_num_dets()
  nscan = get_num_scans()
  !if(info%id == 0) call write_stat_headers(info)


  ! First scan through scans: Build up detector-scan-statistics
  allocate(dtmp(ndet,NUM_DET_STATS),stmp(STAT_NUM),ftmp(ndet,5))
  allocate(det_stats(ndet,DET_STAT_MAX,nscan), det_buffer(ndet,DET_STAT_MAX,nscan))
  allocate(scan_stats(STAT_MAX,nscan), scan_buffer(STAT_MAX,nscan))
  allocate(exists(nscan),existbuf(nscan))
  allocate(scanfreqs(nscan),sftmp(nscan))
  det_stats = 0.d0
  scan_stats   = 0.d0
  exists      = .false.
  sftmp       = 0.d0
  if (apply_scan_det_cuts) then
     call init_task_list(tasks, lockfile, nscan, MPI_COMM_WORLD)
     do while(get_next_task(tasks, snum))     
        call get_scan_info(snum, scan)                   
        write(*,fmt="(i3,a,i4,a)") info%id, " scanning scan ", scan%sid, " (" // trim(itoa(snum)) // "/" // &
             & trim(itoa(get_num_scans())) // ")"
        inquire(file=trim(scan%l3file), exist=exists(snum))
        if (.not. exists(snum)) cycle
        call open_hdf_file(scan%l3file, file, 'r')
        call read_hdf(file, 'stats',       stmp)
        call read_hdf(file, 'det_stats', dtmp(:,:))
        call read_hdf(file, 'sigma0',      dtmp(:,DET_STAT_SIGMA0))
        call read_hdf(file, 'alpha',       dtmp(:,DET_STAT_ALPHA))
        call read_hdf(file, 'fknee',       dtmp(:,DET_STAT_FKNEE))
        call read_hdf(file, 'filter_par',  ftmp)
        dtmp(:,DET_STAT_NU_LOW)  = ftmp(:,FILTER_LOW_NU)
        dtmp(:,DET_STAT_NU_HIGH) = ftmp(:,FILTER_HIGH_NU_SCAN)
        dtmp(:,DET_STAT_AZORDER) = ftmp(:,FILTER_AZ_ORDER)
        call read_hdf(file, 'scanfreq',   sftmp(snum))
        call close_hdf_file(file)
        scan_stats(1:STAT_MAX,snum) = stmp
        det_stats(:,1:DET_STAT_MAX,snum) = dtmp
        call dmem("scanning " // trim(itoa(scan%sid)), 2)
     end do
     call mpi_allreduce(det_stats, det_buffer,     size(det_buffer),     mpi_st, MPI_SUM, mpi_comm_world, ierr)
     call mpi_allreduce(scan_stats,   scan_buffer, size(scan_buffer), mpi_st, MPI_SUM, mpi_comm_world, ierr)
     call mpi_allreduce(exists,      existbuf,   size(exists),     mpi_logical,          MPI_LOR, mpi_comm_world, ierr)
     call mpi_allreduce(sftmp,       scanfreqs,  size(scanfreqs),  mpi_st, MPI_SUM, mpi_comm_world, ierr)
     det_stats = det_buffer
     scan_stats   = scan_buffer
     exists      = existbuf
     call free_task_list(tasks)
  end if
  deallocate(det_buffer, scan_buffer, existbuf, stmp, dtmp, sftmp)
  call dmem("scan done")

  ! Initialize the validation module, including det specific values
  call initialize_scan_validate_mod(parfile); call dmem("validate mod")

  ! Initialize stats structure
  allocate(stats(ndet,nscan))
 ! call initialize_validation_structs(stats)

  ! ********************************************************************
  !                           Apply cuts 
  ! ********************************************************************
  call dmem("cuts start")

  ! Check for dead dets
  do i = 1, ndet
     if (.not. is_alive(i)) then
        det = (i-1)/get_num_dets()
!        stats(i,:)%status = ior(stats(i,:)%status, status_baddi)
     end if
  end do
  call dmem("dead dets")

  ! Missing files
  do i = 1, nscan
     if(.not. exists(i)) stats(:,i)%status = ior(stats(:,i)%status, status_nofile)
  end do
  call dmem("missing files")

  ! Apply static cuts. Since the purpose of the static cuts is to cut
  ! extra scanes than the normal automatic cuts, it does not make sense
  ! that every scan-det needs to be mentioned in the input accept list
  ! in order not to be rejected. The default for missing scanes in this
  ! list is therefore to accept them (and subject them to the normal
  ! set of automatic cuts).
  call initialize_accept_list(acceptlist_in, alist, REJECTED_NONE)
  do i = 1, nscan
     call get_scan_info(i, scan)                   
     do j = 1, ndet
        det  = modulo(j-1,get_num_dets())
        !if ((.not. is_accepted(alist,scan%sid,mod,di) .or. scan%sid < firstscan .or. &
        !     & scan%sid > lastscan) .and. stats(j,i)%status == status_ok) then
        !   stats(j,i)%status = ior(stats(j,i)%status, status_static)
        !end if
     end do
  end do
  call deallocate_acceptlist(alist)
  call dmem("input cuts")

  ! Other simple static cuts:
  !call validate_static(status_static, stats)
  !call dmem("static cuts")

  ! Obvious non-varying cuts
  !call validate_elementary(det_stats, status_elementary, stats)
  !call validate_gain(det_stats, status_gain, stats)

  !call dmem("elementary cuts")

  !call validate_jump(det_stats(:,DET_STAT_JUMP,:), status_jump, stats) 
  ! Apply APEX cuts
  !call validate_APEX(scan_stats, status_apex, stats)
  !call dmem("simple cuts")

  ! For gc/gb: Skip all but the completely obvious cuts to get a baseline acceptlist
  !if (baseline) goto 41

  ! Apply filter cuts
  !call validate_filter_params(det_stats(:,DET_STAT_NU_LOW,:), &
  !     & det_stats(:,DET_STAT_NU_HIGH,:), scanfreqs, status_filter, stats)

  ! Make scan-det cuts only if we have more than N SCANs
  if (nscan > 1 .and. apply_scan_det_cuts) then
     call initialize_det_scan_cuts(det_stats)
     call dmem("init det cuts")
     !call validate_sigma0(det_stats(:,DET_STAT_SIGMA0,:),    status_sigma0, stats)
     !call dmem("sigma0 cut")
     !call validate_fknee(det_stats(:,DET_STAT_FKNEE,:),      status_fknee,  stats)
     !call dmem("fknee cut")
  end if

  ! Proces all scans
41  call init_task_list(tasks, lockfile, nscan, MPI_COMM_WORLD)
  do while(get_next_task(tasks, snum))     
     if(.not. exists(snum)) cycle
     call dmem('scan begin')
     call get_scan_info(snum, scan)

     ! Do not cut on demod for strong objects
     bright_object = trim(scan%object) == 'moon'

     write(*,fmt="(i3,a,i4,a)") info%id, " proscansing scan ", scan%sid, " (" // trim(itoa(snum)) // "/" // &
          & trim(itoa(get_num_scans())) // ")"

     call read_l3_file(scan%l3file, data)            ; call dmem('read_l3_file')
     call calc_fourier(data, ffts, powspecs)        ; call dmem("fourier")
     nfreq = size(powspecs,1)
     ! TMR
     !if (use_templates) then
     !   call get_templates(templates, nfreq, data%samprate)
     !   powspecs = powspecs/templates ! powspecs(nf,ndet)
     !end if

     allocate(tods(size(data%tod,1),size(data%tod,2)))
     allocate(F_filter(nfreq),N_fft(nfreq))

     !if (baseline) goto 42     

     ! Is it faster/better to do this one det at a time?
     !call compute_filtered_tod(data, ffts, tods, use_filter_lowpass, use_filter_highpass)    ; call dmem("filter_tod")
     ! Check statistics for each det
     do i = 1, ndet
        det = (i-1)/get_num_dets()

        ! Initialize statistic structures
        !call update_current_stat(mod, di, stats(i,snum), data, det_stats(:,:,snum))
        !if (.not. is_active(stats(i,snum)%status)) cycle

        ! Apply FFT cuts
        !call get_N_filter(data%samprate, data%sigma0(i), data%fknee(i), &
        ! & data%alpha(i), 1.d0, N_fft, use_filter_highpass)

        ! Now using individually tuned filter parameters
        F_filter = 1.d0
        !if(use_filter_lowpass)  call apodize_filter_fft(.true., size(tods,1), data%samprate, &
        !     & data%filter_par(i,FILTER_LOW_NU), data%filter_par(i,FILTER_LOW_ALPHA), .false., F_filter)
        !if(use_filter_highpass) call apodize_filter_fft(.true.,size(tods,1), data%samprate, &
        !     & data%filter_par(i,FILTER_HIGH_NU_SCAN)*data%scanfreq, data%filter_par(i,FILTER_HIGH_ALPHA), .true., F_filter)

        !call validate_fft(is_polarization_module(mod), data%samprate, &
        ! & data%scanfreq, F_filter, N_fft, powspecs(:,i), stats(i,snum))

        ! Apply TOD cuts
        !if (use_templates) then
        !   call validate_tod(is_polarization_module(mod), (sum(F_filter**2*N_fft*templates(:,i))) / &
        !        & (nfreq), data%orig_point(1,:), tods(:,i), stats(i,snum))
        !else
        !   call validate_tod(is_polarization_module(mod), (sum(F_filter**2*N_fft)) / &
        !        & (nfreq), data%orig_point(1,:), tods(:,i), stats(i,snum))
        !end if
     end do

     ! FIXME: Validate_corr must be generalized
     !call validate_corr(data%corr, status_corr, stats(:,snum)); call dmem("validate corr")
     ! If we have a bright object, ignore certain cuts
     !if(bright_object) stats(:,snum)%status = iand(stats(:,snum)%status, not(mask_bright))
     !call validate_global_scan(data, stats(:,snum));           ; call dmem("global")
!42   !call output_statistic_file(scan%sid, stats(:,snum), info) ; call dmem('output_stat_file')
     !call output_accept_file(scan%sid, stats(:,snum), info)    ; call dmem('output_accept_file')

     call free_lx_struct(data)
     deallocate(ffts, powspecs, tods, F_filter, N_fft)
     if(allocated(templates)) deallocate(templates)
  end do
  
  call close_shared_ofile(info%statfile) 
  call close_shared_ofile(info%acceptfile) 

  ! Generate summary statistics
  allocate(my_status(ndet,nscan), status_tot(ndet,nscan))
  do i = 1, nscan
     do j = 1, ndet
        my_status(j,i) = stats(j,i)%status
     end do
  end do

  call mpi_reduce(my_status, status_tot, size(my_status), mpi_integer, MPI_MAX, 0, mpi_comm_world, ierr)

  if (info%id == 0) then
     !call output_season_summary(status_tot, summary_file)
  end if
  deallocate(my_status, status_tot)

  deallocate(stats, det_stats, scanfreqs)
  call free_task_list(tasks)
  call dmem("finished")
  call mpi_finalize(ierr)

  if (info%id == 0) then
     write(*,*)
     write(*,*) 'scan_validate completed sucscansfully'
     write(*,*)
  end if


contains
!!$
!!$  ! Routines for statistics output
!!$  subroutine output_statistic_file(sid, stats, info)
!!$    implicit none
!!$    type(info_struct),                     intent(in) :: info
!!$    integer(i4b),                          intent(in) :: sid
!!$    type(validation_struct), dimension(:), intent(in) :: stats
!!$
!!$    integer(i4b)        :: i, j, accept, ndet, snum, mod, di
!!$    character(len=10000) :: line
!!$    character(len=10)    :: object
!!$    type(quiet_scan_info) :: scan
!!$
!!$    ndet        = size(stats)
!!$    call get_scan_info(lookup_scan(sid), scan)
!!$    object      = scan%object
!!$
!!$    do i = 1, ndet
!!$       mod = (i-1)/get_num_dets()
!!$       di  = modulo(i-1,get_num_dets())
!!$       write(line,fmt='(i8,a13,i4,i6,z9,i10,f7.2,e12.4,4f10.3,i6,15f12.2)') &
!!$            & sid, adjustr(object), mod, di, stats(i)%status, stats(i)%status, &
!!$            & stats(i)%scan_accept_ratio, stats(i)%sigma0, stats(i)%alpha, &
!!$            & min(stats(i)%fknee,100.d0), stats(i)%nu_high, &
!!$            & stats(i)%nu_low, nint(stats(i)%azorder), stats(i)%full_fft_chisq, &
!!$            & stats(i)%scanfreq_fft_chisq, stats(i)%band_fft_chisq(1:6), &
!!$            & stats(i)%tod_chisq, stats(i)%tod_absmax, stats(i)%tod_chi_az, &
!!$            & stats(i)%typeb, 1d4*stats(i)%wthr1, 1d4*stats(i)%wthr2, stats(i)%jump
!!$       call write_shared_ofile(info%statfile, trim(line))
!!$    end do
!!$    call flush_shared_ofile(info%statfile)
!!$  end subroutine output_statistic_file
!!$
!!$  subroutine output_accept_file(sid, stats, info)
!!$    implicit none
!!$    type(info_struct),                     intent(in) :: info
!!$    integer(i4b),                          intent(in) :: sid
!!$    type(validation_struct), dimension(:), intent(in) :: stats
!!$
!!$    integer(i4b)        :: i, j, mod, di
!!$    character(len=10000) :: line
!!$
!!$    ! Write line to file
!!$    do i = 1, size(stats)
!!$       mod = (i-1)/get_num_dets()
!!$       di  = modulo(i-1,get_num_dets())
!!$       if (stats(i)%status == status_ok) then
!!$          write(line,fmt='(i7,i5,i2,i3)') sid, mod, di, 1
!!$       else 
!!$          write(line,fmt='(i7,i5,i2,i3)') sid, mod, di, 0
!!$       end if
!!$       call write_shared_ofile(info%acceptfile, trim(line))
!!$    end do
!!$    call flush_shared_ofile(info%acceptfile)
!!$  end subroutine output_accept_file
!!$
!!$
!!$  subroutine write_stat_headers(info)
!!$    implicit none
!!$    type(info_struct) :: info
!!$    integer(i4b) :: unit
!!$    call write_shared_ofile(info%statfile, '#    SCAN_id    Object  Mod   Di  StatHex StatDec' // &
!!$         & '   R      Sigma0      Alpha     Fknee   nu_high  nu_low   az_order   0_12.5Hz' //&
!!$         & '       Scan        0_0.2Hz       CMB       10_12.5Hz     1.2Hz    FFT_spike' // &
!!$         & '1.0Hz    TOD_chisq  TOD_abs    chi_az   type_B    Wthr_10s   Wthr_30s   Jump')
!!$    call flush_shared_ofile(info%statfile)
!!$  end subroutine write_stat_headers
!!$
!!$  subroutine compute_filtered_tod(data, ffts, tods, lowpass, highpass)
!!$    implicit none
!!$
!!$    type(lx_struct)                           :: data
!!$    complex(spc), dimension(:,:), intent(in)  :: ffts
!!$    real(spc),    dimension(:,:), intent(out) :: tods
!!$    logical(lgt),                 intent(in)  :: lowpass, highpass
!!$    
!!$    integer(i4b) :: i, ndet, nsamp, nfreq
!!$    real(dp),     allocatable, dimension(:)   :: F, sqrt_invN
!!$    complex(spc), allocatable, dimension(:,:) :: ffts_int
!!$
!!$    nsamp   = size(tods,1)
!!$    ndet     = size(tods,2)
!!$    nfreq   = size(ffts,1)
!!$
!!$    allocate(F(nfreq), sqrt_invN(nfreq), ffts_int(nfreq,ndet))
!!$
!!$    do i = 1, ndet
!!$       F = 1.d0
!!$       if(lowpass)  call apodize_filter_fft(.true., nsamp, data%samprate, &
!!$            & data%filter_par(i,FILTER_LOW_NU), data%filter_par(i,FILTER_LOW_ALPHA), .false., F)
!!$       if(highpass) call apodize_filter_fft(.true., nsamp, data%samprate, &
!!$            & data%filter_par(i,FILTER_HIGH_NU_SCAN)*data%scanfreq, data%filter_par(i,FILTER_HIGH_ALPHA), .true., F)
!!$
!!$       ! TMR: Why do we do this?
!!$       sqrt_invN     = 1.d0
!!$       ffts_int(:,i) = ffts(:,i) * F * sqrt_invN
!!$    end do
!!$    call fft_multi(tods, ffts_int, -1)
!!$
!!$    deallocate(F, sqrt_invN, ffts_int)
!!$
!!$  end subroutine compute_filtered_tod
!!$
!!$  subroutine output_season_summary(stats, summary_file)
!!$    implicit none
!!$
!!$    integer(i4b),            dimension(:,:), intent(in) :: stats
!!$    character(len=*),                        intent(in) :: summary_file
!!$
!!$    integer(i4b) :: i, j, k, l, m, unit, nscan, typ, nper
!!$    integer(i8b) :: fields
!!$    logical(lgt) :: pol, mask(size(stats,1))
!!$    integer(i4b) :: cut_count(num_status), cum_count(num_status), numtot, numtot_pruned
!!$    type(filter_params) :: filteropts
!!$
!!$    nscan      = size(stats,2)
!!$    call get_default_filter(filteropts)
!!$
!!$    unit = getlun()
!!$    open(unit, file=trim(summary_file))
!!$    write(unit,*) '******** SCAN VALIDATION SUMMARY ***********'
!!$    write(unit,*)
!!$
!!$    do typ = 1, 2
!!$       pol = typ == 1
!!$
!!$       write(unit,*)
!!$       if (pol) then
!!$          write(unit,*) '----- Polarization -----'
!!$          mask = abs(quiet_dets%stokes(2)) > abs(quiet_dets%stokes(1))
!!$       else
!!$          write(unit,*) '----- Temperature -----'
!!$          mask = abs(quiet_dets%stokes(2)) <= abs(quiet_dets%stokes(1))
!!$       end if
!!$       nper = count(mask)
!!$
!!$       do i = 1, size(patches)
!!$          ! What scan-dets in this patch are cut by which cuts?
!!$          cut_count = 0
!!$          cum_count = 0
!!$          numtot    = 0
!!$
!!$          do j = 1, nscan
!!$             call get_scan_info(j,scan)
!!$             if(patches(i)%name /= scan%object) cycle
!!$             numtot = numtot  + nper
!!$             fields = 0
!!$             do k = 1, num_status
!!$                fields = ior(fields,ishft(1,k-1))
!!$                cut_count(k) = cut_count(k) + count(mask .and. btest(stats(:,j),k-1))
!!$                cum_count(k) = cum_count(k) + count(mask .and. iand (stats(:,j),fields)==0)
!!$             end do
!!$          end do
!!$
!!$          ! Calculating accept rates relative to proper data only - leaving out 3 first cuts (bad dets,
!!$          ! missing files and static cuts (various garbage)) from the accounting 
!!$          numtot_pruned = cum_count(3)
!!$
!!$          ! Unexposed patch
!!$          if(numtot == 0) cycle
!!$
!!$          ! And output
!!$          write(unit,*)
!!$          write(unit,*) trim(patches(i)%name)
!!$
!!$          1 format('    Total number of SCAN-dets ', a20,' = ',i8,', ',f7.2,'%')
!!$          2 format('    Number of SCAN-dets after ', a20,' = ',i8,', ',f7.2,'%',', ',f7.2,'%')
!!$          3 format('    Total number of viable SCAN-dets ', a20,' = ',i8,', ',f7.2,'%')
!!$          write(unit,1) '', numtot, 100.0
!!$          do j = 1, 2
!!$             write(unit,2), status_desc(j), cum_count(j), 100.0*cum_count(j)/numtot, 100.0*cut_count(j)/numtot
!!$          end do
!!$
!!$          write(unit,3) '', numtot_pruned, 100.0
!!$          do j = 3, num_status
!!$             write(unit,2), status_desc(j), cum_count(j), 100.0*cum_count(j)/numtot_pruned, 100.0*cut_count(j)/numtot_pruned
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    write(unit,*)
!!$    write(unit,*) ' ----- Hex key -----'
!!$    do i = 1, num_status
!!$       write(unit,'(z6,a,a)') ishft(1,i-1), " for ", trim(status_desc(i))
!!$    end do
!!$
!!$    ! Output summary of cut parameters
!!$    write(unit,*)
!!$    write(unit,*)
!!$    write(unit,*) ' ----- Filter parameters -----'
!!$    write(unit,*)
!!$    write(unit,*) '  APPLY_LOWPASS_FILTER               = ', filteropts%apply_lowpass
!!$    write(unit,*) '  APPLY_HIGHPASS_FILTER              = ', filteropts%apply_highpass
!!$    write(unit,*) '  NU_HIGHPASS (scanfreq)             = ', real(filteropts%nu_highpass,sp)
!!$    write(unit,*) '  ALPHA_HIGHPASS                     = ', real(filteropts%alpha_highpass,sp)
!!$    write(unit,*) '  NU_LOWPASS (Hz)                    = ', real(filteropts%nu_lowpass,sp)
!!$    write(unit,*) '  ALPHA_LOWPASS                      = ', real(filteropts%alpha_lowpass,sp)
!!$    write(unit,*)
!!$    write(unit,*)
!!$    write(unit,*) ' ----- Cut definitions -----'
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 1 -- dead dets'
!!$    write(unit,*) '    Taken from focalplane definition file'
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 2 -- missing file cuts'
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 3 -- static cuts'
!!$    write(unit,*) '    Taken from input accept list'
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 4 -- elementary cuts'
!!$    write(unit,*) '    Sanity checks'
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 5 -- correlation cuts'
!!$    write(unit,*) '    Are the correlations invertible?'
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 6 -- APEX cuts'
!!$    write(unit,*) '    APEX_MAX_PWV (mm PWV)             = ', real(apex_max_pwv,sp)
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 7 -- Type-B cuts'
!!$    write(unit,*) '    TYPE_B_CHISQ_THRESHOLD (sigma)    = ', real(typeb_threshold,sp)
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 8 -- Weather cuts'
!!$    write(unit,*) '    WEATHER_CHISQ_THRESHOLD (sigma)   = ', real(weather_threshold,sp)
!!$    write(unit,*) '    WEATHER_SCAN_FRACTION (fraction)   = ', real(weather_scan_fraction,sp)
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 9 -- sigma0 cuts'
!!$    write(unit,*) '    SIGMA0_THRESHOLD (sigma)          = ', real(sigma0_threshold,sp)
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 10 -- jump cuts'
!!$    write(unit,*) '    JUMP_FINDER_FRACTIONAL_JUMP_THRESHO= ', real(jump_finder_fractional_jump_threshold,sp)
!!$    write(unit,*) 
!!$    write(unit,*) ' Level 11 -- tod cuts'
!!$    write(unit,*) '    TOD_CHISQ_THRESHOLD (sigma)       = ', real(tod_chisq_threshold,sp)
!!$    write(unit,*) '    TOD_ABSMAX_THRESHOLD (sigma)      = ', real(tod_absmax_threshold,sp)
!!$    write(unit,*) '    TOD_AZ_BINSIZE (degrees)          = ', real(tod_az_binsize,sp)
!!$    write(unit,*) '    TOD_AZ_MAX_CHISQ (sigma)          = ', real(tod_az_max_chisq,sp)
!!$    write(unit,*)
!!$    write(unit,*) ' Level 12 -- fft cuts'
!!$    write(unit,*) '    FFT_CHISQ_THRESHOLD (CMB; sigma)  = ', real(fft_chisq_threshold,sp)
!!$    write(unit,*) '    FFT_CHISQ_1_2HZ_THRESHOLD (sigma) = ', real(fft_chisq_1_2Hz_threshold,sp)
!!$    write(unit,*) '    FFT_CHISQ_SPIKE_THRESHOLD (sigma) = ', real(fft_spike_threshold,sp)
!!$    write(unit,*) '    FFT_CHISQ_SCAN_THRESHOLD (sigma)  = ', real(fft_chisq_scan_threshold,sp)
!!$    write(unit,*) '    FFT_CHISQ_LOW_THRESHOLD (sigma)   = ', real(fft_chisq_low_threshold,sp)
!!$    write(unit,*) '    FFT_CHISQ_HIGH_THRESHOLD (sigma)  = ', real(fft_chisq_high_threshold,sp)
!!$    write(unit,*) ' Level 13 -- sun cuts'
!!$    write(unit,*) '    Based on sidelobe map files'
!!$    write(unit,*)
!!$    write(unit,*) ' Level 14 -- accept ratio cut'
!!$    write(unit,*) '    Rejects if the det fraction rejected'
!!$    write(unit,*) '    due to a cause which should logically'
!!$    write(unit,*) '    affect all dets, exceeds:'
!!$    write(unit,*) '    SCAN_MIN_DET_ACCEPT_RATIO        = ', real(min_det_accept_ratio,sp)
!!$    write(unit,*) 
!!$    close(unit)
!!$  end subroutine output_season_summary
!!$
  subroutine calc_fourier(data, ffts, powspecs)
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: ndet, i, nf
    complex(spc), dimension(:,:), allocatable :: ffts
    real(sp),     dimension(:,:), allocatable :: powspecs
    ndet   = size(data%tod,2)
    nf    = size(data%tod,1)/2+1
    !allocate(ffts(nf,ndet), powspecs(nf, ndet))
    !call fft_multi(data%tod, ffts, 1)
    !do i = 1, ndet
    !   call extract_powspec(ffts(:,i), powspecs(:,i))
    !end do
  end subroutine;

end program scan_validate
