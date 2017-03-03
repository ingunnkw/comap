! CES validation tool: Run through each ces, and determine whether a diode-CES is acceptable for map making
program ces_validate
  use quiet_utils
  use quiet_fileutils
  use quiet_system_mod
  use quiet_fft_mod
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_filter_mod
  use quiet_mpi_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use quiet_shared_output_mod
  use quiet_acceptlist_mod
  use quiet_lx_mod
  use ces_validate_mod
  use quiet_patch_mod
  implicit none

  type info_struct
     type(shared_ofile) :: statfile, acceptfile
     integer(i4b)       :: id, nproc
     logical(lgt)       :: sigma0, alpha, f_knee, corr, validate
  end type info_struct

  character(len=512)   :: parfile, odir, lockfile, outfile, donefile, acceptlist_in, acceptlist_out
  character(len=512)   :: statfile, runlist_file, l3dir, summary_file
  integer(i4b)         :: ierr, cid, ndi, i, j, nces, cnum, di, nfreq, unit, mod, debug
  integer(i4b)         :: firstces, lastces
  logical(lgt)         :: accepted, apply_ces_diode_cuts, bright_object, use_templates
  logical(lgt)         :: use_filter_highpass, use_filter_lowpass, baseline ! TMR
  type(task_list)      :: tasks
  type(info_struct)    :: info
  type(acceptlist)     :: alist
  type(quiet_runlist)  :: runlist 
  type(quiet_ces_info) :: ces
  type(lx_struct)      :: data
  type(hdf_file)       :: file
  type(validation_struct), dimension(:,:),   allocatable :: stats
  integer(i4b), dimension(:,:),   allocatable :: my_status, status_tot
  real(sp),     dimension(:,:),   allocatable :: powspecs, tods, templates
  real(dp),     dimension(:),     allocatable :: F_filter, N_fft, scanfreqs, sftmp 
  real(dp),     dimension(:,:),   allocatable :: ftmp
  real(sp),     dimension(:),     allocatable :: stmp
  real(sp),     dimension(:,:),   allocatable :: ces_stats, ces_buffer, dtmp
  real(sp),     dimension(:,:,:), allocatable :: diode_stats, buffer
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
  call get_parameter(0, parfile, 'APPLY_CES_DIODE_CUTS',  par_lgt=apply_ces_diode_cuts)
  call get_parameter(0, parfile, 'DEBUG',                 par_int=debug)
  call get_parameter(0, parfile, 'FIRST_CES',             par_int=firstces)
  call get_parameter(0, parfile, 'LAST_CES',              par_int=lastces)
  call get_parameter(0, parfile, 'USE_TEMPLATES',         par_lgt=use_templates)
  call get_parameter(0, parfile, 'APPLY_HIGHPASS_FILTER', par_lgt=use_filter_highpass)
  call get_parameter(0, parfile, 'APPLY_LOWPASS_FILTER',  par_lgt=use_filter_lowpass)
  call get_parameter(0, parfile, 'MAKE_BASELINE_ACCLIST', par_lgt=baseline)

  if (baseline .and. (info%id == 0)) then
     write(*,*) 'Building baseline acceptlist. All but the initial obvious cuts disabled.'
     write(*,*) '***Note: summary.txt and ces_statistics.txt will not be valid.'
  end if

  acceptlist_out = trim(odir) // '/accept.txt'
  summary_file   = trim(odir) // '/summary.txt'
  lockfile       = trim(odir) // "/lock.dat"
  statfile       = trim(odir) // "/ces_statistics.txt"
  call dset(level=debug)

  call dmem("init")

  call mkdirs(trim(lockfile), .true.)
  call open_shared_ofile(info%statfile, trim(statfile), MPI_COMM_WORLD)
  call open_shared_ofile(info%acceptfile, trim(acceptlist_out), MPI_COMM_WORLD)

  call initialize_quiet_patch_mod(parfile);   call dmem("patch mod")
  call initialize_ces_mod(parfile);           call dmem("ces mod")
  call initialize_module_mod(parfile);        call dmem("module mod")
  call initialize_filter_mod(parfile);        call dmem("filter mod")
  call initialize_quiet_apex_mod(parfile);    call dmem("apex mod")

  ndi  = get_num_diodes() * get_num_modules()
  nces = get_num_ces()
  if(info%id == 0) call write_stat_headers(info)

  ! First scan through CES's: Build up diode-CES-statistics
  allocate(dtmp(ndi,NUM_DIODE_STATS),stmp(STAT_NUM),ftmp(ndi,5))
  allocate(diode_stats(ndi,DIODE_STAT_MAX,nces), buffer(ndi,DIODE_STAT_MAX,nces))
  allocate(ces_stats(STAT_MAX,nces), ces_buffer(STAT_MAX,nces))
  allocate(exists(nces),existbuf(nces))
  allocate(scanfreqs(nces),sftmp(nces))
  diode_stats = 0.d0
  ces_stats   = 0.d0
  exists      = .false.
  sftmp       = 0.d0
  if (apply_ces_diode_cuts) then
     call init_task_list(tasks, lockfile, nces, MPI_COMM_WORLD)
     do while(get_next_task(tasks, cnum))     
        call get_ces_info(cnum, ces)                   
        write(*,fmt="(i3,a,i4,a)") info%id, " scanning ces ", ces%cid, " (" // trim(itoa(cnum)) // "/" // &
             & trim(itoa(get_num_ces())) // ")"
        inquire(file=trim(ces%l3file), exist=exists(cnum))
        if (.not. exists(cnum)) cycle
        call open_hdf_file(ces%l3file, file, 'r')
        call read_hdf(file, 'stats',       stmp)
        call read_hdf(file, 'diode_stats', dtmp(:,:))
        call read_hdf(file, 'sigma0',      dtmp(:,DIODE_STAT_SIGMA0))
        call read_hdf(file, 'alpha',       dtmp(:,DIODE_STAT_ALPHA))
        call read_hdf(file, 'fknee',       dtmp(:,DIODE_STAT_FKNEE))
        call read_hdf(file, 'filter_par',  ftmp)
        dtmp(:,DIODE_STAT_NU_LOW)  = ftmp(:,FILTER_LOW_NU)
        dtmp(:,DIODE_STAT_NU_HIGH) = ftmp(:,FILTER_HIGH_NU_SCAN)
        dtmp(:,DIODE_STAT_AZORDER) = ftmp(:,FILTER_AZ_ORDER)
        call read_hdf(file, 'scanfreq',   sftmp(cnum))
        call close_hdf_file(file)
        ces_stats(1:STAT_MAX,cnum) = stmp
        diode_stats(:,1:DIODE_STAT_MAX,cnum) = dtmp
        call dmem("scanning " // trim(itoa(ces%cid)), 2)
     end do
     call mpi_allreduce(diode_stats, buffer,     size(buffer),     mpi_st, MPI_SUM, mpi_comm_world, ierr)
     call mpi_allreduce(ces_stats,   ces_buffer, size(ces_buffer), mpi_st, MPI_SUM, mpi_comm_world, ierr)
     call mpi_allreduce(exists,      existbuf,   size(exists),     mpi_logical,          MPI_LOR, mpi_comm_world, ierr)
     call mpi_allreduce(sftmp,       scanfreqs,  size(scanfreqs),  mpi_st, MPI_SUM, mpi_comm_world, ierr)
     diode_stats = buffer
     ces_stats   = ces_buffer
     exists      = existbuf
     call free_task_list(tasks)
  end if
  deallocate(buffer, ces_buffer, existbuf, stmp, dtmp, sftmp)
  call dmem("scan done")

  ! Initialize the validation module, including diode specific values
  call initialize_ces_validate_mod(parfile); call dmem("validate mod")

  ! Initialize stats structure
  allocate(stats(ndi,nces))
  call initialize_validation_structs(stats)

  ! ********************************************************************
  !                           Apply cuts 
  ! ********************************************************************
  call dmem("cuts start")

  ! Check for dead diodes
  do i = 1, ndi
     if (.not. is_alive(i)) then
        mod = (i-1)/get_num_diodes()
        di  = modulo(i-1,get_num_diodes())
        stats(i,:)%status = ior(stats(i,:)%status, status_baddi)
     end if
  end do
  call dmem("dead diodes")

  ! Missing files
  do i = 1, nces
     if(.not. exists(i)) stats(:,i)%status = ior(stats(:,i)%status, status_nofile)
  end do
  call dmem("missing files")

  ! Apply static cuts. Since the purpose of the static cuts is to cut
  ! extra ceses than the normal automatic cuts, it does not make sense
  ! that every ces-diode needs to be mentioned in the input accept list
  ! in order not to be rejected. The default for missing ceses in this
  ! list is therefore to accept them (and subject them to the normal
  ! set of automatic cuts).
  call initialize_accept_list(acceptlist_in, alist, REJECTED_NONE)
  do i = 1, nces
     call get_ces_info(i, ces)                   
     do j = 1, ndi
        mod = (j-1)/get_num_diodes()
        di  = modulo(j-1,get_num_diodes())
        if ((.not. is_accepted(alist,ces%cid,mod,di) .or. ces%cid < firstces .or. &
             & ces%cid > lastces) .and. stats(j,i)%status == status_ok) then
           stats(j,i)%status = ior(stats(j,i)%status, status_static)
        end if
     end do
  end do
  call deallocate_acceptlist(alist)
  call dmem("input cuts")

  ! Other simple static cuts:
  call validate_static(status_static, stats)
  call dmem("static cuts")

  ! Obvious non-varying cuts
  call validate_elementary(diode_stats, status_elementary, stats)
  call validate_gain(diode_stats, status_gain, stats)

  call dmem("elementary cuts")

  !call validate_jump(diode_stats(:,DIODE_STAT_JUMP,:), status_jump, stats) 
  ! Apply APEX cuts
  call validate_APEX(ces_stats, status_apex, stats)
  call dmem("simple cuts")

  ! For gc/gb: Skip all but the completely obvious cuts to get a baseline acceptlist
  if (baseline) goto 41

  ! Apply filter cuts
  call validate_filter_params(diode_stats(:,DIODE_STAT_NU_LOW,:), &
       & diode_stats(:,DIODE_STAT_NU_HIGH,:), scanfreqs, status_filter, stats)

  ! Make CES-diode cuts only if we have more than N CESs
  if (nces > 1 .and. apply_ces_diode_cuts) then
     call initialize_diode_ces_cuts(diode_stats)
     call dmem("init diode cuts")
     call validate_type_B_chisq(diode_stats(:,DIODE_STAT_TYPEB,:), status_typeb, stats)
     call dmem("type-B cut")
     call validate_weather(diode_stats(:,DIODE_STAT_WEATHER1,:), DIODE_STAT_WEATHER1, status_weather, stats)  ! 10 seconds
     call validate_weather(diode_stats(:,DIODE_STAT_WEATHER2,:), DIODE_STAT_WEATHER2, status_weather, stats)  ! 30 seconds
     call dmem("weather cut")
     !call validate_jump(diode_stats(:,DIODE_STAT_JUMP,:),        status_jump,   stats)
     !call dmem("jump cut")
     call validate_sigma0(diode_stats(:,DIODE_STAT_SIGMA0,:),    status_sigma0, stats)
     call dmem("sigma0 cut")
     call validate_fknee(diode_stats(:,DIODE_STAT_FKNEE,:),      status_fknee,  stats)
     call dmem("fknee cut")
     call validate_fbias(diode_stats,                            status_fbias,  stats)
     call dmem("fbias cut")
     call validate_sun(diode_stats(:,DIODE_STAT_SIDELOBE_HIT,:), status_sun,    stats)
     call dmem("sun cut")
  end if

  ! Process all CES's
41  call init_task_list(tasks, lockfile, nces, MPI_COMM_WORLD)
  do while(get_next_task(tasks, cnum))     
     if(.not. exists(cnum)) cycle
     call dmem('ces begin')
     call get_ces_info(cnum, ces)

     ! Do not cut on demod for strong objects
     bright_object = trim(ces%object) == 'moon'

     write(*,fmt="(i3,a,i4,a)") info%id, " processing ces ", ces%cid, " (" // trim(itoa(cnum)) // "/" // &
          & trim(itoa(get_num_ces())) // ")"

     call read_l3_file(ces%l3file, data)            ; call dmem('read_l3_file')
     call calc_fourier(data, ffts, powspecs)        ; call dmem("fourier")
     nfreq = size(powspecs,1)
     ! TMR
     if (use_templates) then
        call get_templates(templates, nfreq, data%samprate)
        powspecs = powspecs/templates ! powspecs(nf,ndi)
     end if

     allocate(tods(size(data%tod,1),size(data%tod,2)))
     allocate(F_filter(nfreq),N_fft(nfreq))

     if (baseline) goto 42     

     ! Is it faster/better to do this one diode at a time?
     call compute_filtered_tod(data, ffts, tods, use_filter_lowpass, use_filter_highpass)    ; call dmem("filter_tod")
     ! Check statistics for each diode
     do i = 1, ndi
        mod = (i-1)/get_num_diodes()
        di  = modulo(i-1,get_num_diodes())

        ! Initialize statistic structures
        call update_current_stat(mod, di, stats(i,cnum), data, diode_stats(:,:,cnum))
        if (.not. is_active(stats(i,cnum)%status)) cycle

        ! Apply FFT cuts
        call get_N_filter(data%samprate, data%sigma0(i), data%fknee(i), &
         & data%alpha(i), 1.d0, N_fft, use_filter_highpass)

        ! Now using individually tuned filter parameters
        F_filter = 1.d0
        if(use_filter_lowpass)  call apodize_filter_fft(.true., size(tods,1), data%samprate, &
             & data%filter_par(i,FILTER_LOW_NU), data%filter_par(i,FILTER_LOW_ALPHA), .false., F_filter)
        if(use_filter_highpass) call apodize_filter_fft(.true.,size(tods,1), data%samprate, &
             & data%filter_par(i,FILTER_HIGH_NU_SCAN)*data%scanfreq, data%filter_par(i,FILTER_HIGH_ALPHA), .true., F_filter)

        call validate_fft(is_polarization_module(mod), data%samprate, &
         & data%scanfreq, F_filter, N_fft, powspecs(:,i), stats(i,cnum))

        ! Apply TOD cuts
        if (use_templates) then
           call validate_tod(is_polarization_module(mod), (sum(F_filter**2*N_fft*templates(:,i))) / &
                & (nfreq), data%orig_point(1,:), tods(:,i), stats(i,cnum))
        else
           call validate_tod(is_polarization_module(mod), (sum(F_filter**2*N_fft)) / &
                & (nfreq), data%orig_point(1,:), tods(:,i), stats(i,cnum))
        end if
     end do

     ! FIXME: Validate_corr must be generalized
     !call validate_corr(data%corr, status_corr, stats(:,cnum)); call dmem("validate corr")
     ! If we have a bright object, ignore certain cuts
     if(bright_object) stats(:,cnum)%status = iand(stats(:,cnum)%status, not(mask_bright))
     call validate_global_ces(data, stats(:,cnum));           ; call dmem("global")
42   call output_statistic_file(ces%cid, stats(:,cnum), info) ; call dmem('output_stat_file')
     call output_accept_file(ces%cid, stats(:,cnum), info)    ; call dmem('output_accept_file')

     call free_lx_struct(data)
     deallocate(ffts, powspecs, tods, F_filter, N_fft)
     if(allocated(templates)) deallocate(templates)
  end do
  
  call close_shared_ofile(info%statfile) 
  call close_shared_ofile(info%acceptfile) 

  ! Generate summary statistics
  allocate(my_status(ndi,nces), status_tot(ndi,nces))
  do i = 1, nces
     do j = 1, ndi
        my_status(j,i) = stats(j,i)%status
     end do
  end do

  call mpi_reduce(my_status, status_tot, size(my_status), mpi_integer, MPI_MAX, 0, mpi_comm_world, ierr)

  if (info%id == 0) then
     call output_season_summary(status_tot, summary_file)
  end if
  deallocate(my_status, status_tot)

  deallocate(stats, diode_stats, scanfreqs)
  call free_task_list(tasks)
  call dmem("finished")
  call mpi_finalize(ierr)

  if (info%id == 0) then
     write(*,*)
     write(*,*) 'ces_validate completed successfully'
     write(*,*)
  end if


contains

  ! Routines for statistics output
  subroutine output_statistic_file(cid, stats, info)
    implicit none
    type(info_struct),                     intent(in) :: info
    integer(i4b),                          intent(in) :: cid
    type(validation_struct), dimension(:), intent(in) :: stats

    integer(i4b)        :: i, j, accept, ndi, cnum, mod, di
    character(len=10000) :: line
    character(len=10)    :: object
    type(quiet_ces_info) :: ces

    ndi        = size(stats)
    call get_ces_info(lookup_ces(cid), ces)
    object      = ces%object

    do i = 1, ndi
       mod = (i-1)/get_num_diodes()
       di  = modulo(i-1,get_num_diodes())
       write(line,fmt='(i8,a13,i4,i6,z9,i10,f7.2,e12.4,4f10.3,i6,15f12.2)') &
            & cid, adjustr(object), mod, di, stats(i)%status, stats(i)%status, &
            & stats(i)%ces_accept_ratio, stats(i)%sigma0, stats(i)%alpha, &
            & min(stats(i)%fknee,100.d0), stats(i)%nu_high, &
            & stats(i)%nu_low, nint(stats(i)%azorder), stats(i)%full_fft_chisq, &
            & stats(i)%scanfreq_fft_chisq, stats(i)%band_fft_chisq(1:6), &
            & stats(i)%tod_chisq, stats(i)%tod_absmax, stats(i)%tod_chi_az, &
            & stats(i)%typeb, 1d4*stats(i)%wthr1, 1d4*stats(i)%wthr2, stats(i)%jump
       call write_shared_ofile(info%statfile, trim(line))
    end do
    call flush_shared_ofile(info%statfile)
  end subroutine output_statistic_file

  subroutine output_accept_file(cid, stats, info)
    implicit none
    type(info_struct),                     intent(in) :: info
    integer(i4b),                          intent(in) :: cid
    type(validation_struct), dimension(:), intent(in) :: stats

    integer(i4b)        :: i, j, mod, di
    character(len=10000) :: line

    ! Write line to file
    do i = 1, size(stats)
       mod = (i-1)/get_num_diodes()
       di  = modulo(i-1,get_num_diodes())
       if (stats(i)%status == status_ok) then
          write(line,fmt='(i7,i5,i2,i3)') cid, mod, di, 1
       else 
          write(line,fmt='(i7,i5,i2,i3)') cid, mod, di, 0
       end if
       call write_shared_ofile(info%acceptfile, trim(line))
    end do
    call flush_shared_ofile(info%acceptfile)
  end subroutine output_accept_file


  subroutine write_stat_headers(info)
    implicit none
    type(info_struct) :: info
    integer(i4b) :: unit
    call write_shared_ofile(info%statfile, '#    CES_id    Object  Mod   Di  StatHex StatDec' // &
         & '   R      Sigma0      Alpha     Fknee   nu_high  nu_low   az_order   0_12.5Hz' //&
         & '       Scan        0_0.2Hz       CMB       10_12.5Hz     1.2Hz    FFT_spike' // &
         & '1.0Hz    TOD_chisq  TOD_abs    chi_az   type_B    Wthr_10s   Wthr_30s   Jump')
    call flush_shared_ofile(info%statfile)
  end subroutine write_stat_headers

  subroutine compute_filtered_tod(data, ffts, tods, lowpass, highpass)
    implicit none

    type(lx_struct)                           :: data
    complex(spc), dimension(:,:), intent(in)  :: ffts
    real(spc),    dimension(:,:), intent(out) :: tods
    logical(lgt),                 intent(in)  :: lowpass, highpass
    
    integer(i4b) :: i, ndi, nsamp, nfreq
    real(dp),     allocatable, dimension(:)   :: F, sqrt_invN
    complex(spc), allocatable, dimension(:,:) :: ffts_int

    nsamp   = size(tods,1)
    ndi     = size(tods,2)
    nfreq   = size(ffts,1)

    allocate(F(nfreq), sqrt_invN(nfreq), ffts_int(nfreq,ndi))

    do i = 1, ndi
       F = 1.d0
       if(lowpass)  call apodize_filter_fft(.true., nsamp, data%samprate, &
            & data%filter_par(i,FILTER_LOW_NU), data%filter_par(i,FILTER_LOW_ALPHA), .false., F)
       if(highpass) call apodize_filter_fft(.true., nsamp, data%samprate, &
            & data%filter_par(i,FILTER_HIGH_NU_SCAN)*data%scanfreq, data%filter_par(i,FILTER_HIGH_ALPHA), .true., F)

       ! TMR: Why do we do this?
       sqrt_invN     = 1.d0
       ffts_int(:,i) = ffts(:,i) * F * sqrt_invN
    end do
    call fft_multi(tods, ffts_int, -1)

    deallocate(F, sqrt_invN, ffts_int)

  end subroutine compute_filtered_tod

  subroutine output_season_summary(stats, summary_file)
    implicit none

    integer(i4b),            dimension(:,:), intent(in) :: stats
    character(len=*),                        intent(in) :: summary_file

    integer(i4b) :: i, j, k, l, m, unit, nces, typ, nper
    integer(i8b) :: fields
    logical(lgt) :: pol, mask(size(stats,1))
    integer(i4b) :: cut_count(num_status), cum_count(num_status), numtot, numtot_pruned
    type(filter_params) :: filteropts

    nces      = size(stats,2)
    call get_default_filter(filteropts)

    unit = getlun()
    open(unit, file=trim(summary_file))
    write(unit,*) '******** CES VALIDATION SUMMARY ***********'
    write(unit,*)

    do typ = 1, 2
       pol = typ == 1

       write(unit,*)
       if (pol) then
          write(unit,*) '----- Polarization -----'
          mask = abs(quiet_diodes%stokes(2)) > abs(quiet_diodes%stokes(1))
       else
          write(unit,*) '----- Temperature -----'
          mask = abs(quiet_diodes%stokes(2)) <= abs(quiet_diodes%stokes(1))
       end if
       nper = count(mask)

       do i = 1, size(patches)
          ! What ces-diodes in this patch are cut by which cuts?
          cut_count = 0
          cum_count = 0
          numtot    = 0

          do j = 1, nces
             call get_ces_info(j,ces)
             if(patches(i)%name /= ces%object) cycle
             numtot = numtot  + nper
             fields = 0
             do k = 1, num_status
                fields = ior(fields,ishft(1,k-1))
                cut_count(k) = cut_count(k) + count(mask .and. btest(stats(:,j),k-1))
                cum_count(k) = cum_count(k) + count(mask .and. iand (stats(:,j),fields)==0)
             end do
          end do

          ! Calculating accept rates relative to proper data only - leaving out 3 first cuts (bad diodes,
          ! missing files and static cuts (various garbage)) from the accounting 
          numtot_pruned = cum_count(3)

          ! Unexposed patch
          if(numtot == 0) cycle

          ! And output
          write(unit,*)
          write(unit,*) trim(patches(i)%name)

          1 format('    Total number of CES-diodes ', a20,' = ',i8,', ',f7.2,'%')
          2 format('    Number of CES-diodes after ', a20,' = ',i8,', ',f7.2,'%',', ',f7.2,'%')
          3 format('    Total number of viable CES-diodes ', a20,' = ',i8,', ',f7.2,'%')
          write(unit,1) '', numtot, 100.0
          do j = 1, 2
             write(unit,2), status_desc(j), cum_count(j), 100.0*cum_count(j)/numtot, 100.0*cut_count(j)/numtot
          end do

          write(unit,3) '', numtot_pruned, 100.0
          do j = 3, num_status
             write(unit,2), status_desc(j), cum_count(j), 100.0*cum_count(j)/numtot_pruned, 100.0*cut_count(j)/numtot_pruned
          end do
       end do
    end do

    write(unit,*)
    write(unit,*) ' ----- Hex key -----'
    do i = 1, num_status
       write(unit,'(z6,a,a)') ishft(1,i-1), " for ", trim(status_desc(i))
    end do

    ! Output summary of cut parameters
    write(unit,*)
    write(unit,*)
    write(unit,*) ' ----- Filter parameters -----'
    write(unit,*)
    write(unit,*) '  APPLY_LOWPASS_FILTER               = ', filteropts%apply_lowpass
    write(unit,*) '  APPLY_HIGHPASS_FILTER              = ', filteropts%apply_highpass
    write(unit,*) '  NU_HIGHPASS (scanfreq)             = ', real(filteropts%nu_highpass,sp)
    write(unit,*) '  ALPHA_HIGHPASS                     = ', real(filteropts%alpha_highpass,sp)
    write(unit,*) '  NU_LOWPASS (Hz)                    = ', real(filteropts%nu_lowpass,sp)
    write(unit,*) '  ALPHA_LOWPASS                      = ', real(filteropts%alpha_lowpass,sp)
    write(unit,*)
    write(unit,*)
    write(unit,*) ' ----- Cut definitions -----'
    write(unit,*) 
    write(unit,*) ' Level 1 -- dead diodes'
    write(unit,*) '    Taken from focalplane definition file'
    write(unit,*) 
    write(unit,*) ' Level 2 -- missing file cuts'
    write(unit,*) 
    write(unit,*) ' Level 3 -- static cuts'
    write(unit,*) '    Taken from input accept list'
    write(unit,*) 
    write(unit,*) ' Level 4 -- elementary cuts'
    write(unit,*) '    Sanity checks'
    write(unit,*) 
    write(unit,*) ' Level 5 -- correlation cuts'
    write(unit,*) '    Are the correlations invertible?'
    write(unit,*) 
    write(unit,*) ' Level 6 -- APEX cuts'
    write(unit,*) '    APEX_MAX_PWV (mm PWV)             = ', real(apex_max_pwv,sp)
    write(unit,*) 
    write(unit,*) ' Level 7 -- Type-B cuts'
    write(unit,*) '    TYPE_B_CHISQ_THRESHOLD (sigma)    = ', real(typeb_threshold,sp)
    write(unit,*) 
    write(unit,*) ' Level 8 -- Weather cuts'
    write(unit,*) '    WEATHER_CHISQ_THRESHOLD (sigma)   = ', real(weather_threshold,sp)
    write(unit,*) '    WEATHER_CES_FRACTION (fraction)   = ', real(weather_ces_fraction,sp)
    write(unit,*) 
    write(unit,*) ' Level 9 -- sigma0 cuts'
    write(unit,*) '    SIGMA0_THRESHOLD (sigma)          = ', real(sigma0_threshold,sp)
    write(unit,*) 
    write(unit,*) ' Level 10 -- jump cuts'
    write(unit,*) '    JUMP_FINDER_FRACTIONAL_JUMP_THRESHO= ', real(jump_finder_fractional_jump_threshold,sp)
    write(unit,*) 
    write(unit,*) ' Level 11 -- tod cuts'
    write(unit,*) '    TOD_CHISQ_THRESHOLD (sigma)       = ', real(tod_chisq_threshold,sp)
    write(unit,*) '    TOD_ABSMAX_THRESHOLD (sigma)      = ', real(tod_absmax_threshold,sp)
    write(unit,*) '    TOD_AZ_BINSIZE (degrees)          = ', real(tod_az_binsize,sp)
    write(unit,*) '    TOD_AZ_MAX_CHISQ (sigma)          = ', real(tod_az_max_chisq,sp)
    write(unit,*)
    write(unit,*) ' Level 12 -- fft cuts'
    write(unit,*) '    FFT_CHISQ_THRESHOLD (CMB; sigma)  = ', real(fft_chisq_threshold,sp)
    write(unit,*) '    FFT_CHISQ_1_2HZ_THRESHOLD (sigma) = ', real(fft_chisq_1_2Hz_threshold,sp)
    write(unit,*) '    FFT_CHISQ_SPIKE_THRESHOLD (sigma) = ', real(fft_spike_threshold,sp)
    write(unit,*) '    FFT_CHISQ_SCAN_THRESHOLD (sigma)  = ', real(fft_chisq_scan_threshold,sp)
    write(unit,*) '    FFT_CHISQ_LOW_THRESHOLD (sigma)   = ', real(fft_chisq_low_threshold,sp)
    write(unit,*) '    FFT_CHISQ_HIGH_THRESHOLD (sigma)  = ', real(fft_chisq_high_threshold,sp)
    write(unit,*) ' Level 13 -- sun cuts'
    write(unit,*) '    Based on sidelobe map files'
    write(unit,*)
    write(unit,*) ' Level 14 -- accept ratio cut'
    write(unit,*) '    Rejects if the diode fraction rejected'
    write(unit,*) '    due to a cause which should logically'
    write(unit,*) '    affect all diodes, exceeds:'
    write(unit,*) '    CES_MIN_DIODE_ACCEPT_RATIO        = ', real(min_diode_accept_ratio,sp)
    write(unit,*) 
    close(unit)
  end subroutine output_season_summary

  subroutine calc_fourier(data, ffts, powspecs)
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: ndi, i, nf
    complex(spc), dimension(:,:), allocatable :: ffts
    real(sp),     dimension(:,:), allocatable :: powspecs
    ndi   = size(data%tod,2)
    nf    = size(data%tod,1)/2+1
    allocate(ffts(nf,ndi), powspecs(nf, ndi))
    call fft_multi(data%tod, ffts, 1)
    do i = 1, ndi
       call extract_powspec(ffts(:,i), powspecs(:,i))
    end do
  end subroutine

end program ces_validate
