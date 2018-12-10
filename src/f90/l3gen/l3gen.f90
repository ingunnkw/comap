! l3gen: Produce proper L3 files from L2 files
program l3gen
  use quiet_utils
  use comap_defs
  use quiet_fileutils
  use quiet_system_mod
  use comap_scan_mod
  use comap_detector_mod
  use quiet_mpi_mod
  use quiet_mpi_utils
  use quiet_task_mod
  use comap_lx_mod
  !use quiet_fft_mod
  use comap_pointing_mod
  use comap_noise_estimation_mod
  use powell_mod
  use comap_gain_mod
  !use quiet_stat_mod
  !use quiet_patch_detect_mod
  !use quiet_sidelobe_mod
  use quiet_nr_mod
  !use quiet_filter_mod
  use quiet_hdf_mod
  use spline_1d_mod
  use comap_patch_mod
  implicit none

  type info_struct
     integer(i4b)       :: id, nproc
  end type info_struct

  character(len=512)    :: parfile, odir, outfile, tsys_loc, freqmaskfile
  character(len=512)    :: lockfile, tmpfile, point_objs, tilt_file, offset_mask_file, coord_out
  integer(i4b)          :: ierr, snum, nmod, i, j, isys, osys, mod, nside_l3, debug
  integer(i4b)          :: num_corr_bins, nscan, numfreq, unit
  real(dp)              :: fix_highpass_freq_scan, fix_lowpass_freq, t1, t2, scanmask_width
  real(dp)              :: scanfreq_min, scanfreq_max
  logical(lgt)          :: reprocess, exist, scanmask, inter_module, no_filters, use_templates, found
  type(task_list)       :: tasks
  type(info_struct)     :: info
  type(comap_scan_info) :: scan
  type(lx_struct)       :: data, data_l2_fullres
  type(patch_info)      :: pinfo
  real(sp),     dimension(:,:,:,:), allocatable :: powspecs
  complex(spc), dimension(:,:,:,:), allocatable :: ffts
  logical(lgt), dimension(:,:),     allocatable :: mask

  call getarg(1, parfile)
  call get_parameter(0, parfile, 'OUTPUT_DIR',           par_string=odir)
  call get_parameter(0, parfile, 'TSYS_LOC',             par_string=tsys_loc)
  call get_parameter(unit, parfile, 'FREQUENCY_MASK',    par_string=freqmaskfile)
  call get_parameter(unit, parfile, 'NUMFREQ',           par_int=numfreq)

  call get_parameter(0, parfile, 'REPROCESS_ALL_FILES',  par_lgt=reprocess)
  call get_parameter(0, parfile, 'DEBUG',                par_int=debug)
  call get_parameter(0, parfile, 'L3_FAST',              par_lgt=no_filters, desc=&
       & "Set this to .false.")
  call get_parameter(0, parfile, 'FIX_HIGHPASS_FREQ',    par_dp=fix_highpass_freq_scan, desc=&
       & "Set this to negative to fit")
  call get_parameter(0, parfile, 'FIX_LOWPASS_FREQ',     par_dp=fix_lowpass_freq, desc=&
       & "Set this to negative to fit")
  call get_parameter(0, parfile, 'L2_COORDINATE_SYSTEM', par_string=coord_out, desc=&
       & "Coordinate system for L2 and L3 files; galactic or celestial")
  call get_parameter(0, parfile, 'APPLY_SCANMASK',      par_lgt=scanmask, desc=&
   & "Whether to mask out freqs near scanfreq when estimating noise parameters. Should normally be .true.")
  call get_parameter(0, parfile, 'SCANMASK_WIDTH',      par_dp=scanmask_width, desc=&
   & "Window to cut around scan frequency and its harmonics during noise estimation, in Hz.")
  call get_parameter(0, parfile, 'SCANFREQ_MIN',      par_dp=scanfreq_min, desc=&
   & "Smallest allowed scanning frequency in Hz.")
  call get_parameter(0, parfile, 'SCANFREQ_MAX',      par_dp=scanfreq_max, desc=&
   & "Largest allowed scanning frequency in Hz.")

  isys = COORD_CEL
  if (trim(coord_out) == 'galactic') then
     osys = COORD_GAL
  else
     osys = COORD_CEL
  end if

  call wall_time(t1)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, info%id, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, info%nproc, ierr)
  call dset(id=info%id,level=debug)
  !call print_host_mapping

  lockfile = trim(odir) // "/lock.dat"
  call mkdirs(trim(lockfile), .true.)

  call initialize_scan_mod(parfile);             call dmem("scan mod")
  call initialize_noise_estimation_mod(parfile); call dmem("noise mod")
  call initialize_comap_pointing_mod(parfile);   call dmem("pointing mod")
  call initialize_gain_mod(parfile);             call dmem("gain mod")
  call initialize_comap_patch_mod(parfile)
  !call initialize_patch_detect_mod(parfile);     call dmem("patch detect mod")
  !call initialize_filter_mod(parfile);           call dmem("filter mod")

  ! Process all CES's
  !call init_task_list(tasks, lockfile, get_num_scans(), MPI_COMM_WORLD)
  nscan = get_num_scans()
  do snum = 1+info%id, nscan, info%nproc
 !do while(get_next_task(tasks, snum))
     call get_scan_info(snum, scan)
     found = get_patch_info(scan%object, pinfo) 
     tmpfile = trim(scan%l3file) // ".part"
     inquire(file=tmpfile,exist=exist)
     if(exist) then
        write(*,fmt="(i3,a,2i5)") info%id, " found incomplete run:", snum, scan%sid
        call rm(tmpfile)
     end if
     inquire(file=scan%l3file,exist=exist)
     if(exist .and. .not. reprocess) then
        write(*,fmt="(i3,a,2i5,a)") info%id, " skipping already finished scan:", snum, scan%sid
        cycle
     end if
     inquire(file=scan%l2file,exist=exist)
     if(.not. exist) then
        write(*,fmt="(i3,a,2i5,a)") info%id, " data missing for scan:", snum, scan%sid, " " // trim(scan%l2file)
        cycle
     end if
     write(*,fmt="(i3,a,i4,a)") info%id, " processing scan ", scan%sid, " (" // trim(itoa(snum)) // "/" // &
          & trim(itoa(get_num_scans())) // ")"

     call dmem("scan start")

     ! Read L2 data
     call read_L2_file(scan%l2file, data)     ; call dmem("read l2")

     ! Initialize frequency mask
     !call initialize_frequency_mask(freqmaskfile, numfreq, data)

     ! Process it
     !call calc_weather_template(data)         ; call dmem("weather")
     call calc_point (data,  isys, osys)      ; call dmem("point")
     !call calc_objrel(data,  point_objs)      ; call dmem("objrel")
     !call calc_pixels(data, nside_l3)         ; call dmem("pixels")
     if (trim(pinfo%type) == 'gal' .or. trim(pinfo%type) == 'cosmo') then
        call apply_gain_cal(tsys_loc,data)
        call calc_scanfreq(data);                ; call dmem("scanfreq")
        !call apply_az_filter(data)               ; call dmem("az_filter")
        call calc_fourier(data, ffts, powspecs)  ; call dmem("fourier")
        call fit_noise(data, powspecs, snum)     ; call dmem("noise")
        !call calc_gain(data)                     ; call dmem("gain")
        if (allocated(ffts)) deallocate(ffts)
        if (allocated(powspecs)) deallocate(powspecs)
     end if
     !call calc_diode_stats(data, powspecs)    ; call dmem("diode_stats")
     !call calc_stats(data)                    ; call dmem("stats")
     !allocate(data%filter_par(size(data%tod,2),NUM_FILTER_PAR))
     !data%filter_par = -1
     !call calc_az_filter_par(data)            ; call dmem("az_par")
     !call calc_bandpass_filter_par(data, powspecs)     ; call dmem("bandpass")

     ! Output L3 data
     call mkdirs(tmpfile, .true.)
     call write_L3_file(tmpfile, data)        ; call dmem("write l3")

     call mkdirs(trim(scan%l3file), .true.)
     call mv(tmpfile, scan%l3file)

     call free_lx_struct(data)
  end do

  call free_task_list(tasks)
  write(*,fmt="(i3,a)") info%id, " finished"
  
  call mpi_finalize(ierr)
  call wall_time(t2)

  write(*,*) 'Time elapsed: ', t2-t1
contains

  subroutine initialize_frequency_mask(freqmaskfile, nfreq, data)
    implicit none
    character(len=*),                                intent(in)    :: freqmaskfile
    integer(i4b),                                    intent(in)    :: nfreq
    type(Lx_struct),                                 intent(inout) :: data


    integer(i4b) :: i, j, k, nfreq_full, nsb, ndet, unit, det, sb, freq, dfreq, ierr
    logical(lgt) :: first
    character(len=1024) :: line, val, equal

    ndet = get_num_dets()
    nsb  = get_num_sideband()
    unit = getlun()
    open(unit, file=trim(freqmaskfile), recl=1024)
    first = .true.
    do while (.true.)
       read(unit,'(a)', end=99) line
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle
       if (first) then
          read(line,*) nfreq_full
          allocate(data%freqmask(nfreq,nsb,ndet))
          allocate(data%freqmask_full(nfreq_full,nsb,ndet))
          data%freqmask_full = 1.d0
          first = .false.
       else
          read(line,*) freq, sb, det
          if (det == 0 .and. sb == 0 .and. freq == 0) then
             write(*,*) 'ERROR: All frequencies removed by freqmask!'
             stop
          else if (det == 0 .and. sb == 0) then
             data%freqmask_full(freq,:,:) = 0.d0
          else if (det == 0 .and. freq == 0) then
             data%freqmask_full(:,sb,:) = 0.d0
          else if (sb == 0 .and. freq == 0) then
             data%freqmask_full(:,:,det) = 0.d0
          else if (det == 0) then
             data%freqmask_full(freq,sb,:) = 0.d0
          else if (sb == 0) then
             data%freqmask_full(freq,:,det) = 0.d0
          else if (freq == 0) then
             data%freqmask_full(:,sb,det) = 0.d0
          else
             data%freqmask_full(freq,sb,det) = 0.d0
          end if
          if (all(data%freqmask_full == 0)) then
             write(*,*) 'ERROR: All frequencies removed by freqmask!'
             stop
             stop
          end if
       end if
    end do
99  close(unit)

    do k = 1, ndet
       if (.not. is_alive(k)) data%freqmask_full(:,:,k) = 0.d0
    end do

    dfreq = nfreq_full/nfreq
    if (dfreq == 1) then
       data%freqmask = data%freqmask_full
    else
       data%freqmask = 1.
       do k = 1, ndet
          do j = 1, nsb
             do i = 1, nfreq
                if (all(data%freqmask_full((i-1)*dfreq+1:i*dfreq,j,k) == 0.d0)) data%freqmask(i,j,k) = 0.0
             end do
          end do
       end do
    end if

  end subroutine initialize_frequency_mask

  subroutine apply_gain_cal(tsys_file, data)
    implicit none
    character(len=*),            intent(in)       :: tsys_file
    type(Lx_struct),             intent(inout)    :: data
    type(hdf_file)                                :: file
    real(dp), dimension(:,:,:,:), allocatable     :: tsys_fullres
    integer(i4b)                                  :: nfreq, nfreq_fullres, nsb, ndet, i, j, k, l, mjd_index, nsamp, dnu
    integer(i4b)                                  :: nsamp_gain(7), num_bin, n
    real(dp)                                      :: mjd_high,w, sum_w_t, sum_w, t1, t2, tsys
    real(dp), dimension(:), allocatable           :: time

    nsamp = size(data%tod,1)
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)
    allocate(data%gain(1,nfreq,nsb,ndet), data%time_gain(1))
    data%time_gain(1) = data%mjd_start
    if (nsamp == 0) return

    ! 1) get tsys-values
    call open_hdf_file(tsys_file, file, "r")
    call read_hdf(file, "nfreq", nfreq_fullres)
    call get_size_hdf(file, "MJD", nsamp_gain)

    !allocate(tsys_fullres(nsamp_gain(1), nfreq_fullres, nsb, ndet))
    allocate(tsys_fullres(ndet, nsb, nfreq_fullres, nsamp_gain(1)))
    allocate(time(nsamp_gain(1)))

    call read_hdf(file, "MJD", time)
    call read_hdf(file, "tsys", tsys_fullres)
    call close_hdf_file(file)
    ! finding closest time-value
    mjd_index = max(locate(time, data%mjd_start),1)
    write(*,*) 'time = ', time(mjd_index), data%mjd_start, mjd_index
    !stop
!!$    mjd_start = data%mjd_start
!!$    mjd_high = 1.d10
!!$    mjd_index = 10000000
!!$    do i = 1, nsamp(1)
!!$       if (abs(mjd_start-time(i)) < mjd_high) then
!!$          mjd_high = abs(mjd_start-time(i))
!!$          mjd_index = i
!!$       end if
!!$    end do

    dnu   = nfreq_fullres/nfreq
    do i=1, ndet
       if (.not. is_alive(i)) cycle
       do j=1, nsb
          do k=1, nfreq
             if (sum(data%freqmask_full((k-1)*dnu+1:k*dnu,j,i)) == 0) then
                ! No live frequencies; return zero
                tsys = 0.d0
                data%tod(:,k,j,i) = 0
             else
                do l=(k-1)*dnu+1,k*dnu
                   ! Remove NaN's
                   if (isnan(tsys_fullres(i, j,l,mjd_index))) tsys_fullres(i, j, l,mjd_index) = 0.0
                end do
                tsys = sum(tsys_fullres(i, j, (k-1)*dnu+1:k*dnu,mjd_index) * &
                     & 1.d0/data%var_fullres((k-1)*dnu+1:k*dnu,j,i)*data%freqmask_full((k-1)*dnu+1:k*dnu,j,i)) / &
                     & sum(1.d0/data%var_fullres((k-1)*dnu+1:k*dnu,j,i)*data%freqmask_full((k-1)*dnu+1:k*dnu,j,i))
                data%gain(1,k,j,i) = tsys
                data%tod(:,k,j,i) = data%tod(:,k,j,i)*tsys
             end if
          end do
       end do
    end do
  deallocate(tsys_fullres)
  deallocate(time)

  end subroutine apply_gain_cal


  subroutine apply_az_filter(data)
    implicit none
    type(lx_struct) :: data

    integer(i4b) :: i, j, k, l, det, freq, sb, n, i_high, i_low, n_scan, n_tod
    real(dp)     :: a, b
    real(dp), allocatable, dimension(:) :: slopes
    
    n_tod  = size(data%tod,1)
    n_scan = 0.5 * data%samprate / data%scanfreq(1) 

    open(58,file='tod.dat')
    do i = 1, n_tod-1
       write(58,*) i, data%tod(i,30,4,3), data%point_tel(1,i+1,1)-data%point_tel(1,i,1)
    end do
    close(58)
    call mpi_finalize(ierr)
    stop

    i_low = 1
    open(58,file='az1.dat')
    open(59,file='az2.dat')
    do while (i_low <= n_tod)
       i_high = i_low+1
       do while (i_high <= n_tod-1 .and. &
            & ((data%point_tel(1,i_high-1,1)-data%point_tel(1,i_high,1))*&
            & (data%point_tel(1,i_high,1)-data%point_tel(1,i_high+1,1)) > 0 .or. i_high-i_low < n_scan))
          i_high = i_high + 1
       end do
       n = i_high - i_low + 1

       allocate(slopes(n*(n+1)/2))
       
       do j = i_low, i_high
          write(58,*) j, data%tod(j,30,4,3)-1
       end do

!       do det = 1, size(data%tod,4)
!          do sb = 1, size(data%tod,3)
       do det = 3, 3
          do sb = 4, 4
             write(*,*) det, sb
!             do freq = 1, size(data%tod,2)
             do freq = 30, 30

                ! Subtract median-evaluated linear azimuth function
                l = 0
                do j = i_low, i_high-1
                   do k = j+1, i_high
                      l = l+1
                      if (abs(data%point_tel(1,j,1)-data%point_tel(1,k,1)) / abs(data%point_tel(1,j,1)) < 1.d-6) then
                         slopes(l) = 1.d30
                      else
                         slopes(l) = (data%tod(j,freq,sb,det)-data%tod(k,freq,sb,det)) / &
                              & (data%point_tel(1,j,1) - data%point_tel(1,k,1))
                      end if
                   end do
                end do
                a = median(slopes)
                b = median(data%tod(i_low:i_high,freq,sb,det)-a*data%point_tel(1,i_low:i_high,1))
                data%tod(i_low:i_high,freq,sb,det) = data%tod(i_low:i_high,freq,sb,det) - &
                     & (a*data%point_tel(1,i_low:i_high,1)+b)

             end do
          end do
       end do

       do j = i_low, i_high
          write(59,*) j, data%tod(j,30,4,3)
       end do
       i_low = i_high+1

       deallocate(slopes)
    end do
    close(58)
    close(59)
    call mpi_finalize(ierr)
    stop
    

  end subroutine apply_az_filter

  subroutine calc_weather_template(data)
    implicit none

    type(lx_struct) :: data


    integer(i4b) :: i, j, k, l, ndet, nsamp, nfreq, nsb, n
    real(dp)     :: sigma0, alpha, fknee, nu, chisq, g
    real(dp),     allocatable, dimension(:) :: dt 
    complex(dpc), allocatable, dimension(:) :: dv

    nsamp  = size(data%tod,1)
    nfreq  = size(data%tod,2)
    nsb    = size(data%tod,3)
    ndet   = size(data%tod,4)
    n      = nsamp+1

    allocate(data%weather_temp(nsamp,nsb,ndet), data%rel_gain(nfreq,nsb,ndet))
    allocate(dt(2*nsamp), dv(0:n-1))
    do i = 2, ndet
       do j = 1, nsb

          ! Normalize gain with a low-pass filter
          data%weather_temp(:,j,i) = 0.d0
          open(58,file='tods.dat')
          do k = 1, nfreq
             dt(1:nsamp) = data%tod(:,k,j,i)
             dt(2*nsamp:nsamp+1:-1) = dt(1:nsamp)
             call fft(dt, dv, 1)
             dv(21:n-1) = 0.d0 ! Remove highest modes
             call fft(dt, dv, -1)
             data%weather_temp(:,j,i) = data%weather_temp(:,j,i) + data%tod(:,k,j,i) / dt(1:nsamp)
             do l = 1, nsamp
                write(58,*) l, data%tod(l,k,j,i)/dt(l)
                data%tod(l,k,j,i) = data%tod(l,k,j,i)/dt(l)
             end do
             write(58,*) 
          end do
          close(58)
          data%weather_temp(:,j,i) = data%weather_temp(:,j,i) / nfreq
          
          ! Compute mean TOD over all frequencies
          !g = mean(data%tod(:,15,j,i))/mean(data%tod(:,17,j,i))
             !do k = 1, nsamp
             !dt(k) = mean(data%tod(k,:,j,i))
             !dt(k) = data%tod(k,15,j,i)-data%tod(k,17,j,i)*g
          !end do
          dt(1:nsamp)            = data%weather_temp(:,j,i)
          dt(2*nsamp:nsamp+1:-1) = dt(1:nsamp)

          open(58,file='weather_raw.dat')
          do k = 1, nsamp
             write(58,*) k, dt(k)
          end do
          close(58)


          ! Fit 1/f noise profile
          call fit_1overf_profile(data%samprate, data%scanfreq, scanmask_width, sigma0, &
               & alpha, fknee, chisq_out=chisq, apply_scanmask=scanmask, tod=dt)
          write(*,fmt='(3i8,f10.2,3f8.3)') info%id, i, j, sigma0, alpha, fknee, chisq



          ! Wiener filter time stream with frequency-averaged noise parameters
          call fft(dt, dv, 1)
          do k = 1, n-1 ! Frequency 0 is unchanged
             nu    = ind2freq(k+1, data%samprate, n)
             dv(k) = (1.d0 / (1.d0 + (nu/fknee)**(-alpha)))**5 * dv(k)
          end do
          call fft(dt, dv, -1)

          open(58,file='weather_wf.dat')
          do k = 1, nsamp
             write(58,*) k, dt(k)
          end do
          close(58)

          do k = 1, nfreq
             dt(1:nsamp) = data%tod(:,k,j,i)-data%weather_temp(:,j,i)
             dt(2*nsamp:nsamp+1:-1) = dt(1:nsamp)
             call fit_1overf_profile(data%samprate, data%scanfreq, scanmask_width, sigma0, &
                  & alpha, fknee, chisq_out=chisq, apply_scanmask=scanmask, tod=dt)
             write(*,fmt='(3i8,f10.2,3f8.3)') info%id, i, j, sigma0, alpha, fknee, chisq
          end do

          call mpi_finalize(ierr)
          stop
       end do
    end do
    
  end subroutine calc_weather_template


  subroutine calc_point(data, isys, osys)
    implicit none
    type(lx_struct) :: data
    integer(i4b) :: i, j, nsamp, mod, nmod, isys, osys, ndet
    real(dp)     :: op(3), np(3), mat(3,3)
    nsamp = size(data%tod,1)
    ndet  = size(data%tod,4)
    allocate(data%point(3,nsamp,ndet))
    if (nsamp == 0) return
    do j = 1, ndet
       do i = 1, nsamp
          op = data%point_cel(:,i,j)
          np = op
!      call swap_coordinate_convention(op(1), op(2), op(3), isys)
!      call coord_convert(isys, op(1), op(2), op(3), osys, np(1), np(2), np(3), &
!       & mjd=data%time(i), euler=mat)
          data%point(:,i,j) = np
       end do
    end do
    data%coord_sys = osys

    ! Find map extent
    do j = 1, ndet
       call make_angles_safe(data%point(1,:,j),real(360.d0,sp))
    end do
    data%point_lim(1) = minval(data%point(1,:,:))
    data%point_lim(2) = maxval(data%point(1,:,:))
    data%point_lim(3) = minval(data%point(2,:,:))
    data%point_lim(4) = maxval(data%point(2,:,:))
    do j = 1, ndet
       data%point(1,:,j)   = mod(data%point(1,:,j),360.)
    end do
  end subroutine

  subroutine calc_scanfreq(data)
    implicit none
    type(lx_struct) :: data
    if (size(data%tod,1) == 0) return
    call get_scanfreq(data%point_tel(1,:,1), data%samprate, data%scanfreq(1))
    call get_scanfreq(data%point_tel(2,:,1), data%samprate, data%scanfreq(2))
  end subroutine calc_scanfreq

  subroutine fit_noise(data, powspecs, snum)
    implicit none
    type(lx_struct) :: data
    integer(i4b) :: ndet, nsb, nfreq, i, j, k, l, snum, p
    character(len=4) :: myid_text
    real(sp)     :: powspecs(:,:,:,:)
    real(dp)     :: chisq, nsamp, scale
    type(comap_scan_info) :: scan
    nsamp  = size(data%tod,1)
    nfreq  = size(data%tod,2)
    nsb    = size(data%tod,3)
    ndet   = size(data%tod,4)
    allocate(data%sigma0(nfreq,nsb,ndet), data%alpha(nfreq,nsb,ndet), data%fknee(nfreq,nsb,ndet))
    if (nsamp == 0) then
       data%sigma0 = 0.d0
       data%alpha  = 0.d0
       data%fknee  = 0.d0
       return
    end if
    call get_scan_info(snum, scan)
    call int2string(info%id, myid_text)
    scale = 1.d0 / sqrt(abs(data%nu(1,1,1)-data%nu(2,1,1))*1d9/data%samprate)
    do i = 1, ndet
!    do i = 2, 2 
       if (.not. is_alive(i)) then
          data%sigma0(:,:,i) = 0.d0
          data%alpha(:,:,i)  = 0.d0
          data%fknee(:,:,i)  = 0.d0
          cycle
       end if
       do k = 1, nsb
!       do k = 1, 1
          do j = 1, nfreq      
!          do j = 56, 56
             if (data%gain(1,j,k,i) == 0) then
                if (j == nfreq/2) write(*,fmt='(a,i4,i3,i6,a)') scan%id, i, k, j, ' -- zero Tsys'
                cycle
             end if
             if (.false.) then
                ! White noise 
                data%sigma0(j,k,i) = sqrt(variance(data%tod(:,j,k,i)))
                data%alpha(j,k,i)  = -10.d0
                data%fknee(j,k,i)  = 1d-6
                chisq = (sum((data%tod(:,j,k,i)/data%sigma0(j,k,i))**2)-nsamp)/sqrt(2.d0*nsamp)
                write(*,fmt='(4i8,e10.2,3f10.3)') info%id, i, j, k, data%sigma0(j,k,i)
             else
                call fit_1overf_profile(data%samprate, data%scanfreq, scanmask_width, data%sigma0(j,k,i), &
                     & data%alpha(j,k,i), data%fknee(j,k,i), tod_ps=real(powspecs(:,j,k,i),dp), limits=[0.01d0,5.d0], &
                     & snum=scan%sid, frequency=j, detector=i, chisq_out=chisq, apply_scanmask=scanmask, fit_par=[.true.,.false.])
                if (j == nfreq/2) write(*,fmt='(a,i4,i3,i6,f10.3,2f8.3,f8.1,f8.3,a)') scan%id, i, k, j, &
                     & data%sigma0(j,k,i)/scale/data%gain(1,j,k,i), data%alpha(j,k,i), &
                     & data%fknee(j,k,i), data%gain(1,j,k,i), chisq, '   '//trim(scan%object)
!!$                write(*,*) scan%id, i, k, j, &
!!$                     & data%sigma0(j,k,i)/scale/data%gain(1,j,k,i), data%alpha(j,k,i), &
!!$                     & data%fknee(j,k,i), data%gain(1,j,k,i), chisq, '   '//trim(scan%object)
!!$                stop
             end if
          end do
       end do
    end do

    if (data%polyorder >= 0) then
       
       p = data%polyorder
       allocate(data%sigma0_poly(0:p,nsb,ndet), data%alpha_poly(0:p,nsb,ndet), data%fknee_poly(0:p,nsb,ndet))
       do i = 1, ndet
          if (.true. .or. .not. is_alive(i)) then
             data%sigma0_poly(:,:,i) = 0.d0
             data%alpha_poly(:,:,i)  = 0.d0
             data%fknee_poly(:,:,i)  = 0.d0
             cycle
          end if
          do k = 1, nsb
             do j = 0, p
                if (.false.) then
                   data%sigma0_poly(j,k,i) = sqrt(variance(data%tod_poly(:,j,k,i)))
                   data%alpha_poly(j,k,i)  = -10.d0
                   data%fknee_poly(j,k,i)  = 1d-6
                   chisq = (sum((data%tod_poly(:,j,k,i)/data%sigma0_poly(j,k,i))**2)-nsamp)/sqrt(2.d0*nsamp)
                   write(*,fmt='(4i8,e10.2,3f8.3)') info%id, i, j, k, data%sigma0_poly(j,k,i)
                else
                   call fit_1overf_profile(data%samprate, data%scanfreq, scanmask_width, data%sigma0_poly(j,k,i), &
                        & data%alpha_poly(j,k,i), data%fknee_poly(j,k,i), tod=real(data%tod_poly(:,j,k,i),dp), &
                        & snum=scan%sid, frequency=j, detector=i, chisq_out=chisq, apply_scanmask=scanmask)
                   write(*,fmt='(a,i4,i3,i6,3f8.3,f8.3,a)') scan%id, i, k, j, data%sigma0_poly(j,k,i)/scale, data%alpha_poly(j,k,i), &
                     & data%fknee_poly(j,k,i), chisq, '   '//trim(scan%object)//' (no poly)'
                end if
             end do
          end do
       end do

    end if

  end subroutine fit_noise


  subroutine get_scanfreq(point, samprate, scanfreq)
    implicit none
    real(sp)                                :: point(:)
    real(dp)                                :: scanfreq, samprate
    real(sp),     dimension(:), allocatable :: tod, pow
    complex(spc), dimension(:), allocatable :: ft
    integer(i4b)        :: n, m, i, i_min, i_max
    n = size(point) ! min(10000, size(az))
    m = n/2+1
    allocate(tod(n), ft(m), pow(m))
    tod = point(1:n) - mean(point(1:n))
    call fft(tod, ft, 1)
    pow = ft*conjg(ft)
    i_min = freq2ind(scanfreq_min, samprate, m)
    i_max = freq2ind(scanfreq_max, samprate, m)
    i = maxloc(pow(i_min:i_max),1)+1+i_min
    scanfreq = ind2freq(i, samprate, m)
    deallocate(tod, ft, pow)
  end subroutine

  subroutine calc_gain(data)
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: n, ndet, nsb, nfreq, i, j, k, l, m, tmin, tmax, parfile_time
    real(dp)        :: g, a, sigma0, chisq, tsys, tau, dnu, const
    real(dp), dimension(:), allocatable :: el, dat

    ! Constants used for calibration with radiometer equation
    ! g = sigma0*sqrt()/T_0
    tsys = 35.d6
    tau = 1./data%samprate
    dnu = data%nu(2,1,1)-data%nu(1,1,1)
    dnu = 2.d9/1024.
    const = sqrt(tau*dnu)/tsys
    write(*,*) tsys, 'T_sys'    
    write(*,*) tau, 'tau'
    write(*,*) dnu, 'delta_nu'
    write(*,*) const, 'const'


!    m=a*data%samprate/data%scanfreq(2)      ! Number of tod-samples per gain estimate
    m   = (size(data%time))
    write(*,*) 'Number of samples per gain estimate    =', m
!    n   = (size(data%time)+m-1)/m           ! Number of gain samples
    n   = (size(data%time))/m               ! Number of gain samples
    
    write(*,*) 'Number of gain estimates for this scan =', n
    write(*,*) 'n*m                                    =', n*m
    write(*,*) 'Total number of samples                =', size(data%time)
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)
    write(*,*) nfreq, '= nfreq', nsb, '= nsb', ndet, '= ndet'
    write(*,*) '---------------------------------------------------------'
    allocate(data%time_gain(n), data%gain(n,nfreq,nsb,ndet))
    allocate(el(m), dat(m))
    data%time_gain = data%time(::m) ! OBS may reallocate length of time_gain!!!
    !open(13,file='gain.dat')
    !open(14,file='chisq.dat')
    do k = 1, ndet
       do l = 1, nsb
          do j = 1, nfreq
             sigma0 = data%sigma0(j,l,k)
             !if (sigma0 == 0) write(*,*) 'sigma0 =', real(sigma0,sp), k, '= det', l, '=sb', j, '= freq'
             g = sigma0*const
             !write(*,*) 'gain =', g
             data%gain(1,j,l,k) = g
!          do i = 1, n
!             tmin = (i-1)*m+1
!             tmax = i*m
!             el  = data%point_tel(2,tmin:tmax)
!             dat = data%tod(tmin:tmax,j,k)
!             write(*,*) tmin, tmax, 'min max'
!             call estimate_gain(el,dat,g,sigma0,chisq)
!             data%gain(i,j,k) = g
!             write(13,*) g
!             write(13,*) (g-5.6d10)/5.6d10*100
!             write(14,*) chisq
!          end do
          end do
       end do
    end do
    !close(13)
    !close(14)
    !deallocate(el, dat)
  end subroutine

!!$  subroutine calc_stats(data)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    integer(i4b)    :: i, j, k, m, n
!!$    real(dp)        :: mjd
!!$    data%stats = NaN
!!$    mjd        = (data%time(1)+data%time(size(data%time)))/2
!!$    data%stats(STAT_MJD)             = mjd
!!$    data%stats(STAT_LST)             = mjd2lst(mjd, QUIET_GEODETIC_LONGITUDE)
!!$    data%stats(STAT_AZ)              = average_ang(real(data%orig_point(1,:),dp))
!!$    data%stats(STAT_EL)              = average_ang(real(data%orig_point(2,:),dp))
!!$    data%stats(STAT_DK)              = average_ang(real(data%orig_point(3,:),dp))
!!$    data%stats(STAT_PWV)             = mean(real(data%hk%apex%value(:,APEX_PWV),dp))
!!$    data%stats(STAT_PWV_CHANGE)      = sqrt(variance(real(data%hk%apex%value(:,APEX_PWV),dp)))
!!$    data%stats(STAT_HUMIDITY)        = mean(real(data%hk%apex%value(:,APEX_HUMIDITY),dp))
!!$    data%stats(STAT_HUMIDITY_CHANGE) = sqrt(variance(real(data%hk%apex%value(:,APEX_HUMIDITY),dp)))
!!$    data%stats(STAT_WIND)            = mean(real(data%hk%apex%value(:,APEX_WIND_SPEED),dp))
!!$    data%stats(STAT_T_AMBIENT)       = mean(real(data%hk%apex%value(:,APEX_TEMPERATURE),dp))
!!$    data%stats(STAT_TENC)            = mean(real(data%hk%encl%value(:,P3T09),dp))
!!$    data%stats(STAT_TENC_CHANGE)     = sqrt(variance(real(data%hk%encl%value(:,P3T09),dp)))
!!$    data%stats(STAT_CRYO)            = mean(real(data%hk%cryo%value(:,P2T19),dp))
!!$    data%stats(STAT_CRYO_CHANGE)     = sqrt(variance(real(data%hk%cryo%value(:,P2T19),dp)))
!!$    !data%stats(STAT_BIAS_ELECTRONICS_TEMERATURE and CHANGE) = NaN
!!$    !data%stats(STAT_SUN_SIDELOBE and MOON) = NaN
!!$    !data%stats(STAT_SUN_DISTANCE and MOON) = NaN
!!$  end subroutine calc_stats
!!$
!!$  subroutine calc_diode_stats(data, powspecs)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    real(sp)        :: powspecs(:,:)
!!$    integer(i4b)    :: i, j, k, m, n, ndi, obj, pix
!!$    logical(lgt)    :: alive
!!$    real(dp)        :: mjd, avgpt(3), phi, theta, psi, va(3), v(3), b2h(3,3), mat(3,3)
!!$    real(dp)        :: t1, p1
!!$    logical(lgt), dimension(:), allocatable :: sidelobe_hit
!!$
!!$    ndi = size(data%tod(1,:))
!!$    allocate(data%diode_stats(ndi,NUM_DIODE_STATS))
!!$    data%diode_stats = NaN
!!$    do i = 1, ndi
!!$       data%diode_stats(i,:) = 0.d0
!!$       if (is_alive(i)) then
!!$          !print *, "In here is ", i
!!$          call compute_typeB_chisq(data%samprate, data%sigma0(i), data%tod(:,i), data%tp(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_TYPEB))
!!$          call compute_weather_stat(data%samprate, 10.d0, data%tp(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_WEATHER1))  ! 10 sec
!!$          call compute_weather_stat(data%samprate, 30.d0, data%tp(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_WEATHER2))  ! 30 sec
!!$          call jump_finder(1001, i, data%tod(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_JUMP)) ! check for jumps
!!$          !call tp_rms_finder(1000, i, data%tp(:,i), data%diode_stats(i,DIODE_STAT_TP_RMS))
!!$          call compute_single_fft_chisq(data%samprate, 9.9d0, 10.1d0, data%sigma0(i), powspecs(:,i), &
!!$               & data%diode_stats(i,DIODE_STAT_10HZ))
!!$          call compute_single_fft_chisq(data%samprate, 1.19d0, 1.21d0, &
!!$               & data%sigma0(i), powspecs(:,i), data%diode_stats(i,DIODE_STAT_1_2HZ))
!!$          call compute_single_fft_chisq(data%samprate, data%scanfreq-0.001d0, data%scanfreq+0.001d0, &
!!$               & data%sigma0(i), powspecs(:,i), data%diode_stats(i,DIODE_STAT_SSS))
!!$          ! Check biases for non-gaussianity/outliers
!!$          !call compute_bias_stat(data%hk%bias%value, data%diode_stats(i,DIODE_STAT_BIAS)) 
!!$          data%diode_stats(i,DIODE_STAT_SIGMA0)   = data%sigma0(i)
!!$          data%diode_stats(i,DIODE_STAT_ALPHA)    = data%alpha(i)
!!$          data%diode_stats(i,DIODE_STAT_FKNEE)    = data%fknee(i)
!!$          data%diode_stats(i,DIODE_STAT_GAIN)     = real(sum(data%gain(:,i),1)/size(data%gain(:,i),1),dp)
!!$       end if
!!$    end do
!!$    ! We also want some sidelobe stats: Sidelobe elevation and
!!$    ! sidelobe sun hit. The sun is assumed to be the first object
!!$    ! in the object-relative coordinates.
!!$
!!$    ! Find the sidelobe hits
!!$    allocate(sidelobe_hit(0:size(quiet_horns)-1))
!!$    sidelobe_hit = .false.
!!$    data%diode_stats(:,DIODE_STAT_SIDELOBE_HIT) = 0
!!$    obj = 1
!!$    do j = 1, size(data%point_objrel,2)
!!$       pix = ang2pix(sidelobes%nside, sidelobes%order, &
!!$        & real(data%point_objrel(2,j,obj),dp), real(data%point_objrel(1,j,obj),dp))
!!$       k = lookup_sidelobe_range(sidelobes, data%time(j))
!!$       sidelobe_hit = sidelobe_hit .or. sidelobes%masks(pix,:,k)
!!$    end do
!!$    where(sidelobe_hit(quiet_diodes%horn)) data%diode_stats(:,DIODE_STAT_SIDELOBE_HIT) = 1
!!$    deallocate(sidelobe_hit)
!!$
!!$    ! Find sidelobe horizontal pointing
!!$    call calc_average_pointing(data, avgpt)
!!$    call swap_coordinate_convention(avgpt(1), avgpt(2), avgpt(3))
!!$    k   = lookup_sidelobe_range(sidelobes, data%time(1))
!!$    b2h = rot_boresight2hor(data%time(1), avgpt(1), avgpt(2), avgpt(3), -1, -1)
!!$    data%diode_stats(:,DIODE_STAT_SIDELOBE_AZ) = 0
!!$    data%diode_stats(:,DIODE_STAT_SIDELOBE_EL) = 0
!!$    do mod = 0, size(quiet_horns)-1
!!$       ! For each pixel in the sidelobe, transform to telescope coordinates
!!$       if(isnan(sidelobes%centers(1,mod,k))) then
!!$          phi = nan; theta = nan
!!$       else
!!$          call vec2ang(sidelobes%centers(:,mod,k), theta, phi)
!!$          mat = matmul(b2h,angles2rot(phi, theta, 0d0))
!!$          call rot2angles(mat, phi, theta, psi)
!!$          call swap_coordinate_convention(phi, theta, psi)
!!$       end if
!!$       data%diode_stats(quiet_horns(mod)%diodes,DIODE_STAT_SIDELOBE_EL) = theta
!!$       data%diode_stats(quiet_horns(mod)%diodes,DIODE_STAT_SIDELOBE_AZ) = phi
!!$!call vec2ang(sidelobes%centers(:,mod,k), t1, p1)
!!$!write(*,'(i4,6f9.2)') mod, avgpt(1:2)*RTOD, p1*RTOD, t1*RTOD, phi*RTOD, theta*RTOD
!!$    end do
!!$  end subroutine calc_diode_stats
!!$
!!$
!!$  ! Find the parameters for the bandpass filter. At the lower end,
!!$  ! we need to worry about the scanning frequency and its multiples.
!!$  ! We want to go as low as possible while maintaining a good chisquare,
!!$  ! and also want to ensure that the filter has finished falling by the
!!$  ! time we reach a harmonic with nasty stuff in it. We will therefore
!!$  ! go in steps 0.5f, 1.5f, 2.5f, 3.5f, ...
!!$  ! We take the spikefilter into account, so that we do not cut because
!!$  ! of spikes which will be filtered anyway.
!!$  subroutine calc_bandpass_filter_par(data, powspecs)
!!$    implicit none
!!$    type(lx_struct)              :: data
!!$    real(sp)                     :: powspecs(:,:)
!!$    integer(i4b)                 :: i, j, k, numbin, ndi, numsamp, nfreq, ind1, ind2, n
!!$    integer(i4b)                 :: delta, n_smooth, delta_high, nd, delta_harm, nharm, indh
!!$    integer(i4b)                 :: nlim
!!$    integer(i4b),allocatable     :: pos(:,:), mask(:)
!!$    real(dp)                     :: dnu, chisq0, nu_min, nu_max, nu, cum, sigma_lim
!!$    real(dp)                     :: acceptlim, highlow_boundary
!!$    real(dp), allocatable, dimension(:) :: chisq, N_fft, chisq_cumdev, chisq_harm, filter, chisq_scan, spikefilter
!!$    type(filter_params) :: filter_opts
!!$!real(dp), allocatable :: cumdevs(:,:), harms(:,:)
!!$!type(hdf_file)        :: hfile
!!$
!!$    if(no_filters) return
!!$    ndi        = size(data%tod,2)
!!$    numsamp    = size(data%tod,1)
!!$    nfreq      = size(powspecs,1)
!!$    dnu        = ind2freq(2, data%samprate, nfreq)
!!$    delta_harm = 19
!!$    sigma_lim  = 4
!!$    acceptlim  = 0.75
!!$    delta_high = 100 ! Smoothing scale for high-frequency chi-squares. Should be even
!!$    ! We do not want the filters to overlap, so we restrict highpass to
!!$    ! be below highlow_boundary and lowpass to be above this.
!!$    highlow_boundary = 2.0
!!$    nlim             = freq2ind(highlow_boundary, data%samprate, nfreq)
!!$    nharm            = nlim/freq2ind(data%scanfreq, data%samprate, nfreq)-1
!!$    allocate(chisq_cumdev(nharm), chisq(nfreq), N_fft(nfreq), chisq_harm(nharm))
!!$    allocate(pos(ndi,2))
!!$    pos = -1
!!$!allocate(harms(nharm,ndi),cumdevs(nharm,ndi))
!!$    
!!$    ! Load the filter to check for spike filtering
!!$    call get_default_filter(filter_opts)
!!$    allocate(mask(nfreq))
!!$    mask = 1
!!$
!!$    if(filter_opts%apply_spike) then
!!$       allocate(spikefilter(nfreq))
!!$       spikefilter = 1.d0
!!$       call spike_filter_fft(.false., filter_opts%spikefreqs, filter_opts%num_spikes, numsamp, data%samprate, nfreq, spikefilter)
!!$       where(spikefilter < 1.d0) mask = 0 ! Mask skips all frequencies with spikes, if any
!!$    end if
!!$
!!$    do i = 1, ndi
!!$       if (data%sigma0(i) == 0) cycle
!!$       call get_inv_noise_filter_fft(.false., numsamp, data%samprate, &
!!$        & data%sigma0(i), data%fknee(i), data%alpha(i), N_fft, filter_opts%apply_highpass)
!!$       ! This is based on complex numbers, so each entry is a chisquare
!!$       ! with 2 degrees of freedom.
!!$
!!$       ! Just mask out the spikefiltered frequencies everywhere.
!!$       chisq(2:nfreq) = powspecs(2:nfreq,i) / N_fft(2:nfreq)
!!$       chisq = chisq*mask ! Masking spikefilter frequencies - if apply_spike=.false., this doesn't do anything
!!$
!!$       cum = 0
!!$       do j = nharm, 1, -1
!!$          ind1 = freq2ind((j-0.5)*data%scanfreq, data%samprate, nfreq)
!!$          ind2 = freq2ind((j+0.5)*data%scanfreq, data%samprate, nfreq)
!!$          indh = freq2ind(      j*data%scanfreq, data%samprate, nfreq)
!!$!          n   = size(chisq(ind1:nlim)) ! Version without spikefilter
!!$          n= sum(mask(ind1:nlim)) ! If there are spikefiltered frequencies, they are not counted.
!!$
!!$          if(j == nharm) then
!!$             cum = sum(chisq(ind1:nlim))   ! And I set them to zero in the chisq array, so they do not count here either
!!$          else
!!$             cum = cum + sum(chisq(ind1:ind2-1))
!!$          end if
!!$          ! Number of sigmas away from expected chisquare in gaussian approx.
!!$          ! DET ER SÅNN SOM DETTE HER ALTSÅ!
!!$          chisq_cumdev(j) = (cum - n)/sqrt(1d0*n)
!!$          ! Find the same thing just around the current harmonic
!!$          ind1  = max(1,min(nfreq,indh-delta_harm))
!!$          ind2  = max(1,min(nfreq,indh+delta_harm))
!!$
!!$          !n     = ind2-ind1+1     ! Version without spikefilter
!!$          n = sum(mask(ind1:ind2)) ! taking mask into account again
!!$          chisq_harm  (j) = (sum(chisq(ind1:ind2))-n)/sqrt(1d0*n)
!!$       end do
!!$
!!$       ! Accept the first position where the cum chisq is ok and all
!!$       ! following harm chisqs are also ok
!!$       do j = nharm, 2, -1
!!$          if(chisq_harm(j) > sigma_lim) exit
!!$       end do
!!$       do k = j, nharm-1
!!$          if(chisq_cumdev(k) < sigma_lim) exit
!!$       end do
!!$       pos(i,1) = k 
!!$!harms(:,i) = chisq_harm
!!$!cumdevs(:,i) = chisq_cumdev
!!$
!!$       ! Then handle the high frequency part of the filter
!!$       ind1   = nlim
!!$       ind2   = nfreq ! nint(12.5 / dnu)
!!$       
!!$       ! Search first for thin spikes
!!$       j = ind1
!!$       do while (chisq(j) < 15.d0 .and. j < ind2)
!!$          j = j+1
!!$       end do
!!$
!!$       ! Search for extended features
!!$       do k = max(ind1,delta_high/2+1), min(ind2,nfreq-delta_high/2)
!!$          chisq0 = sum(chisq(k-delta_high/2:k+delta_high/2))
!!$!          chisq0 = (chisq0 - (delta_high+1)) / sqrt(delta_high+1.d0) ! Version wthout spikefilter
!!$          n = sum(mask(k-delta_high/2:k+delta_high/2)) 
!!$          chisq0 = (chisq0 - n) / sqrt(1.d0*n)         
!!$          if (chisq0 > 5.d0) exit
!!$       end do
!!$       pos(i,2) = min(j,k)
!!$
!!$    end do
!!$    ! The first index (0.5*scanfreq) is extra suspicious. We only
!!$    ! accept it if most of the others also accept it.
!!$    if(count(pos(:,1)==1) < ndi*acceptlim) where(pos(:,1)==1) pos(:,1) = 2
!!$!call open_hdf_file("foo.hdf", hfile, "w")
!!$!call write_hdf(hfile, "harm", harms)
!!$!call write_hdf(hfile, "cum",  cumdevs)
!!$!call close_hdf_file(hfile)
!!$!deallocate(harms, cumdevs)
!!$!stop
!!$
!!$    ! Translate pos into filter parameters
!!$    do i = 1, ndi
!!$       if(data%sigma0(i) == 0) cycle
!!$       if (fix_highpass_freq_scan > 0.d0) then
!!$          data%filter_par(i,FILTER_HIGH_NU_SCAN) = fix_highpass_freq_scan
!!$       else
!!$          data%filter_par(i,FILTER_HIGH_NU_SCAN)= (pos(i,1)-0.5)
!!$       end if
!!$       data%filter_par(i,FILTER_HIGH_ALPHA) = -30
!!$       if (fix_lowpass_freq > 0.d0) then
!!$          data%filter_par(i,FILTER_LOW_NU)     = fix_lowpass_freq
!!$       else
!!$          data%filter_par(i,FILTER_LOW_NU)     = ind2freq(pos(i,2), data%samprate, nfreq) - 0.2d0
!!$       end if
!!$       data%filter_par(i,FILTER_LOW_ALPHA)  = -300
!!$    end do
!!$
!!$    deallocate(chisq_cumdev, chisq, N_fft, chisq_harm, pos, mask)
!!$    if(allocated(spikefilter)) deallocate(spikefilter)
!!$
!!$  end subroutine calc_bandpass_filter_par

  subroutine calc_fourier(data, ffts, powspecs)
    implicit none
    type(lx_struct) :: data
    integer(i4b)    :: nfreq, nsb, ndet, i, j, k, nf
    complex(spc), dimension(:,:,:,:), allocatable :: ffts
    real(sp),     dimension(:,:,:,:), allocatable :: powspecs
    if (size(data%tod,1) == 0) return
    nfreq = size(data%tod,2)
    nsb   = size(data%tod,3)
    ndet  = size(data%tod,4)
    nf    = size(data%tod,1)/2+1
    allocate(ffts(nf,nfreq,nsb,ndet), powspecs(nf,nfreq,nsb,ndet))
    do i = 1, ndet
       do k = 1, nsb
          call fft_multi(data%tod(:,:,k,i), ffts(:,:,k,i), 1)
       end do
    end do
    do i = 1, ndet
       do k = 1, nsb
          do j = 1, nfreq
             call extract_powspec(ffts(:,j,k,i), powspecs(:,j,k,i))
          end do
       end do
    end do
  end subroutine
!!$
!!$  function ptrans_rot(phi1, theta1, phi2, theta2) result(angle)
!!$    implicit none
!!$    real(dp)           :: phi1, phi2, theta1, theta2, c2a, s2a, r1(3), r2(3), angle
!!$    call ang2vec(theta1,phi1,r1)
!!$    call ang2vec(theta2,phi2,r2)
!!$    call qu_transport_rot(r1,r2,c2a,s2a)
!!$    angle = atan2(s2a,c2a)/2
!!$  end function
!!$
!!$  subroutine calc_average_pointing(data, avgpt)
!!$    implicit none
!!$    type(lx_struct) :: data
!!$    real(dp)        :: avgpt(:), v(3), va(3), pt(2), foo
!!$    integer(i4b)    :: i
!!$    va  = 0
!!$    foo = 0
!!$    do i = 1, size(data%orig_point,2)
!!$       pt = data%orig_point(:,i)
!!$       call swap_coordinate_convention(pt(1), pt(2), foo)
!!$       call ang2vec(pt(2), pt(1), v)
!!$       va = va + v
!!$    end do
!!$    call vec2ang(va, avgpt(2), avgpt(1))
!!$    call swap_coordinate_convention(avgpt(1), avgpt(2), foo)
!!$    avgpt(3) = average_ang(real(data%orig_point(3,:),dp))
!!$  end subroutine

end program l3gen
