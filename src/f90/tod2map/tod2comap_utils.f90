module tod2comap_utils
  use comap_lx_mod
  use cholesky_decomposition_mod
  use rngmod
  implicit none

  type tod_type
     real(dp)     :: samprate, Tsys
     integer(i4b) :: nsamp, ndet, nfreq, nsb, nsim
     real(dp)     :: fmin, fmax, df, mean_el, mean_az

     integer(i4b), allocatable, dimension(:)   :: feeds                ! active feeds
     real(sp), allocatable, dimension(:,:,:)   :: freqmask             ! (freq, sb, det)
     real(dp), allocatable, dimension(:)       :: t                    ! (time) 
     real(dp), allocatable, dimension(:,:,:,:) :: d, d_long, d_raw, g  ! (time, freq,  sb, det)
     real(dp), allocatable, dimension(:,:,:,:) :: d_sim, d_raw_sim    ! (time, freq,  sb, det)
     real(dp), allocatable, dimension(:,:,:)   :: rms                    ! (freq,  sb, det)
     real(dp), allocatable, dimension(:,:,:)   :: rms_sim                ! (freq,  sb, det)   
     real(dp), allocatable, dimension(:,:,:)   :: sigma0, fknee, alpha,f ! (freq, sb, det)
     real(dp), allocatable, dimension(:,:)     :: pixel               ! (sb, freq) or (time, det)
     real(dp), allocatable, dimension(:,:,:)   :: point, point_tel     ! (det, 3, time)

  end type tod_type

contains

  subroutine get_tod(l2file, tod, parfile)
    implicit none
    character(len=*), intent(in)    :: l2file, parfile
    type(tod_type),   intent(inout) :: tod

    integer(i4b) :: i, j, k, l
    real(dp)     :: nu_cut
    logical(lgt) :: hifreq
    type(lx_struct) :: data

    nu_cut = 0.1d0

    ! Read data
    call read_l2_file(l2file, data)
    call free_tod_type(tod)

    call get_parameter(0, parfile, 'APPLY_HIGHPASS_FILTER', par_lgt=hifreq)
 
    tod%samprate = 50.d0 !data%samprate
    tod%nsamp = size(data%time)
    tod%nfreq = size(data%nu,1)
    tod%nsb   = size(data%tod,3)
    tod%ndet  = size(data%tod,4)
    
    !write(*,*) tod%nsamp, tod%nfreq, tod%nsb, tod%ndet
    allocate( tod%t(tod%nsamp), tod%f(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%point(3,tod%nsamp,tod%ndet), tod%pixel(tod%nsamp, 19), &
         & tod%point_tel(3,tod%nsamp, tod%ndet), &
         & tod%d_raw(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%d(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%g(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%rms(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%feeds(tod%ndet) )

    allocate( tod%sigma0(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%fknee(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%alpha(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%freqmask(tod%nfreq, tod%nsb, tod%ndet))

    tod%t = data%time; tod%f = data%nu
    tod%point = data%point_cel ! call make_angles_safe(tod%point(1,:),maxang)
    tod%point_tel = data%point_tel
    !tod%g     = data%gain
    tod%feeds  = data%pixels
    tod%sigma0 = data%sigma0
    tod%fknee  = data%fknee
    tod%alpha  = data%alpha
    tod%freqmask = data%freqmask
    tod%mean_el = mean(data%point_tel(2,:,1)) ! Mean boresight
    tod%mean_az = mean(data%point_tel(1,:,1)) ! Mean boresight
    !write(*,*) shape(data%point), tod%nsamp
    tod%pixel = -200

    do k = 1, tod%ndet
       do l = 1, tod%nsb
          do j = 1, tod%nfreq
             if (tod%freqmask(j,l,k) == 0) cycle
             !do j = 6, 6
             do i = 1, tod%nsamp
                tod%d_raw(i,j,l,k) = data%tod(i,j,l,k)
             end do
             ! Apply high pass filter
             tod%d(:,j,l,k) = tod%d_raw(:,j,l,k) !- 1.d0 ! remove at some point
             !tod%d(:,j,l,k) = tod%d(:,j,l,k) - mean(tod%d(:,j,l,k))
             if (hifreq) call hp_filter(nu_cut, tod%d(:,j,l,k),tod%samprate)

             ! Estimate RMS
             tod%rms(j,l,k) = sqrt(variance(tod%d(:,j,l,k)))
          end do
       end do
    end do

    call free_lx_struct(data)

  end subroutine get_tod

  subroutine get_sim(l2file, tod, parfile, rng_handle)
    implicit none
    character(len=*), intent(in)    :: l2file, parfile
    type(tod_type),   intent(inout) :: tod
    type(planck_rng), intent(inout) :: rng_handle

    integer(i4b) :: i, j, k, l
    real(dp)     :: nu_cut
    logical(lgt) :: hifreq
    type(lx_struct) :: data

    nu_cut = 0.1d0

    ! Read data
    call read_l2_file(l2file, data)
    call free_tod_type(tod)

    call get_parameter(0, parfile, 'APPLY_HIGHPASS_FILTER', par_lgt=hifreq)
 
    tod%samprate = 50.d0 !data%samprate
    tod%nsamp = size(data%time)
    tod%nfreq = size(data%nu,1)
    tod%nsb   = size(data%tod,3)
    tod%ndet  = size(data%tod,4)
    
    !write(*,*) tod%nsamp, tod%nfreq, tod%nsb, tod%ndet
    allocate( tod%t(tod%nsamp), tod%f(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%point(3,tod%nsamp,tod%ndet), tod%pixel(tod%nsamp, 19), &
         & tod%point_tel(3,tod%nsamp, tod%ndet), &
         & tod%g(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%d_raw_sim(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%d_sim(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%rms_sim(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%feeds(tod%ndet) )

    allocate( tod%sigma0(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%fknee(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%alpha(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%freqmask(tod%nfreq, tod%nsb, tod%ndet))

    ! Get simulated data
    call simulate_tod(data, tod%d_raw_sim, parfile, rng_handle)

    tod%t = data%time; tod%f = data%nu
    tod%point = data%point_cel ! call make_angles_safe(tod%point(1,:),maxang)
    tod%point_tel = data%point_tel
    !tod%g     = data%gain
    tod%feeds  = data%pixels
    tod%sigma0 = data%sigma0
    tod%fknee  = data%fknee
    tod%alpha  = data%alpha
    tod%freqmask = data%freqmask
    tod%mean_el = mean(data%point_tel(2,:,1)) ! Mean boresight
    tod%mean_az = mean(data%point_tel(1,:,1)) ! Mean boresight
    !write(*,*) shape(data%point), tod%nsamp
    tod%pixel = -200

    do k = 1, tod%ndet
       do l = 1, tod%nsb
          do j = 1, tod%nfreq
             if (tod%freqmask(j,l,k) == 0) cycle             
             ! Apply high pass filter
             tod%d_sim(:,j,l,k) = tod%d_raw_sim(:,j,l,k)
             !tod%d(:,j,l,k) = tod%d(:,j,l,k) - mean(tod%d(:,j,l,k))
             if (hifreq) call hp_filter(nu_cut, tod%d_sim(:,j,l,k),tod%samprate)
             
             ! Estimate RMS
             tod%rms_sim(j,l,k) = sqrt(variance(tod%d_sim(:,j,l,k)))

          end do
       end do
    end do

    call free_lx_struct(data)

  end subroutine get_sim



  subroutine hp_filter(nu_cut, tod_d, samp_rate)
    implicit none   
    real(dp),                    intent(in)    :: nu_cut, samp_rate
    real(dp),     dimension(1:), intent(inout) :: tod_d
    
    integer(i4b) :: n, ind_cut
    complex(dpc), allocatable, dimension(:) :: tod_fft

    n = size(tod_d)/2 + 1
    allocate(tod_fft(n))
    ind_cut = freq2ind(nu_cut, samp_rate, n) ! n or nsamps??????
    call fft(tod_d, tod_fft, 1)
    tod_fft(1:ind_cut) = 0.d0
    call fft(tod_d, tod_fft, -1)
    
    deallocate(tod_fft)

  end subroutine hp_filter


  subroutine simulate_tod(data, tod_sim, parfile, rng_handle)
    implicit none 
  
    type(Lx_struct),  intent(in)                    :: data
    character(len=*), intent(in)                    :: parfile
    real(dp), dimension(:,:,:,:), intent(inout)     :: tod_sim 
    type(planck_rng), intent(inout)                 :: rng_handle
    type(hdf_file)                                  :: file
  
    real(dp),     allocatable, dimension(:,:)       :: cholesky_of_corr      ! Cholesky decomp of correlation matrix       
    real(dp),     allocatable, dimension(:)         :: z                     ! Array with gaussian distributed random values
    real(dp),     allocatable, dimension(:,:,:,:)   :: data_corr             ! Correlation matrix of data 
    real(dp),     allocatable, dimension(:,:,:,:)   :: cholesky_of_data_corr ! Cholesky decomposition of correlatio nmatrix of data
    real(dp),     allocatable, dimension(:,:)       :: data_current
    real(dp),     allocatable, dimension(:,:)       :: A
    real(dp),     allocatable, dimension(:,:)       :: B

    character(len=512)   :: corrmatrixfile

    integer(i4b) :: n_freq, n_samples, n_bands, n_feeds, i, j, k, l, m, o, p, x, y, n_sim, brute_force
    real(dp)     :: s, x_bar, y_bar, x_std, y_std

    call get_parameter(0, parfile, 'N_NOISE_SIMULATIONS',       par_int=n_sim) 
    call get_parameter(0, parfile, 'CORR_MATRIX_LOC',           par_string=corrmatrixfile)
    call get_parameter(0, parfile, 'BRUTE_FORCE_SIM',           par_int=brute_force)

    n_samples    = size(data%tod,1)      ! Number of time samples
    n_freq       = size(data%freqmask,1) ! Number of frequency channels   
    n_bands      = size(data%tod,3)      ! Number of side-bands  
    n_feeds      = size(data%tod,4)      ! Number of feeds     

    allocate(cholesky_of_corr(n_freq, n_freq))
    allocate(z(n_freq))

    allocate(data_corr(n_freq, n_freq, n_bands, n_feeds))
    allocate(cholesky_of_data_corr(n_freq, n_freq, n_bands, n_feeds))

    allocate(data_current(n_samples, n_freq))
    allocate(A(n_freq, n_freq))
    allocate(B(n_freq, n_freq))

    ! ------------ SIMULATIONS FROM PREDEFINED CORR MATRICES -----------    
    if (brute_force == 0) then
       ! Reading in cholesky decomposition of correlation matrix       
       call open_hdf_file(corrmatrixfile, file, "r")
       if (data%polyorder == 1) then
          call read_hdf(file, "cholesky1", cholesky_of_corr)   ! 1. order polynom
       else if (data%polyorder == 0) then
          call read_hdf(file, "cholesky0", cholesky_of_corr)   ! 0. order polynom  
       else if (data%polyorder == -1) then
          call read_hdf(file, "cholesky_1", cholesky_of_corr)   ! Polyfilter turned off 
       end if
       call close_hdf_file(file)

       ! Transpose cholesky decomposition because Python and Fortran have opposite array indexing 
       cholesky_of_corr = transpose(cholesky_of_corr)

       do k=1, n_feeds
          do j=1, n_bands
             do i=1, n_samples
                
                do o=1, n_freq
                   z(o) = rand_gauss(rng_handle)
                end do

                tod_sim(i,:,j,k) = matmul(cholesky_of_corr,z) * data%sigma0(:,j,k) 

             end do
          end do
       end do

    ! --------- BRUTE-FORCE SIMULATIONS FROM DATA CORR MATRICES ----------- 
    else

       ! Calculating correlation matrix for all bands and feeds 

       do j=1, n_feeds
          if (.not. is_alive(j)) cycle
          do k=1, n_bands
             data_current = data%tod(:,:,k,j)
             do x=1, n_freq
                if (data%freqmask(x,k,j) == 0.d0) cycle
                do y=1, n_freq
                   if (data%freqmask(y,k,j) == 0.d0) cycle
                   x_bar = sum(data_current(:,x))/n_samples
                   y_bar = sum(data_current(:,y))/n_samples

                   s = 0.d0
                   do i=1, n_samples
                      s = s + (data_current(i,x) - x_bar)*(data_current(i,y) - y_bar)
                   end do

                   x_std = sqrt( sum((data_current(:,x) - x_bar)**2) / (n_samples-1) )
                   y_std = sqrt( sum((data_current(:,y) - y_bar)**2) / (n_samples-1) )

                   data_corr(x,y,k,j) = s/(n_samples - 1) / sqrt( x_std**2 * y_std**2 )
                end do
             end do
          end do
       end do

       ! Calculating cholesky decomposition of the correlation matrix for all bands and feeds 
       do j=1, n_feeds
          if (.not. is_alive(j)) cycle
          do k=1, n_bands
             A = data_corr(:,:,k,j)
             call cholesky_decomposition(A, n_freq, B)
             cholesky_of_data_corr(:,:,k,j) = B
          end do
       end do
       ! Calculating simulated data using the cholesky decompositions       
       do k=1, n_feeds
          do j=1, n_bands
             do i=1, n_samples
                
                do o=1, n_freq
                   z(o) = rand_gauss(rng_handle)
                end do
                
                tod_sim(i,:,j,k) = matmul(cholesky_of_data_corr(:,:,j,k), z) * data%sigma0(:,j,k)

             end do
          end do
       end do
    end if

  end subroutine simulate_tod


  subroutine free_tod_type(tod)
    implicit none
    type(tod_type), intent(inout) :: tod 
    if (allocated(tod%t))         deallocate(tod%t)
    if (allocated(tod%f))         deallocate(tod%f)
    if (allocated(tod%d))         deallocate(tod%d)
    if (allocated(tod%d_raw))     deallocate(tod%d_raw)
    if (allocated(tod%g))         deallocate(tod%g)
    if (allocated(tod%rms))       deallocate(tod%rms)
    if (allocated(tod%point))     deallocate(tod%point)
    if (allocated(tod%point_tel)) deallocate(tod%point_tel)
    if (allocated(tod%sigma0))    deallocate(tod%sigma0)
    if (allocated(tod%fknee))     deallocate(tod%fknee)
    if (allocated(tod%alpha))     deallocate(tod%alpha)
    if (allocated(tod%pixel))     deallocate(tod%pixel)
    if (allocated(tod%freqmask))  deallocate(tod%freqmask)
    if (allocated(tod%feeds))     deallocate(tod%feeds)

    ! Simulated data 
    if (allocated(tod%d_sim))     deallocate(tod%d_sim)
    if (allocated(tod%d_raw_sim)) deallocate(tod%d_raw_sim) 
    if (allocated(tod%rms_sim))   deallocate(tod%rms_sim)


  end subroutine free_tod_type

end module
