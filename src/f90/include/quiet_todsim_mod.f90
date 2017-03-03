module quiet_todsim_mod
  use quiet_gain_mod
  use quiet_filter_mod
  use rngmod
  use quiet_fft_mod
  use quiet_assembly_mod
  use math_tools
  use quiet_utils
  use quiet_module_mod
  use quiet_fileutils
  use ziggurat
  implicit none

  type multimap
     integer(i4b) :: n, nnoise_per_sim, ncomp, nmap, nside, npix, order
     integer(i4b), dimension(:),     allocatable :: pixels, map2mask
     real(sp),     dimension(:,:,:), allocatable :: maps
  end type

  type sim_diode
     real(dp),     dimension(:,:,:), allocatable :: phase ! (comp,horn,samp)
     integer(i4b), dimension(:,:),   allocatable :: pix, ind
  end type

  type sim_cache
     type(sim_diode), dimension(:), allocatable :: diode_info
     integer(i4b),    dimension(:), allocatable :: comps
  end type

  ! Data structures for simulation
  type(multimap), private  :: maps
  logical(lgt),   private  :: do_add_signal, do_add_oneoverf_noise, do_add_white_noise, do_add_weather
  logical(lgt),   private  :: do_add_diode_corr

  ! Auxilliary variables
  integer(i4b),       private :: num_skies, num_sim_per_sky, oversampling
  integer(i4b),       private :: order_tot, nside_tot
  logical(lgt),       private :: initialized = .false., do_sim
  character(len=512), private :: target_name

contains

  ! Initalization routine
  subroutine initialize_todsim_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    character(len=512) :: mapdir

    if(initialized) return
    call get_parameter(0, parfile, 'ANALYZE_SIMULATED_DATA', par_lgt=do_sim)
    if(do_sim) then
       call get_parameter(0, parfile, 'ADD_SIGNAL',         par_lgt=do_add_signal)
       call get_parameter(0, parfile, 'ADD_ONEOVERF_NOISE', par_lgt=do_add_oneoverf_noise)
       call get_parameter(0, parfile, 'ADD_WHITE_NOISE',    par_lgt=do_add_white_noise)
       call get_parameter(0, parfile, 'USE_TEMPLATES',      par_lgt=do_add_weather)
       call get_parameter(0, parfile, 'ADD_DIODE_CORR',     par_lgt=do_add_diode_corr)
       call get_parameter(0, parfile, 'INPUT_SKY_MAPS',     par_string=mapdir)
       call get_parameter(0, parfile, 'NSIDE_OUT',          par_int=nside_tot)
       call get_parameter(0, parfile, 'NUM_SKY_MAPS',       par_int=num_skies)
       call get_parameter(0, parfile, 'NUM_SIM_PER_SKY',    par_int=num_sim_per_sky)
       call get_parameter(0, parfile, 'SIM_OVERSAMPLING',   par_int=oversampling, &
        & desc="How many times the double demod freq (50 Hz) to use when sampling the simulation maps.")
       call get_parameter(0, parfile, 'TARGET_NAME',        par_string=target_name)
       order_tot    = NEST

       if(do_add_signal) call read_multimap(mapdir, maps)
    end if
    initialized = .true.
  end subroutine initialize_todsim_mod

  function get_num_sim() result(n)
     implicit none
     integer(i4b) :: n
     if(.not. initialized .or. .not. do_sim) then
        n = 0
     else
        n = num_skies*num_sim_per_sky
     end if
  end function

  ! To be able to use fft_multi, we need to process more than
  ! one simulation at a time. Therefore, which is a range,
  ! and tod has one extra dimension. On the other hand, doing
  ! all simulations at once will not fit in memory, so one still
  ! needs to loop over chunks of simulations, for example 32 at a time.
  ! Note that the tod which is sent in here is assumed to be just
  ! the part relevant for the range specified in which.
  subroutine get_tod_simulation(assembly, cache, which, tod, off)
    implicit none
    type(quiet_assembly)  :: assembly
    type(sim_cache)       :: cache
    integer(i4b)          :: i, j, which(2), unit, o
    integer(i4b), optional:: off
    real(dp)              :: tod(:,:,:)  ! (numsamp, num_diodes, sim_loc)
    o = 1; if(present(off)) o = off
    tod = 0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
    !$OMP DO
    do i = 1, which(2)-which(1)+1
       j = max(1,(which(1)+i+off-3)/num_sim_per_sky+1)
       if(do_add_signal)      call add_signal     (assembly, cache, j, tod(:,:,i))
    end do
    !$OMP END PARALLEL
    if(do_add_oneoverf_noise .or. do_add_white_noise .or. do_add_weather) &
     & call add_oneoverf_noise(assembly, which, tod)
  end subroutine get_tod_simulation

  ! ***************************************************************
  !            Routines for adding individual components
  ! ***************************************************************

  subroutine add_signal(assembly, cache, which, tod)
    implicit none
    type(quiet_assembly) :: assembly
    type(sim_cache)      :: cache
    integer(i4b)         :: which, n, m, ndi, d, i, j
    real(dp)             :: tod(:,:)  ! (numsamp, num_diodes)
    real(dp),allocatable :: tod_full(:,:)
    if(size(cache%diode_info) == 0) return
    m   = assembly%decimation*oversampling
    n   = size(cache%diode_info(1)%ind,2)
    ndi = assembly%num_diodes
    allocate(tod_full(n,ndi))

    tod_full = 0
    do d = 1, ndi
       do i = 1, n
          do j = 1, size(cache%diode_info(d)%phase,2)
             tod_full(i,d) = tod_full(i,d) + sum(maps%maps(cache%diode_info(d)%ind(j,i), &
                  & cache%comps, which) * cache%diode_info(d)%phase(:,j,i))
          end do
       end do
    end do

    ! And finally decimate this
    do d = 1, ndi
       do i = 1, assembly%numsamp-1
          tod(i,d) = sum(tod_full((i-1)*m+1:i*m,d))/m
       end do
       tod(assembly%numsamp,d) = tod_full(n,d)
    end do
    deallocate(tod_full)
  end subroutine add_signal

  ! Adds oneoverf+white noise with frequency-dependent correlations
  subroutine add_oneoverf_noise(assembly, which, tod)
    implicit none
    type(quiet_assembly)  :: assembly
    integer(i4b)          :: which(2)
    real(dp)              :: tod(:,:,:) ! (numsamp, num_diodes,sims)

    real(dp)              :: srate, time, nu, t1, t2, fact, white, oneoverf, foo
    integer(i4b)          :: i,j,k, numsamp, fft_numsamp, ndi, d, m, nsim, nbin, sim
    integer(i4b),     dimension(:),     allocatable :: binmap(:)
    complex(spc),     dimension(:,:,:), allocatable :: noise_fft
    real(sp),         dimension(:,:,:), allocatable :: noise
    real(sp),         dimension(:,:),   allocatable :: wtmp
    real(dp),         dimension(:,:),   allocatable :: weathertemps
    real(dp),         dimension(:,:,:), allocatable :: L
    real(dp),         dimension(:),     allocatable :: amp
    type(planck_rng), dimension(:,:),   allocatable :: rngs

    if(.not. do_add_white_noise .and. .not. do_add_oneoverf_noise) return

    !Initialize parameters
    srate       = assembly%samprate ! QUIET sampling rate is 50 Hz, but might be downsampled
    numsamp     = assembly%numsamp
    time        = assembly%time(1)
    fft_numsamp = (numsamp/2+1)
    ndi         = assembly%num_diodes
    nsim        = which(2)-which(1)+1
    allocate(noise(numsamp,ndi,nsim),L(ndi,ndi,fft_numsamp))
    allocate(noise_fft(fft_numsamp,ndi,nsim),rngs(ndi,nsim))
    
    white    = 1; if(.not. do_add_white_noise)    white    = 0
    oneoverf = 1; if(.not. do_add_oneoverf_noise) oneoverf = 0

    ! TMR adding support for weather templates
    allocate(weathertemps(fft_numsamp,ndi))
    weathertemps = 1.d0
    if(do_add_weather) then
       call get_templates(wtmp,fft_numsamp,srate,assembly%diodes_abs)
       weathertemps = real(wtmp,dp)
    end if

    ! Precompute cholesky matrices and indexing
    k = 1
    do k = 1, fft_numsamp
        call eigen_pow(assembly%corr(k,:,:),0.5d0,L(:,:,k))
    end do

    ! Set up rngs (in order to get deterministic results, even with
    ! OMP)
    do sim = 1, nsim
       do i = 1, ndi
          call rand_init(rngs(i,sim), nint(time*2000), which(1), assembly%diodes_abs(i),sim)
       end do
    end do

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,i,nu,amp,sim,d)
    !$OMP DO COLLAPSE(2)
    ! Generate white noise in fourier
    do sim = 1, nsim
       do d = 1, ndi
          noise_fft(1,d,sim) = 0
          do i = 2, fft_numsamp-1
             noise_fft(i,d,sim) = cmplx(rand_gauss(rngs(d,sim)),rand_gauss(rngs(d,sim)))/2**0.5
          end do
          noise_fft(fft_numsamp,d,sim) = rand_gauss(rngs(d,sim))
       end do
    end do

    ! Then apply correlations and oneoverf
    allocate(amp(ndi))
    !$OMP DO
    do i = 2, fft_numsamp
       nu  = ind2freq(i, srate, fft_numsamp)
!       amp = assembly%sigma0 * (white + oneoverf*(nu/assembly%fknee)**(assembly%alpha))**0.5

       ! TMR adding weather
       amp = assembly%sigma0 * (white + oneoverf*(nu/assembly%fknee)**(assembly%alpha))**0.5 * weathertemps(i,:)**0.5
       do sim = 1, nsim
          if(do_add_diode_corr) noise_fft(i,:,sim) = matmul(L(:,:,i),noise_fft(i,:,sim))
          noise_fft(i,:,sim) = noise_fft(i,:,sim) * amp
       end do
    end do
    deallocate(amp)
    !$OMP END PARALLEL


    call fft_multi(noise, noise_fft, -1)

    !$OMP PARALLEL WORKSHARE
    tod = tod + noise
    !$OMP END PARALLEL WORKSHARE

    deallocate(noise_fft, noise, L, rngs, weathertemps)
    if(allocated(wtmp)) deallocate(wtmp)
  end subroutine add_oneoverf_noise

  subroutine test_oneoverf_noise
    implicit none
    type(quiet_assembly)    :: assembly
    integer(i4b), parameter :: n = 100000, ndi = 2, nf = n/2+1
    integer(i4b)            :: i
    real(dp)                :: tod(n,ndi,1)
    allocate(assembly%time(n), assembly%corr(nf,ndi,ndi))
    allocate(assembly%diodes_abs(ndi))
    allocate(assembly%fknee(ndi), assembly%sigma0(ndi), assembly%alpha(ndi))
    assembly%numsamp    = n
    assembly%samprate   = 25
    assembly%time       = 55500+irange(n)/24./60/60/assembly%samprate
    assembly%num_diodes = ndi
    assembly%diodes_abs = [1,2]
    assembly%sigma0     = [10.0, 1.0]
    assembly%fknee      = [ 0.1, 1.0]
    assembly%alpha      = [-2.0,-4.0]
    do i = 1, nf/2
       assembly%corr(i,:,:) = reshape([1.0,-0.3,-0.3,1.0],[2,2])
    end do
    do i = nf/2+1, nf
       assembly%corr(i,:,:) = reshape([1.0, 0.6, 0.6,1.0],[2,2])
    end do
    tod = 0
    call add_oneoverf_noise(assembly, [1,1], tod)
    open(44,file="tod.txt")
    do i = 1, n
       write(44,'(4e15.7)') (assembly%time(i)-assembly%time(1))*24*60*60, tod(i,:,1)
    end do
    close(44)
    call deallocate_quiet_assembly(assembly)
  end subroutine





  ! To avoid needlessly repeating calculations between each simulation,
  ! we use a sim_helper to store the interpolated pointing and gain
  ! as the phase and pix vectors.
  subroutine init_sim_cache(assembly, comps, cache)
    implicit none
    type(quiet_assembly),    intent(in)    :: assembly
    integer(i4b),            intent(in)    :: comps(:)
    type(sim_cache),         intent(inout) :: cache
    integer(i4b), parameter, dimension(4)  :: stokes_map = [1,2,2,3]
    integer(i4b)         :: i, j, k, l, h, c, d, m, n, nh, nc, ndi, ndec, adi, mod
    real(dp)             :: gain(2), ivec(3), ipsi, igain, a2t, x, tmp
    real(dp), allocatable:: vec(:,:,:), psi(:,:)
    call free_sim_cache(cache)
    if(.not. do_add_signal) return
    ndec = assembly%decimation * oversampling
    n    = (assembly%numsamp-1)*ndec
    ndi  = assembly%num_diodes
    nc   = size(comps)

    allocate(cache%comps(nc),cache%diode_info(ndi))
    cache%comps = comps
    do d = 1, ndi
       nh  = size(assembly%diode_info(d)%gamp)
       adi = assembly%diodes_abs(d)
       mod = assembly%diodes(d,1)
       a2t = ant2thermo(get_diode_freq(adi))
       allocate(cache%diode_info(d)%phase(nc,nh,n))
       allocate(cache%diode_info(d)%ind(nh,n))
       allocate(cache%diode_info(d)%pix(nh,n))
       allocate(vec(3,nh,2),psi(2,nh))
       do j = 1, assembly%numsamp-1
          ! Establish the boundaries for this interpolation
          gain = assembly%diode_info(d)%gain(j:j+1)*1d-9/a2t ! gain -> V/uK and antenna

          do h = 1, nh
             call ang2vec(assembly%diode_info(d)%point(2,j,  h),assembly%diode_info(d)%point(1,j,  h),vec(:,h,1))
             call ang2vec(assembly%diode_info(d)%point(2,j+1,h),assembly%diode_info(d)%point(1,j+1,h),vec(:,h,2))
             psi(:,h) = assembly%diode_info(d)%point(3,j:j+1,h)
             call make_angles_safe(psi(:,h), 2*pi)
          end do

          ! Actually perform the interpolation
          do k = 1, ndec
             l     = (j-1)*ndec+k
             x     = real(k-1,dp)/ndec
             igain = sum(gain*[1-x,x])
             do h = 1, nh
                ivec = vec(:,h,1)*(1-x)+vec(:,h,2)*x
                cache%diode_info(d)%pix(h,l) = vec2pix(maps%nside, maps%order, ivec)
                cache%diode_info(d)%ind(h,l) = maps%map2mask(cache%diode_info(d)%pix(h,l))
                ipsi = sum(psi(:,h)*[1-x,x])
                do c = 1, size(comps)
                   select case(comps(c))
                      case(1); tmp = 1
                      case(2); tmp = cos(2d0 * ipsi)
                      case(3); tmp = sin(2d0 * ipsi)
                      case(4); tmp = 1
                   end select
                   cache%diode_info(d)%phase(c,h,l) = &
                    & quiet_diodes(adi)%stokes(stokes_map(comps(c))) * &
                    & quiet_horns(mod)%amps(h) * igain * tmp
                end do
             end do
          end do
       end do
       deallocate(vec,psi)
    end do
  end subroutine init_sim_cache

  subroutine free_sim_diode(di)
    implicit none
    type(sim_diode) :: di
    if(allocated(di%phase)) deallocate(di%phase)
    if(allocated(di%pix))   deallocate(di%pix)
    if(allocated(di%ind))   deallocate(di%ind)
  end subroutine free_sim_diode

  subroutine free_sim_cache(cache)
     implicit none
     type(sim_cache) :: cache
     integer(i4b)    :: i
     if(allocated(cache%comps)) deallocate(cache%comps)
     if(allocated(cache%diode_info)) then
        do i = 1, size(cache%diode_info)
           call free_sim_diode(cache%diode_info(i))
        end do
        deallocate(cache%diode_info)
     end if
   end subroutine free_sim_cache

  ! ***************************************************************
  !            Routines for handling the multimaps
  ! ***************************************************************

  subroutine read_multimap(topdir, mmap)
    implicit none
    character(len=*)                                :: topdir
    type(multimap)                                  :: mmap
    character(len=512)                              :: line, dirname
    character(len=512), dimension(:),   allocatable :: mapnames
    integer(i4b),       dimension(:),   allocatable :: pixels
    real(dp),           dimension(:,:), allocatable :: vals
    integer(i4b)                                    :: unit, i, nside, order, n
    logical(lgt)                                    :: exist
    call free_multimap(mmap)

    dirname = trim(topdir) // "/" // trim(target_name)
    inquire(file=trim(dirname) //"/info.txt",exist=exist)
    if(.not. exist) dirname = trim(topdir) // "/default"

    ! Parse the info file
    unit = getlun()
    open(unit,file=trim(dirname) // "/info.txt")
    read(unit,*) n
    call assert(n >= num_skies, "Too few sky realizations available!")
    mmap%nmap           = num_skies
    mmap%nnoise_per_sim = num_sim_per_sky
    allocate(mapnames(mmap%nmap))
    do i = 1, mmap%nmap
       read(unit,fmt="(a)") line
       mapnames(i) = trim(dirname) // "/" // trim(line)
    end do
    close(unit)
    call assert(mmap%nmap >= 1, "No maps in multimap: " // trim(dirname))

    ! Read the first map to get the pixels etc.
    write(*,*) trim(mapnames(1))
    call read_map(vals, mmap%pixels, mmap%nside, mmap%order, trim(mapnames(1)))
    call assert(mmap%nside >= nside_tot, "Nside in sim maps smaller than NSIDE_OUT!")
  !  call assert(mmap%order == order_tot, "Order in sim maps inconsistent with run order (NEST)!")
    mmap%npix  = 12*mmap%nside**2
    mmap%n     = size(vals, 1)
    mmap%ncomp = size(vals, 2)
    allocate(mmap%maps(0:mmap%n, mmap%ncomp, mmap%nmap))
    mmap%maps(1:,:,1) = vals
    deallocate(vals)

    ! Then all the rest
    do i = 2, mmap%nmap
       write(*,*) trim(mapnames(i))
       call read_map(vals, pixels, nside, order, trim(mapnames(i)))
       call assert(size(vals,1) == mmap%n .and. size(vals,2) == mmap%ncomp &
        & .and. nside == mmap%nside .and. order == mmap%order, &
        & "Error in read multimap: Map " // trim(itoa(i)) // " is inconsistent!")
       mmap%maps(1:,:,i) = vals
       deallocate(vals, pixels)
    end do

    ! Default value is 0
    mmap%maps(0,:,:) = 0

    ! If input maps are ringed, then covert to nest
    if (mmap%order ==1) then 
       call convert_ring2nest_sparse(mmap%nside, mmap%pixels)
       mmap%order = 2
    end if

    ! And bulid map2mask
    allocate(mmap%map2mask(0:mmap%npix-1))
    mmap%map2mask = 0
    do i = 1, mmap%n; mmap%map2mask(mmap%pixels(i)) = i; end do

  end subroutine read_multimap

  subroutine free_multimap(mmap)
    implicit none
    type(multimap) :: mmap
    if(allocated(mmap%pixels))   deallocate(mmap%pixels)
    if(allocated(mmap%map2mask)) deallocate(mmap%map2mask)
    if(allocated(mmap%maps))     deallocate(mmap%maps)
  end subroutine

end module quiet_todsim_mod
