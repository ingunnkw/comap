module tod2comap_mapmaker
  use quiet_utils
  use tod2comap_utils
  use quiet_assembly_mod
  use quiet_filter_mod
  use quiet_todsim_mod
  use quiet_nr_mod
  implicit none

  ! The data used for computing a map (asside from what is in the
  ! assembly structure). Contains both the parts pol and tmp needs.
  ! Only the relevant parts are allocated.

  ! Unified version
  ! phase: (ncomp,nhorn,nsamp,ndi)

  type phase_info
     real(cov_type), dimension(:,:,:), allocatable :: val ! (ncomp,nhorn,nsamp)
  end type

  type mapinfo
     ! common part
     integer(i4b),     dimension(:),       allocatable :: map2mask
     real(dp),         dimension(:,:),     allocatable :: N_d
     type(phase_info), dimension(:),       allocatable :: phase
     ! TOD template variables
     real(rhs_type), dimension(:,:,:),     allocatable :: PNT, PNT_sq !(ntemp,npix,ncomp)
     real(dp),       dimension(:,:,:),     allocatable :: T           !(nsamp,ndi,ntemp)
     real(dp),       dimension(:,:),       allocatable :: Cinv, Cinv_sq
     integer(i4b) :: ntemp
  end type mapinfo

  logical(lgt), private :: initialized = .false.
  logical(lgt), private :: sim_only, filter_azimuth, use_precomp_azfilter, apply_spike, filter_spikes = .false., fftconv=.true.
  real(dp),     private :: az_filter_param

  integer(i4b), private :: bench_prepare, bench_rhs, bench_az, bench_azrhs, bench_sim, bench_div, bench_cov, bench_azcov
  integer(i4b)          :: bench_cov1, bench_cov2, bench_cov3, bench_cov4, bench_cov5

contains

  subroutine initialize_map_maker(parfile, info)
    implicit none
    character(len=*), intent(in) :: parfile
    type(common_info)            :: info
    if(initialized) return
    call initialize_todsim_mod(parfile)
    call get_parameter(0, parfile, "REPLACE_DATA_WITH_SIM", par_lgt=sim_only, desc=&
         & "Normally, the first component of data%map is the real data, and the rest" // &
         & "are simulations. With this to true, all are simulations.")
    call get_parameter(0, parfile, "APPLY_AZIMUTH_FILTER", par_lgt=filter_azimuth)
    call get_parameter(0, parfile, 'USE_PRECOMPUTED_FILTER', par_lgt=use_precomp_azfilter)
    if(.not. use_precomp_azfilter) then
       call get_parameter(0, parfile, 'AZIMUTH_FILTER_ORDER', par_dp=az_filter_param)
    end if
    call get_parameter(0, parfile, "APPLY_SPIKEFILTER", par_lgt=apply_spike)
    if (apply_spike) then
       call get_parameter(0, parfile, 'USE_WOODBURY_SPIKEFILTER', par_lgt=filter_spikes)
    end if
    call get_parameter(0, parfile, 'FFT_CONVOLUTION', par_lgt=fftconv, desc=&
         & "Use fourier transforms to apply the filters when possible. Rhs " // &
         & "evaluation will be slow if you disable this.")
    call init_map_maker_bench(info)

    initialized  = .true.
  end subroutine

  subroutine init_map_maker_bench(info)
    implicit none
    type(common_info)            :: info
    bench_prepare = bench_init(info%bench, "prepare")
    bench_rhs     = bench_init(info%bench, "rhs")
    bench_az      = bench_init(info%bench, "az")
    bench_azrhs   = bench_init(info%bench, "azrhs")
    bench_sim     = bench_init(info%bench, "sim")
    bench_div     = bench_init(info%bench, "div")
    bench_cov     = bench_init(info%bench, "cov")
    bench_azcov   = bench_init(info%bench, "azcov")

    bench_cov1   = bench_init(info%bench, "cov_rbuf")
    bench_cov2   = bench_init(info%bench, "cov_fill1")
    bench_cov3   = bench_init(info%bench, "cov_fill2")
    bench_cov4   = bench_init(info%bench, "cov_sample")
  end subroutine

  subroutine finalize_map_maker
  end subroutine

  ! It turned out that everything here only really ever dealt with a single
  ! set. The assignment to other sets can happen outside.
  subroutine process_assembly_data(assembly, map2mask, info, out)
    implicit none
    type(quiet_assembly), intent(inout) :: assembly
    integer(i4b),         intent(in)    :: map2mask(0:)
    type(common_info),    intent(in)    :: info
    type(mapdata),        intent(inout) :: out
    type(mapinfo)                       :: data
    integer(i4b)                        :: i

    call update_mon(info, "proc_ass_data_start")

    call bench_start(info%bench, bench_prepare)
    call prepare_data(assembly, map2mask, info, data, out);          call update_mon(info, "prepare_data")

!    write(*,*) info%myid, info%cid, "correlation length:", assembly%noise_corr_length
    call dmem("correlation length: " // trim(itoa(assembly%noise_corr_length)),2)
    call bench_stop (info%bench, bench_prepare)
    call bench_start(info%bench, bench_rhs)
    if(.not. sim_only) call build_rhs(assembly, info, data, data%N_d, out%rhs(:,:,1)); call update_mon(info, "build_rhs")
    call bench_stop (info%bench, bench_rhs)

    call bench_start(info%bench, bench_az)
    call process_azimuth_filter(assembly, info, data, out);          call update_mon(info, "build_az")
    call bench_stop (info%bench, bench_az)

    call bench_start(info%bench, bench_azrhs)
    if (.not. sim_only) call apply_azimuth_filter_to_rhs(assembly, info, data, data%N_d, out%rhs(:,:,1)); &
         & call update_mon(info, "apply_az")
    call bench_stop (info%bench, bench_azrhs)

    call bench_start(info%bench, bench_sim)
    call build_sim(assembly, info, data, out);                       call update_mon(info, "build_sim")
    call bench_stop (info%bench, bench_sim)

    call bench_start(info%bench, bench_div)
    call build_div(assembly, info, data, out);                       call update_mon(info, "build_div")
    call bench_stop (info%bench, bench_div)

    call bench_start(info%bench, bench_cov)
    call build_cov(assembly, info, data, out);                       call update_mon(info, "build_cov")
!    call build_cov_test(assembly, info, data, out);                       call update_mon(info, "build_cov")
    call bench_stop (info%bench, bench_cov)

    call bench_start(info%bench, bench_azcov)
    call correct_az_cov(assembly, info, data, out);                  call update_mon(info, "correct_az_cov")
    call bench_stop (info%bench, bench_azcov)
!    call test_mapmaker(out)

    call build_dir(assembly, info, data, out)
    call build_ang(assembly, info, data, out)
    call free_mapinfo(data)

  end subroutine process_assembly_data

  subroutine free_mapinfo(minfo)
    implicit none
    type(mapinfo) :: minfo
    integer(i4b)  :: i
    if(allocated(minfo%map2mask))      deallocate(minfo%map2mask)
    if(allocated(minfo%N_d))           deallocate(minfo%N_d)
    if(allocated(minfo%PNT))           deallocate(minfo%PNT)
    if(allocated(minfo%PNT_sq))        deallocate(minfo%PNT_sq)
    if(allocated(minfo%T))             deallocate(minfo%T)
    if(allocated(minfo%Cinv))          deallocate(minfo%Cinv)
    if(allocated(minfo%Cinv_sq))       deallocate(minfo%Cinv_sq)
    if(allocated(minfo%phase)) then
       do i = 1, size(minfo%phase)
          if(allocated(minfo%phase(i)%val)) deallocate(minfo%phase(i)%val)
       end do
       deallocate(minfo%phase)
    end if
  end subroutine


  ! Prepares data structures for mapmaking.
  subroutine prepare_data(assembly, map2mask, info, data, out)
    implicit none
    type(quiet_assembly),intent(inout) :: assembly
    type(common_info),   intent(in)    :: info
    integer(i4b),        intent(in)    :: map2mask(0:)
    type(mapinfo),       intent(out)   :: data
    type(mapdata)                      :: out
    integer(i4b) :: i
    call dmem("prepare_data_start",3)
    call setup_mapinfo(assembly, info, map2mask, data, out)
    call dmem("setup_mapinfo",3)
    call initialize_filter_assembly(assembly)
    call dmem("initialize_filter_assembly: "//trim(itoa(assembly%noise_corr_length)),3)
    if(fftconv) then
       call apply_filter_fft_matrix(assembly%N_corr_fft, assembly%tod, data%N_d)
       call dmem("apply_filter_fft_matrix",3)
    else
       call apply_filter_real_matrix(assembly%N_corr_real, assembly%tod, data%N_d)
       call dmem("apply_filter_real_matrix",3)
    end if

! ingunn her
!    open(58,file='filtered_tod.dat')
!    do i = 1, size(assembly%tod(:,1))
!       write(58,*) i, data%N_d(i,1)
!    end do
!    close(58)
!stop
    !data%N_d = assembly%tod

  end subroutine

  subroutine build_sim(assembly, info, data, out)
    implicit none
    type(quiet_assembly),     intent(in) :: assembly
    type(common_info),        intent(in) :: info
    type(mapinfo),            intent(in) :: data
    type(mapdata)                        :: out
    type(sim_cache)                      :: cache
    integer(i4b)                         :: i, j, k, off, m, n, nsim, a, b, ngroup
    real(dp), allocatable                :: tod(:,:,:), Nd(:,:,:)
!type(hdf_file) :: hfile
    off = 1; if(sim_only) off = 0
    nsim   = out%nmap-off
    ngroup = 16
    if(nsim == 0) return
    call init_sim_cache(assembly, info%comps, cache)
    do i = 0, (nsim-1)/ngroup
       a = i*ngroup+1
       b = min(nsim,(i+1)*ngroup)
       allocate(tod(assembly%numsamp,assembly%num_diodes,b-a+1))
       allocate(Nd(assembly%numsamp, assembly%num_diodes,b-a+1))
       call get_tod_simulation(assembly, cache, [a,b], tod, off=off)
!call open_hdf_file("foo.hdf", hfile, "w")
!call write_hdf(hfile, "tod", tod)
!call write_hdf(hfile, "time", assembly%time)
!call close_hdf_file(hfile)
!write(*,*) assembly%diodes_abs
       if(fftconv) then
          call apply_filter_fft_multi(assembly%N_corr_fft, tod, Nd) 
       else
          call apply_filter_real_multi(assembly%N_corr_real, tod, Nd) 
       end if
       do j = 1, b-a+1
          call build_rhs(assembly, info, data, Nd(:,:,j), out%rhs(:,:,a+j-1+off))
          call apply_azimuth_filter_to_rhs(assembly, info, data, Nd(:,:,j), out%rhs(:,:,a+j-1+off))
       end do
       deallocate(tod,Nd)
    end do
    call free_sim_cache(cache)
  end subroutine

  subroutine setup_mapinfo(assembly, info, map2mask, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    integer(i4b)         :: map2mask(0:)
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out
    type(filter_params)  :: filter
    integer(i4b)         :: nsamp, ndi, i, j, h, npix, ntemp, di, ngroups, mod, adi, n
    integer(i4b)         :: az_order, nspikes
    real(dp)             :: g, s, a2t
    real(dp),     allocatable :: rphase(:,:), fact(:)
    integer(i4b), allocatable :: stok(:)
    integer(i4b), parameter, dimension(4) :: stokes_map = [1,2,2,3]

    nsamp    = assembly%numsamp
    ndi      = assembly%num_diodes
    npix     = size(map2mask)
    n        = out%n
    az_order = maxval(assembly%az_order)
    if(.not. use_precomp_azfilter) az_order = az_filter_param
    nspikes = 0
    if(filter_spikes) then
       call get_default_filter(filter)
       nspikes = filter%num_spikes
    end if
    ntemp      = ndi*(az_order + nspikes*2)

    call free_mapinfo(data)
    allocate(data%map2mask(0:npix-1))
    allocate(data%N_d(nsamp, ndi))
    allocate(data%phase(ndi))
    allocate(data%PNT(ntemp,n,info%ncomp))
    allocate(data%T(nsamp,ndi,ntemp))
    allocate(data%Cinv(ntemp,ntemp))
    if (associated(out%cov)) allocate(data%PNT_sq(ntemp,n,info%ncomp))
    if (associated(out%cov)) allocate(data%Cinv_sq(ntemp,ntemp))
    data%ntemp    = ntemp
    data%map2mask = map2mask

    allocate(stok(info%ncomp), rphase(info%ncomp,nsamp), fact(info%ncomp))
    stok = stokes_map(info%comps)
    do di = 1, ndi
       adi     = assembly%diodes_abs(di)
       mod     = quiet_diodes(adi)%horn
       a2t     = ant2thermo(get_diode_freq(adi))
       ngroups = size(quiet_horns(mod)%groups)
       allocate(data%phase(di)%val(info%ncomp, ngroups, nsamp))
       do h = 1, ngroups
          do i = 1, info%ncomp
             select case(info%comps(i))
                case(T); rphase(i,:) = 1
                case(Q); rphase(i,:) = cos(2d0 * assembly%diode_info(di)%point(3,:,h))
                case(U); rphase(i,:) = sin(2d0 * assembly%diode_info(di)%point(3,:,h))
                case(V); rphase(i,:) = 1
             end select
          end do

          ! Weight of components from this horn, mV/K (ant.) -> V/uK (thermo)
          fact = quiet_diodes(adi)%stokes(stok) * quiet_horns(mod)%amps(h) * 1d-9 / a2t
          do j = 1, nsamp
             data%phase(di)%val(:,h,j) = rphase(:,j) * assembly%diode_info(di)%gain(j) * fact
          end do
       end do
    end do
    deallocate(stok, rphase, fact)
  end subroutine setup_mapinfo

  subroutine build_rhs(assembly, info, data, from, to)
    implicit none
    type(quiet_assembly),     intent(in) :: assembly
    type(common_info),        intent(in) :: info
    type(mapinfo),            intent(in) :: data
    real(dp),                 intent(in) :: from(:,:)
    real(rhs_type),           intent(out):: to(:,:)

    integer(i4b)              :: d, h, i, c, p

    do d = 1, assembly%num_diodes
       do i = 1, assembly%numsamp
          do h = 1, size(assembly%diode_info(d)%pix,1)
             p = data%map2mask(assembly%diode_info(d)%pix(h,i))
             if(p <= 0) cycle
             do c = 1, info%ncomp
                to(p,c) = to(p,c) + data%phase(d)%val(c,h,i) * from(i,d)
             end do
          end do
       end do
    end do
  end subroutine

  ! We ignore time correlations, but still include diode
  ! correlations. The per-sample inverse covariance matrix is
  ! assembly%inv_N_eff(ndi,ndi). We will project this down
  ! to comppix space with phase(ndi,ncomp*pix_hit) such that
  ! inv_N_qu = phase' * inv_N_eff * phase. Multiple horns
  ! would normally lead to pixel correlations, but since we
  ! ignore those, the horn loop just assigns the (ncomp,ncomp)
  ! matrix repeatedly to each horn.
  subroutine build_div(assembly, info, data, out)
    implicit none
    type(quiet_assembly)                        :: assembly
    type(common_info)                           :: info
    type(mapinfo)                               :: data
    type(mapdata)                               :: out

    integer(i4b)                                :: i, j, k, g, p
    integer(i4b),   allocatable, dimension(:)   :: horns, d
    real(dp),       allocatable, dimension(:,:) :: icov
    real(rhs_type), allocatable, dimension(:,:) :: phase

    if(.not. associated(out%div)) return
    ! Three effects lead to pixel correlations:
    !  1. Time correlations: Simulated with N_eff
    !  2. Multiple horns with their own timestreams:
    !       These can be decoupled and handled one by one,
    !       ignoring their correlations.
    !  3. Differential horns:
    !       A single time-stream comes from several points
    !       on the sky. Calculate for one horn, then
    !       distribute to each.

    ! First divide into groups of compatible pointing.
    call uniqi(quiet_diodes(assembly%diodes_abs)%horn, horns)
    do j = 1, size(horns)
       call wherei(quiet_diodes(assembly%diodes_abs)%horn == horns(j), d)
       ! Get the icov for the diodes in this horn
       allocate(icov(size(d),size(d)), phase(info%ncomp,size(d)))
       icov = assembly%N_eff(d,d)
       call invert_matrix(icov, cholesky=.true.)
       do i = 1, assembly%numsamp
          ! Loop over signal contributions for each horn. Since we
          ! already spit into horns, all diodes have the same pointing.
          do g = 1, size(quiet_horns(horns(j))%groups)
             p = data%map2mask(assembly%diode_info(d(1))%pix(g,i))
             if(p <= 0) cycle
             ! Copy needed because :%: is not allowed
             do k = 1, size(d)
                phase(:,k) = data%phase(d(k))%val(:,g,i)
             end do
             out%div(p,:,:) = out%div(p,:,:) + &
              & matmul(phase, matmul(icov, transpose(phase)))
          end do
       end do
       deallocate(icov, phase, d)
    end do
    deallocate(horns)
  end subroutine

  subroutine build_cov_test2(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out

    integer(i4b)         :: i1, i2, d1, d2, h1, h2, p1, p2, c1, c2, i, j, k, ndi, nsamp, n
    integer(i4b)         :: omp_get_thread_num, bsize, nblock, b1, b2, hmax, j1, j2
    integer(i4b)         :: r1(2), r2(2), l1, l2, locmax
    integer(i4b)         :: ngroups(assembly%num_diodes), corrlen, col
    real(cov_type)       :: N_phase(info%ncomp), N_phase_sq(info%ncomp)
    real(cov_type), allocatable                   :: buf(:,:,:,:,:)
    real(cov_type), allocatable, dimension(:,:,:) :: corr, corr_sq
    integer(i4b),   allocatable                   :: loc(:,:,:,:), loc2glob(:,:), tpix(:)
    integer(i4b),   allocatable                   :: upix(:), locn(:)

    if(.not. associated(out%cov)) return

call bench_start(info%bench, bench_cov3)
    ndi     = assembly%num_diodes
    nsamp   = assembly%numsamp
    corrlen = assembly%noise_corr_length
    bsize   = 32
    nblock  = (nsamp+bsize-1)/bsize
    do i = 1, ndi
       ngroups(i) = size(data%phase(i)%val,2)
    end do
    hmax = maxval(ngroups)

    allocate(corr(ndi,ndi,-corrlen:corrlen), corr_sq(ndi,ndi,-corrlen:corrlen))
    do i = -corrlen, corrlen
       corr(:,:,i)    = assembly%N_corr_real(i,:,:)
       corr_sq(:,:,i) = assembly%N_corr_F_sq_real(i,:,:)
    end do

    ! Set up mapping for block-local pixels
    allocate(loc(nblock,bsize,ndi,hmax),loc2glob(nblock,bsize*hmax*ndi), tpix(bsize*hmax*ndi))
    allocate(locn(nblock))
    do b1 = 1, nblock
       ! gather pixels O(nsamp)
       k = 0
       do i1 = (b1-1)*bsize+1, min(b1*bsize,nsamp)
          do d1 = 1, ndi
             do h1 = 1, ngroups(d1)
                k = k+1
                tpix(k) = data%map2mask(assembly%diode_info(d1)%pix(h1,i1))
             end do
          end do
       end do
       ! find unique O(nsamp log(bsize))
       call uniqi(tpix(:k),upix)
       locn(b1) = size(upix)
       loc2glob(b1,:size(upix)) = upix
       ! assign mapping O(nsamp*bsize)
       k = 0
       j = 0
       do i1 = (b1-1)*bsize+1, min(b1*bsize,nsamp)
          j = j+1
          do d1 = 1, ndi
             do h1 = 1, ngroups(d1)
                k = k+1
                p1 = data%map2mask(assembly%diode_info(d1)%pix(h1,i1))
                do i = 1, size(upix)
                   if(upix(i) == p1) exit
                end do
                loc(b1,j,d1,h1) = i
             end do
          end do
       end do
       deallocate(upix)
    end do
    deallocate(tpix)
    locmax = maxval(locn)
call bench_stop(info%bench, bench_cov3)

    !$OMP PARALLEL PRIVATE(b1,b2,i1,i2,j1,j2,d1,d2,h1,h2,p1,p2,c1,c2,k,N_phase,N_phase_sq,buf,col,r1,r2,l1,l2)
    allocate(buf(info%ncomp,info%ncomp,locmax,locmax,2))
    !$OMP DO SCHEDULE(guided)
    do b1 = 1, nblock
       r1 = [ (b1-1)*bsize+1, min(b1*bsize,nsamp) ]
       ! loop through all correlated blocks
       do b2 = (max(r1(1)-corrlen,1)-1)/bsize+1, (min(r1(2)+corrlen,nsamp)-1)/bsize+1
          r2 = [ (b2-1)*bsize+1, min(b2*bsize,nsamp) ]
          buf = 0
          j1 = 0
call bench_start(info%bench, bench_cov1)
          do i1 = r1(1), r1(2)
             j1 = j1+1
             do d1 = 1, ndi
                do h1 = 1, ngroups(d1)
                   l1 = loc(b1,j1,d1,h1)
                   j2 = 0
                   do i2 = r2(1), r2(2) ! can be shortened slightly
                      j2 = j2+1
                      k  = i2-i1
                      do d2 = 1, ndi
                         N_phase   = data%phase(d1)%val(:,h1,i1) * corr(d1,d2,k)
                         N_phase_sq= data%phase(d1)%val(:,h1,i1) * corr_sq(d1,d2,k)
                         do h2 = 1, ngroups(d2)
                            l2 = loc(b2,j2,d2,h2)
                            do c1 = 1, info%ncomp
                               do c2 = 1, c1
                                  buf(c2,c1,l2,l1,1) = buf(c2,c1,l2,l1,1) + &
                                   & N_phase(c1) * data%phase(d2)%val(c2,h2,i2)
                               end do
                            end do
                            do c1 = 1, info%ncomp
                               do c2 = 1, c1
                                  buf(c2,c1,l2,l1,2) = buf(c2,c1,l2,l1,2) + &
                                   & N_phase_sq(c1) * data%phase(d2)%val(c2,h2,i2)
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
call bench_stop(info%bench, bench_cov1)
call bench_start(info%bench, bench_cov2)
          do l1 = 1, locn(b1)
             p1 = loc2glob(b1,l1)
             do l2 = 1, locn(b2)
                p2 = loc2glob(b2,l2)
                if(p2 <= p1) then
                   do c1 = 1, info%ncomp
                      do c2 = 1, c1
                         !$OMP ATOMIC
                         out%cov(p2,c2,p1,c1) = out%cov(p2,c2,p1,c1) + buf(c2,c1,l2,l1,1)
                      end do
                   end do
                else
                   do c1 = 1, info%ncomp
                      do c2 = 1, c1
                         !$OMP ATOMIC
                         out%cov(p2+1,c2,p1,c1) = out%cov(p2+1,c2,p1,c1) + buf(c2,c1,l2,l1,2)
                      end do
                   end do
                end if
             end do
          end do
call bench_stop(info%bench, bench_cov2)
       end do
    end do
    !$OMP END DO

call bench_start(info%bench, bench_cov4)
    ! And finally fill in the other symmetric part of the components
    !$OMP DO
    do p1 = 1, out%n
       do p2 = 1, out%n+1
          do c1 = 1, info%ncomp
             do c2 = 1, c1-1
                out%cov(p2,c1,p1,c2) = out%cov(p2,c2,p1,c1)
             end do
          end do
       end do
    end do
call bench_stop(info%bench, bench_cov4)

    deallocate(buf)
    !$OMP END PARALLEL
    deallocate(loc,loc2glob,locn,corr)
  end subroutine

  ! Like build_cov, but assumes that the assembly is homogeneous.
  ! This means that each diode has the same properties: The same
  ! components and the same horns. This is the case for our
  ! current assemblies, but not for all possible ones.
  !
  ! This is faster than build_cov with one thread, but scales
  ! poorly, and with 4 threads is is significantly slower.
  subroutine build_cov_test(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out

    integer(i4b)         :: i1,i2,d,h,p,i,j,k,c,ndi,nsamp,ncomp,nhorn,npc
    integer(i4b)         :: omp_get_thread_num, ia,ib,ic
    integer(i4b)         :: ngroups(assembly%num_diodes), corrlen, col
    real(cov_type)       :: N_phase(info%ncomp), N_phase_sq(info%ncomp)
    real(cov_type), allocatable :: phase(:,:,:), tphase(:,:,:), invN(:,:,:,:),lcov(:,:), tcov(:,:)
    real(cov_type), allocatable :: colbuf(:,:,:), colbuf2(:,:,:,:,:)
    integer(i4b),   allocatable :: pixels(:,:), pixcomps(:,:)

    if(.not. associated(out%cov)) return

    ndi     = assembly%num_diodes
    nsamp   = assembly%numsamp
    corrlen = assembly%noise_corr_length
    ncomp   = info%ncomp
    nhorn   = size(data%phase(1)%val,2)
    npc     = ncomp*nhorn

    ! Precompute as much as possible before the loop into pix,comp-flattened
    ! arrays in memory-sensible order.
    allocate(phase(npc,ndi,nsamp), pixels(nhorn,nsamp), invN(ndi,ndi,2,-corrlen:corrlen))
    allocate(pixcomps(npc,nsamp), tphase(ndi,npc,nsamp))
    do i1 = 1, nsamp
       pixels(:,i1)  = data%map2mask(assembly%diode_info(1)%pix(:,i1))
       do d = 1, ndi
          phase(:,d,i1) = reshape(transpose(data%phase(d)%val(:,:,i1)),[npc])
       end do
       tphase(:,:,i1) = transpose(phase(:,:,i1))
       do h = 1, nhorn
          do c = 1, ncomp
             pixcomps(h+nhorn*(c-1),i1) = pixels(h,i1)+out%n*(c-1)
          end do
       end do
    end do
    do k = -corrlen, corrlen
       invN(:,:,1,k) = assembly%N_corr_real(k,:,:)
       invN(:,:,2,k) = assembly%N_corr_F_sq_real(k,:,:)
    end do

    !$OMP PARALLEL PRIVATE(lcov,colbuf,colbuf2,i1,i2,i,j,k,h,p,tcov,ia,ib,ic)
    allocate(lcov(npc,npc), colbuf2(out%n,ncomp,nhorn,npc,2))
    allocate(tcov(npc,ndi))
    !$OMP DO SCHEDULE(guided)
    do i1 = 1, nsamp
       call bench_start(info%bench, bench_cov4)
       if(all(pixels(:,i1) <= 0 .or. pixels(:,i1) < out%col_start &
        & .or. pixels(:,i1) > out%col_stop)) cycle
       call bench_start(info%bench, bench_cov1)
       colbuf2 = 0
       do j = 1, 2
          do i2 = max(i1-corrlen,1), min(i1+corrlen, nsamp)
             i = i1-i2
             !lcov(:,:,j) = matmul(phase(:,:,i1),matmul(invN(:,:,j,i),transpose(phase(:,:,i2))))
             !tcov = matmul(invN(:,:,j,i),tphase(:,:,i2))
             !lcov(:,:,j) = matmul(phase(:,:,i1),tcov)
             !tcov = (invN*tphase)' = tphase'*invN' = tphase'*invN
             do ia = 1, size(tcov,2)
             do ib = 1, size(tcov,1)
                tcov(ib,ia) = sum(tphase(:,ib,i2)*invN(:,ia,j,i))
             end do
             end do
             !and lcov = phase * invN * phase' = phase * tcov'
             do ia = 1, size(lcov,2)
             do ib = 1, size(lcov,2)
                lcov(ib,ia) = sum(phase(:,ib,i1)*tcov(:,ia))
             end do
             end do
             colbuf2(pixels(:,i2),:,:,:,j) = colbuf2(pixels(:,i2),:,:,:,j) + &
              & reshape(lcov(:,:),[nhorn,ncomp,nhorn,ncomp])
          end do
       end do
       call bench_stop(info%bench, bench_cov1)
       !! Ok, we have this column. Fill in the areas in the actual cov
       !do h = 1, nhorn
       !   p = pixels(h,i1)
       !   if(p <= 0) cycle
       !   out%cov(:p  ,:,p,:) = out%cov(:p,  :,p,:) + colbuf(1:p,:,h,:,1)
       !   out%cov(p+1:,:,p,:) = out%cov(p+1:,:,p,:) + colbuf(p: ,:,h,:,2)
       !end do
       call bench_start(info%bench, bench_cov2)
       do h = 1, nhorn
          p = pixels(h,i1)
          if(p <= 0) cycle
          do k = 1, p
             if(all(colbuf2(k,:,h,:,1)==0)) cycle
             !$OMP CRITICAL
             out%cov(k,:,p,:) = out%cov(k,  :,p,:) + colbuf2(k,:,h,:,1)
             !$OMP END CRITICAL
          end do
          do k = p, out%n
             if(all(colbuf2(k,:,h,:,2)==0)) cycle
             !$OMP CRITICAL
             out%cov(k+1,:,p,:) = out%cov(k+1,:,p,:) + colbuf2(k ,:,h,:,2)
             !$OMP END CRITICAL
          end do
       end do
       call bench_stop(info%bench, bench_cov2)
       call bench_stop(info%bench, bench_cov4)
    end do
    !$OMP END DO
    deallocate(lcov,colbuf2,tcov)
    !$OMP END PARALLEL
    deallocate(phase, invN, pixels, pixcomps, tphase)
  end subroutine

  ! build_cov_pol and build_cov_temp could in theory be combined
  ! into a single routine, with support for both multiple components
  ! and multiple horns at the same time. This would make it more
  ! general and easier to maintain.
  subroutine build_cov(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out

    integer(i4b)         :: i1, i2, d1, d2, h1, h2, p1, p2, c1, c2, i, j, k, ndi, nsamp, n
    integer(i4b)         :: omp_get_thread_num
    integer(i4b)         :: ngroups(assembly%num_diodes), corrlen, col
    real(dp)             :: t1, t2
    real(cov_type)       :: N_phase(info%ncomp), N_phase_sq(info%ncomp)
    real(cov_type), allocatable, dimension(:,:,:) :: rbuffer, rbuffer_sq, corr, corr_sq

    if(.not. associated(out%cov)) return

    call wall_time(t1)

    ndi     = assembly%num_diodes
    nsamp   = assembly%numsamp
    corrlen = assembly%noise_corr_length
    do i = 1, ndi
       ngroups(i) = size(data%phase(i)%val,2) 
    end do

    ! Use optimized version if possible. We need there to only
    ! be a single horn involved, and the components to be Q and U only.
    !  Do al diodes in the assembly belong to the same module?
    if(all(assembly%diodes(:,1) == assembly%diodes(1,1)) .and. &
     ! Does that module only get data from one horn?
     & size(quiet_horns(assembly%diodes(1,1))%groups) == 1 .and. &
     ! Do we have two components only?
     & info%ncomp == 2) then
     ! Are those two components Q and U?
       if(all(info%comps == [Q,U])) then
          call build_cov_pol(assembly, info, data, out)
          return
       end if
    end if

    if (all(info%comps == [T])) then
       call build_cov_temp(assembly, info, data, out)
       return
    end if

    allocate(corr(ndi,ndi,-corrlen:corrlen), corr_sq(ndi,ndi,-corrlen:corrlen))
    do i = -corrlen, corrlen
       corr(:,:,i)    = assembly%N_corr_real(i,:,:)
       corr_sq(:,:,i) = assembly%N_corr_F_sq_real(i,:,:)
    end do

    !$OMP PARALLEL PRIVATE(i1,i2,d1,d2,h1,h2,p1,p2,c1,c2,k,N_phase,N_phase_sq,rbuffer,rbuffer_sq,col)
    allocate(rbuffer   (info%ncomp,info%ncomp,out%col_start:out%col_stop))
    allocate(rbuffer_sq(info%ncomp,info%ncomp,out%col_start:out%col_stop))
    !$OMP DO SCHEDULE(guided)
    ! For each sample-diode-group...
    do i1 = 1, nsamp
       if (mod(i1,1000) == 0) write(*,*) i1, nsamp, corrlen, info%cid
       call bench_start(info%bench, bench_cov4)
       do d1 = 1, ndi
          do h1 = 1, ngroups(d1)
             ! ... find its contribution to its row in the covariance matrix ...
             rbuffer    = 0
             rbuffer_sq = 0
             p1 = data%map2mask(assembly%diode_info(d1)%pix(h1,i1))

             ! ... by looping over every correlated sample-diode-group ...
             call bench_start(info%bench, bench_cov1)
             do i2 = max(i1-corrlen,1), min(i1+corrlen, nsamp)
                k = i2-i1
                do d2 = 1, ndi
                   N_phase   = data%phase(d1)%val(:,h1,i1) * corr(d1,d2,k)
                   N_phase_sq= data%phase(d1)%val(:,h1,i1) * corr_sq(d1,d2,k)
                   do h2 = 1, ngroups(d2)
                      ! ... and computing its contribution to each polarization component
                      p2 = data%map2mask(assembly%diode_info(d2)%pix(h2,i2))
                      if(.not. (p2 >= out%col_start .and. p2 <= out%col_stop)) cycle
                      do c1 = 1, info%ncomp
                         do c2 = 1, info%ncomp
                            ! 8 levels of loops! Yummy!
                            rbuffer(c1,c2,p2)    = rbuffer(c1,c2,p2)    + &
                             & N_phase(c1) * data%phase(d2)%val(c2,h2,i2)
                            rbuffer_sq(c1,c2,p2) = rbuffer_sq(c1,c2,p2) + &
                             & N_phase_sq(c1) * data%phase(d2)%val(c2,h2,i2)
                         end do
                      end do
                   end do
                end do
             end do
             call bench_stop(info%bench, bench_cov1)

             ! We now have the contribution to this row, so do an atomic update
             ! (to avoid clobbering)
             call bench_start(info%bench, bench_cov2)
             do col = out%col_start, out%col_stop
                if (rbuffer(1,1,col) == 0.d0 .or. p1 > col .or. p1 <= 0) cycle
                do c2 = 1, info%ncomp
                   do c1 = 1, info%ncomp
                      !$OMP ATOMIC
                      out%cov(p1,c1,col,c2) = out%cov(p1,c1,col,c2) + rbuffer(c1,c2,col)
                   end do
                end do
             end do
             call bench_stop(info%bench, bench_cov2)

            call bench_start(info%bench, bench_cov3)
             do col = out%col_start, out%col_stop
                if (rbuffer_sq(1,1,col) == 0.d0 .or. p1 < col .or. p1 <= 0) cycle
                do c2 = 1, info%ncomp
                   do c1 = 1, info%ncomp
                      !$OMP ATOMIC
                      out%cov(p1+1,c1,col,c2) = out%cov(p1+1,c1,col,c2) + rbuffer_sq(c1,c2,col)
                   end do
                end do
             end do
             call bench_stop(info%bench, bench_cov3)
          end do
       end do
       call bench_stop(info%bench, bench_cov4)
    end do
    !$OMP END DO
    deallocate(rbuffer,rbuffer_sq)
    !$OMP END PARALLEL

!!$    call wall_time(t2)
!!$    write(*,*) 'wall time = ', t2-t1
!!$    write(*,*) 'element   = ', count(out%cov == 0.d0)
!!$    write(*,*) 'sum       = ', sum(abs(out%cov))
!!$    stop

  end subroutine build_cov


  subroutine build_cov_pol(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out

    real(dp)             :: r
    integer(i4b)         :: i1, i2, d1, d2, p1, p2, c1, c2, i, j, k, ndi, nsamp, n
    integer(i4b)         :: omp_get_thread_num
    integer(i4b)         :: corrlen, col
    real(cov_type)       :: N_phase(info%ncomp), N_phase_sq(info%ncomp)
    real(cov_type), allocatable, dimension(:,:,:) :: rbuffer, rbuffer_sq, corr, corr_sq
    real(dp), allocatable, dimension(:,:) :: phase

    if(.not. associated(out%cov)) return

    ndi     = assembly%num_diodes
    nsamp   = assembly%numsamp
    corrlen = assembly%noise_corr_length

    allocate(corr(ndi,ndi,-corrlen:corrlen), corr_sq(ndi,ndi,-corrlen:corrlen))
    do i = -corrlen, corrlen
       corr(:,:,i)    = assembly%N_corr_real(i,:,:)
       corr_sq(:,:,i) = assembly%N_corr_F_sq_real(i,:,:)
    end do

    !$OMP PARALLEL PRIVATE(i1,i2,d1,d2,p1,p2,c1,c2,k,N_phase,N_phase_sq,rbuffer,rbuffer_sq,col)
    allocate(rbuffer   (info%ncomp,info%ncomp,out%col_start:out%col_stop))
    allocate(rbuffer_sq(info%ncomp,info%ncomp,out%col_start:out%col_stop))
    !$OMP DO SCHEDULE(guided)
    ! For each sample-diode-group...
    do i1 = 1, nsamp
       call bench_start(info%bench, bench_cov4)
       rbuffer    = 0.d0
       rbuffer_sq = 0.d0
       p1 = data%map2mask(assembly%diode_info(1)%pix(1,i1))
       
       call bench_start(info%bench, bench_cov1)
       do i2 = max(i1-corrlen,1), min(i1+corrlen, nsamp)
          p2 = data%map2mask(assembly%diode_info(1)%pix(1,i2))
          if(.not. (p2 >= out%col_start .and. p2 <= out%col_stop)) cycle
          k = i2-i1
          do d2 = 1, ndi
             do d1 = 1, ndi
                rbuffer(1,1,p2)    = rbuffer(1,1,p2)    + data%phase(d1)%val(1,1,i1) * corr(d1,d2,k)    * data%phase(d2)%val(1,1,i2)
                rbuffer(2,1,p2)    = rbuffer(2,1,p2)    + data%phase(d1)%val(2,1,i1) * corr(d1,d2,k)    * data%phase(d2)%val(1,1,i2)
                rbuffer(1,2,p2)    = rbuffer(1,2,p2)    + data%phase(d1)%val(1,1,i1) * corr(d1,d2,k)    * data%phase(d2)%val(2,1,i2)
                rbuffer(2,2,p2)    = rbuffer(2,2,p2)    + data%phase(d1)%val(2,1,i1) * corr(d1,d2,k)    * data%phase(d2)%val(2,1,i2)
                rbuffer_sq(1,1,p2) = rbuffer_sq(1,1,p2) + data%phase(d1)%val(1,1,i1) * corr_sq(d1,d2,k) * data%phase(d2)%val(1,1,i2)
                rbuffer_sq(2,1,p2) = rbuffer_sq(2,1,p2) + data%phase(d1)%val(2,1,i1) * corr_sq(d1,d2,k) * data%phase(d2)%val(1,1,i2)
                rbuffer_sq(1,2,p2) = rbuffer_sq(1,2,p2) + data%phase(d1)%val(1,1,i1) * corr_sq(d1,d2,k) * data%phase(d2)%val(2,1,i2)
                rbuffer_sq(2,2,p2) = rbuffer_sq(2,2,p2) + data%phase(d1)%val(2,1,i1) * corr_sq(d1,d2,k) * data%phase(d2)%val(2,1,i2)
             end do
          end do
       end do
       call bench_stop(info%bench, bench_cov1)
          
       ! We now have the contribution to this row, so do an atomic update
       ! (to avoid clobbering)
       call bench_start(info%bench, bench_cov2)
       do col = out%col_start, out%col_stop
          if (rbuffer(1,1,col) == 0.d0 .or. p1 > col .or. p1 <= 0) cycle
          do c2 = 1, info%ncomp
             do c1 = 1, info%ncomp
                !$OMP ATOMIC
                out%cov(p1,c1,col,c2) = out%cov(p1,c1,col,c2) + rbuffer(c1,c2,col)
             end do
          end do
       end do
       call bench_stop(info%bench, bench_cov2)
       
       call bench_start(info%bench, bench_cov3)
       do col = out%col_start, out%col_stop
          if (rbuffer_sq(1,1,col) == 0.d0 .or. p1 < col .or. p1 <= 0) cycle
          do c2 = 1, info%ncomp
             do c1 = 1, info%ncomp
                !$OMP ATOMIC
                out%cov(p1+1,c1,col,c2) = out%cov(p1+1,c1,col,c2) + rbuffer_sq(c1,c2,col)
             end do
          end do
       end do
       call bench_stop(info%bench, bench_cov3)
       
       call bench_stop(info%bench, bench_cov4)
    end do
    !$OMP END DO
    deallocate(rbuffer,rbuffer_sq)
    !$OMP END PARALLEL

  end subroutine build_cov_pol

  subroutine build_cov_temp(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out

    integer(i4b)         :: i1, i2, d1, d2, p1, p2, c1, c2, i, j, k, ndi, nsamp, n, h1, h2
    integer(i4b)         :: omp_get_thread_num
    integer(i4b)         :: corrlen, col, ngroups(assembly%num_diodes)
    real(dp)             :: t1, t2
    real(cov_type)       :: N_phase(info%ncomp), N_phase_sq(info%ncomp)
    real(cov_type), allocatable, dimension(:)     :: rbuffer, rbuffer_sq
    real(cov_type), allocatable, dimension(:,:,:) :: corr, corr_sq

    if(.not. associated(out%cov)) return

    call wall_time(t1)

    ndi     = assembly%num_diodes
    nsamp   = assembly%numsamp
    corrlen = assembly%noise_corr_length

    allocate(corr(ndi,ndi,-corrlen:corrlen), corr_sq(ndi,ndi,-corrlen:corrlen))
    do i = -corrlen, corrlen
       corr(:,:,i)    = assembly%N_corr_real(i,:,:)
       corr_sq(:,:,i) = assembly%N_corr_F_sq_real(i,:,:)
    end do

    !$OMP PARALLEL PRIVATE(i1,i2,d1,d2,p1,p2,c1,c2,h1,h2,k,N_phase,N_phase_sq,rbuffer,rbuffer_sq,col)
    allocate(rbuffer   (out%col_start:out%col_stop))
    allocate(rbuffer_sq(out%col_start:out%col_stop))
    !$OMP DO SCHEDULE(guided)
    ! For each sample-diode-group...
    do i1 = 1, nsamp
       !if (mod(i1,1000) == 0) write(*,*) i1, nsamp, corrlen, info%cid
       call bench_start(info%bench, bench_cov4)
       do h1 = 1, 2
          ! ... find its contribution to its row in the covariance matrix ...
          rbuffer    = 0
          rbuffer_sq = 0
          p1 = data%map2mask(assembly%diode_info(1)%pix(h1,i1))

          ! ... by looping over every correlated sample-diode-group ...
          call bench_start(info%bench, bench_cov1)
          do h2 = 1, 2
             do i2 = max(i1-corrlen,1), min(i1+corrlen, nsamp)
                k = i2-i1
                p2 = data%map2mask(assembly%diode_info(1)%pix(h2,i2))
                if(.not. (p2 >= out%col_start .and. p2 <= out%col_stop)) cycle
                do d2 = 1, ndi
                   do d1 = 1, ndi
                      rbuffer(p2)    = rbuffer(p2)    + data%phase(d1)%val(1,h1,i1) * corr(d1,d2,k)    * data%phase(d2)%val(1,h2,i2)
                      rbuffer_sq(p2) = rbuffer_sq(p2) + data%phase(d1)%val(1,h1,i1) * corr_sq(d1,d2,k) * data%phase(d2)%val(1,h2,i2)
                   end do
                end do
             end do
             call bench_stop(info%bench, bench_cov1)
          end do
             
          ! We now have the contribution to this row, so do an atomic update
          ! (to avoid clobbering)
          call bench_start(info%bench, bench_cov2)
          do col = out%col_start, out%col_stop
             if (rbuffer(col) == 0.d0 .or. p1 > col .or. p1 <= 0) cycle
             !$OMP ATOMIC
             out%cov(p1,1,col,1) = out%cov(p1,1,col,1) + rbuffer(col)
          end do
          call bench_stop(info%bench, bench_cov2)
          
          call bench_start(info%bench, bench_cov3)
          do col = out%col_start, out%col_stop
             if (rbuffer_sq(col) == 0.d0 .or. p1 < col .or. p1 <= 0) cycle
             !$OMP ATOMIC
             out%cov(p1+1,1,col,1) = out%cov(p1+1,1,col,1) + rbuffer_sq(col)
          end do
          call bench_stop(info%bench, bench_cov3)
       end do
       call bench_stop(info%bench, bench_cov4)
    end do
    !$OMP END DO
    deallocate(rbuffer,rbuffer_sq)
    !$OMP END PARALLEL

!!$    call wall_time(t2)
!!$    write(*,*) 'wall time = ', t2-t1
!!$    write(*,*) 'element   = ', count(out%cov == 0.d0)
!!$    write(*,*) 'sum       = ', sum(abs(out%cov))
!!$    stop

  end subroutine build_cov_temp

  subroutine build_dir(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out
    integer(i4b) :: i, d, numpix, num_samples, num_diodes, p
    real(dp)     :: theta, phi, psi, maz, mel, foo, ang, N_corr

    if(.not. associated(out%dir)) return
    ! The diode angle is almost, but not quite useable as it is -
    ! we just need to subtract the diode angle and deck angle.
    ! To do this, we transform a psi=0 from TELE to HOR - the
    ! result angle will be just what we must subtract.
    d = 1; i = 1
    call pix2ang_nest(assembly%nside, assembly%diode_info(d)%pix(1,i), theta, phi)
    call coord_convert(assembly%coordinate_system, phi, theta, 0d0, &
      & COORD_TELE, maz, mel, foo, mod=assembly%diodes(d,1), &
      & diode=assembly%diodes(d,2), mjd=assembly%time(i))
    call coord_convert(COORD_TELE, maz, mel, 0d0, COORD_HOR, &
      & phi, theta, psi, mod=assembly%diodes(d,1), &
      & diode=assembly%diodes(d,2), mjd=assembly%time(i))
    num_diodes  = assembly%num_diodes
    num_samples = assembly%numsamp
    N_corr      = 0
    do i = 1, num_diodes
       N_corr = N_corr + assembly%sigma0(i)**(-2)
    end do
    do i = 1, num_samples
       p = data%map2mask(assembly%diode_info(d)%pix(1,i))
       if(p <= 0) cycle
       ang = assembly%diode_info(d)%point(3,i,1)-psi
       out%dir(p,1) = out%dir(p,1) + cos(2*ang)*N_corr
       out%dir(p,2) = out%dir(p,2) + sin(2*ang)*N_corr
       out%dir(p,3) = out%dir(p,3) + N_corr
    end do
  end subroutine

  subroutine build_ang(assembly, info, data, out)
    implicit none
    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out
    integer(i4b) :: i, d, num_diodes, num_samples, p
    real(dp)     :: N_corr

    if(.not. associated(out%ang)) return
    num_diodes  = assembly%num_diodes
    num_samples = assembly%numsamp

    N_corr = 0
    do i = 1, num_diodes
       N_corr = N_corr + assembly%sigma0(i)**(-2)
    end do
    do d = 1, num_diodes
       do i = 1, num_samples
          p = data%map2mask(assembly%diode_info(d)%pix(1,i))
          if(p <= 0) cycle
          out%ang(p,1) = out%ang(p,1) + sin(2*assembly%diode_info(d)%point(3,i,1))*N_corr
          out%ang(p,2) = out%ang(p,2) + cos(2*assembly%diode_info(d)%point(3,i,1))*N_corr
          out%ang(p,3) = out%ang(p,3) + N_corr
       end do
    end do
  end subroutine

  ! This is for making necessary corrections to the mapmaking component matrices
  ! due to azimuth and spike filtering. The routine name is now slightly inaccurate...
  subroutine process_azimuth_filter(assembly, info, data, out)
    implicit none

    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out
    type(filter_params)  :: filter

    integer(i4b)  :: i, j, k, ntemp, ndi, numsamp, bin, ncomp, npix, d1, az_order, nspikes
    logical(lgt)  :: positive
    real(dp)      :: az_min, az_max, daz, x, t1, t2, time_min, time_max, dt, eff_freq
    real(dp), allocatable, dimension(:,:,:) :: invN_T, invN_T_sq
    real(dp), allocatable, dimension(:)     :: amp, P, az, time, spikes, myeigvals !TMR
    real(dp), allocatable, dimension(:,:)   :: Q, myeigvecs ! TMR

    if ((.not. filter_azimuth) .and. (.not. filter_spikes)) return

    ! ntemp was defined in setup_mapdata.
    ntemp = data%ntemp
    if (ntemp == 0) return

    call wall_time(t1)
    call update_mon(info, "az_init1")
    call get_default_filter(filter)
    ! Hm. It doesn't feel quite right to have the list of spikes in filter_params, but filter_struct is apparently
    ! no longer in use, and putting it in the assembly structure makes for an awful lot of unneccesary copies.

    ndi        = assembly%num_diodes
    numsamp    = assembly%numsamp
    npix       = out%n
    ncomp      = out%ncomp
    nspikes    = 0
    az_order   = 0

    if (filter_azimuth) then
       az_order   = maxval(assembly%az_order)
       if(.not. use_precomp_azfilter) az_order = az_filter_param
       allocate(az(numsamp))
       az = assembly%orig_point(1,:)
       call make_angles_safe(az, 2*pi)
       az_min     = minval(az)
       az_max     = maxval(az)
       daz        = az_max - az_min
    end if

    if (filter_spikes) then
       nspikes = filter%num_spikes
       allocate(time(numsamp))
       allocate(spikes(nspikes))
       time=assembly%time * 86400 ! converting to seconds here!
       spikes=filter%spikefreqs
       time_min = time(1)
       time_max = time(numsamp)
       dt    = time_max - time_min
    end if

!    ntemp      = ndi*(az_order + nspikes*2)
    ! Why do we do this here? data%ntemp was assigned in setup_mapinfo.
 !   data%ntemp = ntemp

    ! Setup of azimuth filter templates
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, x, P)
    !$OMP WORKSHARE
    data%T = 0
    !$OMP END WORKSHARE
    if (filter_azimuth) then
       allocate(P(az_order))
       !$OMP DO
       do i = 1, numsamp
          x = (az(i)-az_min)/daz
          do j = 1, az_order
             P(j) = cos(1*pi*j*x)
          end do
          do j = 1, ndi
             do k = 1, az_order
                data%T(i,j,(j-1)*az_order+k) = P(k)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(P)
    end if
    !$OMP END PARALLEL

    ! Setup of spike filter templates
    if (filter_spikes) then

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, x, Q)
       allocate(Q(nspikes,2))
       !$OMP DO
       do i = 1,numsamp
          x = time(i)-time_min ! time variable in seconds
          do j=1,nspikes
!TMR: testing changing the frequency a little to make periodic function
             eff_freq = floor(spikes(j)*dt) / dt
             Q(j,1) = cos(2*pi*eff_freq*x)
             Q(j,2) = sin(2*pi*eff_freq*x)
          end do
          do j=1,ndi
             do k=1,nspikes
                 data%T(i,j,ndi*az_order + 2*(j-1)*nspikes + 2*k - 1) = Q(k,1)
                 data%T(i,j,ndi*az_order + 2*(j-1)*nspikes + 2*k)     = Q(k,2)
             end do
          end do
       end do
       !$OMP END DO
       deallocate(Q)
       !$OMP END PARALLEL
    end if

!write(*,*) spikes(1), eff_freq

    ! Compute inverse template covariance
    allocate(invN_T(numsamp,ndi,ntemp))
    if(fftconv) then
       call apply_filter_fft_multi(assembly%N_corr_fft,      data%T, invN_T)
    else
       call apply_filter_real_multi(assembly%N_corr_real,    data%T, invN_T)
    end if

!!$open(54,file='template2.txt',action='write')
!!$do i=1,ndi
!!$   do j=1,numsamp
!!$      write(54,*) data%T(j,i,1), invN_T(j,i,1)
!!$   end do
!!$end do
!!$close(54)

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, d1, j)
    !$OMP DO
    do i = 1, ntemp
       if(i<=ndi*az_order) then 
          d1 = (i-1)/az_order + 1 ! az_order may be zero, but then the condition is never true
       else
          d1 = (i-1-ndi*az_order)/(2*nspikes) + 1
       end if
       do j = i, ntemp
          data%Cinv(i,j)    = sum(data%T(:,d1,i)*invN_T(:,d1,j)) 
          data%Cinv(j,i)    = data%Cinv(i,j)
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

!!$write(*,*) 'Cinv: before inversion'
!!$do i=1,ntemp
!!$   write(*,'(4g15.6)') data%Cinv(i,1:4)
!!$end do
!!$

    call invert_matrix(data%Cinv)

!!$!Checking whats up with the correlations
!!$allocate(myeigvals(ntemp))
!!$allocate(myeigvecs(ntemp,ntemp))
!!$
!!$write(*,*) 'Cinv:'
!!$do i=1,ntemp
!!$   write(*,'(4g15.6)') data%Cinv(i,1), data%Cinv(i,2), data%Cinv(i,3), data%Cinv(i,4)
!!$end do
!!$
!!$
!!$call get_eigen_decomposition(info%myid,data%Cinv,myeigvals,myeigvecs)
!!$write(*,*) 'eigvals:'
!!$do i =1,ntemp
!!$   write(*,*) info%myid, myeigvals(i)
!!$end do
!!$write(*,*) 'eigvecs:'
!!$do i=1,ntemp
!!$   write(*,'(4f15.6)') myeigvecs(i,1), myeigvecs(i,2), myeigvecs(i,3), myeigvecs(i,4)
!!$end do

 
    ! Compute P * N_inv * T 
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
    !$OMP WORKSHARE
    data%PNT = 0
    !$OMP END WORKSHARE
    !$OMP DO
    do i = 1, ntemp
       call build_rhs(assembly, info, data, invN_T(:,:,i), data%PNT(i,:,:))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    deallocate(invN_T)

    if (allocated(data%PNT_sq)) then

       ! Compute inverse template covariance for F_sq
       allocate(invN_T_sq(numsamp,ndi,ntemp))
       if(fftconv) then
          call apply_filter_fft_multi(assembly%N_corr_F_sq_fft, data%T, invN_T_sq)
       else
          call apply_filter_real_multi(assembly%N_corr_F_sq_real, data%T, invN_T_sq)
       end if

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, d1, j)
       !$OMP DO
       do i = 1, ntemp
          if(i<=ndi*az_order) then 
             d1 = (i-1)/az_order + 1 
          else
             d1 = (i-1-ndi*az_order)/(2*nspikes) + 1
          end if
          do j = i, ntemp
             data%Cinv_sq(i,j) = sum(data%T(:,d1,i)*invN_T_sq(:,d1,j))          
             data%Cinv_sq(j,i) = data%Cinv_sq(i,j)
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       call invert_matrix(data%Cinv_sq)
       
       ! Compute P * N_inv * T 
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
       !$OMP WORKSHARE
       data%PNT_sq = 0
       !$OMP END WORKSHARE
       !$OMP DO
       do i = 1, ntemp
          call build_rhs(assembly, info, data, invN_T_sq(:,:,i), data%PNT_sq(i,:,:))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       deallocate(invN_T_sq)

       call wall_time(t2)
       !write(*,*) 'CPU time for az filter = ', t2-t1
       call update_mon(info, "az_init2")

    end if

    if(allocated(az)) deallocate(az)
    if(allocated(time)) deallocate(time)
    if(allocated(spikes)) deallocate(spikes)
  end subroutine process_azimuth_filter

  subroutine apply_azimuth_filter_to_rhs(assembly, info, data, N_d, rhs)
    implicit none

    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    real(dp), dimension(:,:), intent(in)    :: N_d
    real(rhs_type), dimension(:,:), intent(inout) :: rhs

    integer(i4b) :: i, d
    real(dp)     :: t1, t2
    real(dp), allocatable, dimension(:)     :: amp

    if ((.not. filter_azimuth) .and. (.not. filter_spikes)) return

    ! Apply azimuth corrections to rhs
    allocate(amp(data%ntemp))
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
    !$OMP DO
    do i = 1, data%ntemp
       amp(i) = sum(data%T(:,:,i)*N_d)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if(any(isnan(rhs))) write(*,*) info%myid, 'found nan in rhs array in apply_az, before applied' !TMR
    amp = matmul(data%Cinv, amp)
    do i = 1, data%ntemp
       rhs = rhs - data%PNT(i,:,:) * amp(i)
    end do
    if(any(isnan(rhs))) write(*,*) info%myid, 'found nan in rhs array in apply_az, after applied' !TMR
    deallocate(amp)

  end subroutine apply_azimuth_filter_to_rhs

  subroutine correct_az_cov(assembly, info, data, out)
    implicit none

    type(quiet_assembly) :: assembly
    type(common_info)    :: info
    type(mapinfo)        :: data
    type(mapdata)        :: out

    integer(i4b) :: p, q, i1, i2, c1, c2
    real(dp), allocatable, dimension(:,:,:) :: PNTC, PNTC_sq

    if(.not. associated(out%cov) .or. ((.not. filter_azimuth) .and. (.not. filter_spikes))) return
    if(data%ntemp == 0) return
    call update_mon(info, "az_cov1")

    ! Compute variance weighted basis vectors
    allocate(PNTC(data%ntemp,out%col_start:out%col_stop,info%ncomp))
    allocate(PNTC_sq(data%ntemp,out%col_start:out%col_stop,info%ncomp))
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(q, c2)
    !$OMP DO
    do q = out%col_start, out%col_stop
       do c2 = 1, info%ncomp
          PNTC(:,q,c2)    = matmul(data%Cinv,    data%PNT(:,q,c2))
          PNTC_sq(:,q,c2) = matmul(data%Cinv_sq, data%PNT_sq(:,q,c2))
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Apply Woodbury correction, -P^T invN^-1 T^T C^-1 T invN P
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(p, q, c1, c2)
    do c2 = 1, info%ncomp
       !$OMP DO
       do q = out%col_start, out%col_stop
          do c1 = 1, info%ncomp
             do p = 1, q
                out%cov(p,c1,q,c2) = out%cov(p,c1,q,c2) - sum(data%PNT(:,p,c1) * PNTC(:,q,c2))
             end do
             do p = q+1, out%n+1
                out%cov(p,c1,q,c2) = out%cov(p,c1,q,c2) - sum(data%PNT_sq(:,p-1,c1) * PNTC_sq(:,q,c2))
             end do
          end do
       end do
       !$OMP END DO
    end do
    !$OMP END PARALLEL
    deallocate(PNTC, PNTC_sq)
    call update_mon(info, "az_cov2")

  end subroutine correct_az_cov

  ! This very slow version of build_rhs_pol takes in an *UNFILTERED*
  ! tod, and builds rhs using assembly%N_corr_real, just like
  ! build_cov does. The purpose of this routine is to test
  ! for inconsistencies between rhs and cov due to inconsistencies
  ! between the real-space and fourier version of the filters.
  subroutine build_rhs_debug(assembly, info, data, from, to)
    implicit none
    type(quiet_assembly),     intent(in) :: assembly
    type(common_info),        intent(in) :: info
    type(mapinfo),            intent(in) :: data
    real(dp),                 intent(in) :: from(:,:)
    real(rhs_type),           intent(out):: to(:,:)
    integer(i4b)              :: d, h, i, j, c, p, ndi, n
    real(dp)                  :: foo(assembly%num_diodes)
    ndi = assembly%num_diodes
    n   = assembly%numsamp
    do i = 1, n
       foo = 0
       do j = max(i-assembly%noise_corr_length,1), min(i+assembly%noise_corr_length, n)
          foo = foo + matmul(assembly%N_corr_real(i-j,:,:), from(j,:))
       end do
       do d = 1, ndi
          do h = 1, size(assembly%diode_info(d)%pix,1)
             p = data%map2mask(assembly%diode_info(d)%pix(h,i))
             if(p <= 0) cycle
             to(p,:) = to(p,:) + data%phase(d)%val(:,h,i)*foo(d)
          end do
       end do
    end do
  end subroutine build_rhs_debug

  !! This is used in conjunctino with build_cov_pol_uvec to build
  !! P e_i.
  !subroutine build_unit_tod_pol(assembly, info, data, pix, comp, tod)
  !  implicit none
  !  type(quiet_assembly),     intent(in) :: assembly
  !  type(common_info),        intent(in) :: info
  !  type(mapinfo),            intent(in) :: data
  !  integer(i4b),             intent(in) :: pix, comp
  !  real(dp),                 intent(out):: tod(:,:)
  !  integer(i4b) :: d
  !  tod = 0
  !  do d = 1, size(tod,2)
  !     where(data%map2mask(assembly%pix(1,:,d)) == pix) tod(:,d) = data%phase(comp,:,d)
  !  end do
  !end subroutine

  !! iN = P'iNtFP = P'iNtF PI = [build_rhs(P e_i)]. iN can
  !! be built by using just build_rhs repeatedly. This is
  !! not optimal, but ensures that there is no inconsistency
  !! between rhs and cov. This prototype implementation only
  !! calculates the F part, not the F^2 part. Thus, the cov
  !! that comes out from this routine cannot be used
  !! for serious work, but is enough to solve for the final
  !! map.
  !subroutine build_cov_pol_uvec(assembly, info, data, out)
  !  implicit none
  !  type(quiet_assembly) :: assembly
  !  type(common_info)    :: info
  !  type(mapinfo)        :: data
  !  type(mapdata)        :: out
  !  integer(i4b)         :: i, j, k, m, n, c, nsamp, ndi
  !  real(dp),       allocatable :: uvec(:,:)
  !  real(rhs_type), allocatable :: covcol(:,:)
  !  nsamp = size(assembly%pix,2)
  !  ndi   = assembly%num_diodes
  !  n     = size(out%cov,3)
  !  allocate(uvec(nsamp,ndi),covcol(n,out%ncomp))
  !  do c = 1, out%ncomp
  !     do i = 1, n
  !        call build_unit_tod_pol(assembly, info, data, i, c, uvec)
  !        call apply_filter_fft_matrix(assembly%N_corr_fft, uvec, uvec)
  !        call build_rhs_pol(assembly, info, data, out, 1, tod=uvec, rhs=covcol)
  !        out%cov(1:n,:,i,c) = covcol
  !     end do
  !  end do
  !  deallocate(uvec,covcol)
  !end subroutine

  subroutine write_cov_indata(assembly, info, data, out, fname)
    implicit none
    type(quiet_assembly)                        :: assembly
    type(common_info)                           :: info
    type(mapinfo)                               :: data
    type(mapdata)                               :: out
    character(len=*)                            :: fname
    type(hdf_file)                              :: file
    integer(i4b)                                :: i
    character(len=512)                          :: pre
    call open_hdf_file(fname, file, "w")
    call create_hdf_group(file, "assembly")
    call write_hdf(file, "assembly/numsamp",           assembly%numsamp)
    call write_hdf(file, "assembly/num_diodes",        assembly%num_diodes)
    call write_hdf(file, "assembly/noise_corr_length", assembly%noise_corr_length)
    call write_hdf(file, "assembly/N_corr_real",       assembly%N_corr_real)
    call write_hdf(file, "assembly/N_corr_F_sq_real",  assembly%N_corr_F_sq_real)
    call write_hdf(file, "assembly/diode_info_n", size(assembly%diode_info))
    do i = 1, size(assembly%diode_info)
       pre = "assembly/diode_info" // trim(itoa(i))
       call create_hdf_group(file, pre)
       call write_hdf(file, trim(pre) // "/pix", assembly%diode_info(i)%pix)
    end do
    call create_hdf_group(file, "info")
    call write_hdf(file, "info/ncomp", info%ncomp)
    call create_hdf_group(file, "data")
    call write_hdf(file, "data/phase_n", size(data%phase))
    do i = 1, size(data%phase)
       pre = "data/phase" // trim(itoa(i))
       call create_hdf_group(file, pre)
       call write_hdf(file, trim(pre) // "/val", data%phase(i)%val)
    end do
    call write_hdf(file, "data/map2mask", data%map2mask)
    call create_hdf_group(file, "out")
    call write_hdf(file, "out/n", int(out%n))
    call write_hdf(file, "out/col_start", int(out%col_start))
    call write_hdf(file, "out/col_stop", int(out%col_stop))
  end subroutine write_cov_indata

  subroutine read_cov_indata(assembly, info, data, out, fname)
    implicit none
    type(quiet_assembly)                        :: assembly
    type(common_info)                           :: info
    type(mapinfo)                               :: data
    type(mapdata)                               :: out
    character(len=*)                            :: fname
    type(hdf_file)                              :: file
    integer(i4b)                                :: i, n, ext(7)
    character(len=512)                          :: pre
    call open_hdf_file(fname, file, "r")
    call read_hdf(file, "assembly/numsamp",           assembly%numsamp)
    call read_hdf(file, "assembly/num_diodes",        assembly%num_diodes)
    call read_hdf(file, "assembly/noise_corr_length", assembly%noise_corr_length)
    call get_size_hdf(file, "assembly/N_corr_real",ext)
    n = (ext(1)-1)/2
    allocate(assembly%N_corr_real(-n:n,ext(2),ext(3)))
    allocate(assembly%N_corr_F_sq_real(-n:n,ext(2),ext(3)))
    call read_hdf(file, "assembly/N_corr_real",       assembly%N_corr_real)
    call read_hdf(file, "assembly/N_corr_F_sq_real",  assembly%N_corr_F_sq_real)
    call read_hdf(file, "assembly/diode_info_n", n)
    allocate(assembly%diode_info(n))
    do i = 1, size(assembly%diode_info)
       pre = "assembly/diode_info" // trim(itoa(i))
       call read_alloc_hdf(file, trim(pre) // "/pix", assembly%diode_info(i)%pix)
    end do
    call read_hdf(file, "info/ncomp", info%ncomp)
    call read_hdf(file, "data/phase_n", n)
    allocate(data%phase(n))
    do i = 1, size(data%phase)
       pre = "data/phase" // trim(itoa(i))
       call read_alloc_hdf(file, trim(pre) // "/val", data%phase(i)%val)
    end do
    call get_size_hdf(file, "data/map2mask", ext)
    allocate(data%map2mask(0:ext(1)-1))
    call read_hdf(file, "data/map2mask", data%map2mask)
    call read_hdf(file, "out/n", n); out%n = n
    call read_hdf(file, "out/col_start", n); out%col_start = n
    call read_hdf(file, "out/col_stop", n); out%col_stop = n
  end subroutine

  subroutine test_mapmaker(out)
    implicit none

    type(mapdata) :: out
    integer(i4b) :: i, j, k, l, n, ntot, unit, ncomp
    real(dp), allocatable, dimension(:,:) :: invcov, V
    real(dp), allocatable, dimension(:)   :: W

    n     = out%n
    ncomp = out%ncomp
    ntot  = n * ncomp
    allocate(invcov(ntot,ntot), V(ntot,ntot), W(ntot))
    do i = 1, ncomp
       do j = 1, ncomp
          invcov((i-1)*n+1:i*n,(j-1)*n+1:j*n) = out%cov(1:n,i,1:n,j)          
       end do
    end do

    do k = 1, ncomp
       do l = 1, ncomp
          do i = 1, n
             do j = i, n
                invcov((l-1)*n+j,(k-1)*n+i) = invcov((k-1)*n+i,(l-1)*n+j)
             end do
          end do
       end do
    end do

    call get_eigen_decomposition(0, invcov, W, V)
    write(*,*) 'min/max eigenval = ', W(1), W(ntot)
    
    unit = getlun()
    open(unit,file='invcov_eigenvals.dat')
    do i = 1, ntot
       write(unit,*) W(i)
    end do
    close(unit)
    stop

  end subroutine test_mapmaker

!  subroutine test_cov
!    ndi     = assembly%num_diodes
!    nsamp   = assembly%numsamp
!    corrlen = assembly%noise_corr_length
!    ngroups(i) = size(data%phase(i)%val,2) 
!    if(all(assembly%diodes(:,1) == assembly%diodes(1,1)) .and. &
!     & size(quiet_horns(assembly%diodes(1,1))%groups) == 1 .and. &
!     & info%ncomp == 2) then
!       if(all(info%comps == [Q,U])) then
!          call build_cov_pol(assembly, info, data, out)
!       corr(:,:,i)    = assembly%N_corr_real(i,:,:)
!       corr_sq(:,:,i) = assembly%N_corr_F_sq_real(i,:,:)
!    allocate(rbuffer   (info%ncomp,info%ncomp,out%col_start:out%col_stop))
!    allocate(rbuffer_sq(info%ncomp,info%ncomp,out%col_start:out%col_stop))
!       call bench_start(info%bench, bench_cov4)
!             p1 = data%map2mask(assembly%diode_info(d1)%pix(h1,i1))
!             call bench_start(info%bench, bench_cov1)
!                   N_phase   = data%phase(d1)%val(:,h1,i1) * corr(d1,d2,k)
!                   N_phase_sq= data%phase(d1)%val(:,h1,i1) * corr_sq(d1,d2,k)
!                      p2 = data%map2mask(assembly%diode_info(d2)%pix(h2,i2))
!                      if(.not. (p2 >= out%col_start .and. p2 <= out%col_stop)) cycle
!                      do c1 = 1, info%ncomp
!                         do c2 = 1, info%ncomp
!                             & N_phase(c1) * data%phase(d2)%val(c2,h2,i2)
!                             & N_phase_sq(c1) * data%phase(d2)%val(c2,h2,i2)
!             call bench_stop(info%bench, bench_cov1)
!             call bench_start(info%bench, bench_cov2)
!             do col = out%col_start, out%col_stop
!  end subroutine

  subroutine initialize_map_debug
    implicit none
    use_precomp_azfilter = .true.
  end subroutine
end module
