module map2cl_data_mod
  use healpix_types
  use quiet_utils
  use map2cl_utils
  use fitstools
  use quiet_fileutils
  use rngmod
  use alm_tools
  use math_tools
  use scalautils
  use quiet_postutils
  implicit none

  type map_struct
     real(dp), allocatable, dimension(:,:) :: data
  end type map_struct

  type data_set
     integer(i4b) :: nside, npix, nmaps, ntot, n_accept, num_templates
     logical(lgt) :: subtract_cl_fg
     real(dp),       allocatable, dimension(:,:) :: vec
     real(dp),       allocatable, dimension(:,:) :: maps, beam, cl_fg
     type(scalamat)                              :: N_cov, L_noise, inv_L_noise
     integer(i4b),   pointer,     dimension(:)   :: map2mask, mask2map     
  end type data_set

  integer(i4b) :: lmax, polar, nbin, numreal, nfield
  integer(i4b) :: verbosity, num_data_set
  integer(i4b), private :: numsim
  logical(lgt)          :: nullify_noise_matrix

  real(dp),       allocatable, dimension(:,:) :: cl_fid, cl_fid_BB
  integer(i4b),   pointer,     dimension(:,:) :: bins

  real(dp)                                  :: noise_level_transfunc
  type(planck_rng),   private               :: rng_handle
  logical(lgt),       private               :: add_signal, add_noise
  character(len=2),   private               :: units
  character(len=256), private               :: basis
  logical(lgt), dimension(6)                :: enable_spec
  integer(i4b)                              :: num_lowl_junk_bins, num_highl_junk_bins
  integer(i4b)                              :: first_bin_print, last_bin_print

  type(data_set), allocatable, dimension(:) :: data_sets

contains

  subroutine initialize_data_mod(parfile, rng_handle_in, info)
    implicit none

    character(len=*), intent(in) :: parfile
    type(planck_rng), intent(in) :: rng_handle_in
    type(scinfo),     intent(in) :: info

    integer(i4b)       :: unit, i, j, k, l, nside_in, ordering, seed, polarization, m, lmax_bin
    integer(i8b)       :: n
    real(dp)           :: t1, t2
    character(len=256) :: binfile, filename, cl_BB_file
    logical(lgt)       :: inv, exist1, exist2
    real(dp),           pointer, dimension(:,:)   :: cl_in
    character(len=80),           dimension(1:180) :: header
    character(len=256)                            :: clfile, noiseformat
    character(len=256)                            :: outtext, parname
    type(scalamat)                                :: mat

    unit       = getlun()
    polar      = 1
    rng_handle = rng_handle_in

    ! Read parameters from file
    call get_parameter(0, parfile, 'BINFILE',              par_string=binfile)
    call get_parameter(0, parfile, 'NUM_CMBFILES',         par_int=numreal)
    call get_parameter(0, parfile, 'NOISEFORMAT',          par_string=noiseformat)
    call get_parameter(0, parfile, 'FID_CL_FILE',          par_string=clfile)
    call get_parameter(0, parfile, 'FID_CL_BB_FILE',       par_string=cl_BB_file)
    call get_parameter(0, parfile, 'SEED',                 par_int=seed)
    call get_parameter(0, parfile, 'ADD_SIGNAL_TO_SIMS',   par_lgt=add_signal)
    call get_parameter(0, parfile, 'ADD_NOISE_TO_SIMS',    par_lgt=add_noise)
    call get_parameter(0, parfile, 'VERBOSITY',            par_int=verbosity)
    call get_parameter(0, parfile, 'NUMSIM',               par_int=numsim)
    call get_parameter(0, parfile, 'BASIS',                par_string=basis)
    call get_parameter(0, parfile, 'NUM_LOWL_JUNK_BINS',   par_int=num_lowl_junk_bins)
    call get_parameter(0, parfile, 'NUM_HIGHL_JUNK_BINS',  par_int=num_highl_junk_bins)
    call get_parameter(0, parfile, 'FIRST_BIN_TO_PRINT',   par_int=first_bin_print)
    call get_parameter(0, parfile, 'LAST_BIN_TO_PRINT',    par_int=last_bin_print)
    call get_parameter(0, parfile, 'NULLIFY_NOISE_MATRIX', par_lgt=nullify_noise_matrix)
    call get_parameter(0, parfile, 'UNITS',                par_string=units)
    call get_parameter(0, parfile, 'NOISE_LEVEL_TRANSFUNC',par_dp=noise_level_transfunc)
    call get_parameter(0, parfile, 'LMAX',                 par_int=lmax)

    call get_parameter(0, parfile, 'INCLUDE_TT',          par_lgt=enable_spec(1))
    call get_parameter(0, parfile, 'INCLUDE_TE',          par_lgt=enable_spec(2))
    call get_parameter(0, parfile, 'INCLUDE_TB',          par_lgt=enable_spec(3))
    call get_parameter(0, parfile, 'INCLUDE_EE',          par_lgt=enable_spec(4))
    call get_parameter(0, parfile, 'INCLUDE_EB',          par_lgt=enable_spec(5))
    call get_parameter(0, parfile, 'INCLUDE_BB',          par_lgt=enable_spec(6))

    call get_parameter(0, parfile, 'NUM_DATA_SETS',       par_int=num_data_set)
    
    if (enable_spec(1) .and. (enable_spec(4) .or. enable_spec(6))) then
       nfield = 3
    else if (enable_spec(4) .or. enable_spec(6)) then
       nfield = 2
    else
       nfield = 1
    end if

    ! Set up bin information
    call read_bins(binfile, bins, lmax_bin)
    if (lmax_bin .NE. lmax) then
       if (info%myid == 0) write(*,*) 'WARNING: lmax_binfile= ', lmax_bin, ' is not equal to lmax=',lmax
    end if
    nbin = size(bins(:,1))

    ! Read fiducial power spectrum
    allocate(cl_fid(0:lmax,6))
    allocate(cl_fid_BB(0:lmax,6))
    call read_powspec_map2cl(clfile, cl_fid)
    
!    do l = 2, lmax
!       cl_fid(l,1) = 1000.d0  / (real(l*(l+1),dp)/2.d0/pi)
!       cl_fid(l,2) = 50.d0    / (real(l*(l+1),dp)/2.d0/pi)
!       cl_fid(l,3) = -10.d0   / (real(l*(l+1),dp)/2.d0/pi)
!       cl_fid(l,4) = 30.d0    / (real(l*(l+1),dp)/2.d0/pi)
!       cl_fid(l,5) =  5d0     / (real(l*(l+1),dp)/2.d0/pi)
!       cl_fid(l,6) = 10.d0    / (real(l*(l+1),dp)/2.d0/pi)
!    end do

    call read_powspec_map2cl(cl_BB_file, cl_fid_BB)
    cl_fid_BB(:,1:5) = 0.d0

    ! Set up data sets
    allocate(data_sets(num_data_set))
    do k = 1, num_data_set
       call initialize_data_set(info, parfile, k, data_sets(k))
    end do

  end subroutine initialize_data_mod


  subroutine initialize_data_set(info, parfile, id, d)
    implicit none

    character(len=*), intent(in)    :: parfile
    integer(i4b),     intent(in)    :: id
    type(data_set),   intent(inout) :: d
    type(scinfo),     intent(in)    :: info

    integer(i4b)       :: unit, i, j, l, nside_in, ordering, seed, polarization, m
    integer(i8b)       :: n
    real(dp)           :: t1, t2, alpha
    character(len=256) :: binfile, cmbfile, maskfile, filename, cl_BB_file, pixwinfile
    logical(lgt)       :: inv, exist1, exist2, inc_temp, inc_pol
    real(dp), allocatable, dimension(:,:)   :: map_in, mask_in
    real(dp), pointer,     dimension(:,:)   :: pixwin, cl_in
    real(dp), allocatable, dimension(:)     :: W_noise
    real(dp), allocatable, dimension(:,:)   :: t
    character(len=4)                        :: nside_text, map_text
    character(len=2)                        :: temp_text, set_text
    character(len=80),     dimension(1:180) :: header
    character(len=256)                      :: clfile, beamfile, noisefile, noiseformat, basis_prefix
    character(len=256)                      :: sqrtnoisefile, invsqrtnoisefile
    character(len=256)                      :: outtext,parname
    type(scalamat)                          :: mat

    unit       = getlun()

    call int2string(id, set_text)

    call get_parameter(0, parfile, 'INCLUDE_TEMPERATURE',    par_lgt=inc_temp)
    call get_parameter(0, parfile, 'INCLUDE_POLARIZATION',   par_lgt=inc_pol)
    call get_parameter(0, parfile, 'BASIS_AND_DERIV_PREFIX', par_string=basis_prefix)

    ! Read mask
    parname = 'MASKFILE' // set_text
    call get_parameter(0, parfile, parname, par_string=filename)
    call read_map(mask_in, ordering, filename, nside=d%nside, nmap=d%nmaps)
    if (inc_temp) then
       if (ordering /= 2) then
          call convert_ring2nest(d%nside, mask_in(:,1))
       end if
       call get_map2mask_from_mask(mask_in(:,1), d%map2mask, d%mask2map)
    else
       if (ordering /= 2) then
          call convert_ring2nest(d%nside, mask_in(:,2))
       end if
       call get_map2mask_from_mask(mask_in(:,2), d%map2mask, d%mask2map)
    end if
    deallocate(mask_in)    
    d%npix     = 12*d%nside**2
    d%n_accept = size(d%mask2map)
    d%ntot     = 0
    if (inc_temp) d%ntot = d%ntot +   d%n_accept ! T  map 
    if (inc_pol)  d%ntot = d%ntot + 2*d%n_accept ! QU map 

    ! Compute pixel vectors
    allocate(d%vec(3,d%n_accept))
    do i = 1, d%n_accept
       call pix2vec_nest(d%nside, d%mask2map(i), d%vec(:,i))
    end do

    ! Check for foreground C_l template
    parname = 'SUBTRACT_CL_FG_TEMPLATE' // set_text
    call get_parameter(0, parfile, parname, par_lgt=d%subtract_cl_fg)
    if (d%subtract_cl_fg) then
       parname = 'CL_FG_TEMPLATE' // set_text
       call get_parameter(0, parfile, parname, par_string=clfile)
       allocate(d%cl_fg(0:lmax,6))
       call read_powspec_map2cl(clfile, d%cl_fg)
    end if

    ! Read maps
    allocate(d%maps(d%ntot,numreal))
    do i = 1, numreal
       call int2string(i, map_text)
       parname = 'CMBFILE' // set_text // '_' // map_text
       call get_parameter(0, parfile, parname, par_string=cmbfile)
       call read_map(map_in, ordering, cmbfile, nside=nside_in, nmap=d%nmaps)
       if (units == 'mK') map_in = map_in * 1.d3
       if (nside_in /= d%nside) then
          write(*,*) 'Error: Mask has different Nside than map.'
          stop
       end if
       if (ordering /= 2) then
          do j = 1, d%nmaps
             call convert_ring2nest(d%nside, map_in(:,j))
          end do
       end if

       call map2array_prealloc(map_in, d%map2mask, nfield, d%maps(:,i))
       deallocate(map_in)
    end do

    ! Read noise covariance matrix for Scalapack
    filename = trim(basis_prefix) // '_set' // set_text// '_N_trans.unf'
    inquire(file=trim(filename), exist=exist1)
    if (.not. exist1) then
       if (verbosity > 0 .and. info%myid == 0) write(*,*) 'Reading noise covariance matrix'
       call wall_time(t1)
       parname = 'NOISEFILE' // set_text
       call get_parameter(0, parfile, parname, par_string=filename)
       call read_cov_sc(unit, info, filename, m, ordering, polarization, d%N_cov, inv)
       if (units == 'mK') d%N_cov%data = d%N_cov%data * 1.d6
       if (noise_level_transfunc > 0.d0) d%N_cov%data = d%N_cov%data * noise_level_transfunc**2
       !    if (nullify_noise_matrix) d%N_cov%data = 0.d0
       call wall_time(t2)
       if (verbosity > 0 .and. info%myid == 0) write(*,*) 'Wall time read_cov_sc = ', real(t2-t1,sp)
       if (ordering /= 2 .and. info%myid == 0) then
          write(*,*) 'Error: map2cl expects the noise covariance matrix to be in NESTED ordering'
          stop
       end if

       ! Add template outer-products to noise covariance matrix
       parname = 'NUM_TEMPLATES' // set_text
       call get_parameter(0, parfile, parname, par_int=d%num_templates)
       parname = 'TEMPLATE_AMPLITUDE' // set_text
       call get_parameter(0, parfile, parname, par_dp=alpha)
       do i = 1, d%num_templates
          allocate(t(d%ntot,1))
          if (info%myid == 0) then
             call int2string(i, temp_text)
             parname = 'TEMPLATE' // set_text // '_' // temp_text
             call get_parameter(0, parfile, parname, par_string=filename)
             call read_map(map_in, ordering, filename, nside=nside_in, nmap=d%nmaps)
             if (units == 'mK') map_in = map_in * 1.d3
             if (nside_in /= d%nside) then
                write(*,*) 'Error: Template has different Nside than map.'
                stop
             end if
             if (ordering /= 2) then
                do j = 1, d%nmaps
                   call convert_ring2nest(d%nside, map_in(:,j))
                end do
             end if
             call map2array_prealloc(map_in, d%map2mask, nfield, t(:,1))
             deallocate(map_in)
          end if
          
          call sc_alloc(mat, d%ntot, 1, info)
          call sc_set(mat, t, 0)
          call sc_rank1_update(d%N_cov, alpha, mat, mat)
          call sc_dealloc(mat)
          deallocate(t)
       end do
    end if

    ! Get square root of noise covariance matrix
    filename = trim(basis_prefix) // '_set' // set_text// '_basis.unf'
    inquire(file=trim(filename), exist=exist1)
    if (numsim /= 0 .or. (trim(basis) /= 'pixel' .and. .not. exist1)) then
       call wall_time(t1)
       outtext = 'N_cov'
       call sc_alloc(d%L_noise, m, m, info)
       call sc_alloc(d%inv_L_noise, m, m, info)
       parname = 'SQRT_NOISEFILE' // set_text
       call get_parameter(0, parfile, parname, par_string=sqrtnoisefile)
       parname = 'INV_SQRT_NOISEFILE' // set_text
       call get_parameter(0, parfile, parname, par_string=invsqrtnoisefile)
       inquire(file=trim(sqrtnoisefile), exist=exist1)
       inquire(file=trim(invsqrtnoisefile), exist=exist2)
!       if (nullify_noise_matrix) then
!          d%L_noise%data     = 0.d0
!          d%inv_L_noise%data = 0.d0
!       else if (exist1 .and. exist2 .and. d%num_templates == 0) then
       if (exist1 .and. exist2 .and. d%num_templates == 0) then
          if (numsim > 0 .or. noise_level_transfunc>0) then
             call read_cov_sc(unit, info, sqrtnoisefile, m, ordering, &
                  & polarization, d%L_noise, inv)
             d%L_noise%data = d%L_noise%data * noise_level_transfunc
          end if
          if (units == 'mK') d%L_noise%data = d%L_noise%data * 1.d3
          call read_cov_sc(unit, info, invsqrtnoisefile, m, ordering, polarization, d%inv_L_noise, inv)
          if (units == 'mK') d%inv_L_noise%data = d%inv_L_noise%data * 1.d-3
          if (noise_level_transfunc > 0.d0) d%inv_L_noise%data = d%inv_L_noise%data / noise_level_transfunc
       else
          call sc_alloc(mat, m, m, info)
          call sc_copy(d%N_cov, mat)
          call sqrt_matrix_sc(unit, outtext, mat, d%L_noise, d%inv_L_noise)
          call sc_dealloc(mat)
       end if
       call wall_time(t2)
       if (verbosity > 0 .and. info%myid == 0) write(*,*) 'Wall time L_noise = ', real(t2-t1,sp)
    end if

    ! Read beam (and pixel window)
    parname = 'BEAMFILE' // set_text
    call get_parameter(0, parfile, parname, par_string=filename)
    call get_parameter(0, parfile, 'PIXEL_WINDOW', par_string=pixwinfile)
    allocate(d%beam(0:lmax, d%nmaps))
    call read_beam(filename, d%beam)
    if (trim(pixwinfile) /= 'none') then
       call read_pixwin(d%nside, d%nmaps, pixwin, filename=pixwinfile)
       d%beam = d%beam * pixwin(0:lmax,1:d%nmaps)
       deallocate(pixwin)
    end if

  end subroutine initialize_data_set

  subroutine deallocate_data_set(d)
    implicit none

    type(data_set) :: d

    call sc_dealloc(d%N_cov)
    call sc_dealloc(d%L_noise)
    call sc_dealloc(d%inv_L_noise)

  end subroutine deallocate_data_set

  subroutine cleanup_data_mod
    implicit none

    integer(i4b) :: i

    do i = 1, num_data_set
       call deallocate_data_set(data_sets(i))
    end do

  end subroutine cleanup_data_mod


  subroutine get_simulation(info, maps)
    implicit none

    type(scinfo),                    intent(in)            :: info
    type(map_struct), dimension(1:), intent(inout)         :: maps

    integer(i4b) :: i, l, m, ierr, d
    real(dp)     :: cl
    real(dp),          allocatable, dimension(:,:)   :: s, noise
    complex(dpc),      allocatable, dimension(:,:,:) :: alms, alms2
    character(len=80),              dimension(1:180) :: header
    type(scalamat)                                   :: eta, N

    if (info%myid == 0 .and. add_signal) then
       allocate(alms(3,0:lmax,0:lmax))
       call create_alms(cl_fid, rng_handle, alms)
    end if

    do d = 1, num_data_set
       if (info%myid == 0 .and. add_signal) then
          allocate(alms2(data_sets(d)%nmaps, 0:lmax,0:lmax))
          allocate(s(0:data_sets(d)%npix-1,data_sets(d)%nmaps))
          do i = 1, data_sets(d)%nmaps
             do l = 0, lmax
                alms2(i,l,:) = alms(i,l,:) * data_sets(d)%beam(l,i)
             end do
          end do
          if (data_sets(d)%nmaps == 1) then
             call alm2map(data_sets(d)%nside, lmax, lmax, alms2, s(:,1))
          else
             call alm2map(data_sets(d)%nside, lmax, lmax, alms2, s)
          end if
          do i = 1, data_sets(d)%nmaps
             call convert_ring2nest(data_sets(d)%nside, s(:,i))
          end do
          call map2array_prealloc(s, data_sets(d)%map2mask, nfield, maps(d)%data(:,1))
          deallocate(alms2)
          deallocate(s)
       else
          maps(d)%data = 0.d0
       end if

       if (add_noise) call add_noise_to_data(info, d, maps(d)%data(:,1))
!       if (add_noise) then
!          allocate(noise(data_sets(d)%ntot,1))
!          call sc_alloc(eta, data_sets(d)%ntot, 1, info)
!          call sc_alloc(N, data_sets(d)%ntot, 1, info)
!          if (eta%cloc > 0) then
!             do i = 1, eta%rloc
!                eta%data(i,1) = rand_gauss(rng_handle)
!             end do
!          end if
!          call sc_matmul(data_sets(d)%L_noise, eta, N)
!          call sc_get(N, noise, 0)
!          if (info%myid == 0) then
!             maps(d)%data = maps(d)%data + noise
!          end if
!          call sc_dealloc(eta)
!          call sc_dealloc(N)
!          deallocate(noise)
!       end if

    end do

    if (allocated(alms)) deallocate(alms)

  end subroutine get_simulation

  subroutine add_noise_to_data(info, d, map)
    implicit none

    type(scinfo),            intent(in)    :: info
    integer(i4b),            intent(in)    :: d
    real(dp), dimension(1:), intent(inout) :: map

    integer(i4b) :: i, m
    real(dp), allocatable, dimension(:,:) :: noise
    type(scalamat)                        :: eta, N

    m = size(map)

    allocate(noise(m,1))
    call sc_alloc(eta, m, 1, info)
    call sc_alloc(N,   m, 1, info)
    if (eta%cloc > 0) then
       do i = 1, eta%rloc
          eta%data(i,1) = rand_gauss(rng_handle)
       end do
    end if
    call sc_matmul(data_sets(d)%L_noise, eta, N)
    call sc_get(N, noise, 0)
    if (info%myid == 0) then
       map = map + noise(:,1)
    end if
    call sc_dealloc(eta)
    call sc_dealloc(N)
    deallocate(noise)

  end subroutine add_noise_to_data

  subroutine create_alms(cls, rng_handle, alms)
    implicit none

    real(dp),         dimension(0:,1:),    intent(in)    :: cls
    type(planck_rng),                      intent(inout) :: rng_handle
    complex(dpc),     dimension(1:,0:,0:), intent(out)   :: alms

    integer(i4b) :: i, l, m
    real(dp)     :: cl
    real(dp), dimension(3,3) :: S_cov, L_cov
    real(dp), dimension(3)   :: eta1, eta2
    
    S_cov = 0.d0
    do l = 2, lmax
       call convert_cls2covmat(cls(l,:), S_cov)
       if (all(S_cov == 0.d0)) then
          alms(:,l,:) = 0.d0
          cycle
       else
          call cholesky_decompose_with_mask(S_cov, L_cov)
       end if
       do m = 0, l
          do i = 1, 3
             eta1(i) = rand_gauss(rng_handle)
             eta2(i) = rand_gauss(rng_handle)
          end do
          eta1 = matmul(L_cov, eta1)
          eta2 = matmul(L_cov, eta2)
          if (m == 0) then
             alms(:,l,m) = cmplx(eta1, 0.d0)
          else
             alms(:,l,m) = cmplx(eta1, eta2) / sqrt(2.d0)
          end if
       end do
    end do

  end subroutine create_alms

  subroutine convert_cls2covmat(cls, S_cov)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: cls
    real(dp), dimension(1:,1:), intent(out) :: S_cov

    integer(i4b) :: i, j, k, n

    n = size(S_cov(1,:))

    k = 1
    do i = 1, n
       do j = i, n
          S_cov(i,j) = cls(k)
          S_cov(j,i) = cls(k)
          k          = k+1
       end do
    end do

  end subroutine convert_cls2covmat

!!$  subroutine test_noise_matrix(L, N)
!!$    implicit none
!!$
!!$    real(dp), dimension(:,:), intent(in) :: L, N
!!$
!!$    integer(i4b) :: i, j, numsim, m
!!$    real(dp), allocatable, dimension(:,:) :: N_sim, eta
!!$
!!$    numsim = 1000
!!$    m      = size(L,1)
!!$    allocate(eta(m,1))
!!$    allocate(N_sim(m,m))
!!$
!!$    N_sim = 0.d0
!!$    do i = 1, numsim
!!$       write(*,*) i, numsim
!!$       do j = 1, m
!!$          eta(j,1) = rand_gauss(rng_handle)
!!$       end do
!!$       eta(:,1) = matmul(L_noise, eta(:,1))
!!$!       write(*,*) eta(:,1)
!!$!       stop
!!$       N_sim    = N_sim + matmul(eta, transpose(eta))
!!$    end do
!!$    N_sim = N_sim / real(numsim,dp)
!!$
!!$    write(*,*) N_sim(1,1), N(1,1)
!!$    write(*,*) 'Maximum difference = ', maxval(abs(N_sim-N) / abs(N_sim+N))
!!$    write(*,*) 'Minimum difference = ', minval(abs(N_sim-N) / abs(N_sim+N))
!!$
!!$    open(58,file='test.dat')
!!$    do i = 1, m
!!$       j = i
!!$       write(58,*) abs(N_sim(i,j)+N(i,j)), abs(N_sim(i,j)-N(i,j)) / abs(N_sim(i,j)+N(i,j))
!!$    end do
!!$    close(58)
!!$
!!$    deallocate(N_sim)
!!$    deallocate(eta)
!!$    stop
!!$
!!$  end subroutine test_noise_matrix

end module map2cl_data_mod
