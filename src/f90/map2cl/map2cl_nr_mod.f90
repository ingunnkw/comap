module map2cl_nr_mod
  use spline_1D_mod
  use math_tools
  use map2cl_data_mod
  use quiet_nr_mod
  implicit none

  ! Memory cost for module variables = (2*npar+3) * N^2
  logical(lgt),                                  private :: precompute_derivatives, spline_corrfuncs
  logical(lgt),                                  private :: use_cholesky, precompute_pls
  logical(lgt),                                  private :: use_precomputed_basis
  integer(i4b),                                  private :: nspline, max_num_iter, cl_init_mode
  integer(i4b),   allocatable, dimension(:)              :: np
  logical(lgt),                                  private :: change_basis, estimate_errors_at_zero_cl
  real(dp),                                      private :: spline_limit, cl_init
  real(dp),                                      private :: cl_amp_basis
  real(dp),       allocatable, dimension(:),     private :: cos_max_dist
  real(dp),                                      private :: conv_crit, eigenmode_threshold
  real(dp),       allocatable, dimension(:),     private :: ln_det_C
  real(dp),                                      private :: gradient_thresh
  type(planck_rng),                              private :: rng_handle
  type(scalamat), allocatable, dimension(:),     private :: C_inv, data, B, D_minus_C
  type(scalamat), allocatable, dimension(:)              :: N_trans
  type(scalamat), allocatable, dimension(:,:),   private :: dC_dx
  type(powspec)                                          :: x_template
  character(len=512),                            private :: cl_init_file

  integer(i4b), parameter, private :: CONSTANT = 1
  integer(i4b), parameter, private :: FILE     = 2
  integer(i4b), parameter, private :: PCL      = 3

  type pls_struct
     real(dp), allocatable, dimension(:,:,:) :: pls
  end type pls_struct

contains

  subroutine initialize_nr_mod(parfile, info, root, rng_handle_in)
    implicit none

    character(len=*), intent(in)    :: parfile
    type(scinfo),     intent(in)    :: info
    logical(lgt),     intent(in)    :: root
    type(planck_rng), intent(inout) :: rng_handle_in

    character(len=2)   :: set_text
    character(len=4)   :: p_text
    character(len=256) :: basis_def, basis_prefix, basis_file, filename
    integer(i4b)       :: i, j, seed, npar, unit
    logical(lgt)       :: compute_ML_spectrum, exist
    type(scalamat) :: S
    type(powspec)  :: x

    unit = getlun()

    ! Read parameters
    call get_parameter(0, parfile, 'PRECOMPUTE_DERIVATIVES',       par_lgt=precompute_derivatives)
    call get_parameter(0, parfile, 'SPLINE_CORRFUNCS',             par_lgt=spline_corrfuncs)
    call get_parameter(0, parfile, 'SPLINE_LIMIT',                 par_dp=spline_limit)
    call get_parameter(0, parfile, 'CONVERGENCE_CRITERION',        par_dp=conv_crit)
    call get_parameter(0, parfile, 'MAX_NUM_NR_ITERATIONS',        par_int=max_num_iter)
    call get_parameter(0, parfile, 'BASIS',                        par_string=basis_def)
    call get_parameter(0, parfile, 'BASIS_EIGENMODE_THRESHOLD',    par_dp=eigenmode_threshold)
    call get_parameter(0, parfile, 'BASIS_SPECTRUM_AMPLITUDE',     par_dp=cl_amp_basis)
    call get_parameter(0, parfile, 'BASIS_AND_DERIV_PREFIX',       par_string=basis_prefix)
    call get_parameter(0, parfile, 'USE_PRECOMP_BASIS_AND_DERIVS', par_lgt=use_precomputed_basis)
    call get_parameter(0, parfile, 'COMPUTE_ML_SPECTRUM',          par_lgt=compute_ML_spectrum)
    call get_parameter(0, parfile, 'USE_CHOLESKY',                 par_lgt=use_cholesky)
    call get_parameter(0, parfile, 'CL_INIT_MODE',                 par_int=cl_init_mode)
    if (cl_init_mode == CONSTANT) then
       call get_parameter(0, parfile, 'CL_INIT',                   par_dp=cl_init)
    else if (cl_init_mode == FILE) then
       call get_parameter(0, parfile, 'CL_INIT_FILE',              par_string=cl_init_file)
    else if (cl_init_mode == PCL) then

    else
    end if
!    call get_parameter(0, parfile, 'GRADIENT_THRESH',              par_dp=gradient_thresh)
    call get_parameter(0, parfile, 'ESTIMATE_ERRORS_AT_ZERO_CL',   par_lgt=estimate_errors_at_zero_cl)
    spline_limit = cos(spline_limit * pi / 180.d0) ! Do exact calculations for theta < spline_limit
    gradient_thresh = 10.d0

    ! Initialize random number generator
    seed = int(1.d8*rand_uni(rng_handle_in))
    call rand_init(rng_handle, seed)

    ! Current convergence criterion: Delta chi-square must be smaller than specified 
    ! fraction times the number of parametersx
    call allocate_powspec(x, enable_spec, bins, first_bin_print, last_bin_print, num_lowl_junk_bins, num_highl_junk_bins, &
         & num_data_set)
    call copy_powspec(x, x_template)
    npar      = x%npar
    conv_crit = conv_crit * npar 

    if (spline_corrfuncs) then
       allocate(cos_max_dist(num_data_set))
       call get_parameter(0, parfile, 'NUM_SPLINE_KNOTS', par_int=nspline)
       do i = 1, num_data_set
          call get_max_dist(data_sets(i)%vec, cos_max_dist(i))
       end do
    end if

    ! Compute new basis
    allocate(B(num_data_set))
    allocate(np(num_data_set))
    do i = 1, num_data_set
       call int2string(i, set_text)
       basis_file = trim(basis_prefix) // '_set' // set_text// '_basis.unf'
       call setup_basis(info, root, basis_def, basis_file, i, np(i), B(i), change_basis)
    end do

    ! Prepare data
    allocate(data(num_data_set))
    allocate(N_trans(num_data_set))
    allocate(C_inv(num_data_set))
    allocate(D_minus_C(num_data_set))
    allocate(ln_det_C(num_data_set))
    do i = 1, num_data_set
       call sc_alloc(data(i),      np(i),     1, info)
       call sc_alloc(N_trans(i),   np(i), np(i), info)
       call sc_alloc(C_inv(i),     np(i), np(i), info)
       call sc_alloc(D_minus_C(i), np(i), np(i), info)
       call int2string(i, set_text)
       filename = trim(basis_prefix) // '_set' // set_text// '_N_trans.unf'
       inquire(file=trim(filename), exist=exist)
       if (exist) then
          if (info%myid == 0) write(*,*) 'Reading transformed noise covariance from ', trim(filename)
          call read_genmat_sc(unit, info, filename, N_trans(i))
       else
          call change_mat_basis(info, data_sets(i)%N_cov, B(i), N_trans(i))
          if (root) call write_genmat_sc(unit, info, filename, N_trans(i))
       end if
    end do

    if (precompute_derivatives .and. compute_ML_spectrum) then
       allocate(dC_dx(npar,num_data_set))
       do j = 1, num_data_set
          call int2string(j, set_text)
          do i = 1, npar
             call sc_alloc(dC_dx(i,j), np(j), np(j), info)
             call int2string(i, p_text)
             filename = trim(basis_prefix) // '_set' // set_text // '_deriv' // p_text // '.unf'
             
             if (x%ind2set(i) == 0 .or. x%ind2set(i) == j) then
                call get_derivative_matrix(info, root, x%bins(x%ind2bin(i),:), &
                     & x%ind2spec(i), j, filename, dC_dx(i,j))
             else
                dC_dx(i,j)%data = 0.d0
             end if
          end do
       end do
    end if

  end subroutine initialize_nr_mod

 ! Memory cost ~ 0
  subroutine compute_max_lnL_spectrum(info, realization, map_in, x, dx, lnL, chisq, status)
    implicit none

    type(scinfo),                   intent(in)    :: info
    integer(i4b),                   intent(in)    :: realization
    type(map_struct), dimension(:), intent(in)    :: map_in
    type(powspec),                  intent(inout) :: x, dx
    real(dp),                       intent(out)   :: lnL, chisq
    integer(i4b),                   intent(out)   :: status

    logical(lgt) :: converged, pos_def
    integer(i4b) :: iter, i, j, ierr, max_alpha_iter
    real(dp)     :: s, t1, t2, chisq_bf, lnL_bf, alpha
    type(powspec)                         :: x_bf, delta_x
    real(dp), allocatable, dimension(:,:)   :: mat, cl_fix
    real(dp), allocatable, dimension(:,:,:) :: cl_in
    real(dp), allocatable, dimension(:)     :: dCdx
    character(len=4)   :: itertext, realtext
    character(len=128) :: filename
    type(scalamat)     :: map

    allocate(dCdx(x%npar))

    call initialize_map(info, map_in)

    ! Search for maximum likelihood spectrum
    if (cl_init_mode == CONSTANT) then
       do i = 1, x%npar
!          if (x%ind2spec(i) == 1) then
!             x%coeff(i) = 2000.
!          else if (x%ind2spec(i) == 4) then
!             x%coeff(i) = 1.d0
!          else
!             x%coeff(i) = 1.d0
!          end if

          if (x%ind2spec(i) == 1 .or. x%ind2spec(i) == 4 .or. x%ind2spec(i) == 6) then
             x%coeff(i) = cl_init
          else
             x%coeff(i) = 0.d0
          end if
       end do
    else if (cl_init_mode == FILE) then
       allocate(cl_fix(0:lmax,6), cl_in(0:lmax,6,num_data_set))
       cl_fix = 0.d0
       call input_powspec(cl_fix, cl_in, filename=cl_init_file)
       call convert_cls2x(cl_in, x)
       deallocate(cl_fix, cl_in)
    else if (cl_init_mode == PCL) then
       write(*,*) 'PCL INIT MODE NOT IMPLEMENTED YET'
       stop
    else
       write(*,*) 'Unknown initialization mode = ', cl_init_mode
       stop
    end if
    dx%coeff = 0.d0
    call copy_powspec(x, x_bf)
    call copy_powspec(x, delta_x)

!    if (info%myid == 0) then
!       do i = 1, x%npar
!          write(*,*) i, x%coeff(i), x%ind2spec(i), x%ind2bin(i), x%ind2set(i)
!       end do
!    end if
!    call mpi_finalize(ierr)
!    stop

    converged = .false.
    iter      = 0
    s         = 1.d0
    t1        = 0.d0
    t2        = 0.d0
    s         = 0.d0
    alpha     = 1.d0
    max_alpha_iter = 20

    ! Initialize covariance matrix
    call wall_time(t1)
    call compute_covar_mat(info, x=x)
    call wall_time(t2)
    if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time compute_covar_mat = ', t2-t1
    call compute_lnL_and_chisq(info, lnL_bf, chisq_bf)
    if (verbosity > 0 .and. info%myid == 0) then
       write(*,fmt='(a, i4, a,i4,a,f12.2,a,f12.2)') 'Real = ', realization, ', Iter = ', iter, &
            & ' -- chisq(sigma) = ', (chisq-sum(np))/sqrt(2.d0*sum(np)), ' -- lnL = ', lnL
    end if

    ! Do a Newton-Raphson search 
    do while (.not. converged)

       iter = iter+1
       if (iter > max_num_iter) then
          x%coeff  =  1.d30
          dx%coeff =  1.d30
          lnL      = -1.d30
          chisq    =  1.d30
          status   = 1
          exit
       end if

       call wall_time(t1)
       call newton_raphson_step(info, x, dx, delta_x, dCdx_out=dCdx)
       call wall_time(t2)       

       ! Initialize covariance matrix
       pos_def = .false.
       do j = 1, max_alpha_iter
          call wall_time(t1)
          x%coeff = x_bf%coeff + alpha * delta_x%coeff
          call compute_covar_mat(info, x=x, positive_definite=pos_def)
          call wall_time(t2)
          if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time compute_covar_mat = ', t2-t1
          if (pos_def) then
             call compute_lnL_and_chisq(info, lnL, chisq)
          else
             lnL   = lnL_bf   - 100
             chisq = chisq_bf + 100
          end if
          s = lnL-lnL_bf
!          if (s < 0.d0 .or. .not. pos_def) then
          if (.not. pos_def) then
             alpha = 0.5d0*alpha
             if (verbosity > 1 .and. info%myid == 0) then
                write(*,fmt='(a, i4, a,I4,a,f8.2,a,l2,a,f8.2)') 'Real = ', realization, ', Iter = ', iter, &
                     & ' -- alpha = ', alpha, ', pos_def = ', pos_def, ', s = ', s
             end if
          else
             exit
          end if
       end do
       ! If we have been a tight valley, try to move faster again
       if (alpha < 1.d0) alpha = 2.d0*alpha

       if (j >= max_alpha_iter) then
          if (info%myid == 0) then
             write(*,fmt='(a, i4, a,I4,a,f8.2,a,l2)') 'Real = ', realization, ', Iter = ', iter, &
                  & ' did not converge'
          end if
          x%coeff  =  1.d30
          dx%coeff =  1.d30
          lnL      = -1.d30
          chisq    =  1.d30
          status   = 1
          exit
       end if

       ! Check for convergence
       converged = s > 0.d0 .and. s < conv_crit .and. pos_def .and. &
            & (alpha == 1.d0 .or. (alpha >= 0.25d0 .and. any(abs(dCdx) > 10.d0)))
       if (s > 0.d0 .and. pos_def) then
          lnL_bf     = lnL
          chisq_bf   = chisq
          x_bf%coeff = x%coeff
       end if

       if (verbosity > 0 .and. info%myid == 0) then
          if (abs(s) > 1.d10) then
             write(*,fmt='(a, i4, a,I4,a,f5.2,a,f9.1)') 'Real = ', realization, ', Iter = ', iter, &
                  & ' -- dlnL > 1.d10  (', conv_crit, ')  -- time = ', t2-t1             
          else
             write(*,fmt='(a, i4, a,I4,a,f10.2,a,f5.2,a,f9.1)') 'Real = ', realization, ', Iter = ', iter, &
                  & ' -- dlnL = ', real(s,sp), ' (', conv_crit, ')  -- time = ', t2-t1
          end if
          
          write(*,fmt='(a, i4, a,I4,a,f8.2,a,f12.2)') 'Real = ', realization, ', Iter = ', iter, &
               & ' -- chisq(sigma) = ', (chisq-sum(np))/sqrt(2.d0*sum(np)), ' -- lnL = ', lnL
       end if


       if (verbosity > 1 .and. info%myid == 0) then
          ! Output power spectrum
          call int2string(iter, itertext)
          call int2string(realization, realtext)
          filename = 'cls_real' // realtext // '_iter' // itertext // '.dat'
          call output_powspec(filename, .false., x, dx)
       end if
    end do

    ! Do one last call to get the error bars right
    status = 0
    if (estimate_errors_at_zero_cl) then
       ! Set signal power spectrum to zero; useful for null-tests, which should have zero signal
       x%coeff = 0.d0
    else
       ! Compute Fisher matrix at maximum likelihood point
       call copy_powspec(x_bf, x)
    end if
    call compute_covar_mat(info, x=x) 
    call newton_raphson_step(info, x, dx, delta_x, dCdx_out=dCdx)
    call compute_lnL_and_chisq(info, lnL, chisq)

    ! Return maximum likelihood spectrum
    call copy_powspec(x_bf, x)

    ! If the derivatives are non-zero, we compute asymmetric error bars
    ! FIXME: Criterion should be made unitless
    !gradient_thresh=10.0d0
!    if (any(abs(dCdx) > gradient_thresh)) then
!       if (verbosity > 8 .and. info%myid == 0) then
!          write(*,*) 'Real = ', realization, 'Iter = ', iter, ' dCdx fails non-zero test'
!       endif
!       call compute_asymm_error(info, x, dx)
!    end if

    if (verbosity > 0 .and. info%myid == 0) then
       ! Output power spectrum
       call int2string(realization, realtext)
       filename = 'cls_real' // realtext // '.dat'
       call output_powspec(filename, .false., x, dx)
    end if

    deallocate(dCdx)
    call deallocate_powspec(x_bf)
    call deallocate_powspec(delta_x)
    
  end subroutine compute_max_lnL_spectrum

  subroutine initialize_map(info, map_in)
    implicit none

    type(scinfo),                     intent(in)  :: info
    type(map_struct), dimension(:),   intent(in)  :: map_in

    integer(i4b)       :: d
    type(scalamat)     :: map

    ! Prepare data
    do d = 1, num_data_set
       call sc_alloc(map, data_sets(d)%ntot, 1, info)
       call sc_set(map, map_in(d)%data, 0)
       if (change_basis) then
          call sc_matmul(B(d), map, data(d), transa='T')
       else
          call sc_copy(map, data(d))
       end if
       call sc_dealloc(map)
    end do


  end subroutine initialize_map


  ! Memory cost = (2*npar+2) * N^2
  subroutine newton_raphson_step(info, x, sigma, dx, dCdx_out)
    implicit none

    type(scinfo),           intent(in)    :: info
    type(powspec),          intent(in)    :: x
    type(powspec),          intent(inout) :: sigma, dx
    real(dp), dimension(:), intent(out), optional :: dCdx_out

    integer(i4b) :: i, j, k, ierr, d
    real(dp)     :: tr, t1, t2, t3, t4, ddot, my_trace
    logical(lgt) :: root
    real(dp),       allocatable, dimension(:)     :: tr_D_tot, W
    real(dp),       allocatable, dimension(:,:)   :: tr_D, F_tot
    real(dp),       allocatable, dimension(:,:,:) :: F
    type(scalamat)                              :: D_mat, invC_dCdx_invC

    root = info%myid == 0

    allocate(F(x%npar,x%npar,num_data_set))
    allocate(W(x%npar))
    allocate(tr_D(x%npar,num_data_set))
    allocate(F_tot(x%npar,x%npar))
    allocate(tr_D_tot(x%npar))

    call wall_time(t1)
    do d = 1, num_data_set

       call sc_alloc(invC_dCdx_invC, np(d), np(d), info)
       call sc_alloc(D_mat, np(d), np(d), info)

       ! Precompute derivative products
       do i = 1, x%npar
!          if (verbosity > 3 .and. root) write(*,*) 'Computing inv_C * dCdx,  p = ', i
          if (precompute_derivatives) then
             if (use_cholesky) then
                call sc_matmul(C_inv(d), dC_dx(i,d), D_mat, symm='l')
             else
                call sc_matmul(C_inv(d), dC_dx(i,d), D_mat)
             end if
          else
             call get_derivative_matrix(info, root, x%bins(x%ind2bin(i),:), &
                  & x%ind2spec(i), d, 'none', invC_dCdx_invC)
             if (use_cholesky) then
                call sc_matmul(C_inv(d), invC_dCdx_invC, D_mat, symm='l')
             else
                call sc_matmul(C_inv(d), invC_dCdx_invC, D_mat)
             end if
          end if
          if (use_cholesky) then
             call sc_matmul(C_inv(d), D_mat, invC_dCdx_invC, symm='r')
          else
             call sc_matmul(D_mat, C_inv(d), invC_dCdx_invC)
          end if
          tr_D(i,d) = sc_trace_mat_symm_prod(D_minus_C(d), invC_dCdx_invC)       
          do j = i, x%npar
             if (precompute_derivatives) then
                F(i,j,d) = sc_trace_mat_symm_prod(dC_dx(j,d), invC_dCdx_invC)       
             else
                call get_derivative_matrix(info, root, x%bins(x%ind2bin(j),:), &
                     & x%ind2spec(j), d, 'none', D_mat)
                F(i,j,d) = sc_trace_mat_symm_prod(D_mat, invC_dCdx_invC)       
             end if
             F(j,i,d) = F(i,j,d)
          end do
       end do

       call sc_dealloc(invC_dCdx_invC)
       call sc_dealloc(D_mat)

    end do
    call wall_time(t2)
    if (verbosity > 2 .and. root) write(*,*) 'CPU time derivatives = ', t2-t1

    ! Compute Fisher matrix
    F_tot    = 0.d0
    tr_D_tot = 0.d0
    do i = 1, num_data_set
       F_tot    = F_tot + F(:,:,i)
       tr_D_tot = tr_D_tot + tr_D(:,i)
    end do
    if (present(dCdx_out)) dCdx_out = tr_D_tot

!    call wall_time(t1)
    F_tot = 0.5d0 * F_tot
    call invert_matrix(F_tot)
!    call wall_time(t2)
!    if (verbosity > 2 .and. root) write(*,*) 'CPU time Fisher matrix = ', t2-t1

!    call get_eigenvalues(F_tot,W)
!    if (verbosity > 2 .and. root) write(*,*) 'Min eigenval Fisher = ', minval(W)
!    if (verbosity > 2 .and. root) write(*,*) 'dC_dx = ', real(tr_D_tot,sp)

    do i = 1, x%npar
       sigma%coeff(i) = sqrt(F_tot(i,i))
    end do

    ! Find update vector
    dx%coeff = 0.5d0 * matmul(F_tot, tr_D_tot)
    call mpi_bcast(dx%coeff, size(dx%coeff), MPI_DOUBLE_PRECISION, 0, info%comm, ierr)
 !   call wall_time(t2)

!    if (info%myid == 0) then
!       write(*,*) 'dx    = ', real(dx%coeff,sp)
!       write(*,*) 'dCdx  = ', real(dCdx_out,sp)
!       write(*,*) 'sigma = ', real(sigma%coeff,sp)
!    end if

 !   if (verbosity > 2 .and. root) write(*,*) 'CPU time update = ', t2-t1

    deallocate(F, tr_D, F_tot, tr_D_tot, W)

  end subroutine newton_raphson_step

  subroutine compute_lnL_and_chisq(info, lnL, chisq)
    implicit none

    type(scinfo), intent(in)  :: info
    real(dp),     intent(out) :: lnL, chisq

    integer(i4b)   :: ierr, d
    real(dp)       :: chisq_sub
    type(scalamat) :: C_inv_d

    chisq = 0.d0
    lnL   = 0.d0
    do d = 1, num_data_set

       call sc_alloc(C_inv_d, np(d), 1, info)
       if (use_cholesky) then
          call sc_matmul(C_inv(d), data(d), C_inv_d, symm='l')
       else
          call sc_matmul(C_inv(d), data(d), C_inv_d)
       end if
       chisq_sub = sc_dotprod(data(d), C_inv_D)
       chisq     = chisq + chisq_sub
       lnL       = -0.5d0*(chisq_sub + ln_det_C(d))

       call sc_dealloc(C_inv_d)

    end do

    call mpi_bcast(chisq, 1, MPI_DOUBLE_PRECISION, 0, info%comm, ierr)
    call mpi_bcast(lnL,   1, MPI_DOUBLE_PRECISION, 0, info%comm, ierr)

  end subroutine compute_lnL_and_chisq

  subroutine compute_lnL_and_chisq_full(info, cls, lnL, chisq, C_in)
    implicit none

    type(scinfo),                      intent(in)            :: info
    real(dp),     dimension(0:,1:,1:), intent(in)            :: cls
    real(dp),                          intent(out)           :: lnL, chisq
    type(scalamat), dimension(:),      intent(in), optional  :: C_in

    integer(i4b)   :: ierr, status, d, i
    real(dp)       :: ln_det, chisq_sub
    type(scalamat) :: C_inv_d
    type(scalamat) :: C_cov

    chisq = 0.d0
    lnL   = 0.d0
    do d = 1, num_data_set

       call sc_alloc(C_cov, np(d), np(d), info)

       if (present(C_in)) then
          call sc_copy(C_in(d), C_cov)
       else
          ! Compute covariance matrix
          call compute_S_mat(info, d, cls(:,:,d), C_cov, .true.)    
          if (size(C_cov%data) > 0) C_cov%data = C_cov%data + N_trans(d)%data
       end if
       
       ! Factorize matrix
       call sc_cholesky_decompose(C_cov, status=status)
       if (status /= 0) then
          if (info%myid == 0) then
!             write(*,*) 'Component ', d, ', not positive definite'
          end if
          lnL   = -1.d30
          chisq = 1.d30
          call sc_dealloc(C_cov)
          return
       end if

       ! Compute chisquare and determinant
       call sc_alloc(C_inv_d, np(d), 1, info)
       call sc_copy(data(d), C_inv_d)
       call sc_cholesky_solve(C_cov, C_inv_d)
       chisq_sub  = sc_dotprod(data(d), C_inv_D)
       ln_det     = sc_cholesky_logdet(C_cov)
       chisq      = chisq + chisq_sub
       lnL        = lnL -0.5d0*(chisq_sub + ln_det)
!       if (info%myid == 0) then
!          write(*,*) 'chisq_sub = ', chisq_sub
!          write(*,*) 'ln_det    = ', ln_det
!          write(*,*) 'chisq     = ', chisq
!          write(*,*) 'lnL       = ', lnL
!       end if
       call sc_dealloc(C_inv_d)
       call sc_dealloc(C_cov)

    end do

    call mpi_bcast(chisq, 1, MPI_DOUBLE_PRECISION, 0, info%comm, ierr)
    call mpi_bcast(lnL,   1, MPI_DOUBLE_PRECISION, 0, info%comm, ierr)

  end subroutine compute_lnL_and_chisq_full

  ! Memory cost = 2 * N^2
  subroutine compute_covar_mat(info, x, cls, positive_definite)
    implicit none

    type(scinfo),                      intent(in)            :: info
    type(powspec),                     intent(in),  optional :: x
    real(dp),     dimension(0:,1:,1:), intent(in),  optional :: cls
    logical(lgt),                      intent(out), optional :: positive_definite

    real(dp)     :: t1, t2, my_W_max, W_max, my_ln_det_C
    integer(i4b) :: i, j, l, d, status, ierr, my_num_neg_eigenvals, num_neg_eigenvals
    logical(lgt) :: pos_def
    type(scalamat) :: V, C_cov, W
    real(dp),     allocatable, dimension(:,:,:)   :: cls_int

    ! Set up power spectrum
    allocate(cls_int(0:lmax,6, num_data_set))
    if (present(x)) then
       call convert_x2cls(x, cls_int)
    else if (present(cls)) then
       cls_int(0:1,:,:) = 0.d0
       do l = 2, lmax
          cls_int(l,:,:) = cls(l,:,:) * 2.d0*pi / real(l*(l+1),dp)
       end do
    end if
    
    do d = 1, num_data_set

       call sc_alloc(V,     np(d), np(d), info)
       call sc_alloc(W,     np(d),     1, info)
       call sc_alloc(C_cov, np(d), np(d), info)
       
       ! Compute signal covariance matrix
       call wall_time(t1)
       call compute_S_mat(info, d, cls_int(:,:,d), C_cov, .true.)
       call wall_time(t2)
       if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time compute_S_mat = ', t2-t1

       ! Add noise covariance
       if (size(C_cov%data) > 0 .and. .not. nullify_noise_matrix) C_cov%data = C_cov%data + N_trans(d)%data

       ! Set up centered data matrix
       call sc_matmul(data(d), data(d), D_minus_C(d), transb='t')
       D_minus_C(d)%data = D_minus_C(d)%data - C_cov%data

       ! Compute inverse covariance matrix and log-determinant
       if (use_cholesky) then

          ! Factorize matrix
          call wall_time(t1)
          call sc_cholesky_decompose(C_cov, uplo='l', status=status)
          call wall_time(t2)
          if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time cholesky C = ', t2-t1
          
          pos_def = status == 0
          if (present(positive_definite)) positive_definite = pos_def
          if (pos_def) then
             
             ! Compute determinant
             ln_det_C(d) = sc_cholesky_logdet(C_cov)
             
             ! Invert matrix
             call wall_time(t1)
             call sc_cholesky_invert_matrix(C_cov,uplo='l',L=.true.)
             call sc_copy(C_cov, C_inv(d))
             call wall_time(t2)
             if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time inverse C = ', t2-t1
             
          else
             call sc_dealloc(V)
             call sc_dealloc(W)
             call sc_dealloc(C_cov)
             deallocate(cls_int)
             return
          end if

       else

          ! Factorize matrix
          call wall_time(t1)
          call sc_eigenvalue_decompose(C_cov, W, 'l')
          call wall_time(t2)
          if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time eigen_decomp C = ', t2-t1
       
          ! Set non-positive eigenvalues to zero
          my_ln_det_C = 0.d0
          if (size(W%data) > 0) then
             call wall_time(t1)
             my_W_max = maxval(W%data)
             call mpi_allreduce(my_W_max, W_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, info%comm, ierr)
             
             my_num_neg_eigenvals = count(W%data < -1.d-12*W_max)
             call mpi_allreduce(count(W%data <= 0.d0), i, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
             if (info%myid == 0) write(*,*) 'Number of negative eigenvals = ', i, ' of ', np(d)
!             open(58,file='eigenvals.dat')
!             do i = 1, W%rloc
!                write(58,*) i, W%data(i,1)
!             end do
!             close(58)
             do i = 1, W%rloc
                if (W%data(i,1) / W_max > 1.d-20) then
                   my_ln_det_C = my_ln_det_C + log(W%data(i,1))
                   W%data(i,1) = 1.d0 / W%data(i,1)
                else
                   W%data(i,1) = 0.d0
                end if
             end do
          else
             my_W_max             = -1.d30
             my_num_neg_eigenvals = 0
             call mpi_allreduce(my_W_max, W_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, info%comm, ierr)
             call mpi_allreduce(0,            i, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
          end if
          call mpi_allreduce(my_ln_det_C, ln_det_C(d), 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)
       
          if (present(positive_definite)) then
             call mpi_allreduce(my_num_neg_eigenvals, num_neg_eigenvals, 1, MPI_INTEGER, MPI_SUM, info%comm, ierr)
             positive_definite = num_neg_eigenvals == 0
             if (.not. positive_definite) then
                call sc_dealloc(V)
                call sc_dealloc(W)
                call sc_dealloc(C_cov)
                deallocate(cls_int)
                return
             end if
          end if
       
          ! Compute inverse covariance matrix
          call sc_copy(C_cov, V)
          call sc_matmul_diag(V, W)
          call sc_matmul(V, C_cov, C_inv(d), transb='t')
          call wall_time(t2)
          if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time inverse C = ', t2-t1
          
       end if

       call sc_dealloc(V)
       call sc_dealloc(W)
       call sc_dealloc(C_cov)

    end do

    deallocate(cls_int)

  end subroutine compute_covar_mat


  subroutine compute_corrfuncs(z, d, lmin, lmax, cls, C)
    implicit none

    real(dp),                            intent(in)           :: z
    integer(i4b),                        intent(in)           :: d
    integer(i4b),                        intent(in)           :: lmin, lmax
    real(dp),     dimension(0:,1:),      intent(in)           :: cls
    real(dp),     dimension(6),          intent(out)          :: C

    integer(i4b) :: l
    real(dp)     :: Q, R
    real(dp), allocatable, dimension(:,:,:) :: pls

    allocate(pls(0:lmax,-1:1,-1:1))

    call compute_legendre_polynomial(z, pls)
    !call compute_legendre_polynomial_QQ_UU(z, pls)

    C = 0.d0
    if (data_sets(d)%nmaps == 1) then
       do l = lmin, lmax
          C(1) = C(1) + real(2*l+1,dp) * data_sets(d)%beam(l,1)**2 * cls(l,1) * pls(l,0,0)
       end do
    else
       do l = lmin, lmax
          Q = 0.5d0 * (pls(l,1,1) + pls(l,1,-1))
          R = 0.5d0 * (pls(l,1,1) - pls(l,1,-1))
          
          C(1) = C(1) + real(2*l+1,dp) * data_sets(d)%beam(l,1)**2 * cls(l,1) * pls(l,0,0)
          C(2) = C(2) - real(2*l+1,dp) * data_sets(d)%beam(l,1)*data_sets(d)%beam(l,2) * &
               & cls(l,2) * pls(l,0,1)
          C(3) = C(3) - real(2*l+1,dp) * data_sets(d)%beam(l,1) * data_sets(d)%beam(l,3) * &
               & cls(l,3) * pls(l,0,1)
          C(4) = C(4) + real(2*l+1,dp) * (data_sets(d)%beam(l,2)**2 * cls(l,4) * Q + &
               & data_sets(d)%beam(l,3)**2 * cls(l,6) * R)
          C(5) = C(5) + real(2*l+1,dp) *  data_sets(d)%beam(l,2)*data_sets(d)%beam(l,3) * &
               & cls(l,5) * (Q - R)
          C(6) = C(6) + real(2*l+1,dp) * (data_sets(d)%beam(l,2)**2 * cls(l,4) * R + &
               & data_sets(d)%beam(l,3)**2 * cls(l,6) * Q)
       end do
    end if
    C = C / (4.d0*pi)

    deallocate(pls)

  end subroutine compute_corrfuncs



 ! Memory cost = N_tot^2
  subroutine compute_S_mat(info, d, cls, S_out, apply_fg_corr, change_basis_int, skip_EB_corr)
    implicit none

    type(scinfo),                   intent(in)             :: info
    integer(i4b),                   intent(in)             :: d
    real(dp),     dimension(0:,1:), intent(in)             :: cls
    type(scalamat),                 intent(inout)          :: S_out
    logical(lgt),                   intent(in)             :: apply_fg_corr 
    logical(lgt),                   intent(in),   optional :: change_basis_int, skip_EB_corr

    logical(lgt) :: cb
    integer(i4b) :: i, j, k, l, row, col, lmin, lmax, p, q, nc, np, counter, nspec, ntot
    real(dp)     :: z, t1, t2
    real(dp),       allocatable, dimension(:)     :: x
    real(dp),       allocatable, dimension(:,:,:) :: C
    real(dp),       allocatable, dimension(:,:)   :: t2a_1, t2a_2, corr1, corr2, cls_int
    real(dp),                    dimension(6)     :: C_lin
    type(scalamat), allocatable, dimension(:,:)   :: S
    type(scalamat)                                :: S_glob

    nc   = nfield
    np   = data_sets(d)%n_accept
    ntot = data_sets(d)%ntot

    if (all(cls == 0.d0)) then
       S_out%data = 0.d0
       return
    end if

    allocate(cls_int(0:size(cls(:,1))-1,6))
    cls_int = cls

    ! Apply foreground corrections, if any
    if (data_sets(d)%subtract_cl_fg .and. apply_fg_corr) then
       if (present(skip_EB_corr)) then
          if (skip_EB_corr) then
             cls_int(:,1:4) = cls_int(:,1:4) + data_sets(d)%cl_fg(:,1:4)
             cls_int(:,6)   = cls_int(:,6) + data_sets(d)%cl_fg(:,6)
          else
             cls_int = cls_int + data_sets(d)%cl_fg
          end if
       else
          cls_int = cls_int + data_sets(d)%cl_fg
       end if
    end if

    lmin = 0
    do while (all(cls_int(lmin,:) == 0.d0))
       lmin = lmin+1
    end do
    lmax = size(cls_int(:,1))-1
    do while (all(cls_int(lmax,:) == 0.d0))
       lmax = lmax-1
    end do

    call wall_time(t1)
    if (spline_corrfuncs) then
       ! Spline correlation functions
       allocate(C(nspline,2,6), x(nspline))
       do i = 1, nspline
          x(i) = cos_max_dist(d) + (1.d0-cos_max_dist(d)) * real(i-1,dp) / real(nspline-1,dp)
          call compute_corrfuncs(x(i), d, lmin, lmax, cls_int, C(i,1,:))
       end do
       do i = 1, 6
          call spline(x, C(:,1,i), 1.d30, 1.d30, C(:,2,i)) ! QQ
       end do
    end if
    call wall_time(t2)
!    if (info%myid == 0) write(*,*) 'spline CPU = ', t2-t1

!    if (info%myid == 0) then
!       open(58,file='corr.dat')
!       do i = 1, 6
!          do j = 1, nspline
!             write(58,*) x(j), C(j,1,i)
!          end do
!          write(58,*)
!       end do
!       close(58)
!    end if
    

    ! Loop over all pixel pairs
    allocate(S(nc,nc))
    allocate(t2a_1(nc,nc), t2a_2(nc,nc), corr1(nc,nc), corr2(nc,nc))
    do i = 1, nc
       do j = 1, nc
          call sc_alloc(S(i,j), np, np, info)
       end do
    end do
    call wall_time(t1)
    do p = 1, S(1,1)%rloc
       i = S(1,1)%rmap(p)
       do q = 1, S(1,1)%cloc
          j = S(1,1)%cmap(q)
          if (i /= j) then
             z = sum(data_sets(d)%vec(:,i)*data_sets(d)%vec(:,j))
          else
             z = 1.d0 !1.d-1
          end if

          ! Get rotation angles
          call compute_rotation_angle(data_sets(d)%vec(:,i), data_sets(d)%vec(:,j), t2a_1, t2a_2)
!          call compute_rotation_angle_old(data_sets(d)%vec(:,i), data_sets(d)%vec(:,j), t2a)

          ! Get rotationally invariant correlation functions
          if (spline_corrfuncs .and. z <= spline_limit) then
             do k = 1, 6
                C_lin(k) = splint_uniform_grid(x, C(:,1,k), C(:,2,k), z)
             end do
          else
             call compute_corrfuncs(z, d, lmin, lmax, cls_int, C_lin)
             !if (p == q) write(*,*) p, z, d, lmin, lmax, real(C_lin,sp)
          end if
          call convert_lin2mat(C_lin, corr1)

          corr2 = matmul(t2a_1, matmul(corr1, transpose(t2a_2)))

!          write(*,*)
          do row = 1, nc
!             write(*,*) real(corr2(row,:),sp)
             do col = 1, nc
                S(row,col)%data(p,q) = corr2(row,col)
             end do
          end do

       end do
    end do
    call wall_time(t2)
!    if (info%myid == 0) write(*,*) 'build CPU = ', t2-t1

    ! Copy information into a global matrix
    call wall_time(t1)
    call sc_alloc(S_glob, ntot, ntot, info)    
    do i = 1, nc
       do j = 1, nc
          call sc_copy(S(i,j), S_glob, rect=(/(i-1)*np+1,(j-1)*np+1,i*np,j*np/), which=2)
          call sc_dealloc(S(i,j))
       end do
    end do


    if (present(change_basis_int)) then
       cb = change_basis_int
    else
       cb = change_basis
    end if
    if (cb) then
       call change_mat_basis(info, S_glob, B(d), S_out)
    else
       call sc_copy(S_glob, S_out)
    end if
    call wall_time(t2)
!    if (info%myid == 0) write(*,*) 'copy CPU = ', t2-t1

    if (spline_corrfuncs) deallocate(C, x)
    call sc_dealloc(S_glob)
    deallocate(S)
    deallocate(t2a_1, t2a_2, corr1, corr2, cls_int)

  end subroutine compute_S_mat



  subroutine get_derivative_matrix(info, root, bin, spec, d, filename, dC_dx)
    implicit none

    integer(i4b),                       intent(in)    :: spec, d
    integer(i4b),     dimension(2),     intent(in)    :: bin
    type(scinfo),                       intent(in)    :: info
    logical(lgt),                       intent(in)    :: root
    character(len=*),                   intent(in)    :: filename
    type(scalamat),                     intent(inout) :: dC_dx
    
    integer(i4b) :: i, j, b, l, s, unit
    logical(lgt) :: exist
    real(dp)     :: t1, t2
    character(len=4)   :: p_text
    real(dp), allocatable, dimension(:,:)   :: dCl


    unit = getlun()

    allocate(dCl(0:lmax,6))

    call wall_time(t1)
    inquire(file=trim(filename), exist=exist)
    if (exist .and. use_precomputed_basis) then
       if (info%myid == 0) write(*,*) 'Reading derivative matrix from ', trim(filename)
       call read_genmat_sc(unit, info, filename, dC_dx)
    else
          
       ! Set up derivative spectrum
       dCl = 0.d0
       do l = bin(1), bin(2)
          dCl(l,spec) = 2.d0*pi / real(l*(l+1),dp) 
       end do

       ! Compute covariance matrix
       call compute_S_mat(info, d, dCl, dC_dx, .false.)

       ! Write derivatives to disk
       if (root .and. trim(filename) /= 'none') then
          if (info%myid == 0) write(*,*) 'Writing derivative matrix to ', trim(filename)
          call write_genmat_sc(unit, info, filename, dC_dx)
       end if
    end if

    call wall_time(t2)
!    if (verbosity > 2 .and. info%myid == 0) write(*,*) 'CPU time precompute derivatives = ', t2-t1

    deallocate(dCl)

  end subroutine get_derivative_matrix


  subroutine setup_basis(info, root, basis_def, basis_file, d, np, B, change_basis)
    implicit none

    type(scinfo),     intent(in)  :: info
    character(len=*), intent(in)  :: basis_def, basis_file
    logical(lgt),     intent(in)  :: root
    integer(i4b),     intent(in)  :: d
    integer(i4b),     intent(out) :: np
    type(scalamat),   intent(out) :: B
    logical(lgt),     intent(out) :: change_basis

    integer(i4b) :: i, j, l, status, ierr, unit, ntot
    logical(lgt) :: exist
    real(dp), allocatable, dimension(:,:) :: cls, W_glob, S_old
    type(scalamat)                          :: S, V, W

    unit = getlun()
    ntot = data_sets(d)%ntot

    change_basis = trim(basis_def) /= 'pixel'
    if (trim(basis_def) == 'pixel') then
       ! Set up output variables
       np           = data_sets(d)%ntot
    else if (trim(basis_def) == 'S_eigenmodes' .or. trim(basis_def) == 'karhunen_loeve' .or. &
         & trim(basis_def) == 'S/N_eigenmodes') then

       inquire(file=trim(basis_file), exist=exist)
       if (exist .and. use_precomputed_basis) then
          ! Read precomputed basis
          if (info%myid == 0) write(*,*) 'Reading basis vectors from ', trim(basis_file)
          call read_genmat_sc(unit, info, basis_file, B)
          np = B%cglob
       else

          ! Compute signal covariance matrix for fiducal C_l = 1 spectrum
          allocate(cls(0:lmax,6))
          allocate(W_glob(ntot,1))
          call sc_alloc(S, ntot, ntot, info)
          call sc_alloc(V, ntot, ntot, info)
          call sc_alloc(W, ntot,    1, info)
          cls = 0.d0
          do l = 2, lmax
             if (trim(basis_def) == 'S_eigenmodes') then
                do i = 1, 6
                   if (enable_spec(i)) cls(l,i) = 2.d0*pi / real(l*(l+1),dp) 
                end do
             else if (trim(basis_def) == 'karhunen_loeve' .or. trim(basis_def) == 'S/N_eigenmodes') then
                do i = 1, 6
                   if (enable_spec(i)) then
                      if (cl_fid(l,i) == 0.d0 .and. i == 6) then
                         cls(l,i) = cl_fid(l,4) * cl_amp_basis
                      else
                         cls(l,i) = cl_fid(l,i) * cl_amp_basis
                      end if
                   end if
                end do
             end if
          end do
          call compute_S_mat(info, d, cls, S, .false., change_basis_int=.false.)
          deallocate(cls)

          if (trim(basis_def) == 'karhunen_loeve' .or. trim(basis_def) == 'S/N_eigenmodes') then
             ! Pre-whiten with noise covariance
             call sc_matmul(data_sets(d)%inv_L_noise, S, V)
             call sc_matmul(V, data_sets(d)%inv_L_noise, S)
          end if

          ! Decompose into eigenmodes
          call sc_eigenvalue_decompose(S, W, 'l')
          call sc_get(W, W_glob, 0)

          if (info%myid == 0) then
             
             if (verbosity > 0) then
                open(58,file='eigenvals_for_new_basis.dat')
                do l = 1, ntot
                   write(58,*) l, abs(W_glob(l,1)/W_glob(ntot,1)), W_glob(l,1)
                end do
                close(58)
             end if
             
             ! Pick out modes larger than given threshold
             np = count(W_glob/W_glob(ntot,1) > eigenmode_threshold)

          end if
          call mpi_bcast(np, 1, MPI_INTEGER, 0, info%comm, ierr)

          call sc_alloc(B, ntot, np, info)
          call sc_copy(S, B, rect=(/1, ntot-np+1, ntot, ntot/), which=1)

          if (root) then
             ! Write precomputed basis
             if (info%myid == 0) write(*,*) 'Writing basis vectors to ', trim(basis_file)
             call write_genmat_sc(unit, info, basis_file, B)
          end if

          call sc_dealloc(V)
          call sc_dealloc(W)
          call sc_dealloc(S)
          deallocate(W_glob)

       end if

    else
       write(*,*) 'Unknown basis = ', trim(basis_def)
       write(*,*) 'Valid choices are {S_eigenmodes}'
       stop
    end if

  end subroutine setup_basis


  ! A -- ntot x ntot ; pixel space
  ! B -- ntot x np   ; projection operator
  ! C -- np x np     ; reduced space
  subroutine change_mat_basis(info, A, B, C)
    implicit none

    type(scinfo),    intent(in)    :: info
    type(scalamat),  intent(in)    :: A, B
    type(scalamat),  intent(inout) :: C

    integer(i4b) :: ntot, np
    type(scalamat) :: V

    ! C = B^T * A * B
    if (change_basis) then
       ntot = B%rglob 
       np   = B%cglob
       call sc_alloc(V, ntot, np, info)
       call sc_matmul(A, B, V)
       call sc_matmul(B, V, C, transa='T')
       call sc_dealloc(V)
    else
       call sc_copy(A, C)
    end if

  end subroutine change_mat_basis

  subroutine compute_asymm_error(info, x, dx)
    implicit none

    type(scinfo),           intent(in)    :: info
    type(powspec),          intent(in)    :: x
    type(powspec),          intent(inout) :: dx

    integer(i4b) :: d, i, j, k, spec, bin, set, lmin_fit, lmax_fit
    real(dp)     :: upper, lower, lnL_max, lnL, step, chisq, mu, sigma, dCl
    real(dp),       allocatable, dimension(:,:,:) :: cls, cls_fix
    type(scalamat), allocatable, dimension(:)     :: C_fix, C_cov, C_fit
    real(dp),       allocatable, dimension(:)     :: cl, P

    step = 0.1

    allocate(cls_fix(0:lmax,6,num_data_set))
    call convert_x2cls(x, cls_fix)

    do i = 1, x%npar
       spec = x%ind2spec(i)
       set  = x%ind2set(i)
       bin  = x%ind2bin(i)

       ! Only bother with CMB bins
       if (set == 0) then

          call plot_single_slice(info, cls_fix, set, spec, bin, cl, P)
          dCl   = cl(2)-cl(1)
          mu    = sum(P*cl) * dCl
          sigma = sqrt(sum(P * (cl-mu)**2) * dCl)

!          if (info%myid == 0 .and. i ==3) then
!             open(58,file='slice.dat')
!             do k = 1, size(cl)
!                write(58,*) cl(k), P(k)
!             end do
!             close(58)
!          end if

!          write(*,*) 'Sigma for par ', i, ' = ', sigma
          dx%coeff(i) = sigma
          
          deallocate(cl, P)

       end if
    end do

    deallocate(cls_fix)

  end subroutine compute_asymm_error


  subroutine plot_single_slice(info, cls_fix, set, spec, bin, cl, P)
    implicit none

    integer(i4b),                                       intent(in)  :: set, spec, bin
    type(scinfo),                                       intent(in)  :: info
    real(dp),                      dimension(0:,1:,1:), intent(in)  :: cls_fix
    real(dp),         allocatable, dimension(:),        intent(out) :: cl, P

    integer(i4b) :: b, d, i, j, l, r_num, lmax_fit, lmin_fit, pos
    real(dp)     :: r_min, r_max, chisq, dr, r_fid, left, right, dx
    real(dp),                    dimension(1000)  :: x, lnL
    real(dp),       allocatable, dimension(:,:,:) :: cls
    type(scalamat), allocatable, dimension(:)     :: C_fix, C_cov, C_fit

    ! Read parameters
    dx = 0.1d0

    allocate(cls(0:lmax,6,num_data_set))
    allocate(C_cov(num_data_set), C_fix(num_data_set), C_fit(num_data_set))
    do d = 1, num_data_set
       call sc_alloc(C_cov(d), np(d), np(d), info)
       call sc_alloc(C_fix(d), np(d), np(d), info)
       call sc_alloc(C_fit(d), np(d), np(d), info)
    end do

    lmin_fit = bins(bin,1)
    lmax_fit = bins(bin,2)

    ! Set up fixed part of covariance matrix
    cls = 0.d0
    cls(2:lmax,:,:)               = cls_fix(2:lmax,:,:)  
    cls(lmin_fit:lmax_fit,spec,:) = 0.d0
    do d = 1, num_data_set
       call compute_S_mat(info, d, cls(:,:,d), C_fix(d), .true.)    
    end do
             
    ! Set up covariance matrix for usable range
    cls = 0.d0
    do d = 1, num_data_set
       cls(lmin_fit:lmax_fit,spec,d) = cls_fix(lmin_fit:lmax_fit,spec,d) 
       call compute_S_mat(info, d, cls(:,:,d), C_fit(d), .false.)    
    end do

    ! Initialize search at given point
    pos    = 1
    x(pos) = 1.d0
    do d = 1, num_data_set
       C_cov(d)%data = x(pos) * C_fit(d)%data + C_fix(d)%data + N_trans(d)%data
    end do
    call compute_lnL_and_chisq_full(info, cls, lnL(pos), chisq, C_in=C_cov)

    ! Search to the right 
    do while (pos < 1000)
       pos    = pos+1
       x(pos) = minval(x(1:pos-1)) - dx
       do d = 1, num_data_set
          C_cov(d)%data = x(pos) * C_fit(d)%data + C_fix(d)%data + N_trans(d)%data
       end do
       call compute_lnL_and_chisq_full(info, cls, lnL(pos), chisq, C_in=C_cov)
       if (maxval(lnL(1:pos)) - lnL(pos) > 5.d0) exit
    end do

    ! Search to the right 
    do while (pos < 1000)
       pos    = pos+1
       x(pos) = maxval(x(1:pos-1)) + dx
       do d = 1, num_data_set
          C_cov(d)%data = x(pos) * C_fit(d)%data + C_fix(d)%data + N_trans(d)%data
       end do
       call compute_lnL_and_chisq_full(info, cls, lnL(pos), chisq, C_in=C_cov)
       if (maxval(lnL(1:pos)) - lnL(pos) > 5.d0) exit
    end do
    
    ! Convert from relative to absolute C_l
    x(1:pos) = x(1:pos) * cls_fix(lmin_fit,spec,1) * (lmin_fit*(lmin_fit+1)/(2.d0*pi))

    ! Sort values to get them in the right order
    call QuickSort_dp_dist(lnL(1:pos), x(1:pos))

    ! Convert to probability
    lnL(1:pos) = exp(lnL(1:pos) - maxval(lnL(1:pos)))
    lnL(1:pos) = lnL(1:pos) / (sum(lnL(1:pos)) * abs(x(2)-x(1)))
  
    ! Output slice
    allocate(cl(pos),P(pos))
    cl = x(1:pos)
    P  = lnL(1:pos)

    deallocate(cls)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
       call sc_dealloc(C_fix(d))
       call sc_dealloc(C_fit(d))
    end do
    deallocate(C_cov, C_fix, C_fit)
    
  end subroutine
  

end module map2cl_nr_mod
