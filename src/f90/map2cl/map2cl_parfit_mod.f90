module map2cl_parfit_mod
  use map2cl_nr_mod
  use map2cl_data_mod
  implicit none

contains


  subroutine fit_qn_model(info, parfile, cls_fix, map, iter, q_out, sigma_q, n_out, sigma_n)
    implicit none

    type(scinfo),                          intent(in)  :: info
    character(len=*),                      intent(in)  :: parfile
    real(dp),         dimension(0:,1:,1:), intent(in)  :: cls_fix
    type(map_struct), dimension(:),        intent(in)  :: map
    integer(i4b),                          intent(in)  :: iter
    real(dp),                              intent(out) :: q_out, sigma_q, n_out, sigma_n

    integer(i4b) :: d, i, j, l, n_num, q_num, lmax_fit, lmin_fit
    real(dp)     :: q_min, q_max, n_min, n_max, l_pivot, chisq, dq, dn
    character(len=4)   :: iter_text
    character(len=128) :: filename
    real(dp),       allocatable, dimension(:)     :: q, n, ls, P_q, P_n
    real(dp),       allocatable, dimension(:,:)   :: lnL
    real(dp),       allocatable, dimension(:,:,:) :: cls
    type(scalamat), allocatable, dimension(:)     :: C_fix, C_cov, C_fit

    ! Set up data
    call initialize_map(info, map)    

    ! Read parameters and set up grid
    call get_parameter(0, parfile, 'NUM_Q_POINTS', par_int=q_num)
    call get_parameter(0, parfile, 'NUM_N_POINTS', par_int=n_num)
    call get_parameter(0, parfile, 'N_MIN',        par_dp=n_min)
    call get_parameter(0, parfile, 'N_MAX',        par_dp=n_max)
    call get_parameter(0, parfile, 'Q_MIN',        par_dp=q_min)
    call get_parameter(0, parfile, 'Q_MAX',        par_dp=q_max)
    call get_parameter(0, parfile, 'L_PIVOT',      par_dp=l_pivot)
    call get_parameter(0, parfile, 'LMIN_FOR_FIT', par_int=lmin_fit)
    call get_parameter(0, parfile, 'LMAX_FOR_FIT', par_int=lmax_fit)
    if (q_num <= 1) then
       dq = 1.d0
    else
       dq = (q_max-q_min) / real(q_num-1,dp)
    end if
    if (n_num == 1) then
       dn = 1.d0
    else
       dn = (n_max-n_min) / real(n_num-1,dp)
    end if
    call int2string(iter, iter_text)
    
    allocate(q(q_num), n(n_num), lnL(q_num, n_num), cls(0:lmax,6,num_data_set), ls(0:lmax))
    allocate(P_q(q_num), P_n(n_num))

    do i = 1, q_num; q(i) = q_min + dq * real(i-1,dp); end do
    do i = 1, n_num; n(i) = n_min + dn * real(i-1,dp); end do
    do l = 0, lmax; ls(l) = real(l,dp); end do

    allocate(C_cov(num_data_set), C_fix(num_data_set), C_fit(num_data_set))
    do i = 1, num_data_set
       call sc_alloc(C_cov(i), np(i), np(i), info)
       call sc_alloc(C_fix(i), np(i), np(i), info)
       call sc_alloc(C_fit(i), np(i), np(i), info)
    end do

    ! Set up fixed part of covariance matrix
    cls = 0.d0
    cls(2:lmin_fit-1,:,:)    = cls_fix(2:lmin_fit-1,:,:)
    cls(lmax_fit+1:lmax,:,:) = cls_fix(lmax_fit+1:lmax,:,:)
    do d = 1, num_data_set
       call compute_S_mat(info, d, cls(:,:,d), C_fix(d), .true.)    
    end do

    ! Loop through grid points and compute log-likelihood
    cls = 0.d0
    do j = 1, n_num
       ! Set up covariance matrix
       do d = 1, num_data_set
          cls(lmin_fit:lmax_fit,4,d) = &
               & cl_fid(lmin_fit:lmax_fit,4) * (ls(lmin_fit:lmax_fit) / l_pivot)**n(j) 
          call compute_S_mat(info, d, cls(:,:,d), C_fit(d), .false.)    
       end do
       do i = 1, q_num
          do d = 1, num_data_set
             C_cov(d)%data = q(i) * C_fit(d)%data + C_fix(d)%data + N_trans(d)%data
          end do
          call compute_lnL_and_chisq_full(info, cls, lnL(i,j), chisq, C_in=C_cov)
          if (lnL(i,j) == -1.d30) then
             if (info%myid == 0) write(*,fmt='(a,i5,a,f5.2,a,f5.2,a)') 'Real = ', iter, &
                  & ', q = ', q(i), ', n = ', n(j), ' not positive definite'
          else
             if (info%myid == 0) write(*,fmt='(a,i5,a,f5.2,a,f5.2,a,f10.2,a,f10.2)') 'Real = ', iter, &
                  & ', q = ', q(i), ', n = ', n(j), ', lnL = ', lnL(i,j), ', chisq = ', chisq
          end if
       end do
    end do
    
    if (info%myid == 0) then

       ! Output "chi-square"
       filename = 'qn_2d_real' // iter_text // '.dat'
       call output_2d_prob(filename, q, n, lnL)
       
       ! Convert to probability
       lnL = exp(lnL - maxval(lnL))
       lnL = lnL / (sum(lnL) * dq*dn)
       
       ! Generate marginals and posterior mean and standard deviations
       do i = 1, q_num
          P_q(i) = sum(lnL(i,:))
       end do
       P_q     = P_q / sum(P_q * dq)
       q_out   = sum(q * P_q) * dq
       sigma_q = sqrt(sum((q-q_out)**2 * P_q) * dq)

       filename = 'q_1d_real' // iter_text // '.dat'
       call output_1d_prob(filename, q, P_q)
       
       do i = 1, n_num
          P_n(i) = sum(lnL(:,i))
       end do
       P_n     = P_n / sum(P_n * dn)
       n_out   = sum(n * P_n) * dn
       sigma_n = sqrt(sum((n-n_out)**2 * P_n) * dn)
       
       filename = 'n_1d_real' // iter_text // '.dat'
       call output_1d_prob(filename, n, P_n)
    end if

    deallocate(q, n, lnL, cls, ls, P_q, P_n)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
       call sc_dealloc(C_fix(d))
       call sc_dealloc(C_fit(d))
    end do
    deallocate(C_cov, C_fix, C_fit)
    
  end subroutine


  subroutine fit_r(info, parfile, cls_fix, map, iter, r_out, sigma_r)
    implicit none

    type(scinfo),                          intent(in)  :: info
    character(len=*),                      intent(in)  :: parfile
    real(dp),         dimension(0:,1:,1:), intent(in)  :: cls_fix
    type(map_struct), dimension(:),        intent(in)  :: map
    integer(i4b),                          intent(in)  :: iter
    real(dp),                              intent(out) :: r_out, sigma_r

    integer(i4b) :: d, i, j, l, r_num, lmax_fit, lmin_fit
    real(dp)     :: r_min, r_max, chisq, dr, r_fid
    character(len=4)   :: iter_text
    character(len=128) :: filename
    real(dp),       allocatable, dimension(:)     :: r, lnL
    real(dp),       allocatable, dimension(:,:,:) :: cls
    type(scalamat), allocatable, dimension(:)     :: C_fix, C_cov, C_fit

    ! Set up data
    call initialize_map(info, map)    

    ! Read parameters and set up grid
    call get_parameter(0, parfile, 'NUM_R_POINTS', par_int=r_num)
    call get_parameter(0, parfile, 'R_MIN',        par_dp=r_min)
    call get_parameter(0, parfile, 'R_MAX',        par_dp=r_max)
    call get_parameter(0, parfile, 'R_FIDUCIAL',   par_dp=r_fid)
    call get_parameter(0, parfile, 'LMIN_FOR_FIT', par_int=lmin_fit)
    call get_parameter(0, parfile, 'LMAX_FOR_FIT', par_int=lmax_fit)
    if (r_num <= 1) then
       dr = 1.d0
    else
       dr = (r_max-r_min) / real(r_num-1,dp)
    end if
    call int2string(iter, iter_text)
    
    allocate(r(r_num), lnL(r_num), cls(0:lmax,6,num_data_set))
    do i = 1, r_num; r(i) = r_min + dr * real(i-1,dp); end do

    allocate(C_cov(num_data_set), C_fix(num_data_set), C_fit(num_data_set))
    do d = 1, num_data_set
       call sc_alloc(C_cov(d), np(d), np(d), info)
       call sc_alloc(C_fix(d), np(d), np(d), info)
       call sc_alloc(C_fit(d), np(d), np(d), info)
    end do

!    if (info%myid == 0) then
!       do l = 
!    end if

    ! Set up fixed part of covariance matrix
    cls = 0.d0
    cls(2:lmax,4:6,:)           = cls_fix(2:lmax,4:6,:)    ! Fix EE, EB and BB spectra outside rane
!    cls(2:lmax,6,:)            = cls_fix(2:lmax,6,:)    ! Fix BB spectrum outside range
    cls(lmin_fit:lmax_fit,5:6,:) = 0.d0
    do d = 1, num_data_set
       call compute_S_mat(info, d, cls(:,:,d), C_fix(d), .true.)    
    end do

    ! Set up covariance matrix for usable range
    cls = 0.d0
    do d = 1, num_data_set
       cls(lmin_fit:lmax_fit,6,d) = cl_fid_BB(lmin_fit:lmax_fit,6) 
       call compute_S_mat(info, d, cls(:,:,d), C_fit(d), .false.)    
       C_fit(d)%data = C_fit(d)%data / r_fid       ! Normalize to r = 1
    end do

    ! Loop through grid points and compute log-likelihood
    do i = 1, r_num
       do d = 1, num_data_set
          C_cov(d)%data = r(i) * C_fit(d)%data + C_fix(d)%data + N_trans(d)%data
       end do
       call compute_lnL_and_chisq_full(info, cls, lnL(i), chisq, C_in=C_cov)
       if (lnL(i) == -1.d30) then
          if (info%myid == 0) write(*,fmt='(a,i5,a,f5.2,a)') 'Real = ', iter, &
               & ', r = ', r(i), '   not positive definite'
       else
          if (info%myid == 0) write(*,fmt='(a,i5,a,f7.2,a,f10.2,a,f10.2)') 'Real = ', iter, &
               & ', r = ', r(i), ', lnL = ', lnL(i), ', chisq = ', chisq
       end if
    end do
    
    if (info%myid == 0) then

       ! Convert to probability
       lnL = exp(lnL - maxval(lnL))
       lnL = lnL / (sum(lnL) * dr)
       
       r_out   = sum(r * lnL) * dr
       sigma_r = sqrt(sum((r-r_out)**2 * lnL) * dr)

       filename = 'r_1d_real' // iter_text // '.dat'
       call output_1d_prob(filename, r, lnL)
       
    end if

    deallocate(r, lnL, cls)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
       call sc_dealloc(C_fix(d))
       call sc_dealloc(C_fit(d))
    end do
    deallocate(C_cov, C_fix, C_fit)
    
  end subroutine


  subroutine fit_r2(info, root, parfile, cls_fix, map, iter, r_out, sigma_r)
    implicit none

    type(scinfo),                          intent(in)  :: info
    logical(lgt),                          intent(in)  :: root
    character(len=*),                      intent(in)  :: parfile
    real(dp),         dimension(0:,1:,1:), intent(in)  :: cls_fix
    type(map_struct), dimension(:),        intent(in)  :: map
    integer(i4b),                          intent(in)  :: iter
    real(dp),                              intent(out) :: r_out, sigma_r

    integer(i4b) :: d, i, j, k, l, r_num, lmax_fit, lmin_fit, n_EE, n_BB, unit
    real(dp)     :: r_min, r_max, chisq, dr, r_fid, dEE, dBB, marg_EE_min, marg_EE_max, marg_BB_min, marg_BB_max
    character(len=4)   :: iter_text
    character(len=128) :: filename
    logical(lgt)       :: exist
    real(dp),       allocatable, dimension(:)     :: r, lnL_marg, EEs, BBs
    real(dp),       allocatable, dimension(:,:)   :: lnL_marg2
    real(dp),       allocatable, dimension(:,:,:) :: cls, lnL
    type(scalamat), allocatable, dimension(:)     :: C_fix, C_cov, C_fit
    type(scalamat), allocatable, dimension(:,:)   :: C_marg

    ! Set up data
    call initialize_map(info, map)    

    ! Read parameters and set up grid
    call get_parameter(0, parfile, 'NUM_R_POINTS', par_int=r_num)
    call get_parameter(0, parfile, 'R_NUM_EE_MARG',   par_int=n_EE)
    call get_parameter(0, parfile, 'R_NUM_BB_MARG',   par_int=n_BB)
    call get_parameter(0, parfile, 'R_MARG_EE_MIN', par_dp=marg_EE_min)
    call get_parameter(0, parfile, 'R_MARG_EE_MAX', par_dp=marg_EE_max)
    call get_parameter(0, parfile, 'R_MARG_BB_MIN', par_dp=marg_BB_min)
    call get_parameter(0, parfile, 'R_MARG_BB_MAX', par_dp=marg_BB_max)
    call get_parameter(0, parfile, 'R_MIN',        par_dp=r_min)
    call get_parameter(0, parfile, 'R_MAX',        par_dp=r_max)
    call get_parameter(0, parfile, 'R_FIDUCIAL',   par_dp=r_fid)
    call get_parameter(0, parfile, 'LMIN_FOR_FIT', par_int=lmin_fit)
    call get_parameter(0, parfile, 'LMAX_FOR_FIT', par_int=lmax_fit)
    unit = getlun()
    if (n_EE <= 1) then
       dEE = 0.d0
    else
       dEE = (marg_EE_max-marg_EE_min) / (n_EE-1)
    end if
    if (n_BB <= 1) then
       dBB = 0.d0
    else
       dBB = (marg_BB_max-marg_BB_min) / (n_BB-1)
    end if
    if (r_num <= 1) then
       dr = 1.d0
    else
       dr = (r_max-r_min) / real(r_num-1,dp)
    end if
    call int2string(iter, iter_text)
    
    allocate(r(r_num), lnL(r_num, n_EE, n_BB), cls(0:lmax,6,num_data_set))
    do i = 1, r_num; r(i) = r_min + dr * real(i-1,dp); end do

    allocate(C_cov(num_data_set), C_fix(num_data_set), C_fit(num_data_set), C_marg(num_data_set,2))
    do d = 1, num_data_set
       call sc_alloc(C_cov(d), np(d), np(d), info)
       call sc_alloc(C_fix(d), np(d), np(d), info)
       call sc_alloc(C_fit(d), np(d), np(d), info)
       call sc_alloc(C_marg(d,1), np(d), np(d), info)
       call sc_alloc(C_marg(d,2), np(d), np(d), info)
    end do

    ! Set up fixed part of covariance matrix
    if (info%myid == 0) write(*,*) 'Computing fixed covariance matrix'
    filename = 'precomp_S_fix.unf'
    inquire(file=trim(filename), exist=exist)
    if (exist) then
       if (info%myid == 0) write(*,*) 'Reading transformed noise covariance from ', trim(filename)
       do d = 1, num_data_set
          call read_genmat_sc(unit, info, filename, C_fix(d))
       end do
    else
       do d = 1, num_data_set
          cls = 0.d0
          cls(lmin_fit:lmax,4,d)   = cls_fix(lmin_fit:lmax,4,d) ! Fix EE and BB spectra at high l's
          cls(lmin_fit:lmax,5,d)   = cls_fix(lmin_fit:lmax,5,d) ! Fix EE and BB spectra at high l's
          cls(lmax_fit+1:lmax,6,d) = cls_fix(lmax_fit+1:lmax,6,d)
          call compute_S_mat(info, d, cls(:,:,d), C_fix(d), .true.)    
          if (info%myid == root) write(*,*) 'Writing matrix = ', trim(filename)
          if (root) call write_genmat_sc(unit, info, filename, C_fix(d))
       end do
    end if

    ! Set up EE marginalization part of covariance matrix
    if (info%myid == 0) write(*,*) 'Computing EE low-l covariance matrix'
    filename = 'precomp_S_EE_lowl.unf'
    inquire(file=trim(filename), exist=exist)
    if (exist) then
       if (info%myid == 0) write(*,*) 'Reading transformed noise covariance from ', trim(filename)
       do d = 1, num_data_set
          call read_genmat_sc(unit, info, filename, C_marg(d,1))
       end do
    else
       cls = 0.d0
       do l = 2, lmin_fit-1
          cls(l,4,:) = 2.d0*pi / real(l*(l+1),dp)    ! Set up unit EE matrix at low l's
       end do
       do d = 1, num_data_set
          call compute_S_mat(info, d, cls(:,:,d), C_marg(d,1), .false.)    
          if (info%myid == root) write(*,*) 'Writing matrix = ', trim(filename)
          if (root) call write_genmat_sc(unit, info, filename, C_marg(d,1))
       end do
    end if

    ! Set up BB marginalization part of covariance matrix
    if (info%myid == 0) write(*,*) 'Computing BB low-l covariance matrix'
    filename = 'precomp_S_BB_lowl.unf'
    inquire(file=trim(filename), exist=exist)
    if (exist) then
       if (info%myid == 0) write(*,*) 'Reading transformed noise covariance from ', trim(filename)
       do d = 1, num_data_set
          call read_genmat_sc(unit, info, filename, C_marg(d,2))
       end do
    else
       cls = 0.d0
       do l = 2, lmin_fit-1
          cls(l,6,:) = 2.d0*pi / real(l*(l+1),dp)    ! Set up unit BB matrix at low l's
       end do
       do d = 1, num_data_set
          call compute_S_mat(info, d, cls(:,:,d), C_marg(d,2), .false.)    
          if (info%myid == root) write(*,*) 'Writing matrix = ', trim(filename)
          if (root) call write_genmat_sc(unit, info, filename, C_marg(d,2))
       end do
    end if

    ! Set up covariance matrix for usable range
    if (info%myid == 0) write(*,*) 'Computing active covariance matrix'
    filename = 'precomp_S_fit.unf'
    inquire(file=trim(filename), exist=exist)
    if (exist) then
       if (info%myid == 0) write(*,*) 'Reading transformed noise covariance from ', trim(filename)
       do d = 1, num_data_set
          call read_genmat_sc(unit, info, filename, C_fit(d))
       end do
    else
       cls = 0.d0
       do d = 1, num_data_set
          cls(lmin_fit:lmax_fit,6,d) = cl_fid_BB(lmin_fit:lmax_fit,6) 
          call compute_S_mat(info, d, cls(:,:,d), C_fit(d), .false.)    
          C_fit(d)%data = C_fit(d)%data / r_fid       ! Normalize to r = 1
          if (info%myid == root) write(*,*) 'Writing matrix = ', trim(filename)
          if (root) call write_genmat_sc(unit, info, filename, C_fit(d))
       end do
    end if

    ! Loop through grid points and compute log-likelihood
    lnL = -1.d30
main: do i = r_num, 1, -1
       do k = n_BB, 1, -1
          do j = n_EE, 1, -1
             do d = 1, num_data_set
                C_cov(d)%data = r(i) * C_fit(d)%data + C_fix(d)%data + N_trans(d)%data + &
                     & (marg_EE_min + (j-1)*dEE) * C_marg(d,1)%data + &
                     & (marg_BB_min + (k-1)*dBB) * C_marg(d,2)%data 
             end do
             call compute_lnL_and_chisq_full(info, cls, lnL(i,j,k), chisq, C_in=C_cov)
             if (lnL(i,j,k) == -1.d30) then
                if (info%myid == 0) write(*,fmt='(a,i5,a,f8.2,a,f8.2,a,f8.2,a,f10.2,a,f10.2)') 'Real = ', iter, &
                     & ', EE = ', (marg_EE_min + (j-1)*dEE), ', BB = ', (marg_BB_min + (k-1)*dBB), ', r = ', r(i), ' not positive definite'
                if (k == n_BB .and. j == n_EE) exit main
                exit
             else
                if (info%myid == 0) write(*,fmt='(a,i5,a,f8.2,a,f8.2,a,f8.2,a,f10.2,a,f10.2)') 'Real = ', iter, &
                     & ', EE = ', (marg_EE_min + (j-1)*dEE), ', BB = ', (marg_BB_min + (k-1)*dBB), ', r = ', r(i), ', lnL = ', lnL(i,j,k), ', chisq = ', chisq
             end if
          end do
       end do
       !if (all(lnL(i,:,:) == -1.d30)) exit
    end do main
    
    if (info%myid == 0) then

       ! Convert to probability
       lnL = exp(lnL - maxval(lnL))

       allocate(lnL_marg(r_num))
       do i = 1, r_num
          lnL_marg(i) = sum(lnL(i,:,:))
       end do

       lnL_marg = lnL_marg / (sum(lnL_marg) * dr)
       
       r_out   = sum(r * lnL_marg) * dr
       sigma_r = sqrt(sum((r-r_out)**2 * lnL_marg) * dr)

       filename = 'r_1d_real' // iter_text // '.dat'
       call output_1d_prob(filename, r, lnL_marg)
       deallocate(lnL_marg)

       allocate(lnL_marg2(n_EE, n_BB), EEs(n_EE), BBs(n_BB))
       do i = 1, n_EE
          EEs(i) = marg_EE_min + (i-1)*dEE
       end do
       do i = 1, n_BB
          BBs(i) = marg_BB_min + (i-1)*dBB
       end do
       lnL_marg2 = 0.d0
       do i = 1, n_EE
          do j = 1, n_BB
             lnL_marg2(i,j) = sum(lnL(:,i,j))
          end do
       end do
       where (lnL_marg2 == 0.d0) lnL_marg2 = 1.d-10
       lnL_marg2 = lnL_marg2 / maxval(lnL_marg2)
       lnL_marg2 = -2.d0*log(lnL_marg2)

       filename = 'EE_BB_2d_real' // iter_text // '.dat'
       call output_2d_prob(filename, EEs, BBs, lnL_marg2)

       ! Output raw data file
       open(unit,file='r_EE_BB.dat', recl=1024)
       do k = 1, n_BB
          do j = 1, n_EE
             do i = 1, r_num
                write(unit,fmt='(4f16.5)') BBs(k), EEs(j), r(i), lnL(i,j,k)
             end do
          end do
       end do
       close(unit)

       deallocate(lnL_marg2, EEs, BBs)
    end if

    deallocate(r, lnL, cls)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
       call sc_dealloc(C_fix(d))
       call sc_dealloc(C_fit(d))
       call sc_dealloc(C_marg(d,1))
       call sc_dealloc(C_marg(d,2))
    end do
    deallocate(C_cov, C_fix, C_fit, C_marg)
    
  end subroutine


  subroutine fit_tab_amplitude(info, parfile, map, realization)

    type(scinfo),                          intent(in)  :: info
    character(len=*),                      intent(in)  :: parfile
    type(map_struct), dimension(:),        intent(in)  :: map
    integer(i4b),                          intent(in)  :: realization

    integer(i4b) :: d, i, j, l, r_num, unit, unit2, n, lmin_tab, lmax_tab
    real(dp)     :: r_min, r_max, chisq, dr, r_fid, a, C(6)
    character(len=4)   :: real_text
    character(len=256) :: directory, filename, fid_cl_tab_filename
    real(dp),       allocatable, dimension(:)     :: amp, lnL
    real(dp),       allocatable, dimension(:,:,:) :: cls
    type(scalamat), allocatable, dimension(:)     :: C_cov
    real(dp)            :: TT, TE, EE, BB

    ! Set up data
    call initialize_map(info, map)    

    ! Read parameters and set up grid
    call get_parameter(0, parfile, 'CL_DIRECTORY', par_string=directory)
    call get_parameter(0, parfile, 'LMIN_TAB', par_int=lmin_tab)
    call get_parameter(0, parfile, 'LMAX_TAB', par_int=lmax_tab)
    call get_parameter(0, parfile, 'FID_CL_TAB_FILE', par_string=fid_cl_tab_filename)

    if (lmin_tab < 2) then
       if (info%myid == 0) write(*,*) 'WARNING: lmin_tab should be 2 or greater, lmin_tab= ', lmin_tab
    end if

    unit = getlun()
    filename = trim(directory) // '/info.txt'
    open(unit, file=trim(filename))
    read(unit,*) n
    
    allocate(amp(n), lnL(n), cls(0:lmax,6,num_data_set))

    allocate(C_cov(num_data_set))
    do d = 1, num_data_set
       call sc_alloc(C_cov(d), np(d), np(d), info)
    end do

    ! Read in cls values below lmin_tab and above lmax_tab
    cls = 0.d0
    unit2 = getlun()
    open(unit2,file=trim(fid_cl_tab_filename))
    if (lmin_tab > 2) then
       do l = 2, lmin_tab-1
          read(unit2,*) i, TT, EE, BB, TE
          cls(l,1,:) = TT / (real(l*(l+1),dp) / (2.d0*pi))
          cls(l,2,:) = TE / (real(l*(l+1),dp) / (2.d0*pi))
          cls(l,4,:) = EE / (real(l*(l+1),dp) / (2.d0*pi))
          cls(l,6,:) = BB / (real(l*(l+1),dp) / (2.d0*pi))
       end do
    end if
    do l = lmin_tab, lmax_tab
       read(unit2,*) i, TT, EE, BB, TE
    end do
    do l = lmax_tab+1, lmax
       read(unit2,*) i, TT, EE, BB, TE
       cls(l,1,:) = TT / (real(l*(l+1),dp) / (2.d0*pi))
       cls(l,2,:) = TE / (real(l*(l+1),dp) / (2.d0*pi))
       cls(l,4,:) = EE / (real(l*(l+1),dp) / (2.d0*pi))
       cls(l,6,:) = BB / (real(l*(l+1),dp) / (2.d0*pi))
    end do
    close(unit2)

    ! Loop through grid points and compute log-likelihood
    do i = 1, n
       read(unit,*) amp(i), filename
       filename = trim(directory) // '/' // trim(filename)
       call input_powspec_simple(filename, cls(:,:,1), lmin_tab, lmax_tab)
       do d = 2, num_data_set
          cls(:,:,d) = cls(:,:,1)
       end do
       call compute_lnL_and_chisq_full(info, cls, lnL(i), chisq)
       if (lnL(i) == -1.d30) then
          if (info%myid == 0) write(*,fmt='(a,i5,a,f5.2,a)') 'Real = ', realization, &
               & ', amp = ', amp(i), '   not positive definite'
       else
          if (info%myid == 0) write(*,fmt='(a,i5,a,f7.2,a,f10.2,a,f10.2)') 'Real = ', realization, &
               & ', amp = ', amp(i), ', lnL = ', lnL(i), ', chisq = ', chisq
       end if
    end do
    close(unit)

    if (info%myid == 0) then
       ! Convert to probability
       lnL = exp(lnL - maxval(lnL))
       filename = 'tabamp_1d_real' // real_text // '.dat'
       call output_1d_prob(filename, amp, lnL)
    end if

    deallocate(amp, lnL, cls)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
    end do
    deallocate(C_cov)

  end subroutine fit_tab_amplitude


  subroutine plot_likelihood_slices(info, parfile, cls_fix, map, iter)
    implicit none

    type(scinfo),                          intent(in)  :: info
    character(len=*),                      intent(in)  :: parfile
    real(dp),         dimension(0:,1:,1:), intent(in)  :: cls_fix
    type(map_struct), dimension(:),        intent(in)  :: map
    integer(i4b),                          intent(in)  :: iter

    integer(i4b) :: b, d, i, j, k, l, r_num, lmax_fit, lmin_fit, pos
    real(dp)     :: r_min, r_max, chisq, dr, r_fid, left, right, dx
    character(len=4)   :: iter_text
    character(len=2)   :: spec_text, bin_text
    character(len=128) :: filename
    real(dp),                    dimension(1000)  :: x, lnL
    real(dp),       allocatable, dimension(:,:,:) :: cls
    type(scalamat), allocatable, dimension(:)     :: C_fix, C_cov, C_fit

    ! Set up data
    call initialize_map(info, map)    

    ! Read parameters
    call int2string(iter, iter_text)
    call get_parameter(0, parfile, 'SLICE_STEP_SIZE',   par_dp=dx)

    allocate(cls(0:lmax,6,num_data_set))
    allocate(C_cov(num_data_set), C_fix(num_data_set), C_fit(num_data_set))
    do d = 1, num_data_set
       call sc_alloc(C_cov(d), np(d), np(d), info)
       call sc_alloc(C_fix(d), np(d), np(d), info)
       call sc_alloc(C_fit(d), np(d), np(d), info)
    end do

!    if (info%myid == 0) then
!       do l = 0, lmax
!          write(*,*) l, real(cls_fix(l,:,1),sp)
!       end do
!    end if
!    stop

    do k = 1, 6
       if (enable_spec(k)) then
          do b = num_lowl_junk_bins+1, nbin-num_highl_junk_bins

             lmin_fit = bins(b,1)
             lmax_fit = bins(b,2)

             ! Set up fixed part of covariance matrix
             cls = 0.d0
             cls(2:lmax,:,:)            = cls_fix(2:lmax,:,:)  
             cls(lmin_fit:lmax_fit,k,:) = 0.d0
             do d = 1, num_data_set
                call compute_S_mat(info, d, cls(:,:,d), C_fix(d), .true.)    
             end do
             
             ! Set up covariance matrix for usable range
             cls = 0.d0
             do d = 1, num_data_set
                cls(lmin_fit:lmax_fit,k,d) = cls_fix(lmin_fit:lmax_fit,k,d) 
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
                if (lnL(pos) == -1.d30) then
                   if (info%myid == 0) write(*,fmt='(a,i5,a,f5.2,a)') 'Real = ', iter, &
                        & ', x = ', x(pos), '   not positive definite'
                else
                   if (info%myid == 0) write(*,fmt='(a,i5,a,f10.2,a,f12.2,a,f12.2)') 'Real = ', iter, &
                        & ', x = ', x(pos), ', lnL = ', lnL(pos), ', chisq = ', chisq
                end if
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
                if (lnL(pos) == -1.d30) then
                   if (info%myid == 0) write(*,fmt='(a,i5,a,f5.2,a)') 'Real = ', iter, &
                        & ', x = ', x(pos), '   not positive definite'
                else                   
                   if (info%myid == 0) write(*,fmt='(a,i5,a,f10.2,a,f12.2,a,f12.2)') 'Real = ', iter, &
                        & ', x = ', x(pos), ', lnL = ', lnL(pos), ', chisq = ', chisq
                end if
                if (maxval(lnL(1:pos)) - lnL(pos) > 5.d0) exit
             end do

             if (info%myid == 0) then

                ! Convert from relative to absolute C_l
                x(1:pos) = x(1:pos) * cls_fix(lmin_fit,k,1) * (lmin_fit*(lmin_fit+1)/(2.d0*pi))

                ! Sort values to get them in the right order
                call QuickSort_dp_dist(lnL(1:pos), x(1:pos))

                ! Convert to probability
                lnL(1:pos) = exp(lnL(1:pos) - maxval(lnL(1:pos)))
                lnL(1:pos) = lnL(1:pos) / (sum(lnL(1:pos)) * abs(x(2)-x(1)))
                
                call int2string(b, bin_text)          
                call int2string(k, spec_text)          
                filename = 'slice_1d_real' // iter_text // '_s' // spec_text // '_b' // bin_text // '.dat'
                call output_1d_prob(filename, x(1:pos), lnL(1:pos))
                
             end if

          end do
       end if
    end do

    deallocate(cls)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
       call sc_dealloc(C_fix(d))
       call sc_dealloc(C_fit(d))
    end do
    deallocate(C_cov, C_fix, C_fit)
    
  end subroutine


  subroutine compute_chisq_vs_null(info, parfile, cls_fix, map, iter)
    implicit none

    type(scinfo),                          intent(in)  :: info
    character(len=*),                      intent(in)  :: parfile
    real(dp),         dimension(0:,1:,1:), intent(in)  :: cls_fix
    type(map_struct), dimension(:),        intent(in)  :: map
    integer(i4b),                          intent(in)  :: iter

    integer(i4b) :: b, d, i, j, k, l, r_num, lmax_fit, lmin_fit, pos, unit, b1, b2, l1, l2
    real(dp)     :: r_min, r_max, chisq, dr, r_fid, left, right, dx, lnL0, lnL, nu
    logical(lgt) :: posdef
    character(len=4)   :: iter_text
    character(len=2)   :: spec_text, bin_text
    character(len=128) :: filename
    real(dp),       allocatable, dimension(:,:,:) :: cls
    type(scalamat), allocatable, dimension(:)     :: C_cov


    ! Set up data
    call initialize_map(info, map)    

    ! Read parameters
    call int2string(iter, iter_text)
    nu       = count(enable_spec)
    unit     = getlun()
    filename = 'chisq_null_real' // iter_text // '.dat'

    allocate(cls(0:lmax,6,num_data_set))
    allocate(C_cov(num_data_set))
    do d = 1, num_data_set
       call sc_alloc(C_cov(d), np(d), np(d), info)
    end do

    ! Compute lnL for default spectrum
    l1 = bins(num_lowl_junk_bins+1,1)
    l2 = bins(nbin-num_highl_junk_bins,2)

    cls = 0.d0
!    cls(2:lmax,:,:)            = cls_fix(2:lmax,:,:)  
    cls(l1:l2,:,:)            = cls_fix(l1:l2,:,:)  
!    if (info%myid == 0) then
!       do l = 0, lmax
!          write(*,*) l, real(cls(l,4:6,1) * l*(l+1)/2./pi, sp)
!       end do
!    end if


!    cls(:,:,:) = 0.d0
!    do d = 1, num_data_set
!       call compute_S_mat(info, d, cls(:,:,d), C_cov(d), .true.)    
!    end do
!    call compute_lnL_and_chisq_full(info, cls, lnL0, chisq, C_in=C_cov)
    call compute_lnL_and_chisq_full(info, cls, lnL0, chisq)
    if (info%myid == 0) write(*,*) 'lnL0 = ', lnL0

!    if (.false.) then
    if (lnL0 == -1.d30) then
       l1 = 2
       l2 = bins(nbin-num_highl_junk_bins,2)

       cls = 0.d0
       cls(l1:l2,:,:)            = cls_fix(l1:l2,:,:)  
       call compute_lnL_and_chisq_full(info, cls, lnL0, chisq)
       if (info%myid == 0) write(*,*) 'lnL0 = ', lnL0

    end if

    ! Compute lnL for each bin set to zero
    if (info%myid == 0) then
       open(unit, file=trim(filename), recl=1024)
       write(unit,*) '#  l_bin    l_1   l_2       delta lnL     chi-square     nu     PTE'
    end if
!    do b = 1,0 
    do b = num_lowl_junk_bins+1, nbin-num_highl_junk_bins

       posdef = .false.
       i      = 0
       do while (.not. posdef)

          b1       = max(b-i, num_lowl_junk_bins+1)
          b2       = min(b+i, nbin-num_highl_junk_bins)
          lmin_fit = bins(b1,1)
          lmax_fit = bins(b2,2)

          ! Set up covariance matrix for usable range
          cls = 0.d0
!          cls(2:lmax,:,:)            = cls_fix(2:lmax,:,:)
          cls(l1:l2,:,:)            = cls_fix(l1:l2,:,:)  
          cls(lmin_fit:lmax_fit,:,:) = 0.d0
!       do d = 1, num_data_set
!          call compute_S_mat(info, d, cls(:,:,d), C_cov(d), .false.)    
!       end do
!       call compute_lnL_and_chisq_full(info, cls, lnL, chisq, C_in=C_cov)
          call compute_lnL_and_chisq_full(info, cls, lnL, chisq)
          if ((lnL == -1.d30 .or. lnL < lnL0) .and. i < 2) then
             i = i+1
             cycle
          else if (lnL == -1.d30 .and. i == 2) then
             ! Return zero degrees of freedom
             lnL = lnL0
             b2  = 0
             b1  = 1 
          end if
          posdef = .true.

          if (info%myid == 0) then
             ! Output change in likelihood, "chi-square" and PTE
             write(*,fmt='(f10.3,i6,i6,f12.4,f12.4,i8,f12.4,f12.4)') 0.5*(lmin_fit+lmax_fit), lmin_fit, lmax_fit, &
                  & lnL0-lnL, 2*(lnL0-lnL), nint(nu*(b2-b1+1)), 0.d0, lnL
             write(unit,fmt='(f10.3,i6,i6,f12.4,f12.4,i8,f12.4)') 0.5*(lmin_fit+lmax_fit), lmin_fit, lmax_fit, &
                  & lnL0-lnL, 2*(lnL0-lnL), nint(nu*(b2-b1+1)), 0.d0
          end if
       end do
    end do

    ! Set up covariance matrix for full range
    b1       = num_lowl_junk_bins+1
    b2       = nbin-num_highl_junk_bins
    lmin_fit = bins(b1,1)
    lmax_fit = bins(b2,2)
    
    cls = 0.d0
!    cls(2:lmax,:,:)            = cls_fix(2:lmax,:,:)
    cls(l1:l2,:,:)            = cls_fix(l1:l2,:,:)  
    cls(lmin_fit:lmax_fit,:,:) = 0.d0
    call compute_lnL_and_chisq_full(info, cls, lnL, chisq)

    if (info%myid == 0) then
       ! Output change in likelihood, "chi-square" and PTE
       write(*,fmt='(f10.3,i6,i6,e16.4,e16.4,i8,e16.4,e16.4)') 0.5*(lmin_fit+lmax_fit), lmin_fit, lmax_fit, &
            & lnL0-lnL, 2*(lnL0-lnL), nint(nu*(b2-b1+1)), 0.d0, lnL
       write(unit,fmt='(f10.3,i6,i6,e16.4,e16.4,i8,e16.4)') 0.5*(lmin_fit+lmax_fit), lmin_fit, lmax_fit, &
            & lnL0-lnL, 2*(lnL0-lnL), nint(nu*(b2-b1+1)), 0.d0
    end if


    ! Set up covariance matrix for full range, this time setting all bins to zero
!    b1       = num_lowl_junk_bins+1
!    b2       = nbin-num_highl_junk_bins
!    lmin_fit = bins(b1,1)
!    lmax_fit = bins(b2,2)
!    
!    cls = 0.d0
!    call compute_lnL_and_chisq_full(info, cls, lnL, chisq)

!    if (info%myid == 0) then
!       ! Output change in likelihood, "chi-square" and PTE
!       write(*,fmt='(f10.3,i6,i6,e16.4,e16.4,i8,e16.4,e16.4)') 0.5*(lmin_fit+lmax_fit), lmin_fit, lmax_fit, &
!            & lnL0-lnL, 2*(lnL0-lnL), nint(nu*(b2-b1+1)), 0.d0, lnL
!       write(unit,fmt='(f10.3,i6,i6,e16.4,e16.4,i8,e16.4)') 0.5*(lmin_fit+lmax_fit), lmin_fit, lmax_fit, &
!            & lnL0-lnL, 2*(lnL0-lnL), nint(nu*(b2-b1+1)), 0.d0
!    end if

    if (info%myid == 0) close(unit)

    deallocate(cls)
    do d = 1, num_data_set
       call sc_dealloc(C_cov(d))
    end do
    deallocate(C_cov)
    
  end subroutine


end module map2cl_parfit_mod
