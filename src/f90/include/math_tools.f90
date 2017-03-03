module math_tools
  use healpix_types

   type cg_search
      real(dp), dimension(:), allocatable :: r, p, b, z, x
      real(dp)                            :: rz, rz0, err
      integer(i4b)                        :: i, n
   end type

   type bcg_search
      real(dp), dimension(:), allocatable :: r, rh, p, b, z, x, t, y, s, v
      real(dp)                            :: alpha, w, rho, err, bb
      integer(i4b)                        :: i, n
   end type

  interface invert_matrix
     module procedure invert_matrix_dpc, invert_matrix_dp, invert_matrix_sp
  end interface

  interface invert_matrix_with_mask
     module procedure invert_matrix_with_mask_dpc, invert_matrix_with_mask_dp
  end interface

  interface solve_linear_system
     module procedure solve_linear_system_dp, solve_linear_system_eigen, solve_linear_system_eigen_with_mask
  end interface

  interface qnorm
     module procedure ppnd7, ppnd16
  end interface

  interface get_diag
     module procedure get_diag_dp, get_diag_sp
  end interface

contains

  subroutine solve_linear_system_dp(A, X, B, status)
    
    real(dp), dimension(1:,1:), intent(in)  :: A
    real(dp), dimension(1:),    intent(in)  :: B
    real(dp), dimension(1:),    intent(out) :: X
    integer(i4b), optional,     intent(out) :: status

    integer(i4b) :: N, nrhs, lda, ldb, info
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(dp),     allocatable, dimension(:,:) :: b_int, A_int

    if(present(status)) status = 0
    N    = size(X)
    nrhs = 1
    lda  = n
    ldb  = n
    info = 0

    allocate(A_int(N,N))
    allocate(b_int(N,1))
    allocate(ipiv(N))

    A_int      = A
    b_int(:,1) = B

    call dgesv(N, nrhs, A_int, lda, ipiv, b_int, ldb, info)
    if (info /= 0) then
       if(present(status)) then
          status = info
          return
       else
          write(*,*) 'Error in solution of real system. Info = ', info
          stop
       end if
    end if

    X = b_int(:,1)

    deallocate(ipiv)
    deallocate(A_int)
    deallocate(b_int)
  end subroutine solve_linear_system_dp


  subroutine solve_linear_system_eigen(myid, A, X, B, cutoff, eigenvals)
    
    integer(i4b),               intent(in)            :: myid
    real(dp),                   intent(in)            :: cutoff
    real(dp), dimension(1:,1:), intent(in)            :: A
    real(dp), dimension(1:),    intent(in)            :: B
    real(dp), dimension(1:),    intent(out)           :: X
    real(dp), dimension(1:),    intent(out), optional :: eigenvals

    integer(i8b)     :: n
    real(dp)         :: cutoff_int
    real(dp),     allocatable, dimension(:,:) :: A_int
    real(dp),     allocatable, dimension(:)   :: W

    n      = size(X)
       
    ! Perform eigenvalue decomposition
    allocate(A_int(n,n))
    allocate(W(n))
    A_int = A
    call get_eigen_decomposition(myid, A, W, A_int)

    if (present(eigenvals)) then
       ! Output eigenvalues
       eigenvals = W
    end if
       
    ! Invert matrix
    cutoff_int = W(n) * cutoff
    where (W > cutoff_int)
       W = 1.d0 / W
    elsewhere
       W = 0.d0
    end where

    ! Solve equations
    X = matmul(A_int, W * matmul(transpose(A_int), B))

    deallocate(A_int)
    deallocate(W)

  end subroutine solve_linear_system_eigen


  subroutine solve_linear_system_eigen_with_mask(myid, A, X, B, mask, cutoff, eigenvals)
    
    integer(i4b),                   intent(in)            :: myid
    real(dp),                       intent(in)            :: cutoff
    real(dp),     dimension(1:,1:), intent(in)            :: A
    real(dp),     dimension(1:),    intent(in)            :: B
    real(dp),     dimension(1:),    intent(out)           :: X
    logical(lgt), dimension(1:),    intent(in)            :: mask
    real(dp),     dimension(1:),    intent(out), optional :: eigenvals

    integer(i4b)     :: i, j, n, m, ip, jp
    real(dp)         :: cutoff_int
    real(dp),     allocatable, dimension(:,:) :: A_int, V
    real(dp),     allocatable, dimension(:)   :: X_int, W, rhs

    n      = size(X)
    m      = count(mask)
       
    ! Perform eigenvalue decomposition
    allocate(A_int(m,m))
    allocate(V(m,m))
    allocate(W(m))
    allocate(rhs(m))
    allocate(X_int(m))
    
    ip = 0
    do i = 1, n
       if (mask(i)) then
          ip      = ip+1
          rhs(ip) = B(i)

          jp = 0
          do j = 1, n
             if (mask(j)) then
                jp = jp+1
                A_int(ip,jp) = A(i,j)
             end if
          end do
       end if
    end do
       
    call get_eigen_decomposition(myid, A_int, W, V)

    if (present(eigenvals)) then
       ! Output eigenvalues
       eigenvals = 0.d0
       eigenvals(1:m) = W
    end if
       
    ! Invert matrix
    cutoff_int = maxval(W) * cutoff
    where (W > cutoff_int)
       W = 1.d0 / W
    elsewhere
       W = 0.d0
    end where

    ! Solve equations
    X_int = matmul(V, W * matmul(transpose(V), rhs))

    ip = 0
    X  = 0.d0
    do i = 1, n
       if (mask(i)) then
          ip   = ip+1
          X(i) = X_int(ip)
       end if
    end do

    deallocate(A_int)
    deallocate(W)
    deallocate(V)
    deallocate(X_int)
    deallocate(rhs)

  end subroutine solve_linear_system_eigen_with_mask



  subroutine get_eigen_decomposition(myid, matrix, eigenvals, eigenvectors)
    
    integer(i4b),               intent(in)            :: myid
    real(dp), dimension(1:,1:), intent(in)            :: matrix
    real(dp), dimension(1:),    intent(out)           :: eigenvals
    real(dp), dimension(1:,1:), intent(out)           :: eigenvectors

    integer(i8b)     :: i, n, liwork, lwork, lda, ldb, info
    character(len=1) :: job, uplo
    real(dp)         :: cutoff_int
    real(dp),     allocatable, dimension(:,:) :: A_int
    real(dp),     allocatable, dimension(:)   :: W, work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    job    = 'v'
    uplo   = 'l'
    n      = size(eigenvals)
    lda    = n
    ldb    = n
    liwork = 5*n + 3
    lwork  = 2*n**2 + 6*n + 1
    info   = 0
      
!    write(*,*) myid, "Eigenvalue work array has size", lwork*8, " B"
    ! Perform eigenvalue decomposition
    allocate(work(lwork))
    allocate(iwork(liwork))
!    write(*,120) myid, n
!120 format('   Myid = ', I4, ' -- number of elements in eigen-decomposition  = ', I6)
    call cpu_time(t1)
    eigenvectors = matrix
    call dsyevd(job, uplo, n, eigenvectors, lda, eigenvals, work, lwork, iwork, liwork, info)
    call cpu_time(t2)
!    if (info /= 0) write(*,*) 'dsyevd info = ', info
!    write(*,130) myid, real(t2-t1,sp)
!130 format('   Myid = ', I4, ' -- CPU time for eigenvalue decomposition    = ', F8.2, ' sec')

    deallocate(work)
    deallocate(iwork)

  end subroutine get_eigen_decomposition


  subroutine get_eigenvalues(A, eigenvals)
    
    real(dp), dimension(1:,1:), intent(in)            :: A
    real(dp), dimension(1:),    intent(out)           :: eigenvals

    integer(i4b)     :: n, liwork, lwork, lda, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:,:) :: A_copy
    real(dp),     allocatable, dimension(:)   :: W, work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    n      = size(eigenvals)

    if (n == 1) then

       eigenvals(1) = A(1,1)

    else if (n == 2) then

       eigenvals(1) = 0.5d0*(A(1,1)+A(2,2) + sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))
       eigenvals(2) = 0.5d0*(A(1,1)+A(2,2) - sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))

    else

       job    = 'n'
       uplo   = 'l'
       lda    = n
       liwork = 1
       lwork  = 2*n+1
       info   = 0

       ! Perform eigenvalue decomposition
       allocate(work(lwork))
       allocate(iwork(liwork))
       allocate(A_copy(n,n))
       A_copy = A
       call dsyevd(job, uplo, n, A_copy, lda, eigenvals, work, lwork, iwork, liwork, info)
       if (info /= 0) write(*,*) 'get_eigenvalues -- dsyevd info = ', info
       
       deallocate(work)
       deallocate(iwork)
       deallocate(A_copy)

    end if

  end subroutine get_eigenvalues

  subroutine invert_singular_matrix(matrix, threshold)
    implicit none

    real(dp), dimension(1:,1:), intent(inout)         :: matrix
    real(dp),                   intent(in)            :: threshold
    
    real(dp), allocatable, dimension(:)               :: eigenvals
    real(dp), allocatable, dimension(:,:)             :: eigenvectors, matrix2
    integer(i4b)                                      :: myid, i, n    
    real(dp)                                          :: maxeigenval  

    myid = 1000
    n = size(matrix(1,:))
    allocate(matrix2(n,n))
    allocate(eigenvals(n))
    allocate(eigenvectors(n,n))
    call get_eigen_decomposition(myid, matrix, eigenvals, eigenvectors)
    maxeigenval = maxval(eigenvals)
    do i = 1, n
       if (eigenvals(i) == 0.d0) then 
          cycle
       else if (eigenvals(i) > threshold * maxeigenval .or. threshold ==0.d0) then
          eigenvals(i) = 1/eigenvals(i)
       else
          eigenvals(i) = 0.d0
       end if
    end do

    matrix = transpose(eigenvectors)
    do i = 1, n
       matrix(i,:) = matrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, matrix, n, 0.d0, matrix2, n)
    matrix = matrix2

    deallocate(eigenvals)
    deallocate(eigenvectors)
    deallocate(matrix2)

  end subroutine invert_singular_matrix


  subroutine invert_singular_matrix_top(matrix, threshold)
    implicit none

    real(dp), dimension(1:,1:), intent(inout)         :: matrix
    real(dp),                   intent(in)            :: threshold
    
    real(dp), allocatable, dimension(:)               :: eigenvals
    real(dp), allocatable, dimension(:,:)             :: eigenvectors, matrix2
    integer(i4b)                                      :: myid, i, n    
    real(dp)                                          :: mineigenval  

    myid = 1000
    n = size(matrix(1,:))
    allocate(matrix2(n,n))
    allocate(eigenvals(n))
    allocate(eigenvectors(n,n))
    call get_eigen_decomposition(myid, matrix, eigenvals, eigenvectors)
    mineigenval = minval(eigenvals)
    do i = 1, n
       if (eigenvals(i) == 0.d0) then 
          cycle
       else if (eigenvals(i) < threshold * mineigenval) then
          eigenvals(i) = 1/eigenvals(i)
       else
          eigenvals(i) = 0.d0
       end if
    end do

    matrix = transpose(eigenvectors)
    do i = 1, n
       matrix(i,:) = matrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, matrix, n, 0.d0, matrix2, n)
    matrix = matrix2

    deallocate(eigenvals)
    deallocate(eigenvectors)
    deallocate(matrix2)

  end subroutine invert_singular_matrix_top
  


  subroutine invert_matrix_dpc(matrix)
    implicit none

    complex(dpc), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    complex(dpc), allocatable, dimension(:)   :: work

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    allocate(ipiv(n))
    allocate(work(n))

    call ZGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'ZGETRF: LU factorization failed. Info = ', info
       stop
    else

       call ZGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'ZGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_dpc


  subroutine invert_matrix_dp(matrix, cholesky, status)
    implicit none

    real(dp), dimension(1:,1:), intent(inout)         :: matrix
    logical(lgt),               intent(in),  optional :: cholesky
    integer(i4b),               intent(out), optional :: status

    integer(i4b)     :: i, j, n, lda, info, lwork
    logical(lgt)     :: use_cholesky
    character(len=1) :: uplo
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(dp),     allocatable, dimension(:)   :: work

    if(present(status)) status = 0
    use_cholesky = .false.; if (present(cholesky)) use_cholesky = cholesky
    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    uplo  = 'l'
    allocate(ipiv(n))
    allocate(work(n))

    if (use_cholesky) then
       call DPOTRF(uplo, n, matrix, lda, info)
    else
       call DGETRF(n, n, matrix, lda, ipiv, info)
    end if
    if (info /= 0) then
       if(present(status)) then
          status = info
          return
       end if
       write(*,*) 'DGETRF: Factorization failed. Info = ', info
       stop
    else

       if (use_cholesky) then
          call DPOTRI(uplo, n, matrix, lda, info)
       else
          call DGETRI(n, matrix, lda, ipiv, work, lwork, info)
       end if

       if (info /= 0) then
          if(present(status)) then
             status = info
             return
          end if
          write(*,*) 'DGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    if (use_cholesky) then
       do i = 1, n
          do j = i+1, n
             matrix(i,j) = matrix(j,i)
          end do
       end do
    end if

    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_dp

  subroutine invert_matrix_sp(matrix)
    implicit none

    real(sp), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(sp),     allocatable, dimension(:)   :: work

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    allocate(ipiv(n))
    allocate(work(n))

    call SGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'SGETRF: LU factorization failed. Info = ', info
       stop
    else

       call SGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'SGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_sp


  subroutine invert_matrix_with_mask_dpc(matrix, cholesky)
    implicit none

    complex(dpc), dimension(1:,1:), intent(inout) :: matrix
    logical(lgt),               intent(in), optional :: cholesky

    integer(i4b)     :: i, j, n, lda, info, lwork
    logical(lgt)     :: use_cholesky
    character(len=1) :: uplo
    integer(i4b), allocatable, dimension(:)   :: ipiv
    complex(dpc), allocatable, dimension(:)   :: work
    logical(lgt), allocatable, dimension(:)   :: mask

    use_cholesky = .false.; if (present(cholesky)) use_cholesky = cholesky
    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    uplo  = 'l'
    allocate(ipiv(n))
    allocate(work(n))
    allocate(mask(n))

    mask = .true.
    do i = 1, n
       if (abs(matrix(i,i)) <= 0.d0) then
          mask(i)     = .false.
          matrix(i,i) = cmplx(1.d0,0.d0,dp)
       end if
    end do

    if (use_cholesky) then
       call ZPOTRF(uplo, n, matrix, lda, info)
    else
       call ZGETRF(n, n, matrix, lda, ipiv, info)
    end if
    if (info /= 0) then
       write(*,*) 'ZGETRF/ZPOTRF: LU factorization failed. Info = ', info
       stop
    else

       if (use_cholesky) then
          call ZPOTRI(uplo, n, matrix, lda, info)
       else
          call ZGETRI(n, matrix, lda, ipiv, work, lwork, info)
       end if

       if (info /= 0) then
          write(*,*) 'ZGETRI/ZPOTRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    do i = 1, n
       if (.not. mask(i)) then
          matrix(i,i) = cmplx(0.d0,0.d0,dp)
       end if
    end do

    if (use_cholesky) then
       do i = 1, n
          do j = i+1, n
             matrix(i,j) = conjg(matrix(j,i))
          end do
       end do
    end if

    deallocate(mask)
    deallocate(work)
    deallocate(ipiv)

  end subroutine invert_matrix_with_mask_dpc
  

  subroutine invert_matrix_with_mask_dp(matrix, cholesky)
    implicit none

    real(dp), dimension(1:,1:), intent(inout) :: matrix
    logical(lgt),               intent(in), optional :: cholesky

    integer(i4b)     :: i, j, n, lda, info, lwork
    logical(lgt)     :: use_cholesky
    character(len=1) :: uplo
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(dp),     allocatable, dimension(:)   :: work
    logical(lgt), allocatable, dimension(:)   :: mask

    use_cholesky = .false.; if (present(cholesky)) use_cholesky = cholesky
    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    uplo  = 'l'
    allocate(ipiv(n))
    allocate(work(n))
    allocate(mask(n))

    mask = .true.
    do i = 1, n
       if (abs(matrix(i,i)) <= 0.d0) then
          mask(i)     = .false.
          matrix(i,i) = 1.d0
       end if
    end do

    if (use_cholesky) then
       call DPOTRF(uplo, n, matrix, lda, info)
    else
       call DGETRF(n, n, matrix, lda, ipiv, info)
    end if
    if (info /= 0) then
       write(*,*) 'DGETRF: LU factorization failed. Info = ', info
       stop
    else

       if (use_cholesky) then
          call DPOTRI(uplo, n, matrix, lda, info)
       else
          call DGETRI(n, matrix, lda, ipiv, work, lwork, info)
       end if

       if (info /= 0) then
          write(*,*) 'DGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    do i = 1, n
       if (.not. mask(i)) then
          matrix(i,i) = 0.d0
       end if
    end do

    if (use_cholesky) then
       do i = 1, n
          do j = i+1, n
             matrix(i,j) = matrix(j,i)
          end do
       end do
    end if

    deallocate(work)
    deallocate(ipiv)
    deallocate(mask)

  end subroutine invert_matrix_with_mask_dp

  subroutine invert_matrix_with_mask_sp(matrix)
    implicit none

    real(sp), dimension(1:,1:), intent(inout) :: matrix

    integer(i4b)     :: i, n, lda, info, lwork
    integer(i4b), allocatable, dimension(:)   :: ipiv
    real(sp),     allocatable, dimension(:)   :: work
    logical(lgt), allocatable, dimension(:)   :: mask

    n     = size(matrix(1,:))
    lda   = n
    lwork = n
    info  = 0
    allocate(ipiv(n))
    allocate(work(n))
    allocate(mask(n))

    mask = .true.
    do i = 1, n
       if (abs(matrix(i,i)) <= 0.0) then
          mask(i)     = .false.
          matrix(i,i) = 1.0
       end if
    end do

    call SGETRF(n, n, matrix, lda, ipiv, info)
    if (info /= 0) then
       write(*,*) 'SGETRF: LU factorization failed. Info = ', info
       stop
    else

       call SGETRI(n, matrix, lda, ipiv, work, lwork, info)

       if (info /= 0) then
          write(*,*) 'SGETRI: Inversion failed. Info = ', info
          stop
       end if

    end if

    do i = 1, n
       if (.not. mask(i)) then
          matrix(i,i) = 0.0
       end if
    end do

    deallocate(work)
    deallocate(ipiv)
    deallocate(mask)

  end subroutine invert_matrix_with_mask_sp




  !------------------------------------------------------------------
  ! Subroutines for inverting a matrix. Based on 
  ! Bo Einarssons F90-manual.
  !------------------------------------------------------------------

  subroutine solve_system(A, X, B)
    implicit none

    complex(dpc), dimension (:, :)               :: A
    complex(dpc), dimension (:)                  :: X
    complex(dpc), dimension (:)                  :: B

    complex(dpc), dimension(size(B), size(B)+1)  :: m
    integer, dimension (1)                       :: max_loc
    complex(dpc), dimension(size(B)+1)           :: temp_row
    integer                                      :: N, K, I 
    N = size (B)
    m (1:N, 1:N) = A
    m (1:N, N+1) = B 
    
    do K = 1, N - 1

       max_loc = maxloc (abs (m (K:N, K)))
       if ( max_loc(1) /= 1 ) then
          temp_row (K:N+1 ) = m (K, K:N+1)
          m (K, K:N+1)= m (K-1+max_loc(1), K:N+1)
          m (K-1+max_loc(1), K:N+1) = temp_row( K:N+1)
       end if

       temp_row (K+1:N) = m (K+1:N, K) / m (K, K)
       do I = K+1, N
          m (I, K+1:N+1) = m (I, K+1:N+1) - &
               temp_row (I) * m (K, K+1:N+1)
       end do
       m (K+1:N, K) = cmplx(0.d0,0.d0)

    end do 

    do K = N, 1, -1
       X (K) = ( m (K, N+1) - &
            sum (m (K, K+1:N) * X (K+1:N)) ) / m (K, K)
    end do

  end subroutine 


  !------------------------------------------------------------------
  ! Subroutines for inverting a matrix. Based on 
  ! Bo Einarssons F90-manual.
  !------------------------------------------------------------------

  subroutine solve_system_real(A, X, B)
    
    real(dp), dimension (:, :)               :: A
    real(dp), dimension (:)                  :: X
    real(dp), dimension (:)                  :: B

    real(dp), dimension(size(B), size(B)+1)  :: m
    integer, dimension (1)                       :: max_loc
    real(dp), dimension(size(B)+1)           :: temp_row
    integer                                      :: N, K, I 

    write(*,*)"WARNING! solve_system_real() should not be used"
    write(*,*)"use the _much_ faster solve_linear_system() !WARNING"

    N = size (B)
    m (1:N, 1:N) = A
    m (1:N, N+1) = B 
    
    do K = 1, N - 1

       max_loc = maxloc (abs (m (K:N, K)))
       if ( max_loc(1) /= 1 ) then
          temp_row (K:N+1 ) = m (K, K:N+1)
          m (K, K:N+1)= m (K-1+max_loc(1), K:N+1)
          m (K-1+max_loc(1), K:N+1) = temp_row( K:N+1)
       end if

       temp_row (K+1:N) = m (K+1:N, K) / m (K, K)
       do I = K+1, N
          m (I, K+1:N+1) = m (I, K+1:N+1) - &
               temp_row (I) * m (K, K+1:N+1)
       end do
       m (K+1:N, K) = 0.d0

    end do 

    do K = N, 1, -1
       X (K) = ( m (K, N+1) - &
            sum (m (K, K+1:N) * X (K+1:N)) ) / m (K, K)
    end do

  end subroutine 

  subroutine invert_matrix_real(m, len)
    implicit none

    real(dp), dimension(:,:)     :: m
    integer(i4b)                     :: len

    integer(i4b)        :: i, j

    real(dp), dimension(len, len)     :: temp_m
    real(dp), dimension(len)          :: temp_b

    do i = 1, len
       do j = 1, len
          if (j == i) then
             temp_b(j) = 1.d0
          else
             temp_b(j) = 0.d0
          end if
       end do

       call solve_system_real(m, temp_m(:, i), temp_b)
    end do

    m = temp_m

  end subroutine invert_matrix_real


  subroutine cholesky_decompose(A, L, status)
    implicit none
    
    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: L
    integer(i4b), optional,   intent(out) :: status
    
    integer(i4b) :: N, i, j, k, stat
    real(dp) :: temp
    real(dp), allocatable, dimension(:) :: temp_row
    
    if(present(status)) status = 0
    N = size(A(1,:))
    
    L = A
    call dpotrf( 'L', N, L, N, stat )

    if (stat /= 0) then
       if(present(status)) then
          status = stat
          return
       end if
       write(*,*) 'Cholesky decomposition failed. stat = ', stat
       stop
    end if

    do i = 1, N
       do j = i+1, N
          L(i,j) = 0.d0
       end do
    end do

  end subroutine cholesky_decompose


  subroutine cholesky_solve_lapack(L, x, b)
    implicit none
    
    real(dp), dimension(:,:), intent(in)  :: L
    real(dp), dimension(:),   intent(in)  :: b
    real(dp), dimension(:),   intent(out) :: x
    
    integer(i4b) :: N, i, j, k, stat
    real(dp) :: temp
    real(dp), allocatable, dimension(:) :: temp_row
    
    N = size(L(1,:))
    
    x = b
    call dpotrs( 'L', N, 1, L, N, x, N, stat)

  end subroutine cholesky_solve_lapack


  subroutine cholesky_decompose_with_mask(A, L)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: A
    real(dp), dimension(1:,1:), intent(out) :: L

    real(dp)         :: max_value
    integer(i4b)     :: i, j, n, stat
    logical(lgt), allocatable, dimension(:)   :: mask

    max_value = maxval(abs(A))
    n         = size(A(1,:))
    allocate(mask(n))

    L = A

    mask = .true.
    do i = 1, n
       if (abs(A(i,i))/max_value <= 1.d-12 .or. max_value == 0d0) then
          mask(i)     = .false.
          L(i,:) = 0.d0
          L(:,i) = 0.d0
          L(i,i) = 1.d0
       end if
    end do

    call dpotrf( 'L', n, L, n, stat )

    if (stat /= 0) then
       write(*,*) 'Cholesky decomposition failed. stat = ', stat
       stop
    end if

    do i = 1, N
       do j = i+1, N
          L(i,j) = 0.d0
       end do
    end do

    do i = 1, n
       if (.not. mask(i)) then
          L(i,:) = 0.d0
          L(:,i) = 0.d0
       end if
    end do

    deallocate(mask)

  end subroutine cholesky_decompose_with_mask



  subroutine cholesky_solve(L, b, x)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: L
    real(dp), dimension(:),   intent(in)  :: b
    real(dp), dimension(:),   intent(out) :: x

    integer(i4b) :: N, i, j, k
    real(dp) :: temp
    real(dp), allocatable, dimension(:) :: temp_row

    N = size(L(1,:))

    x = 0.d0

    do j = 1, N
       
       temp = 0.d0
       do i = 1, j-1
          temp = temp + L(j,i) * x(i)
       end do

       x(j) = (b(j) - temp) / L(j,j)

    end do

  end subroutine cholesky_solve

  subroutine matmul_symm(side, A, B, C)
    implicit none

    character(len=1),                 intent(in)  :: side
    real(dp),         dimension(:,:), intent(in)  :: A, B
    real(dp),         dimension(:,:), intent(out) :: C

    integer(i4b) :: n

    n = size(A(:,1))

    call dsymm(side, 'l', n, n, 1.d0, A, n, B, n, 0.d0, C, n )

  end subroutine matmul_symm

  subroutine matmul_gen(A, B, C)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: A, B
    real(dp), dimension(:,:), intent(out) :: C
    integer(i4b) :: m, n, k

    m = size(A,1)
    n = size(B,2)
    k = size(A,2)

    call dgemm('N','N',m,n,k,1.d0,A,m,B,k,0.d0,C,m)

  end subroutine matmul_gen


  !  THIS routine returns the value of normalised P_lm(theta) such that
  !  2.*PI*Integral_-1^+1 dx P_lm(x)^2 = 1., where x = cos(theta)
  !  
  !  modified P_lm generating routine from HEALPix, by K.M.G. 4, Sept. 2000
  
  !=======================================================================
  subroutine comp_normalised_Plm(nlmax, m, theta, plm)
    !=======================================================================
    !nlmax (integer) = maximum l
    !theta in radians (double precision)
    !lambda (double precision) = modulus of the complex spherical harmonics
    !  contains lambda(l,m) with l,m in [0,nlmax]
    !  contains lambda(l-m) with l-m in [0,nlmax-m]
    !=======================================================================
    IMPLICIT none
    !
    REAL(dp), PARAMETER:: max_dp  = HUGE(1.0d0)
    REAL(DP), PARAMETER :: PI = 3.141592653589793238462643383279502884197d0
    
    INTEGER(I4B), INTENT(IN)  :: m
    INTEGER(I4B), INTENT(IN)  :: nlmax
    REAL(DP),     INTENT(IN)  :: theta
    
    
    REAL(DP) :: lambda(0:nlmax)
    REAL(DP),     DIMENSION(0:nlmax), INTENT(OUT) :: plm
    
    INTEGER(I4B) :: nmmax
    INTEGER(I4B) l, ith, indl, mm               !, m ...  alm related
    
    REAL(DP) sq4pi_inv
    REAL(DP) cth, sth
    REAL(DP) a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) f2m, fm2, fl2
    
    Character(LEN=7), PARAMETER :: code = 'ALM2MAP'
    INTEGER(I4B) :: status
    
    REAL(DP), PARAMETER :: bignorm = 1.d-20*max_dp
    !=======================================================================
    
    
    !      write(*,*)'   PI   =    ',PI
    
    nmmax = nlmax
    
    LAMBDA = 0.0d0
    
    !     --------------------------------------------
    sq4pi_inv = 1.D0 / SQRT(4.D0*PI)
    
    cth = COS(theta)
    sth = SIN(theta)
    
    plm=0.d0
    
    !      write(*,*)cth,sth
    
    !        lambda_mm tends to go down when m increases (risk of underflow)
    !        lambda_lm tends to go up   when l increases (risk of overflow)
    
    lam_mm = sq4pi_inv * bignorm ! lambda_0_0 * big number --> to avoid underflow
    
    !      do m = 0, nmmax
    
    fm2 = DFLOAT(m) **2
    
    !           ---------- l = m ----------
    par_lm = 1.d0  ! = (-1)^(l+m)
    if (m .ge. 1) then ! lambda_0_0 for m>0
       do mm = 1, m
          f2m = 2.d0 * mm
          lam_mm = - lam_mm * sth * DSQRT( (f2m+1.d0)/ f2m )
       enddo
    endif
    
    lam_lm = lam_mm / bignorm ! actual lambda_mm
    
    LAMBDA(M) = LAM_LM
    
    
    !           ---------- l > m ----------
    lam_0 = 0.d0
    lam_1 = 1.d0 / bignorm    ! small number --> to avoid overflow
    fl2 = DFLOAT(m+1) **2
    a_rec = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
    lam_2 = cth * lam_1 * a_rec
    do l = m+1, nmmax
       par_lm = - par_lm  ! = (-1)^(l+m)
       
       lam_lm = lam_2 * lam_mm ! actual lambda_lm (small and big numbers cancel out)
       
       !            lambda(l,m) = lam_lm
       !            lambda(l-m) = lam_lm
       
       LAMBDA(L) = LAM_LM
       
       lam_0 = lam_1 / a_rec
       lam_1 = lam_2
       fl2 = DFLOAT(l+1) **2
       a_rec = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       lam_2 = (cth * lam_1 - lam_0) * a_rec
    enddo
    
    !      enddo
    
    
    plm = lambda
    
    return
  end subroutine comp_normalised_Plm


  subroutine convert_sigma2frac(frac, sigma)
    implicit none

    real(sp), intent(in)  :: sigma
    real(sp), intent(out) :: frac

    real(dp) :: int_frac, int_sigma

    int_sigma = real(sigma,dp)

    int_frac = corr_erf(int_sigma / sqrt(2.d0))

    frac = (1. + real(int_frac,sp)) / 2.

  end subroutine convert_sigma2frac

  ! Solve the equation erf(x/sqrt(2)) - y = 0 using bisection
  subroutine convert_fract2sigma(sigma, fract)
    implicit none

    real(sp), intent(in)  :: fract
    real(sp), intent(out) :: sigma

    integer(i4b)   :: maxit, j
    real(dp)       :: dx, fr, fl, fm, xl, xr, xm

    maxit = 40

    if (fract > 0.999999426) then
       sigma = 5.
    else if (fract < 1e-7) then
       sigma = 0.
    else

       xl = 0.
       xr = 10.
       xm = (xl+xm)/2.

       fl = -fract
       fr = corr_erf(xr)-fract

       do j = 1, maxit
          if (abs(xl-xr) < 1e-7) exit

          xm = (xl+xr)/2.
          fm = corr_erf(xm)-fract

          if (fm*fl < 0.0) then
             xr = xm
             fr = fm
          else
             xl = xm
             fl = fm
          end if

       end do

       if (j == maxit) then
          write(*,*) 'ERROR: Too many iterations in the fract2sigma search'
       end if

       sigma = sqrt(2.) * xm

    end if

  end subroutine convert_fract2sigma

  ! Computes the error function. Borrowed from Numerical Recipes
  real(dp) function corr_erf(x)
    implicit none

    real(dp), intent(in) :: x

    real(dp) :: ans, z, t

    z = abs(x)
    t=1./(1.+0.5*z)
    ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+&
         & t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+&
         & t*(-0.82215223+t*0.17087277)))))))))
    if (x >= 0.) then
       corr_erf = 1.-ans
    else
       corr_erf = ans - 1.
    end if

  end function corr_erf

  real(dp) function gammln(xx)
    implicit none

    real(dp), intent(in) :: xx

    real(dp) :: x, tmp
    real(dp) :: stp = 2.5066282746310005d0
    integer(i4b) :: i

    real(dp), dimension(6) :: coef

    coef(1) = 76.18009172947146d0
    coef(2) = -86.50532032941677d0
    coef(3) = 24.01409824083091d0
    coef(4) = -1.231739572450155d0
    coef(5) = 0.001208650973866179d0
    coef(6) = -0.000005395239384953d0

    x = xx
    tmp = x + 5.5d0
    tmp = (x+0.5d0) * log(tmp) - tmp
    gammln = 0.d0
    do i = 1, 6
       gammln = gammln + coef(i)/(x+real(i,dp))
    end do
    gammln = tmp+log(stp*(gammln + 1.000000000190015d0)/x)

  end function gammln

  function get_identity(n) result(M)
    implicit none
    integer(i4b)             :: n, i
    real(dp)                 :: M(n,n)
    M = 0.d0
    do i = 1, n
       M(i,i) = 1.d0
    end do
  end function get_identity

  function get_diag_dp(A) result(d)
    implicit none
    real(dp)     :: A(:,:), d(size(A,1))
    integer(i4b) :: i
    do i = 1, size(d)
       d(i) = A(i,i)
    end do
  end function
  function get_diag_sp(A) result(d)
    implicit none
    real(sp)     :: A(:,:), d(size(A,1))
    integer(i4b) :: i
    do i = 1, size(d)
       d(i) = A(i,i)
    end do
  end function

  function diag(d) result(A)
    implicit none
    real(dp)     :: d(:), A(size(d),size(d))
    integer(i4b) :: i
    A = 0
    do i = 1, size(d)
       A(i,i) = d(i)
    end do
  end function

  ! The following routines are very similar to the ones at the top, but
  ! have two differences: They are silent, and they take a status
  ! argument which is nonzero if anything goes wrong.
  subroutine solve_system_eiglim(A, b, x, cutoff, status)
    implicit none
    real(dp), intent(in)  :: A(:,:), b(:), cutoff
    real(dp), intent(out) :: x(:)
    real(dp)              :: c
    real(dp), allocatable :: V(:,:), W(:)
    integer(i4b) :: status, n
    n = size(A,1)
    status = 0
    if(n == 1) then
       x(1) = b(1)/A(1,1)
       return
    end if
    allocate(V(n,n), W(n))
    call eigen_decomp(A, W, V, status)
    if(status == 0) then
       c = maxval(W)*cutoff
       if(any(W < c)) status = 2
       where(W < c)
          W = 0
       elsewhere
          W = 1/W
       end where
       x = matmul(V, W * matmul(transpose(V), b))
    else
       status = 1
    end if
    deallocate(V, W)
  end subroutine

  subroutine eigen_decomp(matrix, eigenvals, eigenvectors, status)
    implicit none
    real(dp)                  :: matrix(:,:), eigenvals(:), eigenvectors(:,:)
    real(dp),     allocatable :: work(:)
    integer(i4b), allocatable :: iwork(:)
    integer(i4b)              :: n, status
    integer(i4b)              :: m
    n = size(matrix,1)
    m = 2*n**2+6*n+1
    allocate(work(m), iwork(5*n+3))
    eigenvectors = matrix
    status = 0
    call dsyevd('v', 'l', n, eigenvectors, n, eigenvals, work, m, &
      & iwork, size(iwork), status)
    deallocate(work, iwork)
  end subroutine eigen_decomp

  subroutine cholesky(A, status)
    implicit none
    real(dp)     :: A(:,:)
    integer(i4b) :: n, status
    n = size(A,1)
    call dpotrf('L', n, A, n, status)
  end subroutine

  subroutine invert_sym(A, status)
    implicit none
    real(dp)     :: A(:,:)
    integer(i4b) :: n, i, j, status
    n = size(A,1)
    call dpotrf('L', n, A, n, status)
    if(status /= 0) return
    call dpotri('L', n, A, n, status)
  end subroutine invert_sym

  function trace(mat) result(res)
    implicit none
    real(dp)     :: mat(:,:), res
    integer(i4b) :: i
    res = 0
    do i = 1, size(mat,1); res = res + mat(i,i); end do
  end function

  subroutine compute_hermitian_power(alpha, cutoff, A, B, W_out)
    implicit none

    real(dp),                 intent(in)  :: alpha, cutoff
    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: B
    real(dp), dimension(:),   intent(out), optional :: W_out

    real(dp) :: c
    real(dp), allocatable :: V(:,:), W(:)
    integer(i4b) :: status, n, i
    n = size(A,1)
    allocate(V(n,n), W(n))
    call eigen_decomp(A, W, V, status)
    if (present(W_out)) W_out = W
    if(status == 0) then
       c = maxval(W)*cutoff
       where(W < c)
          W = 0
       elsewhere
          W = W**(0.5d0*alpha)
       end where
       do i = 1, n
          V(:,i) = V(:,i) * W(i)
       end do
       call dgemm('N','T',n,n,n,1.d0,V,n,V,n,0.d0,B,n)
    else
       status = 1
       write(*,*) 'compute_hermitian_power -- cannot eigen decompose matrix'
       stop
    end if
    deallocate(V, W)

  end subroutine compute_hermitian_power

  subroutine fit_polynomial(x, y, a, status)
    implicit none

    real(dp), dimension(1:),  intent(in) :: x, y
    real(dp), dimension(0:),  intent(out) :: a
    integer(i4b), optional,   intent(out) :: status

    integer(i4b) :: i, j
    real(dp), dimension(0:size(a)-1)              :: b
    real(dp), dimension(0:size(a)-1, 0:size(a)-1) :: C

    do i = 0, size(a)-1
       b(i) = sum(y * x**i)
       do j = i, size(a)-1
          C(i,j) = sum(x**(i+j))
          C(j,i) = C(i,j)
       end do
    end do
    call solve_linear_system(C, a, b, status)

  end subroutine fit_polynomial

  ! xin and xout must be sorted.
  subroutine lin_interpol(xin, yin, xout, yout)
    implicit none
    real(dp)     :: xin(:), yin(:), xout(:), yout(:), x
    integer(i4b) :: i, j, k, m, n

    if(size(xin) <= 0) then
       ! yout = nan ! Should return nan, but nan not available without pulling in utils
       return
    elseif(size(xin) == 1) then
       yout = yin(1)
       return
    end if

    j = 1
    do i = 1, size(xout)
       if(xout(i) <= xin(1)) then
          yout(i) = yin(1)
       elseif(xout(i) >= xin(size(xin))) then
          yout(i) = yin(size(xin))
       else
          ! A [j-1,j] range with i inside it exists
          do while(xin(j) < xout(i))
             j = j+1
          end do
          if(xin(j) == xin(j-1)) then
             x = 0.5d0
          else
             x = (xout(i)-xin(j-1))/(xin(j)-xin(j-1))
          end if
          yout(i) = yin(j-1) + (yin(j)-yin(j-1))*x
       end if
    end do
  end subroutine

  ! CG search routines. Typical usage sketch, for
  ! system Ax=b
  !  call cg_init(cg, b)
  !  do
  !     call cg_pre(cg, matmul(iM, cg%r))
  !     call cg_mul(cg, matmul(A,  cg%p))
  !     if(maxval(abs(cg%r)) < 0.0001) exit
  !  end do
  !  call cg_free(cg)
  !
  ! The second argument of cg_pre is optional, so
  ! you don't need to pass it if you don't use preconditioning.
  ! But cg_pre must always be called before cg_mul anyway.
  !
  ! In practice, you most likely wouldn't use matmul,
  ! but some faster method.
  subroutine cg_init(cg, b)
    implicit none
    type(cg_search)    :: cg
    real(dp)           :: b(:)
    integer(i4b)       :: n
    call cg_free(cg)
    n    = size(b)
    cg%n = n
    cg%i = 0
    allocate(cg%r(n), cg%p(n), cg%x(n), cg%b(n), cg%z(n))
    cg%x  = 0
    cg%b  = b
    cg%r  = b
  end subroutine

  subroutine cg_pre(cg, iMr)
    implicit none
    type(cg_search)    :: cg
    real(dp), optional :: iMr(:)
    real(dp)           :: rz, beta
    if(present(iMr)) then
       cg%z = iMr
    else
       cg%z = cg%r
    end if
    if(cg%i == 0) then
       cg%rz  = dot_product(cg%r, cg%z)
       cg%rz0 = cg%rz
       cg%p = cg%z
    else
       rz    = dot_product(cg%r, cg%z)
       beta  = rz/cg%rz
       cg%rz = rz
       cg%p  = cg%z + cg%p*beta
    end if
  end subroutine

  subroutine cg_mul(cg, Ap)
    implicit none
    type(cg_search) :: cg
    real(dp)        :: Ap(:), alpha
    alpha = cg%rz/dot_product(cg%p, Ap)
    cg%x  = cg%x + cg%p*alpha
    cg%r  = cg%r - Ap  *alpha
    cg%err= sqrt(cg%rz/cg%rz0)
    cg%i  = cg%i + 1
  end subroutine

  subroutine cg_free(cg)
    implicit none
    type(cg_search) :: cg
    if(allocated(cg%r)) deallocate(cg%r)
    if(allocated(cg%p)) deallocate(cg%p)
    if(allocated(cg%x)) deallocate(cg%x)
    if(allocated(cg%b)) deallocate(cg%b)
    if(allocated(cg%z)) deallocate(cg%z)
  end subroutine

  ! Stabilized biconjugate gradient. This is pretty cumbersome,
  ! requireing not just pre and mul, but 3 different pre calls and two
  ! different mul calls. Usage:
  !
  ! call bcg_init(bcg, b)
  ! do
  !    call bcg1(bcg, iK*bcg%p)
  !    call bcg2(bcg,  A*bcg%y)
  !    call bcg3(bcg, iK*bcg%s)
  !    call bcg4(bcg,  A*bcg%z)
  !    call bcg5(bcg, iK*bcg%t)
  !    exit if bcg%x is good enough
  ! end do
  ! call bcg_free(bcg)
  !
  ! Here the equation system is Ax=b, and iK ~ 1/A is the preconditioner.
  subroutine bcg_init(bcg, b)
    implicit none
    type(bcg_search) :: bcg
    real(dp)         :: b(:)
    integer(i4b)     :: n
    call bcg_free(bcg)
    n = size(b)
    allocate(bcg%r(n),bcg%rh(n),bcg%p(n),bcg%b(n),bcg%z(n))
    allocate(bcg%x(n),bcg%t(n), bcg%y(n),bcg%s(n),bcg%v(n))
    bcg%b     = b
    bcg%r     = b
    bcg%rh    = b
    bcg%alpha = 1
    bcg%w     = 1
    bcg%rho   = dot_product(bcg%r,bcg%rh)
    bcg%p     = bcg%r
    bcg%bb    = bcg%rho
    bcg%x     = 0
    bcg%err   = 1
  end subroutine

  subroutine bcg_free(bcg)
    implicit none
    type(bcg_search) :: bcg
    if(allocated(bcg%r))  deallocate(bcg%r)
    if(allocated(bcg%rh)) deallocate(bcg%rh)
    if(allocated(bcg%p))  deallocate(bcg%p)
    if(allocated(bcg%b))  deallocate(bcg%b)
    if(allocated(bcg%z))  deallocate(bcg%z)
    if(allocated(bcg%x))  deallocate(bcg%x)
    if(allocated(bcg%t))  deallocate(bcg%t)
    if(allocated(bcg%y))  deallocate(bcg%y)
    if(allocated(bcg%s))  deallocate(bcg%s)
    if(allocated(bcg%v))  deallocate(bcg%v)
  end subroutine

  subroutine bcg1(bcg, Kp)
    implicit none
    type(bcg_search) :: bcg
    real(dp)         :: Kp(:)
    bcg%y = Kp
  end subroutine

  subroutine bcg2(bcg, Ay)
    implicit none
    type(bcg_search) :: bcg
    real(dp)         :: Ay(:)
    bcg%v     = Ay
    bcg%alpha = bcg%rho/dot_product(bcg%rh,bcg%v)
    bcg%s     = bcg%r - bcg%alpha*bcg%v
  end subroutine

  subroutine bcg3(bcg, Ks)
    implicit none
    type(bcg_search) :: bcg
    real(dp)         :: Ks(:)
    bcg%z = Ks
  end subroutine

  subroutine bcg4(bcg, Az)
    implicit none
    type(bcg_search) :: bcg
    real(dp)         :: Az(:)
    bcg%t = Az
  end subroutine

  subroutine bcg5(bcg, Kt)
    implicit none
    type(bcg_search) :: bcg
    real(dp)         :: Kt(:), beta, rho
    bcg%w   = dot_product(Kt,bcg%z)/dot_product(Kt,Kt)
    bcg%x   = bcg%x + bcg%alpha*bcg%y + bcg%w*bcg%z
    ! Update accuracy thing here
    bcg%r   = bcg%s - bcg%w*bcg%t
    rho     = dot_product(bcg%rh,bcg%r)
    beta    = rho/bcg%rho*bcg%alpha/bcg%w
    bcg%rho = rho
    bcg%p   = bcg%r + beta*(bcg%p-bcg%w*bcg%v)
    bcg%err = sqrt(sum(bcg%r**2)/bcg%bb)
  end subroutine

  subroutine get_legendre_polynomials(x, P)
    implicit none

    real(dp),                intent(in)  :: x
    real(dp), dimension(0:), intent(out) :: P

    integer(i4b) :: i, n

    n = size(P)-1

    if (n > -1) P(0) = 1.d0
    if (n > 0)  P(1) = x

    do i = 1, n-1
       P(i+1) = ((2.d0*real(i,dp)+1.d0) * x * P(i) - real(i,dp) * P(i-1)) / real(i+1,dp)
    end do

  end subroutine get_legendre_polynomials

  ! Compute mat**pow by using eigenvalue decomposition. Eigenvalues
  ! too small compared to the maximum one are set to zero to avoid
  ! numerical problems
  subroutine eigen_pow(mat, pow, res, neff, cutoff)
    implicit none
    real(dp),    intent(in)  :: mat(:,:), pow
    real(dp),    intent(out) :: res(:,:)
    real(dp),    intent(in),  optional :: cutoff
    integer(i4b),intent(out), optional :: neff
    integer(i4b)             :: i, j, k, n
    real(dp)                 :: mval, cut
    real(dp),    allocatable :: eigvals(:), eigvecs(:,:), tmp(:,:)
    cut = 1d-12; if(present(cutoff)) cut = cutoff
    n = size(mat,1)
    allocate(eigvals(n),eigvecs(n,n))
    call get_eigen_decomposition(0, mat, eigvals, eigvecs)
    mval = eigvals(n)
    if(all(eigvals==0)) mval = 1
    where(eigvals/mval < cut)
       eigvals = 0
    elsewhere
       eigvals = eigvals**pow
    end where
    if(present(neff)) neff = count(eigvals /= 0)
    allocate(tmp(n,n))
    do i = 1, n
       tmp(:,i) = eigvecs(:,i)*eigvals(i)
    end do
    call dgemm('N', 'T', n, n, n, 1d0, tmp, n, eigvecs, n, 0d0, res, n)
    deallocate(eigvals, eigvecs, tmp)
  end subroutine

  !ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-
  !484.
  !
  !Produces the normal deviate Z corresponding to a given lower
  !tail area of P; Z is accurate to about 1 part in 10**7.
  !
  !The hash sums below are the sums of the mantissas of the
  !coefficients.   They are included for use in checking
  !transcription.
  !
  real function ppnd7 (p, ifault)
    implicit none
    integer, optional :: ifault
    real zero, one, half, split1, split2, const1, const2, a0, a1, &
     &      a2, a3, b1, b2, b3, c0, c1, c2, c3, d1, d2, e0, e1, e2, &
     &      e3, f1, f2, p, q, r
    parameter (zero = 0.0, one = 1.0, half = 0.5, &
     &      split1 = 0.425, split2 = 5.0, &
     &      const1 = 0.180625, const2 = 1.6)
    !
    !coefficients for p close to 0.5
    !
    parameter (a0 = 3.3871327179e+00, a1 = 5.0434271938e+01, &
     &         a2 = 1.5929113202e+02, a3 = 5.9109374720e+01, &
     &         b1 = 1.7895169469e+01, b2 = 7.8757757664e+01, &
     &         b3 = 6.7187563600e+01)
    !hash sum ab    32.31845 77772
    !
    !coefficients for p not close to 0, 0.5 or 1.
    !
    parameter (c0 = 1.4234372777e+00, c1 = 2.7568153900e+00, &
     &         c2 = 1.3067284816e+00, c3 = 1.7023821103e-01, &
     &         d1 = 7.3700164250e-01, d2 = 1.2021132975e-01)
    !hash sum cd    15.76149 29821
    !
    !coefficients for p near 0 or 1.
    !
    parameter (e0 = 6.6579051150e+00, e1 = 3.0812263860e+00, &
     &         e2 = 4.2868294337e-01, e3 = 1.7337203997e-02, &
     &         f1 = 2.4197894225e-01, f2 = 1.2258202635e-02)
    !hash sum ef    19.40529 10204
    !
    if(present(ifault)) ifault = 0
    q = p - half
    if (abs(q) .le. split1) then
       r     = const1 - q * q
       ppnd7 = q * (((a3 * r + a2) * r + a1) * r + a0) / &
        & (((b3 * r + b2) * r + b1) * r + one)
    else
       if (q .lt. zero) then
          r = p
       else
          r = one - p
       end if
       if (r .le. zero) then
          if(present(ifault)) ifault = 1
          ppnd7 = zero
          return
       end if
       r = sqrt(-log(r))
       if (r .le. split2) then
          r     = r - const2
          ppnd7 = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one)
       else
          r     = r - split2
          ppnd7 = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one)
       end if
       if (q .lt. zero) ppnd7 = - ppnd7
    end if
  end function

  !algorithm as241  appl. statist. (1988) vol. 37, no. 3
  !
  !produces the normal deviate z corresponding to a given lower
  !tail area of p; z is accurate to about 1 part in 10**16.
  !
  !the hash sums below are the sums of the mantissas of the
  !coefficients.   they are included for use in checking
  !transcription.
  double precision function ppnd16 (p, ifault)
    integer, optional :: ifault
    double precision zero, one, half, split1, split2, const1, &
     & const2, a0, a1, a2, a3, a4, a5, a6, a7, b1, b2, b3, &
     &     b4, b5, b6, b7, &
     & c0, c1, c2, c3, c4, c5, c6, c7, d1, d2, d3, d4, d5, &
     & d6, d7, e0, e1, e2, e3, e4, e5, e6, e7, f1, f2, f3, &
     & f4, f5, f6, f7, p, q, r
    parameter (zero = 0.d0, one = 1.d0, half = 0.5d0, &
     & split1 = 0.425d0, split2 = 5.d0, &
     & const1 = 0.180625d0, const2 = 1.6d0)
    !coefficients for p close to 0.5
    parameter (a0 = 3.3871328727963666080d0, &
     &         a1 = 1.3314166789178437745d+2, &
     &         a2 = 1.9715909503065514427d+3, &
     &         a3 = 1.3731693765509461125d+4, &
     &         a4 = 4.5921953931549871457d+4, &
     &         a5 = 6.7265770927008700853d+4, &
     &         a6 = 3.3430575583588128105d+4, &
     &         a7 = 2.5090809287301226727d+3, &
     &         b1 = 4.2313330701600911252d+1, &
     &         b2 = 6.8718700749205790830d+2, &
     &         b3 = 5.3941960214247511077d+3, &
     &         b4 = 2.1213794301586595867d+4, &
     &         b5 = 3.9307895800092710610d+4, &
     &         b6 = 2.8729085735721942674d+4, &
     &         b7 = 5.2264952788528545610d+3)
    !hash sum ab    55.88319 28806 14901 4439
    !
    !coefficients for p not close to 0, 0.5 or 1.
    !
    parameter (c0 = 1.42343711074968357734d0, &
     &         c1 = 4.63033784615654529590d0, &
     &         c2 = 5.76949722146069140550d0, &
     &         c3 = 3.64784832476320460504d0, &
     &         c4 = 1.27045825245236838258d0, &
     &         c5 = 2.41780725177450611770d-1, &
     &         c6 = 2.27238449892691845833d-2, &
     &         c7 = 7.74545014278341407640d-4, &
     &         d1 = 2.05319162663775882187d0, &
     &         d2 = 1.67638483018380384940d0, &
     &         d3 = 6.89767334985100004550d-1, &
     &         d4 = 1.48103976427480074590d-1, &
     &         d5 = 1.51986665636164571966d-2, &
     &         d6 = 5.47593808499534494600d-4, &
     &         d7 = 1.05075007164441684324d-9)
    !hash sum cd    49.33206 50330 16102 89036
    !
    !coefficients for p near 0 or 1.
    !
    parameter (e0 = 6.65790464350110377720d0, &
     &         e1 = 5.46378491116411436990d0, &
     &         e2 = 1.78482653991729133580d0, &
     &         e3 = 2.96560571828504891230d-1, &
     &         e4 = 2.65321895265761230930d-2, &
     &         e5 = 1.24266094738807843860d-3, &
     &         e6 = 2.71155556874348757815d-5, &
     &         e7 = 2.01033439929228813265d-7, &
     &         f1 = 5.99832206555887937690d-1, &
     &         f2 = 1.36929880922735805310d-1, &
     &         f3 = 1.48753612908506148525d-2, &
     &         f4 = 7.86869131145613259100d-4, &
     &         f5 = 1.84631831751005468180d-5, &
     &         f6 = 1.42151175831644588870d-7, &
     &         f7 = 2.04426310338993978564d-15)
    !hash sum ef    47.52583 31754 92896 71629
    !
    if(present(ifault)) ifault = 0
    q = p - half
    if (abs(q) .le. split1) then
       r = const1 - q * q
       ppnd16 = q * (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3) &
        & * r + a2) * r + a1) * r + a0) / &
        &   (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3) &
        & * r + b2) * r + b1) * r + one)
    else
       if (q .lt. zero) then
          r = p
       else
          r = one - p
       end if
       if (r .le. zero) then
          if(present(ifault)) ifault = 1
          ppnd16 = zero
          return
       end if
       r = sqrt(-log(r))
       if (r .le. split2) then
          r = r - const2
          ppnd16 = (((((((c7 * r + c6) * r + c5) * r + c4) * r + c3) &
           & * r + c2) * r + c1) * r + c0) / &
           &  (((((((d7 * r + d6) * r + d5) * r + d4) * r + d3) &
           & * r + d2) * r + d1) * r + one)
       else
          r = r - split2
          ppnd16 = (((((((e7 * r + e6) * r + e5) * r + e4) * r + e3) &
           & * r + e2) * r + e1) * r + e0) / &
           &  (((((((f7 * r + f6) * r + f5) * r + f4) * r + f3) &
           & * r + f2) * r + f1) * r + one)
       end if
       if (q .lt. zero) ppnd16 = - ppnd16
    end if
  end function

end module math_tools
