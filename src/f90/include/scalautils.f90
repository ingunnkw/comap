module scalautils
  use healpix_types
  use quiet_mpi_mod
  use quiet_utils
  use quiet_postutils
  use scalawrap
  implicit none

contains

  !---------------------------------------------------------------------------------
  ! Read eq set from file
  !-------------------------------------------------------------------------------

  subroutine read_eqs_combo_sc(unit, info, filename, n, ordering, polarisation, amat, bmat, rhs, npix, map2mask)
    implicit none

    character(len=*),                         intent(in)  :: filename
    integer(i4b),                             intent(in)  :: unit
    type(scinfo),                             intent(in)  :: info
    integer(i4b),                             intent(out) :: ordering, polarisation, n, npix
    type(scalamat),                           intent(out) :: amat, bmat, rhs
    integer(i4b), pointer, dimension(:)                   :: map2mask
    
    real(dp)      :: t1, t2
    integer(i4b)  :: nrhs, i, j

    if (info%myid==0) write(*,*) 'Reading from ', trim(filename)

    open(unit, file = trim(filename), form="unformatted")
    read(unit) n
    if (info%myid==0) write(*,*) n,'= n'
    read(unit) ordering
    if (info%myid==0) write(*,*) ordering, '= ordering'
    read(unit) polarisation
    if (info%myid==0) write(*,*) polarisation, '= polarisation'
    call sc_alloc(amat, n, n, info)
    call sc_alloc(bmat, n, n, info)
    call sc_read_combo(unit, amat, bmat)
    read(unit) nrhs
    if (info%myid==0) write(*,*) nrhs, '= nrhs'
    call sc_alloc(rhs, n, nrhs, info)
    call sc_read(unit, rhs)
    read(unit) npix
    if (info%myid==0) write(*,*) int(sqrt(real(npix/12))), '= nside'
    allocate(map2mask(0:npix-1))
    read(unit) map2mask
    where(map2mask==0) map2mask=-1
    close(unit)
        
  end subroutine read_eqs_combo_sc
  
  subroutine read_eqs_sc(unit, info, filename, n, ordering, polarisation, invcov, rhs, npix, map2mask)
    implicit none

    character(len=*),                         intent(in)  :: filename
    integer(i4b),                             intent(in)  :: unit
    type(scinfo),                             intent(in)  :: info
    integer(i4b),                             intent(out) :: ordering, polarisation, n, npix
    type(scalamat),                           intent(out) :: invcov, rhs
    integer(i4b), pointer, dimension(:)                   :: map2mask
    
    real(dp)      :: t1, t2
    integer(i4b)  :: nrhs

    if (info%myid==0) write(*,*) 'Reading from ', trim(filename)

    !    call wall_time(t1)

    open(unit, file = trim(filename), form="unformatted")
    read(unit) n
    if (info%myid==0) write(*,*) n,'= n'
    read(unit) ordering
    if (info%myid==0) write(*,*) ordering, '= ordering'
    read(unit) polarisation
    if (info%myid==0) write(*,*) polarisation, '= polarisation'
    call sc_alloc(invcov, n, n, info)
    call sc_read(unit, invcov)!, broadcast=.true.)
!    call sc_alloc(rhs, n, 1, info)
!    call sc_read(unit, rhs)!, broadcast=.true.)
    read(unit) nrhs
    if (info%myid==0) write(*,*) nrhs, '= nrhs'
    call sc_alloc(rhs, n, nrhs, info)
    call sc_read(unit, rhs)
    read(unit) npix
    if (info%myid==0) write(*,*) int(sqrt(real(npix/12))), '= nside'
    allocate(map2mask(0:npix-1))
    read(unit) map2mask
    close(unit)

    !    call mpi_barrier(info%comm)
    !    call wall_time(t2)
    !    if (invcov%info%myid==0) write(*,*) 'Time spent reading from file:', t2-t1 
        
  end subroutine read_eqs_sc

  ! Eats its inputs like anything else
  subroutine write_eqs_sc(unit, info, filename, n, ordering, polarisation, invcov, rhs, npix, map2mask)
    implicit none

    character(len=*),                         intent(in)    :: filename
    integer(i4b),                             intent(in)    :: unit
    type(scinfo),                             intent(in)    :: info
    integer(i4b),                             intent(in)    :: ordering, polarisation, n, npix
    type(scalamat),                           intent(inout) :: invcov, rhs
    integer(i4b), dimension(:),               intent(in)    :: map2mask
    real(dp) :: t1, t2

    open(unit, file = trim(filename), form="unformatted")
    if(info%myid == 0) then
       write(unit) n
       write(unit) ordering
       write(unit) polarisation
    end if
    call sc_write(unit, invcov)
    call sc_write(unit, rhs)
    if(info%myid == 0) then
       write(unit) npix
       write(unit) map2mask
    end if
    close(unit)

    call sc_dealloc(invcov)
    call sc_dealloc(rhs)
  end subroutine write_eqs_sc

  !---------------------------------------------------------------------------------
  ! Read cov matrix from file
  !-------------------------------------------------------------------------------
  subroutine read_cov_sc(unit, info, filename, n, ordering, polarisation, invcov, &
    & isinv)
    implicit none

    character(len=*),                         intent(in)  :: filename
    integer(i4b),                             intent(in)  :: unit
    type(scinfo),                             intent(in)  :: info
    integer(i4b),                             intent(out) :: ordering, polarisation, n
    type(scalamat),                           intent(out) :: invcov
    logical(lgt),                    optional,intent(out) :: isinv

    
    real(dp) :: t1, t2

    open(unit, file = trim(filename), form="unformatted")
    read(unit) n
    if (info%myid==0) write(*,*) n,'= n'
    read(unit) ordering
    if (info%myid==0) write(*,*) ordering, '= ordering'
    read(unit) polarisation
    if (info%myid==0) write(*,*) polarisation, '= polarisation'
    call sc_alloc(invcov, n, n, info)
    call sc_read(unit, invcov)!, broadcast=.true.)
    if (present(isinv)) read(unit) isinv
    close(unit)
  end subroutine read_cov_sc

  !---------------------------------------------------------------------------------
  ! Read eigendec of cov matrix from file
  !-------------------------------------------------------------------------------
  subroutine read_eig_sc(unit, info, filename, n, ordering, polarisation, eigmatrix, eigvals, isinv)
    implicit none

    character(len=*),                         intent(in)  :: filename
    integer(i4b),                             intent(in)  :: unit
    type(scinfo),                             intent(in)  :: info
    integer(i4b),                             intent(out) :: ordering, polarisation, n
    type(scalamat),                           intent(out) :: eigmatrix, eigvals
    logical(lgt),                    optional,intent(out) :: isinv

    open(unit, file = trim(filename), form="unformatted")
    read(unit) n
    if (info%myid==0) write(*,*) n,'= n'
    read(unit) ordering
    if (info%myid==0) write(*,*) ordering, '= ordering'
    read(unit) polarisation
    if (info%myid==0) write(*,*) polarisation, '= polarisation'
    call sc_alloc(eigmatrix, n, n, info)
    call sc_read(unit, eigmatrix)
    call sc_alloc(eigvals, n, 1, info)
    call sc_read(unit, eigvals)
    if (present(isinv)) read(unit) isinv
    close(unit)
  end subroutine

  !---------------------------------------------------------------------------------
  ! Read general matrix from file
  !-------------------------------------------------------------------------------
  subroutine read_genmat_sc(unit, info, filename, matrix)
    implicit none

    character(len=*),  intent(in)    :: filename
    integer(i4b),      intent(in)    :: unit
    type(scinfo),      intent(in)    :: info
    type(scalamat),    intent(inout) :: matrix

    integer(i4b) :: m, n
    real(dp) :: t1, t2

    open(unit, file = trim(filename), form="unformatted")
    read(unit) n
    read(unit) m
    if (info%myid==0) write(*,*) n,'= n', m, ' = m'
    call sc_dealloc(matrix)
    call sc_alloc(matrix, n, m, info)
    call sc_read(unit, matrix)
    close(unit)
  end subroutine read_genmat_sc

  !---------------------------------------------------------------------------------
  ! Write general matrix from file
  !-------------------------------------------------------------------------------
  subroutine write_genmat_sc(unit, info, filename, matrix)
    implicit none

    character(len=*),  intent(in)    :: filename
    integer(i4b),      intent(in)    :: unit
    type(scinfo),      intent(in)    :: info
    type(scalamat),    intent(in)    :: matrix

    integer(i4b) :: m, n
    real(dp) :: t1, t2

    n = matrix%rglob
    m = matrix%cglob

    if (info%myid==0) then
       open(unit, file=trim(filename), form='unformatted')
       write(unit) n
       write(unit) m
    end if
    call sc_write(unit,matrix)
    if (info%myid==0) then
       close(unit)
    end if

  end subroutine write_genmat_sc

  !---------------------------------------------------------------------------------
  ! Write covmatrix to file
  !-------------------------------------------------------------------------------

  subroutine write_covmatrix_sc(unit, info, filename, ordering, pol, matrix, inv)
    integer(i4b),                       intent(in) :: unit, ordering, pol
    type(scinfo),                       intent(in) :: info
    character(len=*),                   intent(in) :: filename
    logical(lgt),                       intent(in) :: inv
    type(scalamat),                     intent(in) :: matrix

    integer(i4b) :: n, m, i, ierr

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    m = matrix%rglob
    n = matrix%cglob

    if (n/= m) then
       write(*,*) ' cov matrix not n x n matrix. quiting'
       call mpi_finalize(ierr)       
       stop
    end if

    if (info%myid==0) then
       open(unit, file=trim(filename), form='unformatted')
       write(unit) n
       write(unit) ordering
       write(unit) pol
    end if
    call sc_write(unit,matrix)
    if (info%myid==0) then
       write(unit) inv
       close(unit)
    end if

  end subroutine write_covmatrix_sc

  !---------------------------------------------------------------------------------
  ! Write eigendec of covmatrix to file
  !-------------------------------------------------------------------------------

  subroutine write_eig_sc(unit, info, filename, ordering, pol, eigmatrix, eigvals, inv)
    integer(i4b),                       intent(in) :: unit, ordering, pol
    type(scinfo),                       intent(in) :: info
    character(len=*),                   intent(in) :: filename
    logical(lgt),                       intent(in) :: inv
    type(scalamat),                     intent(in) :: eigmatrix, eigvals

    integer(i4b) :: n, m, p, i, ierr

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    m = eigmatrix%rglob
    n = eigmatrix%cglob
    p = eigvals%rglob
    i = eigvals%cglob
    if (n/= m) then
       write(*,*) n, m,' cov matrix not n x n matrix. quiting'
       call mpi_finalize(ierr)       
       stop
    end if
    if (p/= m) then
       write(*,*) m,p,' eigmatrix and eigvals not same size. quiting'
       call mpi_finalize(ierr)       
       stop
    end if
    if (i /= 1) then
       write(*,*) i, ' eigvals is not a vector. quiting'
       call mpi_finalize(ierr)       
       stop
    end if

    if (info%myid==0) then
       open(unit, file=trim(filename), form='unformatted')
       write(unit) n
       write(unit) ordering
       write(unit) pol
    end if
    call sc_write(unit, eigmatrix)
    call sc_write(unit, eigvals)
    if (info%myid==0) then
       write(unit) inv
       close(unit)
    end if

  end subroutine

  !---------------------------------------------------------------------------------
  ! Multiply a vector from both sides with matrix of eigenvectors
  !-------------------------------------------------------------------------------

  subroutine eigmult_sc(eigmatrix, vector, out)
    implicit none

    type(scalamat),           intent(in)  :: eigmatrix, vector
    type(scalamat),           intent(out) :: out

    type(scinfo)                          :: info
    type(scalamat)                        :: matrix
    integer(i4b)                          :: i, j, m, n, ierr
    real(dp)                              :: t1, t2

    info = eigmatrix%info
    m = eigmatrix%rglob
    n = eigmatrix%cglob

    if (n/= m) then
       write(*,*) 'eigenvector matrix not n x n matrix. quiting'
       call mpi_finalize(ierr)       
       stop
    end if

    ! out = eigmatrix * diag(vector) * (eigmatrix)^T
    call sc_alloc(out,      n, n, info)    
    call sc_alloc(matrix,   n, n, info)
do j=1,matrix%cloc
   do i = 1, matrix%rloc
      matrix%data(i,j) = eigmatrix%data(i,j)
   end do
end do
    call sc_matmul_diag(matrix, vector)   
    call cpu_time(t1) 
    call sc_matmul(matrix, eigmatrix, out,    transb='t')
    call cpu_time(t2)
    if (info%myid==0) write(*,*) 'Time spent on matrix multiplication:', t2-t1 

    call sc_dealloc(matrix)
  end subroutine eigmult_sc

  !---------------------------------------------------------------------------------
  ! Solve eq set: Find map and cov from invcov and rhs
  !-------------------------------------------------------------------------------

  subroutine solve_eqs_sc(info, unit, n, invcov, rhs, map, cov, bmat, outprefix, outtext, kept)
    implicit none

    type(scinfo),                          intent(in)    :: info
    character(len=*), optional,            intent(in)    :: outprefix
    character(len=*), optional,            intent(inout) :: outtext
    integer(i4b),                          intent(in)    :: unit, n
    type(scalamat),                        intent(inout) :: invcov   ! destroyed
    type(scalamat),                        intent(inout) :: rhs      ! destroyed
    type(scalamat),                        intent(out)   :: map, cov 
    type(scalamat),   optional                           :: bmat     ! destroyed
    integer(i4b),     optional,            intent(inout) :: kept


    type(scalamat)                        :: eig
    integer(i4b)                          :: i,j, nrhs
    real(dp)                              :: maxeigenval, t1, t2
    real(dp), allocatable, dimension(:,:) :: eigenvals, mat

    nrhs = rhs%cglob

    ! decompose inverse covariance matrix into eigenvectors
    call cpu_time(t1)
    call sc_alloc(eig, n, 1, info)
    call sc_eigenvalue_decompose(invcov,eig) !eigenvectors put in invcov, eigenvalues in eig
    call cpu_time(t2)
    if (info%myid==0) then
       write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1 
       allocate(eigenvals(n,1))
    end if

    ! Getting eigenvalues for invcov into eigenvals
    call sc_get(eig, eigenvals, 0)
    ! check and invert eigenvalues
    outtext = 'inveigenval_' // outtext
    if (info%myid==0) call invert_eigenvals(info%myid, eigenvals(:,1), unit, outprefix, outtext, kept=kept)
    ! Putting eigenvalues for noninverse cov into eig
    call sc_set(eig, eigenvals, 0)
    if (info%myid==0) deallocate(eigenvals)

    ! finding map = matmul(eigenvectors, matmul(transpose(eigenvectors), rhs1)*eigenvals)
    call sc_alloc(map, n, nrhs, info)
    call sc_matmul(invcov, rhs, map,    transa='t')
    !rhs%data=map%data*eig%data     ! now destroying rhs
    call sc_matmul_diag(map, eig, rhs, flip=.true.) ! now destroying rhs
    call sc_matmul(invcov, rhs, map)

    ! finding non-inverse covariance matrix
    call eigmult_sc(invcov, eig, cov)

    ! find cov for full mapmaking eqs
    if (present(bmat)) then
       call sc_matmul(cov, bmat, invcov)
       call sc_matmul(invcov, cov, bmat)
       do j=1, cov%cloc
          do i = 1, cov%rloc
             cov%data(i,j) = bmat%data(i,j)
          end do
       end do
       call sc_dealloc(bmat)
    end if

    ! Clean up
    call sc_dealloc(eig)
    call sc_dealloc(invcov)
    call sc_dealloc(rhs)

  end subroutine solve_eqs_sc

  !---------------------------------------------------------------------------------
  ! Solve eq set: Find map and (cov and/or sqrt of invcov/cov) from invcov and rhs
  !-------------------------------------------------------------------------------

  subroutine solve_eqs_withsqrt_sc(invcov, rhs, map, cov, sqrt, sqrt_input, unit, outprefix, outtext)
    implicit none

    character(len=*), optional,            intent(in)    :: outprefix
    character(len=*), optional,            intent(inout) :: outtext
    integer(i4b),     optional,            intent(in)    :: unit
    logical(lgt),     optional,            intent(in)    :: sqrt_input
    type(scalamat),                        intent(inout) :: invcov, rhs ! destroyed
    type(scalamat),                        intent(out)   :: map 
    type(scalamat),   optional,            intent(out)   :: sqrt, cov 

    type(scalamat)                        :: eig, sqrteig, matrix
    type(scinfo)                          :: info
    integer(i4b)                          :: i,j,n
    real(dp)                              :: maxeigenval, t1, t2
    real(dp), allocatable, dimension(:,:) :: eigenvals, sqrteigenvals

    n = invcov%rglob
    info = invcov%info

    ! decompose inverse covariance matrix into eigenvectors
    call cpu_time(t1)
    call sc_alloc(eig, n, 1, info)
    call sc_alloc(sqrteig, n, 1, info)
    call sc_eigenvalue_decompose(invcov,eig) !eigenvectors put in invcov, eigenvalues in eig
    call cpu_time(t2)
    if (info%myid==0) then
       write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1 
       allocate(eigenvals(n,1))
       allocate(sqrteigenvals(n,1))
    end if

    ! Getting eigenvalues for invcov into eigenvals
    call sc_get(eig, eigenvals, 0)
    ! check and invert eigenvalues
    outtext = 'inveigenval_' // outtext
    if (info%myid==0) then
       call invert_andsqrt_eigenvals(info%myid, eigenvals(:,1), sqrteigenvals(:,1), sqrt_input, unit, outprefix, outtext)
    end if
    ! Putting eigenvalues for noninverse cov back into eig
    call sc_set(eig, eigenvals, 0)
    ! Putting sqrt of eigenvalues into sqrteig
    call sc_set(sqrteig, sqrteigenvals, 0)
    if (info%myid==0) then
       deallocate(eigenvals)
       deallocate(sqrteigenvals)
    end if

    ! finding map = matmul(eigenvectors, matmul(transpose(eigenvectors), rhs1)*eigenvals)
    call sc_alloc(map, n, 1, info)
    call sc_matmul(invcov, rhs, map,    transa='t')
    rhs%data=map%data*eig%data     ! now destroying rhs
    call sc_matmul(invcov, rhs, map)

    ! finding sqrt of invcov
    if (present(sqrt)) call eigmult_sc(invcov, sqrteig, sqrt)

    ! finding non-inverse covariance matrix
    if (present(cov)) call eigmult_sc(invcov, eig, cov)

    ! Clean up
    call sc_dealloc(eig)
    call sc_dealloc(sqrteig)
    call sc_dealloc(invcov)
    call sc_dealloc(rhs)

  end subroutine solve_eqs_withsqrt_sc


  !---------------------------------------------------------------------------------
  ! Invert matrix and find sqrt of either input or output matrix or both
  !-------------------------------------------------------------------------------

  subroutine sqrtandinv_matrix_sc(unit, outprefix, outtext, inmatrix, outmatrix, sqrtmatrix, sqrt_input, matrix3)
    implicit none

    character(len=*)                                     :: outtext
    character(len=*),                      intent(in)    :: outprefix
    integer(i4b),                          intent(in)    :: unit
    type(scalamat),                        intent(inout) :: inmatrix ! destroyed
    type(scalamat),                        intent(out)   :: outmatrix, sqrtmatrix 
    type(scalamat),              optional, intent(out)   :: matrix3
    logical(lgt),                          intent(in)    :: sqrt_input

    type(scalamat)                        :: eig
    type(scinfo)                          :: info
    integer(i4b)                          :: i,j,n
    real(dp)                              :: maxeigenval, t1, t2
    real(dp), allocatable, dimension(:,:) :: eigenvals, sqrteigenvals, sqrteigenvals2

    n = inmatrix%rglob
    info = inmatrix%info

    ! decompose input matrix into eigenvectors
    call cpu_time(t1)
    call sc_alloc(eig, n, 1, info)
    call sc_eigenvalue_decompose(inmatrix,eig) !eigenvectors put in inmatrix, eigenvalues in eig
    call cpu_time(t2)
    if (info%myid==0) then
       write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1 
       allocate(eigenvals(n,1))
       allocate(sqrteigenvals(n,1))
       if (present(matrix3)) allocate(sqrteigenvals2(n,1))
    end if
    ! Getting eigenvalues for inmatrix into eigenvals
    call sc_get(eig, eigenvals, 0)

    ! check and invert and get sqrt of eigenvalues
    if (info%myid==0) then
       if (present(matrix3)) then
          call invert_andsqrt_eigenvals(info%myid, eigenvals(:,1), sqrteigenvals(:,1), sqrt_input, unit, outtext=outtext, sqrteigenvals2=sqrteigenvals2(:,1))
       else
          call invert_andsqrt_eigenvals(info%myid, eigenvals(:,1), sqrteigenvals(:,1), sqrt_input, unit, outtext=outtext)
       end if
    end if

    ! finding inverse of input matrix
    call sc_set(eig, eigenvals, 0)
    call eigmult_sc(inmatrix, eig, outmatrix)
    if (info%myid==0) deallocate(eigenvals)

    ! finding sqrt of matrix
    call sc_set(eig, sqrteigenvals, 0)
    call eigmult_sc(inmatrix, eig, sqrtmatrix)
    if (info%myid==0) deallocate(sqrteigenvals)

    ! finding sqrt of other matrix
    if (present(matrix3)) then
       call sc_set(eig, sqrteigenvals2, 0)
       call eigmult_sc(inmatrix, eig, matrix3)
       if (info%myid==0) deallocate(sqrteigenvals2)
    end if

    ! Clean up
    call sc_dealloc(eig)
    call sc_dealloc(inmatrix)

  end subroutine 

  !---------------------------------------------------------------------------------
  ! Find sqrt of a matrix
  !-------------------------------------------------------------------------------

  subroutine sqrt_matrix_sc(unit, outtext, invcov, sqrtinvcov, sqrtcov)
    implicit none

    character(len=*), intent(inout)         :: outtext
    integer(i4b),     intent(in)            :: unit
    type(scalamat),   intent(inout)         :: invcov ! destroyed
    type(scalamat),   intent(out)           :: sqrtinvcov
    type(scalamat),   intent(out), optional :: sqrtcov                     

    type(scalamat)                        :: eig
    type(scinfo)                          :: info
    integer(i4b)                          :: i,j,n
    real(dp)                              :: maxeigenval, t1, t2
    real(dp), allocatable, dimension(:,:) :: eigenvals

    n = invcov%rglob
    info = invcov%info

    ! decompose inverse covariance matrix into eigenvectors
    call cpu_time(t1)
    call sc_alloc(eig, n, 1, info)
    call sc_eigenvalue_decompose(invcov,eig) !eigenvectors put in invcov, eigenvalues in eig
    call cpu_time(t2)
    if (info%myid==0) then
       write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1 
       allocate(eigenvals(n,1))
    end if

    ! Getting eigenvalues for invcov into eigenvals
    call sc_get(eig, eigenvals, 0)
    ! check and invert eigenvalues
    outtext = 'inveigenval_' // outtext
    if (info%myid==0) then
       call sqrt_eigenvals(info%myid, outtext, eigenvals(:,1))
    end if
    ! Putting sqrt of eigenvalues back into eig
    call sc_set(eig, eigenvals, 0)

    ! finding sqrt of inverse covariance matrix
    call eigmult_sc(invcov, eig, sqrtinvcov)

    ! Output inverse of sqrt as well, if requested
    if (present(sqrtcov)) then
       ! check and invert eigenvalues
       outtext = 'eigenval_' // outtext
       if (info%myid==0) then
          where (eigenvals > 0.d0)
             eigenvals = 1.d0 / eigenvals
          end where
       end if
       ! Putting sqrt of eigenvalues back into eig
       call sc_set(eig, eigenvals, 0)
       
       ! finding sqrt of inverse covariance matrix
       call eigmult_sc(invcov, eig, sqrtcov)       
    end if

    ! Clean up
    if (info%myid==0) deallocate(eigenvals)
    call sc_dealloc(eig)
    call sc_dealloc(invcov)

  end subroutine sqrt_matrix_sc


  !---------------------------------------------------------------------
  ! Removing large and uncommon pixels from map and cov using in2red
  !----------------------------------------------------------------------

  subroutine inmap2red_sc(info, in2red, inmap, incov, redmap, redcov, pol, inmat, redmat)
   implicit none

    type(scinfo),                   intent(in)    :: info
    integer(i4b),    dimension(1:), intent(in)    :: in2red
    type(scalamat),                 intent(inout) :: inmap, incov   ! destroyed
    type(scalamat),                 intent(out)   :: redmap, redcov ! allocated
    integer(i4b),                   intent(in)    :: pol
    type(scalamat),  optional,      intent(inout) :: inmat  ! destroyed
    type(scalamat),  optional,      intent(out)   :: redmat ! alloacted

    integer(i4b)                             :: n, nmaps
    integer(i4b), allocatable, dimension(:)  :: polvec

    n = size(in2red)
    nmaps = inmap%cglob

    ! Removing bad pixels from matrices
    allocate(polvec(pol*n))
    polvec(1:n)     = in2red
    if (pol == 2) polvec(n+1:2*n) = in2red + (inmap%rglob)/2
    call sc_alloc(redmap, pol*n, nmaps,   info)      
    call sc_alloc(redcov, pol*n, pol*n, info)
    call sc_copy(inmap, redmap, rows=polvec)
    call sc_copy(incov, redcov, rows=polvec, cols=polvec)
    call sc_dealloc(inmap)
    call sc_dealloc(incov)
    if (present(redmat)) then 
       call sc_alloc(redmat, pol*n, pol*n, info)
       call sc_copy(inmat, redmat, rows=polvec, cols=polvec)
       call sc_dealloc(inmat)
    end if
    deallocate(polvec)

  end subroutine inmap2red_sc

  !---------------------------------------------------------------------
  ! Calculating chi_square = mapdiff * invcov * mapdiff
  !----------------------------------------------------------------------

  subroutine chi_square_sc(map, invcov, chi_square)
   implicit none

    type(scalamat),            intent(in)    :: map, invcov
    real(dp),  dimension(1:),  intent(out)   :: chi_square 

    integer(i4b)                             :: n, nrhs, i, j, rect(4)
    type(scalamat)                           :: vector
    real(dp),  dimension(:,:), allocatable   :: loc

    n    = map%rglob
    nrhs = map%cglob

    call sc_alloc(vector, n, nrhs, map%info)
    call sc_matmul(invcov, map, vector)
!    chi_square = sc_dotprod(map, vector)
    do j = 1, map%cloc
       do i = 1, map%rloc
          vector%data(i,j) =  vector%data(i,j)*map%data(i,j)
       end do
    end do
    if (map%info%myid == 0) allocate(loc(n,1))
    do i = 1, nrhs
       rect = (/1,i, n,i/)
       call sc_get(vector, loc, 0, rect)
       if (map%info%myid == 0) chi_square(i) = sum(loc)
    end do

    call sc_dealloc(vector) 
    if (map%info%myid == 0) deallocate(loc)
    
  end subroutine chi_square_sc
  !---------------------------------------------------------------------
  ! Calculating chi_square = mapdiff * invcov * mapdiff per eigenmode
  !----------------------------------------------------------------------

  subroutine chi_square_eigen_sc(map, eig, eigvals, chi_square, outprefix)
   implicit none

    type(scalamat),            intent(in)    :: map, eig, eigvals
    real(dp),  dimension(1:),  intent(out)   :: chi_square
    character(len=*),optional, intent(in)    :: outprefix 

    integer(i4b)                             :: n, nrhs, i, j, rect(4)
    type(scalamat)                           :: vector
    real(dp),  dimension(:,:), allocatable   :: loc

    n    = map%rglob
    nrhs = map%cglob

    call sc_alloc(vector, n, nrhs, map%info)
    call sc_matmul(eig, map, vector, transa='t')
    do j = 1, map%cloc   
       do i = 1, map%rloc
          vector%data(i,j) =  vector%data(i,j)**2 *eigvals%data(i,1)
       end do
    end do
    if (map%info%myid == 0) allocate(loc(n,1))
    do i = 1, nrhs
       rect = (/1,i, n,i/)
       call sc_get(vector, loc, 0, rect)
       if (map%info%myid == 0) then
          chi_square(i) = sum(loc)
          if (i ==1 .and. present(outprefix)) then
             open(15, file=trim(outprefix) // '_chisq_per_eigenmode.dat')
             do j = 1, n
                write(15,*) loc(j,1)
             end do
             close(15)
          end if
       end if
    end do

    call sc_dealloc(vector) 
    if (map%info%myid == 0) deallocate(loc)
    
  end subroutine

  !---------------------------------------------------------------------
  ! Calculating chi_square = mapdiff * invcov * mapdiff
  !----------------------------------------------------------------------

  subroutine chi_square_single_sc(map, invcov, chi_square)
   implicit none

    type(scalamat),            intent(in)    :: map, invcov
    real(dp),                  intent(out)   :: chi_square 

    integer(i4b)                             :: n, nrhs
    type(scalamat)                           :: vector

    n    = map%rglob
    nrhs = 1

    call sc_alloc(vector, n, nrhs, map%info)
    call sc_matmul(invcov, map, vector)
    chi_square = sc_dotprod(map, vector)
    call sc_dealloc(vector) 
        
  end subroutine chi_square_single_sc

  !---------------------------------------------------------------------
  ! Finding summap, diffmap and invcov for weighted maps
  !----------------------------------------------------------------------

  subroutine weighted_maps_sc(info, map1, map2, invcov1, invcov2, cov, mapsum, mapdiff)
   implicit none

    type(scinfo),   intent(in)    :: info
    type(scalamat), intent(inout) :: map1, map2, invcov1, invcov2 !is destroyed
    type(scalamat), intent(in)    :: cov
    type(scalamat), intent(out)   :: mapsum, mapdiff

    type(scalamat)                :: vec1, vec2
    integer(i4b)                  :: n, i, j, nrhs 

    n    = map1%rglob
    nrhs = map1%cglob
    
    call sc_alloc(mapdiff, n, nrhs, info)
    call sc_alloc(mapsum,  n, nrhs, info)
    call sc_alloc(vec1,    n, nrhs, info)
    call sc_alloc(vec2,    n, nrhs, info)

    ! Computing diff_map and sum_map
    call sc_matmul(invcov1, map1, vec1)
    call sc_matmul(invcov2, map2, vec2)

do j=1,map2%cloc
  do i = 1, map1%rloc
    map1%data(i,j) = vec1%data(i,j) + vec2%data(i,j)
    map2%data(i,j) = vec1%data(i,j) - vec2%data(i,j)
  end do
end do
    call sc_matmul(cov, map1, mapsum)
    call sc_matmul(cov, map2, mapdiff)
    
    ! Clean up
    call sc_dealloc(vec1)
    call sc_dealloc(vec2)
    call sc_dealloc(map1)
    call sc_dealloc(map2)
    call sc_dealloc(invcov1)
    call sc_dealloc(invcov2)

  end subroutine weighted_maps_sc

  !---------------------------------------------------------------------------------
  ! Invert matrix using eigendecomp
  !-------------------------------------------------------------------------------

  subroutine invert_matrix_eigen_sc(A, unit, outprefix, outtext, kept, acc)
    implicit none

    type(scalamat),                  intent(inout) :: A ! destroyed
    integer(i4b),          optional, intent(in)    :: unit   
    character(len=*),      optional, intent(in)    :: outtext, outprefix
    integer(i4b),          optional, intent(out)   :: kept ! # of kept eigenvalues   
    real(dp),              optional, intent(in)    :: acc  ! threshold for eigenvalue removal 

    type(scalamat)                        :: eig, matrix
    type(scinfo)                          :: info
    integer(i4b)                          :: i, j, n, ierr
    real(dp)                              :: maxeigenval, t1, t2
    real(dp), allocatable, dimension(:,:) :: eigenvals
    
    n=A%rglob
    info=A%info

    ! decompose A into eigenvectors
    call cpu_time(t1)
    call sc_alloc(eig, n, 1, A%info)
    call sc_eigenvalue_decompose(A,eig) !eigenvectors put in A, eigenvalues in eig
    call cpu_time(t2)
    if (A%info%myid==0) then
       write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1 
       allocate(eigenvals(n,1))
    end if
    ! Getting eigenvalues for A into eigenvals
    call sc_get(eig, eigenvals, 0)
    ! check and invert eigenvalues
    if (A%info%myid==0) call invert_eigenvals(info%myid, eigenvals(:,1), unit, outprefix, outtext, kept, acc)

    ! Putting eigenvalues for inverse A back into eig
    call sc_set(eig, eigenvals, 0)
    if (A%info%myid==0) then
       deallocate(eigenvals)
    end if
    ! finding inverse of A
    call eigmult_sc(A, eig, matrix)
    ! And putting it back into A
do j=1,A%cloc
   do i = 1, A%rloc
    A%data(i,j)=matrix%data(i,j)
 end do
end do

    ! Clean up
    call sc_dealloc(eig)
    call sc_dealloc(matrix)

  end subroutine invert_matrix_eigen_sc

  !---------------------------------------------------------------------------------
  ! Invert or sqrt matrix from eigendecomp
  !-------------------------------------------------------------------------------

  subroutine matrix_from_eigen_sc(eig, eigvals, unit, mat, invmat, sqrtmat, sqrtinv, invmatvals, sqrtmatvals, sqrtinvvals, outprefix, outtext, kept, acc, cov, ignore)
    implicit none

    type(scalamat)                                 :: eig, eigvals
    type(scalamat),        optional                :: mat, invmat, sqrtmat, sqrtinv
    type(scalamat),        optional                :: invmatvals, sqrtmatvals, sqrtinvvals
    integer(i4b),          optional, intent(in)    :: unit   
    character(len=*),      optional, intent(in)    :: outtext, outprefix
    integer(i4b),          optional, intent(out)   :: kept ! # of kept eigenvalues   
    real(dp),              optional, intent(in)    :: acc  ! threshold for eigenvalue removal 
    logical(lgt),          optional, intent(in)    :: cov
    integer(i4b),          optional, intent(in)    :: ignore ! # of previously cut eigenvalues   

    type(scinfo)                          :: info
    integer(i4b)                          :: i, j, n
    real(dp), allocatable, dimension(:,:) :: eigenvals, invei, sqrtei, sqrtinvei

    n=eig%rglob
    info=eig%info

    ! Getting eigenvalues into local root eigenvals
    if (info%myid==0) then 
       allocate(eigenvals(n,1))
       allocate(invei(n,1))
       allocate(sqrtei(n,1))
       allocate(sqrtinvei(n,1))
    end if
    call sc_get(eigvals, eigenvals, 0)
    ! check and invert eigenvalues
    if (info%myid==0) call modify_eigenvals(info%myid, eigenvals(:,1), invei(:,1), sqrtei(:,1), sqrtinvei(:,1), unit, outprefix, outtext, kept, acc, cov, ignore)


    ! Output invmat
    if (present(invmat)) then
       call sc_set(eigvals, invei, 0)
       call eigmult_sc(eig, eigvals, invmat)
       if (present(invmatvals)) call sc_copy(eigvals, invmatvals)
    end if
    ! Output sqrtmat
    if (present(sqrtmat)) then
       call sc_set(eigvals, sqrtei, 0)
       call eigmult_sc(eig, eigvals, sqrtmat)
       if (present(sqrtmatvals)) call sc_copy(eigvals, sqrtmatvals)
    end if
    ! Output sqrtinv
    if (present(sqrtinv)) then
       call sc_set(eigvals, sqrtinvei, 0)
       call eigmult_sc(eig, eigvals, sqrtinv)
       if (present(sqrtinvvals)) call sc_copy(eigvals, sqrtinvvals)
    end if
    ! Put original (but checked and cut) eigenvals back into eigvals
    call sc_set(eigvals, eigenvals, 0)
    ! Output mat
    if (present(mat)) then
       call eigmult_sc(eig, eigvals, mat)
    end if

    ! Clean up
    if (info%myid==0) deallocate(eigenvals, invei, sqrtei, sqrtinvei)

  end subroutine matrix_from_eigen_sc

  !---------------------------------------------------------------------------------
  ! Solve eq set fast: Find map and cov from invcov and rhs
  !-------------------------------------------------------------------------------

  subroutine fast_solve_eqs_sc(n, invcov, rhs, map, cov, bmat)
    implicit none

    integer(i4b),                intent(in)    :: n
    type(scalamat),              intent(inout) :: invcov, rhs ! destroyed
    type(scalamat),              intent(out)   :: map, cov 
    type(scalamat),   optional,  intent(inout) :: bmat ! destroyed

    integer(i4b)    :: i, j, nrhs
    real(dp)        :: t1, t2

    nrhs=rhs%cglob

    call sc_alloc(cov, n, n, invcov%info)    
    call sc_alloc(map, n, nrhs, invcov%info)
    call cpu_time(t1)
    do j=1, cov%cloc
       do i = 1, cov%rloc
          cov%data(i,j) = invcov%data(i,j)
       end do
    end do
    call cpu_time(t2)
    if (invcov%info%myid==0) write(*,*) 'Time spent on matrix = matrix:', t2-t1 
    call sc_invert(cov)
    call cpu_time(t1)
    if (invcov%info%myid==0) write(*,*) 'Time spent on matrix invertion:', t1-t2 
    call sc_matmul(cov, rhs, map)
    call cpu_time(t2)
    if (invcov%info%myid==0) write(*,*) 'Time spent on matrix multiplication:', t2-t1 

    ! find cov for full mapmaking eqs
    if (present(bmat)) then
       call sc_matmul(cov, bmat, invcov)
       call sc_matmul(invcov, cov, bmat)
       do j=1, cov%cloc
          do i = 1, cov%rloc
             cov%data(i,j) = bmat%data(i,j)
          end do
       end do
       call sc_dealloc(bmat)
    end if

    ! Clean up
    call sc_dealloc(invcov)
    call sc_dealloc(rhs)

  end subroutine fast_solve_eqs_sc

  ! Performs the matrix vector product res = mat*vec, where only mat
  ! is a distributed matrix, and vec and res are normal arrays.
  ! The input vec must be present on all processors, and the result is
  ! distributed to everyone.
  subroutine mul_smat_vec(mat, vec, res)
    implicit none
    type(scalamat) :: mat, svec, sres
    real(dp)       :: vec(:), res(:)
    real(dp), dimension(:,:), allocatable :: tmp
    integer(i4b)   :: n, m, i
    n = mat%cglob; m = mat%rglob
    call sc_alloc(svec, n, 1, mat%info)
    call sc_alloc(sres, m, 1, mat%info)
    call sc_allset(svec,reshape(vec,(/n,1/)))
    call sc_matmul(mat, svec, sres)
    allocate(tmp(m,1))
    call sc_allget(sres,tmp)
    res = tmp(:,1)
    deallocate(tmp)
    call sc_dealloc(svec)
    call sc_dealloc(sres)
  end subroutine

end module scalautils
