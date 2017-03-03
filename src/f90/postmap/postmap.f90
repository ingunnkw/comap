program postmap
  use healpix_types
  use pix_tools
  use fitstools
  use quiet_fileutils
  use math_tools
  use quiet_postutils
  use quiet_utils
  implicit none

  include "mpif.h"

  integer(i4b)       :: iargc
  integer(i4b)       :: i, j, k
  integer(i4b)       :: unit, myid, numprocs, ierr, root
  character(len=30)  :: kommando

  type pixel_set
     integer(i4b), dimension(:), allocatable :: pixels
  end type

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0
  unit        = 21+myid

  if (iargc() == 0) then
     call give_user_info
  else if (myid == root) then
     write(*,*) '-------------------- QUIET map post-processor  --------------------'
     write(*,*)
  end if

  ! Get name of main kommando
  call getarg(1,kommando)

  if (kommando == 'matmul') then
     write(*,*) 'Matrix multiplication'
     call do_matmul(unit)
  else if (kommando == 'proj0') then
     write(*,*) 'Project out null-space'
     call project_out_nullspace(unit)
  else if (kommando == 'finalmap') then
     write(*,*) 'Make final map and diff_map and covariance matrix from two subsets'
     call make_finalmap(unit)
  else if (kommando == 'solve') then
     write(*,*) 'Solve an equation system in the format produced by tod2map'
     call solve_equation(unit)
  else if (kommando == 'cutoff' .or. kommando == 'tempcutoff' .or. kommando == 'tempcut') then
     write(*,*) 'Project out eigenmodes lower than cutoff for temperature data'
     call project_out_eigenmodes(unit, .true.)
  else if (kommando == 'polcutoff' .or. kommando== 'polcut') then
     write(*,*) 'Project out eigenmodes lower than cutoff for polarisation data'
     call project_out_eigenmodes(unit, .false.)
  else if (kommando == 'lowcut' .or. kommando == 'lcut') then
     write(*,*) 'Project out modes lower or equal to l_max'
     call project_out_lowmodes(unit)
  else if (kommando == 'rms') then
     write(*,*) 'Produce an rms map from a covariance matrix and map2mask'
     call output_rms(unit)
  else if (kommando == 'eqn_sum') then
     write(*,*) 'Sum equation sets with possibly different pixel sets'
     call do_sum_equations
  else if (kommando == 'kernel') then
     write(*,*) 'Find coupling kernel for given mask'
     call coupling_kernel(unit)
  else if (kommando == 'eqn_compare') then
     write(*,*) 'Check how different two equation sets are'
     call eqn_compare
  else if (kommando == 'print_eigenvalues') then
     write(*,*) 'Output eigenvalues'
     call print_eigenvalues
  else if (kommando == 'edge') then
     write(*,*) 'Mark the edge of a mask'
     call mask_edge(unit)
  else if (kommando == 'makemask') then
     write(*,*) 'Make mask from center coordinates and radius'
     call make_mask(unit)
  else if (kommando == 'findpatch') then
     write(*,*) 'Find patch with lowest rms'
     call findpatch
  else
     call give_user_info
  end if

  ! And exit
  call mpi_finalize(ierr)
  if (myid == root) then 
     write(*,*)
     write(*,*) '-------------------- QUIET map post-processor completed  ----------'
  end if
  
contains

  !-----------------------------------------------------------------------------------------------
  ! subroutine make_finalmap
  !-----------------------------------------------------------------------------------------------
  
  subroutine make_finalmap(unit)
    implicit none
    
    integer(i4b),              intent(in) :: unit

    character(len=256)                    :: sub1file, sub2file, outprefix, outfile, outtext
    integer(i8b)                          :: npix1, npix2, n, n1, n2, num
    integer(i4b)                          :: n_in, n1_in, n2_in, npix_in, npix2_in, int_in
    integer(i4b)                          :: i, j, k, nmaps, ordering, ordering2
    integer(i4b)                          :: polarisation, polarisation2
    integer(i4b)                          :: nmaps_mask1, nside_mask1, npix_mask1, ordering_mask1
    integer(i4b)                          :: nmaps_mask2, nside_mask2, npix_mask2, ordering_mask2
    integer(i4b)                          :: dummy, cutoff, band
    real(dp)                              :: nullval, maxeigenval, chi_square
    real(dp)                              :: t1, t2, a2t, frequency
    logical(lgt)                          :: anynull, temperature, temperature2, asym1, asym2
    
    real(dp), allocatable, dimension(:,:) :: outmap, matrix1, matrix2 
    real(dp), allocatable, dimension(:,:) :: matrixsum, eigenvectors
    real(dp), allocatable, dimension(:,:) :: invmatrixsum, somematrix, cov1, cov2, rhs1in, rhs2in
    real(dp), allocatable, dimension(:)   :: eigenvals, inveigenvals
    real(dp), allocatable, dimension(:)   :: mapsum, mapdiff, chivec
    real(dp), allocatable, dimension(:)   :: map1, map2, rhs1, rhs2
    real(dp), pointer,  dimension(:)      :: redmap1, redmap2
    real(dp), pointer,  dimension(:,:)    :: redcov1, redcov2
    integer(i4b), pointer,  dimension(:)  :: in2red1, in2red2, healvector1
    integer(i4b), allocatable,dimension(:):: vector1, vector2, nhits
    integer(i4b), allocatable,dimension(:):: mask1, mask2, vector
    
    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'finalmap takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, sub1file)
       call getarg(3, sub2file)
       call getarg(4, outprefix)
    end if

    ! set overall parameters
    nmaps = 3
    band  = 1                                            !OBS Q-band
    if (band ==1) then
       frequency = 43.d0
       if (myid==root)  write(*,*) 'Q-band frequency in GHz is', frequency
    else
       write(*,*) 'Need frequency for Wband'
       stop
    end if

    ! read egn_file1
    ! call read_equ_set_old(unit, sub1file, n1, matrix1, rhs1, temperature, npix1, mask1)
    call read_equ_set(unit, sub1file, n1, ordering, polarisation, matrix1, rhs1in, npix1, mask1)
    write(*,*) n1, '= n for covariance matrix 1'
    write(*,*) int(sqrt(real(npix1/12))), '= nside for map 1'

    ! read eqn_file2
    ! call read_equ_set_old(unit, sub2file, n2, matrix2, rhs2, temperature2, npix2, mask2)
    call read_equ_set(unit, sub2file, n2, ordering2, polarisation2, matrix2, rhs2in, npix2, mask2)
    write(*,*) n2, '= n for covariance matrix 2'
    write(*,*) int(sqrt(real(npix2/12))), '= nside for map 2'

    !For now just using first rhs component (real data)
    allocate(rhs1(n1))
    allocate(rhs2(n2))
    rhs1=rhs1in(:,1)
    rhs2=rhs2in(:,1)
    deallocate(rhs1in,rhs2in)

    ! call check_input
    npix_in = npix1
    npix2_in = npix2
    call check_input(myid, ordering, ordering2, polarisation, polarisation2, npix_in, npix2_in, temperature)
    call check_sym_matrix(matrix1, asym1)
    call check_sym_matrix(matrix2, asym2)
    if (myid==root) then
       if (asym1) write(*,*) 'Covariance matrix 1 is not symmetric.'
       if (asym2) write(*,*) 'Covariance matrix 2 is not symmetric.'
       if (asym1 .or. asym2) write(*,*) 'Covariance matrices are not symmetric. Quiting'
    end if
    if (asym1 .or. asym2) then
       call mpi_finalize(ierr)
       stop
    end if

    ! Solve equation set 1
    allocate(cov1(n1, n1))
    allocate(map1(n1))
    outtext = 'input_map1'
    call solve_eq(myid, unit, outprefix, outtext, matrix1, rhs1, map1, cov1)
    deallocate(matrix1)
    deallocate(rhs1)
    ! Writing input map1 to file
    allocate(outmap(0:npix1-1, nmaps))
    call map2healmap(outmap, map1, mask1, polarisation, npix_in)
    outfile = trim(outprefix) // '_input_map1.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Input map 1 written to file = ', trim(outfile)  
    
    ! Solve equation set 2
    allocate(cov2(n2, n2))
    allocate(map2(n2))
    outtext = 'input_map2'
    call solve_eq(myid, unit, outprefix, outtext, matrix2, rhs2, map2, cov2)
    deallocate(matrix2)
    deallocate(rhs2)
    ! Writing input map1 to file
    call map2healmap(outmap, map2, mask2, polarisation, npix_in)
    outfile = trim(outprefix) // '_input_map2.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Input map 2 written to file = ', trim(outfile)  
    
    ! Removing large and uncommon pixels from map2mask
    n1_in = n1
    n2_in = n2
    npix_in = npix1
    call mask2common(myid, mask1, mask2, in2red1, in2red2, healvector1, n1_in, n2_in, n_in, npix_in)
    n = n_in  
    if (polarisation == 2) n = 2*n
    deallocate(mask1)
    deallocate(mask2)
    !.. and from maps and cov matrices
    call inmap2red(in2red1, map1, cov1, redmap1, redcov1, polarisation)
    deallocate(in2red1)
    deallocate(map1)
    deallocate(cov1)
    call inmap2red(in2red2, map2, cov2, redmap2, redcov2, polarisation)
    deallocate(in2red2)
    deallocate(map2)
    deallocate(cov2)
    
    ! Writing map1 to file
    call map2outmap(outmap, redmap1, healvector1, polarisation)
    outfile = trim(outprefix) // '_red_map1.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Reduced input map1 written to file = ', trim(outfile)
    
    ! Writing map2 to file
    call map2outmap(outmap, redmap2, healvector1, polarisation)
    outfile = trim(outprefix) // '_red_map2.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Reduced input map2 written to file = ', trim(outfile)
    
    ! Writing true diffmap to file
    call map2outmap(outmap, (redmap1 - redmap2)/2, healvector1, polarisation)
    outfile = trim(outprefix) // '_true_diffmap.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* True diffmap written to file = ', trim(outfile)
    
    ! Writing direct sum map to file
    call map2outmap(outmap, (redmap1 + redmap2)/2, healvector1, polarisation)
    outfile = trim(outprefix) // '_direct_summap.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Direct sum map written to file = ', trim(outfile)
    
    ! Writing map2mask to file
    call map2mask2outmap(outmap, healvector1, polarisation)   
    outfile = trim(outprefix) // '_map2mask.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Map2mask written to file = ', trim(outfile)
    
    ! Writing mask to file
    call mask2outmap(outmap, healvector1, polarisation)   
    outfile = trim(outprefix) // '_mask.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Mask written to file = ', trim(outfile)
    
    
    
    ! Need more working space..
    allocate(matrixsum(n, n)) 
    allocate(somematrix(n, n))
    allocate(eigenvectors(n,n)) 
    allocate(eigenvals(n))
    allocate(inveigenvals(n))
    allocate(chivec(n))
    allocate(mapsum(n))
    allocate(mapdiff(n))
    
    ! cov matrix for true diffmap is sum of redcov1 and redcov2
    matrixsum = (redcov1 + redcov2)/4
    
    !call invert_singular_matrix(matrixsum, 0.d0)
    call get_eigen_decomposition(myid, matrixsum, eigenvals, eigenvectors)
    ! check and invert eigenvals
    outtext = 'eigenvals_directsum'
    call invert_eigenvals(myid, eigenvals, unit, outprefix, outtext)
    
    ! write eigenvectors for inverse cov matrix for true diffmap to file
    outfile = trim(outprefix) // '_inv_N_eigen_diffmap.unf'
    open(unit, file=trim(outfile), form='unformatted')
    int_in = n
    write(unit) int_in
    write(unit) eigenvectors
    write(unit) eigenvals
    close(unit)
    write(*,*) '* Eigendecomp for inv cov matrix for true diffmap written to file = ', trim(outfile)
    
    ! invert direct total cov matrix
    somematrix = transpose(eigenvectors)
    do i = 1, n
       somematrix(i,:) = somematrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, matrixsum, n)
    
    ! write inverse cov matrix for true diffmap to file
    outfile = trim(outprefix) // '_inv_N_diffmap.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, matrixsum, .true.)
    write(*,*) '* Inv cov matrix for true diffmap written to file = ', trim(outfile)
    
    ! chi_square check for direct diffmap
    mapdiff = (redmap1 - redmap2)/2
    chi_square = sum(mapdiff * matmul(matrixsum, mapdiff))
    write(*,*) '-----> Chi_squared direct diffmap ==', chi_square
    
    !-------------------------------------------------------------------------------
    ! Making weighted maps
    !-------------------------------------------------------------------------------

    ! Decomposing inverse final covariance matrix into eigenvectors
    call invert_singular_matrix(redcov1, 0.d0)
    call invert_singular_matrix(redcov2, 0.d0)
    matrixsum = redcov1 + redcov2   !inverse finalcov matrix
    call get_eigen_decomposition(myid, matrixsum, inveigenvals, eigenvectors) 
    
    ! Writing eigenvectors for inv final cov matrix to file
    outfile = trim(outprefix) // '_inv_N_eigen.unf'
    open(unit, file=trim(outfile), form='unformatted')
    int_in = n
    write(unit) int_in
    write(unit) eigenvectors
    write(unit) inveigenvals
    close(unit)
    write(*,*) '* Eigendecomposition for inv cov matrix written to file = ', trim(outfile)
    
    ! Writing inverse total covariance matrix to file
    outfile = trim(outprefix) // '_inv_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, matrixsum, .true.)    
    write(*,*) '* Final inv cov matrix written to file = ', trim(outfile)
    
    somematrix = transpose(eigenvectors)
    do i = 1, n
       somematrix(i,:) = somematrix(i,:) * sqrt(inveigenvals(i))
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, matrixsum, n)
    
    ! Writing sqrt of inverse total covariance matrix to file
    outfile = trim(outprefix) // '_sqrt_inv_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, matrixsum, .true.)    
    write(*,*) '* Sqrt of inv cov matrix written to file = ', trim(outfile)
    
    
    ! check and invert eigenvals
    outtext = 'inveigenvals_weightedsum'
    eigenvals = inveigenvals
    call invert_eigenvals(myid, eigenvals, unit, outprefix, outtext)
    
    somematrix = transpose(eigenvectors)
    do i = 1, n
       somematrix(i,:) = somematrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, matrixsum, n)
    
    
    ! Computing diff_map and sum_map
    redmap1 = matmul(redcov1, redmap1)
    redmap2 = matmul(redcov2, redmap2)
    mapsum  = matmul(matrixsum, (redmap1 + redmap2))
    mapdiff = matmul(matrixsum, (redmap1 - redmap2))
    deallocate(redmap1)
    deallocate(redmap2)
    deallocate(redcov1)
    deallocate(redcov2)
    

    ! computing chi_squared
    chivec = matmul(transpose(eigenvectors), mapdiff)
    chivec = chivec * chivec * inveigenvals
    chi_square = sum(chivec)
    write(*,*) '-----> Chi_squared weighted diffmap =', chi_square
    
    ! Output chi_square for each eigenmode
    outfile = trim(outprefix) // '_chisq_weighted_diffmap.txt'
    open(unit, file=trim(outfile))
    do i = 1, n
       write(unit,*) i, chivec(i)
    end do
    close(unit)
    
    ! Writing final map to file
    call map2outmap(outmap, mapsum, healvector1, polarisation)
    outfile = trim(outprefix) // '_map.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Final map written to file = ', trim(outfile)
    
    ! Writing weighted diffmap to file
    call map2outmap(outmap, mapdiff, healvector1, polarisation)
    outfile = trim(outprefix) // '_weighted_diffmap.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Weighted diffmap written to file = ', trim(outfile)
    
    !Writing rms map to outmap
    outmap = -1.6375d30
    if (temperature) then
       do i = 1, n
          outmap(healvector1(i),1) = sqrt(matrixsum(i,i))
       end do
    else
       do i = 1, n/2
          outmap(healvector1(i),2) = sqrt(matrixsum(i,i))
       end do
       do i = 1,n/2
          outmap(healvector1(i),3) = sqrt(matrixsum(i+n/2,i+n/2))
       end do
    end if
    ! Writing rmsmap to file
    outfile = trim(outprefix) // '_rms_map.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Noise rms map written to file = ', trim(outfile)
    
    ! Writing total covariance matrix to file
    outfile = trim(outprefix) // '_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, matrixsum, .false.)
    write(*,*) '* Final cov matrix written to file = ', trim(outfile)
  
    ! clean up
    deallocate(outmap)
    deallocate(eigenvals)
    deallocate(eigenvectors)
    deallocate(somematrix)
    deallocate(matrixsum)
    deallocate(mapsum)
    deallocate(mapdiff)
    deallocate(chivec)
    
  end subroutine make_finalmap

  !-----------------------------------------------------------------------------------------------
  ! subroutine project_out_eigenmodes
  !-----------------------------------------------------------------------------------------------
  
  subroutine project_out_eigenmodes(unit, temperature)
    implicit none

    integer(i4b),              intent(in) :: unit
    logical(lgt),              intent(in) :: temperature
    

    character(len=256)                    :: infile, outfile, map2maskfile, matrixfile
    character(len=256)                    :: cutoffname, outprefix, arg
    integer(i8b)                          :: n
    integer(i4b)                          :: i, int_in, nmaps, nside, npix, ordering, cutoff
    integer(i4b)                          :: nmaps_mask, nside_mask, npix_mask, ordering_mask
    integer(i4b)                          :: polarisation
    real(dp)                              :: nullval, dummy, maxeigenval, t1, t2
    logical                               :: anynull, mismatch, maponly
    
    real(dp), allocatable, dimension(:,:) :: inmap, outmap, mask, eigenvectors, somematrix, newcov
    real(dp), allocatable, dimension(:)   :: vector, eigenvals
    
    ! Get parameters
    if (iargc() < 6) then
       write(*,*) 'cutoff takes 5 parameters'
       call give_user_info
    else 
       call getarg(2, map2maskfile)
       call getarg(3, matrixfile)
       call getarg(4, infile)
       call getarg(5, outprefix)
       call getarg(6, cutoffname)
    end if
    maponly = .false.
    if(iargc() > 6) then
       call getarg(7, arg)
       if(arg == "maponly") maponly = .true.
    end if

    read(cutoffname,*) cutoff

    if (temperature) then
       write(*,*) 'Cutting ',trim(cutoffname),' eigenmodes from temperature data'
    else
       write(*,*) 'Cutting ',trim(cutoffname),' eigenmodes from polarisation data'
    end if

    ! read matrix file = inv N eigenvectors
    open(unit, file=trim(matrixfile), form='unformatted')
    read(unit) int_in
    n = int_in
    write(*,*) n, '= n for matrix'
    allocate(eigenvectors(n,n)) 
    read(unit) eigenvectors
    allocate(eigenvals(n))
    read(unit) eigenvals
    close(unit)
    
    !read map
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    allocate(inmap(0:npix-1,nmaps))
    call read_bintab (infile, inmap, npix, nmaps, nullval, anynull)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 
    
    !Read map2mask
    npix_mask = getsize_fits(map2maskfile, nmaps=nmaps_mask, ordering=ordering_mask, nside=nside_mask)
    
    ! Checking that map2mask and input map is of same size
    mismatch = .false.
    if (nside    /= nside_mask) mismatch = .true.
    if (nmaps    /= nmaps_mask) mismatch = .true.
    if (ordering /= ordering_mask) mismatch = .true.
    if (mismatch) then
       if (myid==root) then 
          write(*,*) nside_mask, '= nside,', nmaps_mask, '= nmaps,',ordering_mask,'= ordering' 
          write(*,*) 'Input map and map2mask is not of sime size. Quiting'
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
    allocate(mask(0:npix-1,nmaps))
    call read_bintab (map2maskfile, mask, npix, nmaps, nullval, anynull)
    
    ! Building n-dim vector from map2mask and inmap
    allocate(vector(n))
    vector = -1.6375d30
    if (temperature) then
       do i = 0, npix-1
          dummy = mask(i,1)
          if (dummy /= -1.d0) vector(int(dummy,i4b)) = inmap(i,1)  ! for temperature 
       end do
    else
       do i = 0, npix-1
          dummy = mask(i,2)
          if (dummy /= -1.d0) then
             vector(nint(dummy,i4b))     = inmap(i,2) 
             vector(nint(dummy,i4b)+n/2) = inmap(i,3) 
          end if
       end do
    end if
    
    !Multiplying above vector by matrix of eigenvectors to get new vector
    vector = matmul(transpose(eigenvectors),vector)
    maxeigenval=maxval(abs(eigenvals))
    
    !do i = 1, n
    !   if (abs(eigenvals(i))<1d-12 * maxeigenval) vector(i) = 0
    !end do
    do i = 1, n
       if (eigenvals(i)<1d-12) then 
          write(*,*) i, eigenvals(i), 'removed'
          vector(i) = 0.d0
          eigenvals(i) = 0.d0
       end if
    end do
    
    write(*,*) eigenvals(1), eigenvals(cutoff), eigenvals(n), 'eigenvals'
    if (eigenvals(1) < eigenvals(n)) then
       write(*,*) 'cutting beginning'
       do i = 1, cutoff
          vector(i) = 0.d0
          eigenvals(i) = 0.d0
       end do
    else
       write(*,*) 'cutting end'
       do i = n-cutoff+1, n
          vector(i) = 0.d0
          eigenvals(i) = 0.d0
       end do
    end if
    
    vector = matmul(eigenvectors,vector)
    
    !Writing new vector to outmap
    allocate(outmap(0:npix-1,nmaps))
    outmap = -1.6375d30
    if (temperature) then
       do i = 0, npix-1
          dummy = mask(i,1)
          if (dummy /= -1.d0) outmap(i,1) = vector(nint(dummy,i4b))
       end do
    else
       do i = 0, npix-1
          dummy = mask(i,2)
          if (dummy /= -1.d0) then
             outmap(i,2) = vector(nint(dummy,i4b))
             outmap(i,3) = vector(nint(dummy,i4b)+n/2)
          end if
       end do
    end if
    
    ! Writing outmap to file
    outfile = trim(outprefix) // '_cut' //trim(cutoffname)// '_map.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) 'Map written to file = ', trim(outfile)
    if(maponly) goto 1
    
    allocate(somematrix(n,n)) 
    allocate(newcov(n,n)) 
    
    ! Removing cut eigenmodes from inv cov matrix
    call cpu_time(t1)
    somematrix = transpose(eigenvectors)
    do i = 1, n
       somematrix(i,:) = somematrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, newcov, n)
    call cpu_time(t2)
    write(*,*) 'Time spent on matrix multiplication:', t2-t1 
    
    ! Writing newinverse total covariance matrix to file
    outfile = trim(outprefix) // '_cut' //trim(cutoffname)// '_inv_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, newcov, .true.)
    write(*,*) 'Inv cov matrix written to file = ', trim(outfile)
    
    ! Corresponding sqrt of  inv cov matrix
    call cpu_time(t1)
    somematrix = transpose(eigenvectors)
    do i = 1, n
       somematrix(i,:) = somematrix(i,:) * sqrt(eigenvals(i))
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, newcov, n)
    call cpu_time(t2)
    write(*,*) 'Time spent on matrix multiplication:', t2-t1 
    
    ! Writing sqrt of inv total covariance matrix to file
    outfile = trim(outprefix) // '_cut' //trim(cutoffname)// '_sqrt_inv_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, newcov, .true.)
    write(*,*) 'Sqrt of inv cov matrix written to file = ', trim(outfile)
    
1   if(allocated(mask))         deallocate(mask)
    if(allocated(somematrix))   deallocate(somematrix)
    if(allocated(newcov))       deallocate(newcov)
    if(allocated(inmap))        deallocate(inmap)
    if(allocated(outmap))       deallocate(outmap)
    if(allocated(vector))       deallocate(vector)
    if(allocated(eigenvals))    deallocate(eigenvals)
    if(allocated(eigenvectors)) deallocate(eigenvectors)
    
  end subroutine project_out_eigenmodes
  
  !-----------------------------------------------------------------------------------------------
  ! subroutine project_out_lowmodes
  !-----------------------------------------------------------------------------------------------
  
  subroutine project_out_lowmodes(unit)
    implicit none
    
    integer(i4b),              intent(in) :: unit
    
    character(len=256)                    :: infile, outprefix, map2maskfile, matrixfile
    character(len=256)                    :: outfile, cutoffname, nocovname
    integer(i8b)                          :: n, nmodes
    integer(i4b)                          :: i, j, k, l, m, p, int_in
    integer(i4b)                          :: nmaps, nside, npix, ordering, lmax, polarisation
    integer(i4b)                          :: nmaps_mask, nside_mask, npix_mask, ordering_mask
    real(dp)                              :: nullval, dummy, maxeigenval, chi_square
    logical(lgt)                          :: anynull, mismatch, inv, temperature, nocov
    
    real(dp), allocatable, dimension(:,:) :: inmap, outmap, map2mask, eigenvectors, fmap, dmap
    real(dp), allocatable, dimension(:,:) :: invcov, somematrix, newmatrix, finf, invfinf, invnf
    real(dp), allocatable, dimension(:)   :: vector, eigenvals, nmap
    
    nocov = .false.
    
    ! Get parameters
    if (iargc() /= 6 .and. iargc() /= 7) then
       write(*,*) 'lcut takes 5(6) parameters'
       call give_user_info
    else 
       call getarg(2, map2maskfile)
       call getarg(3, matrixfile)
       call getarg(4, infile)
       call getarg(5, outprefix)
       call getarg(6, cutoffname)
       if (iargc() == 7) call getarg(7, nocovname)
    end if
    
    if (trim(nocovname)== 'nocov') nocov=.true.
    read(cutoffname,*) lmax
    outprefix = trim(outprefix) // '_lcut' //trim(cutoffname)    
    
    if (nocov) then
       if (myid==0)  write(*,*) 'Running with covariance matrix put to unity.'
       if (myid==0)  write(*,*) 'Just for making pretty pictures, not for further analysis'
    end if
    
    ! read matrix file = inv cov matrix
    call read_covmatrix(unit, matrixfile, ordering, polarisation, invcov, inv, n)
    
    if (polarisation == 1) then
       temperature = .true.
    else if (polarisation == 2) then
       temperature = .false.
    else
       if (myid==root) then 
          write(*,*) "Neither temperature nor polarisation data. quiting"
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
    if (temperature) then
       nmodes = (lmax+1)**2
       write(*,fmt='(a,i5,a,a,a)') 'Projecting out up to',nmodes,' modes with l <= ',trim(cutoffname),' for temperature' 
    else
       nmodes = 2*(lmax+1)**2
       write(*,fmt='(a,i5,a,a,a)') 'Projecting out up to',nmodes,' modes with l <= ',trim(cutoffname),' for polarisation' 
    end if
    
    if (nocov) then
       do i = 1,n
          do j = 1,n
             if (i == j) then 
                invcov(i,j)=1.d0
             else
                invcov(i,j)=0.d0
             end if
          end do
       end do
    else
       ! Checking that inverse cov matrix is symmetric
       do i = 1, n
          do j = i, n
             if ( abs((invcov(i,j)-invcov(j,i))/(invcov(i,j)+invcov(j,i))) > 1d-5 .and. &
                  & abs(invcov(i,j))/sqrt(invcov(i,i)*invcov(j,j)) > 1.d-8) then
                write(*,*) i,j, real(invcov(i,j),sp), real(invcov(j,i),sp), &
                     & real(abs(invcov(i,j))/sqrt(invcov(i,i)*invcov(j,j)),sp)
                if (myid==root) then 
                   write(*,*) 'Covariance matrix is not symmetric. Quiting'
                end if
                call mpi_finalize(ierr)
                stop
             end if
          end do
       end do
    end if
    
    !read map
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    allocate(inmap(0:npix-1,nmaps))
    call read_bintab (infile, inmap, npix, nmaps, nullval, anynull)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 
    
    ! check size of map2mask
    npix_mask = getsize_fits(map2maskfile, nmaps=nmaps_mask, ordering=ordering_mask, nside=nside_mask)
    ! Checking that map2mask and input map is of same size
    mismatch = .false.
    if (nside    /= nside_mask) mismatch = .true.
    if (nmaps    /= nmaps_mask) mismatch = .true.
    if (ordering /= ordering_mask) mismatch = .true.
    if (mismatch) then
       if (myid==root) then 
          write(*,*) nside_mask, '= nside,', nmaps_mask, '= nmaps,',ordering_mask,'= ordering' 
          write(*,*) 'Input map and map2mask is not of same size. Quiting'
       end if
       call mpi_finalize(ierr)
       stop
    end if
    !Read map2mask
    allocate(map2mask(0:npix-1,nmaps))
    call read_bintab (map2maskfile, map2mask, npix, nmaps, nullval, anynull)
    

! nmap = n version of inmap
!allocate(dmap(n, 1))
!call healmap2map(dmap(:,1), inmap, map2mask, polarisation)
!dmap(1:n/2,1) = dmap(1:n/2,1) + 2000.d0

!allocate(fmap(n,1))
!fmap = 1.d0
!dummy = sum(fmap*matmul(invcov,dmap))
!chi_square = sum(fmap*matmul(invcov,fmap))
!dummy = dummy/chi_square
!write(*,*) dummy, sum(dmap)/n
!stop

    ! putting the unwanted modes into fmap 
    allocate(fmap(n, nmodes))
    call get_lowmodes(fmap, lmax, n, nmodes, nmaps, nside, ordering, temperature, map2mask)

    ! invnf = invcov*fmap
    allocate(invnf(n, nmodes))
    call DGEMM('n', 'n', n, nmodes, n, 1.d0, invcov, n, fmap, n, 0.d0, invnf, n)
    
    ! finf = fmap^T*invcov*fmap
    allocate(finf(nmodes, nmodes))
    call DGEMM('t', 'n', nmodes, nmodes, n, 1.d0, fmap, n, invnf, n, 0.d0, finf, nmodes)
    
    ! inverting finf
    allocate(eigenvectors(nmodes,nmodes))
    allocate(eigenvals(nmodes))
    call get_eigen_decomposition(myid, finf, eigenvals, eigenvectors)
    call invert_eigenvals(myid, eigenvals, unit, outprefix, kept=k, acc=1d-4)
    write(*,*) nmodes-k, 'of', nmodes, 'eigenvalues removed.'
    write(*,*) 'Removing a total of', k, 'modes out of',n
    ! finally inverting finf
    allocate(somematrix(nmodes, nmodes))
    allocate(invfinf(nmodes, nmodes))
    somematrix = transpose(eigenvectors)
    do i = 1, nmodes
       somematrix(i,:) = somematrix(i,:) * eigenvals(i)
    end do
    call DGEMM('N', 'N', nmodes, nmodes, nmodes, 1.d0, eigenvectors, nmodes, somematrix, nmodes, 0.d0, invfinf, nmodes)
    deallocate(somematrix)
    deallocate(eigenvectors)
    deallocate(eigenvals)
    
    ! nmap = n version of inmap
    allocate(nmap(n))
    call healmap2map(nmap, inmap, map2mask, polarisation)
        
    ! chi_square for inmap
    chi_square = sum(nmap * matmul(invcov, nmap))
    call write_chi_log(unit, chi_square, int(n,i4b), outprefix, helptext='before')
    
    ! inmap = inmap - f* (f^T*N^-1*f)^-1 * (f^T*N^-1*inmap)
    allocate(newmatrix(nmodes, n))
    call DGEMM('N', 'T', nmodes, n, nmodes, 1.d0, invfinf, nmodes, invnf, n, 0.d0, newmatrix, nmodes)
    allocate(somematrix(nmodes, 1))
    call DGEMM('N', 'N', nmodes, 1, n, 1.d0, newmatrix, nmodes, nmap, n, 0.d0, somematrix, nmodes)
    deallocate(newmatrix)
    allocate(vector(n))
    call DGEMM('N', 'N', n, 1, nmodes, 1.d0, fmap, n, somematrix, nmodes, 0.d0,vector, n)
    nmap = nmap  - vector
    deallocate(somematrix)
    
    !Writing modified map to outmap
    allocate(outmap(0:npix-1,nmaps))
    outmap = -1.6375d30
    if (temperature) then
       do i = 0, npix-1
          dummy = map2mask(i,1)
          if (dummy /= -1.d0) outmap(i,1) = nmap(nint(dummy,i4b))
       end do
    else
       do i = 0, npix-1
          dummy = map2mask(i,2)
          if (dummy /= -1.d0) then
             outmap(i,2) = nmap(nint(dummy,i4b))
             outmap(i,3) = nmap(nint(dummy,i4b)+n/2)
          end if
       end do
    end if
    
    ! Writing outmap to file
    outfile = trim(outprefix) // '_map.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Map written to file = ', trim(outfile)
    
    !Writing removed modes to outmap
    outmap = -1.6375d30
    if (temperature) then
       do i = 0, npix-1
          dummy = map2mask(i,1)
          if (dummy /= -1.d0) outmap(i,1) = vector(nint(dummy,i4b))
       end do
    else
       do i = 0, npix-1
          dummy = map2mask(i,2)
          if (dummy /= -1.d0) then
             outmap(i,2) = vector(nint(dummy,i4b))
             outmap(i,3) = vector(nint(dummy,i4b)+n/2)
          end if
       end do
    end if
    
    ! Writing outmap to file
    outfile = trim(outprefix) // '_removed_modes_map.fits'
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Removed part of map written to file = ', trim(outfile)
    
    deallocate(vector)
    deallocate(inmap)
    deallocate(fmap)
    deallocate(outmap)
    
    if (nocov) return

write(*,*) 'mangler log'
    
    ! invcov = incov - N^-1 * f * (f^T*N^-1*f)^-1 * f^T * N^-1
    allocate(somematrix(n, nmodes))
    call DGEMM('n', 'n', n, nmodes, nmodes, 1.d0, invnf, n, invfinf, nmodes, 0.d0, somematrix, n)
    allocate(newmatrix(n,n))
    call DGEMM('n', 't', n, n, nmodes, 1.d0, somematrix, n, invnf, n, 0.d0, newmatrix, n)
    invcov = invcov - newmatrix
    deallocate(newmatrix)  
    deallocate(somematrix)
    deallocate(invnf)
    deallocate(invfinf)
    
    ! chi_square for outmap
    chi_square = sum(nmap * matmul(invcov, nmap))
    call write_chi_log(unit, chi_square, int(n-k,i4b), outprefix, helptext='after')

    if (trim(nocovname)== 'onlymap') return
    
    ! Writing new inverse covariance matrix to file
    outfile = trim(outprefix) // '_inv_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, invcov, .true.)
    write(*,*) '* Inv cov matrix written to file = ', trim(outfile)
    
    ! Decomposing new inv cov matrix into eigenvectors
    allocate(eigenvectors(n, n))
    allocate(eigenvals(n))
    call get_eigen_decomposition(myid, invcov, eigenvals, eigenvectors)
    allocate(somematrix(n, n))
    ! Find sqrt of inverse cov
    somematrix = transpose(eigenvectors)
    maxeigenval = maxval(abs(eigenvals))
    do i = 1, n
       if (eigenvals(i) > 1d-12 * maxeigenval) then
          somematrix(i,:) = somematrix(i,:) * sqrt(eigenvals(i))         
       else
          somematrix(i,:) = 0.d0
       end if
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, invcov, n)
    
    ! Writing sqrt of new inverse covariance matrix to file
    outfile = trim(outprefix) // '_sqrt_inv_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, invcov, .true.)
    write(*,*) '* Sqrt of inv cov matrix written to file = ', trim(outfile)

    ! Find non-inverse cov
    somematrix = transpose(eigenvectors)
    do i = 1, n
       if (eigenvals(i) > 1d-12 * maxeigenval) then
          somematrix(i,:)  = somematrix(i,:) * 1/(eigenvals(i))
       else
          somematrix(i,:) = 0.d0
       end if
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, invcov, n)

    ! Writing non-inverse covariance matrix to file
    outfile = trim(outprefix) // '_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, invcov, .false.)
    write(*,*) '* Cov matrix written to file = ', trim(outfile)

    ! Find sqrt of non-inverse cov
    somematrix = transpose(eigenvectors)
    do i = 1, n
       if (eigenvals(i) > 1d-12 * maxeigenval) then
          somematrix(i,:)  = somematrix(i,:) * 1/sqrt(eigenvals(i))
       else
          somematrix(i,:) = 0.d0
       end if
    end do
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, somematrix, n, 0.d0, invcov, n)

    ! Writing sqrt of non-inverse covariance matrix to file
    outfile = trim(outprefix) // '_sqrt_N.unf'
    call write_covmatrix(unit, outfile, ordering, polarisation, invcov, .false.)
    write(*,*) '* Sqrt of cov matrix written to file = ', trim(outfile)
    
    deallocate(somematrix)
    deallocate(eigenvectors)
    deallocate(eigenvals)
    deallocate(invcov)
    deallocate(map2mask)
    
  end subroutine project_out_lowmodes

  !-----------------------------------------------------------------------------------------------
  ! subroutine coupling_kernel
  !-----------------------------------------------------------------------------------------------

  subroutine coupling_kernel(unit)
    implicit none
    
    integer(i4b),              intent(in) :: unit
    
    character(len=256)                    :: infile, outprefix, outfile, cutoffname
    character(len=3)                      :: filnum
    integer(i4b)                          :: i, j, k, l, m, p, int_in, nestpix, ringpix, pix
    integer(i4b)                          :: nmaps, nside, npix, ordering, lmax, npol
    integer(i4b)                          :: ndim, ierr, lmax3    
    real(dp)                              :: nullval, tall, norm, l3max, l3min
    logical(lgt)                          :: anynull, mismatch, inv, temperature, nocov
  
    real(dp),     allocatable, dimension(:,:)   :: mask, kernel, dw8
    real(dp),     allocatable, dimension(:)     :: power, cgc
    complex(dpc), allocatable, dimension(:,:,:) :: alm
    real(dp),                  dimension(2)     :: z

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'kernel takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, infile)
       call getarg(3, outprefix)
       call getarg(4, cutoffname)
    end if
    
    read(cutoffname,*) lmax
    write(*,fmt='(a,a,a)') 'Calculating coupling kernel matrix for lmax = ',trim(cutoffname),' for temperature' 
    npol = 1
    
    !read mask
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab (infile, mask, npix, nmaps, nullval, anynull)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 
    if (ordering ==2) then                       ! map2alm requieres ringed
       call convert_nest2ring(nside ,mask)
       ordering = 1
    end if

    ! finding alm's for given skymap (mask)
    lmax3 = 2*lmax
    allocate(dw8(1:2*nside, 1:npol))
    allocate(alm(1:npol, 0:lmax3, 0:lmax3))
    dw8 = 1.d0
    z = 0.d0
    call map2alm(nside, lmax3, lmax3, mask(:, nmaps), alm, z, dw8)
    deallocate(dw8)
    deallocate(mask)

    ! and power spectrum
    allocate(power(0:lmax3))
    power = 0.d0
    do i = 0, lmax3
       do m = 0, i
          if (m == 0) then
             power(i) = power(i) + abs(alm(npol, i, m))**2             
          else
             power(i) = power(i) + 2*abs(alm(npol, i, m))**2
          end if
!write(*,*) i,m,alm(npol,i,m)       
       end do
       power(i) = power(i)/(2.d0*i+1.d0)
    end do
    deallocate(alm)

    ! Output powerspectrum
    outfile = trim(outprefix) // '_powerspec.txt'
    open(unit, file=trim(outfile))
    do i = 0, lmax3
       write(unit,*) i, power(i)
    end do
    close(unit)

    ! kernel matrix
    allocate(kernel(0:lmax, 0:lmax))
    kernel = 0.d0
    do j = 0, lmax
       write(*,*) j 
       do i = 0, lmax
          ndim = i+j - abs(i-j)+1
          allocate(cgc(ndim))
          call DRC3JJ (real(i,dp), real(j,dp), 0.d0, 0.d0,l3min, l3max, cgc,  ndim, ierr)

!write(*,*) i,j
!do m = 1, ndim
!write(*,*) m, cgc(m)
!end do
!write(*,*)
          do k = abs(i-j), abs(i-j) + ndim -1
             kernel(i,j) = kernel(i,j) + (2.d0*k+1.d0)*power(k)*cgc(k-abs(i-j)+1)**2
!write(*,*) i,j,k,kernel(i,j)
!write(*,*) k-abs(i-j)+1,cgc(k-abs(i-j)+1), cgc(k-abs(i-j)+1)**2
!write(*,*)
          end do
          deallocate(cgc)
       end do
       kernel(:,j)=kernel(:,j)* (2.d0*j+1.d0)
    end do
    kernel = kernel/(4*pi)

    deallocate(power)
 
   ! Write to file
    do j = 0, lmax
       call int2string(j, filnum)
       outfile=trim(outprefix)//'_kernel_'//filnum//'.txt'
       open(unit, file=trim(outfile))
       do i = 0, lmax
          write(unit,*) i, kernel(i,j)
       end do
       close(unit)
       write(*,*) 'written to file = ', trim(outfile)
    end do

    tall=maxval(kernel(:,25))
write(*,*), tall, tall/2.7d0
    outfile=trim(outprefix)//'_e.txt'
    open(unit, file=trim(outfile))
    do i = 0, lmax
       write(unit,*) i, tall/2.7d0
    end do
    close(unit)
    write(*,*) 'written to file = ', trim(outfile)

    outfile=trim(outprefix)//'_sum_halflmax.txt'
    open(unit, file=trim(outfile))
    do i = 0, lmax
       tall = 0.d0
       do j = 0, lmax/2
          tall = tall + kernel(i,j)
       end do
    write(unit,*) i, tall/2.7d0
    end do
    close(unit)
    write(*,*) 'written to file = ', trim(outfile)

    outfile=trim(outprefix)//'_prod_5.txt'
    open(unit, file=trim(outfile))
    do i = 0, lmax
       tall = 1.d0
       do j = 0, 5
          norm=maxval(kernel(:,j))
          tall = tall * (1.d0 - kernel(i,j)/norm)
       end do
    write(unit,*) i, 1.d0 - tall
    end do
    close(unit)
    write(*,*) 'written to file = ', trim(outfile)

    outfile=trim(outprefix)//'_prod_10.txt'
    open(unit, file=trim(outfile))
    do i = 0, lmax
       tall = 1.d0
       do j = 0, 10
          norm=maxval(kernel(:,j))
          tall = tall * (1.d0 - kernel(i,j)/norm)
       end do
    write(unit,*) i, 1.d0 - tall
    end do
    close(unit)
    write(*,*) 'written to file = ', trim(outfile)

    outfile=trim(outprefix)//'_prod_25.txt'
    open(unit, file=trim(outfile))
    do i = 0, lmax
       tall = 1.d0
       do j = 0, 25
          norm=maxval(kernel(:,j))
          tall = tall * (1.d0 - kernel(i,j)/norm)
       end do
    write(unit,*) i, 1.d0 - tall
    end do
    close(unit)
    write(*,*) 'written to file = ', trim(outfile)



    deallocate(kernel)

  end subroutine coupling_kernel

!-----------------------------------------------------------------------------------------------
! subroutine project_out_nullspace
!-----------------------------------------------------------------------------------------------

subroutine project_out_nullspace(unit)
  implicit none

  integer(i4b),              intent(in) :: unit

  character(len=256)                    :: infile, outfile, map2maskfile, matrixfile
  integer(i8b)                          :: n
  integer(i4b)                          :: i, int_in, nmaps, nside, npix, ordering
  integer(i4b)                          :: nmaps_mask, nside_mask, npix_mask, ordering_mask
  real(dp)                              :: nullval, dummy, maxeigenval
  logical                               :: anynull, mismatch

  real(dp), allocatable, dimension(:,:) :: inmap, outmap, matrix, mask, eigenvectors
  real(dp), allocatable, dimension(:)   :: vector, eigenvals
    
  ! Get parameters
  if (iargc() /= 5) then
     write(*,*) 'proj0 takes 4 parameters'
     call give_user_info
  else 
     call getarg(2, map2maskfile)
     call getarg(3, matrixfile)
     call getarg(4, infile)
     call getarg(5, outfile)
  end if

  ! read matrix file
  open(unit, file=trim(matrixfile), form='unformatted')
  read(unit) int_in
  n = int_in
  write(*,*) n, '= n for matrix'
  allocate(matrix(n,n))
  read(unit) matrix
  close(unit)

  !read map
  npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
  allocate(inmap(0:npix-1,nmaps))
  call read_bintab (infile, inmap, npix, nmaps, nullval, anynull)
  if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 

  !Read map2mask
  npix_mask = getsize_fits(map2maskfile, nmaps=nmaps_mask, ordering=ordering_mask, nside=nside_mask)
  
  ! Checking that map2mask and input map is of same size
  mismatch = .false.
  if (nside    /= nside_mask) mismatch = .true.
  if (nmaps    /= nmaps_mask) mismatch = .true.
  if (ordering /= ordering_mask) mismatch = .true.
  if (mismatch) then
     if (myid==root) then 
        write(*,*) nside_mask, '= nside,', nmaps_mask, '= nmaps,',ordering_mask,'= ordering' 
        write(*,*) 'Input map and map2mask is not of sime size. Quiting'
     end if
     call mpi_finalize(ierr)
     stop
  end if

  allocate(mask(0:npix-1,nmaps))
  call read_bintab (map2maskfile, mask, npix, nmaps, nullval, anynull)

  ! Building n-dim vector from map2mask and inmap
  allocate(vector(n))
  vector = -1.6375d30
  do i = 0, npix-1
     dummy = mask(i,1)
     if (dummy /= -1.d0) vector(int(dummy,i4b)) = inmap(i,1)  ! for temperature 
  end do
  
  !Decomposing matrix into eigenvectors
  allocate(eigenvals(n))
  allocate(eigenvectors(n,n))
  call get_eigen_decomposition(myid, matrix, eigenvals, eigenvectors)

  !Multiplying vector by matrix of eigenvectors to get new vector
  vector = matmul(transpose(eigenvectors),vector)
  maxeigenval=maxval(abs(eigenvals))
  do i = 1, n
     if (abs(eigenvals(i))<1d-12 * maxeigenval) vector(i) = 0
  end do
  vector = matmul(eigenvectors,vector)

  !Writing new vector to outmap
  allocate(outmap(0:npix-1,nmaps))
  outmap = -1.6375d30
  do i = 0, npix-1
     dummy = mask(i,1)
     if (dummy /= -1.d0) outmap(i,1) = vector(nint(dummy,i4b))
  end do

  ! Writing outmap to file
  call write_map(outmap, ordering, outfile)

  deallocate(matrix)
  deallocate(mask)
  deallocate(inmap)
  deallocate(outmap)
  deallocate(vector)
  deallocate(eigenvals)
  deallocate(eigenvectors)

end subroutine project_out_nullspace

!-----------------------------------------------------------------------------------------------
! subroutine do_matmul
!-----------------------------------------------------------------------------------------------

subroutine do_matmul(unit)
  implicit none

  integer(i4b),              intent(in) :: unit

  character(len=256)                    :: infile, outfile, map2maskfile, matrixfile
  integer(i8b)                          :: n
  integer(i4b)                          :: i, int_in, nmaps, nside, npix, ordering
  integer(i4b)                          :: nmaps_mask, nside_mask, npix_mask, ordering_mask
  real(dp)                              :: nullval, dummy
  logical                               :: anynull, mismatch

  real(dp), allocatable, dimension(:,:) :: inmap, outmap, matrix, mask
  real(dp), allocatable, dimension(:)   :: vector
    
  ! Get parameters
  if (iargc() /= 5) then
     write(*,*) 'matmul takes 4 parameters'
     call give_user_info
  else 
     call getarg(2, map2maskfile)
     call getarg(3, matrixfile)
     call getarg(4, infile)
     call getarg(5, outfile)
  end if

!  call generate_testdata

  ! read matrix file
  open(unit, file=trim(matrixfile), form='unformatted')
  read(unit) int_in
  n = int_in
  write(*,*) n, '= n for matrix'
  allocate(matrix(n,n))
  read(unit) matrix
  close(unit)

  !read map
  npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
  allocate(inmap(0:npix-1,nmaps))
  call read_bintab (infile, inmap, npix, nmaps, nullval, anynull)
  if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 

  !Read map2mask
  npix_mask = getsize_fits(map2maskfile, nmaps=nmaps_mask, ordering=ordering_mask, nside=nside_mask)
  
  ! Checking that map2mask and input map ia of same size
  mismatch = .false.
  if (nside    /= nside_mask) mismatch = .true.
  if (nmaps    /= nmaps_mask) mismatch = .true.
  if (ordering /= ordering_mask) mismatch = .true.
  if (mismatch) then
     if (myid==root) then 
        write(*,*) nside_mask, '= nside,', nmaps_mask, '= nmaps,',ordering_mask,'= ordering' 
        write(*,*) 'Input map and map2mask is not of sime size. Quiting'
     end if
     call mpi_finalize(ierr)
     stop
  end if

  allocate(mask(0:npix-1,nmaps))
  call read_bintab (map2maskfile, mask, npix, nmaps, nullval, anynull)

  ! Building n-dim vector from map2mask and inmap
  allocate(vector(n))
  vector = -1.6375d30
  do i = 0, npix-1
     dummy = mask(i,1)
     if (dummy /= -1.d0) vector(int(dummy,i4b)) = inmap(i,1)  ! for temperature 
  end do

  !Multiplying vector by matrix to get new vector
  vector = matmul(matrix,vector)


  !Writing new vector to outmap
  allocate(outmap(0:npix-1,nmaps))
  outmap = -1.6375d30
  do i = 0, npix-1
     dummy = mask(i,1)
     if (dummy /= -1.d0) outmap(i,1) = vector(nint(dummy,i4b))
  end do

  ! Writing outmap to file
  call write_map(outmap, ordering, outfile)

  deallocate(matrix)
  deallocate(mask)
  deallocate(inmap)
  deallocate(outmap)
  deallocate(vector)

end subroutine do_matmul

!-----------------------------------------------------------------------------------------------
! subroutine solve_equation
!-----------------------------------------------------------------------------------------------

subroutine solve_equation(unit)
  implicit none

  integer(i4b),              intent(in) :: unit

  character(len=256)                    :: sub1file, outprefix, outfile
  integer(i8b)                          :: npix1, n1, mem
  integer(i4b)                          :: i, j, k, nmaps, ordering, int_in, polarisation
  real(dp)                              :: maxeigenval, adiag, mdiag
  logical(lgt)                          :: mismatch, temperature

  real(dp), allocatable, dimension(:,:) :: outmap, matrix1, eigenvectors, rhs1in
  real(dp), allocatable, dimension(:)   :: eigenvals, inveigenvals
  real(dp), allocatable, dimension(:)   :: map1, rhs1
  integer(i4b), allocatable,dimension(:):: mask1

  ! Get parameters
  if (iargc() /= 3) then
     write(*,*) 'solve takes 2 parameters'
     call give_user_info
  else 
     call getarg(2, sub1file)
     call getarg(3, outprefix)
  end if

  mem = 0

  write(*,*) "i8b = ", i8b

  ! read egn_file1
  call read_equ_set(unit, sub1file, n1, ordering, polarisation, matrix1, rhs1in, npix1, mask1)
  write(*,*) n1, '= n for covariance matrix 1'
  write(*,*) int(sqrt(real(npix1/12))), '= nside for map 1'

  !For now just using first rhs component (real data)
  allocate(rhs1(n1))
  rhs1=rhs1in(:,1)
  deallocate(rhs1in)

  if (polarisation == 1) then
     temperature = .true.
  else if (polarisation == 2) then
     temperature = .false.
  else
     write(*,*) 'Unknown polarisation flag = ', polarisation
     stop
  end if

!  open(unit, file=trim(sub1file), form='unformatted')
!  read(unit) int_in
!  n1 = int_in
!  write(*,*) n1, '= n for covariance matrix 1'
!  mem = mem + n1*n1*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
!  allocate(matrix1(n1,n1))                ! inverse covariance matrix
!  read(unit) matrix1
!  allocate(rhs1(n1))                      !rhs
!  mem = mem+n1*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
!  read(unit) rhs1
!  read(unit) temperature
!  read(unit) int_in
!  npix1 = int_in
!  write(*,*) int(sqrt(real(npix1/12))), '= nside for map 1'
!  allocate(mask1(0:npix1-1))              ! map2mask
!  mem = mem+npix1*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
!  read(unit) mask1
!  close(unit)

  if (temperature) then
     write(*,*) 'Processing temperature data'
  else
     write(*,*) 'Processing polarisation data'
  end if

  ! Checking that input matrices are symmetric
  write(*,*) "Asym: ", asymmetry(matrix1)
  mismatch = .false.
  outfile = trim(outprefix) // '_asym.txt'
  open(unit, file=trim(outfile))
  do i = 1, n1
     do j = i, n1
        if ( abs((matrix1(i,j)-matrix1(j,i))/(matrix1(i,j)+matrix1(j,i))) > 1d-5 ) then 
           write(unit,fmt="(2i6,4e18.10,i6)") i,j, matrix1(i,j), matrix1(j,i), adiag, mdiag, 1
           mismatch = .true.
        end if
     end do
  end do
  close(unit)
  if (mismatch) then
     if (myid==root) then 
        write(*,*) 'Covariance matrix is not symmetric. Quiting'
     end if
     call mpi_finalize(ierr)
     stop
  end if

  ! Setting overall parameters
  nmaps    = 3
  ordering = 2              !obs 2=nest, what tod2map gives out

  !solving to find input maps

  ! decompose inverse covariance matrix1 into eigenvectors
  allocate(eigenvals(n1))
  mem = mem + n1*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
  allocate(inveigenvals(n1))
  mem = mem + n1*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
  allocate(eigenvectors(n1,n1))
  mem = mem + n1**2*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
  mem = mem + 2*n1**2*8; write(*,fmt="(a,i12,a)") "mem ", mem, " B"
  call get_eigen_decomposition(myid, matrix1, inveigenvals, eigenvectors)
  deallocate(matrix1)

  ! Output inveigenvals input matrix1
  outfile = trim(outprefix) // '_inveigenvals.dat'
  open(unit, file=trim(outfile))
  do i = 1, n1
     write(unit,*) i, inveigenvals(i)
  end do
  close(unit)

  ! Checking for negative eigenvalues
  maxeigenval=maxval(abs(inveigenvals))
  if (minval(inveigenvals) <-1d-12 * maxeigenval) then
     if (myid==root) then 
        write(*,*) minval(inveigenvals), '= minval(inveigenvals)'
        write(*,*) 'Covariance matrix is not positive definite!'
     end if
     call mpi_finalize(ierr)
     stop
  end if

  ! Writing eigenvectors for inv final cov matrix to file
  outfile = trim(outprefix) // '_inv_N_eigen.unf'
  open(unit, file=trim(outfile), form='unformatted')
  int_in = n1
  write(unit) int_in
  write(unit) eigenvectors
  write(unit) inveigenvals
  close(unit)

  ! finding eigenvals for noninverse cov matrix
  maxeigenval = maxval(abs(inveigenvals))
  write(*,*) maxeigenval, '= maxeigenval inverse covariance matrix 1'
  do i = 1, n1
     if (inveigenvals(i) > 1d-12 * maxeigenval) then
        eigenvals(i) = 1/inveigenvals(i)
     else
        write(*,*) i, inveigenvals(i), 'removed'
        eigenvals(i) = 0.d0 !1d8*maxeigenval
        inveigenvals(i) = 0.d0
     end if
  end do
  ! finding map
  allocate(map1(n1))
  map1 = matmul(eigenvectors, matmul(transpose(eigenvectors), rhs1)*eigenvals)

  deallocate(eigenvals)
  deallocate(inveigenvals)
  deallocate(eigenvectors)
  deallocate(rhs1)

  ! Write the result map to file
  allocate(outmap(0:npix1-1,nmaps))
  outmap = -1.6375d30
  if (temperature) then
     do i = 0, npix1-1
        if (mask1(i) /= -1) outmap(i,1) = map1(mask1(i))
     end do
  else
     do i = 0, npix1-1
        if (mask1(i) /= -1) then
           outmap(i,2) = map1(mask1(i))
           outmap(i,3) = map1(mask1(i)+n1/2)
        end if
     end do
  end if
  outfile = trim(outprefix) // '_map.fits'
  call write_map(outmap, ordering, outfile)

  ! Write map2mask to file 
  outmap = 0
  if (temperature) then
     outmap(:,1) = mask1
  else
     outmap(:,2) = mask1
     outmap(:,3) = mask1
  end if
  ! Writing input map1 to file
  outfile = trim(outprefix) // '_map2mask.fits'
  call write_map(outmap, ordering, outfile)

  deallocate(outmap)
  deallocate(map1)
  deallocate(mask1)
end subroutine solve_equation

!-----------------------------------------------------------------------------------------------
! subroutine output_rms
!-----------------------------------------------------------------------------------------------

subroutine output_rms(unit)
  implicit none
  integer(i4b)       :: unit, mask_order, nmaps
  integer(i4b)       :: nside, order, ncomp, npix, i, j, k, l, m, n, offset, polint
  character(len=512) :: arg, maskname, covname, outname
  real(dp), dimension(:,:), allocatable :: mask, map
  real(dp), dimension(:),   pointer     :: covdiag
  real(dp)           :: badval
  logical(lgt)       :: anynull, polarization, inv

  if(iargc() < 4) call give_user_info
  call getarg(2, maskname)
  call getarg(3, covname)
  call getarg(4, outname)

  ! Get the mask and covariance matrix
  npix = getsize_fits(maskname, nmaps=nmaps, ordering=mask_order, nside=nside)
  allocate(mask(0:npix-1,nmaps))
  call read_bintab(maskname, mask, npix, nmaps, badval, anynull)
  call read_matrix_diag(covname, covdiag, n, order, polarization, inv)

  if(order /= mask_order) then
     write(*,*) "mask and covmat have different ordering!"
     stop
  end if

  if(inv) then
     write(*,*) "Need a normal covariance matrix, not an inverse one!"
!     stop
  end if

  if(polarization) then; ncomp = 2; offset=1; else; ncomp = 1; offset = 0; end if
  allocate(map(0:npix-1,3))
  map = hpx_dbadval
  do i = 1, ncomp
     do j = 0, npix-1
        if(mask(j,i+offset) <= 0) cycle
        k = nint(mask(j,i+offset))+n/2*(i-1)
        map(j,i+offset) = sqrt(covdiag(k))
     end do
  end do

  call write_map(map, order, outname)
  deallocate(mask, covdiag, map)
end subroutine

subroutine read_matrix_diag(name, output, n, ordering, polarization, inv)
  implicit none
  integer(i4b)           :: unit, n, ordering, polint, i
  character(len=*)       :: name
  real(dp), dimension(:), pointer, intent(out) :: output
  real(dp), dimension(:), allocatable          :: column
  logical(lgt)           :: polarization
  logical(lgt), optional, intent(out)          :: inv

  unit = getlun()
  open(unit, file=name, form="unformatted")
  read(unit) n
  read(unit) ordering
  read(unit) polint
  allocate(output(n))
  allocate(column(n))
  do i = 1, n
     read(unit) column
     output(i) = column(i)
  end do
  deallocate(column)
  polarization = polint > 1
  read(unit) inv
  close(unit)
end subroutine

subroutine do_sum_equations
  implicit none
  character(len=512), dimension(:), allocatable :: inames
  character(len=512)                            :: oname
  integer(i4b)                                  :: i, n

  n = iargc()-2
  if(n < 1) call give_user_info
  allocate(inames(n))
  do i = 1, n
     call getarg(i+1, inames(i))
  end do
  call getarg(n+2, oname)
  call sum_equations(inames, oname)
  deallocate(inames)
end subroutine

! Given a list of equation files, compute the union pixel set
! and the mappings from the individual pixel sets to the union set.
! Then allocate the rhs, read individual rhses and populate. After that,
! do the same for the columns of the cov.
!
! Multiresolution pixels are pruned. Another way of handling them
! would be to expand them into maximally correlated single pixels before
! adding stuff together.
!
! This way of doing this uses very little memory (as opposed to the
! naive implementation, which would load all the matrices into memory
! at the same time), but has two weaknesses:
! 1. It assumes that the pixels are sorted (which the should be)
! 2. It effectively uses lots of random access, since it reads in
!    parallel from severl files (though each file is read sequentially).
subroutine sum_equations(infiles, outfile)
  implicit none
  character(len=*) :: infiles(:), outfile
  integer(i4b)     :: i, j, k, l, m, n, minpix, maxpix, a, b, fi, nrhs
  integer(i4b)     :: nfile, ounit, ncomp
  type(pixel_set), dimension(:), allocatable :: inpix, infull
  integer(i4b),    dimension(:), allocatable :: opix, units, pixmap
  integer(i4b),    dimension(:), allocatable :: ns, pols, orders, npixs, ms, nrhss
  integer(i4b),    dimension(:), allocatable :: lastcol, map2mask
  real(dp),        dimension(:), allocatable :: irhs, orhs, icol, ocol

  nfile = size(infiles)
  maxpix = 0
  minpix = 2000000000
  allocate(inpix(nfile), infull(nfile))
  ! Get pixels form file. inpix maps from the pruned set to global pixel.
  ! infull maps from the full set to either 0 or the corresponding inpix entry
  do i = 1, nfile
     write(*,*) "Reading pixels from " // trim(infiles(i))
     call get_pixels(infiles(i), inpix(i))
     minpix = min(minpix, minval(inpix(i)%pixels))
     maxpix = max(maxpix, maxval(inpix(i)%pixels))
  end do
  ! Make the translation map from global pixel numbers to opix by
  ! first marking all used pixels, and then assigning places for them
  ! in the output list by increasing order.
  allocate(pixmap(0:maxpix-minpix))
  write(*,*) "Computing final pixel space"
  pixmap = 0
  do i = 1, nfile
     pixmap(inpix(i)%pixels-minpix) = 1
  end do
  j = 1
  do i = 0, size(pixmap)-1
     if(pixmap(i) == 0) cycle
     pixmap(i) = j
     j = j+1
  end do
  ! We now have the final pixel space:
  m = count(pixmap /= 0)
  allocate(opix(m))
  j = 1
  do i = 0, size(pixmap)-1
     if(pixmap(i) == 0) cycle
     opix(j) = i+minpix
     j = j+1
  end do

  ! Pixel stuff done, Time to do the actual work: Open all the files,
  ! read elementwise from each, combine as appropriate and write out
  write(*,*) "Writing to " // trim(outfile)
  write(*,*) "Writing header"
  allocate(units(nfile))
  do i = 1, nfile
     units(i) = getlun()
     open(units(i),file=infiles(i),form="unformatted")
  end do
  ounit = getlun()
  open(ounit,file=outfile,form="unformatted")

  allocate(ns(nfile), pols(nfile), orders(nfile), npixs(nfile), ms(nfile), nrhss(nfile))
  do i = 1, nfile
     read(units(i)) ns(i)
     read(units(i)) orders(i)
     read(units(i)) pols(i)
  end do
  call assert(all(orders == orders(1)),"Inconsistent orderings in input files!")
  call assert(all(pols   == pols(1)),  "Inconsistent polarization status in input files!")

  ncomp = pols(1)
  ms = ns/ncomp
  n  = m*ncomp
  write(ounit) n
  write(ounit) orders(1)
  write(ounit) pols(1)

  ! Get the covariance matrix: Loop through the columns of the output
  ! matrix. For each column, find the corresponding column in each
  ! inpux matrix, ASSUMING THAT THE PIXELS ARE SORTED. This low-memory
  ! implementation is critically dependent on that. Skip the column
  ! from files that don't contain it, otherwise do as for rhs.
  write(*,*) "Writing inverse covariance matrix"
  allocate(ocol(n+1))
  allocate(lastcol(nfile))
  do i = 0, ncomp-1
     lastcol = 1
     do j = 1, m
        a = j ! Column of output matrix is really j+m*i, but we treat the
              ! matrix as two halves with the same pixel pattern.
        ocol = 0
        do fi = 1, nfile
           if(lastcol(fi) > ms(fi)) cycle ! Done with this file for now
           ! Cycle if we're not ready yet (i.e. if this file doesn't contain
           ! that column).
           if(pixmap(inpix(fi)%pixels(lastcol(fi))-minpix) > a) cycle
           ! Phew! We are at the right position!
           allocate(icol(ns(fi)+1))
           read(units(fi)) icol
           lastcol(fi) = lastcol(fi)+1
           ! We need to know which column we are at in input and output for masked_add_col
           ! to be able to handle the extra diagonal element. Output column is
           ! j+m*i, input column is lastcol(fi)-1. The extra diagonal is at index col+1.
           call masked_add_col(ocol, icol, m, ms(fi), pixmap(inpix(fi)%pixels-minpix), ncomp, &
            & lastcol(fi), j+m*i+1)
           deallocate(icol)
        end do
        write(ounit) ocol
     end do
     ! Skip any left-over columns.
     do fi = 1, nfile
        do while(lastcol(fi) <= ms(fi))
           read(units(fi))
           lastcol(fi) = lastcol(fi) + 1
        end do
     end do
  end do
  deallocate(ocol)

  ! Handle the right hand side
  write(*,*) "Writing RHS"
  do i = 1, nfile
     read(units(i)) nrhss(i)
  end do
  nrhs = nrhss(1)
  call assert(all(nrhss == nrhs), "Inconsistent number of right hand sides!")
  write(ounit) nrhs
  allocate(orhs(n))
  orhs = 0
  do j = 1, nrhs
     do i = 1, nfile
        allocate(irhs(ns(i)))
        read(units(i)) irhs
        call masked_add_col(orhs, irhs, m, ms(i), pixmap(inpix(i)%pixels-minpix), ncomp,0,0)
        deallocate(irhs)
     end do
     write(ounit) orhs
  end do
  deallocate(orhs)

  write(*,*) "Writing pixel information"
  do i = 1, nfile; read(units(i)) npixs(i); end do
  call assert(all(npixs == npixs(1)), "Inconsistent npix in input files!")
  write(ounit) npixs(1)

  ! And finally output the new map2mask
  allocate(map2mask(0:npixs(1)-1))
  map2mask = -1
  do i = 1, m
     map2mask(opix(i)) = i
  end do
  write(ounit) map2mask

  close(ounit)
  do i = 1, nfile; close(units(i)); end do

  ! And free up memory
  do i = 1, nfile
     deallocate(inpix(i)%pixels)
  end do
  deallocate(units, inpix, map2mask, opix, pixmap, nrhss)
  write(*,*) "Done"
end subroutine

subroutine get_pixels(file, pixels)
  implicit none
  character(len=*) :: file
  type(pixel_set)  :: pixels
  integer(i4b)     :: i, j, k, l, m, n, npix, unit, pol, ncomp, nmap
  integer(i4b), dimension(:), allocatable :: map2mask, counts
  unit = getlun()
  open(unit,file=file, form="unformatted")
  read(unit) n
  read(unit) ! order
  read(unit) pol
  do i = 1, n; read(unit); end do ! cov
  read(unit) nmap
  do i = 1, nmap; read(unit); end do !rhs
  read(unit) npix
  allocate(map2mask(0:npix-1))
  read(unit) map2mask
  close(unit)

  ncomp = pol

  allocate(counts(n/ncomp))
  counts = 0
  do i = 0, npix-1
     if(map2mask(i) <= 0) cycle
     counts(map2mask(i)) = counts(map2mask(i)) + 1
  end do

  if(any(counts /= 1)) then
      write(*,*) "Error: Multiresolution pixels detected. Postmap's eqn_sum does not correctly handle these."
      stop
  end if

  ! Extract the pixel information
  m = count(counts == 1)
  allocate(pixels%pixels(m))
  j = 1
  do i = 0, npix-1
     if(map2mask(i) <= 0) cycle
     pixels%pixels(j) = i
     j = j+1
  end do
  deallocate(counts, map2mask)
end subroutine

subroutine masked_add_col(ocol, icol, om, im, map, ncomp, igap, ogap)
  implicit none
  integer(i4b) :: om, im, map(:), i, j, k, ip, op, ncomp,  ioff, ooff, igap, ogap
  real(dp)     :: ocol(:), icol(:)
  ! First do all the normal operations, i.e. excluding the extra diagonal element
  do j = 0, ncomp-1
     do k = 1, im
        ip = k      + im*j
        op = map(k) + om*j
        ! We will now add from logical position input(ip) to output(op).
        ! But due to the strange format, there will be a gap at ioff and ooff:
        ! between ip and ip+1 there is an extra gap.
        ioff = 0; if(igap > 0 .and. ip >= igap) ioff = 1
        ooff = 0; if(ogap > 0 .and. op >= ogap) ooff = 1
        ocol(op+ooff) = ocol(op+ooff) + icol(ip+ioff)
     end do
  end do
  ! Ok, now the extra diagonal remains. Both input and output necessarily have the
  ! diagonal, so this is unproblematic
  if(igap > 0 .and. ogap > 0) ocol(ogap) = ocol(ogap) + icol(igap)
end subroutine

subroutine eqn_compare
  implicit none
  integer(i4b)       :: unit
  character(len=512) :: file1, file2
  integer(i8b)       :: n1, n2, npix, npix2, i
  integer(i4b)       :: order, order2, pol, pol2, nc, k, nrhs1, nrhs2
  integer(i4b), allocatable, dimension(:)   :: count1, count2
  integer(i4b), allocatable, dimension(:)   :: m1, m2
  integer(i4b), allocatable, dimension(:,:) :: ind
  real(dp),     allocatable, dimension(:,:) :: icov1, icov2, rhs1in, rhs2in 
  real(dp),     allocatable, dimension(:)   :: rhs1, rhs2
  real(dp)     :: rhsdev, icovdev, rhscorr, icovcorr
  unit = getlun()
  call getarg(2, file1)
  call getarg(3, file2)
  call read_equ_set(unit, file1, n1, order,  pol,  icov1, rhs1in, npix,  m1)
  call read_equ_set(unit, file2, n2, order2, pol2, icov2, rhs2in, npix2, m2)
  nrhs1 = size(rhs1in(1,:))
  nrhs2 = size(rhs2in(1,:))
  if(order2 /= order .or. pol2 /= pol .or. npix /= npix2 .or. nrhs1 /= nrhs2) then
     write(*,*) "Basic parameter mismatch."
     write(*,fmt="(a,i5,a,i2,a,i2,a,i6)") "Set 1: n: ", n1,", order: ", &
       & order, ", pol: ", pol, ", npix: ", npix
     write(*,fmt="(a,i5,a,i2,a,i2,a,i6)") "Set 2: n: ", n2,", order: ", &
       & order2,", pol: ", pol2,", npix: ", npix2
     stop
  end if
  if(any(m1 /= m2)) then
     write(*,*) "Pixel sets differ: " // trim(itoa(count((m1 > 0 .and. m2 <= 0) .or. (m1 <= 0 .and. m2 > 0)))) // " pixels disagree"
  end if
  ! Find the set of pixels to compare. These are pixels that are in both,
  ! and that are not multiresolution.
  ! First determine which are multires.
  allocate(count1(n1), count2(n2))
  count1 = 0;       count2 = 0
  do i = 0, npix-1; if(m1(i) > 0) count1(m1(i)) = count1(m1(i)) + 1; end do
  do i = 0, npix-1; if(m2(i) > 0) count2(m2(i)) = count2(m2(i)) + 1; end do

  ! Then build up the index arrays.
  nc = 0
  do i = 0, npix-1
     if(m1(i) > 0 .and. m2(i) > 0) then
        if(count1(m1(i)) == 1 .and. count2(m2(i)) == 1) nc = nc+1
     end if
  end do
  allocate(ind(nc,2))
  k = 1
  do i = 0, npix-1
     if(m1(i) > 0 .and. m2(i) > 0) then
        if(count1(m1(i)) == 1 .and. count2(m2(i)) == 1) then
           ind(k,1) = m1(i)
           ind(k,2) = m2(i)
           k = k+1
        end if
     end if
  end do

  allocate(rhs1(n1))
  allocate(rhs2(n2))
  do i =1, nrhs1
     rhs1=rhs1in(:,i)
     rhs2=rhs2in(:,i)
     write(*,*) 'simulation:', i-1
     rhsdev  = sqrt(sum((rhs1(ind(:,1))-rhs2(ind(:,2)))**2)/nc)
     icovdev = sqrt(sum((icov1(ind(:,1),ind(:,1))-icov2(ind(:,2),ind(:,2)))**2)/nc/nc)
!     write(*,fmt="(a,e10.5,a,e10.5)") "difference rms: rhs: ", rhsdev, " icov: ", icovdev

     rhscorr = sum(rhs1(ind(:,1))*rhs2(ind(:,2)))/sqrt(sum(rhs1(ind(:,1))**2)*sum(rhs2(ind(:,2))**2))
     icovcorr = sum(icov1(ind(:,1),ind(:,1))*icov2(ind(:,2),ind(:,2)))/ &
          & sqrt(sum(icov1(ind(:,1),ind(:,1))**2)*sum(icov2(ind(:,2),ind(:,2))**2))
!     write(*,fmt="(a,e10.5,a,e10.5)") "deviation:      rhs: ", 1-rhscorr, " icov: ", 1-icovcorr
     write(*,*)
  end do

  deallocate(rhs1,rhs2,rhs1in, rhs2in)

end subroutine


  subroutine print_eigenvalues
    implicit none
    
    
    character(len=256)                    :: matrixfile
    character(len=256)                    :: outfile
    integer(i8b)                          :: n
    integer(i4b)                          :: i, ordering, polarisation
    logical(lgt)                          :: inv
    
    real(dp), allocatable, dimension(:,:) :: cov
    real(dp), allocatable, dimension(:)   :: W

    unit = getlun()
    
    ! Get parameters
    if (iargc() /= 3) then
       write(*,*) 'print_eigenvalues takes 2 parameters'
       call give_user_info
    else 
       call getarg(2, matrixfile)
       call getarg(3, outfile)
    end if
    
    ! read matrix file
    call read_covmatrix(unit, matrixfile, ordering, polarisation, cov, inv, n)

    ! Compute eigenvalues
    allocate(W(n))
    call get_eigenvalues(cov, W)

    ! Output eigenvalues
    open(unit,file=trim(outfile))
    do i = 1, n
       write(unit,*) i, W(i)
    end do
    close(unit)
    
    deallocate(cov)
    deallocate(W)
    
  end subroutine print_eigenvalues

  !-----------------------------------------------------------------------------------------------
  ! subroutine mask_edge
  !-----------------------------------------------------------------------------------------------

  subroutine mask_edge(unit)
    implicit none
    
    integer(i4b),              intent(in) :: unit
    
    character(len=256)                    :: infile, outprefix, outfile, radiusname
    character(len=3)                      :: filnum
    integer(i4b)                          :: i, j, k, pix, numpix
    integer(i4b)                          :: nmaps, nside, npix, ordering
!    integer(i4b)                          :: ndim, ierr, lmax3    
    real(dp)                              :: nullval, vec(3), radius
    logical(lgt)                          :: anynull, mismatch, inv, temperature, nocov
  
    real(dp),     allocatable, dimension(:,:)   :: mask
    integer(i4b), allocatable, dimension(:)     :: pixlist, healvec

    ! Get parameters
    if (iargc() /= 3 .and. iargc() /= 4) then
       write(*,*) 'mask_edge takes 2(3) parameters'
       call give_user_info
    else 
       call getarg(2, infile)
       call getarg(3, outprefix)
    end if
    if (iargc() == 4) then 
       call getarg(4, radiusname)
       read(radiusname,*) radius
       radius = radius*pi/180.d0 ! in radians
    else
       radius = 1.d0*pi/180.d0 ! in radians
    end if
    write(*, fmt='(a,f5.2,a)') ' Using line thickness =', radius*180.d0/pi, ' degrees'

    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab (infile, mask, npix, nmaps, nullval, anynull)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 
    write(*,*)'npix =', npix
    allocate(pixlist(1:npix))
    allocate(healvec(1:npix))

    do j = 1, nmaps
       write(*,*) j, 'of',nmaps,'nmaps'
       k = 0
       do i = 0, npix-1
          if (mask(i, nmaps)==1.d0) then
             k = k+1
             healvec(k) =  i
          end if
       end do
       write(*,*) k, '= number of masked out pixels'
       do i = 1, k
          if (ordering==1) then
             call pix2vec_ring(nside, healvec(i), vec)          
             call query_disc(nside, vec, radius, pixlist, numpix)
          else if (ordering==2) then
             call pix2vec_nest(nside, healvec(i), vec)          
             call query_disc(nside, vec, radius, pixlist, numpix, nest=1)
          end if
!          if ( count(mask(pixlist(1:numpix),j)==0.d0)>0.3*numpix) mask(healvec(i),j)=2.d0
          if ( any(mask(pixlist(1:numpix),j)==0.d0) .and.  any(mask(pixlist(1:numpix),j)==1.d0)) mask(healvec(i),j)=2.d0
       end do
       do i = 0, npix-1
          if (mask(i,j)==2.d0) then
             mask(i,j) = 0.d0
          else
             mask(i,j) = 1.d0
          end if
       end do
    end do
    where(mask==1.d0) mask=1.d0
    where(mask==0.d0) mask=0.d0

    outfile = trim(outprefix) // '_maskedge.fits'
    call write_map(mask, ordering, outfile)
    write(*,*) '* Edge of mask written to file = ', trim(outfile)  

  end subroutine mask_edge

  !-----------------------------------------------------------------------------------------------
  ! subroutine make_mask
  !-----------------------------------------------------------------------------------------------

  subroutine make_mask(unit)
    implicit none
    
    integer(i4b),              intent(in) :: unit
    
    character(len=256)                    :: outprefix, outfile
    character(len=256)                    :: raname, decname, radiusname, linename
    integer(i4b)                          :: i, j, k, pix, numpix
    integer(i4b)                          :: nmaps, nside, npix, ordering
    real(dp)                              :: nullval, vec(3), radius, ra, dec, line
    logical(lgt)                          :: anynull, fillall
  
    real(dp),     allocatable, dimension(:,:)   :: mask
    integer(i4b), allocatable, dimension(:)     :: pixlist, healvec

    ! Get parameters
    if (iargc() /= 5 .and. iargc() /= 6) then
       write(*,*) 'makemask takes 4(5) parameters'
       call give_user_info
    else 
       call getarg(2, raname)
       call getarg(3, decname)
       call getarg(4, radiusname)
       call getarg(5, outprefix)
    end if
    read(raname,*) ra
    read(decname,*) dec
    read(radiusname,*) radius
    write(*, fmt='(a,f5.2,a,f6.2,a,f6.2,a)') ' Making mask of radius', radius, &
         & ' degrees, centered around (ra =',ra, ', dec =', dec, ')'
    radius = radius*pi/180.d0 ! in radians
    ra     = ra*pi/180.d0     ! in radians
    dec    = dec*pi/180.d0    ! in radians
   
   if (iargc() == 6) then
      call getarg(6, linename)
       read(linename,*) line
       line = line*pi/180.d0 ! in radians
       fillall = .false.
       write(*, fmt='(a,f5.2,a)') ' Using line thickness =', radius*180.d0/pi, ' degrees'
    else
       fillall = .true.
    end if

    nside = 512
    npix  = 12*nside**2
    ordering = 1
    nmaps = 1
    allocate(mask(0:npix-1,nmaps))
    allocate(pixlist(0:npix))
!    allocate(healvec(1:npix))

    call ang2vec(pi/2-dec, ra, vec)
    if (ordering==1) then
       call query_disc(nside, vec, radius, pixlist, numpix)
    else if (ordering==2) then
       call query_disc(nside, vec, radius, pixlist, numpix, nest=1)
    end if
    
    mask = 0.d0
    do i = 0, npix-1
       do j = 0, numpix-1
          if (pixlist(j)==i) mask(i,:) = 1.d0
       end do
    end do
    
    if (fillall) then
       outfile = trim(outprefix) // '_mask.fits'
    else
       outfile = trim(outprefix) // '_maskedge.fits'
    end if
    call write_map(mask, ordering, outfile)
    write(*,*) '* Mask written to file = ', trim(outfile)  

  end subroutine make_mask

!-----------------------------------------------------------------------------------------------
! subroutine findpatch
!-----------------------------------------------------------------------------------------------

  subroutine findpatch
    implicit none

    character(len=256)                    :: mapfile, outprefix, outfile, nside_in, radius_in 
    character(len=256)                    :: maskfile
    integer(i4b)                          :: i, j, nside_search, npix_search, numpix
    integer(i4b)                          :: ordering, nside, nmaps, npix, nummask
    integer(i4b)                          :: ordering2, nside2, nmaps2
    real(dp)                              :: radius, vec(3), ms, colat, long
    real(dp)                              :: healnan=-1.6375d30
    
    real(dp),     allocatable, dimension(:,:) :: inmap, inmask, mask
    real(dp),     allocatable, dimension(:)   :: rms, mean
    integer(i4b), allocatable, dimension(:)   :: pixlist, minpix
    logical(lgt), allocatable, dimension(:)   :: lmask

    ! Get parameters
    if (iargc() /= 6 ) then
       write(*,*) 'findpatch takes 5 parameters: map, mask, nside for search centers, radius, outprefix'
       call give_user_info
    else 
       call getarg(2, mapfile)
       call getarg(3, maskfile)
       call getarg(4, nside_in)
       call getarg(5, radius_in)
       call getarg(6, outprefix)
    end if
    read(nside_in,*) nside_search
    read(radius_in,*) radius
    npix_search = 12*nside_search**2
    write(*, fmt='(a,f5.2,a,i4)') ' Searching for lowest rms patch using circles of radius =', &
         & radius, ' degrees, centered around each healpix center for nside = ', nside_search
    radius = radius*pi/180.d0 ! in radians
    ! Read (and allocate) map
    call read_map(inmap, ordering, mapfile, nside=nside, nmap=nmaps)
    if (myid==0) write(*,*) nside,'= nside', nmaps,'= nmaps', ordering,'= ordering'
    if (nmaps==1) then
       if (myid==0) write(*,*) 'Running at temperature data' 
    else if (nmaps==3) then
       if (myid==0) write(*,*) 'Running at polarisation data' 
    else
       write(*,*) nmaps, '= nmaps. Unknown number. Quiting'
    end if
    if (ordering==2) then
       call convert_nest2ring(nside, inmap)
       ordering = 1
       if (myid==0) write(*,*) 'Converting input map from nest to ring' 
    end if
    npix = 12*nside**2
    call read_map(inmask, ordering2, maskfile, nside=nside2, nmap=nmaps2)
    if (myid==0) write(*,*) nside2,'= nside', nmaps2,'= nmaps', ordering2,'= ordering'
    call assert(nside == nside2, "Nside mismatch between map and mask")
    call assert(nmaps == nmaps2, "Nmaps mismatch between map and mask")
    if (ordering2==2) then
       call convert_nest2ring(nside, inmask)
       ordering2 = 1
       if (myid==0) write(*,*) 'Converting input mask from nest to ring' 
    end if

    ! Get ready to work
    allocate(pixlist(0:nside**2-1))
    allocate(rms(0:npix_search-1))
    allocate(mean(0:npix_search-1))
    allocate(mask(0:npix-1,1:1))
    mask= 0.d0
    rms = -healnan
    ! Loop over search centers
    do i = 0+myid, npix_search-1, numprocs
       if (mod(i,1000) == 0) write(*,*) i, npix_search
       call pix2vec_ring(nside_search, i, vec)
       call query_disc(nside, vec, radius, pixlist, numpix)
       ! Calculate rms for each patch
       mean(i) = sum(inmap(pixlist(0:numpix-1),nmaps))/real(numpix,dp)
       ms = 0.d0
       do j = 0, numpix-1
          ms = ms + (inmap(pixlist(j),nmaps)-mean(i))**2
       end do
       rms(i) = sqrt(ms/real(numpix,dp))
!       if (any(inmask(pixlist(0:numpix-1),nmaps) == 0.d0)) rms(i) = -healnan
       if (sum(inmask(pixlist(0:numpix-1),nmaps)) < real(numpix,dp)/2.d0) rms(i) = -healnan
       if (rms(i)<200. .and. rms(i)>0.d0) mask(pixlist(0:numpix-1),1) = 1.d0
    end do
    write(*,*) 'Unmasked part in %', sum(mask)/real(npix,dp)*100.d0

    ! Write rmsmask to file
    outfile = trim(outprefix) // '_rmsmask.fits'
    call write_map(mask, ordering, outfile)
    write(*,*) '* rmsmask written to file = ', trim(outfile)  
   
    ! Write ampmask to file
    mask = 0.d0
    where (inmap < 1000.d0) mask=1.d0
    outfile = trim(outprefix) // '_ampmask.fits'
    call write_map(mask, ordering, outfile)
    write(*,*) '* ampmask written to file = ', trim(outfile)  
    deallocate(mask)

    ! Write result map to file
    allocate(mask(0:npix_search-1, 1))
    mask(:,1) = rms
    where(mask==-healnan) mask=healnan
    outfile = trim(outprefix) // '_rms.fits'
    call write_map(mask, ordering, outfile)
    write(*,*) '* rms written to file = ', trim(outfile)  
    deallocate(mask)

    !Find best patches
    allocate(lmask(0:npix_search-1))
    lmask = .true.
!    lmask(npix_search/2:npix_search-1) = .false.
!    lmask(0:npix_search/2) = .false.
    nummask = 1
    allocate(minpix(nummask))
    do i = 1, nummask
       minpix(i) = minloc(rms, 1, lmask) - 1   !OBS minloc does not take start point 0 into account!!!
       call pix2ang_ring(nside_search, minpix(i), colat, long)
       write(*,fmt='(i2,a,i5,a,i5,a,f6.2,a,f6.2,a,f7.2,a,f7.2)') i, ': pix=', minpix(i),'(',&
            & npix_search,') rms=',rms(minpix(i)),'  mean=', mean(minpix(i)), '   lat=', &
            & 90.d0-colat*180.d0/pi, '   lon=', long*180.d0/pi
       lmask(minpix(i)) = .false.
    end do

    ! Make mask
    allocate(mask(0:npix-1, 1))
    mask = 0.d0
    do i = 1, 1!nummask
       call pix2vec_ring(nside_search, minpix(i), vec)
       call query_disc(nside, vec, radius, pixlist, numpix)
       mask(pixlist(0:numpix-1),1) = 1.d0
    end do

    ! Write patch map to file
    mask = inmap*mask
    outfile = trim(outprefix) // '_patch.fits'
    call write_map(mask, ordering, outfile)
    write(*,*) '* patch map written to file = ', trim(outfile)  

    ! Clean up
    deallocate(inmap,rms,mask,lmask)
    deallocate(pixlist)

  end subroutine findpatch

!-----------------------------------------------------------------------------------------------
! subroutine give_user_info
!-----------------------------------------------------------------------------------------------

subroutine give_user_info
  implicit none

  if (myid == root) then
     write(*,*) 'Usage: mpirun -n N postmap command params'
     write(*,*) 
     write(*,*) 'command = matmul   : Matrix multiplication'
     write(*,*) ' matmul [map2mask] [matrix file] [inmap] [outmap]'
     write(*,*) 
     write(*,*) 'command = proj0    : Project out null-space'
     write(*,*) ' proj0 [map2mask] [inv N file] [inmap] [outmap]'
     write(*,*) 
     write(*,*) 'command = finalmap : Make final map, cov. matrix etc from 2 subsets'
     write(*,*) ' finalmap [eqns.unf 1] [eqns.unf 2] [outprefix]'
     write(*,*) 
     write(*,*) 'command = solve    : Solve equation system as output from tod2map.'
     write(*,*) ' solve eqnsystem.unf outprefix'
     write(*,*) 
     write(*,*) 'command = cutoff   : Project out eigenmodes lower than cutoff'
     write(*,*) ' cutoff [map2mask] [inv_N_eigen] [inmap] [outprefix] cutoff'
     write(*,*) 
     write(*,*) 'command = polcut   : Project out eigenmodes lower than cutoff for pol'
     write(*,*) ' polcut [map2mask] [inv_N_eigen] [inmap] [outprefix] cutoff'
     write(*,*)    
     write(*,*) 'command = lcut     : Project out eigenmodes lower than lmax'
     write(*,*) ' lcut map2mask inv_N.unf inmap outprefix lmax (nocov/onlymap)'
     write(*,*)    
     write(*,*) 'command = rms      : Produce rms map from covar and map2mask'
     write(*,*) ' rms map2mask N.unf outname'
     write(*,*)    
     write(*,*) 'command = eqn_sum  : Sum several equation sets. Need not have same pixel set.'
     write(*,*) ' eqn_sum in1.unf [in2.unf [in3.unf ...]] out.unf'
     write(*,*)    
     write(*,*) 'command = kernel  : Output coupling kernel matrix for given mask'
     write(*,*) ' kernel mask.fits outprefix lmax'
     write(*,*)
     write(*,*) 'command = eqn_compare  : Compare two equation sets'
     write(*,*) ' eqn_compare eqn1.unf eqn2.unf'
     write(*,*)     
     write(*,*) 'OBS: Only postmap lcut, rms and eqn_sum are known to work and still usefull.'
     write(*,*) '     Finalmap and solve are replaced by scalapost versions'
     write(*,*) '     Eqn_compare might need to be updated to new eq set format'
  end if
  call mpi_finalize(ierr)
  stop

end subroutine give_user_info


end program postmap
