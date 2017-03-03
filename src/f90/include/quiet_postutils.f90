module quiet_postutils
  use healpix_types
  use math_tools
  use quiet_fileutils
  use pix_tools
  use alm_tools
  implicit none

contains

  !---------------------------------------------------------------------
  !  putting the unwanted modes into fmap 
  !-----------------------------------------------------------------------
  subroutine get_lowmodes(fmap, lmax, n, nmodes, nmaps, nside, ordering, temperature, map2mask)
    implicit none

    integer(i4b),                          intent(in) :: lmax, nmaps, nside, ordering
    integer(i8b),                          intent(in) :: n, nmodes
    logical(lgt),                          intent(in) :: temperature
    real(dp), dimension(:,:), allocatable, intent(in) :: map2mask
    
    integer(i4b)                                :: i, j, k, l, m, p, pix, npix
    real(dp)                                    :: dummy
    real(dp),     dimension(:,:),   allocatable :: fmap, map
    complex(dpc), dimension(:,:,:), allocatable :: alm

    npix = 12*nside**2
    allocate(map(0:npix-1,nmaps)) 
    if (temperature) then
       i=0
       do l = 0, lmax
          write(*,*) 'l =',l
          allocate(alm(nmaps,0:l,0:l))
          do m = 0, l
             if (m==0) then 
                alm = 0.d0
                i = i+1
                alm(1,l,m) = cmplx(1.d0, 0.d0)
                call alm2map(nside, l, l, alm, map(:,1))
                if (ordering ==2) then                       ! alm2map always gives out ringed
                   call convert_ring2nest(nside ,map)
                end if
                do j = 0, npix -1
                   dummy = map2mask(j,1)
                   if (dummy /= -1.d0) then
                      pix = int(dummy,i4b)
                      fmap(pix, i) = map(j,1)   ! obs temp
                   end if
                end do
             else 
                alm = 0.d0
                i = i+1
                alm(1,l,m) = cmplx(1.d0, 0.d0)
                call alm2map(nside, l, l, alm, map(:,1))
                if (ordering ==2) then
                   call convert_ring2nest(nside ,map)
                end if
                do j = 0, npix -1
                   dummy = map2mask(j,1)
                   if (dummy /= -1.d0) fmap(int(dummy,i4b), i) = map(j,1)   ! obs temp
                end do
                alm = 0.d0
                i = i+1
                alm(1,l,m) = cmplx(0.d0 ,1.d0)
                call alm2map(nside, l, l, alm, map(:,1))
                if (ordering ==2) then
                   call convert_ring2nest(nside ,map)
                end if
                do j = 0, npix -1
                   dummy = map2mask(j,1)
                   if (dummy /= -1.d0) fmap(int(dummy,i4b), i) = map(j,1)   ! obs temp
                end do
             end if
          end do
          deallocate(alm)
       end do
    else ! if polarisation
       i=0
       fmap = 0.d0
       do p=2,3
          do l = 0, lmax
             write(*,*) 'l =',l
             allocate(alm(nmaps,0:l,0:l))
             do m = 0, l
                if (m==0) then
                   alm = 0.d0
                   i = i+1
                   alm(p,l,m) = cmplx(1.d0, 0.d0)
                   call alm2map(nside, l, l, alm, map)
                   if (ordering ==2) then                       ! alm2mask always gives out ringed
                      call convert_ring2nest(nside ,map)
                   end if
                   do j = 0, npix -1
                      dummy = map2mask(j,2)
                      if (dummy /= -1.d0) then
                         pix = int(dummy,i4b)
                         fmap(pix, i) = map(j,2)
                         fmap(pix+n/2, i) = map(j,3)
                      end if
                   end do
                else 
                   alm = 0.d0
                   i = i+1
                   alm(p,l,m) = cmplx(1.d0, 0.d0)
                   call alm2map(nside, l, l, alm, map)
                   if (ordering ==2) then
                      call convert_ring2nest(nside ,map)
                   end if
                   do j = 0, npix -1
                      dummy = map2mask(j,2)
                      if (dummy /= -1.d0) then
                         pix = int(dummy,i4b)
                         fmap(pix, i) = map(j,2)
                         fmap(pix+n/2, i) = map(j,3)
                      end if
                   end do
                   alm = 0.d0
                   i = i+1
                   alm(p,l,m) = cmplx(0.d0 ,1.d0)
                   call alm2map(nside, l, l, alm, map)
                   if (ordering ==2) then
                      call convert_ring2nest(nside ,map)
                   end if
                   do j = 0, npix -1
                      dummy = map2mask(j,2)
                      if (dummy /= -1.d0) then
                         pix = int(dummy,i4b)
                         fmap(pix, i) = map(j,2)
                         fmap(pix+n/2, i) = map(j,3)
                      end if
                   end do
                end if
             end do
             deallocate(alm)
          end do
       end do
    end if
    deallocate(map)
    
  end subroutine get_lowmodes
  
  !---------------------------------------------------------------------
  ! Check that input data ere consistent
  !-----------------------------------------------------------------------
  subroutine check_oneinput(myid, ordering, pol)
    implicit none

    integer(i4b),   intent(in) :: myid, ordering, pol

    integer(i4b)               :: ierr, root=0


    if (ordering /= 1 .and. ordering /= 2) then
       if (myid==root) write(*,*) "Ordering is neither ring nor nested. Quiting"
       call mpi_finalize(ierr)
       stop
    end if

    if (pol /= 1 .and. pol /= 2) then
       if (myid==root) write(*,*) "Neither temperature nor polarisation. Quiting"
       call mpi_finalize(ierr)
       stop
    end if

  end subroutine check_oneinput

  !---------------------------------------------------------------------
  ! Check that input data ere consistent
  !-----------------------------------------------------------------------
  subroutine check_input(myid, ordering1, ordering2, pol1, pol2, npix1, npix2, temperature)
    implicit none

    integer(i4b)               :: myid, ordering1, ordering2, pol1, pol2, npix1, npix2
    logical(lgt)               :: temperature

    integer(i4b)               :: n1, n2, ierr, root=0
    logical(lgt)               :: mismatch, temperature2

    ! Check that we have temperature or polarisation data
    if (pol1 == 1) then
       temperature = .true.
    else if (pol1 == 2) then
       temperature = .false.
    else
       if (myid==root) then 
          write(*,*) "Input1 is neither temperature nor polarisation data. Quiting"
       end if
       call mpi_finalize(ierr)
       stop
    end if

    if (pol2 == 1) then
       temperature2 = .true.
    else if (pol2 == 2) then
       temperature2 = .false.
    else
       if (myid==root) then 
          write(*,*) "Imput2 is neither temperature nor polarisation data. Quiting"
       end if
       call mpi_finalize(ierr)
       stop
    end if

    ! check that both are temperature or both are polarisation  
    if (temperature .neqv. temperature2) then
       if (myid==root) then 
          write(*,*) 'Both maps must be either temperature or polarization. Quiting'
       end if
       call mpi_finalize(ierr)
       stop
    end if

    ! Check that both input maps is of same pixel size
    if (npix1    /= npix2) then
       if (myid==root) then 
          write(*,*) 'Input map_1 and map_2 is not of same resolution. Quiting'
       end if
       call mpi_finalize(ierr)
       stop
    end if

    ! Finished and everything ok, writing message
    if (myid == 0) then
       if (temperature) then
          write(*,*) 'Running at temperature data'
       else
          write(*,*) 'Running at polarisation data'
       end if
    end if

  end subroutine check_input

  !---------------------------------------------------------------------
  ! Check that matrix is symmetric
  !-----------------------------------------------------------------------
  subroutine check_sym_matrix(matrix, asym)
    implicit none

    real(dp), dimension(1:,1:) :: matrix
    logical(lgt)               :: asym

    integer(i4b)               :: n, i, j

    ! Checking that input matrices are symmetric
    asym = .false.
    n = size(matrix(1,:))
    do i = 1, n
       do j = i, n
          if ( abs((matrix(i,j)-matrix(j,i))/(matrix(i,j)+matrix(j,i))) > 1d-8 ) then 
             write(*,*) i,j, matrix(i,j), matrix(j,i)
             write(*,*) i,j, matrix(i,i), matrix(j,j)
             asym = .true.
          end if
       end do
    end do

  end subroutine check_sym_matrix

  !------------------------------------------------------------------------------
  ! solving equationset of inv cov and rhs to find map and cov
  !------------------------------------------------------------------------------
  subroutine solve_eq(myid, unit, outprefix, outtext, matrix, rhs, map, cov)
    implicit none

    integer(i4b),               intent(in)    :: myid, unit
    character(len=*),           intent(in)    :: outprefix
    character(len=*),           intent(inout) :: outtext
    real(dp), dimension(1:,1:), intent(inout) :: matrix    !invcov, is overwritten
    real(dp), dimension(1:),    intent(in)    :: rhs
    real(dp), dimension(1:,1:), intent(out)   :: cov
    real(dp), dimension(1:),    intent(out)   :: map
    
    integer(i4b)                              :: n, i, j
    real(dp)                                  :: t1, t2    
    real(dp), allocatable, dimension(:)       :: eigenvals
    real(dp), allocatable, dimension(:,:)     :: eigenvectors

    n = size(rhs)
    allocate(eigenvals(n))
    allocate(eigenvectors(n,n))

    ! decompose inverse covariance matrix into eigenvectors
    call get_eigen_decomposition(myid, matrix, eigenvals, eigenvectors)
    ! check and invert eigenvalues
    outtext = 'inveigenval_' // outtext
    call invert_eigenvals(myid, eigenvals, unit, outprefix, outtext)
    
    ! finding map
    map = matmul(eigenvectors, matmul(transpose(eigenvectors), rhs)*eigenvals)
    
    ! finding non-inverse covariance matrix
    matrix = transpose(eigenvectors)
    do i = 1, n
       matrix(i,:) = matrix(i,:) * eigenvals(i)
    end do
    call cpu_time(t1)  
    call DGEMM('N', 'N', n, n, n, 1.d0, eigenvectors, n, matrix, n, 0.d0, cov, n)
    call cpu_time(t2)
    write(*,*) 'Time spent on matrix multiplication:', t2-t1 

    ! Clean up
    deallocate(eigenvals)
    deallocate(eigenvectors)

  end subroutine solve_eq

  !----------------------------------------------------------------
  ! checks and inverts eigenvalues
  !----------------------------------------------------------------
  subroutine invert_eigenvals(myid, eigenvals, unit, outprefix, outtext, kept, acc, extravec)
    implicit none

    integer(i4b),                    intent(in)    :: myid
    real(dp),         dimension(1:), intent(inout) :: eigenvals
    integer(i4b),          optional, intent(in)    :: unit   
    character(len=*),      optional, intent(in)    :: outtext, outprefix
    integer(i4b),          optional, intent(out)   :: kept
    real(dp),              optional, intent(in)    :: acc   
    real(dp),              optional, dimension(:)  :: extravec

    integer(i4b)            :: n, i, ierr
    character(len=256)      :: outfile
    real(dp)                :: maxeigenval

    n = size(eigenvals)

    ! Output inveigenvals input invcov
    if (present(unit) .and. present(outprefix) .and. present(outtext)) then
       outfile = trim(outprefix) // '_' // trim(outtext) // '_before.txt'
       open(unit, file=trim(outfile))
       do i = 1, n
          write(unit,*) i, eigenvals(i)
       end do
       close(unit)
    end if

    ! Checking for negative eigenvalues
    if(present(outtext)) then
       call check_eigenvals(myid, outtext, eigenvals, kept, acc, extravec=extravec)
    else
       call check_eigenvals(myid, "", eigenvals, kept, acc, extravec=extravec)
    end if

    ! inverting eigenvals
    maxeigenval = maxval(abs(eigenvals))
    do i = 1, n
       if (eigenvals(i) > 1d-12 * maxeigenval) then
          eigenvals(i) = 1/eigenvals(i)
       else
          eigenvals(i)    = 0.d0 
       end if
    end do

    ! Output inveigenvals input non-inv covariance matrix
    if (present(unit) .and. present(outprefix) .and. present(outtext)) then
       outfile = trim(outprefix) // '_inv_' //trim(outtext) // '_after.txt'
       open(unit, file=trim(outfile))
       do i = 1, n
          write(unit,*) i, eigenvals(i)
       end do
       close(unit)
    end if

  end subroutine invert_eigenvals

  !----------------------------------------------------------------
  ! checks and inverts and take sqrt of eigenvalues
  !----------------------------------------------------------------
  subroutine modify_eigenvals(myid, eigenvals, invei, sqrtei, sqrtinvei, unit, outprefix, outtext, kept, acc, cov, ignore)
    implicit none

    integer(i4b),                    intent(in)    :: myid
    real(dp),         dimension(1:), intent(inout) :: eigenvals, invei, sqrtei, sqrtinvei
    integer(i4b),          optional, intent(in)    :: unit   
    character(len=*),      optional, intent(in)    :: outtext, outprefix
    integer(i4b),          optional, intent(out)   :: kept
    real(dp),              optional, intent(in)    :: acc   
    logical(lgt),          optional, intent(in)    :: cov   
    integer(i4b),          optional, intent(in)    :: ignore

    integer(i4b)            :: n, i, ierr
    character(len=256)      :: outfile
    real(dp)                :: maxeigenval, accuracy

    n = size(eigenvals)
    if (present(acc)) then
       accuracy = acc
    else
       accuracy = 1d-12
    end if

    ! Output inveigenvals
    if (present(unit) .and. present(outprefix) .and. present(outtext)) then
       outfile = trim(outprefix) // '_' // trim(outtext) // 'eigenvals_before.txt'
       open(unit, file=trim(outfile))
       do i = 1, n
          write(unit,*) i, eigenvals(i)
       end do
       close(unit)
    end if

    ! Checking for negative eigenvalues and below machine precistion
    if(present(outtext)) then
       call check_eigenvals(myid, outtext, eigenvals, kept, acc, cov=cov, ignore=ignore)
    else
       call check_eigenvals(myid, "", eigenvals, kept, acc, cov=cov, ignore=ignore)
    end if

    ! inverting eigenvals etc
    do i = 1, n
       if (eigenvals(i) > 0.d0) then
          sqrtei(i)    = sqrt(eigenvals(i))
          invei(i)     = 1/eigenvals(i)
          sqrtinvei(i) = sqrt(invei(i))
       else
          eigenvals(i) = 0.d0 
          sqrtei(i)    = 0.d0
          invei(i)     = 0.d0
          sqrtinvei(i) = 0.d0
       end if
    end do

    ! Output inveigenvals
    if (present(unit) .and. present(outprefix) .and. present(outtext)) then
       outfile = trim(outprefix) // '_' //trim(outtext) // 'eigenvals__after.txt'
       open(unit, file=trim(outfile))
       do i = 1, n
          write(unit,*) i, eigenvals(i)
       end do
       close(unit)
    end if

  end subroutine modify_eigenvals

  !----------------------------------------------------------------
  ! checks and inverts and take sqrt of eigenvalues
  !----------------------------------------------------------------
  subroutine invert_andsqrt_eigenvals(myid, eigenvals, sqrteigenvals, sqrt_input, unit, outprefix, outtext, sqrteigenvals2)
    implicit none

    integer(i4b),                         intent(in)    :: myid
    real(dp),              dimension(1:), intent(inout) :: eigenvals
    real(dp),         optional,     dimension(1:), intent(out)   :: sqrteigenvals, sqrteigenvals2
    logical(lgt),     optional,                    intent(in)    :: sqrt_input
    integer(i4b),          optional,      intent(in)    :: unit   
    character(len=*),      optional,      intent(in)    :: outtext, outprefix

    integer(i4b)            :: n, i, ierr
    character(len=256)      :: outfile
    real(dp)                :: maxeigenval

    n = size(eigenvals)

    ! Output input inveigenvals 
    if (present(unit) .and. present(outprefix) .and. present(outtext)) then
       outfile = trim(outprefix) // '_' // trim(outtext) // '_before.txt'
       open(unit, file=trim(outfile))
       do i = 1, n
          write(unit,*) i, eigenvals(i)
       end do
       close(unit)
    end if

    ! Checking for negative eigenvalues
    if(present(outtext)) then
       call check_eigenvals(myid, outtext, eigenvals)
    else
       call check_eigenvals(myid, "", eigenvals)
    end if

    ! inverting eigenvals and taking sqrt
    maxeigenval = maxval(abs(eigenvals))
    do i = 1, n
       if (eigenvals(i) > 1d-12 * maxeigenval) then
          if (sqrt_input) sqrteigenvals(i) = sqrt(eigenvals(i))
          if (present(sqrteigenvals2) .and. .not. sqrt_input) sqrteigenvals2(i) = sqrt(eigenvals(i))
          eigenvals(i)  = 1/eigenvals(i)
          if (.not. sqrt_input) sqrteigenvals(i) = sqrt(eigenvals(i))
          if (present(sqrteigenvals2) .and. sqrt_input) sqrteigenvals2(i) = sqrt(eigenvals(i))
       else
!          write(*,*) i, eigenvals(i), 'removed'
          eigenvals(i)     = 0.d0 
          sqrteigenvals(i) = 0.d0 
          if (present(sqrteigenvals2)) sqrteigenvals2(i) = 0.d0 
       end if
    end do

    ! Output inverted eigenvals
    if (present(unit) .and. present(outprefix) .and. present(outtext)) then
       outfile = trim(outprefix) // '_inv_' //trim(outtext) // '_after.txt'
       open(unit, file=trim(outfile))
       do i = 1, n
          write(unit,*) i, eigenvals(i)
       end do
       close(unit)
    end if

  end subroutine invert_andsqrt_eigenvals

  !----------------------------------------------------------------
  ! checks eigenvalues
  !----------------------------------------------------------------
  subroutine check_eigenvals(myid, outtext, eigenvals, kept, acc, extravec, cov, ignore)
    implicit none

    integer(i4b),                         intent(in)    :: myid
    character(len=*),                     intent(in)    :: outtext
    real(dp),              dimension(1:), intent(inout) :: eigenvals
    integer(i4b),               optional, intent(out)   :: kept
    real(dp),                   optional, intent(in)    :: acc   
    real(dp),   optional,  dimension(1:), intent(inout) :: extravec
    logical(lgt),               optional, intent(in)    :: cov   
    integer(i4b),               optional, intent(in)   :: ignore

    integer(i4b)             :: n, i, j, k, ierr, ig
    real(dp)                 :: maxeigenval, mineigenval, accuracy
    real(dp),    allocatable :: copy(:)

    n = size(eigenvals)
    if (present(acc)) then
       accuracy = acc
    else
       accuracy = 1d-12
    end if
    if (present(ignore)) then
       ig = ignore
    else
       ig = 0
    end if

    if (myid==0) then 
       write(*,*) minval(eigenvals(ig+1:n)), '= mineigenval ', trim(outtext)
       write(*,*) maxval(eigenvals(ig+1:n)), '= maxeigenval ', trim(outtext)
    end if

    if (maxval(eigenvals(ig+1:n)) /= maxval(abs(eigenvals(ig+1:n)))) then
       if (myid==0) write(*,*) trim(outtext), 'Covariance matrix is not positive definite. Quiting'
       !       call mpi_finalize(ierr)
       stop
    end if

    ! Putting to zero numbers smaller than accuarcy
    j = 0
    if (present(cov) .and. cov) then   
       allocate(copy(n-ig))
       copy = eigenvals(ig+1:n)
       k = 0
       do i = 1, n-ig
          if (copy(i) < 0.d0) then
             copy(i) = 1d+30
             k = k+1
          end if
       end do
       if (myid==0) write(*,*) trim(outtext), ':', k,'negative eigenvalues put to zero'       
       mineigenval = minval(copy)
       do i = 1, ig
          eigenvals(i) = 0.d0
          if (present(extravec)) extravec(i) = 0.d0
       end do
       do i = ig+1, n
          if (eigenvals(i) >= mineigenval/accuracy .or. eigenvals(i) < 0.d0) then
             eigenvals(i) = 0.d0
             if (present(extravec)) extravec(i) = 0.d0
             j = j+1
          end if
       end do
       if (j>k) then
          if (myid==0) write(*,*) trim(outtext), ':', j-k,'LARGE eigvals put to zero for cov'
       end if
    else ! if inv cov
       do i = 1, ig
          eigenvals(i) = 0.d0
          if (present(extravec)) extravec(i) = 0.d0
       end do
       maxeigenval = maxval(eigenvals)
       do i = ig+1, n
          if (eigenvals(i) <= accuracy * maxeigenval) then
             eigenvals(i) = 0.d0
             if (present(extravec)) extravec(i) = 0.d0
             j = j+1
          end if
       end do
    end if

    if (myid==0) write(*,*) trim(outtext), ':', j+ig,'eigenvalues put to zero'
    if (present(kept)) kept = n-j-ig

    ! Checking for negative eigenvalues
    if (minval(eigenvals) < -1d-8 ) then
       if (myid==0) write(*,*) trim(outtext), 'Covariance matrix is not positive definite. Quiting'
       !       call mpi_finalize(ierr)
       stop
    end if

  end subroutine check_eigenvals

  !----------------------------------------------------------------
  ! checks and take sqrt of eigenvalues
  !----------------------------------------------------------------
  subroutine sqrt_eigenvals(myid, outtext, eigenvals)
    implicit none

    integer(i4b),                           intent(in)  :: myid
    character(len=*),                       intent(in)  :: outtext
    real(dp),              dimension(1:), intent(inout) :: eigenvals

    integer(i4b)            :: n, i, ierr
    real(dp)                :: maxeigenval

    call check_eigenvals(myid, outtext, eigenvals)

    eigenvals = sqrt(eigenvals)

  end subroutine sqrt_eigenvals

  !----------------------------------------------------------------
  ! writes maps to file
  !----------------------------------------------------------------

  subroutine write_maps(unit, map, map2heal, pol, npix, ordering, outprefix, simprefix, maptext)
    implicit none

    character(len=*),                   intent(in)    :: outprefix, simprefix, maptext
    real(dp),         dimension(1:,1:), intent(in)    :: map 
    integer(i4b),     dimension(1:),    intent(in)    :: map2heal
    integer(i4b),                       intent(in)    :: unit, pol, npix, ordering

    real(dp), allocatable, dimension(:,:) :: outmap
    character(len=512)                    :: outfile
    character(len=3)                      :: filnum
    integer(i4b)                          :: i, nmaps, nrhs, n
    
    n = size(map(:,1))
    nrhs = size(map(1,:))
    if (pol == 1) then 
       nmaps = 1
    else if (pol == 2) then
       nmaps = 3
    end if

    allocate(outmap(0:npix-1, nmaps))
    call map2outmap(outmap, map(:,1), map2heal, pol)
    outfile = trim(outprefix) // '_' // trim(maptext)
    call write_map(outmap, ordering, trim(outfile))
    write(*,*) '* Written to file = ', trim(outfile)  
    if (nrhs > 1) then
       do i = 1, nrhs-1
          call int2string(i, filnum)
          call map2outmap(outmap, map(:,i+1), map2heal, pol)
          outfile = trim(simprefix)//'_sim'//filnum// '_' // trim(maptext)
          call write_map(outmap, ordering, outfile)
       end do
    end if
    deallocate(outmap)

  end subroutine write_maps


  !----------------------------------------------------------------
  ! writes maps to file
  !----------------------------------------------------------------

  subroutine write_maps_m2m(unit, map, map2mask, pol, npix, ordering, outprefix, simprefix, maptext)
    implicit none

    character(len=*),                   intent(in)    :: outprefix, simprefix, maptext
    real(dp),         dimension(1:,1:), intent(in)    :: map 
    integer(i4b),     dimension(0:),    intent(in)    :: map2mask
    integer(i4b),                       intent(in)    :: unit, pol, npix, ordering

    real(dp), allocatable, dimension(:,:) :: outmap
    character(len=512)                    :: outfile
    character(len=3)                      :: filnum
    integer(i4b)                          :: i, nmaps, nrhs, n
    
    n = size(map(:,1))
    nrhs = size(map(1,:))
    if (pol == 1) then 
       nmaps = 1
    else if (pol == 2) then
       nmaps = 3
    end if

    allocate(outmap(0:npix-1, nmaps))
    call map2healmap(outmap, map(:,1), map2mask, pol, npix)
    outfile = trim(outprefix) // '_' // trim(maptext)
    call write_map(outmap, ordering, outfile)
    write(*,*) '* Written to file = ', trim(outfile)  
    if (nrhs > 1) then
       do i = 1, nrhs-1
          call int2string(i, filnum)
          call map2healmap(outmap, map(:,i+1), map2mask, pol, npix)
          outfile = trim(simprefix)//'_sim'//filnum// '_' // trim(maptext)
          call write_map(outmap, ordering, outfile)
       end do
    end if
    deallocate(outmap)

  end subroutine write_maps_m2m

  !----------------------------------------------------------------
  ! writes chi_squares for lcut to log and screen
  !----------------------------------------------------------------

  subroutine write_chi_log(unit, chi_square, n, outprefix, helptext)
    implicit none

    character(len=*),           intent(in)    :: outprefix, helptext
    real(dp),                   intent(in)    :: chi_square
    integer(i4b),               intent(in)    :: unit, n

    character(len=512)                        :: outfile

    ! Write to screen
    write(*,*)'--->', chi_square,'= chi_sq ',trim(helptext)
    write(*,*)' ', n,'     = number of modes'
    write(*,*)'    ', chi_square/(n),'= chi_sq per mode'
    write(*,*)'    ',(chi_square-n)/sqrt(2.d0*(n)),'= sigma'
    
    ! Write to logfile
    outfile = trim(outprefix) // '_log.txt'
    if (trim(helptext)=='before') open(unit+1, file=trim(outfile))
    write(unit+1,*)'--->', chi_square,'= chi_sq ',trim(helptext)
    write(unit+1,*)' ', n,'     = number of modes'
    write(unit+1,*)'    ', chi_square/(n),'= chi_sq per mode'
    write(unit+1,*)'    ',(chi_square-n)/sqrt(2.d0*(n)),'= sigma'
    write(unit+1,*)
    if (trim(helptext)=='after') close(unit+1)

  end subroutine write_chi_log

  !----------------------------------------------------------------
  ! writes chi_squares to file and screen
  !----------------------------------------------------------------

  subroutine write_chis(unit, chi_square, n, outprefix, simprefix, maptext, helptext)
    implicit none

    character(len=*),                intent(in)    :: outprefix, simprefix, maptext
    real(dp),         dimension(1:), intent(in)    :: chi_square
    integer(i4b),                    intent(in)    :: unit, n
    character(len=*),      optional, intent(in)    :: helptext

    character(len=512)                    :: outfile, writetext
    character(len=3)                      :: filnum
    integer(i4b)                          :: i, nmaps, nrhs

    if (present(helptext)) then
       writetext = helptext
    else
       writetext = ''
    end if

    write(*,*)'----->', chi_square(1),'= chi_sq ', trim(writetext)
    write(*,fmt='(a,i20,a)')'      ', n,'      = number of modes'
    write(*,*)'      ', chi_square(1)/(n),'= chi_sq per mode'
    write(*,*)'      ',(chi_square(1)-n)/sqrt(2.d0*(n)),'= sigma'   

    outfile = trim(outprefix) // '_' // trim(maptext) // '_log.txt'    
    open(unit, file=trim(outfile))
    do i = 1, size(chi_square)
       if (i==1) then
          write(unit,*) 'Actual ', trim(writetext) 
       else
          write(unit,*) 'Simulation number', i-1
       end if
       write(unit,*)'----->', chi_square(i),'= chi_sq ', trim(writetext)
       write(unit,fmt='(a,i20,a)')'      ', n,'      = number of modes'
       write(unit,*)'      ', chi_square(i)/(n),'= chi_sq per mode'
       write(unit,*)'      ',(chi_square(i)-n)/sqrt(2.d0*(n)),'= sigma'   
       write(unit,*)
    end do
    close(unit)

    if (size(chi_square)>1) then
       outfile = trim(simprefix) // '_' // trim(maptext) // '_sims.txt'
       open(unit, file=trim(outfile))
       do i = 2, size(chi_square)
          write(unit,*) i-1, (chi_square(i)-n)/sqrt(2.d0*(n)), chi_square(i)
       end do
       close(unit)
    end if

  end subroutine write_chis

  !----------------------------------------------------------------
  ! writes a healpix map to a n-map using map2mask
  !----------------------------------------------------------------

  subroutine healmap2map(nmap, inmap, map2mask, pol)
    implicit none

    real(dp),       dimension(1:),    intent(out) :: nmap
    real(dp),       dimension(0:,1:), intent(in)  :: inmap, map2mask
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)                                  :: j, npix, n
    real(dp)                                      :: dummy

    n = size(nmap)
    npix = size(inmap(:,1))
    do j = 0, npix -1
       if (pol==1) then
          dummy = map2mask(j,1)
          if (dummy /= -1.d0) nmap(int(dummy,i4b)) = inmap(j,1)   ! obs temp
       else if (pol==2) then
          dummy = map2mask(j,2)
          if (dummy /= -1.d0) then
             nmap(int(dummy,i4b)) = inmap(j,2)  
             nmap(int(dummy,i4b)+n/2) = inmap(j,3)  
          end if
       end if
    end do

  end subroutine healmap2map

  !----------------------------------------------------------------
  ! writes a n-map to a healpix outmap using map2mask
  !----------------------------------------------------------------

  subroutine map2healmap(outmap, map, map2mask, pol, npix)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    real(dp),       dimension(1:),    intent(in)  :: map
    integer(i4b),   dimension(0:),    intent(in)  :: map2mask
    integer(i4b),                     intent(in)  :: pol, npix
  
    integer(i4b)                                  :: i, n

    n = size(map) 
    
    outmap = -1.6375d30
    if (pol==1) then
       do i = 0, npix-1
          if (map2mask(i) > 0) outmap(i,1) = map(map2mask(i))
       end do
    else if (pol==2) then
       do i = 0, npix-1
          if (map2mask(i) > 0) then
             outmap(i,2) = map(map2mask(i))
             outmap(i,3) = map(map2mask(i)+n/2)
          end if
       end do
    end if

  end subroutine map2healmap

  !----------------------------------------------------------------
  ! writes a n-map to a healpix outmap using map2heal-vector
  !----------------------------------------------------------------
  subroutine map2outmap(outmap, map, map2heal, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    real(dp),       dimension(1:),    intent(in)  :: map
    integer(i4b),   dimension(1:),    intent(in)  :: map2heal
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)                                  :: n

    n = size(map)
    outmap = -1.6375d30
    if (pol==1) then
       outmap(map2heal,1) = map
    else if (pol==2) then
       outmap(map2heal,2) = map(1:n/2)
       outmap(map2heal,3) = map(n/2+1:n)
    end if

  end subroutine map2outmap

  !----------------------------------------------------------------
  ! writes a multidim n-map to a healpix outmap using map2heal-vector
  !----------------------------------------------------------------
  subroutine nmap2outmap(outmap, map, map2heal)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    real(dp),       dimension(1:,1:), intent(in)  :: map
    integer(i4b),   dimension(1:),    intent(in)  :: map2heal
    integer(i4b)                                  :: i

    outmap = -1.6375d30
    do i = 1, size(map(1,:))
       outmap(map2heal,i) = map(:,i)
    end do

  end subroutine nmap2outmap

  !----------------------------------------------------------------
  ! writes a map2mask to a healpix outmap using map2heal-vector
  !----------------------------------------------------------------
  subroutine map2mask2outmap(outmap, map2heal, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    integer(i4b),   dimension(1:),    intent(in)  :: map2heal
    integer(i4b),                     intent(in)  :: pol
    integer(i4b)                                  :: i, j

    !Writing map2mask to outmap
    outmap = -1.d0
    do i = 1, size(outmap,2)
       do j = 1, size(map2heal)
          outmap(map2heal(j),i) = j
       end do
    end do
  end subroutine map2mask2outmap

  !----------------------------------------------------------------
  ! writes a mask to a healpix outmap using map2heal-vector
  !----------------------------------------------------------------
  subroutine mask2outmap(outmap, map2heal, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    integer(i4b),   dimension(1:),    intent(in)  :: map2heal
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)                                  :: n, i

    n = size(map2heal)

    !Writing mask to outmap
    outmap = 0.d0
    if (pol == 1) then
       do i = 1, n
          outmap(map2heal(i),1) = 1.d0
       end do
    else if (pol == 2) then
       do i = 1, n
          outmap(map2heal(i),2) = 1.d0
          outmap(map2heal(i),3) = 1.d0
       end do
    end if
       
  end subroutine mask2outmap

  !----------------------------------------------------------------
  ! writes a mask to a healpix outmap using map2mask
  !----------------------------------------------------------------
  subroutine mask2outmap_viamap2mask(outmap, map2mask, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    integer(i4b),   dimension(0:),    intent(in)  :: map2mask
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)                                  :: npix, i

    npix = size(map2mask)

    !Writing mask to outmap
    outmap = 0.d0
    if (pol == 1) then
       do i = 0, npix-1
          if (map2mask(i) > 0) outmap(i,1) = 1.d0
       end do
    else if (pol == 2) then
       do i = 0, npix-1
          if (map2mask(i) > 0) then
             outmap(i,2) = 1.d0
             outmap(i,3) = 1.d0
          end if
       end do
    end if
       
  end subroutine mask2outmap_viamap2mask

  !---------------------------------------------------------------------------------
  ! Writing rms map to outmap given diag of covariance matrix and map2heal
  !-------------------------------------------------------------------------------

  subroutine rms2outmap(outmap, covdiag, map2heal, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    real(dp),       dimension(1:),    intent(in)  :: covdiag
    integer(i4b),   dimension(1:),    intent(in)  :: map2heal
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)                                  :: n, i

    n = size(map2heal)

    !Writing rms map to outmap
    outmap = -1.6375d30
    if (pol == 1) then
       do i = 1, n
          outmap(map2heal(i),1) = sqrt(covdiag(i))
       end do
    else if (pol == 2) then
       do i = 1, n
          outmap(map2heal(i),2) = sqrt(covdiag(i))
          outmap(map2heal(i),3) = sqrt(covdiag(i+n))
       end do
    end if
       
  end subroutine rms2outmap

  !---------------------------------------------------------------------------------
  ! Writing rms map to outmap given diag of covariance matrix and map2mask
  !-------------------------------------------------------------------------------

  subroutine rms2outmap_viamap2mask(outmap, covdiag, map2mask, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    real(dp),       dimension(1:),    intent(in)  :: covdiag
    integer(i4b),   dimension(0:),    intent(in)  :: map2mask
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)                                  :: npix, i, n

    n = size(covdiag)
    npix = size(map2mask)

    !Writing rms map to outmap
    outmap = -1.6375d30
    if (pol == 1) then
       do i = 0, npix-1
          if (map2mask(i) > 0) outmap(i,1) = sqrt(covdiag(map2mask(i)))
       end do
    else if (pol == 2) then
       do i = 0, npix-1
          if (map2mask(i) > 0) then
             outmap(i,2) = sqrt(covdiag(map2mask(i)))
             outmap(i,3) = sqrt(covdiag(map2mask(i)+n/2))
          end if
       end do
    end if
       
  end subroutine rms2outmap_viamap2mask

  !---------------------------------------------------------------------------------
  ! Make mask given diag of covariance matrix and map2heal
  !-------------------------------------------------------------------------------

  subroutine rms2mask(outmap, covdiag, map2heal, pol)
    implicit none

    real(dp),       dimension(0:,1:), intent(out) :: outmap
    real(dp),       dimension(1:),    intent(in)  :: covdiag
    integer(i4b),   dimension(1:),    intent(in)  :: map2heal
    integer(i4b),                     intent(in)  :: pol
  
    integer(i4b)            :: n, i
    real(dp)                :: minrms, lim

    lim = 2.d0**2

    n = size(map2heal)
    minrms = minval(covdiag)

    !Writing mask to outmap
    outmap = 0.d0
    if (pol == 1) then
       do i = 1, n
          if (covdiag(i) < lim*minrms) outmap(map2heal(i),1) = 1.d0
       end do
    else if (pol == 2) then
       do i = 1, n
          if (covdiag(i) < lim*minrms) outmap(map2heal(i),2) = 1.d0
          if (covdiag(i) < lim*minrms) outmap(map2heal(i),3) = 1.d0
       end do
    end if
       
  end subroutine rms2mask

  !---------------------------------------------------------------------------
  ! Make map2heal from map2mask
  !---------------------------------------------------------------------------
  subroutine map2mask2map2heal(map2mask, map2heal, n)
    implicit none

    integer(i4b),          dimension(0:), intent(in)  :: map2mask
    integer(i4b), pointer, dimension(:)               :: map2heal
    integer(i4b), optional,               intent(in)  :: n        ! not used!

    integer(i4b)                                      :: i, j, npix, m

    npix = size(map2mask)
    m = count(map2mask>0)
    allocate(map2heal(m))
    j=0
    do i = 0, npix-1
       if (map2mask(i) > 0) then
          j = j+1
          map2heal(j)= i
       end if
    end do
    call assert(j == m, 'map2mask2map2heal: something is wrong')
  end subroutine map2mask2map2heal
  !---------------------------------------------------------------------------
  ! Make map2heal from mask
  !---------------------------------------------------------------------------
  subroutine mask2heal(mask, map2heal)
    implicit none

    integer(i4b),              dimension(0:), intent(in)  :: mask
    integer(i4b), allocatable, dimension(:)               :: map2heal

    integer(i4b)                                          :: i, j, npix, n

    npix = size(mask)
    n = count(mask == 1)
    write(*,*) n,'=n', npix,'=npix'
    allocate(map2heal(n))
    j=0
    do i = 0, npix-1
       if (mask(i) > 0) then
          j = j+1
          map2heal(j)= i
       end if
    end do
    call assert(j == n, 'mask2heal: something is wrong')
  end subroutine mask2heal

  !---------------------------------------------------------------------------
  ! Removing (unhit), large, uncommon (and masked out) pixels from map2mask
  !---------------------------------------------------------------------------
  subroutine mask2common(myid, map2mask1, map2mask2, in2red1, in2red2, map2heal, n1, n2, n, npix, diag1, diag2, mask)
    implicit none

    integer(i4b),                          intent(in)    :: n1, n2, npix, myid
    integer(i4b),                          intent(out)   :: n
    integer(i4b),           dimension(0:), intent(inout) :: map2mask1, map2mask2
    integer(i4b), pointer,  dimension(:)                 :: in2red1, in2red2, map2heal
    real(dp),     optional, dimension(:),  intent(in)    :: diag1, diag2
    integer(i4b), optional, dimension(0:), intent(in)    :: mask

    integer(i4b)                             :: i, j, dummy, num, ierr
    integer(i4b), allocatable, dimension(:)  :: nhits, finalmask1, finalmask2, usemask
    
    ! Prepare mask
    allocate(usemask(0:npix-1))
    if (present(mask)) then
       usemask = mask
    else
       usemask = 1
    end if
 
    ! Remove unhit pixels
    if (present(diag1) .and. present(diag2)) then
       if (myid==0) then
          write(*,*) count(diag1 == 0.d0), '= number of unhit pixels in set 1'
          write(*,*) count(diag2 == 0.d0), '= number of unhit pixels in set 2'
       end if
       do i = 0, npix-1
          if (map2mask1(i) > 0) then
             if (diag1(map2mask1(i)) == 0.d0) map2mask1(i) = -1
          end if
          if (map2mask2(i) > 0) then
             if (diag2(map2mask2(i)) == 0.d0) map2mask2(i) = -1
          end if
       end do
    end if

    ! Removing big pixels from final map2mask
    allocate(finalmask1(0:npix-1))             
    allocate(finalmask2(0:npix-1)) 
    finalmask1 = -1
    finalmask2 = -1

    allocate(nhits(n1))
    nhits = 0
    do i = 0, npix-1
       dummy = map2mask1(i)
       if (dummy /= -1) nhits(dummy) = nhits(dummy) + 1
    end do
    do  i = 0, npix-1
       dummy = map2mask1(i)
       if (dummy /= -1) then
          if (nhits(dummy) == 1) finalmask1(i) = map2mask1(i)
       end if
    end do
    deallocate(nhits)

    allocate(nhits(n2))
    nhits = 0
    do i = 0, npix-1
       dummy = map2mask2(i)
       if (dummy /= -1) nhits(dummy) = nhits(dummy) + 1
    end do
    do  i = 0, npix-1
       dummy = map2mask2(i)
       if (dummy /= -1) then
          if (nhits(dummy) == 1) finalmask2(i) = map2mask2(i)
       end if
    end do
    deallocate(nhits)

    do i = 0, npix-1
       if ( finalmask1(i) == -1 .or. finalmask2(i) == -1 .or. usemask(i) /= 1 ) then
          finalmask1(i) = -1
          finalmask2(i) = -1
       end if
    end do

    ! Checking for pixels left
    n = count(finalmask1 /= -1)                 ! ='n/2' for polarisation
    if (myid==0) write(*,*) n, '= number of pixels in final cov matrix'
    if (n < 1) then
       write(*,*) 'map_1 and map_2 (and mask) have no common pixels. Quiting' 
       call mpi_finalize(ierr)       
       stop
    end if

    ! from common pixels to in2red and map2heal
    allocate(in2red1(n)) 
    allocate(in2red2(n))  
    allocate(map2heal(n))
    j=0
    do i = 0, npix-1
       if ( finalmask1(i) /= -1 ) then
          j = j+1       
          in2red1(j)   = finalmask1(i)
          map2heal(j) = i                        ! faktiske healpixlene some er med  
       end if
    end do
    j=0
    do i = 0, npix-1
       if ( finalmask2(i) /= -1 ) then
          j = j+1       
          in2red2(j)   = finalmask2(i)
       end if
    end do
    
    deallocate(finalmask1)
    deallocate(finalmask2)

  end subroutine mask2common

  !---------------------------------------------------------------------------
  ! Removing (unhit), large (and masked out) pixels from map2mask
  !---------------------------------------------------------------------------
  subroutine map2mask2in2red(myid, map2mask, in2red, map2heal, n, npix, ordering, mask, diag)
    implicit none

    integer(i4b),                           intent(in)    :: npix, myid, ordering
    integer(i4b),                           intent(out)   :: n
    integer(i4b),            dimension(0:), intent(inout) :: map2mask
    integer(i4b), optional,  dimension(0:), intent(in)    :: mask
    integer(i4b), pointer,   dimension(:)                 :: in2red, map2heal
    real(dp),     optional,  dimension(:),  intent(in)    :: diag
  
    integer(i4b)                             :: i, j, dummy, num, ierr, nside, listpix(10000), m
    integer(i4b), allocatable, dimension(:)  :: donald, finalmask, nhits
    real(dp)                                 :: k, vec0(3), r

    ! Remove unhit pixels
    if (present(diag)) then
       if (myid==0) write(*,*) count(diag == 0.d0), '= number of unhit pixels'
       do i = 0, npix-1
          if (map2mask(i) > 0) then
             if (diag(map2mask(i)) == 0.d0)  map2mask(i) = -1
          end if
       end do
    end if

    num = count(map2mask > 0)
    ! Removing big pixels from final map2mask
    allocate(finalmask(0:npix-1))
    finalmask = -1
    allocate(nhits(num))
    nhits = 0
    do i = 0, npix-1
       dummy = map2mask(i)
       if (dummy /= -1) nhits(dummy) = nhits(dummy) + 1
    end do
    k=0
    do  i = 0, npix-1
       dummy = map2mask(i)
       if (dummy > 0) then
          if (nhits(dummy) == 1) then 
             finalmask(i) = map2mask(i)
          else
             k = k + 1/real(nhits(dummy), dp)
          end if
       end if
    end do
    deallocate(nhits)
    if (myid==0) write(*,*) int(k, i4b), '= number of large pixels removed'

    ! Removing single isolated pixels -- if there are less than five pixels within
    ! a radius of 4 pixels, remove it
    nside = nint(sqrt(real(npix,sp)/12.))
    r = 4.d0 * 60.d0 / real(nside,dp) * DEG2RAD ! Include a radius of ~4 pixels
    do i = 0, npix-1
       if (finalmask(i) /= -1) then
          if (ordering == 1) then
             call pix2vec_ring(nside, i, vec0)
          else
             call pix2vec_nest(nside, i, vec0)
          end if
          call query_disc(nside, vec0, r, listpix, m, nest=ordering-1)
          if (count(map2mask(listpix(:m)) > 0) < 5) finalmask(i) = -1
       end if
    end do

    if (present(mask)) then
       ! Finding num pixels common for final map2mask and mask
       n = count(finalmask > 0 .and. mask == 1)
       if (n==0) then
          if (myid==0) write(*,*) 'Map2mask and mask have no common pixels in current component. Returning'
          return
       end if
    else
       n = count(finalmask /= -1)
    end if
    if (myid==0) write(*,*) n, '= number of pixels in final cov matrix'

    ! Removing masked out pixels
    allocate(in2red(n))
    allocate(map2heal(n))
    map2mask = -1

    if (present(mask)) then
       j=0
       do i = 0, npix-1
          if (finalmask(i) /= -1 .and. mask(i) == 1) then
             j = j+1       
             in2red(j)   = finalmask(i)
             map2heal(j) = i                        ! faktiske healpixlene some er med  
             map2mask(i) = j
          end if
       end do
    else
       j=0
       do i = 0, npix-1
          if (finalmask(i) /= -1) then
             j = j+1       
             in2red(j)   = finalmask(i)
             map2heal(j) = i                        ! faktiske healpixlene some er med  
             map2mask(i) = j
          end if
       end do
    end if

    deallocate(finalmask)

  end subroutine

  !---------------------------------------------------------------------------
  ! Removing unhit pixels from map2mask
  !---------------------------------------------------------------------------
  subroutine unhit_remove(myid, map2mask, diag, in2red, pol)
    implicit none

    integer(i4b),                           intent(in)    :: myid, pol
    integer(i4b),            dimension(0:), intent(inout) :: map2mask
    integer(i4b), pointer,   dimension(:)                 :: in2red
    real(dp),                dimension(:),  intent(in)    :: diag
  
    integer(i4b)                             :: i, j, dummy, num, ierr, n, npix, k

    num = count(diag == 0.d0) 
    if (myid==0) write(*,*) num, '= number of unhit pixels removed'
    n = (size(diag) - num)
    if (myid==0) write(*,*) n, '= number of pixels left' 
    !if (myid==0) write(*,*) count(map2mask /= -1), '= number of pixels in map2mask to start with'

    npix = size(map2mask)
    ! Remove unhit pixels
    k = 0
    do i = 0, npix-1
       if (map2mask(i) > 0) then
          if (diag(map2mask(i)) == 0.d0) then
!             write(*,*) myid, i, map2mask(i), diag(map2mask(i)) 
             k = k +1
             map2mask(i) = -1
          end if
       end if
    end do
    !if (myid==0) write(*,*) k, '= number of unhit pixels removed from m2m'

    ! Making new in2red and map2mask
    allocate(in2red(n))
    j=0
    do i = 1, size(diag)
       if (diag(i) /= 0.d0) then
          j = j+1       
          in2red(j)   = i
       end if
    end do
    do i = 0, npix-1
       if (map2mask(i) > 0) then
          do j = 1, n
             if (in2red(j) == map2mask(i)) map2mask(i)=j 
          end do
       end if
    end do
    
    !if (myid==0) write(*,*) count(map2mask > 0), '= number of pixels in map2mask'

  end subroutine

  !---------------------------------------------------------------------------
  ! Removing masked out pixels from map2mask
  !---------------------------------------------------------------------------
  subroutine mask2in2red(myid, mask, map2mask, in2red, n, npix)
    implicit none

    integer(i4b),                         intent(in)  :: npix, myid
    integer(i4b),                         intent(out) :: n
    integer(i4b),          dimension(0:), intent(in)  :: mask, map2mask
    integer(i4b), pointer, dimension(:)               :: in2red

    integer(i4b)                             :: i, j, ierr

    n = count(map2mask > 0 .and. mask == 1)
    if(n == 0) then
       write(*,*) 'Map2mask and mask have no common pixels. Quiting'
       call mpi_finalize(ierr)
       stop
    end if

    allocate(in2red(n))
    j=0
    do i = 0, npix-1
       if (map2mask(i) > 0 .and. mask(i) == 1) then
          j = j+1       
          in2red(j) = map2mask(i)
       end if
    end do
  end subroutine mask2in2red

  !---------------------------------------------------------------------
  ! Removing large and uncommon pixels from map and cov using in2red
  !----------------------------------------------------------------------

  subroutine inmap2red(in2red, inmap, incov, redmap, redcov, pol, inmat, redmat)
   implicit none

    integer(i4b),       dimension(1:),    intent(in)  :: in2red
    real(dp),           dimension(1:),    intent(in)  :: inmap
    real(dp),           dimension(1:,1:), intent(in)  :: incov
    real(dp),  pointer, dimension(:),     intent(out) :: redmap
    real(dp),  pointer, dimension(:,:),   intent(out) :: redcov
    integer(i4b),                         intent(in)  :: pol  
    real(dp),  optional,          dimension(1:,1:), intent(in)  :: inmat
    real(dp),  optional, pointer, dimension(:,:),   intent(out) :: redmat

    integer(i4b)                             :: n
    integer(i4b), allocatable, dimension(:)  :: polvector

    n = size(in2red)

    ! Removing bad pixels from matrices
     if (pol == 1) then
       allocate(redmap(n))      
       allocate(redcov(n, n))
       redmap = inmap(in2red)       
       redcov = incov(in2red, in2red)
       if (present(redmat)) then 
          allocate(redmat(n, n))
          redmat = inmat(in2red, in2red)
       end if
    else if (pol == 2) then
       allocate(redmap(2*n))    
       allocate(redcov(2*n, 2*n))
       allocate(polvector(2*n))
       polvector(1:n)     = in2red
       polvector(n+1:2*n) = in2red + size(inmap)/2
       redmap = inmap(polvector)
       redcov = incov(polvector,polvector)
       if (present(redmat)) then 
          allocate(redmat(2*n, 2*n))
          redmat = inmat(polvector,polvector)
       end if
       deallocate(polvector)
    end if

  end subroutine inmap2red

  !---------------------------------------------------------------------
  ! Read and check mask
  !----------------------------------------------------------------------

  subroutine read_and_check_mask(myid, maskfile, mask, nmaps, ordering, npix, pol)
   implicit none

    character(len=*),                     intent(in) :: maskfile
    integer(i4b),                         intent(in) :: nmaps, ordering, npix, pol, myid
    integer(i4b), pointer, dimension(:)              :: mask

    integer(i4b)       :: npix_mask, nmaps_mask, ordering_mask, nside_mask, ierr
    logical(lgt)       :: anynull
    real(dp)           :: nullval
    real(dp),     allocatable, dimension(:,:) :: inmask
 
    ! Read mask
    if (myid==0) write(*,*) 'Reading from ', trim(maskfile)
    call read_map(inmask, ordering_mask, maskfile, nside=nside_mask, nmap=nmaps_mask)
    if (myid==0) write(*,*) nside_mask, '= nside_mask,', nmaps_mask, '= nmaps_mask,',ordering_mask,'= ordering_mask'
    npix_mask=12*nside_mask**2
    ! Check that both input maps is of same pixel size
    if (npix_mask /= npix) then
       if (myid==0) then 
          write(*,*) 'Input map(s) and mask is not of same resolution. Quiting'
       end if
       stop
    end if
    ! Check ordering
    if (ordering_mask /= 1 .and. ordering_mask /= 2) then
       if (myid==0) write(*,*) "Ordering is neither ring nor nested. Quiting"
       stop
    end if
    ! Convert ordering if necessary
    if (ordering_mask  /= ordering) then
       if (ordering == 1) then 
          call convert_nest2ring(nside_mask ,inmask)
          if (myid==0) write(*,*) 'Converting mask from nested to ring'
       else if (ordering == 2) then 
          call convert_ring2nest(nside_mask ,inmask)
          if (myid==0) write(*,*) 'Converting mask from ring to nested'
       end if
    end if
    ! Check nmaps/polarisation
    if (nmaps /= nmaps_mask) then
       if (nmaps_mask == 1 .and. pol == 2) then
          if (myid==0) then 
             write(*,*) 'OBS: Using component 1 mask for polarisation data'
          end if
       else if (nmaps_mask == 3 .and. pol == 1) then
          if (myid==0) then 
             write(*,*) 'OBS: Mask has 3 components. Using component 1 for temperature data'
          end if
          nmaps_mask = 1
       else
          if (myid==0) then 
             write(*,*) 'Confused. nmaps_mask =',nmaps_mask,'polarisation =',pol,' Quiting'
          end if
          stop
       end if
    end if
    ! Put mask in 1-dim mask
    allocate(mask(0:npix-1))
    mask = nint(inmask(:,nmaps_mask))
    deallocate(inmask)
  end subroutine read_and_check_mask

  ! Thesese subroutines assume 3 fields in the map.
  subroutine map2array(map, map2mask, nfield, arr)
    implicit none
    real(dp),     dimension(0:,:) :: map
    integer(i4b), dimension(0:)   :: map2mask
    real(dp),     dimension(:), allocatable :: arr
    integer(i4b), dimension(:), pointer     :: map2heal
    integer(i4b) :: nfield, i, n, off
    n = maxval(map2mask)
    off = 0; if(nfield == 2) off = 1
    if(.not. allocated(arr)) allocate(arr(n*nfield))
    call map2mask2map2heal(map2mask, map2heal, n)
    do i = 1, nfield
       arr((i-1)*n+1:i*n) = map(map2heal, i+off)
    end do
    deallocate(map2heal)
  end subroutine

  subroutine map2array_prealloc(map, map2mask, nfield, arr)
    implicit none
    real(dp),     dimension(0:,:) :: map
    integer(i4b), dimension(0:)   :: map2mask
    real(dp),     dimension(:)    :: arr
    integer(i4b), dimension(:), pointer     :: map2heal
    integer(i4b) :: nfield, i, n, off
    n = maxval(map2mask)
    off = 0; if(nfield == 2) off = 1
    call map2mask2map2heal(map2mask, map2heal, n)
    do i = 1, nfield
       arr((i-1)*n+1:i*n) = map(map2heal, i+off)
    end do
    deallocate(map2heal)
  end subroutine

  subroutine array2map(arr, map2mask, nfield, map)
    implicit none
    real(dp),     dimension(:)    :: arr
    integer(i4b), dimension(0:)   :: map2mask
    real(dp),     dimension(:,:), allocatable :: map
    integer(i4b), dimension(:),   pointer     :: map2heal
    integer(i4b) :: nfield, i, n, off
    n = maxval(map2mask)
    off = 0; if(nfield == 2) off = 1
    if(.not. allocated(map)) allocate(map(0:size(map2mask), 3))
    call map2mask2map2heal(map2mask, map2heal, n)
    map = 0
    do i = 1, nfield
       map(map2heal, i+off) = arr((i-1)*n+1:i*n)
    end do
    deallocate(map2heal)
  end subroutine

  subroutine array2map_prealloc(arr, map2mask, nfield, map)
    implicit none
    real(dp),     dimension(:)    :: arr
    integer(i4b), dimension(0:)   :: map2mask
    real(dp),     dimension(:,:)  :: map
    integer(i4b), dimension(:),   pointer     :: map2heal
    integer(i4b) :: nfield, i, n, off
    n = maxval(map2mask)
    off = 0; if(nfield == 2) off = 1
    call map2mask2map2heal(map2mask, map2heal, n)
    map = 0
    do i = 1, nfield
       map(map2heal, i+off) = arr((i-1)*n+1:i*n)
    end do
    deallocate(map2heal)
  end subroutine

end module quiet_postutils
