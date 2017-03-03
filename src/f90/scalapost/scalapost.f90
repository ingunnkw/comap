program scalapost
  use scalawrap
  use quiet_fileutils
  use quiet_postutils
  use quiet_utils
  use scalautils
  use pix_tools
  use quiet_mapfile_mod
  use quiet_healpix_mod
  use quiet_angcorr_mod
  use alm_tools
  use rngmod
  implicit none

  integer(i4b)       :: unit, root, ierr
  type(scinfo)       :: info
  character(len=30)  :: kommando

  call sc_init(info, active=.true.)
  call dset(id=info%myid)
  unit=42+info%myid
  root = 0

  if (info%myid == root) then
     write(*,*) '-------------------- QUIET map post-processor  --------------------'
     write(*,*)
  end if

  ! Get name of main kommando
  call getarg(1,kommando)

  if (kommando == 'finalmap') then
     if (info%myid == root) write(*,*) 'Make final map and diff_map and covariance matrix from two subsets'
     call make_finalmap_sc(unit, info)
  else if (kommando == 'genmap') then
     if (info%myid == root) write(*,*) 'Make final map and diff_map and covariance matrix from two subsets'
     call make_genmap_sc(unit, info)
  else if (kommando == 'solve') then
     if (info%myid == root) write(*,*) 'Solve an equation system in the format produced by tod2map'
     call solve_equation_sc(unit, info)
  else if (kommando == 'lcut') then
     if (info%myid == root) write(*,*) 'Project out modes lower or equal to l_max'
     call project_out_lowmodes_sc(unit, info)
  else if (kommando == 'cutoff') then
     if (info%myid == root) write(*,*) 'Project out eigenmodes lower than cutoff'
     call project_out_eigenmodes_sc(unit, info)
  else if (kommando == 'sqrt') then
     if (info%myid == root) write(*,*) 'Find sqrt of (inv)cov'
     call get_sqrt_sc(unit, info)
  else if (kommando == 'invertcov') then
     if (info%myid == root) write(*,*) 'Find inverse of (inv)cov'
     call get_inverse_sc(unit, info)
  else if (kommando == 'prune') then
     if (info%myid == root) write(*,*) 'Prune multiresolution pixels. Did not I delete this function?'
     call prune_multires(unit, info)
  else if (kommando == 'convert_C_to_F90') then
     if (info%myid == root) write(*,*) 'Converting matrix from C to F90 format'
     call convert_C_to_f90(unit, info)
  else if (kommando == 'convert_WMAP_to_F90') then
     if (info%myid == root) write(*,*) 'Converting matrix from WMAP to F90 format'
     call convert_WMAP_to_f90(unit, info)
  else if (kommando == 'invert') then
     if (info%myid == root) write(*,*) "Simply invert a matrix"
     call spost_invert_matrix(unit, info)
  else if (kommando == 'mask') then
     if (info%myid == root) write(*,*) "Mask matrix"
     call spost_mask_matrix(unit, info)
  else if (kommando == 'fit_template') then
     if (info%myid == root) write(*,*) "Fit template to map"
     call fit_template(unit, info)
  else if (kommando == 'multiply_by_scalar') then
     if (info%myid == root) write(*,*) 'Multiply matrix by scalar'
     call operate_on_matrix_with_scalar_sc(unit, info, '*')
  else if (kommando == 'add_scalar') then
     if (info%myid == root) write(*,*) 'Add scalar to matrix'
     call operate_on_matrix_with_scalar_sc(unit, info, '+')
  else if (kommando == 'add') then
     if (info%myid == root) write(*,*) "Add two matrices"
     call operate_on_matrices_sc(unit, info, '+')
  else if (kommando == 'fit_template_with_noise ') then
     if (info%myid == root) write(*,*) "Fit noisy template to map"
     call fit_template_with_noise(unit, info)
  else if (kommando == 'fit_gaussian') then
     if (info%myid == root) write(*,*) "Fit 2D Gaussian to a map feature"
     call fit_gaussian(unit, info)
  else if (kommando == 'covcat') then
     if (info%myid == 0) write(*,*) "Block-diagonal concatenation of covariance matrices"
     call concat_cov(unit, info)
  else if (kommando == 'eigendecomp') then
     if (info%myid == 0) write(*,*) "Compute eigenvalue decomposition of matrix"
     call eigendecomp_cmd(info)
  else if (kommando == 'eigencut' .or. kommando == 'eigenswap') then
     call eigencut_cmd(info)
  else if (kommando == 'eigenchi') then
     call eigenchi_cmd(info)
  else if (kommando == 'convert_mask_to_map2mask') then
     call convert_mask_to_map2mask(info)
  else if (kommando == 'weighted_sum') then
     call compute_weighted_sum(info)
  else if (kommando == 'upadd') then
     call upadd(info)
  else if (kommando == 'buildmats') then
     call buildmats(info)
  else if (kommando == 'apply_mask') then
      call apply_mask(info)
  else if (kommando == 'cov2rms') then
      call cov2rms(unit, info)
  else if (kommando == 'rms2cov') then
      call rms2cov(unit, info)
  else if (kommando == 'build_S_mat') then
     call build_S_matrix(info)
  else if (kommando == 'extract_minimaps') then
     call extract_minimaps(info)
  else if (kommando == 'co_add') then
     if (info%myid == 0) write(*,*) "Co-add maps with invN weighing"
     call co_add_maps(unit, info)
  else if (kommando == 'cov2corr') then
     call cov2corr(unit, info)
  else if (kommando == 'test') then
     call test(unit, info)
  else
     call give_user_info_sc(info)
  end if

  ! And exit
  call mpi_finalize(ierr)
  if (info%myid == root) then 
     write(*,*)
     write(*,*) '-------------------- QUIET map post-processor completed  ----------'
  end if

contains

!-----------------------------------------------------------------------------------------------
! subroutine make_finalmap
!-----------------------------------------------------------------------------------------------

  subroutine make_finalmap_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)              :: amat1, amat2, bmat1, bmat2, rhs1, rhs2
    type(scalamat)              :: amat1in, amat2in, bmat1in, bmat2in, rhs1in, rhs2in
    type(scalamat)              :: cov1, cov2, map1, map2
    type(scalamat)              :: redcov1, redcov2, redmap1, redmap2
    type(scalamat)              :: mapsum, mapdiff
    character(len=256)          :: sub1file, sub2file, outprefix, outfile, outtext, maskfile, simprefix
    integer(i4b)                :: ierr, root=0, ignore, ig1, ig2
    integer(i4b)                :: ordering, ordering2, pol, pol2, nrhs
    integer(i4b)                :: n, n1, n2, npix, npix2, i, j, int, nmaps
    logical(lgt)                :: inv, fast, havemask, asym1, asym2
    real(dp)                    :: t1, t2

    integer(i4b), pointer,     dimension(:)   :: in2red1, in2red2, map2heal1, map2heal2, healvec, mask
    integer(i4b), pointer,     dimension(:)   :: map2mask1, map2mask2
    real(dp),     allocatable, dimension(:,:) :: outmap, locmap1, locmap2 
    real(dp),     allocatable, dimension(:)   :: diag1, diag2

    fast = .false.
    havemask = .false.

    ! Get parameters
    if (iargc() /= 4 .and. iargc() /= 5) then
       write(*,*) 'finalmap takes 3(4) parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, sub1file)
       call getarg(3, sub2file)
       call getarg(4, outprefix)
    end if

    if (iargc() == 5) then 
       call getarg(5, maskfile)
       if (trim(maskfile)=='fast') then
          fast = .true.
          if (info%myid==0)  write(*,*) 'Running in fast mode; might not work for bad matrices'
       else
          havemask = .true.
          if (info%myid==0)  write(*,*) 'Will apply mask ', trim(maskfile),' to data'
       end if
    end if

    ! Read equation set 1
    call read_eqs_combo_sc(unit, info, sub1file, n1, ordering, pol, amat1in, bmat1in, rhs1in, npix, map2mask1)
    call assert(count(map2mask1 > 0) >= n1/pol, 'size cov1 and map2mask1 incosistent')
    ! Read equation set 2
    call read_eqs_combo_sc(unit, info, sub2file, n2, ordering2, pol2, amat2in, bmat2in, rhs2in, npix2, map2mask2)
    call assert(count(map2mask2 > 0) >= n2/pol, 'size cov2 and map2mask2 incosistent')

!!$    call check_input(info%myid, ordering, ordering2, pol, pol2, npix, npix2, temperature)
!!$    call check_sym_matrix(invcov1in%data, asym1)
!!$    call check_sym_matrix(invcov2in%data, asym2)
!!$    if (info%myid==root) then
!!$       if (asym1) write(*,*) 'Covariance matrix 1 is not symmetric.'
!!$       if (asym2) write(*,*) 'Covariance matrix 2 is not symmetric.'
!!$       if (asym1 .or. asym2) write(*,*) 'Covariance matrices are not symmetric. Quitting'
!!$     end if
!!$    if (asym1 .or. asym2) then
!!$       call mpi_finalize(ierr)
!!$       stop
!!$    end if
!!$    if (info%myid==root) then
!!$       write(*,*) 'Number of rhs in equation set 1 =',rhs1in%cglob
!!$       write(*,*) 'Number of rhs in equation set 2 =',rhs2in%cglob
!!$    end if
    if (rhs1in%cglob /= rhs2in%cglob) then 
       if (info%myid==root) write(*,*) 'Not an equal number of right-hand-sides. Quiting'
       call mpi_finalize(ierr)
       stop
    end if
    nrhs = rhs1in%cglob
    if (nrhs > 1) call mkdir(trim(outprefix) // '_sims')
    simprefix = trim(outprefix) // '_sims/' // trim(outprefix )

    if (pol==1) then
       if (info%myid==0) write(*,*) 'Running on temperature data'
       nmaps = 1
    else if (pol==2) then
       if (info%myid==0) write(*,*) 'Running on polarisation data'
       nmaps = 3
    else 
       write(*,*) 'Unknown polarisation. Quiting'
       stop
    end if
    ! Read mask
    if (havemask) call read_and_check_mask(info%myid, maskfile, mask, nmaps, ordering, npix, pol)

    ! Removing unhit pixels
    allocate(diag1(n1))
    allocate(diag2(n2))
    call sc_get_diag(amat1in, diag1)
    call sc_get_diag(amat2in, diag2)
    call unhit_remove(info%myid, map2mask1, diag1, in2red1, pol)
    call unhit_remove(info%myid, map2mask2, diag2, in2red2, pol)
    call inmap2red_sc(info, in2red1, rhs1in, amat1in, rhs1, amat1, 1, bmat1in, bmat1)
    call inmap2red_sc(info, in2red2, rhs2in, amat2in, rhs2, amat2, 1, bmat2in, bmat2)
   
    deallocate(diag1)
    deallocate(diag2)
    n1 = size(in2red1)
    n2 = size(in2red2)
    deallocate(in2red1)
    deallocate(in2red2)

    ! Writing input rhs to file
    allocate(locmap1(n1,nrhs))
    call sc_get(rhs1, locmap1, 0)
    !call sc_get(rhs1, locmap1, 0, (/1,1, n1,nrhs/))   
    if (info%myid==0) then
       call write_maps_m2m(unit, locmap1, map2mask1, pol, npix, ordering, outprefix, simprefix, 'input_rhs1.fits')
    end if
    deallocate(locmap1)
    allocate(locmap1(n2,nrhs))
    call sc_get(rhs2, locmap1, 0)
    if (info%myid==0) then
       call write_maps_m2m(unit, locmap1, map2mask2, pol, npix, ordering, outprefix, simprefix, 'input_rhs2.fits')
    end if
    deallocate(locmap1)

    ! Solve equation set 1
    outtext = 'input_map1'
    if (fast) then
       call fast_solve_eqs_sc(n1, amat1, rhs1, map1, cov1, bmat1)
       ig1 = 0
    else 
       call solve_eqs_sc(info, unit, n1, amat1, rhs1, map1, cov1, bmat=bmat1, outtext=outtext, kept=ig1)
       ig1 = n1 - ig1
    end if
    allocate(locmap1(n1, nrhs))
    call sc_get(map1, locmap1, 0)   
    ! Writing input map1 to file
    if (info%myid==0) then
       call write_maps_m2m(unit, locmap1, map2mask1, pol, npix, ordering, outprefix, simprefix, 'input_map1.fits')
    end if

    ! Solve equation set 2
    outtext = 'input_map2'
    if (fast) then
       call fast_solve_eqs_sc(n2, amat2, rhs2, map2, cov2, bmat2)
       ig2 = 0
    else
       call solve_eqs_sc(info, unit, n2, amat2, rhs2, map2, cov2, bmat=bmat2, outtext=outtext, kept=ig2)
       ig2 = n2 - ig2
    end if
    ignore = min(ig1,ig2)
    allocate(locmap2(n2,nrhs))
    call sc_get(map2, locmap2, 0)   
    ! Writing input map2 to file
    if (info%myid==0) then
       call write_maps_m2m(unit, locmap2, map2mask2, pol, npix, ordering, outprefix, simprefix, 'input_map2.fits')
    end if

    ! Removing large and uncommon pixels from map2mask
    if (havemask) then
       call mask2common(info%myid, map2mask1, map2mask2, in2red1, in2red2, healvec, n1, n2, n, npix, mask=mask)
    else
       call mask2common(info%myid, map2mask1, map2mask2, in2red1, in2red2, healvec, n1, n2, n, npix)
    end if
    if (pol == 2) n = 2*n
    if (info%myid==0) write(*,*) "Number of modes is",n
    deallocate(map2mask1)
    deallocate(map2mask2)
    !.. and from maps and cov matrices
    call inmap2red_sc(info, in2red1, map1, cov1, redmap1, redcov1, pol)
    deallocate(in2red1)
    call inmap2red_sc(info, in2red2, map2, cov2, redmap2, redcov2, pol)
    deallocate(in2red2)
    deallocate(locmap1)
    deallocate(locmap2)
    allocate(locmap1(n,nrhs))
    allocate(locmap2(n,nrhs))
    call sc_get(redmap1, locmap1, 0)   
    call sc_get(redmap2, locmap2, 0)  

    if (info%myid==0) then
       ! Writing redused map1 to file
       call write_maps(unit, locmap1, healvec, pol, npix, ordering, outprefix, simprefix, 'red_map1.fits')
    
       ! Writing redused map2 to file
       call write_maps(unit, locmap2, healvec, pol, npix, ordering, outprefix, simprefix, 'red_map2.fits')
   
       ! Writing true diffmap to file
       call write_maps(unit, (locmap1-locmap2)/2, healvec, pol, npix, ordering, outprefix, simprefix, 'true_diffmap.fits')
    
       ! Writing direct sum map to file
       call write_maps(unit, (locmap1+locmap2)/2, healvec, pol, npix, ordering, outprefix, simprefix, 'direct_summap.fits')
    
       ! Note: map2mask is not consistent with an empty temperature layer. TMR
              
       ! Writing map2mask to file
       allocate(outmap(0:npix-1, nmaps))
       call map2mask2outmap(outmap, healvec, pol)   
       outfile = trim(outprefix) // '_map2mask.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Map2mask written to file = ', trim(outfile)
    
       ! Writing mask to file
       call mask2outmap(outmap, healvec, pol)   
       outfile = trim(outprefix) // '_mask.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Mask written to file = ', trim(outfile)

       deallocate(outmap, locmap1, locmap2)
    end if

    ! Produce and output final maps and covs
    call produce_finalmaps_sc(unit, info, redmap1, redmap2, redcov1, redcov2, healvec, outprefix, simprefix, ordering, pol, nmaps, npix, nrhs, fast, ignore)

  end subroutine make_finalmap_sc

!-----------------------------------------------------------------------------------------------
! subroutine make_genmap
!-----------------------------------------------------------------------------------------------

  subroutine make_genmap_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)              :: cov1, cov2, map1, map2
    type(scalamat)              :: mapsum, mapdiff
    character(len=256)          :: map1file, map2file, cov1file, cov2file, maskfile
    character(len=256)          :: outprefix, runmode, outfile, outtext, simprefix
    integer(i4b)                :: ierr, root=0
    integer(i4b)                :: ordering, ordering2, pol, pol2, nrhs, nside, nside2
    integer(i4b)                :: n, n2, npix, i, j, int, nmaps, nmaps2, ignore, ig1, ig2
    logical(lgt)                :: inv, inv2, fast, havemask

    integer(i4b), pointer,     dimension(:)   :: healvec
    real(dp),     allocatable, dimension(:,:) :: outmap, locmap1, locmap2 
    real(dp), allocatable, dimension(:,:) :: map1in, map2in, map2mask


    fast = .false.
    havemask = .false.

    ! Get parameters
    if (iargc() /= 7 .and. iargc() /= 8) then
       write(*,*) 'finalmap takes 6(7) parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, map1file)
       call getarg(3, cov1file)
       call getarg(4, map2file)
       call getarg(5, cov2file)
       call getarg(6, maskfile)
       call getarg(7, outprefix)
    end if
    
    if (iargc() == 8) then 
       call getarg(8, runmode)
       if (trim(runmode)=='fast') then
          fast = .true.
          if (info%myid==0)  write(*,*) 'Running in fast mode; might not work for bad matrices'
       end if
    end if

    ! Read input
    call read_cov_sc(unit, info, cov1file, n, ordering, pol, cov1, inv)
    call read_cov_sc(unit, info, cov2file, n2, ordering2, pol2, cov2, inv2)
    call assert(n       == n2,        "Size mismatch between cov1 and cov2")
    call assert(ordering == ordering2, "Ordering mismatch between cov1 and cov2")
    call assert(pol      == pol2,      "Polarisation mismatch between cov1 and cov2")
    call read_map(map1in, ordering2, map1file, nside=nside, nmap=nmaps)
    call assert(ordering == ordering2, "Ordering mismatch between cov and map1")
    call read_map(map2in, ordering2, map2file, nside=nside2, nmap=nmaps2)
    call assert(ordering == ordering2, "Ordering mismatch between cov and map2")
    call assert(nside    == nside2,    "Nside mismatch between map1 and map2")
    call assert(nmaps    == nmaps2,    "Nmaps mismatch between map1 and map2")
    call read_map(map2mask, ordering2, maskfile, nside=nside2, nmap=nmaps2)
    call assert(ordering == ordering2, "Ordering mismatch between cov and map2mask")
    call assert(nside    == nside2,    "Nside mismatch between maps and map2mask")
    call assert(nmaps    == nmaps2,    "Nmaps mismatch between maps and map2mask")
    npix = 12*nside**2
    ! Invert covariance matrices if needed
    if (inv) then
       if (info%myid==0)  write(*,*) 'Inverting cov1'
       if (fast) then
          call sc_invert(cov1)
          ig1 = 0
       else
          call invert_matrix_eigen_sc(cov1, outtext='cov1', kept=ig1)
          ig1 = n - ig1
       end if
    else
       ig1 = 0
    end if
    if (inv2) then
       if (info%myid==0)  write(*,*) 'Inverting cov2'
       if (fast) then
          call sc_invert(cov2)
          ig2 = 0
       else
          call invert_matrix_eigen_sc(cov2, outtext='cov2', kept=ig2)
          ig2 = n - ig2
       end if
    else
       ig2 = 0
    end if
    ignore = min(ig1,ig2)

    ! making scalaversion of inmaps
    if (info%myid==0) then
       allocate(locmap1(n,1))
       allocate(locmap2(n,1))
       call healmap2map(locmap1(:,1), map1in, map2mask, pol)
       call healmap2map(locmap2(:,1), map2in, map2mask, pol)
    end if
    call sc_alloc(map1, n, 1, info)
    call sc_alloc(map2, n, 1, info)
    call sc_set(map1, locmap1, 0)
    call sc_set(map2, locmap2, 0)
    if (info%myid==0) deallocate(locmap1, locmap2, map1in, map2in)

    if (map1%cglob /= map2%cglob) then 
       if (info%myid==root) write(*,*) 'Not an equal number of right-hand-sides. Quiting'
       call mpi_finalize(ierr)
       stop
    end if
    nrhs = map1%cglob
    if (nrhs > 1) call mkdir(trim(outprefix) // '_sims')
    simprefix = trim(outprefix) // '_sims/' // trim(outprefix )

    if (pol==1) then
       if (info%myid==0) write(*,*) 'Running on temperature data'
       nmaps = 1
    else if (pol==2) then
       if (info%myid==0) write(*,*) 'Running on polarisation data'
       nmaps = 3
    else 
       write(*,*) 'Unknown polarisation. Quiting'
       stop
    end if

    ! Putting maps to local root for outputting
    if (info%myid==0) allocate(locmap1(n, nrhs), locmap2(n2, nrhs))
    call sc_get(map1, locmap1, 0) 
    call sc_get(map2, locmap2, 0)   
    if (info%myid==0) then
       call map2mask2map2heal(nint(map2mask(:,nmaps)), healvec, n/pol)
       ! Writing redused map1 to file
       call write_maps(unit, locmap1, healvec, pol, npix, ordering, outprefix, simprefix, 'input_map1.fits')
    
       ! Writing redused map2 to file
       call write_maps(unit, locmap2, healvec, pol, npix, ordering, outprefix, simprefix, 'input_map2.fits')
   
       ! Writing true diffmap to file
       call write_maps(unit, (locmap1-locmap2)/2, healvec, pol, npix, ordering, outprefix, simprefix, 'true_diffmap.fits')
    
       ! Writing direct sum map to file
       call write_maps(unit, (locmap1+locmap2)/2, healvec, pol, npix, ordering, outprefix, simprefix, 'direct_summap.fits')
    
       ! Writing map2mask to file
       allocate(outmap(0:npix-1, nmaps))
       call map2mask2outmap(outmap, healvec, pol)   
       outfile = trim(outprefix) // '_map2mask.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Map2mask written to file = ', trim(outfile)
    
       ! Writing mask to file
       call mask2outmap(outmap, healvec, pol)   
       outfile = trim(outprefix) // '_mask.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Mask written to file = ', trim(outfile)
    
       deallocate(outmap, locmap1, locmap2)
    end if

    ! Produce and output final maps and covs
    call produce_finalmaps_sc(unit, info, map1, map2, cov1, cov2, healvec, outprefix, simprefix, ordering, pol, nmaps, npix, nrhs, fast, ignore)

  end subroutine make_genmap_sc

!-----------------------------------------------------------------------------------------------
! subroutine produce_finalmaps
!-----------------------------------------------------------------------------------------------

  subroutine produce_finalmaps_sc(unit, info, map1, map2, cov1, cov2, healvec, outprefix, simprefix, ordering, pol, nmaps, npix, nrhs, fast, ignore)
    implicit none

    integer(i4b),          intent(in) :: unit
    type(scinfo),          intent(in) :: info
    integer(i4b), pointer, intent(in) :: healvec(:)
    character(len=256),    intent(in) :: outprefix, simprefix
    integer(i4b),          intent(in) :: ordering, pol, nmaps, npix, nrhs, ignore
    logical(lgt),          intent(in) :: fast
    type(scalamat)                    :: cov1, cov2, map1, map2

    type(scalamat)              :: eigvals, cov, invcov, sqrtinvcov
    type(scalamat)              :: mapsum, mapdiff, inveigvals
    character(len=256)          :: outfile, outtext
    integer(i4b)                :: ierr, root=0
    integer(i4b)                :: i, j, n, ig1, ig2, ig
    logical(lgt)                :: inv
    real(dp)                    :: t1, t2
    real(dp),     allocatable, dimension(:,:) :: outmap, locmap1, locmap2 
    real(dp),     allocatable, dimension(:)   :: diag, chi_square

    n = cov1%cglob

    ! cov matrix for true diffmap is sum of cov1 and cov2
    call sc_alloc(cov, n, n, info)        
    do j=1,cov%cloc
       do i = 1, cov%rloc
          cov%data(i,j) = (cov1%data(i,j) + cov2%data(i,j))/4
       end do
    end do

    ! write cov matrix for true diffmap to file
    outfile = trim(outprefix) // '_N_diffmap.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, cov, .false.)
    if (info%myid==0) write(*,*) '* Cov for true diffmap written to file = ', trim(outfile)

    ! invert cov for true diffmap
    if (fast) then
       call sc_invert(cov)
    else
       call cpu_time(t1)
       call sc_alloc(eigvals, n, 1, cov%info)
       call sc_alloc(inveigvals, n, 1, cov%info)
       call sc_eigenvalue_decompose(cov, eigvals)
       call cpu_time(t2)
       if (info%myid==0) write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1 
       ! Finding invcov
       call matrix_from_eigen_sc(cov, eigvals, invmat=invcov, invmatvals=inveigvals, outtext='cov_truediff', cov=.true., ignore=ignore)
       ! Writing eigendec of invcov for true diffmap to file
       outfile = trim(outprefix) // '_inv_N_eigen_diffmap.unf'
       call write_eig_sc(unit, info, outfile, ordering, pol, cov, inveigvals, .true.)
       if (info%myid==0) write(*,*) '* Eigendec of invcov for true diffmap written to file = ', trim(outfile)
    end if

    ! write inverse cov matrix for true diffmap to file
    outfile = trim(outprefix) // '_inv_N_diffmap.unf'
    if (fast) then
       call write_covmatrix_sc(unit, info, outfile, ordering, pol, cov, .true.)
    else
       call write_covmatrix_sc(unit, info, outfile, ordering, pol, invcov, .true.)
       call sc_dealloc(invcov)
    end if
    if (info%myid==0) write(*,*) '* Invcov for true diffmap written to file = ', trim(outfile)
    
    ! chi_square check for direct diffmap
    call sc_alloc(mapdiff, n, nrhs, info)
    do j=1,mapdiff%cloc
       do i = 1, mapdiff%rloc
          mapdiff%data(i,j) = (map1%data(i,j) - map2%data(i,j))/2
       end do
    end do
    allocate(chi_square(nrhs))
    if (fast) then
       call chi_square_sc(mapdiff, cov, chi_square)
    else
       outtext=trim(outprefix) //'_true_diffmap'
       call chi_square_eigen_sc(mapdiff, cov, inveigvals, chi_square, outtext)
       call sc_dealloc(inveigvals)
    end if
    if (info%myid==0) call write_chis(unit, chi_square, n, outprefix, simprefix, 'chisq', 'direct diffmap')
   
    call sc_dealloc(mapdiff)
    call sc_dealloc(cov)

    !-------------------------------------------------------------------------------
    ! Making weighted maps
    !-------------------------------------------------------------------------------
      
    ! Inverting cov1 and cov2
    if (fast) then
       call sc_invert(cov1)
       call sc_invert(cov2)
       ig1 = 0
       ig2 = 0
    else
       outtext='cov_map1'
       call invert_matrix_eigen_sc(cov1, outtext=outtext, kept=ig1)
       outtext='cov_map2'
       call invert_matrix_eigen_sc(cov2, outtext=outtext, kept=ig2)
       ig1 = n - ig1
       ig2 = n - ig2
    end if
    ig = min(ignore, ig1, ig2)
    ! Invcov for final weighted map is sum of invcov1 and invcov2
    call sc_alloc(invcov, n, n, info)
    do j=1,invcov%cloc
       do i = 1, invcov%rloc
          invcov%data(i,j) = cov1%data(i,j) + cov2%data(i,j)
       end do
    end do
    ! Writing inverse cov matrix for weighted maps to file
    outfile = trim(outprefix) // '_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, invcov, .true.)
    if (info%myid==0) write(*,*) '* Invcov for final map written to file = ', trim(outfile)

    ! Inverting invcov to find cov
    if (fast) then    
       call sc_alloc(cov,    n, n, info)
       do j=1,cov%cloc
          do i = 1, cov%rloc
             cov%data(i,j) = invcov%data(i,j)
          end do
       end do
       call sc_invert(cov)
    else
       call cpu_time(t1)
       call sc_eigenvalue_decompose(invcov, eigvals)
       call cpu_time(t2)
       if (info%myid==0) write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1
       ! Writing eigendec of invcov for weighted maps to file
       outfile = trim(outprefix) // '_inv_N_eigen.unf'
       call write_eig_sc(unit, info, outfile, ordering, pol, invcov, eigvals, .true.)
       if (info%myid==0) write(*,*) '* Eigendec of invcov for final map written to file = ', trim(outfile)
       ! Finding and writing sqrt of inverse cov matrix for weighted maps to file
       call matrix_from_eigen_sc(invcov, eigvals, sqrtmat=cov, outtext='sqrtinvcov_weighted', ignore=ig)       
       outfile = trim(outprefix) // '_sqrt_inv_N.unf'
       call write_covmatrix_sc(unit, info, outfile, ordering, pol, cov, .true.)
       if (info%myid==0) write(*,*) '* Sqrt of invcov for final map written to file = ', trim(outfile)
       call sc_dealloc(cov)
       ! Finding and writing sqrt of cov matrix for weighted maps to file
       call matrix_from_eigen_sc(invcov, eigvals, sqrtinv=cov, outtext='sqrtcov_weighted', ignore=ig)       
       outfile = trim(outprefix) // '_sqrt_N.unf'
       call write_covmatrix_sc(unit, info, outfile, ordering, pol, cov, .false.)
       if (info%myid==0) write(*,*) '* Sqrt of cov for final map written to file = ', trim(outfile)
       call sc_dealloc(cov)
       ! Finding cov
       call matrix_from_eigen_sc(invcov, eigvals, invmat=cov, outtext='cov_weighted', ignore=ig)       
    end if
    ! Writing cov matrix for weighted maps to file
    outfile = trim(outprefix) // '_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, cov, .false.)    
    if (info%myid==0) write(*,*) '* Cov for weighted map written to file = ', trim(outfile)
   
    ! Finding summap and diffmap
    call weighted_maps_sc(info, map1, map2, cov1, cov2, cov, mapsum, mapdiff)
    if (info%myid==0) allocate(locmap1(n,nrhs), locmap2(n,nrhs)) 
    call sc_get(mapsum, locmap1, 0)   
    call sc_get(mapdiff, locmap2, 0)   
    if (info%myid==0) then
       ! Writing final map to file
       call write_maps(unit, locmap1, healvec, pol, npix, ordering, outprefix, simprefix, 'map.fits')
       ! Writing weighted diffmap to file
       call write_maps(unit, locmap2, healvec, pol, npix, ordering, outprefix, simprefix, 'weighted_diffmap.fits')
       deallocate(locmap1, locmap2) 
    end if

    ! chi_square check for weighted diffmap
    if (fast) then
       call chi_square_sc(mapdiff, invcov, chi_square)
    else
       outtext=trim(outprefix) //'_weighted_diffmap'
       call chi_square_eigen_sc(mapdiff, invcov, eigvals, chi_square, outtext)
       call sc_dealloc(eigvals)
    end if
    if (info%myid==0) call write_chis(unit, chi_square, n, outprefix, simprefix, 'chisq_weighted_diffmap', 'weighted diffmap') 

    ! Writing rms map to file
    allocate(diag(n))
    call  sc_get_diag(cov, diag)
    if (info%myid==0) then
       allocate(outmap(0:npix-1, nmaps))
       call rms2outmap(outmap, diag, healvec, pol)
       outfile = trim(outprefix) // '_rms.fits'
       call write_map(outmap, ordering, outfile)
       write(*,*) '* Noise rms map written to file = ', trim(outfile)   
       ! Writing default mask based on rms map to file
       call rms2mask(outmap, diag, healvec, pol)
       outfile = trim(outprefix) // '_rms_mask.fits'
       call write_map(outmap, ordering, outfile)
       write(*,*) '* Mask from rms map written to file = ', trim(outfile)   
       deallocate(outmap)
    end if
    deallocate(diag)

    ! Clean up
    call sc_dealloc(cov)
    call sc_dealloc(invcov)
    call sc_dealloc(mapsum)
    call sc_dealloc(mapdiff)

  end subroutine produce_finalmaps_sc

!-----------------------------------------------------------------------------------------------
! subroutine solve_equation
!-----------------------------------------------------------------------------------------------

  subroutine solve_equation_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)     :: amatin, bmatin, rhsin
    type(scalamat)     :: amat, bmat, rhs, map, cov, redcov, redmap, eigvals
    character(len=256) :: filename, outprefix, outfile, argname, outtext, maskfile, simprefix, rhsfile
    integer(i4b)       :: i, j, int, ordering, pol, nmaps, kept, ignore
    integer(i4b)       :: n, npix, npix_mask, nside, nrhs, n2, ordering2, pol2
    logical(lgt)       :: onlymap, rhs_override
    real(dp)           :: t1, t2
    integer(i4b), pointer,     dimension(:)   :: map2mask, mask, in2red, map2heal
    real(dp),     allocatable, dimension(:,:) :: outmap, locmap
    real(dp),     allocatable, dimension(:)   :: diag
    
    onlymap      = .false.
    rhs_override = .false.
    ! Get parameters
    if (iargc() < 4) then
       write(*,*) 'solve takes at least 3 parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, filename)
       call getarg(3, maskfile)
       call getarg(4, outprefix)
    end if
    i = 5
    do while(i <= iargc())
       call getarg(i, argname)
       select case(argname)
          case("onlymap")
             if (info%myid==0) write(*,*) 'Running in onlymap mode; not computing cov matrices'
             onlymap=.true.
          case("rhs")
             if (info%myid==0) write(*,*) 'Using specified rhs instead of the one in eqn'
             i=i+1; call getarg(i, rhsfile)
             rhs_override = .true.
          case default
             if (info%myid==0) write(*,*) 'Unrecognized option ' // trim(argname)
             call give_user_info_sc(info)
       end select
       i=i+1
    end do

    ! Read equation set
    call read_eqs_combo_sc(unit, info, filename, n, ordering, pol, amatin, bmatin, rhsin, npix, map2mask)
    if(rhs_override) then
       call sc_dealloc(rhsin)
       open(unit, file=rhsfile, form="unformatted", action="read", status="old")
       read(unit) n2
       read(unit) ordering2
       read(unit) pol2
       read(unit) nrhs
       call assert(all([n2,ordering2,pol2] == [n,ordering,pol]), "Rhs and eqn inconsistent! Rhs has n: " // trim(itoa(n2)) // " order: " // trim(itoa(ordering2)) // " pol: " // trim(itoa(pol2)))
       call sc_alloc(rhsin, n, nrhs, info)
       call sc_read(unit, rhsin, type=sp)
       close(unit)
       if(info%myid == 0) write(*,*) "Successfully read " // trim(itoa(nrhs)) // " rhses"
    end if
    call check_oneinput(info%myid, ordering, pol)
    nrhs = rhsin%cglob
    if (nrhs > 1) call mkdir(trim(outprefix) // '_sims')
    simprefix = trim(outprefix) // '_sims/' // trim(outprefix )

    if (pol==1) then
       if (info%myid==0) write(*,*) 'Running on temperature data'
       nmaps = 1
    else if (pol==2) then
       if (info%myid==0) write(*,*) 'Running on polarisation data'
       nmaps = 3
    else 
       write(*,*) 'Unknown polarisation. Quiting'
       stop
    end if
    ! Read mask
    if (trim(maskfile) /= 'nomask') then
       call read_and_check_mask(info%myid, maskfile, mask, nmaps, ordering, npix, pol)
    end if

    ! Removing unhit pixels
    allocate(diag(n))
    call sc_get_diag(amatin, diag)
    call unhit_remove(info%myid, map2mask, diag, in2red, pol)
    call inmap2red_sc(info, in2red, rhsin, amatin, rhs, amat, 1, bmatin, bmat)

    deallocate(diag)
    n = size(in2red)
    deallocate(in2red)

    ! Writing rhs to file
    allocate(locmap(n,nrhs))
    call sc_get(rhs, locmap, 0)
    if (info%myid==0) then
       call write_maps_m2m(unit, locmap, map2mask, pol, npix, ordering, outprefix, simprefix, 'rhs.fits')
    end if
    deallocate(locmap)

    ! Solve equation set
    outtext = 'map '
    if (onlymap) then
       call fast_solve_eqs_sc(n, amat, rhs, map, cov)
       ignore = 0
    else
       call solve_eqs_sc(info, unit, n, amat, rhs, map, cov, bmat=bmat, outtext=outtext, kept=kept)
       ignore = n-kept
    end if

    ! Writing input map to file
    allocate(locmap(n,nrhs))
    call sc_get(map, locmap, 0)
    if (info%myid==0) then
       call write_maps_m2m(unit, locmap, map2mask, pol, npix, ordering, outprefix, simprefix, 'inmap.fits')
    end if
    deallocate(locmap)

    ! Removing large and uncommon pixels from map2mask
    if (trim(maskfile) /= 'nomask') then
       call map2mask2in2red(info%myid, map2mask, in2red, map2heal, n, npix, ordering, mask)
    else
       call map2mask2in2red(info%myid, map2mask, in2red, map2heal, n, npix, ordering)
    end if
    n = pol*n
    if (info%myid==0) write(*,*) "Number of modes is",n
    ! ..and from map and cov
    call inmap2red_sc(info, in2red, map, cov, redmap, redcov, pol)
    deallocate(in2red)

    ! Putting map in locmap
    allocate(locmap(n,nrhs))
    call sc_get(redmap, locmap, 0)
    call sc_dealloc(redmap)
    if (info%myid==0) then
       ! Writing input map to file
       call write_maps(unit, locmap, map2heal, pol, npix, ordering, outprefix, simprefix, 'map.fits')
       
       ! Writing map2mask to file
       allocate(outmap(0:npix-1, nmaps))
       call map2mask2outmap(outmap, map2heal, pol)

       outfile = trim(outprefix) // '_map2mask.fits'
       call write_map(outmap, ordering, outfile)
       write(*,*) '* Map2mask written to file = ', trim(outfile)
    
       ! Writing mask to file
       call mask2outmap(outmap,map2heal, pol)   
       outfile = trim(outprefix) // '_mask.fits'
       call write_map(outmap, ordering, outfile)
       write(*,*) '* Mask written to file = ', trim(outfile)
    end if

    if (onlymap) return

    ! Writing cov matrix to file
    outfile = trim(outprefix) // '_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, redcov, .false.)
    if (info%myid==0) write(*,*) '* Covariance matrix written to file = ', trim(outfile)
    ! Writing rmsmap to file
    allocate(diag(n))
    call  sc_get_diag(redcov, diag)
    if (info%myid==0) then
       call rms2outmap(outmap, diag, map2heal, pol)
       outfile = trim(outprefix) // '_rms.fits'
       call write_map(outmap, ordering, outfile)
       write(*,*) '* Noise rms map written to file = ', trim(outfile)   
       ! Writing default rms-based mask to file
       call rms2mask(outmap, diag, map2heal, pol)
       outfile = trim(outprefix) // '_rms_mask.fits'
       call write_map(outmap, ordering, outfile)
       write(*,*) '* Mask from rms map written to file = ', trim(outfile)   
    end if
    deallocate(diag)

    ! Inverting cov to find invcov and sqrt of invcov
    call sc_alloc(eigvals, n, 1, info)
    call cpu_time(t1)
    call sc_eigenvalue_decompose(redcov, eigvals)
    call cpu_time(t2)
    if (info%myid==0) write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1
    call sc_alloc(amat, n, n, info)
    call sc_alloc(bmat, n, 1, info)
    ! Finding and writing inverse cov matrix to file
    call matrix_from_eigen_sc(redcov, eigvals, invmat=amat, invmatvals=bmat, outtext='invcov', cov=.true., ignore=ignore)  
    outfile = trim(outprefix) // '_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, amat, .true.)
    if (info%myid==0) write(*,*) '* Invcov written to file = ', trim(outfile)
    ! Writing eigendec of invcov to file
    outfile = trim(outprefix) // '_inv_N_eigen.unf'
    call write_eig_sc(unit, info, outfile, ordering, pol, redcov, bmat, .true.)
    if (info%myid==0) write(*,*) '* Eigendec of invcov written to file = ', trim(outfile)
    ! Finding and writing sqrt of inverse cov matrix to file
    call matrix_from_eigen_sc(redcov, eigvals, sqrtinv=amat, outtext='sqrtinv', cov=.true., ignore=ignore)       
    outfile = trim(outprefix) // '_sqrt_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, amat, .true.)
    if (info%myid==0) write(*,*) '* Sqrt of invcov written to file = ', trim(outfile)
    ! Finding and writing sqrt of cov matrix to file
    call matrix_from_eigen_sc(redcov, eigvals, sqrtmat=amat, outtext='sqrtcov', cov=.true., ignore=ignore)  
    outfile = trim(outprefix) // '_sqrt_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, amat, .false.)
    if (info%myid==0) write(*,*) '* Sqrt of cov written to file = ', trim(outfile)

    ! Clean up
    if (allocated(locmap))      deallocate(locmap)
    if (allocated(outmap))      deallocate(outmap)
    deallocate(map2mask)
    call sc_dealloc(amat)
    call sc_dealloc(bmat)
    call sc_dealloc(redcov)
    call sc_dealloc(eigvals)

  end subroutine solve_equation_sc

  !-----------------------------------------------------------------------------------------------
  ! subroutine project_out_lowmodes
  !-----------------------------------------------------------------------------------------------
  
  subroutine project_out_lowmodes_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)                        :: invcov, nmap, fmap, finf, invnf, ifti, mvec, nvec, eigvals
    character(len=256)                    :: infile, outprefix, map2maskfile, matrixfile
    character(len=256)                    :: outfile, cutoffname, nocovname
    integer(i4b)                          :: n, nmodes, ordering, pol, lmax, k, i, j
    integer(i8b)                          :: n_i8b, nmodes_i8b
    integer(i4b)                          :: nmaps, nside, tmporder, tmpnside, tmpnmaps, npix
    real(dp), dimension(1)                :: chi_square
    real(dp)                              :: t1, t2
    logical(lgt)                          :: inv, nocov, temperature
    
    real(dp), allocatable, dimension(:,:) :: inmap, outmap, map2mask, locmap
    
    ! Get parameters
    if (iargc() /= 6 .and. iargc() /= 7) then
       write(*,*) 'lcut takes 5(6) parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, map2maskfile)
       call getarg(3, matrixfile)
       call getarg(4, infile)
       call getarg(5, outprefix)
       call getarg(6, cutoffname)
       if (iargc() == 7) call getarg(7, nocovname)
    end if
    nocov = .false.
    if (trim(nocovname)== 'nocov') then 
       nocov=.true.
       if (info%myid==0)  write(*,*) 'Running with covariance matrix put to unity.'
       if (info%myid==0)  write(*,*) 'Just for making pretty pictures, not for further analysis'
    end if
    read(cutoffname,*) lmax
    outprefix = trim(outprefix) // '_lcut' //trim(cutoffname)    
        
    ! Read covariance matrix and check polarisation
    call read_cov_sc(unit, info, matrixfile, n, ordering, pol, invcov, inv)
    if (pol == 1) then
       temperature = .true.
       nmodes = (lmax+1)**2
       if (info%myid==0) write(*,fmt='(a,i5,a,a,a)') 'Projecting out up to',nmodes,' modes with l <= ',trim(cutoffname),' for temperature'  
    else if (pol == 2) then
       temperature = .false.
       nmodes = 2*(lmax+1)**2
       if (info%myid==0) write(*,fmt='(a,i5,a,a,a)') 'Projecting out up to',nmodes,' modes with l <= ',trim(cutoffname),' for polarisation'  
    else
       if (info%myid==0) write(*,*) "Neither temperature nor polarisation data. quiting"
       call mpi_finalize(ierr)
       stop
    end if

    ! Read map and map2mask
    call read_map(map2mask, tmporder, map2maskfile, nside=nside, nmap=nmaps)
    call assert(tmporder == ordering, "Ordering mismatch between invcov and map2mask!")
    call read_map(inmap, tmporder, infile, nside=tmpnside, nmap=tmpnmaps)
    call assert(tmporder == ordering, "Ordering mismatch between invcov and map!")
    call assert(tmpnside == nside,    "Nside mismatch between map and map2mask!")
    !call assert(tmpnmaps == nmaps,    "Nmaps mismatch between map and map2mask!")
    npix = 12*nside**2

    ! Make unity cov matrix in nocov mode
    if (nocov) then
       do j=1, invcov%cloc
          do i = 1, invcov%rloc
             invcov%data(i,j) = 0.d0
          end do
       end do
       do i = 1, n
          call sc_set_entry(invcov, i, i, 1.d0, 0)
       end do
    end if

    ! nmap = n version of inmap
    if (info%myid==0) then
       allocate(locmap(n,1))
       call healmap2map(locmap(:,1), inmap, map2mask, pol)
    end if
    call sc_alloc(nmap, n, 1, info)
    call sc_set(nmap, locmap, 0)
    if (info%myid==0) deallocate(locmap, inmap)

    ! chi_square for inmap
    call chi_square_sc(nmap, invcov, chi_square)
    if (info%myid==0) call write_chi_log(unit, chi_square(1), n, outprefix, helptext='before')

    ! putting the unwanted modes into fmap
    n_i8b=n
    nmodes_i8b=nmodes   
    if (info%myid==0) then
       allocate(locmap(n, nmodes))
       call get_lowmodes(locmap, lmax, n_i8b, nmodes_i8b, nmaps, nside, ordering, temperature, map2mask)
    end if
    call sc_alloc(fmap, n, nmodes, info)
    call sc_set(fmap, locmap, 0)
    if (info%myid==0) deallocate(locmap)

    ! invnf = invcov*fmap
    call sc_alloc(invnf, n, nmodes, info)
    call sc_matmul(invcov, fmap, invnf)
    
    ! finf = fmap^T*invcov*fmap
    call sc_alloc(finf, nmodes, nmodes, info)
    call sc_matmul(fmap, invnf, finf, transa='t')
      
    ! inverting finf
    call invert_matrix_eigen_sc(finf, unit, kept=k, acc=1d-4)
    if (info%myid==0) then
       write(*,*) nmodes-k, 'of', nmodes, 'eigenvalues removed.'
       write(*,*) 'Removing a total of', k, 'modes out of',n
    end if

    ! ifti = invfinf*invnf^T
    call sc_alloc(ifti, nmodes, n, info)
    call sc_matmul(finf, invnf, ifti, transb='t')
    call sc_dealloc(finf)

    ! inmap = inmap - f* (f^T*N^-1*f)^-1 * (f^T*N^-1*inmap)
    call sc_alloc(mvec, nmodes, 1, info)
    call sc_matmul(ifti, nmap, mvec)
    call sc_alloc(nvec, n, 1, info)
    call sc_matmul(fmap, mvec, nvec)
    do j=1,nmap%cloc
       do i = 1, nmap%rloc
          nmap%data(i,j) = nmap%data(i,j) - nvec%data(i,j)
       end do
    end do
    call sc_dealloc(mvec)
    call sc_dealloc(fmap)
    
    ! write map and removed modes to file
    if (info%myid==0) allocate(locmap(n,1))
    call sc_get(nmap, locmap, 0)
    if (info%myid==0) then
       allocate(outmap(0:npix-1, nmaps))
       call map2healmap(outmap, locmap(:,1), nint(map2mask(:,nmaps)), pol, npix)
       outfile = trim(outprefix) // '_map.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end if
    call sc_get(nvec, locmap, 0)
    if (info%myid==0) then
       call map2healmap(outmap, locmap(:,1), nint(map2mask(:,nmaps)), pol, npix)
       outfile = trim(outprefix) // '_removed_modes_map.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Removed part of map written to file = ', trim(outfile)
       deallocate(outmap, locmap)
    end if
    call sc_dealloc(nvec)

    if (nocov) return

    ! invcov = incov - N^-1 * f * (f^T*N^-1*f)^-1 * f^T * N^-1
    call sc_alloc(finf, n, n, info)
    call sc_matmul(invnf, ifti, finf)
    do j=1,invcov%cloc
       do i = 1, invcov%rloc
          invcov%data(i,j) = invcov%data(i,j) - finf%data(i,j)
       end do
    end do
    call sc_dealloc(finf)
    call sc_dealloc(ifti)
    call sc_dealloc(invnf)

    ! chi_square for outmap
    call chi_square_sc(nmap, invcov, chi_square)
    if (info%myid==0) call write_chi_log(unit, chi_square(1), n-k, outprefix, helptext='after')
    call sc_dealloc(nmap)

    if (trim(nocovname)== 'onlymap') return

    ! Writing new inverse covariance matrix to file
    outfile = trim(outprefix) // '_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, invcov, .true.)
    if (info%myid==0)   write(*,*) '* Inv cov matrix written to file = ', trim(outfile)

    ! Eigendecomposition of invcov
    call sc_alloc(eigvals, n, 1, info)
    call cpu_time(t1)
    call sc_eigenvalue_decompose(invcov, eigvals)
    call cpu_time(t2)
    if (info%myid==0) write(*,*) 'Time spent on eigenvalue decomposition:', t2-t1
    ! Writing eigendec of invcov to file
    outfile = trim(outprefix) // '_inv_N_eigen.unf'
    call write_eig_sc(unit, info, outfile, ordering, pol, invcov, eigvals, .true.)
    if (info%myid==0) write(*,*) '* Eigendec of invcov written to file = ', trim(outfile)

    ! Finding and writing cov matrix to file
    call matrix_from_eigen_sc(invcov, eigvals, invmat=finf, outtext='cov')       
    outfile = trim(outprefix) // '_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, finf, .false.)
    if (info%myid==0) write(*,*) '* Cov written to file = ', trim(outfile)
    call sc_dealloc(finf)
    ! Finding and writing sqrt of inverse cov matrix to file
    call matrix_from_eigen_sc(invcov, eigvals, sqrtmat=finf, outtext='sqrtinv')       
    outfile = trim(outprefix) // '_sqrt_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, finf, .true.)
    if (info%myid==0) write(*,*) '* Sqrt of invcov written to file = ', trim(outfile)
    call sc_dealloc(finf)
    ! Finding and writing sqrt of cov matrix to file
    call matrix_from_eigen_sc(invcov, eigvals, sqrtinv=finf, outtext='sqrtcov')       
    outfile = trim(outprefix) // '_sqrt_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, finf, .false.)
    if (info%myid==0) write(*,*) '* Sqrt of cov written to file = ', trim(outfile)
    call sc_dealloc(finf)
    call sc_dealloc(invcov)
    call sc_dealloc(eigvals)

  end subroutine project_out_lowmodes_sc

  !-----------------------------------------------------------------------------------------------
  ! subroutine project_out_eigenmodes
  !-----------------------------------------------------------------------------------------------
  
  subroutine project_out_eigenmodes_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)                        :: eigmat, eigvals, nmap, eigmap, mat
    character(len=256)                    :: map2maskfile, matrixfile, mapfile, outprefix
    character(len=256)                    :: outfile, cutoffname
    integer(i4b)                          :: n, nmodes, ordering, pol, cutoff, k, i, j
    integer(i8b)                          :: n_i8b, nmodes_i8b
    integer(i4b)                          :: nmaps, nside, tmporder, tmpnside, tmpnmaps, npix
    logical(lgt)                          :: inv
    real(dp)                              :: first_eigval
    
    real(dp), allocatable, dimension(:,:) :: inmap, outmap, map2mask, locmap, eigenvals
    
    ! Get parameters
    if (iargc() /= 6) then
       write(*,*) 'lcut takes 5 parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, map2maskfile)
       call getarg(3, matrixfile)
       call getarg(4, mapfile)
       call getarg(5, outprefix)
       call getarg(6, cutoffname)
    end if
    read(cutoffname,*) cutoff
    outprefix = trim(outprefix) // '_eigcut' //trim(cutoffname)    
        
    ! Read covariance matrix and check polarisation
    call read_eig_sc(unit, info, matrixfile, n, ordering, pol, eigmat, eigvals)
    if (pol == 1) then
       if (info%myid==0) write(*,fmt='(a,a,a)') 'Projecting out ',trim(cutoffname),' eigenmodes for temperature'  
    else if (pol == 2) then
       if (info%myid==0) write(*,fmt='(a,a,a)') 'Projecting out ',trim(cutoffname),' eigenmodes for polarisation'  
    else
       if (info%myid==0) write(*,*) "Neither temperature nor polarisation data. quiting"
       call mpi_finalize(ierr)
       stop
    end if

    ! Read map and map2mask
    call read_map(map2mask, tmporder, map2maskfile, nside=nside, nmap=nmaps)
    call assert(tmporder == ordering, "Ordering mismatch between invcov and map2mask!")
    call read_map(inmap, tmporder, mapfile, nside=tmpnside, nmap=tmpnmaps)
    call assert(tmporder == ordering, "Ordering mismatch between invcov and map!")
    call assert(tmpnside == nside,    "Nside mismatch between map and map2mask!")
    call assert(tmpnmaps == nmaps,    "Nmaps mismatch between map and map2mask!")
    npix = 12*nside**2

    ! nmap = n version of inmap
    if (info%myid==0) then
       allocate(locmap(n,1))
       call healmap2map(locmap(:,1), inmap, map2mask, pol)
    end if
    call sc_alloc(nmap, n, 1, info)
    call sc_set(nmap, locmap, 0)
    if (info%myid==0) deallocate(inmap)

    ! eigenversion of nmap
    call sc_alloc(eigmap, n, 1, info)
    call sc_matmul(eigmat, nmap, eigmap, transa='t')    
    if (info%myid==0) allocate(eigenvals(n,1))
    call sc_get(eigmap, locmap, 0)
    call sc_get(eigvals, eigenvals, 0)
    ! cutting bad eigenmodes
    if (info%myid==0) then 
       call check_eigenvals(info%myid, 'input', eigenvals(:,1), extravec=locmap(:,1))
       i = 1
       do while (eigenvals(i,1) == 0.d0)
          i = i+1
       end do
       write(*,*) eigenvals(i,1), eigenvals(cutoff,1), eigenvals(n,1), 'eigenvals'
       if (eigenvals(i,1) < eigenvals(n,1)) then
          if (info%myid==0) write(*,*) 'cutting beginning'
          do i = 1, cutoff
             locmap(i,1)    = 0.d0
             eigenvals(i,1) = 0.d0
          end do
       else
          if (info%myid==0) write(*,*) 'cutting end'
          do i = n-cutoff+1, n
             locmap(i,1)    = 0.d0
             eigenvals(i,1) = 0.d0
          end do
       end if
    end if
    call sc_set(eigmap, locmap, 0)
    call sc_set(eigvals, eigenvals, 0)
    call sc_matmul(eigmat, eigmap, nmap)
    ! Write eigencut map to file
    call sc_get(nmap, locmap, 0)
    if (info%myid==0) then
       allocate(outmap(0:npix-1, nmaps))
       call map2healmap(outmap, locmap(:,1), nint(map2mask(:,nmaps)), pol, npix)
       outfile = trim(outprefix) // '_map.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! Writing eigendec of invcov to file
    outfile = trim(outprefix) // '_inv_N_eigen.unf'
    call write_eig_sc(unit, info, outfile, ordering, pol, eigmat, eigvals, .true.)
    if (info%myid==0) write(*,*) '* Eigendec of invcov written to file = ', trim(outfile)

    ! Finding and writing new inverse covariance matrix to file
    call eigmult_sc(eigmat, eigvals, mat)
    outfile = trim(outprefix) // '_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat, .true.)
    if (info%myid==0)   write(*,*) '* Inv cov matrix written to file = ', trim(outfile)
    call sc_dealloc(mat)

    ! Finding and writing cov matrix to file
    call matrix_from_eigen_sc(eigmat, eigvals, invmat=mat, outtext='cov')       
    outfile = trim(outprefix) // '_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat, .false.)
    if (info%myid==0) write(*,*) '* Cov written to file = ', trim(outfile)
    call sc_dealloc(mat)

    ! Finding and writing sqrt of inverse cov matrix to file
    call matrix_from_eigen_sc(eigmat, eigvals, sqrtmat=mat, outtext='sqrtinv')       
    outfile = trim(outprefix) // '_sqrt_inv_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat, .true.)
    if (info%myid==0) write(*,*) '* Sqrt of invcov written to file = ', trim(outfile)
    call sc_dealloc(mat)

    ! Finding and writing sqrt of cov matrix to file
    call matrix_from_eigen_sc(eigmat, eigvals, sqrtinv=mat, outtext='sqrtcov')       
    outfile = trim(outprefix) // '_sqrt_N.unf'
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat, .false.)
    if (info%myid==0) write(*,*) '* Sqrt of cov written to file = ', trim(outfile)
    call sc_dealloc(mat)
    call sc_dealloc(eigmat)
    call sc_dealloc(eigvals)

  end subroutine project_out_eigenmodes_sc

!-----------------------------------------------------------------------------------------------
! subroutine get_sqrt_sc
!-----------------------------------------------------------------------------------------------

  subroutine get_sqrt_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)     :: invcov, sqrtinvcov
    character(len=256) :: filename, outprefix, outfile, argname, outtext
    integer(i4b)       :: i, int, ordering, pol
    integer(i4b)       :: n
    real(dp)           :: t1, t2
    logical(lgt)       :: inv
        
    ! Get parameters
    if (iargc() /= 3) then
       write(*,*) 'sqrt takes 2 parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, filename)
       call getarg(3, outprefix)
    end if

    ! Read covariance matrix
    call read_cov_sc(unit, info, filename, n, ordering, pol, invcov, inv)

    ! Get sqrt of invcov
    outtext = outprefix
    call sqrt_matrix_sc(unit, outtext, invcov, sqrtinvcov)

    ! Writing inverse cov matrix for weighted maps to file
    if (inv) then
       outfile = trim(outprefix) // '_sqrt_inv_N.unf'
       call write_covmatrix_sc(unit, info, outfile, ordering, pol, sqrtinvcov, .true.)
       if (info%myid==0) write(*,*) '* Sqrt of inv cov matrix written to file = ', trim(outfile)
    else
       outfile = trim(outprefix) // '_sqrt_N.unf'
       call write_covmatrix_sc(unit, info, outfile, ordering, pol, sqrtinvcov, .false.)
       if (info%myid==0) write(*,*) '* Sqrt of cov matrix written to file = ', trim(outfile)
    end if

    !Clean up
    call sc_dealloc(sqrtinvcov)

  end subroutine get_sqrt_sc

!-----------------------------------------------------------------------------------------------
! subroutine get_inverse_sc
!-----------------------------------------------------------------------------------------------

  subroutine get_inverse_sc(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)     :: mat
    character(len=256) :: filename, outprefix, outfile
    integer(i4b)       :: n, ordering, pol
    logical(lgt)       :: inv
        
    ! Get parameters
    if (iargc() /= 3) then
       write(*,*) 'invert takes 2 parameters'
       call give_user_info_sc(info)
    end if
    call getarg(2, filename)
    call getarg(3, outfile)

    ! Read input matrix
    call read_cov_sc(unit, info, filename, n, ordering, pol, mat, inv)

    ! Get inverse of invcov
    call invert_matrix_eigen_sc(mat, unit, outtext=outprefix)

    ! Writing inverse cov matrix for weighted maps to file
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat, .not. inv)

    !Clean up
    call sc_dealloc(mat)

  end subroutine get_inverse_sc

  subroutine operate_on_matrices_sc(unit, info, operation)
    implicit none

    integer(i4b),     intent(in) :: unit
    type(scinfo),     intent(in) :: info
    character(len=*), intent(in) :: operation

    type(scalamat)     :: mat1, mat2
    character(len=256) :: filename1, filename2, outprefix, outfile
    integer(i4b)       :: n, ordering, pol
    logical(lgt)       :: inv
        
    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'operate takes 4 parameters'
       call give_user_info_sc(info)
    end if
    call getarg(2, filename1)
    call getarg(3, filename2)
    call getarg(4, outfile)

    ! Read input matrices
    call read_cov_sc(unit, info, filename1, n, ordering, pol, mat1, inv)
    call read_cov_sc(unit, info, filename2, n, ordering, pol, mat2, inv)

    if (trim(operation) == '+') then
       mat1%data = mat1%data + mat2%data 
    else if (trim(operation) == '-') then
       mat1%data = mat1%data + mat2%data 
    else if (trim(operation) == '*') then
       mat1%data = mat1%data * mat2%data 
    else if (trim(operation) == '/') then
       mat1%data = mat1%data / mat2%data 
    else
       write(*,*) 'Unknown operation. Exiting.'
       stop
    end if

    ! Writing inverse cov matrix for weighted maps to file
    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat1, inv)

    !Clean up
    call sc_dealloc(mat1)
    call sc_dealloc(mat2)

  end subroutine operate_on_matrices_sc

  subroutine operate_on_matrix_with_scalar_sc(unit, info, operation)
    implicit none

    integer(i4b),     intent(in) :: unit
    type(scinfo),     intent(in) :: info
    character(len=*), intent(in) :: operation

    type(scalamat)     :: mat
    real(dp)            :: scalar
    character(len=256) :: filename, outfile, scalarname
    integer(i4b)       :: n, ordering, pol
    logical(lgt)       :: inv
        
    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'operate takes 4 parameters'
       call give_user_info_sc(info)
    end if
    call getarg(2, filename)
    call getarg(3, scalarname)
    call getarg(4, outfile)

    read(scalarname, *) scalar

    ! Read input matrices
    call read_cov_sc(unit, info, filename, n, ordering, pol, mat, inv)

    if (trim(operation) == '+') then
       mat%data = mat%data + scalar
    else if (trim(operation) == '-') then
       mat%data = mat%data + scalar
    else if (trim(operation) == '*') then
       mat%data = mat%data * scalar
    else if (trim(operation) == '/') then
       mat%data = mat%data / scalar
    else
       write(*,*) 'Unknown operation. Exiting.'
       stop
    end if

    call write_covmatrix_sc(unit, info, outfile, ordering, pol, mat, inv)

    !Clean up
    call sc_dealloc(mat)

   end subroutine

  subroutine compute_weighted_sum(info)
    implicit none

    type(scinfo),     intent(in) :: info

    type(scalamat)     :: mat, totmat, map, tmp, totmap
    character(len=256) :: filename, infofile, outprefix, outfile, map2maskname
    integer(i4b)       :: i, j, k, n, ordering, pol, unit, nummaps, order, nside, nmaps
    logical(lgt)       :: inv
    real(dp), allocatable, dimension(:,:) :: map_in, map2mask, locmap
    integer(i4b), pointer, dimension(:)   :: healvec
    character(len=256), dimension(100)    :: mapname, matname

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'weighted_sum takes 3 parameters'
       call give_user_info_sc(info)
    end if
    call getarg(2, infofile)
    call getarg(3, map2maskname)
    call getarg(4, outprefix)

    call read_map(map2mask, order, map2maskname, nside=nside, nmap=nmaps)

    unit  = getlun()
    open(unit,file=trim(infofile))
    do nummaps = 1, 100
       read(unit,*,end=99) mapname(nummaps), matname(nummaps)
    end do
99  close(unit)
    nummaps = nummaps-1

    ! Read first data set by hand
    call read_cov_sc(unit, info, matname(1), n, ordering, pol, mat, inv)
    call sc_alloc(totmat, n, n, info)
    call sc_copy(mat, totmat)

    if (info%myid==0) then
       allocate(locmap(n,1))
       call read_map(map_in, order, mapname(1), nside=nside, nmap=nmaps)
       call healmap2map(locmap(:,1), map_in, map2mask, pol)
    end if
    call sc_alloc(map, n, 1, info)
    call sc_alloc(tmp, n, 1, info)
    call sc_alloc(totmap, n, 1, info)
    call sc_set(map, locmap, 0)
    call sc_matmul(mat, map, totmap)

    do k = 2, nummaps

       write(*,*) 'Reading ', trim(matname(k))

       ! Add inverse covariance matrices
       call read_cov_sc(unit, info, matname(k), n, ordering, pol, mat, inv)
       do j = 1, mat%cloc
          do i = 1, mat%rloc
             totmat%data(i,j) = totmat%data(i,j) + mat%data(i,j)
          end do
       end do

       ! Add weighted maps
       if (info%myid==0) then
          deallocate(map_in)
          call read_map(map_in, order, mapname(k))
          call healmap2map(locmap(:,1), map_in, map2mask, pol)
       end if
       call sc_set(tmp, locmap, 0)
       call sc_matmul(mat, tmp, map)
       do i = 1, map%rloc
          totmap%data(i,1) = totmap%data(i,1) + map%data(i,1)
       end do
          
    end do

    ! Writing inverse cov matrix for weighted maps to file
    call write_covmatrix_sc(unit, info, trim(outprefix) // '_invN.unf', ordering, pol, totmat, .true.)

    ! Generate properly weighted map
    write(*,*) 'Inverting matrix'
    call sc_invert_matrix(totmat, 'LU')
    call sc_matmul(totmat, totmap, map)
    call sc_get(map, locmap, 0)

    if (info%myid == 0) then
       call map2healmap(map_in, locmap(:,1), nint(map2mask(:,nmaps)), pol, 12*nside**2)
       call write_map(map_in, ordering, trim(outprefix)//'_map.fits')
    end if
    
    !Clean up
    call sc_dealloc(totmat)
    call sc_dealloc(mat)
    
  end subroutine compute_weighted_sum

  !-----------------------------------------------------------------------------------------------
  ! subroutine extract_minimaps
  !-----------------------------------------------------------------------------------------------
  
  subroutine extract_minimaps(info)
    implicit none

    type(scinfo),    intent(in) :: info

    type(scalamat)                        :: invcov, nmap, fmap, finf, invnf, ifti, mvec, nvec, eigvals
    character(len=256)                    :: infile, outprefix, map2maskfile, matrixfile
    character(len=256)                    :: outfile, cutoffname, nocovname, par, minimapfile
    character(len=4)                      :: mapstring
    integer(i4b)                          :: n, nmodes, ordering, pol, lmax, k, i, j, nlist
    integer(i8b)                          :: n_i8b, nmodes_i8b
    integer(i4b)                          :: nmaps, nside, tmporder, tmpnside, tmpnmaps, npix, listpix(100000), pix(100000), mapnum
    real(dp), dimension(1)                :: chi_square
    real(dp)                              :: t1, t2, fwhm, a, b, lat, lon, cvec(3)
    logical(lgt)                          :: inv, nocov, temperature
    logical(lgt), allocatable, dimension(:) :: skip
    real(dp), allocatable, dimension(:,:) :: inmap, outmap, map2mask, locmap
    
    ! Get parameters
    if (iargc() /= 7) then
       write(*,*) 'extract_minimaps takes 7 parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, minimapfile)
       call getarg(3, par)
       read(par,*) fwhm
       call getarg(4, map2maskfile)
       call getarg(5, matrixfile)
       call getarg(6, infile)
       call getarg(7, outprefix)
    end if
        
    ! Read covariance matrix and data
    call read_cov_sc(unit, info, matrixfile, n, ordering, pol, invcov, inv)
    call read_map(map2mask, tmporder, map2maskfile, nside=nside, nmap=nmaps)
    call read_map(inmap, tmporder, infile)
    npix = 12*nside**2
    n    = n/2

    ! Output minimaps and covariance matrix
    mapnum = 0
    open(58,file=trim(minimapfile))
    do while (.true.)
       read(58,*,end=998) a, b, lon, lat
       mapnum = mapnum+1
       write(*,*) 'Writing minimap no.', mapnum, lon, lat
       lon    = lon*DEG2RAD
       lat    = 0.5d0*pi - lat*DEG2RAD 

       call ang2vec(lat, lon, cvec)
       call query_disc(nside, cvec, fwhm/60.d0*DEG2RAD, pix, nlist, nest=1)

       call int2string(mapnum, mapstring)

       outfile = trim(outprefix) // '_id' // mapstring // '_map.dat'
       open(59,file=trim(outfile))
       write(59,*) '# Nside = ', nside, ', Ordering = nest,    npix = ', nlist
       do j = 1, nlist
          write(59,*) pix(j), inmap(pix(j),2), ' Q'
       end do
       write(59,*)
       do j = 1, nlist
          write(59,*) pix(j), inmap(pix(j),3), ' U'
       end do
       close(59)

       allocate(skip(nlist))
       skip = (inmap(pix,2) == hpx_dbadval .or. inmap(pix,3) == hpx_dbadval)

       if (count(skip)>0) write(*,*) 'Found ', trim(itoa(count(skip))), ' masked pixels in minimap. These will not be included in cov file'

       outfile = trim(outprefix) // '_id' // mapstring // '_cov.dat'
       open(59,file=trim(outfile))
       write(59,*) '# Nside = ', nside, ', Ordering = nest,    npix = ', nlist
       do j = 1, nlist
          do k = 1, nlist
             if ((.not. skip(j)) .and. (.not. skip(k))) write(59,*) pix(j), pix(k), invcov%data(map2mask(pix(j),2), map2mask(pix(k),2)), ' QQ'
          end do
       end do
       write(59,*)
       do j = 1, nlist
          do k = 1, nlist
             if ((.not. skip(j)) .and. (.not. skip(k))) write(59,*) pix(j), pix(k), invcov%data(map2mask(pix(j),2), n+map2mask(pix(k),3)), ' QU'
          end do
       end do
       write(59,*)
       do j = 1, nlist
          do k = 1, nlist
             if ((.not. skip(j)) .and. (.not. skip(k))) write(59,*) pix(j), pix(k), invcov%data(n+map2mask(pix(j),3), map2mask(pix(k),2)), ' UQ'
          end do
       end do
       write(59,*)
       do j = 1, nlist
          do k = 1, nlist
             if ((.not. skip(j)) .and. (.not. skip(k))) write(59,*) pix(j), pix(k), invcov%data(n+map2mask(pix(j),3), n+map2mask(pix(k),3)), ' UU'
          end do
       end do
       close(59)
       deallocate(skip)

    end do
998 close(58)
    call sc_dealloc(invcov)

  end subroutine extract_minimaps


  !---------------------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------------------

  subroutine prune_multires(unit, info)
    implicit none

    type(scinfo),    intent(in) :: info
    character(len=512)          :: ifile, ofile
    integer(i4b)                :: nred, unit, i, m, n, pol, npix, ordering

    type(scalamat)              :: cov, icov, redcov, redicov, rhs, redrhs
    type(scalamat)              :: map, redmap
    integer(i4b), dimension(:), pointer :: map2mask, in2red, redmask2map, redmap2mask

    ! Get parameters
    if (iargc() /= 3) then
       write(*,*) 'prune takes 2 parameters'
       call give_user_info_sc(info)
    end if
    call getarg(2, ifile)
    call getarg(3, ofile)

    write(*,*) "Reading equation set"
    call read_eqs_sc(unit, info, ifile, n, ordering, pol, icov, rhs, npix, map2mask)
    write(*,*) "Solving for map and covariance matrix"
    call solve_eqs_sc(info, unit, n, icov, rhs, map, cov)
    write(*,*) "Mapping multiresolution pixels"
    call map2mask2in2red(info%myid, map2mask, in2red, redmask2map, m, npix, ordering)
    write(*,*) "Computing map2mask"
    n = m*pol
    allocate(redmap2mask(0:npix-1))
    redmap2mask = -1
    do i = 1, m
       redmap2mask(redmask2map(i)) = i
    end do
    write(*,*) "Removing multiresolution pixels"
    call inmap2red_sc(info, in2red, map, cov, redmap, redcov, pol)
    write(*,*) "Solving for rhs and inverse covariance matrix"
    call solve_eqs_sc(info, unit, n, redcov, redmap, redrhs, redicov)
    write(*,*) "Writing equation set"
    call write_eqs_sc(unit, info, ofile, n, ordering, pol, redicov, redrhs, npix, redmap2mask)
    write(*,*) "Done"
    deallocate(in2red, map2mask, redmask2map, redmap2mask)
  end subroutine

  subroutine convert_C_to_F90(unit, info)
    implicit none

    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    character(len=256) :: filename, infile, outfile, argname, outtext
    logical(lgt)       :: invcov
    integer(i4b)       :: i, int, ordering, pol
    integer(i4b)       :: n, unit_out 
    real(dp)           :: t1, t2, scale_factor, sigma_t, sigma_p
    real(dp), allocatable, dimension(:) :: M
        
    ! Get parameters
    if (iargc() /= 10) then
       write(*,*) 'convert_C_to_f90 takes 9 parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, infile)
       call getarg(3, outtext)
       read(outtext,*) n
       call getarg(4, outtext)
       read(outtext,*) ordering
       call getarg(5, outtext)
       read(outtext,*) pol
       call getarg(6, outtext)
       read(outtext,*) invcov
       call getarg(7, outtext)
       read(outtext,*) scale_factor
       call getarg(8, outtext)
       read(outtext,*) sigma_t
       call getarg(9, outtext)
       read(outtext,*) sigma_p
       call getarg(10, outfile)
    end if
    write(*,*) pol, '= pol,', real(sigma_t,sp), '= sigma_T,', real(sigma_p,sp), '= sigma_P,'
    write(*,*) n, '= n'

    unit_out = unit+1

    ! Read equation set from binary C file
    if (info%myid==0) write(*,*) '* Reading covariance matrix from file = ', trim(infile)
    allocate(M(n))
    open(unit,file=trim(infile),form='unformatted',access='direct', recl=8*n)
    open(unit_out,file=trim(outfile),form='unformatted')
    write(unit_out) n
    write(unit_out) ordering
    write(unit_out) pol
    if (pol==1) then
       do i = 1, n
          read(unit,rec=i) M
          M = M * scale_factor
          M(i) = M(i) + sigma_t**2
          if (i==1 .or. i==2 .or. i==n-1 .or. i==n .or. M(i)<1e-10) write(*,*) i, real(M(i),sp)
          write(unit_out) M
       end do
    else if (pol==3) then
       do i = 1, n/3
          read(unit,rec=i) M
          M = M * scale_factor
          M(i) = M(i) + sigma_t**2
          if (i==1 .or. i==2 .or. i==n/3-1 .or. i==n/3  .or. M(i)<1e-10) write(*,*) i, real(M(i),sp)
          write(unit_out) M
       end do
       do i = n/3+1, n
          read(unit,rec=i) M
          M = M * scale_factor
          M(i) = M(i) + sigma_p**2
          if (i==n/3+1 .or. i==n/3+2 .or. i==n-1 .or. i==n  .or. M(i)<1e-10) write(*,*) i, real(M(i),sp)
          write(unit_out) M
       end do
    end if
    write(unit_out) invcov
    close(unit_out)

    if (info%myid==0) write(*,*) '* Covariance matrix written to file = ', trim(outfile)

    ! Clean up
    deallocate(M)

  end subroutine convert_C_to_F90
  
  subroutine convert_WMAP_to_F90(unit, info)
    implicit none

    type(scinfo),    intent(in) :: info
    character(len=512)          :: ifile, ofile, remove_masked_name
    integer(i4b)                :: unit, n, ordering, polarization, status
    logical(lgt)                :: inv, remove_masked
    real(sp), pointer, dimension(:,:) :: mat_in
    logical(lgt), allocatable, dimension(:)  :: internal_mask
    type(scalamat)              :: mat
    integer(i4b)                :: num_notmasked, i, j, k, m

    ! Get parameters
    if (iargc() /= 3 .and. iargc() /= 4) then
       write(*,*) 'convert_WMAP_to_F90 takes two or three arguments!'
       call give_user_info_sc(info)
    end if
    call getarg(2, ifile)
    call getarg(3, ofile)
    remove_masked = .false.
    if (iargc() == 4) then
      call getarg(4, remove_masked_name)
      read(remove_masked_name, *) remove_masked
   end if

    ! Read covariance matrix on WMAP format
    print *, remove_masked
    call WMAP_Read_NInv(ifile, status, mat_in)
    if (remove_masked) then
       num_notmasked = size(mat_in, 1)
       allocate(internal_mask(size(mat_in, 1)))
       internal_mask = .True.
      do i = 1, size(mat_in, 1)
         print *, i
         if (mat_in(i, i) == 0) then
            num_notmasked = num_notmasked - 1
            internal_mask(i) = .False.
         end if
      end do
      n = num_notmasked
   else
      n = size(mat_in, 1)
    end if

    ! Set up the scalapost matrix
    call sc_alloc(mat, n, n, info)
    if (remove_masked) then
       k = 1
       do i = 1, size(mat_in, 1)
          m = 1
          if (internal_mask(i)) then
            do j = 1, size(mat_in, 1)
               if (internal_mask(j)) then
                  mat%data(m, k) = mat_in(j, i)
                  m = m + 1
               end if
            end do
            k = k + 1
         end if
      end do
    else
      mat%data = mat_in
   end if

    ! Output matrix
    call write_covmatrix_sc(unit, info, ofile, 2, 2, mat, .true.)
    call sc_dealloc(mat)
    deallocate(mat_in)
  end subroutine convert_WMAP_to_F90

  subroutine spost_invert_matrix(unit, info)
    implicit none

    type(scinfo),    intent(in) :: info
    character(len=512)          :: ifile, ofile, method
    integer(i4b)                :: unit, n, ordering, polarization
    logical(lgt)                :: inv
    type(scalamat)              :: mat

    real(dp), allocatable, dimension(:) :: eigenvals

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'invert takes three arguments!'
       call give_user_info_sc(info)
    end if
    call getarg(2, method)
    call getarg(3, ifile)
    call getarg(4, ofile)

    call read_cov_sc(unit, info, ifile, n, ordering, polarization, mat, inv)
    call sc_invert_matrix(mat, method)
    call write_covmatrix_sc(unit, info, ofile, ordering, polarization, mat, .not. inv)
    call sc_dealloc(mat)
  end subroutine spost_invert_matrix

  subroutine spost_mask_matrix(unit, info)
    implicit none

    type pixset
       integer(i4b), dimension(:), pointer :: pixels
    end type

    type(scinfo),    intent(in) :: info
    character(len=512)          :: imatfile, imaskfile, imap2maskfile, omatfile, omap2maskfile
    integer(i4b)                :: unit, n, ordering, fields, nmaps, nside, tmpnside, tmporder, tmpnmaps, npix, i, j, k, l, m
    integer(i4b), parameter     :: maxfields = 3
    integer(i4b)                :: ns_out(maxfields), ns_in(maxfields), n_out
    logical(lgt)                :: inv
    type(scalamat)              :: omat, imat
    real(dp),     dimension(:,:), allocatable :: imask, imap2mask, omap2mask, tmparr
    integer(i4b), dimension(:),   allocatable :: pixels
    integer(i4b), dimension(:,:), allocatable :: mask, map2mask
    integer(i4b), dimension(:),   pointer     :: map2heal
    type(pixset), dimension(maxfields)        :: fieldpix

    ! Get parameters
    if (iargc() /= 6) then
       write(*,*) 'mask takes 5 arguments!'
       call give_user_info_sc(info)
    end if
    call getarg(2, imatfile)
    call getarg(3, imap2maskfile)
    call getarg(4, imaskfile)
    call getarg(5, omatfile)
    call getarg(6, omap2maskfile)

    ! Ensuring pointer is not associated unless it is
    nullify(map2heal)

    ! Read matrix, map2mask and mask, convert to internal formats, and check.
    call read_cov_sc(unit, info, imatfile, n, ordering, fields, imat, inv)
    if (info%myid==0) write(*,*) '* Covariance matrix read from file = ', trim(imatfile)
    call read_map(imap2mask, tmporder, imap2maskfile, nside=nside, nmap=nmaps)
    call assert(tmporder == ordering, "Ordering mismatch between matrix and map2mask!")
    npix = 12*nside**2
    if(nmaps /= maxfields) then
       call assert(nmaps == 1, "Only 1 or 3 components supported for map2mask!")
       allocate(tmparr(npix,maxfields))
       do i = 1, maxfields
          tmparr(:,i) = imap2mask(:,1)
       end do
       deallocate(imap2mask)
       allocate(imap2mask(npix,maxfields))
       imap2mask = tmparr
       deallocate(tmparr)
    end if
    call read_map(imask, tmporder, imaskfile, nside=tmpnside, nmap=tmpnmaps)
    call assert(tmporder == ordering, "Ordering mismatch between matrix and mask!")
    call assert(tmpnside == nside,    "Nside mismatch between matrix and mask!")
    call assert(tmpnmaps == maxfields, "Expected " // trim(itoa(maxfields)) // "-component mask!")
    allocate(mask(0:npix-1,maxfields))
    allocate(map2mask(0:npix-1,maxfields))
    do i = 1, maxfields
       mask(:,i) = nint(imask(:,i))
       map2mask(:,i) = nint(imap2mask(:,i))
    end do
    deallocate(imask, imap2mask)
    if (info%myid==0) write(*,*) '* Map2mask and mask read from files'

    ! Everything has been read in. Time to build the list of remaining pixels.
    ! There could be a different map2mask and mask for each component. We
    ! translate each independently, and then stitch them together, using the
    ! arrays ns_in, holding the lengths of each field in the input, and ns_out
    ! doing the corresponding for the output.
    do i = 1, maxfields
       ns_in(i) = count(map2mask(:,i) > 0)
       if (ns_in(i) > 0) then
          call map2mask2in2red(info%myid, map2mask(:,i), fieldpix(i)%pixels, map2heal, ns_out(i), npix, ordering, mask(:,i))
       else
          ns_out(i) = 0
       end if
       if (info%myid==0) write(*,*) '* Processing field ', i, ' of ', maxfields, ' ns_out(i) = ', ns_out(i)
       map2mask(:,i) = -1
       do j = 1, ns_out(i)
          map2mask(map2heal(j),i) = j
       end do
       if (associated(map2heal)) deallocate(map2heal)
    end do

    ! NB! When the temperature layer is empty in the mask but not in map2mask this crashes.
    ! However this could be because map2mask did not match the map, so the bug is probably in finalmap
    ! (which made the map) and not here. To be worked out eventually. TMR
    n_out = sum(ns_out)
    allocate(pixels(n_out))
    j = 1
    k = 0
    do i = 1, maxfields
       if(ns_out(i) > 0) pixels(j:j+ns_out(i)-1) = fieldpix(i)%pixels + k
       j = j + ns_out(i)
       k = k + ns_in(i)
    end do
    if (info%myid==0) write(*,*) '* Translation done'
    ! Prune our matrix
    call sc_alloc(omat, n_out, n_out, info)
    call sc_copy(imat, omat, rows=pixels, cols=pixels)
    call sc_dealloc(imat)
    if (info%myid==0) write(*,*) '* Removed rows and columns'
    ! Construct output map2mask by converting to real
    allocate(omap2mask(0:npix-1, maxfields))
    omap2mask = real(map2mask,dp)

    ! And output results
    if(info%myid==0) then
       call write_map(omap2mask(:,1:nmaps), ordering, omap2maskfile)
    end if
    call write_covmatrix_sc(unit, info, omatfile, ordering, fields, omat, inv)
    call sc_dealloc(omat)
    if (info%myid==0) write(*,*) '* File written to file = ', trim(omatfile)
    deallocate(map2mask, mask, omap2mask, pixels)
    do i = 1, size(fieldpix)
       if (associated(fieldpix(i)%pixels)) deallocate(fieldpix(i)%pixels)
    end do
  end subroutine spost_mask_matrix


  subroutine cov2rms(unit, info)
    implicit none

    type(scinfo),    intent(in) :: info
    character(len=512)          :: imatfile, imaskfile, imap2maskfile, omatfile, omap2maskfile, ofile
    integer(i4b)                :: unit, n, ordering, fields, nmaps, nside, tmpnside, tmporder, tmpnmaps, npix, i, j, k, l, m
    integer(i4b), parameter     :: maxfields = 3
    integer(i4b)                :: ns_out(maxfields), ns_in(maxfields), n_out, p, offset
    logical(lgt)                :: inv
    type(scalamat)              :: omat, imat
    integer(i4b), allocatable, dimension(:,:) :: imap2mask
    real(dp),     allocatable, dimension(:,:) :: rms

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'cov2rms takes 3 arguments!'
       call give_user_info_sc(info)
    end if
    call getarg(2, imatfile)
    call getarg(3, imap2maskfile)
    call getarg(4, ofile)

    ! Read matrix, map2mask and mask, convert to internal formats, and check.
    call read_cov_sc(unit, info, imatfile, n, ordering, fields, imat, inv)
    if (info%myid==0) write(*,*) '* Covariance matrix read from file = ', trim(imatfile)
    call read_map(imap2mask, tmporder, imap2maskfile, nside=nside, nmap=nmaps)
    call assert(tmporder == ordering, "Ordering mismatch between matrix and map2mask!")
    npix = 12*nside**2
    if(nmaps /= maxfields) then
       call assert(nmaps == 1, "Only 1 or 3 components supported for map2mask!")
       stop
    end if
    allocate(rms(0:npix-1,nmaps))
    rms = 0.d0
    offset = 0
    do i = 1, nmaps
       do j = 0, npix-1
          if (imap2mask(j,i) > 0) then
             p = offset + imap2mask(j,i)
             rms(j,i) = sqrt(imat%data(p,p))
          end if
       end do
       offset = offset + count(imap2mask(:,i) > 0)
    end do

    ! And output results
    call write_map(rms, ordering, ofile)

  end subroutine


  subroutine rms2cov(unit, info)
    implicit none

    type(scinfo),    intent(in) :: info
    character(len=512)          :: irmsfile, ofile, input
    integer(i4b)                :: unit, n, ordering, nmaps, nside, i, j, k, ncomp
    real(dp)                    :: power
    real(dp),     allocatable, dimension(:,:) :: irms
    real(dp),     allocatable, dimension(:)   :: row

    unit = getlun()

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'rms2cov takes 3 arguments!'
       call give_user_info_sc(info)
    end if
    call getarg(2, irmsfile)
    call getarg(3, input)
    read(input,*) power
    call getarg(4, ofile)
    power = 2*power

    ! Read rms file
    call read_map(irms, ordering, irmsfile, nside=nside, nmap=nmaps)
    n = count(irms > 0.d0)
    ncomp = 0
    do i = 1, nmaps
       if (any(irms(:,i) > 0.d0)) ncomp = ncomp+1
    end do

    ! Write covariance matrix
    write(*,*) '    n        = ', n
    write(*,*) '    ordering = ', ordering
    write(*,*) '    ncomp    = ', ncomp
    open(unit, file=ofile, form="unformatted")
    write(unit) n
    write(unit) ordering
    write(unit) ncomp
    k = 0
    allocate(row(n))
    do j = 1, nmaps
       do i = 0, 12*nside**2-1
          if (irms(i,j) > 0.d0) then
             k   = k+1
             row = 0.d0
             row(k) = irms(i,j)**power
             write(unit) row
          end if
       end do
    end do
    write(unit) power < 0.d0
    close(unit)

    deallocate(irms, row)

  end subroutine

  subroutine cov2corr(unit, info)
    implicit none

    type(scinfo),    intent(in) :: info
    character(len=512)          :: imatfile, imaskfile, imap2maskfile, omatfile, omap2maskfile, ofile, param
    integer(i4b)                :: unit, n, ordering, fields, nmaps, nside, tmpnside, tmporder, tmpnmaps, npix, i, j, k, l, m
    integer(i4b), parameter     :: maxfields = 3
    integer(i4b)                :: ns_out(maxfields), ns_in(maxfields), n_out, p, offset, col
    logical(lgt)                :: inv
    type(scalamat)              :: omat, imat
    integer(i4b), allocatable, dimension(:,:) :: imap2mask

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'cov2corr takes 3 arguments!'
       call give_user_info_sc(info)
    end if
    call getarg(2, imatfile)
    call getarg(3, param)
    read(param,*) col
    call getarg(4, ofile)

    ! Read matrix, map2mask and mask, convert to internal formats, and check.
    call read_cov_sc(unit, info, imatfile, n, ordering, fields, imat, inv)
    if (info%myid==0) write(*,*) '* Covariance matrix read from file = ', trim(imatfile)
    unit = getlun()
    open(unit,file=trim(ofile))
    do i = 1, n
       write(unit,*) i, imat%data(i,col) / sqrt(imat%data(i,i)*imat%data(col,col))
    end do
    close(unit)

  end subroutine


  subroutine sc_invert_matrix(mat, method)
    implicit none
    type(scalamat)   :: mat
    character(Len=*) :: method

    if (method(1:2) == "LU") then
       call sc_invert(mat, symmetric=.false.)
    else if(method(1:3) == "eig") then
       call invert_matrix_eigen_sc(mat)
    else if(method(1:4) == "chol") then
       call sc_invert(mat, symmetric=.true.)
    else
       write(*,*) "Unrecognized method '" // trim(method) // "'!"
       stop
    end if

  end subroutine


  ! Given an inverse covariance matrix, map2mask, a solved map and a
  ! template, find the maximum likelihood amplitude of the template
  subroutine fit_template(unit, info)
    implicit none
    integer(i4b)       :: unit
    type(scinfo)       :: info

    character(len=512) :: icovfile, map2maskfile, mapfile, templatefile
    integer(i4b)       :: i, j, k, m, n, order, order2, nside, fields, nmaps
    integer(i4b)       :: nside2, offset
    logical(lgt)       :: inv
    type(scalamat)     :: icov, mmat, tmat
    integer(i4b), dimension(:),   allocatable :: map2mask
    real(dp),     dimension(:,:), allocatable :: map, template, imap2mask
    real(dp),     dimension(:),   allocatable :: marr, tarr, icmap, ictem
    real(dp) :: amplitude, sigma, var

    if(iargc() < 5) call give_user_info_sc(info)

    call getarg(2, icovfile)
    call getarg(3, map2maskfile)
    call getarg(4, mapfile)
    call getarg(5, templatefile)

    ! Read the data. Since an ordering problem is so easy to fix, we
    ! fix it ourselves.
    call read_cov_sc(unit, info, icovfile, n, order, fields, icov, inv)
    call read_map(imap2mask, order2, map2maskfile, nside=nside, nmap=nmaps)
    call set_ordering(order, order2, imap2mask)
    call read_map(map, order2, mapfile)
    call set_ordering(order, order2, map)
    call read_map(template, order2, templatefile)
    call set_ordering(order, order2, template)

    ! Extract the appropriate arrays from the maps
    offset = 0; if(fields == 2) offset = 1
    allocate(map2mask(0:size(imap2mask,1)-1))
    map2mask = nint(imap2mask(:,1+offset))

    deallocate(imap2mask)
    call map2array(map,      map2mask, fields, marr)
    call map2array(template, map2mask, fields, tarr)
    deallocate(map, template)

    if(any(healnot(marr))) write(*,*) "Warning: Unset pixels in map."
    if(any(healnot(tarr))) write(*,*) "Warning: Unset pixels in template."
    where(healnot(marr) .or. healnot(tarr)) marr = 0
    where(healnot(marr) .or. healnot(tarr)) tarr = 0

    ! The maximul likelihood amplitude is
    ! a = var*(t'*icov*m), var = 1/(t'*icov*t)
    allocate(icmap(n), ictem(n))
    call mul_smat_vec(icov, marr, icmap)
    call mul_smat_vec(icov, tarr, ictem)
    var = 1/(dot_product(tarr, ictem))
    amplitude = dot_product(tarr, icmap) * var
    sigma = sqrt(var)

    if(info%myid == 0) write(*,*) amplitude, sigma
    deallocate(icmap, ictem, tarr, marr, map2mask)
    call sc_dealloc(icov)
  end subroutine

  ! Given an inverse covariance matrix, map2mask, a solved map and a
  ! template, find the maximum likelihood amplitude of the template
  subroutine fit_template_with_noise(unit, info)
    implicit none
    integer(i4b)       :: unit
    type(scinfo)       :: info

    character(len=512) :: covfile, map2maskfile, mapfile, templatefile, param, rmsfile
    integer(i4b)       :: i, j, k, m, n, order, order2, nside, fields, nmaps
    integer(i4b)       :: nside2, offset, numbin, status
    real(dp)           :: alpha_min, alpha_max, dalpha, chisq, lndet, mean, stddev, ln_det_C, my_ln_det_C
    real(dp)           :: W_max, my_W_max
    logical(lgt)       :: inv
    type(scalamat)     :: cov, mmat, tmat, totmat, d, r, t, invC_r, W, tmp
    integer(i4b), dimension(:),   allocatable :: map2mask
    real(dp),     dimension(:,:), allocatable :: map, template, imap2mask, rms
    real(dp),     dimension(:),   allocatable :: marr, tarr, icmap, ictem, rmsarr
    real(dp),     dimension(:),   allocatable :: alpha, P
    real(dp) :: amplitude, sigma, var

    if(iargc() < 5) call give_user_info_sc(info)

    call getarg(2, covfile)
    call getarg(3, map2maskfile)
    call getarg(4, mapfile)
    call getarg(5, templatefile)
    call getarg(6, rmsfile)
    call getarg(7, param)
    read(param,*) alpha_min
    call getarg(8, param)
    read(param,*) alpha_max
    call getarg(9, param)
    read(param,*) numbin
    dalpha = (alpha_max - alpha_min) / real(numbin-1,dp)

    ! Read the data. Since an ordering problem is so easy to fix, we
    ! fix it ourselves.
    call read_cov_sc(unit, info, covfile, n, order, fields, cov, inv)
    call read_map(imap2mask, order2, map2maskfile, nside=nside, nmap=nmaps)
    call set_ordering(order, order2, imap2mask)
    call read_map(map, order2, mapfile)
    call set_ordering(order, order2, map)
    call read_map(template, order2, templatefile)
    call set_ordering(order, order2, template)
    call read_map(rms, order2, rmsfile)
    call set_ordering(order, order2, rms)

    ! Extract the appropriate arrays from the maps
    offset = 0; if(fields == 2) offset = 1
    allocate(map2mask(0:size(imap2mask,1)-1))
    map2mask = nint(imap2mask(:,1+offset))

    deallocate(imap2mask)
    call map2array(map,      map2mask, fields, marr)
    call map2array(template, map2mask, fields, tarr)
    call map2array(rms,      map2mask, fields, rmsarr)
    deallocate(map, template, rms)

    if(any(healnot(marr))) write(*,*) "Warning: Unset pixels in map."
    if(any(healnot(tarr))) write(*,*) "Warning: Unset pixels in template."
    if(any(healnot(tarr))) write(*,*) "Warning: Unset pixels in rmsmap."
    where(healnot(marr) .or. healnot(tarr)) marr = 0
    where(healnot(marr) .or. healnot(tarr)) tarr = 0
    where(healnot(marr) .or. healnot(tarr)) rmsarr = 0

    ! Set up data, template and template noise matrix
    call sc_alloc(tmat,   n, n, info)
    call sc_alloc(totmat, n, n, info)
    call sc_alloc(d,      n, 1, info)
    call sc_alloc(t,      n, 1, info)
    call sc_alloc(r,      n, 1, info)
    call sc_alloc(invC_r, n, 1, info)
    call sc_alloc(W,      n, 1, info)
    call sc_alloc(tmp,    n, 1, info)
    call sc_set(d, reshape(marr,(/n,1/)), 0)
    call sc_set(t, reshape(tarr,(/n,1/)), 0)
    tmat%data = 0.d0
    do i = 1, tmat%rloc
       do j = 1, tmat%cloc
          if (tmat%rmap(i) == tmat%cmap(j)) then
             tmat%data(i,j) = rmsarr(tmat%rmap(i))**2
!             cov%data(i,j) = cov%data(i,j) + 0.0001d0
          end if
       end do
    end do

    ! Compute a grid of alpha
    allocate(alpha(numbin))
    allocate(P(numbin))
    do i = 1, numbin
       ! Set up residual and covariance matrix
       alpha(i)    = alpha_min + real(i-1,dp) * dalpha
       if (size(totmat%data) > 0) totmat%data = cov%data! + alpha(i)**2 * tmat%data
       if (size(r%data) > 0) r%data = d%data - alpha(i) * t%data

       ! Compute likelihood 
!       call sc_cholesky_decompose(totmat,status=status)
       call sc_eigenvalue_decompose(totmat, W, 'l')

       ! Compute determinant
       my_ln_det_C = 0.d0
       if (size(W%data) > 0) then
          my_W_max = maxval(W%data)
          call mpi_allreduce(my_W_max, W_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, info%comm, ierr)
          do j = 1, W%rloc
             if (W%data(j,1) / W_max > 1.d-20) then
                my_ln_det_C = my_ln_det_C + log(W%data(j,1))
                W%data(j,1) = 1.d0 / W%data(j,1)
             else
                W%data(j,1) = 0.d0
             end if
          end do
       else
          my_W_max             = -1.d30
          call mpi_allreduce(my_W_max, W_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, info%comm, ierr)
       end if
       call mpi_allreduce(my_ln_det_C, lndet, 1, MPI_DOUBLE_PRECISION, MPI_SUM, info%comm, ierr)

       call sc_copy(r, invC_r)
       call sc_matmul(totmat, invC_r, tmp, transa='t')
       if (size(W%data) > 0) tmp%data = tmp%data * W%data
       call sc_matmul(totmat, tmp, invC_r)
       chisq = sc_dotprod(r, invC_r)
       P(i) = -0.5d0*(chisq + lndet)

       if (info%myid == 0) write(*,*) 'alpha = ', real(alpha(i),sp), ', lnL = ', real(P(i),sp), ', chisq = ', real(chisq,sp), ', lndet = ', real(lndet,sp)
    end do

    if (info%myid == 0) then
       P      = exp(P-maxval(P))
       P      = P / sum(P * dalpha)
       mean   = sum(alpha * P) * dalpha
       stddev = sqrt(sum((alpha-mean)**2 * P) * dalpha)
       write(*,*) 'alpha = ', mean, ' +/- ', stddev

       open(85,file='alpha.dat')
       do i = 1, numbin
          write(85,*) real(alpha(i),sp), real(P(i),sp)
       end do

    end if

    deallocate(P, alpha)
    call sc_dealloc(cov)
    call sc_dealloc(tmat)
    call sc_dealloc(totmat)
    call sc_dealloc(d)
    call sc_dealloc(t)
    call sc_dealloc(r)
    call sc_dealloc(invC_r)
  end subroutine

  ! Given an inverse covariance matrix, map2mask, a solved map and a
  ! template, fit a 2D gaussian
  subroutine fit_gaussian(unit, info)
    implicit none
    integer(i4b)       :: unit
    type(scinfo)       :: info

    character(len=512) :: icovfile, map2maskfile, mapfile, templatefile
    integer(i4b)       :: i, j, k, m, n, order, order2, nside, fields, nmaps
    integer(i4b)       :: nside2, offset, nstep, npar, seed, naccept, nreject
    logical(lgt)       :: inv
    type(scalamat)     :: icov, mmat, tmat
    real(dp)           :: cvec(3), lnL_curr, lnL_prop
    type(planck_rng)   :: rng
    integer(i4b), dimension(:),   allocatable :: map2mask
    real(dp),     dimension(:,:), allocatable :: map, template, imap2mask, vec, p
    real(dp),     dimension(:),   allocatable :: marr, tarr, icmap, ictem
    real(dp),     dimension(:),   allocatable :: dp, eta, p_curr, p_prop

    if(iargc() < 4) call give_user_info_sc(info)

    call getarg(2, icovfile)
    call getarg(3, map2maskfile)
    call getarg(4, mapfile)

    ! Read the data. Since an ordering problem is so easy to fix, we
    ! fix it ourselves.
    call read_cov_sc(unit, info, icovfile, n, order, fields, icov, inv)
    call read_map(imap2mask, order2, map2maskfile, nside=nside, nmap=nmaps)
    call set_ordering(order, order2, imap2mask)
    call read_map(map, order2, mapfile)
    call set_ordering(order, order2, map)

    ! Extract the appropriate arrays from the maps
    offset = 0; if(fields == 2) offset = 1
    allocate(map2mask(0:size(imap2mask,1)-1))
    map2mask = nint(imap2mask(:,1+offset))

    deallocate(imap2mask)
    call map2array(map, map2mask, fields, marr)
    deallocate(map)

    if(any(healnot(marr))) write(*,*) "Warning: Unset pixels in map."
    where(healnot(marr)) marr = 0

    ! Find position vectors
    allocate(vec(3,n))
    j = 0
    do i = 0, size(map2mask)-1
       if (map2mask(i)>0) then
          j = j+1
          if (order == 1) then
             call pix2vec_ring(nside, i, vec(:,j))
          else
             call pix2vec_nest(nside, i, vec(:,j))
          end if
       end if
    end do
    do i = 1, 3
       cvec(i) = sum(vec(i,:))
    end do
    cvec = cvec / sqrt(sum(cvec**2))

    ! Fit the 2D gaussian using MCMC
    nstep = 100000
    npar  = 5
    allocate(p(nstep,npar), p_curr(npar), p_prop(npar), dp(npar), eta(npar))
    p_curr(1) = maxval(marr)                                   ! Amplitude
    call vec2ang(cvec,p_curr(2),p_curr(3))                     ! Theta, phi 
    p_curr(4) = 10.d0 * DEG2RAD / 60.d0                        ! FWHM
    p_curr(5) = 0.d0                                           ! Offset
    lnL_curr = calc_2D_gauss_likelihood(info, p_curr, marr, vec, icov)

    dp(1) = 30.d0  ! Amplitude
    dp(2) = 1.d0 * DEG2RAD / 60.d0
    dp(3) = 1.d0 * DEG2RAD / 60.d0 
    dp(4) = 1.d0 * DEG2RAD / 60.d0 
    dp(5) = 10.d0
    !dp = 0.5d0 * dp

    seed = 3418947
    call rand_init(rng, seed)

    naccept = 0
    nreject = 0
    open(unit, file='gauss_par.dat', recl=1024)
    do i = 1, nstep
       do j = 1, npar
          p_prop(j) = p_curr(j) + dp(j) * rand_gauss(rng)
       end do
       lnL_prop = calc_2D_gauss_likelihood(info, p_prop, marr, vec, icov)

       if (exp(lnL_prop-lnL_curr) > rand_uni(rng)) then
          p_curr   = p_prop
          lnL_curr = lnL_prop
          naccept  = naccept + 1
       else
          nreject  = nreject + 1
       end if
       p(i,:)   = p_curr
!       write(unit,fmt='(i6,6f10.6)') i, lnL_curr, p(i,:)
!       write(*,fmt='(i6,6f10.6)')    i, lnL_curr, p(i,:)
       if (mod(i,10)==0) write(unit,*) i, (-2.d0*lnL_curr)/n, p(i,1), p(i,2)*RAD2DEG, p(i,3)*RAD2DEG, p(i,4)*RAD2DEG*60.d0, p(i,5)
       if (mod(i,1000)==0) write(*,*)    i, (-2.d0*lnL_curr)/n, p(i,1), p(i,2)*RAD2DEG, p(i,3)*RAD2DEG, p(i,4)*RAD2DEG*60.d0, p(i,5)

       if (i == 3000) then
          do j = 1, npar
             dp(j) = sqrt(variance(p(1:3000,j)))
          end do
       end if

    end do
    close(unit)

    write(*,*) 'accept rate = ', real(naccept,sp) / real(naccept+nreject,sp)

    deallocate(marr, map2mask)
    call sc_dealloc(icov)
  end subroutine

  function calc_2D_gauss_likelihood(info, p, map, vec, icov)
    implicit none

    type(scinfo)                         :: info
    real(dp), dimension(:),   intent(in) :: p, map
    real(dp), dimension(:,:), intent(in) :: vec
    type(scalamat)                       :: icov
    real(dp)                             :: calc_2D_gauss_likelihood

    integer(i4b) :: i, n
    real(dp)     :: cvec(3), r
    real(dp), allocatable, dimension(:) :: g

    call ang2vec(p(2), p(3), cvec)
    n = size(map)
    allocate(g(n))
!    open(58,file='data.dat')
    g = 0.d0
    do i = 1, n
       r    = acos(sum(vec(:,i)*cvec))
       if (r*RAD2DEG*60.d0 < 20.d0) then
          g(i) = p(1) * exp(-0.5d0 * r**2 / (p(4)/sqrt(8.d0*log(2.d0)))**2) + p(5)
       end if
!       write(58,*) map(i), g(i), icov%data(i,i)
!       write(58,*) acos(sum(vec(:,i)*cvec)), p(4)/sqrt(8.d0*log(2.d0))
!       write(58,*)
    end do
!    close(58)
!    stop
    calc_2D_gauss_likelihood = -0.5d0 * sum((map-g)* matmul(icov%data,(map-g))) - 0.5d0*(p(5)/30.d0)**2
!    stop

  end function calc_2D_gauss_likelihood

  subroutine concat_cov(unit, info)
    implicit none
    integer(i4b) :: unit, i
    type(scinfo) :: info
    character(len=512), dimension(:), allocatable :: ifiles
    character(len=512)                            :: ofile
    if(iargc() < 3) then
        if(info%myid == 0) write(*,*) "covcat takes at least 3 arguments!"
        call give_user_info_sc(info)
    end if
    allocate(ifiles(iargc()-2))
    do i = 1, iargc()-2
       call getarg(i+1,ifiles(i))
    end do
    call getarg(iargc(), ofile)
    call concat_cov_helper(info, ifiles, ofile)
    deallocate(ifiles)
  end subroutine

  ! Simple, hugely memory-wasting implemenation: Read in all matrices, allocate
  ! large one, copy over, write and deallocate.
  subroutine concat_cov_helper(info, ifiles, ofile)
    implicit none
    type(scinfo)     :: info
    character(len=*) :: ifiles(:), ofile
    integer(i4b)     :: unit, n_in(size(ifiles)), p_in(size(ifiles)), o_in(size(ifiles)), nfile, n, p
    integer(i4b)     :: i, j
    logical(lgt)     :: inv_in(size(ifiles))
    type(scalamat)   :: imats(size(ifiles)), omat

    nfile = size(ifiles)
    unit = getlun()
    do i = 1, nfile
       call read_cov_sc(unit, info, ifiles(i), n_in(i), o_in(i), p_in(i), imats(i), inv_in(i))
    end do
    call assert(all(n_in/p_in == n_in(1)/p_in(1)), "Inconsistent number of pixels!")
    call assert(all(o_in == o_in(1)), "Inconsistent orderings!")
    call assert(all(inv_in .eqv. inv_in(1)), "Inconsistent inverseness!")
    n = n_in(1)/p_in(1)
    p = sum(p_in)
    call sc_alloc(omat, n*p, n*p, info)
    omat%data = 0
    j = 0
    do i = 1, nfile
       call sc_copy(imats(i), omat, rect=[j+1,j+1,j+n_in(i),j+n_in(i)],which=2)
       call sc_dealloc(imats(i))
       j = j+n_in(i)
    end do
    call write_covmatrix_sc(unit, info, ofile, o_in(1), p, omat, inv_in(1))
    call sc_dealloc(omat)
  end subroutine

  subroutine eigendecomp_cmd(info)
    implicit none
    type(scinfo)          :: info
    character(len=512)    :: ifile, ofile, eigvalfile
    integer(i4b)          :: unit, n, order, ncomp, inv
    type(scalamat)        :: V, M
    real(dp), allocatable :: eigvals(:,:)
    if(iargc() < 3) call give_user_info_sc(info)
    call getarg(2, ifile)
    call getarg(3, ofile)
    if(iargc() > 3) call getarg(4, eigvalfile)

    if(info%myid == 0) call dmem("Start",0)
    unit = getlun()
    open(unit, file=ifile, form="unformatted", status="old", action="read")
    read(unit) n
    read(unit) order
    read(unit) ncomp
    call sc_alloc(V, n, n, info)
    call sc_alloc(M, n, 1, info)
    call sc_read(unit, V)
    close(unit)
    if(info%myid == 0) call dmem("Read " // trim(ifile),0)

    call sc_eigenvalue_decompose(V, M)
    if(info%myid == 0) call dmem("Eigendecomp",0)
    allocate(eigvals(n,1))
    call sc_get(M, eigvals, 0)

    open(unit, file=ofile, form="unformatted")
    if(info%myid == 0) then
       write(unit) n
       write(unit) order
       write(unit) ncomp
    end if
    call sc_write(unit, V)
    call sc_write(unit, M)
    close(unit)
    if(info%myid == 0 .and. iargc() > 3) call dump_matrix(eigvals, eigvalfile)
    if(info%myid == 0) call dmem("Wrote result",0)
    deallocate(eigvals)
    call sc_dealloc(V)
    call sc_dealloc(M)
  end subroutine

  ! Remove the least significant eigenmodes from a map.
  ! Specifically, cut those with inverse eigenvalues
  ! below a specified threshold.
  ! Send in the eigendecomposition of inv_N togehter
  ! with the map to be cut and the fraction
  subroutine eigencut_cmd(info)
    implicit none
    type(scinfo)              :: info
    character(len=512)        :: map2maskfile, invNeigfile, imapfile(2), omapfile(2), arr
    integer(i4b)              :: unit, n, order, order2, ncomp, inv, off, i, cut, npix
    integer(i4b)              :: pos(2), nlim, ilim, nside, mode, order3, mcomp, ni
    integer(i4b), parameter   :: mode_cut = 0, mode_swap = 1
    type(scalamat)            :: V, M, Map, Maps, Emap, Emaps
    real(dp)                  :: tval, mval
    real(dp),     allocatable :: limits(:), maps_loc(:,:), emaps_loc(:,:)
    real(dp),     allocatable :: Mloc(:,:), omap(:,:,:), cmap(:,:,:)
    integer(i4b), allocatable :: pixels(:)
    if(iargc() < 7) call give_user_info_sc(info)
    n = 1
    call getarg(n, arr);          n=n+1
    select case(arr)
       case("eigencut");   mode = mode_cut
       case("eigenswap");  mode = mode_swap
       case default;       call assert(.false., "Bug in eigencut_cmd!")
    end select
    call getarg(n, map2maskfile); n=n+1
    call getarg(n, invNeigfile);  n=n+1
    call getarg(n, imapfile(1));     n=n+1
    if(mode == mode_swap) then; call getarg(n, imapfile(2)); n=n+1; endif
    call getarg(n, omapfile(1));     n=n+1
    call getarg(n, omapfile(2));     n=n+1
    nlim = iargc()-n+1
    allocate(limits(nlim))
    do i = 1, nlim
       call getarg(n, arr); read(arr,*) limits(i); n=n+1
    end do
    ni = 1; if(mode == mode_swap) ni = 2

    if(info%myid == 0) call dmem("Start", 0)
    call read_eigs_scala(info, invNeigfile, M, V, ncomp, order, verbose=.true.)
    call read_pixels(map2maskfile, pixels, nside)
    call read_maps_scala(info, imapfile(1:ni), ncomp, order, nside, pixels, Map, verbose=.true.)
    n = size(pixels)*ncomp

    ! Get spread M to everybody
    allocate(Mloc(n,1))
    call sc_allget(M, Mloc)
    call sc_alloc(Emap,  n, 2,    info)
    call sc_alloc(Maps,  n, 2*nlim, info)
    call sc_alloc(Emaps, n, 2*nlim, info)

    ! Find the eigenspace version of map
    call sc_matmul(V, Map, Emap, transa='T')
    if(info%myid == 0) call dmem("Decomposed map")
    ! And get a local version of it
    allocate(emaps_loc(n, 2*nlim))
    call sc_get(Emap, emaps_loc(:,1:1),           0, [1,1,n,1])
    call sc_get(Emap, emaps_loc(:,nlim+1:nlim+1), 0, [1,2,n,2])
    if(info%myid == 0) then
       emaps_loc(:,2:nlim)  = spread(emaps_loc(:,1),     2,nlim-1)
       emaps_loc(:,nlim+2:) = spread(emaps_loc(:,nlim+1),2,nlim-1)
       ! Swap chosen eigenvalues
       mval = maxval(Mloc(:,1))
       do ilim = 1, nlim
          cut = 0
          do i = 1, n
             if(Mloc(i,1) >= mval*limits(ilim)) cycle
             tval = emaps_loc(i,nlim+ilim)
             emaps_loc(i,nlim+ilim) = emaps_loc(i,ilim)
             emaps_loc(i,ilim) = tval
             cut = cut+1
          end do
          if(mode == mode_cut)  call dmem("Removed " // trim(itoa(cut)) // " eigvals", 0)
          if(mode == mode_swap) call dmem("Swapped " // trim(itoa(cut)) // " eigvals", 0)
       end do
    end if
    call sc_set(Emaps, emaps_loc, 0)

    ! Transform back to normal space
    call sc_matmul(V, Emaps, Maps, transa='N')
    if(info%myid == 0) call dmem("Recomposed maps", 0)

    ! And output it
    allocate(maps_loc(n, 2*nlim))
    call sc_get(Maps, maps_loc, 0)
    off = 1; if(ncomp==2) off = 2
    mcomp = 1; if(ncomp > 1) mcomp = 3
    if(info%myid == 0) then
       allocate(omap(n/ncomp,mcomp,nlim))
       allocate(cmap(n/ncomp,mcomp,nlim))
       omap = hpx_dbadval
       cmap = hpx_dbadval
       omap(:,off:off+ncomp-1,:)    = reshape(maps_loc(:,:nlim),   [n/ncomp,ncomp,nlim])
       cmap(:,off:off+ncomp-1,:)    = reshape(maps_loc(:,nlim+1:), [n/ncomp,ncomp,nlim])
       call write_map(omap, pixels, nside, order, omapfile(1))
       call write_map(cmap, pixels, nside, order, omapfile(2))
       deallocate(omap, cmap)
       call dmem("Wrote result", 0)
    end if

    call sc_dealloc(V)
    call sc_dealloc(M)
    call sc_dealloc(Map)
    call sc_dealloc(Maps)
    call sc_dealloc(Emap)
    call sc_dealloc(Emaps)

    deallocate(Mloc, emaps_loc, maps_loc, pixels)
  end subroutine

  ! Calculate the chisquare per eigenmode of a map
  subroutine eigenchi_cmd(info)
    implicit none
    type(scinfo)              :: info
    character(len=512)        :: map2maskfile, invNeigfile, imapfile, ofile
    integer(i4b)              :: unit, n, order, ncomp, nside, i
    type(scalamat)            :: V, M, Map, Emap
    real(dp),     allocatable :: chisq(:), eigmap(:,:), eigvals(:,:)
    integer(i4b), allocatable :: pixels(:)
    if(iargc() < 5) call give_user_info_sc(info)
    call getarg(2, map2maskfile)
    call getarg(3, invNeigfile)
    call getarg(4, imapfile)
    call getarg(5, ofile)

    if(info%myid == 0) call dmem("Start", 0)
    call read_eigs_scala(info, invNeigfile, M, V, ncomp, order)
    if(info%myid == 0) call dmem("Read matrix", 0)
    call read_pixels(map2maskfile, pixels, nside)
    call read_maps_scala(info,[imapfile], ncomp, order, nside, pixels,Map)
    n = size(pixels)*ncomp
    if(info%myid == 0) call dmem("Read maps", 0)

    ! Eigendecompose map and calculate chisq per egienmode
    call sc_alloc(Emap, n, 1, info)
    call sc_matmul(V, Map, Emap, transa='T')
    if(info%myid == 0) call dmem("Decomposed map", 0)
    allocate(chisq(n), eigvals(n,1), eigmap(n,1))
    call sc_get(M,    eigvals, 0)
    call sc_get(Emap, eigmap,  0)
    chisq = eigmap(:,1)**2*eigvals(:,1)
    if(info%myid == 0) call dmem("Calculated chisq", 0)

    ! And output
    if(info%myid == 0) then
       unit = getlun()
       open(unit, file=ofile)
       do i = 1, n
          write(unit,'(i5,3e15.7)') i, chisq(i), eigvals(i,1), eigmap(i,1)
       end do
       close(unit)
       call dmem("Wrote output",0)
    end if
    call sc_dealloc(Map)
    call sc_dealloc(V)
    call sc_dealloc(M)
    call sc_dealloc(Emap)
    deallocate(eigvals, eigmap, chisq, pixels)
  end subroutine



  subroutine convert_mask_to_map2mask(info)
    implicit none
    type(scinfo)              :: info
    character(len=512)        :: maskfile, map2maskfile

    integer(i4b) :: nside, npix, i, j, k, nmaps, ordering
    real(dp), allocatable, dimension(:,:) :: mask, map2mask

    ! Get parameters
    if (iargc() /= 3) then
       write(*,*) 'convert_mask_to_map2mask takes 3(4) parameters'
       call give_user_info_sc(info)
    else 
       call getarg(2, maskfile)
       call getarg(3, map2maskfile)
    end if

    call read_map(mask, ordering, maskfile)
    npix  = size(mask,1)
    nmaps = size(mask,2)

    allocate(map2mask(0:npix-1,nmaps))
    map2mask = -1.d0
    do i = 1, nmaps
       k = 1
       do j = 0, npix-1
          if (mask(j,i) > 0.d0) then
             map2mask(j,i) = k
             k             = k+1
          end if
       end do
    end do

    call write_map(map2mask, ordering, map2maskfile)

    deallocate(map2mask)

  end subroutine convert_mask_to_map2mask

  subroutine read_pixels(map2maskfile, pixels, nside)
    implicit none
    type(scinfo)                    :: info
    character(len=*), intent(in)    :: map2maskfile
    integer(i4b),     allocatable   :: pixels(:)
    integer(i4b),     intent(out)   :: nside
    integer(i4b)                    :: order, pos(2), npix, n, i, j
    integer(i4b),     allocatable   :: map2mask(:,:)
    if(allocated(pixels)) deallocate(pixels)
    ! The map2mask file is in a stupid format with 3 components,
    ! where the first one may be empty. So read in all three and use
    ! first non-empty column
    call read_map(map2mask, order, map2maskfile)
    pos = maxloc(map2mask)
    npix = size(map2mask,1)
    nside= npix2nside(npix)
    n = maxval(map2mask(:,pos(2)))
    allocate(pixels(n))
    do i = 0, npix-1
       if(map2mask(i,pos(2)) <= 0) cycle
       pixels(map2mask(i,pos(2))) = i
    end do
  end subroutine

  ! Read an eigenmatrixfile
  subroutine read_eigs_scala(info, eigfile, eigvals, eigvecs, ncomp, order, verbose)
    implicit none
    type(scinfo)                    :: info
    character(len=*), intent(in)    :: eigfile
    type(scalamat),   intent(inout) :: eigvals, eigvecs
    integer(i4b),     intent(out)   :: ncomp, order
    logical(lgt),     optional      :: verbose
    logical(lgt)                    :: v
    integer(i4b)                    :: unit, n
    v = .false.; if(present(verbose)) v = verbose
    call sc_dealloc(eigvals)
    call sc_dealloc(eigvecs)
    unit = getlun()
    open(unit, file=eigfile, form="unformatted", action="read", status="old")
    read(unit) n
    read(unit) order
    read(unit) ncomp
    call sc_alloc(eigvecs, n, n, info)
    call sc_alloc(eigvals, n, 1, info)
    call sc_read(unit, eigvecs)
    call sc_read(unit, eigvals)
    close(unit)
    if(v .and. info%myid == 0) call dmem("Read matrix", 0)
  end subroutine

  ! Read a set of maps, extract T, QU or TQU depending on ncomp, and
  ! the pixels corresponding to pixels, and produce a scalamat with
  ! one column for each map. Only reads the first map from each file.
  ! Checks that the nside and order match those given.
  subroutine read_maps_scala(info, fnames, ncomp, order, nside, pixels, maps, verbose)
    implicit none
    type(scinfo)                    :: info
    character(len=*), intent(in)    :: fnames(:)
    integer(i4b),     intent(in)    :: ncomp, order, nside, pixels(:)
    type(scalamat),   intent(inout) :: maps
    logical(lgt),     optional      :: verbose
    integer(i4b)                    :: off, n, nmap, i, order2, npix, npix2
    logical(lgt)                    :: v
    real(dp),         allocatable   :: imap(:,:)
    call sc_dealloc(maps)
    off = 1; if(ncomp == 2) off = 2
    n   = size(pixels)*ncomp
    nmap= size(fnames)
    npix= 12*nside**2
    v   = .false.; if(present(verbose)) v = verbose
    call sc_alloc(maps, n, nmap, info)
    do i = 1, nmap
       call read_map(imap, order2, fnames(i))
       if(v .and. info%myid == 0) call dmem("Read map " // trim(itoa(i)))
       npix2 = size(imap,1)
       call assert(order == order2, "Inconsistent ordering in map " // trim(itoa(i)))
       call assert(npix  == npix2,  "Inconsistent nside for map " // trim(itoa(i)))
       call assert(size(imap,2) >= ncomp, "Too few components in map " // trim(itoa(i)))
       call sc_allset(maps, reshape(imap(pixels,off:off+ncomp-1),[n,1]), rect=[1,i,n,i])
       if(v .and. info%myid == 0) call dmem("Distributed map " // trim(itoa(i)))
       deallocate(imap)
    end do
  end subroutine

  subroutine co_add_maps(unit, info)
    implicit none
    integer(i4b),    intent(in) :: unit
    type(scinfo),    intent(in) :: info

    type(scalamat)     :: cov1,cov2,cov3,map1,map2
    real(dp), dimension(:,:), allocatable :: healmap1, healmap2, locmap, map2mask
    real(dp), dimension(:), allocatable :: totmap1, totmap2
    character(len=256) :: covfile1, covfile2,mapfile1,mapfile2, outmap1, outmap2, outsummap, covsumfile
    integer(i4b)       :: i, npi, n, n2, ordering, ordering2, pol, pol2, nside, nside2, nmaps, nmaps2
    logical(lgt)       :: inv, inv2
    integer(i4b), dimension(:),allocatable :: pixels1, pixels2!, map2mask
        
    ! Get parameters
    if (iargc() /= 9) then
       write(*,*) 'co_add takes 8 parameters'
       call give_user_info_sc(info)
    end if
    call getarg(2, covfile1)
    call getarg(3, covfile2)
    call getarg(4, covsumfile)
    call getarg(5, mapfile1)
    call getarg(6, mapfile2)
    call getarg(7, outmap1)
    call getarg(8, outmap2)
    call getarg(9, outsummap)

    ! Read input
    call read_cov_sc(unit, info, covfile1, n, ordering, pol, cov1, inv)
    call read_cov_sc(unit, info, covfile2, n2, ordering2, pol2, cov2, inv2)
    call assert(n       == n2,        "Size mismatch between cov1 and cov2")
    call assert(ordering == ordering2, "Ordering mismatch between cov1 and cov2")
    call assert(pol      == pol2,      "Polarisation mismatch between cov1 and cov2")
    call read_map(healmap1, pixels1, nside, ordering2, mapfile1)
    call assert(ordering == ordering2, "Ordering mismatch between cov and map1")
    call read_map(healmap2, pixels2,nside, ordering2, mapfile2)
    call assert(ordering == ordering2, "Ordering mismatch between cov and map2")
    write(*,*) "finished reading input"

    write(*,*)size(healmap1,1), size(healmap1,2)
    npi=size(healmap1,1)
    allocate(totmap1(n),totmap2(n))
    totmap1(1:npi)=healmap1(:,2)
    totmap2(1:npi)=healmap2(:,2)
    totmap1(npi+1:n)=healmap1(:,3)
    totmap2(npi+1:n)=healmap2(:,3)

    write(*,*)"A",totmap1(1:20)
    write(*,*)"B",totmap2(1:20)
    write(*,*)"C",cov1%data(1:20,1)
    write(*,*)"D",cov2%data(1:20,1)

    totmap1 = matmul(cov1%data , totmap1)
    totmap2 = matmul(cov2%data , totmap2)
    write(*,*)"E",cov1%data(1,1)

!    !add matrices
!    cov1%data = cov1%data + cov2%data

    !clean up
    call sc_dealloc(cov2)    
    call sc_dealloc(cov1)

    !read inverse-sum-cov
    call read_cov_sc(unit, info, covsumfile, n, ordering, pol, cov3, inv)

    totmap1 = matmul(cov3%data , totmap1)
    totmap2 = matmul(cov3%data , totmap2)

    healmap1(:,2) = totmap1(1:npi)
    healmap2(:,2) = totmap2(1:npi)
    healmap1(:,3) = totmap1(npi+1:n)
    healmap2(:,3) = totmap2(npi+1:n)

    !write out the two maps
    call write_map(healmap1, pixels1,nside,ordering,outmap1)
    call write_map(healmap2, pixels1,nside,ordering,outmap2)

    !make the sum map and write out
    healmap1=healmap1+healmap2
    call write_map(healmap1, pixels1,nside,ordering,outsummap)

    call sc_dealloc(cov3)    

  end subroutine co_add_maps

  subroutine test(unit, info)
    implicit none
    integer(i4b)       :: unit, n
    type(scinfo)       :: info
    type(scalamat)     :: A, B, C
    real(dp)           :: Af(3,3), Vf(3,3), Ef(3,1)

    n = 3
    Af = reshape(&
     & [8, 4, 2, &
     &  4,16, 4, &
     &  2, 4, 7 ], [3, 3])
    call sc_alloc(A, n, n, info, 1)
    call sc_alloc(B, n, 1, info, 1)
    call sc_allset(A, Af)
    call sc_eigenvalue_decompose(A, B)
    call sc_allget(A, Vf)
    call sc_allget(B, Ef)

    if(info%myid == 0) then
       write(*,*) "Vectors:"
       write(*,fmt="(3e15.7)") transpose(Vf)
       write(*,*) "Eigvals:"
       write(*,fmt="(3e15.7)") transpose(Ef)
    end if
  end subroutine

  ! Given two maps with corresponding covariance matrices, coadd them
  ! at the resolution of the high-res map. The high-res map should be
  ! specified first. In order to produce a final result with a
  ! simple beam, per-pixel weighting will not be used, only global
  ! weighting. The maps are assumed to be sparse, with the same
  ! pixels as the covariance matrices!
  subroutine upadd(info)
    implicit none
    type(scinfo)              :: info
    character(len=512)        :: hmfile, lmfile, hcfile, lcfile, omfile, ocfile
    integer(i4b)              :: hnside, lnside, order, order2, unit, hn, ln
    integer(i4b)              :: nc, nc2, o, step, i1, i2, j1, j2, i, j, c1, c2, k
    integer(i4b)              :: h1, h2, l1, l2
    real(dp)                  :: hw, lw, wnorm, hval, lval, g
    type(scalamat)            :: hcov, lcov
    real(dp),     allocatable :: hmap(:,:), lmap(:,:), omap(:,:), var(:), col(:)
    integer(i4b), allocatable :: hpixels(:), lpixels(:), hmap2mask(:), lmap2mask(:)
    call dset(id=info%myid)
    call dmem("start")
    call getarg(2, hmfile)
    call getarg(3, hcfile)
    call getarg(4, lmfile)
    call getarg(5, lcfile)
    call getarg(6, omfile)
    call getarg(7, ocfile)
    call read_map(hmap, hpixels, hnside, order,  hmfile)
    call read_map(lmap, lpixels, lnside, order2, lmfile)
    call assert(order == nest, "Ordering should be nested!")
    call assert(order == order2, "Inconsistent ordering in maps!")
    call assert(hnside >= lnside, "Highest res map should be first!")
    call dmem("read maps")

    ! Compute the mapping from high to low pixels. This is simply an
    ! integer pixel division
    step = (hnside/lnside)**2
    allocate(hmap2mask(0:12*hnside**2-1), lmap2mask(0:12*lnside**2-1))
    hmap2mask = 0
    lmap2mask = 0
    do i = 1, size(hpixels)
       hmap2mask(hpixels(i)) = i
    end do
    do i = 1, size(lpixels)
       lmap2mask(lpixels(i)) = i
    end do
    call assert(all(lmap2mask(hpixels/step)>0), "There must be a corresponding map 2 pixel for every map 1 pixel!")

    unit = getlun()
    call read_cov_sc(unit, info, hcfile, hn, order2, nc,  hcov)
    call assert(hn/nc == size(hpixels), "Map 1 and cov 1 have inconsistent pixels!")
    call assert(order == order2, "Map 1 and cov 1 have inconsistent ordering!")

    ! Read the low-res matrix in as a high-res matrix compatible with hcov
    call sc_alloc(lcov, hn, hn, info)
    open(unit,file=lcfile,form="unformatted",action="read",status="old")
    read(unit) ln
    read(unit) order2
    read(unit) nc2
    allocate(col(ln))
    call assert(ln/nc == size(lpixels), "Map 2 and cov 2 have inconsistent pixels!")
    call assert(order == order2, "Map 2 and cov 2 have inconsistent ordering!")
    call assert(nc == nc2, "Set 1 and 2 have inconsistent number of components!")
    o = 0; if(nc == 2) o = 1
    call assert(size(hmap,2) >= o+nc, "Map 1 has too few components!")
    call assert(size(lmap,2) >= o+nc, "Map 2 has too few components!")

    lcov%data = hcov%data
    do l1 = 1, ln
       ! Loop through low-res columns
       read(unit) col
       j  = modulo(l1-1,ln/nc)+1
       c1 = (l1-1)/(ln/nc)
       ! Loop through high-res columns in those low-res columns
       do j1 = lpixels(j)*step, (lpixels(j)+1)*step-1
          h1 = hmap2mask(j1)
          if(h1 <= 0) cycle
          h1 = h1 + c1*size(hpixels)
          if(lcov%cown(h1) /= lcov%info%gridpos(2)) cycle
          ! Loop through the local high-res rows we want to fill
          do i = 1, lcov%rloc
             h2 = lcov%rmap(i)
             ! Calculate the corresponding low-res position
             k  = modulo(h2-1,hn/nc)+1
             c2 = (h2-1)/(hn/nc)
             l2 = lmap2mask(hpixels(k)/step)
             if(l2 <= 0) cycle
             l2 = l2 + c2*size(lpixels)
             lcov%data(i,lcov%cpos(h1)) = col(l2)
          end do
       end do
    end do
    deallocate(col)
    close(unit)

    call dmem("read covs")

    if(info%myid==0) write(*,*) "Hack! Scaling map by 1000 and ant2thermo (and its cov by 1e6)"
    g = 1e3*ant2thermo(44d0)
    lmap(:,o:) = lmap(:,o:) * g
    lcov%data  = lcov%data  * g**2

    ! Ok, input is finished. Now determine the weights of the maps. Since
    ! we are not doing this per pixel, this will have to be a compromise.
    ! I choose to maximize the singal level in the center of the map, so
    ! I will simply use the 1/minval(var) as the weights.
    allocate(var(hn))
    call sc_get_diag(hcov, var)
    hw = 1/minval(var)
    call sc_get_diag(lcov, var)
    lw = 1/minval(var)/step
    deallocate(var)
    wnorm = hw+lw
    hw = hw/wnorm
    lw = lw/wnorm
    if(info%myid == 0) write(*,'(a,2f10.7)') "Weights: ", hw, lw

    ! Build the merged map. This will have the same pixels as the high
    ! resolution map. Any low res pixels outside of this will be discarded.
    ! We require all the high res pixels to be inside the low res ones.
    allocate(omap(size(hpixels),nc+o))
    omap = hmap
    do j = 1, nc
       do i = 1, size(hpixels)
          k = lmap2mask(hpixels(i)/step)
          if(k <= 0) cycle
          omap(i,j+o) = hmap(i,j+o)*hw + lmap(k,j+o)*lw
       end do
    end do
    if(info%myid == 0) call write_map(omap, hpixels, hnside, order, omfile)
    call dmem("combined maps")

    ! Build the merged covariance matrix. This is
    ! hw**2 * hcov + lw**2 * lcov. When reading. We have already
    ! upgraded lcov while reading, so this can be done directly now
    hcov%data = hcov%data * hw**2 + lcov%data * lw**2

    call dmem("combined covs")
    call write_covmatrix_sc(unit, info, ocfile, order, nc, hcov, .false.)
    call dmem("wrote covs")

    call sc_dealloc(hcov)
    call sc_dealloc(lcov)
  end subroutine

  ! Read in a matrix, and use it to produce other specified
  ! versions of the matrix. Allocates at most 3 matrices
  ! internally at any time, meaning that the memory requirements
  ! are about 3 times the size of one matrix.
  subroutine buildmats(info)
    implicit none
    type(scinfo)       :: info
    character(len=512) :: arg, fname
    type(scalamat)     :: eigvec, eigval, eigval2, mat, mat2
    integer(i4b)       :: ai, i, j, k, m, n, unit, order, nc
    unit = getlun()
    call assert(iargc() >= 5, "Syntax: buildmats type imat type omat [type omat ...]")
    ai = 2
    call getarg(ai, arg)
    call getarg(ai+1, fname)
    if(info%myid == 0) write(*,*) "Reading " // trim(fname)
    select case(arg)
       case("eig")
          call read_eig_sc(unit, info, fname, n, order, nc, eigvec, eigval)
       case("ieig")
          call read_eig_sc(unit, info, fname, n, order, nc, eigvec, eigval)
          where(eigval%data /= 0) eigval%data = 1/eigval%data
       case("cov")
          call read_cov_sc(unit, info, fname, n, order, nc, eigvec)
          call sc_alloc(eigval, n, 1, info)
          call sc_eigenvalue_decompose(eigvec, eigval)
       case("icov")
          call read_cov_sc(unit, info, fname, n, order, nc, eigvec)
          call sc_alloc(eigval, n, 1, info)
          call sc_eigenvalue_decompose(eigvec, eigval)
          where(eigval%data /= 0) eigval%data = 1/eigval%data
       case default
          call assert(.false., "Unrecognized input format '" // trim(arg) // "'!")
    end select
    where(eigval%data < 0) eigval%data = 0
    call sc_alloc(eigval2, n, 1, info)
    call sc_alloc(mat,     n, n, info)
    ! Then loop through the outputs
    ai = ai+2
    do while(ai+1 <= iargc())
       call getarg(ai, arg)
       call getarg(ai+1,fname)
       if(info%myid == 0) write(*,*) "Writing " // trim(arg) // ": " // trim(fname)
       select case(arg)
          case("eig")
             call write_eig_sc(unit, info, fname, order, nc, eigvec, eigval, .false.)
          case("ieig")
             eigval2%data = 0
             where(eigval%data /= 0) eigval2%data = 1/eigval%data
             call write_eig_sc(unit, info, fname, order, nc, eigvec, eigval2, .true.)
          case("cov","mat")
             call sc_alloc(mat2, n, n, info)
             call sc_matmul_diag(eigvec, eigval, mat2)
             call sc_matmul(mat2, eigvec, mat, transb='t')
             call sc_dealloc(mat2)
             call write_covmatrix_sc(unit, info, fname, order, nc, mat, .false.)
          case("icov","imat")
             eigval2%data = 0
             where(eigval%data /= 0) eigval2%data = 1/eigval%data
             call sc_alloc(mat2, n, n, info)
             call sc_matmul_diag(eigvec, eigval2, mat2)
             call sc_matmul(mat2, eigvec, mat, transb='t')
             call sc_dealloc(mat2)
             call write_covmatrix_sc(unit, info, fname, order, nc, mat, .true.)
          case("sqrt")
             eigval2%data = eigval%data**0.5
             call sc_alloc(mat2, n, n, info)
             call sc_matmul_diag(eigvec, eigval, mat2)
             call sc_matmul(mat2, eigvec, mat, transb='t')
             call sc_dealloc(mat2)
             call write_covmatrix_sc(unit, info, fname, order, nc, mat, .false.)
          case("isqrt")
             where(eigval%data /= 0) eigval2%data = eigval%data**(-0.5)
             call sc_alloc(mat2, n, n, info)
             call sc_matmul_diag(eigvec, eigval, mat2)
             call sc_matmul(mat2, eigvec, mat, transb='t')
             call sc_dealloc(mat2)
             call write_covmatrix_sc(unit, info, fname, order, nc, mat, .true.)
          case default
             if(info%myid == 0) write(*,*) "  Ignored unknown type " // trim(arg)
       end select
       ai = ai+2
    end do
    call sc_dealloc(eigvec)
    call sc_dealloc(eigval)
    call sc_dealloc(mat)
  end subroutine

  ! apply_mask isparse icov mask osparse ocov
  subroutine apply_mask(info)
    ! TMR TODO: Rewrite to ask which mask component to use (or just use all?)
    ! Make applying to cov optional.
    type(scinfo)               :: info
    character(len=512)         :: imapfile, icovfile, maskfile, omapfile, ocovfile
    integer(i4b)               :: nside, nside2, order, iunit, ounit, n, nc, order2, m
    integer(i4b)               :: i, j, k, c
    logical(lgt)               :: isinv
    integer(i4b),  allocatable :: ipix(:), opix(:), o2i(:)
    real(dp),      allocatable :: imap(:,:), omap(:,:), mask(:,:), icol(:), ocol(:)
    if(info%myid /= 0) return
    call dmem("start")
    call getarg(2, imapfile)
    call getarg(3, icovfile)
    call getarg(4, maskfile)
    call getarg(5, omapfile)
    call getarg(6, ocovfile)
    call read_map(imap, ipix, nside, order,  imapfile)
    call read_map(mask,              order2, maskfile)
    call dmem("read maps")
    nside2 = npix2nside(size(mask,1))
    call assert(nside == nside2, "Inconsistent nside in imap and mask")
    call assert(order == order2, "Inconsistent order in imap and mask")
    where(mask <= 0); mask = 0
    elsewhere;        mask = 1
    end where

    ! Create a mapping from reduced to full
    ! TMR: This only used T component. Changed to Q, but ideally it should use all to account for possibly different masks for each component. 
    ! also, mask is real, == test unsafe?
    call wherei(mask(ipix,2) == 1, o2i)
    m = size(o2i)
    allocate(opix(m))
    opix = ipix(o2i)

    ! Write output map
    allocate(omap(m,size(imap,2)))
    omap = imap(o2i,:)
    call write_map(omap, opix, nside, order, omapfile)
    call dmem("wrote map")

    ! Now read and write while filtering out unwanted pixels
    iunit = getlun()
    open(iunit,file=icovfile,action="read",status="old",form="unformatted")
    read(iunit) n
    read(iunit) order2
    read(iunit) nc
    call assert(size(ipix)*nc == n, "Inconsistent size of imap and icov")
    call assert(order == order2, "Inconsistent order of imap and icov")
    ounit = getlun()
    open(ounit,file=ocovfile,form="unformatted")
    write(ounit) m*nc
    write(ounit) order
    write(ounit) nc
    allocate(icol(n),ocol(m*nc))
    do i = 1, n
       j = modulo(i-1,n/nc)+1
       if(mask(ipix(j),1) == 0) then
          read(iunit)
       else
          read(iunit) icol
          do c = 0, nc-1
             ocol(1+c*m:m+c*m) = icol(o2i+c*n/nc)
          end do
          write(ounit) ocol
       end if
       if(modulo(i,100)==0) write(*,*) "col", i
    end do
    isinv = .false.
    read(iunit,end=1) isinv
    1 write(ounit) isinv
    close(iunit)
    close(ounit)
    call dmem("wrote cov")
  end subroutine

  subroutine build_S_matrix(info)
    implicit none
    type(scinfo)                  :: info
    character(len=512)            :: psfile, mapfile, bfile, arg, ofile, icfile
    integer(i4b)                  :: unit, i, j, k, m, n, lmax, l, nc, ncomp, ncomp2
    integer(i4b)                  :: nside, order, order2
    logical(lgt)                  :: inv, invert
    integer(i4b), allocatable     :: comps(:), pixels(:)
    real(dp),     allocatable     :: pslin(:,:), psmat(:,:,:), beam(:,:), map(:,:)
    real(dp)                      :: foo(6), bsize, C(3,3)
    type(angcorr_type)            :: cfun
    type(scalamat)                :: cov
    if(info%myid==0) call dmem("start")
    call getarg(2, psfile)
    call getarg(3, bfile)
    call getarg(4, mapfile)
    call getarg(5, arg);    read(arg,*) ncomp
    call getarg(6, arg);    invert = .false.; if(arg=="inv") invert = .true.
    call getarg(7, ofile)

    ! Read the linear power spectrum
    unit = getlun()
    open(unit,file=psfile,action="read",status="old")
    lmax = 0
    i    = 0
    do
       if(i == 0) then
          read(unit,'(a)',end=1) arg
          nc = num_tokens(arg, " ")-1
       else
          read(unit,*,end=1) l
          if(l > lmax) lmax = l
       end if
       i=i+1
    end do
    1 rewind(unit)
    allocate(pslin(0:lmax,nc),psmat(0:lmax,3,3),beam(0:lmax,3))
    pslin = 0
    do
       read(unit,*,end=2) l, foo(:nc)
       pslin(l,:) = foo(:nc)*2*pi/(l*(l+1))
    end do
    2 close(unit)
    if(info%myid==0) call dmem("read powspec")

    ! Read the beam. Simply the gaussian width FWHM arcmin for now
    read(bfile,*) bsize
    bsize = bsize*pi/180/60/sqrt(8*log(2d0))
    do l = 0, lmax
       beam(l,:) = exp(-0.5*l*(l+1)*bsize**2)
    end do
    if(info%myid==0) call dmem("read beam")

    ! Build the effective power spectrum
    call build_powspec_matrix(pslin, beam, psmat)
    if(info%myid==0) call dmem("build powspec mat")

    ! Read the pixel information from a map
    call read_map(map, pixels, nside, order, mapfile)

    ! The components
    allocate(comps(ncomp))
    select case(ncomp)
       case(1); comps = 1
       case(2); comps = [2,3]
       case(3); comps = [1,2,3]
    end select
    if(info%myid==0) call dmem("read pixels")

    ! Set up the covariance matrix
    if(iargc() >= 8) then
       call getarg(8, icfile)
       open(unit,file=icfile,action="read",status="old",form="unformatted")
       read(unit) n
       read(unit) order2
       read(unit) ncomp2
       close(unit)
       call assert(order2 == order, "Inconsitent ordering in map and cov!")
       call assert(ncomp2 == ncomp, "Requested number of comps different from cov!")
       call assert(n/ncomp2 == size(pixels), "Inconsistent pixels!")
       call read_cov_sc(unit, info, icfile, n, order2, ncomp2, cov, inv)
    else
       n = size(pixels)*ncomp
       call sc_alloc(cov, n, n, info)
       cov%data = 0
    end if
    if(info%myid==0) call dmem("init cov")

    ! And calculate the S matrix
    call build_cmb_cov(psmat, pixels, nside, order, comps, cov, old=1d0)
    if(info%myid==0) call dmem("build smat")

    ! Possibly invert it
    if(invert) then
       call sc_invert(cov)
       if(info%myid==0) call dmem("invert")
    end if

    ! And output
    call write_covmatrix_sc(unit, info, ofile, order, ncomp, cov, invert)
    if(info%myid==0) call dmem("write")

    call angcorr_free(cfun)
    deallocate(pslin, psmat, beam)
    if(info%myid==0) call dmem("done")

  end subroutine


  ! Build a matrix powspec from a linear powspec and a beam,
  ! all of which must be allocated already.
  subroutine build_powspec_matrix(powspec, beam, ps)
    implicit none
    real(dp)     :: powspec(:,:), beam(:,:), ps(:,:,:)
    integer(i4b) :: i, j, k
    do i = 1, size(powspec,1)
       call convert_lin2mat2(powspec(i,:), ps(i,:,:))
    end do
    do j = 1, size(ps,3)
      do i = j, size(ps,2)
         ps(:,i,j) = ps(:,i,j) * beam(:,i) * beam(:,j)
         ps(:,j,i) = ps(:,i,j)
      end do
    end do
  end subroutine

  ! Given a power spectrum in matrix form, build its corresponding
  ! signal covariance matrix as a scalamat. This function is suboptimal
  ! if a single core has the whole matrix, as it calls angcorr_get
  ! individually for each component in a pixel pair. This is necessary
  ! because we cannot in general know that the other components are
  ! held by the same cpu (by making the function significantly more
  ! convoluted, it should be possible to avoid this, though).
  subroutine build_cmb_cov(ps, pixels, nside, order, comps, cov, old)
    implicit none
    real(dp)             :: ps(0:,:,:), v1(3), v2(3), C(3,3), old_
    real(dp), optional   :: old
    integer(i4b)         :: pixels(:), nside, order, comps(:)
    integer(i4b)         :: i1, i2, g1, g2, p1, p2, c1, c2, nc, np
    type(scalamat)       :: cov
    type(angcorr_type)   :: cfun
    old_ = 0; if(present(old)) old_ = old
    nc = size(comps)
    np = size(pixels)
    call angcorr_init(cfun, ps, lmin=2, nspline=100000)
    do i2 = 1, cov%cloc
       g2 = cov%cmap(i2)
       c2 = comps((g2-1)/np+1)
       p2 = modulo(g2-1,np)+1
       v2 = pix2vec(nside, order, pixels(p2))
       do i1 = 1, cov%rloc
          g1 = cov%rmap(i1)
          c1 = comps((g1-1)/np+1)
          p1 = modulo(g1-1,np)+1
          v1 = pix2vec(nside, order, pixels(p1))
          call angcorr_get(cfun, v1, v2, C)
          cov%data(i1,i2) = cov%data(i1, i2)*old_ + C(c1,c2)
       end do
    end do
    call angcorr_free(cfun)
  end subroutine

!-----------------------------------------------------------------------------------------------
! subroutine give_user_info
!-----------------------------------------------------------------------------------------------

  subroutine give_user_info_sc(info)
    implicit none

    type(scinfo),    intent(in) :: info

    if (info%myid == root) then
       write(*,*) 'Usage: mpirun/qrun -n N scalapost command params'
       write(*,*)
       write(*,*) 'scalapost finalmap eqns1.unf eqns2.unf outprefix (maskfile.fits/fast)'
       write(*,*)
       write(*,*) "scalapost solve eqns.unf mask.fits/'nomask' outprefix (onlymap)"
       write(*,*)
       write(*,*) 'scalapost sqrt inv_N.unf outprefix'
       write(*,*)
       write(*,*) 'scalapost lcut map2mask inv_N.unf inmap outprefix lmax (nocov/onlymap)'
       write(*,*) 
       write(*,*) 'scalapost cutoff map2mask inv_N_eigen inmap outprefix  cutoff'
       write(*,*)
       write(*,*) 'scalapost genmap map1.fits map2.fits, icov1.unf icov2.unf map2mask.fits outprefix (fast)'
       write(*,*) 
       write(*,*) 'scalapost invert method mat.unf omat.unf'
       write(*,*) '   method can be eigenvalue, cholesky or LU'
       write(*,*)
       write(*,*) 'scalapost mask inmat.unf inmap2mask.fits inmask.fits outmat.unf outmap2mask.unf'
       write(*,*)
       write(*,*) 'scalapost convert_C_to_F90 matrix_in.unf n ordering pol invcov scale_factor reg_T reg_P matrix_out.unf'
       write(*,*)
       write(*,*) 'scalapost convert_WMAP_to_F90 matrix_in.fits matrix_out.unf'
       write(*,*) 
       write(*,*) 'scalapost fit_template icov.unf map2mask.fits map.fits template.fits'
       write(*,*) '  All must have consistent nside and components.'
       write(*,*)
       write(*,*) 'scalapost covcat icov.unf [icov.unf] [icov.unf] [...] ocov.unf'
       write(*,*) '  Perform a block-diagonal concatenation. Useful for getting IQU from I and QU.'
       write(*,*) '  Uses as much memory as the sum of the input and output matrices.'
       write(*,*)
       write(*,*) 'scalapost eigendecomp cov.unf eig.unf'
       write(*,*) '  Compute eigenvalue decomposition of matrix'
       write(*,*)
       write(*,*) 'scalapost eigencut map2mask inv_N_eigen inmap outmap.hdf outcuts.hdf cutval [cutval [...]]'
       write(*,*) '  Remove eigenmodes lower than the cutvals. Outmap.hdf will contain the result, while outcuts.hdf will contain what was removed.'
       write(*,*)
       write(*,*) 'scalapost eigenswap map2mask inv_N_eigen inmap1 inmap2 outmap1.hdf outmap2.hdf cutval [cutval [...]]'
       write(*,*) '  Swap eigenmodes lower than the cutvals between inmap1 and inmap2. The outmaps will be the swapped versions.'
       write(*,*)
       write(*,*) 'scalapost eigenchi map2mask inv_N_eigen map ofile.txt'
       write(*,*) '  Calculate chisquare per eigenmode for map'
       write(*,*) 'scalapost multiply_by_scalar inmatrix value outmatrix'
       write(*,*) '  Multiply matrix by scalar'
       write(*,*) 'scalapost add_scalar inmatrix value outmatrix'
       write(*,*) '  Add scalar to matrix'
       write(*,*) 'scalapost cov2rms covmat_in map2mask rms_out'
       write(*,*) '  Extract square root of diagonal'
       write(*,*) 'scalapost rms2cov rms_in power covmat_out'
       write(*,*) '  Write diagonal covariance matrix; power refers to covariance matrix'
    end if

    call mpi_finalize(ierr)
    stop

  end subroutine give_user_info_sc

end program scalapost

