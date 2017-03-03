module tod2map_cl_mod
  use alm_tools
  use quiet_fileutils
  use tod2map_utils
  use quiet_healpix_mod
  implicit none

  integer(i4b),                              private :: lmin, lmax, dl, numbin
  integer(i4b),                              private :: nside
  real(dp),     pointer,     dimension(:,:), private :: weights
  real(dp),     allocatable, dimension(:,:), private :: beam
  integer(i4b), allocatable, dimension(:,:), private :: bins
    

contains

  subroutine initialize_tod2map_cl_mod(parfile)
    implicit none

    character(len=*), intent(in) :: parfile

    character(len=512) :: beamfile
    real(dp), pointer, dimension(:,:) :: pixwin

    call get_parameter(0, parfile, 'LMAX_PSEUDO_CLS',      par_int=lmax)
    call get_parameter(0, parfile, 'LMIN_PSEUDO_CLS',      par_int=lmin)
    call get_parameter(0, parfile, 'DL_PSEUDO_CLS',        par_int=dl)
!    call get_parameter(0, parfile, 'BEAMFILE',             par_string=beamfile)
    call get_parameter(0, parfile, 'NSIDE_OUT',            par_int=nside)

    ! Set up beam and pixel window
    weights => get_hpix_ringweights(nside)
    !call read_ringweights(nside, polarization, weights)
!    call read_pixwin(nside, ncomp, pixwin)
!    allocate(beam(0:lmax,ncomp))
!    call read_beam(beamfile, beam)
!    beam = beam * pixwin(0:lmax,1:ncomp)

    ! Set up bins
    call generate_bins(lmin, lmax, dl, bins)
    numbin = size(bins,1)

  end subroutine initialize_tod2map_cl_mod


  subroutine output_pseudo_spectra_and_chisq(info, map2mask, dir, target, split, maps1, maps2, weight, &
       & pte, chisq_tot, chi_tot, chisq_max, nbin_tot)
    implicit none

    type(common_info),                     intent(in)  :: info
    character(len=*),                      intent(in)  :: dir, target, split
    integer(i4b),     dimension(0:),       intent(in)  :: map2mask
    real(dp),         dimension(1:,1:,1:), intent(in)  :: maps1, maps2
    real(dp),         dimension(1:,1:),    intent(in)  :: weight
    real(dp),         dimension(1:,1:),    intent(out) :: pte, chisq_tot, chi_tot, chisq_max
    integer(i4b),     dimension(1:),       intent(out) :: nbin_tot

    integer(i4b)       :: i, j, k, l, m, n, nmap, npix, nspec, bin, unit, ierr, pix, ncomp
    real(dp)           :: t1, t2, chisq
    character(len=512) :: filename
    real(dp),                  dimension(2)       :: z
    real(dp),                  dimension(3)       :: vec
    real(dp),     allocatable, dimension(:,:)     :: map, map_full
    real(dp),     allocatable, dimension(:)       :: mu, sigma
    real(dp),     allocatable, dimension(:,:)     :: cls, cls_tot
    complex(dpc), allocatable, dimension(:,:,:)   :: alms
    complex(dpc), allocatable, dimension(:)       :: cl
   
!    write(*,*) info%myid, 'a'
!    call mpi_barrier(mpi_comm_world, ierr)

    ncomp  = size(maps1,2)
    nmap   = size(maps1,3)
    npix   = 12*nside**2
    nspec  = ncomp*(ncomp+1)/2 
    n      = nspec * numbin
    unit   = getlun()

    call assert(ncomp == 3, "output_pseudo_spectra_and_chisq needs maps1 and maps2 to have 3 components")
    
    ! Find rings to include in spherical harmonics transforms
    z = [1.d0, -1d0]
    do i = 0, npix-1
       if (map2mask(i) > 0) then
          call pix2vec_nest(nside, i, vec)
          z(1) = min(z(1), vec(3))
          z(2) = max(z(2), vec(3))
       end if
    end do

!    write(*,*) info%myid, 'b'
!    call mpi_barrier(mpi_comm_world, ierr)

    ! Compute spectra for all maps
    allocate(alms(ncomp,0:lmax,0:lmax), map(size(maps1,1),ncomp), map_full(0:npix-1,ncomp))
    allocate(cls(n,nmap), cl(0:lmax), cls_tot(n,nmap))
    cls = 0.d0
    do i = 1+info%myid, nmap, info%nproc

       where (maps1(:,:,i) /= hpx_dbadval .and. maps2(:,:,i) /= hpx_dbadval)
          map = weight * 0.5d0 * (maps1(:,:,i) - maps2(:,:,i))
       elsewhere
          map = 0.d0
       end where

       ! Expand onto full Healpix sky
       map_full = 0.d0
       do pix = 0, npix-1
          if (map2mask(pix) > 0) then
             map_full(pix,:) = map(map2mask(pix),:)
          end if
       end do

       do j = 1, ncomp
          call convert_nest2ring(nside, map_full(:,j))
       end do

       ! Compute alms
       if (info%ncomp == 1 .and. info%comps(1) == T) then
          call map2alm(nside, lmax, lmax, map_full(:,1), alms, z, weights)
       else
          call map2alm(nside, lmax, lmax, map_full, alms, z, weights)
       end if

       ! Compute binned spectra
       bin = 1
       do j = 1, ncomp
!          do k = j, j     ! removing EB, keping EE and BB
          do k = j, ncomp ! Keeping all components 
             do l = 0, lmax
                cl(l) = alms(j,l,0)*alms(k,l,0)
                do m = 1, l
                   cl(l) = cl(l) + alms(j,l,m)*conjg(alms(k,l,m)) + conjg(alms(j,l,m)) * alms(k,l,m)
                end do
                cl(l) = cl(l) / (2*l+1) * (l*(l+1)) / (2*pi)
             end do
             do l = 1, numbin
                cls(bin,i) = sum(cl(bins(l,1):bins(l,2))) / (bins(l,2)-bins(l,1)+1)
                bin = bin+1
             end do
          end do
       end do

    end do
    call mpi_reduce(cls, cls_tot, size(cls), mpi_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

!    write(*,*) info%myid, 'c'
!    call mpi_barrier(mpi_comm_world, ierr)

    if (info%myid == 0) then
       cls = cls_tot

       ! Compute mean spectrum and errors
       allocate(mu(n), sigma(n))
       mu    = 0.d0
       sigma = 0.d0
       do i = 1, n
          !mu(i)    = mean(cls(i,:)) 
          !sigma(i) = sqrt(variance(cls(i,:)))
          mu(i)    = sum(cls(i,:)) / real(nmap,dp)
          sigma(i) = sqrt(sum((cls(i,:)-mu(i))**2) / real(nmap,dp))
       end do

       nbin_tot = 0
       do j = 1, n
          if (sigma(j) > 0.d0) then
             nbin_tot(1) = nbin_tot(1) + 1
             if ((j-1)/numbin == 3) nbin_tot(2) = nbin_tot(2) + 1 ! EE
             if ((j-1)/numbin == 5) nbin_tot(3) = nbin_tot(3) + 1 ! BB
          end if
       end do
       
       ! Output mean and standard deviation
       filename = trim(dir) // '/cls/mu_sigma_' // trim(target) // '_' // trim(split) // '.dat'
       call mkdirs(filename, .true.)
       open(unit,file=trim(filename), recl=4096)
       bin = 1
       do j = 1, nspec
          do i = 1, numbin
             if (sigma(bin) /= 0.d0) write(unit,*) 0.5*(bins(i,1)+bins(i,2)), mu(bin), sigma(bin)
             bin = bin+1
          end do
          write(unit,*)
       end do
       close(unit)

       ! Compute normalized spectra
       do i = 1, nmap
          do j = 1, n
             if (sigma(j) > 0.d0) then
                cls(j,i) = (cls(j,i) - mu(j)) / sigma(j)
             else
                cls(j,i) = 0.d0
             end if
          end do
       end do

       ! Output power spectrum
       filename = trim(dir) // '/cls/norm_cl_' // trim(target) // '_' // trim(split) // '.dat'
       call mkdirs(filename, .true.)
       open(unit,file=trim(filename), recl=4096)
       bin = 1
       do j = 1, nspec
          do i = 1, numbin
             if (sigma(bin) /= 0.d0) write(unit,*) 0.5*(bins(i,1)+bins(i,2)), cls(bin,1)
             bin = bin+1
          end do
          write(unit,*)
       end do
       close(unit)
       
       ! Output chi's
       filename = trim(dir) // '/cls/chi_' // trim(target) // '_' // trim(split) // '.dat'
       open(unit,file=trim(filename))
       write(unit,*) '# Nmap = ', nmap, ', numbin = ', nbin_tot
       chi_tot = 0.d0
       do i = 1, nmap
          do j = 1, n
             if (sigma(j) > 0.d0) then
                !write(unit,*) j,cls(j,i)
                write(unit,*) cls(j,i)
                chi_tot(i,1) = chi_tot(i,1) + cls(j,i)
                if ((j-1)/numbin == 3) chi_tot(i,2) = chi_tot(i,2) + cls(j,i) ! EE
                if ((j-1)/numbin == 5) chi_tot(i,3) = chi_tot(i,3) + cls(j,i) ! BB
             end if
          end do
          write(unit,*)
       end do
       close(unit)
       
       ! Output chi-squares
       cls       = cls**2
       chisq_tot = 0.d0
       filename = trim(dir) // '/cls/chisq_' // trim(target) // '_' // trim(split) // '.dat'
       open(unit,file=trim(filename))
       write(unit,*) '# Nmap = ', nmap, ', numbin = ', nbin_tot
       do i = 1, nmap
          do j = 1, n
             if (sigma(j) > 0.d0) then
                write(unit,*) cls(j,i)
                chisq_tot(i,1) = chisq_tot(i,1) + cls(j,i)
                if ((j-1)/numbin == 3) chisq_tot(i,2) = chisq_tot(i,2) + cls(j,i) ! EE
                if ((j-1)/numbin == 5) chisq_tot(i,3) = chisq_tot(i,3) + cls(j,i) ! BB
             end if
          end do
          write(unit,*)
       end do
       close(unit)

       ! Return maximum chisquare 
       chisq_max = 0.d0
       do i = 1, nmap
          do j = 1, n
             if (sigma(j) > 0.d0) then
                chisq_max(i,1) = max(chisq_max(i,1), cls(j,i))
                if ((j-1)/numbin == 3) chisq_max(i,2) = max(chisq_max(i,2),cls(j,i)) ! EE
                if ((j-1)/numbin == 5) chisq_max(i,3) = max(chisq_max(i,3),cls(j,i)) ! BB
             end if
          end do
       end do
       
       ! Return PTEs
       pte = 0.d0
       do j = 1, nmap
          m     = 0
          do i = 1, nmap
             if (i /= j .and. sum(cls(:,i)) > sum(cls(:,j))) m = m+1
          end do
          pte(j,1) = real(m,dp) / real(nmap-1,dp)
       end do

       if (nspec > 1) then
          do j = 1, nmap
             m     = 0
             do i = 1, nmap
                if (i /= j .and. sum(cls(3*numbin+1:4*numbin,i)) > sum(cls(3*numbin+1:4*numbin,j))) m = m+1
             end do
             pte(j,2) = real(m,dp) / real(nmap-1,dp) ! EE
          end do

          do j = 1, nmap
             m     = 0
             do i = 1, nmap
                if (i /= j .and. sum(cls(5*numbin+1:6*numbin,i)) > sum(cls(5*numbin+1:6*numbin,j))) m = m+1
             end do
             pte(j,3) = real(m,dp) / real(nmap-1,dp) ! BB
          end do
       end if

       filename = trim(dir) // '/cls/pte_' // trim(target) // '_' // trim(split) // '.dat'
       open(unit,file=trim(filename))
       write(unit,*) '# Nmap = ', nmap
       do i = 1, nmap
          write(unit,*) pte(i,:)
       end do
       close(unit)

       deallocate(mu, sigma)

    end if
       
    deallocate(alms,cls,cls_tot,map_full)

!    write(*,*) info%myid, 'd'
!    call mpi_barrier(mpi_comm_world, ierr)

  end subroutine output_pseudo_spectra_and_chisq

  subroutine generate_bins(lmin, lmax, dl, bins)
    implicit none

    integer(i4b),                              intent(in) :: lmin, lmax, dl
    integer(i4b), allocatable, dimension(:,:)             :: bins

    integer(i4b) :: numbin, i

    numbin = (lmax-lmin+1) / dl
    if (lmax-lmin+1 > numbin*dl) numbin = numbin+1
    allocate(bins(numbin,2))
    do i = 1, numbin
       bins(i,1) = lmin + (i-1)*dl
       bins(i,2) = min(lmin + i*dl - 1, lmax)
    end do

  end subroutine generate_bins


end module tod2map_cl_mod
