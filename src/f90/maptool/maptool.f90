program maptool
  use healpix_types
  use pix_tools
  use fitstools
  use quiet_fileutils
  use math_tools
  use alm_tools
  use quiet_postutils
  use quiet_utils
  use quiet_mapfile_mod
  use powell_mod
  use rngmod
  use quiet_mpi_mod
  use quiet_task_mod
  use quiet_pixspace_mod
  use comap_map_mod

  implicit none

!  include "mpif.h"

  integer(i4b)       :: iargc
  integer(i4b)       :: unit, myid, numprocs, ierr, root
  character(len=30)  :: kommando
  type(planck_rng)   :: rng_handle, rng_handle_cmb, rng_handle_noise

  ! For scatter_sindex
  integer(i4b)                            :: nside, ordering
  real(dp)                                :: cvec(3) 
  real(dp),     allocatable, dimension(:) :: mapdata
  integer(i4b), allocatable, dimension(:) :: pixlist

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

  root = 0
  unit        = 42+myid

  if (iargc() == 0) then
     call give_user_info
  else if (myid == root) then
     write(*,*) '-------------------- QUIET maptool --------------------'
     write(*,*)
  end if

  ! Get name of main kommando
  call getarg(1,kommando)

  if (kommando == 'sindex') then
     if (myid==0) write(*,*) 'Find spectral index from two maps'
     call spectral_index2(unit)
  else if (kommando == 'scatter') then
     if (myid==0) write(*,*) 'Find scatterplot spectral index from two smoothed antenna maps'
     call scatter_sindex(unit)
  else if (kommando == 'leakage') then
     if (myid==0) write(*,*) 'Estimate I->QU leakage parameters'
     call leakage(unit)
  else if (kommando == 'transfer') then
     if (myid==0) write(*,*) 'Estimate transfer function between two maps'
     call transfer(unit)
  else if (kommando == 'smooth') then
     if (myid==0) write(*,*) 'Smooth transfer function'
     call smooth_transfer(unit)
  else if (kommando == 'sims') then
     if (myid==0) write(*,*) 'Make simulations from powspec'
     call make_sims(unit)
  else if (kommando == 'clusim') then
     if (myid==0) write(*,*) 'Make cluster simulations'
     call clustersims(unit)
  else if (kommando == 'beam') then
     if (myid==0) write(*,*) 'Estimate beam from powspecs'
     call estimate_beam(unit)
  else if (kommando == 'amp') then
     if (myid==0) write(*,*) 'Estimate overall amplitude from powspecs'
     call estimate_amp(unit)
  else if (kommando == 'amp2') then
     if (myid==0) write(*,*) 'Estimate overall amplitude from unbinned powspecs (like growth rate..)'
     call estimate_amp2(unit)
  else if (kommando == 'amp3') then
     if (myid==0) write(*,*) 'Estimate which powspec fits the data'
     call estimate_amp3(unit)
  else if (kommando == 'psamp') then
     if (myid==0) write(*,*) 'Estimate unresolved point source amplitudes from multiple frequency powspecs'
     call estimate_psamp(unit)
  else if (kommando == 'bin') then
     if (myid==0) write(*,*) 'Binning file'
     call binning(unit)
  else if (kommando == 'horns') then
     if (myid==0) write(*,*) 'Plot focalplane'
     call plot_focalplane(unit)
  else if (kommando == 'merge') then
     if (myid==0) write(*,*) 'Merge components to one fits file'
     call merge_fits(unit)
  else if (kommando == 'extract') then
     if (myid==0) write(*,*) 'Extract components from fits file'
     call extract_fits(unit)
  else if (kommando == 'bin') then
     if (myid==0) write(*,*) 'bin'
     call bin_spec(unit)
  else if (kommando == 'process_ASCII_table') then
     if (myid==0) write(*,*) 'Print mean and standard deviation of ASCII file columns; first row is data'
     call process_ASCII_table(unit)
  else if (kommando == 'mdmaps') then
     if (myid==0) write(*,*) 'Output (mono-)/dipole maps from list of values'
     call output_mdmaps(unit)
  else if (kommando == 'beta') then
     if (myid==0) write(*,*) 'synthesize maps based on signal-to-noise'
     call make_weighted_map(unit)
  else if (kommando == 'gain') then
     if (myid==0) write(*,*) 'Esitimate gain between map of source and known strength'
     call estimate_gain_from_map(unit)
  else
     call give_user_info
  end if

  ! And exit
  call mpi_finalize(ierr)
  if (myid == root) then 
     write(*,*)
     write(*,*) '-------------------- QUIET maptool finished -----------'
  end if
  
contains

  !-----------------------------------------------------------------------------------------------
  ! subroutine estimate beam
  !-----------------------------------------------------------------------------------------------

  subroutine estimate_beam(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, datafile, outprefix, outfile
    integer(i4b)                :: i, j, l, bin, lmax, numstep, numbin, numcol, specol
    real(dp)                    :: tall, expt, binpow, omega, arr(1), var, mu, sigma
    real(dp)                    :: amp, amax, amin, theta, tmin, tmax

    real(dp),       allocatable, dimension(:,:)   :: power, data
    real(dp),       allocatable, dimension(:)   :: like, amps

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'beam takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, powspecfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
    end if

    lmax  = 1024.d0

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    allocate(power(0:lmax,1:4))
    call read_powspec(powspecfile, power)
!    open(13, file='power.dat', recl=1024)
    do l = 0, lmax
       power(l,:) = power(l,:) *l*(l+1)/(2*pi)
!       write(13,*) l, power(l,:)
    end do
!    close(13)


    ! Read data
    numcol = 5 !15
    specol = 4 !10
    numbin = 16 !19 !20
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
    allocate(data(numbin,numcol))
    open(13, file=trim(datafile))
    do i = 1, numbin
       read(13,*) data(i,:)
    end do
    close(13)

!data(:,4:15)=data(:,4:15)/0.88**2 !ml

    ! Calculate likelihood
    tmin = -30.d0
    tmax = 10.d0
!    amin = 0.5d0
!    amax = 2.0d0
    numstep = 300
    tall = (pi/(180.d0*60.d0*log(2.d0)*sqrt(8.d0)))**2
!    allocate(like(numstep, numstep))
    allocate(like(numstep))
    like = 0.d0
    do j = 1, numstep
!       amp = amin + (real(j,dp)-1)*(amax-amin)/(numstep-1)
!       do i = 1, numstep
          theta = tmin + (real(j,dp)-1)*(tmax-tmin)/(numstep-1)
          if (theta < 0.d0) then
             expt  = -theta**2*tall
          else
             expt  = theta**2*tall
          end if
amp = 1.d0
          do bin = 1, numbin
             binpow = 0.d0
             do l = data(bin,2), data(bin,3)
                binpow = binpow + power(l,2)*amp*exp(-l*(l+1)*expt)
!write(*,*) l, power(l,2)
             end do
             binpow = binpow/(data(bin,3)-data(bin,2)+1)
             like(j) = like(j) + ((data(bin,specol)-binpow)/data(bin,specol+1))**2
!write(*,*) bin, data(bin,specol), data(bin,specol+1), binpow             
!          end do
       end do
    end do
 
    ! Write data
    arr   = minloc(like)
    theta = tmin + (real(arr(1),dp)-1)*(tmax-tmin)/(numstep-1)
!    amp   = amin + (real(arr(1),dp)-1)*(amax-amin)/(numstep-1)
    write(*,*) 'Min like', minval(like), 'for theta =', theta, ',amp =',amp
!    write(*,*) 'Min like', minval(like), ',for amp =',amp
    outfile = trim(outprefix) // '_chisq.dat'
    open(13, file=trim(outfile))
!    write(13,*) numstep, numstep
    do j = 1, numstep
       theta = tmin + (real(j,dp)-1)*(tmax-tmin)/(numstep-1)
!       do i = 1, numstep
!          theta = tmin + (real(i,dp)-1)*(tmax-tmin)/(numstep-1)
          write(13,*) theta, like(j) !exp(-like(i)/2.d0)
 !      end do
    end do
    close(13)
    write(*,*) '* Written to file = ', trim(outfile)

    allocate(amps(numstep))
    do i = 1, numstep
       amps(i) = tmin + (i-1) * (tmax-tmin)/(numstep-1)
    end do
    like = exp(-0.5*(like-minval(like)))
    like = like / (sum(like)*(tmax-tmin)/(numstep-1))
    mu = sum(like*amps) * (tmax-tmin)/(numstep-1)
    sigma = sqrt(sum(like*(amps-sum(like*amps) * (tmax-tmin)/(numstep-1))**2) * (tmax-tmin)/(numstep-1))
    write(*,*) 'Mean   = ', mu, ' +/- ', sigma
    write(*,*) 'Sigma from one = ', (mu-1.d0)/sigma

  end subroutine estimate_beam

  !-----------------------------------------------------------------------------------------------
  ! subroutine estimate amplitude
  !-----------------------------------------------------------------------------------------------

  subroutine estimate_amp(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, datafile, outprefix, outfile
    integer(i4b)                :: i, j, l, bin, lmax, numstep, numbin, numcol, specol
    real(dp)                    :: tall, expt, binpow, omega, arr(1), var, mu, sigma
    real(dp)                    :: amp, amax, amin, theta, tmin, tmax

    real(dp),       allocatable, dimension(:,:)   :: power, data
    real(dp),       allocatable, dimension(:)   :: like, amps

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'amp takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, powspecfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
    end if

    lmax  = 1100.d0

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    allocate(power(0:lmax,1:4))
    call read_powspec(powspecfile, power)
!    open(13, file='power.dat', recl=1024)
    do l = 0, lmax
       power(l,:) = power(l,:) *l*(l+1)/(2*pi)
!       write(13,*) l, power(l,:)
    end do
!    close(13)


    ! Read data
    numcol = 5 !15
    specol = 4 !10
    numbin = 16 !19 !20
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
    allocate(data(numbin,numcol))
    open(13, file=trim(datafile))
    do i = 1, numbin
       read(13,*) data(i,:)
    end do
    close(13)

!data(:,4:15)=data(:,4:15)/0.88**2 !ml

    ! Calculate likelihood
!    tmin = 0.d0
!    tmax = 0.d0
    amin = 0.5d0
    amax = 2.0d0
    numstep = 300
!    tall = (pi/(180.d0*60.d0*log(2.d0)*sqrt(8.d0)))**2
!    allocate(like(numstep, numstep))
    allocate(like(numstep))
    like = 0.d0
    do j = 1, numstep
       amp = amin + (real(j,dp)-1)*(amax-amin)/(numstep-1)
!       do i = 1, numstep
!          theta = tmin + (real(i,dp)-1)*(tmax-tmin)/(numstep-1)
!          if (theta < 0.d0) then
!             expt  = -theta**2*tall
!          else
!             expt  = theta**2*tall
!          end if
!amp = 1.d0
          do bin = 1, numbin
             binpow = 0.d0
             do l = data(bin,2), data(bin,3)
                binpow = binpow + power(l,2)*amp!*exp(-l*(l+1)*expt)
!write(*,*) l, power(l,2)
             end do
             binpow = binpow/(data(bin,3)-data(bin,2)+1)
             like(j) = like(j) + ((data(bin,specol)-binpow)/data(bin,specol+1))**2
!write(*,*) bin, data(bin,specol), data(bin,specol+1), binpow             
!          end do
       end do
    end do
    
    ! Normalizing
!    like = exp(-like/2.d0)
!    mu = 0.d0
!    do j = 1, numstep
!       do i = 1, numstep
!          mu = mu + like(i,j)*(tmax-tmin+1)/numstep *(amax-amin+1)/numstep 
!       end do
!    end do
!    like = like/mu
!    mu = 0.d0
!    do j = 1, numstep
!       amp = amin + (real(j,dp)-1)*(amax-amin)/(numstep-1)
!       do i = 1, numstep
!          theta = (tmin + (real(i,dp)-1)*(tmax-tmin)/(numstep-1))
!          mu = mu + like(i,j)*theta*(tmax-tmin+1)/numstep
!       end do
!    end do
!    var = 0.d0
!    do j = 1, numstep
!       do i = 1, numstep
!          theta = (tmin + (real(i,dp)-1)*(tmax-tmin)/(numstep-1))
!          var = var + like(i,j)*(tmax-tmin+1)/numstep *(theta-mu)**2
!       end do
!    end do
!    write(*,*) mu, sqrt(var), 'mu, var'

    ! Write data
    arr   = minloc(like)
!    theta = tmin + (real(arr(1),dp)-1)*(tmax-tmin)/(numstep-1)
    amp   = amin + (real(arr(1),dp)-1)*(amax-amin)/(numstep-1)
!    write(*,*) 'Min like', minval(like), 'for theta =', theta, ',amp =',amp
    write(*,*) 'Min like', minval(like), ',for amp =',amp
    outfile = trim(outprefix) // '_chisq.dat'
    open(13, file=trim(outfile))
!    write(13,*) numstep, numstep
    do j = 1, numstep
       amp = amin + (real(j,dp)-1)*(amax-amin)/(numstep-1)
!       do i = 1, numstep
!          theta = tmin + (real(i,dp)-1)*(tmax-tmin)/(numstep-1)
          write(13,*) amp, like(j) !exp(-like(i)/2.d0)
 !      end do
    end do
    close(13)
    write(*,*) '* Written to file = ', trim(outfile)

    allocate(amps(numstep))
    do i = 1, numstep
       amps(i) = amin + (i-1) * (amax-amin)/(numstep-1)
    end do
    like = exp(-0.5*(like-minval(like)))
    like = like / (sum(like)*(amax-amin)/(numstep-1))
    mu = sum(like*amps) * (amax-amin)/(numstep-1)
    sigma = sqrt(sum(like*(amps-sum(like*amps) * (amax-amin)/(numstep-1))**2) * (amax-amin)/(numstep-1))
    write(*,*) 'Mean   = ', mu, ' +/- ', sigma
    write(*,*) 'Sigma from one = ', (mu-1.d0)/sigma

  end subroutine estimate_amp
  !-----------------------------------------------------------------------------------------------
  ! subroutine estimate amplitude
  !-----------------------------------------------------------------------------------------------

  subroutine bin_spec(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, datafile, outprefix, outfile
    integer(i4b)                :: i, j, l, bin, lmax, numstep, numbin, numcol, specol,binsize
    real(dp)                    :: tall, expt, binpow, omega, arr(1), var, mu, sigma
    real(dp)                    :: amp, amax, amin, theta, tmin, tmax

    real(dp),       allocatable, dimension(:,:)   :: power, data, bindata
    real(dp),       allocatable, dimension(:)   :: like, amps

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'amp takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, powspecfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
    end if

    lmax  = 1100.d0

    ! Read data
    numcol = 2 !15
    numbin = 492 !19 !20
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
    allocate(data(numbin,numcol))
    open(13, file=trim(datafile))
    do i = 1, numbin
       read(13,*) data(i,:)
!       write(*,*) data(i,1), data(i,2)  
    end do
    close(13)

    binsize = 5
    allocate(bindata(9,numcol))
    bindata=0.d0
    do i = 1,9
       write(*,*) (i-1)*binsize+1,i*binsize, data((i-1)*binsize+1,1), data(i*binsize,1)
!       do j = (i-1)*binsize+1,i*binsize
       bindata(i,1) = sum(data((i-1)*binsize+1:i*binsize,1))/real(binsize,dp)
       bindata(i,2) = sum(data((i-1)*binsize+1:i*binsize,2))/real(binsize,dp)
!          bindata(i,1) = bindata(i,1) + data(j,1)
!          bindata(i,2) = bindata(i,2) + data(j,2)
!          write(*,*) j, bindata(i,1), bindata(i,2)
!       end do
    end do


    outfile = trim(outprefix)
    open(13, file=trim(outfile))
    do i = 1, 9
       write(13,*) bindata(i,:)
    end do
    close(13)
    write(*,*) '* Written to file = ', trim(outfile)


    tall=0.d0
    do i = 21,55
       tall = tall +i
    end do
    write(*,*) tall, tall/binsize


  end subroutine bin_spec

  !-----------------------------------------------------------------------------------------------
  ! subroutine estimate amplitude unbinned
  !-----------------------------------------------------------------------------------------------

  subroutine estimate_amp2(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, datafile, outprefix, outfile
    integer(i4b)                :: i, j, l, bin, lmax, numstep, numbin, numcol, specol, numpow
    real(dp)                    :: tall, expt, binpow, omega, arr(1), var, mu, sigma
    real(dp)                    :: amp, amax, amin, theta, tmin, tmax

    real(dp),       allocatable, dimension(:,:)   :: power, data
    real(dp),       allocatable, dimension(:)   :: like, amps

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'amp takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, powspecfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
    end if

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    open(13, file=trim(powspecfile))
    read(13,*) numpow
    write(*,*) numpow,'= numpow'
    allocate(power(numpow,2))
    do i = 1, numpow
       read(13,*) power(i,:)
    end do
    close(13)

    ! Read data
    numcol = 3
    specol = 2
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
    open(13, file=trim(datafile))
    read(13,*) numbin
    call assert(numpow==numbin, "mismatching number of datapoints in data and spectrum")
    allocate(data(numbin,numcol))
    open(13, file=trim(datafile))
    do i = 1, numbin
       read(13,*) data(i,:)
    end do
    close(13)

    ! Calculate likelihood
    amin = 0.5d0
    amax = 1.5d0
    numstep = 300
    allocate(like(numstep))
    like = 0.d0
    do j = 1, numstep
       amp = amin + (real(j,dp)-1)*(amax-amin)/(numstep-1)
       do bin = 1, numbin
          binpow = power(bin,2)*amp
          like(j) = like(j) + ((data(bin,specol)-binpow)/data(bin,specol+1))**2
       end do
    end do

    ! Write data
    arr   = minloc(like)
    amp   = amin + (real(arr(1),dp)-1)*(amax-amin)/(numstep-1)
    write(*,*) 'Min chisq =', minval(like), ',for amp =',amp
    outfile = trim(outprefix) // '_chisq.dat'
    open(13, file=trim(outfile))
    do j = 1, numstep
       amp = amin + (real(j,dp)-1)*(amax-amin)/(numstep-1)
       write(13,*) amp, like(j) !exp(-like(i)/2.d0)
    end do
    close(13)
    write(*,*) '* Written to file = ', trim(outfile)

    allocate(amps(numstep))
    do i = 1, numstep
       amps(i) = amin + (i-1) * (amax-amin)/(numstep-1)
    end do
    like = exp(-0.5*(like-minval(like)))
    like = like / (sum(like)*(amax-amin)/(numstep-1))
    mu = sum(like*amps) * (amax-amin)/(numstep-1)
    sigma = sqrt(sum(like*(amps-sum(like*amps) * (amax-amin)/(numstep-1))**2) * (amax-amin)/(numstep-1))
    write(*,*) 'Mean   = ', mu, ' +/- ', sigma
    write(*,*) 'Sigma from one = ', (mu-1.d0)/sigma

  end subroutine estimate_amp2

  !-----------------------------------------------------------------------------------------------
  ! subroutine estimate which powspec fits the data
  !-----------------------------------------------------------------------------------------------

  subroutine estimate_amp3(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, datafile, outprefix, outfile, lookupfile
    integer(i4b)                :: i, j, k, bin, lmax, numstep, numbin, numcol, specol
    integer(i4b)                :: numpoints, numspec, begcov, endcov
    integer(i4b)                :: numlook, startlook, dlook
    real(dp)                    :: tall, expt, binpow, omega, arr(1), var, mu, sigma, min_pos
    real(dp)                    :: wval, hm, lcdm_pos, startx, dx, xlcdm, cov_test(3,3), W(3)
    real(dp)                    :: startz, dz, sigma8, wval2, integral
    logical(lgt)                :: lookup
 
    real(dp),       allocatable, dimension(:,:) :: power, data, pow2der, pow, cov, lookuptable
    real(dp),       allocatable, dimension(:)   :: like, amps, zval, chisq, look2der

!    cov_test(1,1:3) = [0.00305163018336723, 0.00139162336885331, -0.000684684592834213]
!    cov_test(2,1:3) = [  0.00139162336885331,   0.00150223937087767,   0.00160252265527017]
!    cov_test(3,1:3) = [-0.000684684592834213,   0.00160252265527017,   0.00437929812791692]
!    call get_eigenvalues(cov_test, W)
!    write(*,*) W
!    stop

    ! Get parameters
    if (iargc() /= 4 .and. iargc() /= 5) then
       write(*,*) 'amp takes 3(4) parameters'
       call give_user_info
    else 
       call getarg(2, powspecfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
       if (iargc() == 5)  then 
          call getarg(5, lookupfile)
          lookup = .true.
       else
          lookup = .false.
       end if
    end if

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    open(13, file=trim(powspecfile))
    read(13,*) numpoints, numspec, startx, dx, xlcdm, startz, dz
    write(*,*) numspec,'= numspec', numpoints,'= length of powspec'
    write(*,*) 'x running from', real(startx,sp), 'to',real(startx+(numspec-1)*dx,sp), 'in steps of', real(dx,sp)
    write(*,*) 'z running from', real(startz,sp), 'to',real(startz+(numpoints-1)*dz,sp), 'in steps of', real(dz,sp)
    allocate(power(numpoints,numspec))
    do i = 1, numpoints
       read(13,*) power(i,:)
    end do
    close(13)

    ! Read lookup table
    if (lookup) then
       if (myid==root) write(*,*) 'Reading from file ', trim(lookupfile)
       open(13, file=trim(lookupfile))
       read(13,*) numlook  !, startlook, dlook
       write(*,*) numlook,'= number of entries in lookup table'
   !    write(*,*) 'lookup running from', real(startlook,sp), 'to',real(startlook+(numlook-1)*dlook,sp), 'in steps of', real(dlook,sp)
       allocate(lookuptable(numlook,2))
       do i = 1, numlook
          read(13,*) lookuptable(i,:)
       end do
       close(13)
    end if

    ! Read data points
    numcol = 3
    specol = 2
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
    open(13, file=trim(datafile))
    read(13,*) numbin, begcov, endcov
    write(*,*) numbin,'= number of datapoints'
    allocate(data(numbin,numcol))
    open(13, file=trim(datafile))
    do i = 1, numbin
       read(13,*) data(i,:)
    end do
    close(13)

    ! Spline powspecs
    allocate(zval(numpoints))
    allocate(pow2der(numpoints,numspec))
    allocate(pow(numbin,numspec))
    do i = 1, numpoints
       zval(i) = startz +real(i-1,dp)*dz
    end do
    do j = 1, numspec
       call spline(zval, power(:,j), 1.d30, 1.d30, pow2der(:,j))
       do i = 1, numbin
          pow(i,j) = splint(zval, power(:,j), pow2der(:,j), data(i,1))
       end do
    end do

    ! Spline lookup
    if (lookup) then
       allocate(look2der(numlook))
       call spline(lookuptable(:,1), lookuptable(:,2), 1.d30, 1.d30, look2der)
    end if

    ! Build covariance matrix
    allocate(cov(numbin,numbin))
!    cov = 0.d0
!    do bin = 1, numbin
!       cov(bin, bin) = data(bin,specol+1)**2
!    end do
!    do j = begcov, endcov
!       do i = j, endcov
!          if (i-j==1) then
!             cov(i,j) = 0.9 *sqrt(data(i,specol+1)*data(j,specol+1))
!             cov(j,i) = cov(i,j)
!          else if (i-j==2) then
!             cov(i,j) = 0.5 *sqrt(data(i,specol+1)*data(j,specol+1))
!             cov(j,i) = cov(i,j)
!          else if (i-j==3) then
!             cov(i,j) = 0.15 *sqrt(data(i,specol+1)*data(j,specol+1))
!             cov(j,i) = cov(i,j)
!          end if
!       end do
!    end do
!write(*,*) ''
!write(*,*) cov(11,11:13)
!write(*,*) cov(12,11:13)
!write(*,*) cov(13,11:13)
    cov = 0.d0
    do bin = 6, numbin
       cov(bin, bin) = data(bin,specol+1)**2
    end do 
    cov(1,1:2) = [0.00305163018336723, -0.000684684592834213] *1.d0
!    cov(2,1:3) = [  0.00139162336885331,   0.00150223937087767,   0.00160252265527017] *1.d0
    cov(2,1:2) = [-0.000684684592834213,   0.00437929812791692] *1.d0
    cov(3,3:5) = [6.400, 2.570, 0.000] *0.001d0
    cov(4,3:5) = [2.570, 3.969, 2.540]  *0.001d0
    cov(5,3:5) = [0.000, 2.540, 5.184]  *0.001d0
    call invert_matrix(cov)

    ! Calculate likelihood
    allocate(like(numspec))
    like = 0.d0
    do k = 1, numspec
       do j = 1, numbin
          do i = 1, numbin
             like(k) = like(k) + (data(i,specol)-pow(i,k))*cov(i,j)*(data(j,specol)-pow(j,k))
          end do
       end do
    end do

    ! Write chisq data
    arr   = minloc(like)
    wval  = startx + real(arr(1)-1,dp)*dx
    write(*,*) 'Min. chisq =', minval(like), 'for x =',wval
    outfile = trim(outprefix) // '_chisq.dat'
    open(13, file=trim(outfile))
    do j = 1, numspec
       wval2  = startx + real(j-1,dp)*dx
       write(13,*) wval2, like(j)
       if (wval2 == xlcdm) lcdm_pos = j
    end do
    close(13)
    write(*,*) 'Lcdm chisq =', like(lcdm_pos), 'for x =', xlcdm
!    write(*,*) lcdm_pos, '= lcdm_pos in powspec file'
!    write(*,*) '* Written to file = ', trim(outfile)

min_pos=arr(1)
allocate(chisq(numbin))
chisq = 0.d0
do i = 1, 2
   do j = 1, 2
      chisq(1) = chisq(1) + (data(i,specol)-pow(i,min_pos))*cov(i,j)*(data(j,specol)-pow(j,min_pos))
   end do
end do
do i = 3, 5
   do j = 3, 5
      chisq(3) = chisq(3) + (data(i,specol)-pow(i,min_pos))*cov(i,j)*(data(j,specol)-pow(j,min_pos))
   end do
end do
do i = 6, numbin
      chisq(i) = chisq(i) + (data(i,specol)-pow(i,min_pos))*cov(i,i)*(data(i,specol)-pow(i,min_pos))
end do
do i = 1, numbin
write(*,*) i, chisq(i)
end do
write(*,*) 'total chisq', sum(chisq)

    outfile = trim(outprefix) // '_bestfit_spec.dat'
    open(13, file=trim(outfile))
    do i = 1, numpoints
       write(13,*) zval(i), power(i,min_pos)
    end do
    close(13)
    outfile = trim(outprefix) // '_lcdm_spec.dat'
    open(13, file=trim(outfile))
    do i = 1, numpoints
       write(13,*) zval(i), power(i,lcdm_pos)
    end do
    close(13)

    ! Write likelihood sigma
    like = exp(-0.5*(like-minval(like)))
    allocate(amps(numspec))
    do i = 1, numspec
       amps(i) = startx + real(i-1,dp)*dx
    end do
    outfile = trim(outprefix) // '_like.dat'
    open(13, file=trim(outfile))
    do j = 1, numspec
       write(13,*) amps(j), like(j)
    end do
    close(13)
!    write(*,*) '* Written to file = ', trim(outfile)

    ! Calculate sigma
    like = like / (sum(like) * dx)
    mu = sum(like*amps) * dx
    sigma = sqrt(sum(like*(amps-mu)**2) * dx)
    write(*,*) 'Mean   = ', mu, ' +/- ', sigma
    write(*,*) 'Sigma from lcdm = ', (mu-xlcdm)/sigma
    write(*,*) 'Min   = ', wval, ' +/- ', sigma
    if (lookup) then
       integral = 0.d0
       do i = 1, numspec-1
          integral = integral + (like(i+1)+like(i))*(lookuptable(i+1,2)-lookuptable(i,2))/2.d0
       end do
write(*,*) 'norm', integral

       like = like/integral
       mu = 0.d0
       do i = 1, numspec-1
          mu = mu + (like(i+1)*lookuptable(i+1,2)+like(i)*lookuptable(i,2))*(lookuptable(i+1,2)-lookuptable(i,2))/2.d0
       end do
       integral = 0.d0
       do i = 1, numspec-1
          integral = integral + (like(i+1)*(lookuptable(i+1,2)-mu)**2+like(i)*&
               & (lookuptable(i,2)-mu)**2)*(lookuptable(i+1,2)-lookuptable(i,2))/2.d0
       end do
       write(*,*) 'Sigma8 mean   = ', mu, ' +/- ', sqrt(integral)
       sigma8 = splint(lookuptable(:,1), lookuptable(:,2), look2der, wval)
       write(*,*) 'Sigma8 min  = ', sigma8, ' +/- ',  sqrt(integral)
       write(*,*) 'sigma upper = ', splint(lookuptable(:,1), lookuptable(:,2), look2der, wval+sigma)-sigma8
       write(*,*) 'sigma lower = ', sigma8 - splint(lookuptable(:,1), lookuptable(:,2), look2der, wval-sigma)
    end if

  end subroutine estimate_amp3

  !-----------------------------------------------------------------------------------------------
  ! subroutine estimate unresolved point source amplitudes
  !-----------------------------------------------------------------------------------------------

  subroutine estimate_psamp(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, datafile, outprefix, outfile
    character(len=3)            :: ftxt
    integer(i4b)                :: i, j, l, bin, lmax, numstep, numbin, numcol, specol, numpow
    integer(i4b)                :: numfreq, freq, uncertcol, arr(2), adloc, arloc, powcol, acol, lcol
    integer(i4b)                :: firstbin, lastbin
    real(dp)                    :: f0, rindex, dindex, teori, arhit, adhit
    real(dp)                    :: armin, armax, dar, arcorr, admin, admax, dad, adcorr
    real(sp)                    :: cmb, radio, dust

    real(dp),       allocatable, dimension(:,:,:) :: data
    real(dp),       allocatable, dimension(:,:)   :: power, chisq
    real(dp),       allocatable, dimension(:)     :: freqs, a2t

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'psamp takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, powspecfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
    end if

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    open(13, file=trim(powspecfile))
    read(13,*) numpow, numcol, powcol, acol, lcol 
    write(*,*) numpow,'= numbins'
    allocate(power(numpow,numcol))
    do i = 1, numpow
       read(13,*) power(i,:)
    end do
    close(13)

    ! Read data
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
    open(13, file=trim(datafile))
    read(13,*) numfreq, numbin, numcol, specol, uncertcol
    call assert(numpow==numbin, "mismatching number of datapoints in data and spectrum")
    allocate(freqs(numfreq))
    allocate(data(numfreq,numbin,numcol))
    do i = 1, numfreq
       read(13,*) freqs(i)
       write(*,*) freqs(i)
    end do
    do i = 1, numfreq
       read(13,*) powspecfile
!       if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
       open(14, file=trim(powspecfile))
       read(14,*) outfile !top text line
       do j = 1, numbin
          read(14,*) data(i,j,:)
       end do
       close(14)
    end do
    close(13)
    allocate(a2t(numfreq))
    do freq = 1, numfreq
       a2t(freq) = ant2thermo(freqs(freq))
    end do


    ! Calculate chisq
    firstbin = 6! 6 !1
    lastbin  = numbin !16 !numbin 
    f0 = 100.d0
    rindex = -4.d0
    dindex =  2.d0
    armin = 0.d0
    armax = 0.001d0
    admin = 0.0d0
    admax = 0.001d0
    numstep = 5000
    dar = (armax-armin)/(real(numstep,dp)-1)
    dad = (admax-admin)/(real(numstep,dp)-1)
    allocate(chisq(numstep,numstep))
    chisq = 0.d0
    do j = 1, numstep
       arcorr = armin + (real(j,dp)-1)*dar
       do i = 1, numstep
          adcorr = admin + (real(i,dp)-1)*dad
          do freq = 1, numfreq
             do bin = firstbin, lastbin
                teori = power(bin,powcol) + power(bin,acol)*a2t(freq)*&
                     & (arcorr*(freqs(freq)/f0)**rindex + adcorr*(freqs(freq)/f0)**dindex)
                chisq(i,j) = chisq(i,j) + ((data(freq,bin,specol)-teori)/data(freq,bin,uncertcol))**2
             end do
          end do
       end do
    end do

    ! Find chisq minimum
    arr   = minloc(chisq)
    adloc = arr(1)
    arloc = arr(2)
    arhit = armin + (real(arloc,dp)-1)*dar
    adhit = admin + (real(adloc,dp)-1)*dad
    write(*,*) 'Min. chisq =', real(minval(chisq),sp), 'for Ar =',real(arhit,sp), 'Ad =', real(adhit,sp)
    write(*,*) 'chisq per dof =', real(minval(chisq),sp)/real(numfreq*(lastbin-firstbin+1),sp)

    ! Write best fit theory spectrum to files
    do freq = 1, numfreq
       call int2string(int(freqs(freq),i4b), ftxt)
       outfile = trim(outprefix) // '_bestfit_spec_freq'// ftxt //'.dat'
       open(13, file=trim(outfile))
       do bin = 1, numbin
          cmb   = power(bin,powcol)
          radio = power(bin,acol)*arhit*(freqs(freq)/f0)**rindex
          dust  = power(bin,acol)*adhit*(freqs(freq)/f0)**dindex
          write(13,*) real(power(bin,lcol),sp), cmb+radio+dust, cmb, radio, dust
       end do
       close(13)
       write(*,*) '* Written to file = ', trim(outfile)
    end do

    ! clean up
    deallocate(data, power, chisq, freqs, a2t)

  end subroutine estimate_psamp

  !-----------------------------------------------------------------------------------------------
  ! subroutine binning
  !-----------------------------------------------------------------------------------------------

  subroutine binning(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, binfile, bintype, outfile
    integer(i4b)                :: i, j, l, bin, lmin, lmax, lbeg, lend, numl
    integer(i4b)                :: powcol, ellcol, begcol, endcol, numbin, numcol, numpow
    real(dp)                    :: ell, pow, ant, weight, nevn
    logical(lgt)                :: trivial

    real(dp),       allocatable, dimension(:,:)   :: power, bins

    ! Get parameters
    if (iargc() /= 5) then
       write(*,*) 'bin takes 4 parameters'
       call give_user_info
    else 
       call getarg(2, bintype)
       call getarg(3, binfile)
       call getarg(4, powspecfile)
       call getarg(5, outfile)
    end if
    trivial = .true.
    if (trim(bintype)=='lweight') then
       trivial = .false.
       write(*,*) 'weighting by 2l+1'
    end if

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    open(13, file=trim(powspecfile))
    read(13,*) lmin, lmax, numcol, powcol, ellcol 
    write(*,*) lmin, lmax,'= lmin, lmax in powspec'
    allocate(power(lmin:lmax,numcol))
    do i = lmin, lmax
       read(13,*) power(i,:)
    end do
    close(13)

    ! Read binfile
    if (myid==root) write(*,*) 'Reading from file ', trim(binfile)
    open(13, file=trim(binfile))
    read(13,*) numbin, numcol, begcol, endcol 
    write(*,*) numbin,'= numbins'
    allocate(bins(numbin,numcol))
    do i = 1, numbin
       read(13,*) bins(i,:)
    end do
    close(13)

    ! bin powspec and calculate l(l+1) for each bin
    open(13, file=trim(outfile))
    do bin = 1, numbin
       lbeg = bins(bin,begcol)
       lend = bins(bin,endcol)
       if (lbeg < lmin) cycle
       if (lend > lmax) exit
       numl = lend-lbeg+1
       ell = 0.d0
       pow = 0.d0
       ant = 0.d0
       nevn = 0.d0
       do l = lbeg, lend
          if (trivial) then
             weight = 1.d0
          else
             weight = 2.d0*real(l,dp)+1.d0
          end if
          ell = ell + real(l,dp)
          pow = pow + weight*power(l,powcol)
          ant = ant + weight*real(l*(l+1),dp)
          nevn = nevn + weight
       end do
       ell = ell/real(numl,dp)
       pow = pow/nevn
       ant = ant/(nevn*2*pi)
       write(13,*) real(ell,sp), int(lbeg,i2b), int(lend,i2b), pow, ant 
    end do
    close(13)
    write(*,*) '* Written to file = ', trim(outfile)

    ! clean up
    deallocate(power, bins)

  end subroutine binning

  !-----------------------------------------------------------------------------------------------
  ! subroutine plot_focalplane
  !-----------------------------------------------------------------------------------------------

  subroutine plot_focalplane(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: hornsfile, datafile, outprefix, outfile
    integer(i4b)                :: i, j, numhorns
    real(dp)                    :: tall

!    real(dp),       allocatable, dimension(:,:) :: power, data
    real(dp),       allocatable, dimension(:)   :: rad, theta, horn

    ! Get parameters
    if (iargc() /= 4) then
       write(*,*) 'horns takes 3 parameters'
       call give_user_info
    else 
       call getarg(2, hornsfile)
       call getarg(3, datafile)
       call getarg(4, outprefix)
    end if

    ! Read horn positions
    if (myid==root) write(*,*) 'Reading from file ', trim(hornsfile)
    open(13, file=trim(hornsfile))
    read(13,*) numhorns
    write(*,*) numhorns,'= numhorns'
    allocate(rad(numhorns))
    allocate(theta(numhorns))
    allocate(horn(numhorns))
    do i = 1, numhorns
       read(13,*) horn(i), tall, tall, rad(i), theta(i)
    end do
    close(13)
    theta = theta*pi/180.d0

    ! Read data
    if (myid==root) write(*,*) 'Reading from file ', trim(datafile)
!    open(13, file=trim(datafile))
!    read(13,*) numbin
!    call assert(numpow==numbin, "mismatching number of datapoinnts in data and spectrum")
!    allocate(data(numbin,numcol))
!    open(13, file=trim(datafile))
!    do i = 1, numbin
!       read(13,*) data(i,:)
!    end do
!    close(13)

     ! Write data
    outfile = trim(outprefix)! // '_chisq.dat'
    open(13, file=trim(outfile))
    do i = 1, numhorns
       write(13,*) rad(i)*cos(theta(i)), rad(i)*sin(theta(i)), horn(i)
       write(13,*) ''
    end do
    close(13)
write(*,*) cos(60.d0*pi/180.d0), 'cos 60'
    write(*,*) '* Written to file = ', trim(outfile)

  end subroutine plot_focalplane

  !---------------------------------------------------------------------------------
  ! subroutine merge
  !---------------------------------------------------------------------------------

  subroutine merge_fits(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: outfile, numfiles_in
    integer(i4b)                :: i, numfiles, npix
    real(dp)                    :: nullval
    logical(lgt)                :: anynull

    integer(i4b),       allocatable, dimension(:)   :: ordering, nside, nmaps
    real(dp),           allocatable, dimension(:,:) :: map
    character(len=512), allocatable, dimension(:)   :: infile

    ! Get parameters
    if (iargc() < 2) then
       write(*,*) 'merge takes 3 or more parameters'
       call give_user_info
    else 
       call getarg(2, numfiles_in)
    end if
    read(numfiles_in,*) numfiles
    if (iargc() /= numfiles+3) then
       write(*,*) 'merge takes numfiles+2 parameters'
       call give_user_info
    else 
       allocate(infile(1:numfiles))
       do i = 1, numfiles
          call getarg(2+i, infile(i))
       end do
       call getarg(numfiles+3, outfile)
    end if

    ! Read maps
    allocate(nmaps(numfiles))
    allocate(nside(numfiles))
    allocate(ordering(numfiles))
    do i = 1, numfiles
       npix = getsize_fits(infile(1), nmaps=nmaps(i), ordering=ordering(i), nside=nside(i))
       if (myid==root) write(*,*) i, nside(i), '= nside,', nmaps(i), '= nmaps,' ,ordering(i),'= ordering'
    end do
    ! Check maps
    do i = 2, numfiles
       if (nmaps(i) /= nmaps(1)) then
          if (myid==root) write(*,*) 'Different nside. Quiting'
          stop
       end if
    end do
    ! Actually read maps
    allocate(map(0:npix-1,sum(nmaps)))
    do i = 1, numfiles
       call read_bintab(infile(i), map(:,i:i+nmaps(i)), npix, nmaps(i), nullval, anynull)
    end do
    ! convert to same ordering if necassary
    do i = 2, numfiles
       if (ordering(i) /= ordering(1)) then
          if (myid==root) write(*,*) 'Different ordering. Converting map', i
          if (ordering(1) == 1) then
             call convert_nest2ring(nside(1), map(:,i:i+nmaps(i)))
          else if (ordering(1) == 2) then
             call convert_ring2nest(nside(1), map(:,i:i+nmaps(i)))
          end if
       end if
    end do
    ! Write map
    call write_map(map, ordering(1), trim(outfile))
    write(*,*) '* Map written to file = ', trim(outfile)

    ! Clean up
    deallocate(map,nmaps,ordering,nside,infile)


  end subroutine merge_fits

  !-----------------------------------------------------------------------------------------------
  ! subroutine extract
  !-----------------------------------------------------------------------------------------------

  subroutine extract_fits(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: outfile, string, infile
    integer(i4b)                :: i, numcomp, npix, ordering, nside, nmaps
    real(dp)                    :: nullval
    logical(lgt)                :: anynull

    integer(i4b),       allocatable, dimension(:)   :: comps
    real(dp),           allocatable, dimension(:,:) :: inmap, outmap

    ! Get parameters
    if (iargc() < 2) then
       write(*,*) 'extract takes 4 or more parameters'
       call give_user_info
    else 
       call getarg(2, string)
    end if
    read(string,*) numcomp
    if (iargc() /= numcomp+4) then
       write(*,*) 'extract takes numcomp+3 parameters'
       call give_user_info
    else 
       allocate(comps(1:numcomp))
       call getarg(3, infile)
       do i = 1, numcomp
          call getarg(3+i, string)
          read(string,*) comps(i)
       end do
       call getarg(numcomp+4, outfile)
    end if

    ! Read map
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,' ,ordering,'= ordering'
    ! Check map
    if (nmaps < maxval(comps)) then
       if (myid==root) write(*,*) 'Input map has fewer components than asked for. Quiting'
       stop
    end if
    ! Actually read map
    allocate(inmap(0:npix-1,nmaps))
    call read_bintab(infile, inmap, npix, nmaps, nullval, anynull)
    
    ! Collect components
    allocate(outmap(0:npix-1,numcomp))
    outmap(:,:) = inmap(:,comps)
    ! Write map
    call write_map(outmap, ordering, trim(outfile))
    write(*,*) '* Map written to file = ', trim(outfile)
    
    ! Clean up
    deallocate(inmap,outmap,comps)

  end subroutine extract_fits

  !-----------------------------------------------------------------------------------------------
  ! function 
  !-----------------------------------------------------------------------------------------------

  function compute_sz_term(nu)   ! nu's in GHz
    implicit none

    real(dp),     intent(in)  :: nu
    real(dp)                  :: compute_sz_term
    real(dp) :: x
    real(dp) :: k_B     = 1.3806503d-23
    real(dp) :: h       = 6.626068d-34
    real(dp) :: c       = 2.99792458d8
    real(dp) :: T_0     = 2.725d0

    x = h*nu*1d9 / (k_B*T_0)
    compute_sz_term = x*(exp(x)+1.d0)/(exp(x)-1.d0) - 4.d0

  end function compute_sz_term

  !-----------------------------------------------------------------------------------------------
  ! subroutine add_cluster_old
  !-----------------------------------------------------------------------------------------------

  subroutine add_cluster_old(map, nside, ordering, clusterpar, at, ak, nu)
    implicit none

    real(dp),                  dimension(0:)     :: map
    real(dp),                  dimension(1:)     :: clusterpar
    integer(i4b)                                 :: nside, ordering, npix, numpix, pix, i
    real(dp)                                     :: vec(3), pos(3), fk, ft, ak, at, nu
    real(dp)                                     :: lon, lat, delta, alpha, beta, zeta, gamma, rc, r, radius 
    integer(i4b), allocatable, dimension(:)      :: pixlist

    lon   = clusterpar(1)
    lat   = clusterpar(2)
    delta = 1.051d0
    alpha = 0.3081d0
    beta  = 4.931d0
    zeta  = 0.d0
    gamma = 1.177d0
    rc    = clusterpar(7)*gamma/clusterpar(4)/3.6d0
    fk    = 1.d0
    ft = compute_sz_term(nu)

    rc     = rc*pi/180.d0      ! in radians
    radius = 5*rc
    lon    = lon*pi/180.d0     ! in radians
    lat    = lat*pi/180.d0     ! in radians
    npix = 12*nside**2
    allocate(pixlist(0:npix))
    call ang2vec(pi/2-lat, lon, vec)
    if (ordering==1) then
       call query_disc(nside, vec, radius, pixlist, numpix)
    else if (ordering==2) then
       call query_disc(nside, vec, radius, pixlist, numpix, nest=1)
    end if
    ! loop through pixels inside radius and calculate distance from center
    do i = 0, numpix-1
       pix = pixlist(i)
       if (ordering==1) then
          call pix2vec_ring(nside, pix, pos)
       else if (ordering==2) then
          call pix2vec_nest(nside, pix, pos)
       end if
       call angdist(pos, vec, r)
       r=r/rc
       map(pix) = map(pix) + (ak*fk+at*ft/r**zeta)*r**(-alpha)*(1+r**delta)**(-beta) 
    end do

    open(44, file='profile_old.dat')
    do i = 0, 10000
       r = real(i,dp)/1000.d0 +0.0001
       write(44, *) r, -(ak*fk+at*ft/r**zeta)*r**(-alpha)*(1+r**delta)**(-beta) 
    end do
    close(44)
    stop

    deallocate(pixlist)

  end subroutine add_cluster_old

  !-----------------------------------------------------------------------------------------------
  ! subroutine add_cluster
  !-----------------------------------------------------------------------------------------------

  subroutine add_cluster(tmap, kmap, nside, ordering, clusterpar, at, ak, kamp, tflux, kflux, beam, rmax)
    implicit none

    real(dp),    dimension(0:,0:)              :: tmap, kmap
    real(dp),    dimension(1:),    intent(in)  :: clusterpar
    integer(i4b),                  intent(in)  :: nside, ordering
    real(dp),                      intent(in)  :: at, ak, kamp
    real(dp),                      intent(out) :: tflux, kflux
    real(dp),    dimension(1:,1:,1:), optional :: beam
    real(dp),    dimension(1:),       optional :: rmax

    integer(i4b)                                :: npix, numpix, pix, kix, i, j, k, nstep, numband
    real(dp)                                    :: vec(3), t1, t2
    real(dp)                                    :: lon, lat, rc, radius, perp
    integer(i4b), allocatable, dimension(:)     :: pixlist
    real(dp),     allocatable, dimension(:)     :: thermal, kinetic, beam_comp
    real(dp),     allocatable, dimension(:,:)   :: pos, norm, kart, tart

    numband = size(tmap(1,:)) - 1 
    npix    = 12*nside**2
    lon     = clusterpar(1)*pi/180.d0              ! in radians
    lat     = clusterpar(2)*pi/180.d0              ! in radians
    rc      = clusterpar(7)/clusterpar(4)/3.6d0    ! R500 for cluster measured in deg
    rc      = rc*pi/180.d0                         ! in radians
    radius  = 7.d0*rc + 1.d0*maxval(rmax)
    allocate(pixlist(0:npix))
    call ang2vec(pi/2.d0-lat, lon, vec)
    if (ordering==1) then
       call query_disc(nside, vec, radius, pixlist, numpix)
    else if (ordering==2) then
       call query_disc(nside, vec, radius, pixlist, numpix, nest=1)
    end if
    if (numpix>500000) then
       write(*,*) myid, numpix,'= numpix'
       call cpu_time(t1)
    end if
    ! loop through pixels inside radius and calculate distance from center
    allocate(thermal(0:numpix-1))
    allocate(kinetic(0:numpix-1))
    allocate(pos(0:numpix-1,3))
    tflux = 0.d0
    kflux = 0.d0
    thermal = 0.d0
    kinetic = 0.d0
    do i = 0, numpix-1
       pix = pixlist(i)
       if (ordering==1) then
          call pix2vec_ring(nside, pix, pos(i,:))
       else if (ordering==2) then
          call pix2vec_nest(nside, pix, pos(i,:))
       end if
       call angdist(pos(i,:), vec, perp)
       if (at/=0.d0) thermal(i)  = calculate_cluster_profile(perp, clusterpar(6), clusterpar(7), &
            & clusterpar(3), clusterpar(4))
       if (ak/=0.d0) kinetic(i) = calculate_cluster_profile(perp, clusterpar(6), clusterpar(7), &
            & clusterpar(3), clusterpar(4), kinetic=.true.) *kamp
       ! map(pix) = map(pix) +at*termic +ak*kinetic!*cos(pi/2.d0-lat)
       tmap(pix,0) = tmap(pix,0) + thermal(i)
       kmap(pix,0) = kmap(pix,0) + kinetic(i)
       tflux = tflux + thermal(i)
       kflux = kflux + kinetic(i)
    end do
    ! output total flux
    tflux = tflux*4.d0*pi/real(npix, dp)
    kflux = kflux*4.d0*pi/real(npix, dp)

    ! convolve with beam in pixel space
    if (present(beam)) then
       allocate(norm(numband, 0:numpix-1))
       allocate(tart(numband, 0:numpix-1))
       allocate(kart(numband, 0:numpix-1))
       allocate(beam_comp(numband))
       norm = 0.d0
       kart = 0.d0
       tart = 0.d0
       do i = 0, numpix-1
          do k = i, numpix-1
             call angdist(pos(i,:), pos(k,:), perp)
             do j = 1, numband
                if (perp<rmax(j)) then
                   beam_comp(j) = splint(beam(:,j,1), beam(:,j,2), beam(:,j,3), perp)
                   norm(j,i) = norm(j,i) + beam_comp(j)
                   tart(j,i) = tart(j,i) + beam_comp(j)*thermal(k)
                   kart(j,i) = kart(j,i) + beam_comp(j)*kinetic(k)
                   if (k /= i) then
                      norm(j,k) = norm(j,k) + beam_comp(j)
                      tart(j,k) = tart(j,k) + beam_comp(j)*thermal(i)
                      kart(j,k) = kart(j,k) + beam_comp(j)*kinetic(i)
                   end if
                end if
             end do
          end do
       end do
       tart = tart/norm
       kart = kart/norm
       do i = 0, numpix-1
          pix = pixlist(i)
          tmap(pix,1:numband) = tmap(pix,1:numband) + tart(:,i)
          kmap(pix,1:numband) = kmap(pix,1:numband) + kart(:,i)
       end do
       deallocate(norm, beam_comp, kart, tart)
    end if
    deallocate(thermal, kinetic, pos)

    if (numpix>500000) then 
       call cpu_time(t2)
       write(*,*) myid, numpix,'= numpix finished in time:', t2-t1
    end if

!    open(44, file='profile_termic.dat')
!    open(45, file='profile_kinetic_m2.dat')
!    do i = 1, 20000
!       perp = real(i,dp)/4000.d0
!       termic  = calculate_cluster_profile(perp, clusterpar(6), clusterpar(7), clusterpar(3), clusterpar(4))
!       kinetic = calculate_cluster_profile(perp, clusterpar(6), clusterpar(7), clusterpar(3), clusterpar(4), gammain=-2.d0)
!       write(44, *) perp, termic
!       write(45, *) perp, kinetic
!    end do
!    close(44)
!    close(45)
!    stop

    deallocate(pixlist)
  end subroutine add_cluster

  !-----------------------------------------------------------------------------------------------
  ! function calculate cluster profile
  !-----------------------------------------------------------------------------------------------

  function calculate_cluster_profile(r, M500, R500, redshift, scale, p0in, c500in, gammain, alphain, betain, kinetic)
    implicit none

    real(dp)                                     :: calculate_cluster_profile
    real(dp),                        intent(in)  :: r, M500, R500, redshift, scale
    real(dp),              optional, intent(in)  :: p0in, c500in, gammain, alphain, betain
    logical(lgt),          optional, intent(in)  :: kinetic

    integer(i4b)                                 :: i, j, nstep
    real(dp)                                     :: p0, c500, gamma, alpha, beta
    real(dp)                                     :: hubble, OmegaM, OmegaL, h70, rc, kintegral, tintegral, rmax
    real(dp)                                     :: mass, hz2, x, y, z, integral, dz, a, b, a_t, a_k, b_t, b_k
    logical(lgt)                                 :: kin

    if (present(kinetic)) then
       kin = kinetic
    else
       kin = .false.
    end if
    !6.628d-6/Mpc = 1.65d-3keV/cm^3 *thompsonscatter/m_electron/c^2
    hubble = 0.7d0
    OmegaM = 0.3d0
    OmegaL = 0.7d0
    h70    = hubble/0.7d0
    p0     = 8.403d0*h70**(-3.d0/2.d0)
    c500   = 1.177d0
    gamma  = 0.3081d0
    alpha  = 1.051d0
    beta   = 5.4905d0
    if (present(p0in))    p0    = p0in*h70**(-3.d0/2.d0)
    if (present(c500in))  c500  = c500in
    if (present(gammain)) gamma = gammain
    if (present(alphain)) alpha = alphain
    if (present(betain))  beta  = betain
    mass   = M500/3.d0*h70
    hz2    = OmegaL + OmegaM*(1+redshift)**3
    rc     = R500/scale/3.6d0    !R500 for cluster measured in deg
    rc     = rc*pi/180.d0        ! in radians
    y = r/rc 
    if (y<0.01d0) then
       nstep = 10000
    else if (y<0.1d0) then
       nstep = 1000
    else
       nstep = 100
    end if
    integral= 0.d0
    rmax = 15.d0
    dz = rmax/real(nstep,dp)
    a = 0.d0
    do j = -nstep, 0 !half the integral
       z = real(j,dp)*dz
       x = sqrt(y**2+z**2)
       b = mass**(0.887d0-0.22d0*(2.d0*x)**3/(1.d0+(2.d0*x)**3)) *(c500*x)**(-gamma)/(1.d0+(c500*x)**alpha)**((beta-gamma)/alpha)
       if (kin) b = b/(1+0.75d0*x)**-1.6d0
       integral = integral +(a+b)*dz!/2.d0
       a = b
    end do
    if (kin) integral = integral*511.d0/11.2d0/R500**2/hubble**2       /1000.d0  !obs
    calculate_cluster_profile = 6.628d-6*hz2**(4.d0/3.d0)*p0*h70**2*integral*R500*2.725d6


    return
    open(44, file='profile_pressure.dat')
    open(45, file='profile_temperature.dat')
    do i = 1, 20000
       x = real(i,dp)/4000.d0
       b_t = mass**(0.887d0-0.22d0*(2.d0*x)**3/(1.d0+(2.d0*x)**3)) *(c500*x)**(-gamma)/(1.d0+(c500*x)**alpha)**((beta-gamma)/alpha)
       b_k = (1+0.75d0*x)**-1.6d0

       tintegral = 6.628d-6*hz2**(4.d0/3.d0)*p0*h70**2*b_t*511/6.6524d-25/3.0856d24 
       kintegral = b_k*11.2d0*R500**2*hubble**2 
       write(44, *) x, tintegral
       write(45, *) x, kintegral
    end do
    close(44)
    close(45)
    call mpi_finalize(i)
    stop

    open(44, file='profile_thermal2.dat')
    open(45, file='profile_kineticnew2.dat')
    do i = 1, 28000
       y = real(i,dp)/4000.d0
       nstep = 3000
       tintegral= 0.d0
       kintegral= 0.d0
       rmax = 15.d0
       dz = rmax/real(nstep,dp)
       a_t = 0.d0
       a_k = 0.d0
       do j = -nstep, 0
          z= real(j,dp)*dz
          x=sqrt(y**2+z**2)
          b_t = mass**(0.887d0-0.22d0*(2.d0*x)**3/(1.d0+(2.d0*x)**3)) *&
               & (c500*x)**(-gamma)/(1.d0+(c500*x)**alpha)**((beta-gamma)/alpha)
          b_k = b_t/(1+0.75d0*x)**-1.6d0
          tintegral = tintegral +(a_t+b_t)*dz!/2.d0
          kintegral = kintegral +(a_k+b_k)*dz!/2.d0
          a_t = b_t
          a_k = b_k
       end do
       tintegral = 6.628d-6*hz2**(4.d0/3.d0)*p0*h70**2*tintegral*R500   *2.d0
       kintegral = 6.628d-6*hz2**(4.d0/3.d0)*p0*h70**2*kintegral*R500 *511.d0/11.2d0/R500**2/hubble**2  /1000.d0
       write(44, *) y, tintegral
       write(45, *) y, kintegral
       write(*,*) i, a_t, a_k
    end do
    close(44)
    close(45)
    call mpi_finalize(i)
    stop

  end function calculate_cluster_profile

  !-----------------------------------------------------------------------------------------------
  ! subroutine clustersims
  !-----------------------------------------------------------------------------------------------

  subroutine clustersims(unit)
    implicit none

    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: clusterfile, outfile, outprefix, powspecfile, noisefile, beamfile
    character(len=512)          :: covfile, noiseseed_in, cmbseed_in, clusterfile2, ak_in, at_in, nside_in
    integer(i4b)                :: ordering, ordering2, nside, nside2, nmaps, nmaps2, npix, npix2
    integer(i4b)                :: dnside, dnside2, dnpix
    integer(i4b)                :: i, j, l, c, numclusters, lmax, cmbseed, noiseseed, numband, nbeam
    real(dp)                    :: t1, t2, at, ak, eta, tflux, kflux, maxsize, beamcutoff, tmp
    logical(lgt)                :: addbeam, addcmb, addnoise, fullsim, onlyclu, addcov
    type(task_list)             :: tasks

    real(dp),          allocatable, dimension(:)     :: beam_rmax, freq, ampt, ampk, kamp
    real(dp),          allocatable, dimension(:,:)   :: flux, kmap, tmap, iflux, ikmap, itmap, cov, Lmat
    real(dp),          allocatable, dimension(:,:)   :: clusterpar, noisemap2, cmbmap, usemap
    real(dp),          allocatable, dimension(:,:)   :: power, beam, map, dmap, inbeam, inpower, noisemap
    real(dp),          allocatable, dimension(:,:,:) :: beam_r
    character(len=3),  allocatable, dimension(:)     :: bandtxt
    integer(i4b),      allocatable, dimension(:)     :: bandmax
 
    ! Get parameters
    if (iargc() /= 13) then
       write(*,*) 'clusim takes 12 parameters'
       call give_user_info
    else 
       call getarg(2, clusterfile)
       call getarg(3, clusterfile2)
       call getarg(4, powspecfile)
       call getarg(5, beamfile)
       call getarg(6, noisefile)
       call getarg(7, covfile)
       call getarg(8, nside_in)
       call getarg(9, outprefix)
       call getarg(10, at_in)
       call getarg(11, ak_in)
       call getarg(12, cmbseed_in)
       call getarg(13, noiseseed_in)
    end if
    read(nside_in,*)     nside
    read(at_in,*)        at
    read(ak_in,*)        ak
    read(cmbseed_in,*)   cmbseed
    read(noiseseed_in,*) noiseseed
    call initialize_random_seeds(MPI_COMM_WORLD, cmbseed,   rng_handle_cmb)
    call initialize_random_seeds(MPI_COMM_WORLD, noiseseed, rng_handle_noise)
    if (trim(clusterfile2)=='fullsim') then
       fullsim = .true.
       onlyclu = .false.
       if (myid==root) write(*,*) 'Making cluster simulations from scratch'
    else if (trim(clusterfile2)=='onlyclu') then 
       fullsim = .true.
       onlyclu = .true.
       if (myid==root) write(*,*) 'Only making and dumping cluster simulations'
    else
       fullsim = .false.
       onlyclu = .false.
       if (myid==root) write(*,*) 'Taking cluster simulations from input map'
    end if
    addbeam  = .true.
    addcmb   = .true.
    addnoise = .true.
    addcov = .true.
    if (trim(beamfile)=='nobeam') then
       addbeam = .false.
       if (myid==root) write(*,*) 'Running without adding beam'
    end if
    if (trim(powspecfile)=='nocmb') then
       addcmb = .false.
       if (myid==root) write(*,*) 'Running without adding cmb'
    end if
    if (trim(noisefile)=='nonoise') then
       addnoise = .false.
       if (myid==root) write(*,*) 'Running without adding noise'
    end if
    if (trim(covfile)=='nocov') then
       addcov = .false.
       if (myid==root) write(*,*) 'Running without adding velocity covariance'
    end if
    ordering = 1 !ring
    nmaps    = 1
    npix     = 12*nside**2
    lmax     = 3*nside
    if (myid==root .and. fullsim) write(*,*) nside, '= nside for simulations;', lmax, '= lmax'

    if (.not. onlyclu) then
       ! Read beams and freqs
       if (addbeam) then 
          if (myid==0) write(*,*) 'Reading from file ', trim(beamfile)
          open(13, file=beamfile)
          read(13,*) numband
          if (myid==root) write(*,*) numband, '= numband'
          allocate(inbeam(0:lmax,1:4))
          allocate(bandtxt(numband))
          allocate(bandmax(numband))
          allocate(freq(numband))
          allocate(beam(0:lmax,1:numband))
          beam = 0.d0
          inbeam = 0.d0
          bandmax = lmax
          do i = 1, numband
             read(13,*) bandtxt(i)
          end do
          do i = 1, numband
             read(13,*) freq(i)
          end do
          do i = 1, numband
             read(13,*) beamfile
             if (myid==root) write(*,*) trim(beamfile)
             call read_beam(beamfile, inbeam)
             do l = 0, lmax
                if (inbeam(l,1) == 0.d0) then
                   bandmax(i) = l-1
                   if (myid==root) write(*,*) bandmax(i), '= lmax beam ', bandtxt(i)
                   exit
                end if
             end do
             beam(0:bandmax(i),i) = inbeam(0:bandmax(i),1)
          end do
          close(13)
       end if
       deallocate(inbeam)
       ! spline beams
       nbeam = 10000     ! spline resolution
       maxsize = 5.d0*pi/180.d0  ! in radians
       beamcutoff = 0.0001
       allocate(beam_r(nbeam,numband,3))
       allocate(beam_rmax(numband))
       beam_r=0.d0
       do i= 1, nbeam
          beam_r(i,:,1) = real(i,dp) * maxsize / real(nbeam,sp) !in arcmin
       end do
       do j = 1, numband
           call beam_to_radial(beam(0:bandmax(j),j), beam_r(:,j,1), beam_r(:,j,2))
           call spline(beam_r(:,j,1), beam_r(:,j,2), 0.d0, 1.d30, beam_r(:,j,3))
       end do
       ! truncate splined beam
       do j = 1, numband
          beam_rmax(j) = 1.d0
          do i = 1, nbeam
             if (beam_r(i,j,2) < beamcutoff*beam_r(1,j,2)) exit
          end do
          call assert(i<nbeam, "beam falls of too slowly, choose larger maxsize")
          !write(*,*) i, 'stop', beam_r(i,j,2)
          beam_rmax(j) = beam_r(i,j,1)
          if (myid==0) write(*,*) bandtxt(j), beam_rmax(j)*180.d0*60.d0/pi, '= beam_rmax in arcmin'
       end do
!       deallocate(beam)

       ! read noise rms maps and convert to ringed
       if (addnoise) then 
          if (myid==root) write(*,*) 'Reading from file ', trim(noisefile)
          open(13, file=noisefile)
          read(13,*) i
          call assert(i == numband, "Numband mismatch between beam and noise")
          read(13,*) noisefile
          call read_map(noisemap2, ordering2, noisefile, nside=dnside, nmap=nmaps2)
          if (myid==root) write(*,*) dnside,'= nside_noise', nmaps2,'= nmaps', ordering2,'= ordering'
          if (ordering2==2) call convert_nest2ring(dnside, noisemap2)
          allocate(noisemap(0:12*dnside**2-1,numband))
          noisemap(:,1) = noisemap2(:,1)
          deallocate(noisemap2)
          do i = 2, numband
             read(13,*) noisefile
             call read_map(noisemap2, ordering2, noisefile, nside=dnside2, nmap=nmaps2)
             if (myid==root) write(*,*) dnside2,'= nside_noise', nmaps2,'= nmaps', ordering2,'= ordering'
             call assert(dnside == dnside2, "Nside mismatch between noise maps")
             if (ordering2==2) call convert_nest2ring(dnside, noisemap2)
             noisemap(:,i) = noisemap2(:,1)
             deallocate(noisemap2)
          end do
          close(13)
       end if
       ! Read powspec
       if (addcmb) then
          if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
          allocate(inpower(0:lmax,1:4))
          inpower = 0.d0
          call read_powspec(powspecfile, inpower)
          if (myid==root) then
             open(13, file='power.dat', recl=1024)
             do l = 0, lmax
                write(13,*) l, inpower(l,:)*l*(l+1)/2/pi
                if (inpower(l,1) ==0.d0 .and. l>1) then
                   lmax = l-1
                   write(*,*) lmax, '= lmax power'
                   exit
                end if
             end do
             close(13)
          end if
          if (addbeam) lmax = min(lmax, maxval(bandmax))
          if (addbeam .and. myid==root) write(*,*) lmax, '= lmax'
          allocate(power(0:lmax,1:1))
          power(0:lmax,1) = inpower(0:lmax,1)
          deallocate(inpower)
       end if
    end if

    ! set cluster amplitudes
    allocate(ampt(0:numband))
    allocate(ampk(0:numband))
    ampt(0) = at
    do i = 1, numband
       ampt(i) = compute_sz_term(freq(i))*at
       if (myid==0) write(*,*) ampt(i), bandtxt(i)
    end do
    ampk(:) = ak 
    if (myid==0 .and. (.not. onlyclu)) then 
       write(*,*) at,'= thermal amplitude'
       write(*,*) ak,'= kinetic amplitude'
    end if

    if (fullsim) then
       ! Read list of clusters
       if (myid==root) write(*,*) 'Reading from file ', trim(clusterfile)
       open(unit,file=trim(clusterfile))
       read(unit, *) numclusters
       if (myid==root) write(*,*) numclusters, '= number of clusters'
       allocate(clusterpar(numclusters, 7))
       do i=1,numclusters
          read(unit, *) clusterpar(i,:)
       end do
       close(unit)
       if (ak/=0.d0) then
          ! generate velocity amplitudes
          allocate(kamp(numclusters))
          do i = 1, size(kamp)
             kamp(i) = rand_gauss(rng_handle_noise)
          end do
          ! Read cov matrix for kinetic sz
          if (addcov) then
             if (myid==root) write(*,*) 'Reading from file ', trim(covfile)
             allocate(cov(numclusters, numclusters))
             allocate(Lmat(numclusters, numclusters))
             cov = 0.d0
             Lmat = 0.d0
             open(unit,file=trim(covfile))
             do i = 1, numclusters
                do j = i, numclusters
                   read(unit,*) tmp, tmp, cov(j,i)
                   cov(i,j) = cov(j,i)
                end do
             end do
             close(unit)
             call cholesky_decompose(cov, Lmat)
             kamp = matmul(Lmat, kamp)
             deallocate(cov, Lmat)
          end if
       end if

       ! Initialize map
       allocate(iflux(1:2, 1:numclusters))
       allocate(itmap(0:npix-1, 0:numband))
       allocate(ikmap(0:npix-1, 0:numband))
       iflux = 0.d0
       itmap = 0.d0
       ikmap = 0.d0
       ! add clusters
       outfile=trim(outprefix)//'_lock.dat'
       call init_task_list(tasks, outfile, numclusters, MPI_COMM_WORLD)
       if (myid==root) write(*,*) myid, 'Ready to simulate clusters'
       call cpu_time(t1)
       do while(get_next_task(tasks, i))
          call add_cluster(itmap, ikmap, nside, ordering, clusterpar(i,:), at, ak, kamp(i), &
               & iflux(1,i), iflux(2,i), beam_r, beam_rmax)
          if (mod(i,1) == 20.d0) write(*,*) myid, 'finished cluster', int(i,i2b), real(iflux(:,i),sp) 
       end do
       deallocate(beam_r, beam_rmax)
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       call cpu_time(t2)
       if (myid==root) write(*,*) myid, 'clustersims', t2-t1 
       if (myid==root) write(*,*) myid, 'all clusters finished. starting to reduce'
       ! Collect to one root map
       call cpu_time(t1)
       allocate(tmap(0:npix-1, 0:numband))
       call mpi_reduce(itmap,  tmap, size(itmap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       deallocate(itmap)
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       call cpu_time(t2)
       if (myid==root) write(*,*) myid, 'reduced 1(3) tmap', t2-t1 
       allocate(kmap(0:npix-1, 0:numband))
       call mpi_reduce(ikmap,  kmap, size(ikmap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       deallocate(ikmap)
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       call cpu_time(t1)
       if (myid==root) write(*,*) myid, 'reduced 2(3) kmap', t1-t2
       allocate(flux(1:2, 1:numclusters))
       call mpi_reduce(iflux,  flux, size(iflux), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       deallocate(iflux)
       call cpu_time(t2)
       if (myid==root) then 
          write(*,*) myid, 'reduced 3(3) flux', t2-t1
          ! Output fluxes
          outfile=trim(outprefix)//'_fluxes.dat'
          open(17, file=outfile)
          open(18, file=trim(outprefix)//'_fluxratio.dat')
          do i = 1, numclusters
             write(17,*) i, flux(1,i), flux(2,i)
             write(18,*) clusterpar(i,3), flux(2,i)/flux(1,i)
          end do
          close(17)
          close(18)
          write(*,*) '* Fluxes written to file = ', trim(outfile)
       end if
       deallocate(flux, clusterpar)
       if (myid==root) then 
          ! Output maps
          outfile=trim(outprefix)//'_thermal.fits'
          call write_map(tmap, ordering, trim(outfile))
          write(*,*) '* Map written to file = ', trim(outfile)
          outfile=trim(outprefix)//'_kinetic.fits'
          call write_map(kmap, ordering, trim(outfile))
          write(*,*) '* Map written to file = ', trim(outfile)
       end if
       if (onlyclu) return

    else ! if not fullsim, read clustersims from file
       if (myid==0) write(*,*) 'Reading from ', trim(clusterfile)
       call read_map(tmap, ordering2, clusterfile, nside=nside, nmap=nmaps2)
       if (myid==root) write(*,*) nside, '= nside,', nmaps2, '= nmaps,',ordering2,'= ordering'       
       if (ordering2==2) call convert_nest2ring(nside, tmap)

       if (myid==0) write(*,*) 'Reading from ', trim(clusterfile2)
       call read_map(kmap, ordering2, clusterfile2, nside=nside2, nmap=nmaps2)
       if (myid==root) write(*,*) nside2, '= nside,', nmaps2, '= nmaps,',ordering2,'= ordering'       
       call assert(nside == nside2, "Nside mismatch between thermal and kinetic map")
       if (ordering2==2) call convert_nest2ring(nside, kmap)
       npix = 12*nside**2
    end if

    if (myid==0) then
       ! Make clustermap
       allocate(map(0:npix-1, 0:numband))
       do j = 0, numband
          map(:,j) = ampt(j)*tmap(:,j) +ampk(j)*kmap(:,j)
       end do
       deallocate(tmap, kmap)
       ! output clustermap
       outfile=trim(outprefix)//'_clustermap.fits'
       call write_map(map, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       
       ! adding cmb and beam
       if (addcmb) then
          allocate(cmbmap(0:npix-1,1))
          cmbmap = 0.d0
          call add_cmbandbeam(cmbmap, power, beam, nside, .false., addcmb)
          if (addbeam) then
             allocate(usemap(0:npix-1,1))
             do j = 1, numband
                write(*,*) j, 'adding..'
                usemap = cmbmap
                call add_cmbandbeam(usemap, power, beam(0:bandmax(j),j:j), nside, addbeam, .false.)
                map(:,j) = map(:,j) + usemap(:,1)
             end do
             deallocate(usemap)
          end if
          deallocate(cmbmap)
          ! Output signalmap
          outfile=trim(outprefix)//'_signalmap.fits'
          call write_map(map, ordering, trim(outfile))
          write(*,*) '* Map written to file = ', trim(outfile)
       end if
       
       ! downgrade and add noise
       if (.not. addnoise) dnside = nside/4
       dnpix  = 12*dnside**2
       allocate(dmap(0:dnpix-1, 0:numband))
       write(*,*) 'Downgrading from', nside,'to', dnside
       call udgrade_ring(map, nside, dmap, dnside)
       ! Output downgraded signalmap
       outfile=trim(outprefix)//'_dsigmap.fits'
       call write_map(dmap, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       if (addnoise) then
          do j = 1, numband
             do i = 0, dnpix -1
                eta = rand_gauss(rng_handle_noise)
                dmap(i,j) = dmap(i,j) + noisemap(i,j)*eta
             end do
          end do
       end if
       ! Output final map
       do j = 1, numband
          outfile=trim(outprefix)//'_'//trim(bandtxt(j))//'_dmap.fits'
          call write_map(dmap(:,j), ordering, trim(outfile)) 
          write(*,*) '* Map written to file = ', trim(outfile)
       end do
    end if
    
    !clean up
    if (allocated(map))   deallocate(map)
    if (allocated(dmap))  deallocate(dmap)
    if (allocated(power)) deallocate(power)
    if (allocated(beam))  deallocate(beam)

  end subroutine clustersims

  !----------------------------------------------------------------------------------------------
  ! subroutine make_sims
  !----------------------------------------------------------------------------------------------

  subroutine make_sims(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: powspecfile, beamfile, maskfile, outdir, outfile
    character(len=512)          :: nside_in, numsim_in, numside_in, lmax_in, pixwin_in
    character(len=3)            :: simtxt
    character(len=4)            :: sidetxt
    integer(i4b)                :: ordering, nside, nside2, nmaps, nmaps2, npix, npix2
    integer(i4b)                :: dnside, dnpix
    integer(i4b)                :: i, j, l, m, s, lmax, pol, numsim, numpatches, ns, numsides
    real(dp)                    :: t1, t2
    logical(lgt)                :: nobeam, applypixwin

    character(len=1024), allocatable, dimension(:) :: masklist
    character(len=8),    allocatable, dimension(:) :: patchtxt
    real(dp),          allocatable, dimension(:,:) :: power, beam, inmask, map, outmap, downmask
    real(dp),          pointer,     dimension(:,:) :: pixwin
    integer(i4b),      allocatable, dimension(:,:) :: masks
    integer(i4b),      allocatable, dimension(:)   :: healvec
 
    !obs sjekk om noise og cmb forskjellig handle
    call initialize_random_seeds(MPI_COMM_WORLD, 735909, rng_handle)  

    ! Get parameters
    if ((iargc() < 9) .or. (iargc() > 10)) then
       if (myid==root) then
          write(*,*) 'sims takes 8(9) parameters'
          write(*,*) ''
          write(*,*) 'sims powspec beam pixwin mask nside numside outdir numsim (lmax)'
          write(*,*) ''
          write(*,*) 'Ex: sims pow.fits beam.fits true nomask 64 1 sims 10 150'
          write(*,*) ''
       end if
    call mpi_finalize(ierr)
    stop
    else if (iargc() == 10) then
       call getarg(10, lmax_in)
       read(lmax_in,*) lmax
    end if
    call getarg(2, powspecfile)
    call getarg(3, beamfile)
    call getarg(4, pixwin_in)
    call getarg(5, maskfile)
    call getarg(6, nside_in)
    call getarg(7, numside_in)
    call getarg(8, outdir)
    call getarg(9, numsim_in)
    read(nside_in,*) nside
    read(numside_in,*) numsides
    read(numsim_in,*) numsim
    read(pixwin_in,*) applypixwin
    nmaps = 3 ! OBS hardcoding polarization
    if (myid==root) write(*,*) nmaps, ' = nmaps' 
    if (myid==root) write(*,*) nside, ' = nside' 
    npix  = 12*nside**2
    if (iargc() == 9) lmax  = 3*nside
    if (myid==root) write(*,*) lmax, ' = lmax' 
    ! Check for beam
    nobeam = .false.
    if (trim(beamfile)=='nobeam') then
       nobeam = .true. 
       if (myid==root) write(*,*) 'Running without beam'
    end if

    ! Read powspec
    if (myid==root) write(*,*) 'Reading from file ', trim(powspecfile)
    allocate(power(0:lmax,1:4))
    call read_powspec(powspecfile, power)
    if (myid==root) then
       open(13, file='power.dat', recl=1024)
       do l = 0, lmax
          write(13,*) l, power(l,:)*l*(l+1)/2/pi
       end do
       close(13)
    end if

    ! Read beam and apply
    if (.not. nobeam) then 
       ! Read beam
       if (myid==0) write(*,*) 'Reading from file ', trim(beamfile)
       allocate(beam(0:lmax,1:4))
       call read_beam(beamfile, beam)
       if (myid==root) then       
          open(13, file='beam.dat', recl=1024)
          do l = 0, lmax
             write(13,*) l, beam(l,:)
          end do
          close(13)
       end if
       ! Convolve power spectrum with beam
       power(:,1:3) = power(:,1:3)*beam(:,1:3)**2
       power(:,4)   = power(:,4)*beam(:,1)*beam(:,2)
       if (myid==root) then
          open(13, file='powerwithbeam.dat', recl=1024)
          do l = 0, lmax
             write(13,*) l, power(l,:)*l*(l+1)/2/pi
          end do
          close(13)
       end if
    end if
    deallocate(beam) ! done with it after convolving the power

    ! Read pixwin and apply
    if (applypixwin) then 
       ! Read pixwin
       if (myid==0) write(*,*) 'Applying pixel window'
       call read_pixwin(nside, nmaps, pixwin)
       if (myid==root) then       
          open(13, file='pixwin.dat', recl=1024)
          do l = 0, lmax
             write(13,*) l, pixwin(l,:)
          end do
          close(13)
       end if
       ! Convolve power spectrum with pixwin
       power(:,1:3) = power(:,1:3)*pixwin(0:lmax,1:3)**2
       power(:,4)   = power(:,4)*pixwin(0:lmax,1)*pixwin(0:lmax,2)
       if (myid==root) then
          open(13, file='powerwithbeamandpixwin.dat', recl=1024)
          do l = 0, lmax
             write(13,*) l, power(l,:)*l*(l+1)/2/pi
          end do
          close(13)
       end if
    end if


    ! Read masklist
    if (myid==0) write(*,*) 'Reading from file ', trim(maskfile)
    numpatches=0
    open(unit,file=trim(maskfile), recl=1024)
    read(unit, *) numpatches
    allocate(patchtxt(0:numpatches))
    patchtxt(0)='fullsky'
    allocate(masklist(numpatches))
    do i = 1, numpatches
       read(unit,*) patchtxt(i)
    end do
!numpatches=1
    do i = 1, numpatches
       read(unit,*) masklist(i)
    end do
    close(unit)
    ! Read masks
    allocate(masks(0:npix-1, 1:numpatches))
    masks = 0
    do i = 1, numpatches
       if (myid==root) write(*,*) patchtxt(i)
       call read_map(inmask, ordering, masklist(i), nside=nside2, nmap=nmaps2)
       if (myid==root) write(*,*) nside2,'= nside_mask', nmaps2,'= nmaps', ordering,'= ordering'
       if (ordering == 1) then
          call convert_ring2nest(nside2, inmask)
          ordering = 2
       end if
       if (nside2/=nside) then
          allocate(downmask(0:npix-1, 1:nmaps2))
          call udgrade_nest(inmask, nside2, downmask, nside)
          deallocate(inmask)
          where (abs(downmask(:,nmaps)-1.d0) < 1d-25) masks(:,i) = 1   ! OBS lik nside
          deallocate(downmask)
       else 
          where (abs(inmask(:,nmaps)-1.d0) < 1d-25) masks(:,i) = 1   ! OBS lik nside
          deallocate(inmask)
       end if
    end do
    deallocate(masklist)

    ! Make output directories
    do ns = 1, numsides
       dnside = nside/(2**(ns-1))
       call int2string(dnside, sidetxt)
       do i = 0, numpatches
          outfile = trim(outdir) // '/' // sidetxt // '/' // trim(patchtxt(i)) // '/sim.fits'
          call mkdirs(outfile, .true.)     
       end do
    end do

    !Make simulations
    allocate(map(0:npix-1, 1:nmaps))
    do s = 1+myid, numsim, numprocs
       call int2string(s, simtxt)
       write(*,*) 'sim', simtxt, myid, '= myid'
       call make_mapsim(map, power, nside)
!       map = real(s, dp)

       ! Output maps
       call int2string(nside, sidetxt)
       outfile = trim(outdir) // '/' // sidetxt // '/' // trim(patchtxt(0)) // '/sim' // simtxt // '.fits'
       !call cpu_time(t1)
       call write_map(map, ordering, trim(outfile))
       !call cpu_time(t2)
       !if (myid==0) write(*,*) 'Time to write fullsky map to disk:', t2-t1
! What is this doing here?
!       deallocate(outmap)
       do i =1, numpatches
          call mask2heal(masks(:,i), healvec)
          allocate(outmap(1:size(healvec), 1:nmaps))
          do j = 1, nmaps
             outmap(:,j) =  map(healvec,j)
          end do
          outfile = trim(outdir) // '/' // sidetxt // '/' // trim(patchtxt(i)) // '/sim' // simtxt // '.fits'
          !call cpu_time(t1)
          call write_map(outmap, healvec, nside, ordering, trim(outfile))
          !call cpu_time(t2)
          !if (myid==0) write(*,*) 'Time to write ',trim(patchtxt(i)), t2-t1 
          deallocate(healvec, outmap)
       end do

       ! Downgrade and output
!       do ns = 2, numsides
!          dnside = nside/(2**(ns-1))
!          dnpix  = 12*dnside**2
!          allocate(dmap(0:dnpix-1, 1:nmaps))
!          call udgrade_nest(map, nside, dmap, dnside)
!          allocate(dmasks(0:dnpix-1, 1:numpatches))
!       end do

    end do

    deallocate(patchtxt,power,map,masks)

  end subroutine make_sims

  !-----------------------------------------------------------------------------------------------
  ! subroutine make_mapsim
  !-----------------------------------------------------------------------------------------------

  subroutine make_mapsim(map, power, nside)
    implicit none
    real(dp),     allocatable, dimension(:,:)   :: map, power
    integer(i4b)                                :: nside

    integer(i4b)                                :: lmax, l, m, n, i
    real(i4b)                                   :: cl(4), t1, t2 
    real(dp),     allocatable, dimension(:)     :: eta, x, y
    real(dp),     allocatable, dimension(:,:)   :: matrix, Lmat
    complex(dpc), allocatable, dimension(:,:,:) :: alm

    n = 3
    lmax = size(power(:,1)) - 1
    allocate(matrix(n,n))
    allocate(Lmat(n,n))
    allocate(eta(n), x(n), y(n))
    allocate(alm(1:n, 0:lmax, 0:lmax))
    !call cpu_time(t1)
    do l = 0, lmax
       matrix = 0.d0
       matrix(1,1) = power(l, 1)
       matrix(2,2) = power(l, 2)
       matrix(3,3) = power(l, 3)
       matrix(1,2) = power(l, 4)
       matrix(2,1) = power(l, 4)
       if (all(matrix==0.d0)) then
          alm(:,l,:) = cmplx(0,0)
       else
          call cholesky_decompose(matrix, Lmat)
          do m = 0, l
             ! Generate gaussian noise
             do i = 1, n
                eta(i) = rand_gauss(rng_handle) 
             end do
             ! .. with the right cov matrix
             x = matmul(Lmat,eta)
             if (m==0) then
                alm(:,l,m) = cmplx(x,0)
             else
                do i = 1, n
                   eta(i) = rand_gauss(rng_handle) 
                end do
                y = matmul(Lmat,eta)
                alm(:,l,m) = cmplx(x,y)/sqrt(2.d0)
             end if
          end do
       end if
    end do
    !call cpu_time(t2)
    !if (myid==0) write(*,*) 'Time almsim:', t2-t1 

    if (myid==0) then
      ! call cpu_time(t1)
       open(13, file='almpower.dat', recl=1024)
       do l = 0, lmax
          cl(1:3) = alm(:,l,0)*conjg(alm(:,l,0))
          cl(4)   = alm(1,l,0)*conjg(alm(2,l,0))
          do m = 1, l
             cl(1:3) = cl(1:3) + 2*alm(:,l,m)*conjg(alm(:,l,m))
             cl(4)   = cl(4) + alm(1,l,m)*conjg(alm(2,l,m)) + alm(2,l,m)*conjg(alm(1,l,m))
          end do
          cl = cl/(2.d0*l+1)
          write(13,*) l, cl*l*(l+1)/2.d0/pi
       end do
       close(13)
       !call cpu_time(t2)
       !write(*,*) 'Time to write almpower to disk:', t2-t1 
    end if
    
    !call cpu_time(t1)
    call alm2map(nside, lmax, lmax, alm, map)
    !call cpu_time(t2)
    !if (myid==0) write(*,*) 'Time alm2map:', t2-t1 
    deallocate(eta,x,y,matrix,Lmat,alm)
    !call cpu_time(t1) 
    call convert_ring2nest(nside, map)
    !call cpu_time(t2) 
    !if (myid==0) write(*,*) 'Time ring2nest:', t2-t1 

  end subroutine make_mapsim

  !-----------------------------------------------------------------------------------------------
  ! subroutine add_cmbandbeam
  !-----------------------------------------------------------------------------------------------

  subroutine add_cmbandbeam(map, power, beam, nside, addbeam, addcmb)
    implicit none
    real(dp),                 dimension(0:,:)   :: map, power, beam
    integer(i4b)                                :: nside
    logical(lgt)                                :: addbeam, addcmb

    integer(i4b)                                :: lmax, l, m, n, i,nmaps
    real(dp)                                    :: t1, t2, eta, x, y, matrix, z(2)
    complex(dpc), allocatable, dimension(:,:,:) :: alm
    real(dp),     pointer,     dimension(:,:)   :: dw8

    ! Finding alm's for given skymap
    call cpu_time(t1)
    nmaps = 1
    if (.not. addbeam) then
       lmax = size(power(:,1)) - 1
    write(*,*) lmax, '= lmax power'
    else if  (.not. addcmb) then
       lmax = size(beam(:,1)) - 1
    write(*,*) lmax, '= lmax beam'
    else
       lmax = min(size(power(:,1)),size(beam(:,1))) - 1
    write(*,*) lmax, '= lmax min'
    end if
    write(*,*) lmax, '= lmax'
    allocate(alm(1:nmaps, 0:lmax, 0:lmax))
    allocate(dw8(1:2*nside, 1:nmaps))
    call read_ringweights(nside, .false., dw8)
    z = 0.d0
    if (nmaps==1) then
       call map2alm(nside, lmax, lmax, map(:,1), alm, z, dw8)
    else
       call map2alm(nside, lmax, lmax, map, alm, z, dw8)
    end if
    deallocate(dw8)
    call cpu_time(t2)
    if (myid==0) write(*,*) 'Time map2alm:', t2-t1 

    ! adding cmb and beam
    call cpu_time(t1)
    do l = 0, lmax
! write(*,*) l, lmax, 'l'
       if (addcmb) then
          matrix = power(l,1)
          if (matrix==0.d0) then
             alm(:,l,:) = alm(:,l,:) + cmplx(0,0)
          else
             do m = 0, l
                ! Generate gaussian noise
                eta = rand_gauss(rng_handle_cmb) 
                ! .. with the right cov matrix
                x = eta*sqrt(matrix)
                if (m==0) then
                   alm(:,l,m) = alm(:,l,m) +  cmplx(x,0)
                else
                   eta = rand_gauss(rng_handle_cmb) 
                   y = eta*sqrt(matrix)
                   alm(:,l,m) = alm(:,l,m) + cmplx(x,y)/sqrt(2.d0)
                end if
             end do
          end if
       end if
       if (addbeam) alm(:,l,:) = alm(:,l,:)*beam(l,1)
    end do
    call cpu_time(t2)
    if (myid==0) write(*,*) 'Time almsim:', t2-t1 

    call cpu_time(t1)
    call alm2map(nside, lmax, lmax, alm, map(:,1))
    call cpu_time(t2)
    if (myid==0) write(*,*) 'Time alm2map:', t2-t1 
    deallocate(alm)

  end subroutine add_cmbandbeam

  !-----------------------------------------------------------------------------------------------
  ! subroutine transfer
  !-----------------------------------------------------------------------------------------------

  subroutine transfer(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: inmapfile, outmapfile, outprefix, outfile, docomponent
    character(len=5)            :: nsidein, nsideout
    integer(i4b)                :: ordering, ordering2, nside, nside2, nmaps, nmaps2, npix, npix2
    integer(i4b)                :: i, l, m, lmax, pol, poldeg, comp, hakk
    real(dp)                    :: z(2), sum

    real(dp),     allocatable, dimension(:)     :: inpower, outpower, x, y, a
    real(dp),     allocatable, dimension(:,:)   :: inmap, outmap, mask, maskin
    real(dp),     pointer,     dimension(:,:)   :: dw8in, dw8out
    complex(dpc), allocatable, dimension(:,:,:) :: almin, almout
    character(len=2),          dimension(3)     :: compname

    ! Get parameters
    if (iargc() /= 5) then
       write(*,*) 'transfer takes 4 parameters'
       call give_user_info
    else 
       call getarg(2, inmapfile)
       call getarg(3, outmapfile)
       call getarg(4, docomponent)
       call getarg(5, outprefix)
    end if

    ! Read maps
    call read_map(inmap, ordering, inmapfile, nside=nside, nmap=nmaps)
    call read_map(outmap, ordering2, outmapfile, nside=nside2, nmap=nmaps2)
    if (nside /= nside2) then
       if (myid==0) then 
          write(*,*) nside, '= nside for inmap'
          write(*,*) nside2, '= nside for outmap'
          write(*,*) 'Will estimate transfer function included pixwin'
       end if
    end if
    if (myid==0) then 
       write(*,*) nmaps, '= nmaps for inmap'
       write(*,*) nmaps2, '= nmaps for outmap'
    end if

    if (trim(docomponent)=='TT') then
       comp = 1
       if (myid==0) write(*,*) 'Making transfer function for TT' 
    else if (trim(docomponent)=='EE') then
       comp = 2
       write(*,*) comp
       if (myid==0) write(*,*) 'Making transfer function for EE' 
    else if (trim(docomponent)=='BB') then
       comp = 3
       if (myid==0) write(*,*) 'Making transfer function for BB' 
    else
       if (myid==0) write(*,*) 'Please specify either TT, EE or BB' 
       call give_user_info()
    end if
    
    if (comp > nmaps) then
       if (myid==0) write(*,*) 'Not enough map components to do polarization. Aborting...'
       call give_user_info()
    end if

    ! Use nmaps according to input. pol = 1 if tt, else pol=2
    if (comp==1) then
       pol = 1
    else
       pol = 2
    end if

    npix = 12*nside**2
    npix2 = 12*nside2**2
    if (ordering /= 1) then
       if (myid==0) write(*,*) 'Converting inmap to ringed'
       call convert_nest2ring(nside, inmap)
       ordering = 1
    end if
    if (ordering2 /= 1) then
       if (myid==0) write(*,*) 'Converting outmap to ringed'
       call convert_nest2ring(nside2, outmap)
       ordering2 = 1
    end if

    ! Put grey pixels to zero before spherical harmonic dec
    ! Assuming all comps have same mask (we always do that anyway)
    if (nside == nside2) then
       where (abs(-1.6375d30 - inmap(:,pol)) < 1d25 .or. abs(-1.6375d30 - outmap(:,pol)) < 1d25)
          ! This is a bit weird depensing on the value of pol, but it works.
          inmap(:,pol)    = 0.d0
          outmap(:,pol)   = 0.d0
          inmap(:,nmaps)  = 0.d0
          outmap(:,nmaps) = 0.d0
       end where
    else
       allocate(mask(0:npix2-1,nmaps))
       allocate(maskin(0:npix-1,nmaps))
       mask = 1.d0

       ! Write map
       outfile = trim(outprefix) // '_outmap_before.fits'
       call write_map(outmap, ordering2,trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)

       where (abs(-1.6375d30 - outmap(:,pol)) < 1.d0)
          mask(:,pol)     = 0.d0
          outmap(:,pol)   = 0.d0
          mask(:,nmaps)   = 0.d0
          outmap(:,nmaps) = 0.d0
       end where

       ! Write map after zeroing the masked pixels
       outfile = trim(outprefix) // '_outmap.fits'
       call write_map(outmap, ordering2,trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Write mask
       outfile = trim(outprefix) // '_outmask.fits'
       call write_map(mask, ordering2,trim(outfile))
       write(*,*) '* Mask written to file = ', trim(outfile)
       
       write(*,*) count(mask/=0.d0), 'out'

       ! Using the same mask on the input map (upgrading to fit nside)
       call udgrade_ring(mask, nside2, maskin, nside) 

       write(*,*) count(mask/=0.d0), count(maskin/=0.d0), 'out, in', count(maskin/=0.d0)/count(mask/=0.d0)
       write(*,*) npix2, npix, 'out, in (npix)', npix/npix2

       ! Write mask
       outfile = trim(outprefix) // '_inmask.fits'
       call write_map(maskin, ordering,trim(outfile))
       write(*,*) '* Mask written to file = ', trim(outfile)
       ! Write map
       outfile = trim(outprefix) // '_inmap_before.fits'
       call write_map(inmap, ordering,trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)

       where (abs(0.d0 - maskin(:,pol)) < 1d-25) 
          inmap(:,pol)   = 0.d0
          inmap(:,nmaps) = 0.d0
       end where
       ! Write map
       outfile = trim(outprefix) // '_inmap.fits'
       call write_map(inmap, ordering,trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       
       deallocate(mask, maskin)
    end if

    ! Spherical harmonic decomposition
    lmax = 3*min(nside, nside2)
    allocate(almin(1:nmaps, 0:lmax, 0:lmax))
    allocate(almout(1:nmaps, 0:lmax, 0:lmax))
    z = 0.d0
    if (nmaps==1) then
       call read_ringweights(nside, .false., dw8in)
       call map2alm(nside, lmax, lmax, inmap(:,1), almin, z, dw8in)
    else
       call read_ringweights(nside, .true., dw8in)
       call map2alm(nside, lmax, lmax, inmap, almin, z, dw8in)
    end if
    z = 0.d0
    if (nmaps==1) then
       call read_ringweights(nside2, .false., dw8out)
       call map2alm(nside2, lmax, lmax, outmap(:,1), almout, z, dw8out)
    else
       call read_ringweights(nside2, .true., dw8out)
       call map2alm(nside2, lmax, lmax, outmap, almout, z, dw8out)
    end if
    deallocate(dw8in, dw8out, inmap, outmap)

    ! Find power spectrum
    allocate(inpower(0:lmax))
    allocate(outpower(0:lmax))
    inpower = 0.d0
    outpower = 0.d0
    do i = 0, lmax
       do m = 0, i
         if (m == 0) then
             inpower(i)  = inpower(i)  + abs(almin(comp, i, m))**2             
             outpower(i) = outpower(i) + abs(almout(comp, i, m))**2             
          else
             inpower(i)  = inpower(i)  + 2*abs(almin(comp, i, m))**2
             outpower(i) = outpower(i) + 2*abs(almout(comp, i, m))**2
          end if
       end do
       inpower(i)  = inpower(i)/(2.d0*i+1.d0)*i*(i+1.d0)/2.d0*pi
       outpower(i) = outpower(i)/(2.d0*i+1.d0)*i*(i+1.d0)/2.d0*pi
    end do
    deallocate(almin, almout)
    outfile = trim(outprefix) // '_power_spec1.dat'
    open(13, file=trim(outfile))
    do l = 2, lmax
       write(13,*) l, inpower(l)
    end do
    close(13)
    outfile = trim(outprefix) // '_power_spec2.dat'
    open(13, file=trim(outfile))
    do l = 2, lmax
       write(13,*) l, outpower(l)
    end do
    close(13)

    ! Find raw transfer function
    if (nside == nside2) then
       call int2string(nside2, nsideout)
       outfile = trim(outprefix) // '_raw_transfer_func_'//nsideout //'.dat'
    else
       call int2string(nside, nsidein)
       call int2string(nside2, nsideout)
       outfile = trim(outprefix)//'_raw_transfer_with_pixwin_'//nsidein//'_to_'//nsideout//'.dat'
    end if
    open(13, file=trim(outfile))
    do l = 2, lmax
       write(13,*) l, sqrt(outpower(l)/inpower(l))
    end do
    close(13)
    write(*,*) '* Raw transfer function written to file = ', trim(outfile)

    ! Smooth transfer function - skip the first 50 ells
    hakk = 50
    poldeg = 5
    allocate(a(0:poldeg))
    allocate(x(1:lmax), y(1:lmax))
    do i = 1, hakk
       x(i) = i
       y(i) = 1.d0
    end do
    do i = hakk + 1, lmax
       x(i) = i
       y(i) = sqrt(outpower(i)/inpower(i))
    end do
    call fit_polynomial(x, y, a)
    if (myid==0) then
       write(*,*) 'Polynomial coefficients:'
       write(*,*) a
    end if
    if (nside == nside2) then
       call int2string(nside2, nsideout)
       outfile = trim(outprefix) // '_transfer_func_'//nsideout //'.dat'
    else
       call int2string(nside, nsidein)
       call int2string(nside2, nsideout)
       outfile = trim(outprefix)//'_transfer_with_pixwin_'//nsidein//'_to_'//nsideout//'.dat'
    end if
    open(13, file=trim(outfile))
    do l = 1, hakk/2
       write(13,*) l, 1.d0
    end do
    do l = hakk/2+1, lmax
       sum = a(0)
       do i = 1, poldeg
          sum = sum + a(i)*real(l,dp)**i
       end do
       write(13,*) l, sum
    end do
    close(13)
    write(*,*) '* Smoothed transfer function written to file = ', trim(outfile)
   

    ! Clean up
    deallocate(inpower, outpower, a, x, y)

  end subroutine transfer


  !-----------------------------------------------------------------------------------------------
  ! subroutine smooth_transfer
  !-----------------------------------------------------------------------------------------------

  subroutine smooth_transfer(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: infile, outprefix, outfile
    integer(i4b)                :: i, l, m, lmax, pol, num, poldeg, comp, hakk
    real(dp)                    :: verdi

    real(dp),     allocatable, dimension(:)     :: x, y, a
    real(dp),     allocatable, dimension(:,:)   :: data


    ! Get parameters
    if (iargc() /= 3) then
       write(*,*) 'smooth takes 2 parameters'
       call give_user_info
    else 
       call getarg(2, infile)
       call getarg(3, outprefix)
    end if

    ! Read txtfile
    lmax = 1150
    allocate(data(0:lmax,2))
    allocate(x(0:lmax), y(0:lmax))
    open(13, file=trim(infile))
    do l = 2, lmax
       read(13,*) data(l,:)
    end do
    x(0) = 0.d0
    y(0) = 1.d0
    x(1) = 1.d0
    y(1) = 1.d0
    do l = 2, lmax
       x(l) = data(l,1)
       y(l) = data(l,2)
    end do

    ! smooth transfer function
    hakk = 25
    do l = 100, lmax-hakk
       data(l,2) = sum(y(l-hakk:l+hakk))/(2.d0*real(hakk,dp)+1.d0)
    end do
write(*,*) 'a1'
    deallocate(x,y)
write(*,*) 'a2'
    ! first bit
    poldeg = 2
    allocate(a(0:poldeg))
    allocate(x(1:3), y(1:3))
    x(1) = 0.d0;   y(1) = 1.d0
    x(2) = 1.d0;   y(2) = 1.d0
    x(3) = 100.d0; y(3) = data(100,2)
write(*,*) 'a3'
    call fit_polynomial(x, y, a)
    write(*,*) a

    ! output
    outfile = trim(outprefix) // '_smoothed_transfer_func.dat'
    open(13, file=trim(outfile))
    do l = 0, 99
       verdi = a(0)
       do i = 1, poldeg
          verdi = verdi + a(i)*real(l,dp)**i
       end do
       write(13,*) l, verdi
    end do
    do l = 100, lmax-hakk
       write(13,*) l, data(l,2)
    end do
    do l = lmax-hakk+1, 1200
       write(13,*) l, data(lmax-hakk,2)
    end do
    close(13)
    write(*,*) '* Smoothed transfer function written to file = ', trim(outfile)
    
    deallocate(x,y,a, data)
    
  end subroutine smooth_transfer

  !-----------------------------------------------------------------------------------------------
  ! subroutine leakage
  !-----------------------------------------------------------------------------------------------

  subroutine leakage(unit)
    implicit none
    
    integer(i4b),       intent(in) :: unit
    
    character(len=256)             :: tempfile, listfile, infreq1, infreq2, outprefix, outfile, rmsfile1, rmsfile2
    character(len=256)             :: beamfile1, beamfile2, maskfile, infwhm, innside, inradius, incutradius, insindex
    integer(i4b)                   :: i, j, k, n, comp, pol, nlist, p
    integer(i4b)                   :: nmaps, nside, npix, ordering, ierr
    integer(i4b)                   :: numdiodes
    integer(i4b)                   :: areanside, areanpix, num
    real(dp)                       :: nullval, freq1, freq2, fwhm, a2t, radius, cutradius, freq_corr, specindex
    logical(lgt)                   :: anynull, mismatch, inv
    real(dp)                       :: healnan=-1.6375d30, root=0
  
    real(dp),       pointer,     dimension(:,:)   :: tmap, map2, mask, rms1, rms2!, beam1, beam2
    real(dp),       allocatable, dimension(:)     :: leak, uncert
!    integer(i4b),   allocatable, dimension(:)     :: map2heal, antimap2heal
    character(len=256), pointer, dimension(:)     :: maplist
    character(len=3)                              :: filnum
    real(dp),                    dimension(3)     :: gc_vec
    integer(i4b),   allocatable, dimension(: )    :: listpix

    ! Get parameters
    if (iargc() /= 12) then
       write(*,*) 'leakage takes 11 parameters'
       call give_user_info
    else 
       call getarg(2, tempfile)
       call getarg(3, listfile)
       call getarg(4, beamfile1)
       call getarg(5, beamfile2)
       call getarg(6, infreq1)
       call getarg(7, infreq2)
       call getarg(8, insindex)
       call getarg(9, infwhm)
       call getarg(10, inradius)
       call getarg(11, incutradius)
!       call getarg(10, rmsfile1)
!       call getarg(11, rmsfile2)
       call getarg(12, outprefix)
    end if
    
    read(infreq1,*) freq1
    read(infreq2,*) freq2
    read(infwhm,*)  fwhm
    read(inradius,*)  radius
    read(incutradius,*)  cutradius
    read(insindex,*)  specindex
    if (myid==root) write(*,*) 'Estimating leakage from ', trim(infreq1), ' GHz map and ', trim(infreq2), ' GHz map'
    if (myid==root) write(*,*) 'using area with radius ', trim(inradius), ' degrees.'
    write(*,*) 'freq1', freq1
    write(*,*) 'freq2', freq2
    write(*,*) 'fwhm ', fwhm
    write(*,*) 'radius ', radius
    write(*,*) 'cutradius ', cutradius
    write(*,*) 'spectral index ', specindex
    cutradius = (180-cutradius)*pi/180.d0   ! radius in radians
    gc_vec = [-1,0,0]                       ! because we want to remove what's OUTSIDE the disk
    freq_corr= (freq1/freq2)**specindex        
    write(*,*) 'freq correction ', freq_corr

    a2t = ant2thermo(freq2)
   
    ! Read maps and mask and convert to nest if necassary
    call read_mapfile(tempfile,  npix, nside, nmaps, ordering, tmap,  nest=.true.)
!    call read_mapfile(maskfile,  npix, nside, nmaps, ordering, mask,  nest=.true., check_nside=nside)
!    call read_mapfile(rmsfile1,  npix, nside, nmaps, ordering, rms1,  nest=.true., check_nside=nside)
!    call read_mapfile(rmsfile2,  npix, nside, nmaps, ordering, rms2,  nest=.true., check_nside=nside)

    ! Always use temperature data
    pol = 1
    nmaps = 1

    ! Write maps
    outfile = trim(outprefix) // '_input_tmap.fits'
    call write_map(tmap, ordering,trim(outfile))
    write(*,*) '* Map written to file = ', trim(outfile)

    ! Put missing pixels to zero
    where (abs(-1.6375d30 - tmap) < 1d25) tmap=0.d0
    ! Convolve temperature map with common gaussian beam
    if (myid==root) write(*,fmt='(a,f5.1,a)') 'Smoothing to common gaussian beam size of ',fwhm,' arc minutes.'
    call beam_convolve(tmap, ordering, nside, nmaps, beamfile1, fwhm)
    ! Restore unhit pixels
    where (tmap == 0.d0) tmap = healnan

    ! Write maps
    outfile = trim(outprefix) // '_beam_convolved_tmap.fits'
    call write_map(tmap, ordering, trim(outfile))
    write(*,*) '* Map written to file = ', trim(outfile)

    ! Read filelist
    call read_filelist(unit, listfile, numdiodes, maplist)
    allocate(leak(numdiodes))
    allocate(uncert(numdiodes))
    allocate(listpix(0:npix-1))
    ! Loop over files
    do i = 1+myid, numdiodes, numprocs
       write(*,*) 'Processing file no.', i, 'out of',numdiodes
       call read_mapfile(trim(maplist(i)), npix, nside, nmaps, ordering, map2, nest=.true., check_nside=nside)
       where (abs(healnan - map2) > 1d25) map2=map2*a2t*1000.d0
       ! Write map
       call int2string(i, filnum)
       outfile = trim(outprefix) // '_input_map' // filnum // '.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Find pixels outside a cutradius deg radius around gc center
       call query_disc(nside, gc_vec, cutradius, listpix, nlist, nest=1)
       do p=0,nlist-1
          map2(listpix(p), :)=0.d0
       end do
        ! Write map
       call int2string(i, filnum)
       outfile = trim(outprefix) // '_input_map' // filnum // '_cut.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Beam convolve
       where (abs(healnan - map2) < 1d25) map2=0.d0
       call beam_convolve(map2, ordering, nside, 1, beamfile2, fwhm)
       where (map2 == 0.d0) map2 = healnan
       ! Write map
       outfile = trim(outprefix) // '_beam_convolved_map' // filnum // '.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       call estimate_leakage(tmap(:,1), map2(:,1), nside, radius, leak(i), uncert(i))
       write(*,*) i, 'leakage=', leak(i)*100.d0
       write(*,*)
    end do

    ! Write leakage coeffisients to file
    outfile = trim(outprefix) // '_leakage.txt'
    open(unit, file=trim(outfile))
    do i = 1, numdiodes
       write(unit,*) 'i2qu_leak  0  0  -1  -1 ', (i-1)/4, modulo(i-1,4), leak(i)*100.d0*freq_corr
       write(*,*) 'i2qu_leak  0  0  -1  -1 ', (i-1)/4, modulo(i-1,4), leak(i)*100.d0*freq_corr
    end do    
    close(unit)
    write(*,*) '* Leakage coeffisients written to file = ', trim(outfile)

    ! Clean up
    deallocate(leak)
    deallocate(uncert)
    deallocate(listpix)

  end subroutine leakage

  !-----------------------------------------------------------------------------------------------
  ! subroutine spectral_index2
  !-----------------------------------------------------------------------------------------------

  subroutine spectral_index2(unit)
    implicit none
    
    integer(i4b),       intent(in) :: unit
    
    character(len=256)             :: infile1, infile2, infreq1, infreq2, outprefix, outfile, rmsfile1, rmsfile2
    character(len=256)             :: beamfile1, beamfile2, maskfile, infwhm, innside, inuncert
    integer(i4b)                   :: i, j, k, n, p, comp, pol, maskcomp, lstop, nband=2
    integer(i4b)                   :: nmaps, nside, npix, ordering, ierr
    integer(i4b)                   :: nmaps_mask, nmaps_rms
    integer(i4b)                   :: areanside, areanpix, num, max_nside, ns, ns_max, max_npix
    real(dp)                       :: nullval, freq1, freq2, fwhm, uncertlimit
    logical(lgt)                   :: anynull, mismatch, inv, multi
    real(dp)                       :: healnan=-1.6375d30, root=0
  
    real(dp),     pointer,     dimension(:,:)   :: map1, map2, rms1, rms2
    real(dp),     allocatable, dimension(:,:)   :: outmap, areamap, uncert, multimap, multiuncert
    real(dp),     allocatable, dimension(:,:)   :: samlekart, samlekart2, samleuncert, samleuncert2 
    real(dp),     allocatable, dimension(:,:)   :: offset, fg, samleoffset, samlefg
    real(dp),     allocatable, dimension(:)     :: ratio, kart1, kart2, var1, var2
    real(dp),     allocatable, dimension(:,:,:) :: kart
    integer(i4b), allocatable, dimension(:)     :: map2heal, antimap2heal, locmask
    integer(i4b), pointer,     dimension(:)     :: mask


    ! Get parameters
    if (iargc() /= 13 .and. iargc() /= 14) then
       write(*,*) 'index takes 12(13) parameters'
       call give_user_info
    else 
       call getarg(2, infile1)
       call getarg(3, infile2)
       call getarg(4, beamfile1)
       call getarg(5, beamfile2)
       call getarg(6, infreq1)
       call getarg(7, infreq2)
       call getarg(8, infwhm)
       call getarg(9, maskfile)
       call getarg(10, rmsfile1)
       call getarg(11, rmsfile2)
       call getarg(12, innside)
       call getarg(13, outprefix)
       if (iargc()==14) then
          multi = .true.
          call getarg(14, inuncert)
          read(inuncert,*) uncertlimit
          write(*,*) 'Multiresolution: uncertlimit =', uncertlimit
       else
          multi = .false.
       end if
    end if

    read(infreq1,*) freq1
    read(infreq2,*) freq2
    read(infwhm,*)  fwhm
    read(innside,*) max_nside
    ns_max = nint(log(real(max_nside,dp))/log(2.d0))+1
    if (max_nside==0) ns_max=0
    if (2**(ns_max-1) /= max_nside) then
       write(*,*) 'nside =',max_nside,' not a factor of 2. Quitng'
       stop
    end if
    if (myid==root) write(*,*) 'Finding spectral index from ', trim(infreq1), ' GHz map and ', trim(infreq2), ' GHz map'
    if (myid==0) then 
       write(*,*) 'freq1', freq1
       write(*,*) 'freq2', freq2
       write(*,*) 'fwhm', fwhm
       write(*,*) 'nside large pixels', max_nside
    end if
  
    ! Read maps and convert to nest if necassary
    call read_mapfile(infile1,  npix, nside, nmaps, ordering, map1, nest=.true.)
    call read_mapfile(infile2,  npix, nside, nmaps, ordering, map2, nest=.true., check_nside=nside, check_nmaps=nmaps)
    ! Checking polariation
    if (nmaps==1) then
       pol = 1
       if (myid==0) write(*,*) 'Running at temperature data' 
    else if (nmaps==3) then
       pol = 2
       if (myid==0) write(*,*) 'Running at polarisation data' 
    else
       write(*,*) nmaps, '= nmaps. Unknown number. Quiting'
    end if
pol=1
nmaps=1
write(*,*) 'OBS OBS hardcoded to run at temperature data!!!!'
    ! Read mask and convert to nest if necessary
    call read_and_check_mask(myid, maskfile, mask, nmaps, ordering, npix, pol)
    ! Read rms maps and convert to nest if necassary
 if (rmsfile1 /='none') then
    call read_mapfile(rmsfile1, npix, nside, nmaps_rms,  ordering, rms1, nest=.true., check_nside=nside)
    if (nmaps_rms<nmaps) then
       if (myid==root) write(*,*) 'We dont have smoothed rms for right component. Quiting'
       call mpi_finalize(ierr)
       stop
    end if    
    call read_mapfile(rmsfile2, npix, nside, nmaps_rms,  ordering, rms2, nest=.true., check_nside=nside)
    if (nmaps_rms<nmaps) then
       if (myid==root) write(*,*) 'We dont have smoothed rms for right component. Quiting'
       call mpi_finalize(ierr)
       stop
    end if    
end if

    ! Write input maps
    if (myid==0) then 
       outfile = trim(outprefix) // '_input_ymap.fits'
       call write_map(map1, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_input_xmap.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! Check for common pixels
    n = count(healok(map1(:, nmaps)) .and. healok(map2(:, nmaps)) .and. mask == 1)
    if (myid==root) write(*, *) 'The two input maps have ', n, 'common pixels'
    if (n==0) return

    ! Save map2heal
    j=0
    k=0
    allocate(map2heal(n))
    allocate(antimap2heal(npix-n))
    if (pol == 1) then
       do i = 0, npix-1
          if (healok(map1(i, nmaps)) .and. healok(map2(i, nmaps)) .and. mask(i)==1) then
             j = j+1
             map2heal(j) = i
          else
             k = k+1
             antimap2heal(k) = i
             mask(i) = 0
          end if
       end do
    else if (pol==2) then
       do i = 0, npix-1
          if (healok(map1(i, 2)) .and. healok(map2(i, 2)) .and. healok(map1(i, 3)) .and. healok(map2(i, 3)) .and. mask(i)==1) then     !!!!OBS only for pol
             j = j+1
             map2heal(j) = i
          else
             k = k+1
             antimap2heal(k) = i
             mask(i) = 0
          end if
       end do
    end if
    ! Put missing pixels to zero
    do j = 1, nmaps
       do i = 0, npix-1
          if (healnot(map1(i, j))) map1(i, j) = 0.d0
          if (healnot(map2(i, j))) map2(i, j) = 0.d0
       end do
    end do

!!$    ! Write maps
!!$    outfile = trim(outprefix) // '_input_map1_after.fits'
!!$    call write_map(map1, ordering, trim(outfile))
!!$    write(*,*) '* Map written to file = ', trim(outfile)
!!$    outfile = trim(outprefix) // '_input_map2_after.fits'
!!$    call write_map(map2, ordering, trim(outfile))
!!$    write(*,*) '* Map written to file = ', trim(outfile)

    ! Convolve both maps with common gaussian beam (skipping if filename is 'none')
    if (beamfile1 =='none') then
       if (myid==root) write(*,fmt='(a)') 'NB: No smoothing of beams since beam filename1 is none.'
    else
       if (myid==root) write(*,fmt='(a,f5.1,a)') 'Smoothing to common gaussian beam size of ',fwhm,' arc minutes.'
       call beam_convolve(map1, ordering, nside, nmaps, beamfile1, fwhm, onlypol=.true.)
       call beam_convolve(map2, ordering, nside, nmaps, beamfile2, fwhm, onlypol=.true.)
    end if

    ! Restoring missing pixels
    map1(antimap2heal,:) = healnan
    map2(antimap2heal,:) = healnan
    deallocate(antimap2heal)

    ! Write beam_convolved maps
    if (myid==0 .and. beamfile1 /= 'none') then 
       outfile = trim(outprefix) // '_beam_conv_map1.fits'
       call write_map(map1, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_beam_conv_map2.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! loop over nsides of large pixels
    do ns = 0, 0!ns_max
       if (all(mask==0)) exit
       if (allocated(kart1)) deallocate(kart1)
       if (allocated(kart2)) deallocate(kart2)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(areamap)) deallocate(areamap)
       if (allocated(offset)) deallocate(offset)
       ! Divide into healpix areas
       areanside=max_nside/(2**ns)
       if (myid==root) write(*,fmt='(a,i3,a)') 'Scatterplot: Choosing regions to be nside =',areanside,' healpix pixels.'
       areanpix = 12*areanside**2
if (max_nside==0) areanpix=1
       num = (npix/areanpix)   ! pixels per large healpix
!       num = (nside/areanside)**2   ! pixels per large healpix
       allocate(areamap(0:areanpix-1, nmaps))
       allocate(offset(0:areanpix-1, nmaps))
       areamap = healnan
       offset  = healnan
       allocate(kart1(pol*num))
       allocate(kart2(pol*num))
       allocate(locmask(num))
       do i = 0+myid, areanpix-1, numprocs
          locmask  = mask(i*num:(i+1)*num-1)
          if (pol==1) then 
             kart1 = map1(i*num:(i+1)*num-1, 1)
             kart2 = map2(i*num:(i+1)*num-1, 1)
          else if (pol==2) then
             kart1(1:num)       = map1(i*num:(i+1)*num-1, 2)
             kart1(num+1:2*num) = map1(i*num:(i+1)*num-1, 3)
             kart2(1:num)       = map2(i*num:(i+1)*num-1, 2)
             kart2(num+1:2*num) = map2(i*num:(i+1)*num-1, 3)
          end if
          call scatterplot(kart1, kart2, locmask, pol, i, areamap(i, :), freq1, freq2, offset(i, :), outprefix)
          mask(i*num:(i+1)*num-1) = locmask
       end do
    end do
    ! Collect to one root map
    allocate(samlekart(0:areanpix-1, nmaps))
    allocate(samleoffset(0:areanpix-1, nmaps))
    where(areamap == healnan) areamap = 0.d0
    where(offset == healnan)  offset  = 0.d0
    call mpi_reduce(areamap, samlekart,   size(areamap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(offset,  samleoffset, size(areamap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
write(*,*) 'scatter areamap samlekart', size(areamap), size(samlekart)

    ! Write spectral index map and offsetainty to file
    if (myid==0) then 
       where(samlekart == 0.d0)   samlekart   = healnan
       where(samleoffset == 0.d0) samleoffset = healnan
       outfile = trim(outprefix) // '_spectral_index_from_scatter.fits'
       outfile = trim(outprefix) // '_slope_from_scatter.fits'
       call write_map(samlekart, ordering, trim(outfile))
       write(*,*) '* Spectral index map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_offset_from_scatter.fits'
       call write_map(samleoffset, ordering, trim(outfile))
       write(*,*) '* Offset in spectral index written to file = ', trim(outfile)
    end if
    
    ! Clean up
    if (allocated(kart1)) deallocate(kart1)
    if (allocated(kart2)) deallocate(kart2)
    if (allocated(areamap)) deallocate(areamap)
    if (allocated(offset)) deallocate(offset)
    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samleoffset)) deallocate(samleoffset)

return

    ! Do the same for maxlike
    ! loop over nsides of large pixels
    if (max_nside==0) then
       max_npix=1
    else
       max_npix=12*max_nside**2
    end if
    allocate(multimap(0:max_npix-1, nmaps))
    allocate(samlekart(0:max_npix-1, nmaps))
    allocate(samlekart2(0:max_npix-1, nmaps))
    allocate(multiuncert(0:max_npix-1, nmaps))
    allocate(samleuncert(0:max_npix-1, nmaps))
    allocate(samleuncert2(0:max_npix-1, nmaps))
    multimap=healnan
    samlekart2=healnan
    multiuncert=healnan
    samleuncert2=healnan
    ! loop over nsides of large pixels
    if (multi) then
       lstop = ns_max
    else
       lstop = 0
    end if
    do ns = 0, lstop
       if (all(mask==0)) exit
       if (allocated(kart)) deallocate(kart)
       if (allocated(var1)) deallocate(var1)
       if (allocated(var2)) deallocate(var2)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(areamap)) deallocate(areamap)
       if (allocated(uncert)) deallocate(uncert)
       if (allocated(offset)) deallocate(offset)
       if (allocated(fg)) deallocate(fg)
       ! Divide into healpix areas
       areanside=max_nside/(2**ns)
       if (myid==root) write(*,fmt='(a,i3,a)') 'Maxlike: Choosing regions to be nside =',areanside,' healpix pixels.'
       areanpix = 12*areanside**2
       if (areanside==0) areanpix=1
       num = (npix/areanpix)   ! pixels per large healpix
       allocate(areamap(0:areanpix-1, nmaps))
       allocate(uncert(0:areanpix-1, nmaps))
       allocate(offset(0:areanpix-1, nmaps))
       allocate(fg(0:npix-1, nmaps))
       allocate(kart(num, nmaps, nband))
       allocate(var1(pol*num))
       allocate(var2(pol*num))
       allocate(locmask(num))
       areamap = healnan
       uncert  = healnan
       offset  = healnan
       uncert  = healnan
       fg      = healnan
       multimap = healnan
       samlekart = healnan
       do i = 0+myid, areanpix-1, numprocs
          locmask  = mask(i*num:(i+1)*num-1)
          kart(:,:,1) = map1(i*num:(i+1)*num-1,:)
          kart(:,:,2) = map2(i*num:(i+1)*num-1,:)
          if (pol==1) then 
             var1  = rms1(i*num:(i+1)*num-1, 1)
             var2  = rms2(i*num:(i+1)*num-1, 1)
          else if (pol==2) then
             var1(1:num)        = rms1(i*num:(i+1)*num-1, 2)
             var1(num+1:2*num)  = rms1(i*num:(i+1)*num-1, 3)
             var2(1:num)        = rms2(i*num:(i+1)*num-1, 2)
             var2(num+1:2*num)  = rms2(i*num:(i+1)*num-1, 3)
          end if
          call maxlikeplot(kart, var1, var2, locmask, pol, i, freq1, freq2, areamap(i,:), &
               & offset(i,:), uncert(i,:), fg(i*num:(i+1)*num-1, :))
          if (myid==root .and. max_nside==0) write(*,*) 'Maxlike: spectral index =',areamap(i,1),' Uncert =',uncert(i,1)
          if (multi) then
             do j = i*(4**ns), i*(4**ns) + 4**ns -1
                if (healnot(multimap(j,1)) .and. abs(uncert(i,1))<uncertlimit) then 
                   multimap(j,:) = areamap(i,:)
                   multiuncert(j,:) = uncert(i,:)
                end if
             end do
             mask(i*num:(i+1)*num-1) = locmask 
          end if
       end do
       where(multimap == healnan) multimap = 0.d0
       where(multiuncert == healnan) multiuncert = 0.d0
       call mpi_allreduce(multimap,samlekart, size(multimap),mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(multiuncert,samleuncert, size(multiuncert),mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
       where(samlekart   == 0.d0) samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       where(samlekart2 == healnan) samlekart2 = samlekart         
       where(samleuncert2 == healnan) samleuncert2 = samleuncert
    end do

    ! Collect to one root map
    allocate(samleoffset(0:areanpix-1, nmaps))
    allocate(samlefg(0:npix-1, nmaps))
    where(areamap == healnan) areamap = 0.d0
    where(uncert == healnan)  uncert  = 0.d0
    where(offset == healnan)  offset  = 0.d0
    where(fg == healnan)      fg  = 0.d0
    if (.not. multi) then
       call mpi_reduce(areamap,samlekart2,   size(areamap),mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       call mpi_reduce(uncert,  samleuncert2, size(uncert),mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    end if

    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samleuncert)) deallocate(samleuncert)
    allocate(samlekart(0:npix-1, nmaps))
    allocate(samleuncert(0:npix-1, nmaps))
    samlekart = healnan
    samleuncert=healnan
    if (ns_max == 0) then
       ns = nint(log(real(nside,dp))/log(2.d0)) +3
    else
       ns = nint(log(real(nside,dp))/log(2.d0)) - nint(log(real(max_nside,dp))/log(2.d0))  
    end if
    ns = 4**ns
    do k = 1, nmaps
       do i = 0, npix-1
          if (mask(i)==1) then
             samlekart(i,k) = samlekart2(i/ns, k)
             samleuncert(i,k) = samleuncert2(i/ns, k)
          end if
       end do
    end do

    call mpi_reduce(offset,  samleoffset, size(offset),  mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(fg,      samlefg,     size(fg),      mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)


    ! Write maxlike spectral index map and uncertainty and offset to file
    if (myid==0) then 
       where(samlekart   == 0.d0) samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       where(samleoffset == 0.d0) samleoffset = healnan
       where(samlefg     == 0.d0) samlefg = healnan
       outfile = trim(outprefix) // '_spectral_index_from_maxlike.fits'
       call write_map(samlekart, ordering, trim(outfile))
       write(*,*) '* Spectral index map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_offset.fits'
       call write_map(samleoffset, ordering, trim(outfile))
       write(*,*) '* Offset map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_uncert.fits'
       call write_map(samleuncert, ordering, trim(outfile))
       write(*,*) '* Spectral index uncert written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_fg.fits'
       call write_map(samlefg, ordering, trim(outfile))
       write(*,*) '* foreground written to file = ', trim(outfile)
       ! Write histogram to file
       outfile = trim(outprefix) // '_histogram.txt'
       open(unit, file=trim(outfile))
       do i = 0, areanpix-1
          if (healok(samleuncert(i,1))) write(unit, *) (areamap(i,1)-(-2.7))/samleuncert(i,1)
       end do
       do j = 1, nmaps
          do i = 0, areanpix-1
             if (healok(samleuncert(i,j))) then 
                offset(i,j) = (samlekart(i,j)-(-2.7))/samleuncert(i,j)
             else
                offset(i,j) = healnan
             end if
          end do
       end do
       outfile = trim(outprefix) // '_delta.fits'
       call write_map(offset, ordering, trim(outfile))
       write(*,*) '* Delta spectral index written to file = ', trim(outfile)
    end if

    ! Clean up
    if (allocated(kart1)) deallocate(kart1)
    if (allocated(kart2)) deallocate(kart2)
    if (allocated(var1)) deallocate(var1)
    if (allocated(var2)) deallocate(var2)
    if (allocated(locmask)) deallocate(locmask)
    if (allocated(areamap)) deallocate(areamap)
    if (allocated(uncert)) deallocate(uncert)
    if (allocated(offset)) deallocate(offset)
    if (allocated(fg)) deallocate(fg)
    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samlekart2)) deallocate(samlekart2)
    if (allocated(samleuncert)) deallocate(samleuncert)
    if (allocated(samleuncert2)) deallocate(samleuncert2)
    if (allocated(samleoffset)) deallocate(samleoffset)
    if (allocated(samlefg)) deallocate(samlefg)
    if (allocated(multimap)) deallocate(multimap)
    if (allocated(multiuncert)) deallocate(multiuncert)


    ! Extract maps
    allocate(kart1(n))
    allocate(kart2(n))
    if (pol==1) then 
       kart1 = map1(map2heal, 1)
       kart2 = map2(map2heal, 1)
    else if (pol==2) then
       kart1 = sqrt(map1(map2heal, 2)**2 + map1(map2heal,3)**2)
       kart2 = sqrt(map2(map2heal, 2)**2 + map2(map2heal,3)**2)
    end if
 
    if (myid==0) then 
       ! Write scatter plot to file
       outfile = trim(outprefix) // '_scatter_plot.txt'
       open(unit, file=trim(outfile))
       if (pol==1) then 
          do i = 1, n
             write(unit,*) kart2(i), kart1(i)
          end do
       else if (pol==2) then
          do i = 1, n
             write(unit,*) map2(map2heal(i), 2),  map1(map2heal(i), 2)
          end do
          write(unit,*)
          do i = 1, n
             write(unit,*) map2(map2heal(i), 3),  map1(map2heal(i), 3)
          end do
       end if
       close(unit)
       write(*,*) 'Scatter plot written to file = ', trim(outfile)
    end if
    ! Starting to clean up
    deallocate(map1)
    deallocate(map2)
    
    ! Calculate naive ratio
    if (myid==0) then 
       allocate(ratio(n))
       ratio = kart1/kart2
       write(*,*) log(sum(ratio)/n)/log(freq1/freq2), '= spectral index from average ratio'
       ratio = log(ratio)/log(freq1/freq2)
       write(*,*) sum(ratio)/n, '= average spectral index'
       ! Writing spectral index map to file
       allocate(outmap(0:npix-1,nmaps))
       outmap = healnan
       outmap(map2heal,1) = ratio
       outfile = trim(outprefix) // '_spectral_index_from_ratio.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Writing ratio-based mask to file
       outmap = 0.d0
       do i = 1, n
          if (ratio(i)>-3.7d0 .and. ratio(i)<-2.3d0) then
             outmap(map2heal(i),:) = 1.d0
          end if
       end do
       outfile = trim(outprefix) // '_mask_from_ratio.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Mask from ratio written to file = ', trim(outfile)
       deallocate(ratio)    
       deallocate(outmap)
    end if

    ! Clean up
    deallocate(map2heal)
    deallocate(kart1)
    deallocate(kart2)
    deallocate(mask)
    deallocate(rms1)
    deallocate(rms2)

  end subroutine spectral_index2

  !-----------------------------------------------------------------------------------------------
  ! subroutine spectral_index
  !-----------------------------------------------------------------------------------------------

  subroutine spectral_index(unit)
    implicit none
    
    integer(i4b),       intent(in) :: unit
    
    character(len=256)             :: infile1, infile2, infreq1, infreq2, outprefix, outfile, rmsfile1, rmsfile2
    character(len=256)             :: beamfile1, beamfile2, maskfile, infwhm, innside, inuncert
    integer(i4b)                   :: i, j, k, n, p, comp, pol, maskcomp, lstop, nband=2
    integer(i4b)                   :: nmaps, nside, npix, ordering, ierr
    integer(i4b)                   :: nmaps_mask, nmaps_rms
    integer(i4b)                   :: areanside, areanpix, num, max_nside, ns, ns_max, max_npix
    real(dp)                       :: nullval, freq1, freq2, fwhm, uncertlimit
    logical(lgt)                   :: anynull, mismatch, inv, multi
    real(dp)                       :: healnan=-1.6375d30, root=0
  
    real(dp),     pointer,     dimension(:,:)   :: map1, map2, rms1, rms2
    real(dp),     allocatable, dimension(:,:)   :: outmap, areamap, uncert, multimap, multiuncert
    real(dp),     allocatable, dimension(:,:)   :: samlekart, samlekart2, samleuncert, samleuncert2 
    real(dp),     allocatable, dimension(:,:)   :: offset, fg, samleoffset, samlefg
    real(dp),     allocatable, dimension(:)     :: ratio, kart1, kart2, var1, var2
    real(dp),     allocatable, dimension(:,:,:) :: kart
    integer(i4b), allocatable, dimension(:)     :: map2heal, antimap2heal, locmask
    integer(i4b), pointer,     dimension(:)     :: mask


    ! Get parameters
    if (iargc() /= 13 .and. iargc() /= 14) then
       write(*,*) 'index takes 12(13) parameters'
       call give_user_info
    else 
       call getarg(2, infile1)
       call getarg(3, infile2)
       call getarg(4, beamfile1)
       call getarg(5, beamfile2)
       call getarg(6, infreq1)
       call getarg(7, infreq2)
       call getarg(8, infwhm)
       call getarg(9, maskfile)
       call getarg(10, rmsfile1)
       call getarg(11, rmsfile2)
       call getarg(12, innside)
       call getarg(13, outprefix)
       if (iargc()==14) then
          multi = .true.
          call getarg(14, inuncert)
          read(inuncert,*) uncertlimit
          write(*,*) 'Multiresolution: uncertlimit =', uncertlimit
       else
          multi = .false.
       end if
    end if

    read(infreq1,*) freq1
    read(infreq2,*) freq2
    read(infwhm,*)  fwhm
    read(innside,*) max_nside
    ns_max = nint(log(real(max_nside,dp))/log(2.d0))+1
    if (max_nside==0) ns_max=0
    if (2**(ns_max-1) /= max_nside) then
       write(*,*) 'nside =',max_nside,' not a factor of 2. Quitng'
       stop
    end if
    if (myid==root) write(*,*) 'Finding spectral index from ', trim(infreq1), ' GHz map and ', trim(infreq2), ' GHz map'
    if (myid==0) then 
       write(*,*) 'freq1', freq1
       write(*,*) 'freq2', freq2
       write(*,*) 'fwhm', fwhm
       write(*,*) 'nside large pixels', max_nside
    end if
  
    ! Read maps and convert to nest if necassary
    call read_mapfile(infile1,  npix, nside, nmaps, ordering, map1, nest=.true.)
    call read_mapfile(infile2,  npix, nside, nmaps, ordering, map2, nest=.true., check_nside=nside, check_nmaps=nmaps)
    ! Checking polariation
    if (nmaps==1) then
       pol = 1
       if (myid==0) write(*,*) 'Running at temperature data' 
    else if (nmaps==3) then
       pol = 2
       if (myid==0) write(*,*) 'Running at polarisation data' 
    else
       write(*,*) nmaps, '= nmaps. Unknown number. Quiting'
    end if
    ! Read mask and convert to nest if necassary
    call read_and_check_mask(myid, maskfile, mask, nmaps, ordering, npix, pol)
    ! Read rms maps and convert to nest if necassary
    call read_mapfile(rmsfile1, npix, nside, nmaps_rms,  ordering, rms1, nest=.true., check_nside=nside)
    if (nmaps_rms<nmaps) then
       if (myid==root) write(*,*) 'We dont have smoothed rms for right component. Quiting'
       call mpi_finalize(ierr)
       stop
    end if    
    call read_mapfile(rmsfile2, npix, nside, nmaps_rms,  ordering, rms2, nest=.true., check_nside=nside)
    if (nmaps_rms<nmaps) then
       if (myid==root) write(*,*) 'We dont have smoothed rms for right component. Quiting'
       call mpi_finalize(ierr)
       stop
    end if    

    ! Write input maps
    if (myid==0) then 
       outfile = trim(outprefix) // '_input_map1.fits'
       call write_map(map1, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_input_map2.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! Check for common pixels
    n = count(healok(map1(:, nmaps)) .and. healok(map2(:, nmaps)) .and. mask == 1)
    if (myid==root) write(*, *) 'The two input maps have ', n, 'common pixels'
    if (n==0) return

    ! Save map2heal
    j=0
    k=0
    allocate(map2heal(n))
    allocate(antimap2heal(npix-n))
    if (pol == 1) then
       do i = 0, npix-1
          if (healok(map1(i, nmaps)) .and. healok(map2(i, nmaps)) .and. mask(i)==1) then
             j = j+1
             map2heal(j) = i
          else
             k = k+1
             antimap2heal(k) = i
             mask(i) = 0
          end if
       end do
    else if (pol==2) then
       do i = 0, npix-1
          if (healok(map1(i, 2)) .and. healok(map2(i, 2)) .and. healok(map1(i, 3)) .and. healok(map2(i, 3)) .and. mask(i)==1) then     !!!!OBS only for pol
             j = j+1
             map2heal(j) = i
          else
             k = k+1
             antimap2heal(k) = i
             mask(i) = 0
          end if
       end do
    end if
    ! Put missing pixels to zero
    do j = 1, nmaps
       do i = 0, npix-1
          if (healnot(map1(i, j))) map1(i, j) = 0.d0
          if (healnot(map2(i, j))) map2(i, j) = 0.d0
       end do
    end do

!!$    ! Write maps
!!$    outfile = trim(outprefix) // '_input_map1_after.fits'
!!$    call write_map(map1, ordering, trim(outfile))
!!$    write(*,*) '* Map written to file = ', trim(outfile)
!!$    outfile = trim(outprefix) // '_input_map2_after.fits'
!!$    call write_map(map2, ordering, trim(outfile))
!!$    write(*,*) '* Map written to file = ', trim(outfile)

    ! Convolve both maps with common gaussian beam (skipping if filename is 'none')
    if (beamfile1 =='none') then
       if (myid==root) write(*,fmt='(a)') 'NB: No smoothing of beams since beam filename1 is none.'
    else
       if (myid==root) write(*,fmt='(a,f5.1,a)') 'Smoothing to common gaussian beam size of ',fwhm,' arc minutes.'
       call beam_convolve(map1, ordering, nside, nmaps, beamfile1, fwhm, onlypol=.true.)
       call beam_convolve(map2, ordering, nside, nmaps, beamfile2, fwhm, onlypol=.true.)
    end if

    ! Restoring missing pixels
    map1(antimap2heal,:) = healnan
    map2(antimap2heal,:) = healnan
    deallocate(antimap2heal)

    ! Write beam_convolved maps
    if (myid==0) then 
       outfile = trim(outprefix) // '_beam_conv_map1.fits'
       call write_map(map1, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_beam_conv_map2.fits'
       call write_map(map2, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! loop over nsides of large pixels
    do ns = 0, 0!ns_max
       if (all(mask==0)) exit
       if (allocated(kart1)) deallocate(kart1)
       if (allocated(kart2)) deallocate(kart2)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(areamap)) deallocate(areamap)
       if (allocated(uncert)) deallocate(uncert)
       ! Divide into healpix areas
       areanside=max_nside/(2**ns)
       if (myid==root) write(*,fmt='(a,i3,a)') 'Scatterplot: Choosing regions to be nside =',areanside,' healpix pixels.'
       areanpix = 12*areanside**2
if (max_nside==0) areanpix=1
       num = (npix/areanpix)   ! pixels per large healpix
!       num = (nside/areanside)**2   ! pixels per large healpix
       allocate(areamap(0:areanpix-1, nmaps))
       allocate(uncert(0:areanpix-1, nmaps))
       areamap = healnan
       uncert  = healnan
       allocate(kart1(pol*num))
       allocate(kart2(pol*num))
       allocate(locmask(num))
       do i = 0+myid, areanpix-1, numprocs
          locmask  = mask(i*num:(i+1)*num-1)
          if (pol==1) then 
             kart1 = map1(i*num:(i+1)*num-1, 1)
             kart2 = map2(i*num:(i+1)*num-1, 1)
          else if (pol==2) then
             kart1(1:num)       = map1(i*num:(i+1)*num-1, 2)
             kart1(num+1:2*num) = map1(i*num:(i+1)*num-1, 3)
             kart2(1:num)       = map2(i*num:(i+1)*num-1, 2)
             kart2(num+1:2*num) = map2(i*num:(i+1)*num-1, 3)
          end if
!          call scatterplot(kart1, kart2, locmask, pol, i, areamap(i, :), freq1, freq2, uncert(i, :) )
          mask(i*num:(i+1)*num-1) = locmask
       end do
    end do
    ! Collect to one root map
    allocate(samlekart(0:areanpix-1, nmaps))
    allocate(samleuncert(0:areanpix-1, nmaps))
    where(areamap == healnan) areamap = 0.d0
    where(uncert == healnan)  uncert  = 0.d0
    call mpi_reduce(areamap, samlekart,   size(areamap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(uncert,  samleuncert, size(areamap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
write(*,*) 'A3 scatter areamap samlekart', size(areamap), size(samlekart)

    ! Write spectral index map and uncertainty to file
    if (myid==0) then 
       where(samlekart == 0.d0)   samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       outfile = trim(outprefix) // '_spectral_index_from_scatter.fits'
       call write_map(samlekart, ordering, trim(outfile))
       write(*,*) '* Spectral index map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_spectral_index_uncertainty_from_scatter.fits'
       call write_map(samleuncert, ordering, trim(outfile))
       write(*,*) '* Uncertainty in spectral index written to file = ', trim(outfile)
    end if
    
    ! Clean up
    if (allocated(kart1)) deallocate(kart1)
    if (allocated(kart2)) deallocate(kart2)
    if (allocated(areamap)) deallocate(areamap)
    if (allocated(uncert)) deallocate(uncert)
    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samleuncert)) deallocate(samleuncert)

    ! Do the same for maxlike
    ! loop over nsides of large pixels
    if (max_nside==0) then
       max_npix=1
    else
       max_npix=12*max_nside**2
    end if
    allocate(multimap(0:max_npix-1, nmaps))
    allocate(samlekart(0:max_npix-1, nmaps))
    allocate(samlekart2(0:max_npix-1, nmaps))
    allocate(multiuncert(0:max_npix-1, nmaps))
    allocate(samleuncert(0:max_npix-1, nmaps))
    allocate(samleuncert2(0:max_npix-1, nmaps))
    multimap=healnan
    samlekart2=healnan
    multiuncert=healnan
    samleuncert2=healnan
    ! loop over nsides of large pixels
    if (multi) then
       lstop = ns_max
    else
       lstop = 0
    end if
    do ns = 0, lstop
       if (all(mask==0)) exit
       if (allocated(kart)) deallocate(kart)
       if (allocated(var1)) deallocate(var1)
       if (allocated(var2)) deallocate(var2)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(areamap)) deallocate(areamap)
       if (allocated(uncert)) deallocate(uncert)
       if (allocated(offset)) deallocate(offset)
       if (allocated(fg)) deallocate(fg)
       ! Divide into healpix areas
       areanside=max_nside/(2**ns)
       if (myid==root) write(*,fmt='(a,i3,a)') 'Maxlike: Choosing regions to be nside =',areanside,' healpix pixels.'
       areanpix = 12*areanside**2
       if (areanside==0) areanpix=1
       num = (npix/areanpix)   ! pixels per large healpix
       allocate(areamap(0:areanpix-1, nmaps))
       allocate(uncert(0:areanpix-1, nmaps))
       allocate(offset(0:areanpix-1, nmaps))
       allocate(fg(0:npix-1, nmaps))
       allocate(kart(num, nmaps, nband))
       allocate(var1(pol*num))
       allocate(var2(pol*num))
       allocate(locmask(num))
       areamap = healnan
       uncert  = healnan
       offset  = healnan
       uncert  = healnan
       fg      = healnan
       multimap = healnan
       samlekart = healnan
       do i = 0+myid, areanpix-1, numprocs
          locmask  = mask(i*num:(i+1)*num-1)
          kart(:,:,1) = map1(i*num:(i+1)*num-1,:)
          kart(:,:,2) = map2(i*num:(i+1)*num-1,:)
          if (pol==1) then 
             var1  = rms1(i*num:(i+1)*num-1, 1)
             var2  = rms2(i*num:(i+1)*num-1, 1)
          else if (pol==2) then
             var1(1:num)        = rms1(i*num:(i+1)*num-1, 2)
             var1(num+1:2*num)  = rms1(i*num:(i+1)*num-1, 3)
             var2(1:num)        = rms2(i*num:(i+1)*num-1, 2)
             var2(num+1:2*num)  = rms2(i*num:(i+1)*num-1, 3)
          end if
          call maxlikeplot(kart, var1, var2, locmask, pol, i, freq1, freq2, areamap(i,:), &
               & offset(i,:), uncert(i,:), fg(i*num:(i+1)*num-1, :))
          if (myid==root .and. max_nside==0) write(*,*) 'Maxlike: spectral index =',areamap(i,1),' Uncert =',uncert(i,1)
          if (multi) then
             do j = i*(4**ns), i*(4**ns) + 4**ns -1
                if (healnot(multimap(j,1)) .and. abs(uncert(i,1))<uncertlimit) then 
                   multimap(j,:) = areamap(i,:)
                   multiuncert(j,:) = uncert(i,:)
                end if
             end do
             mask(i*num:(i+1)*num-1) = locmask 
          end if
       end do
       where(multimap == healnan) multimap = 0.d0
       where(multiuncert == healnan) multiuncert = 0.d0
       call mpi_allreduce(multimap,samlekart, size(multimap),mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(multiuncert,samleuncert, size(multiuncert),mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
       where(samlekart   == 0.d0) samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       where(samlekart2 == healnan) samlekart2 = samlekart         
       where(samleuncert2 == healnan) samleuncert2 = samleuncert
    end do

    ! Collect to one root map
    allocate(samleoffset(0:areanpix-1, nmaps))
    allocate(samlefg(0:npix-1, nmaps))
    where(areamap == healnan) areamap = 0.d0
    where(uncert == healnan)  uncert  = 0.d0
    where(offset == healnan)  offset  = 0.d0
    where(fg == healnan)      fg  = 0.d0
    if (.not. multi) then
       call mpi_reduce(areamap,samlekart2,   size(areamap),mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       call mpi_reduce(uncert,  samleuncert2, size(uncert),mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    end if

    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samleuncert)) deallocate(samleuncert)
    allocate(samlekart(0:npix-1, nmaps))
    allocate(samleuncert(0:npix-1, nmaps))
    samlekart = healnan
    samleuncert=healnan
    if (ns_max == 0) then
       ns = nint(log(real(nside,dp))/log(2.d0)) +3
    else
       ns = nint(log(real(nside,dp))/log(2.d0)) - nint(log(real(max_nside,dp))/log(2.d0))  
    end if
    ns = 4**ns
    do k = 1, nmaps
       do i = 0, npix-1
          if (mask(i)==1) then
             samlekart(i,k) = samlekart2(i/ns, k)
             samleuncert(i,k) = samleuncert2(i/ns, k)
          end if
       end do
    end do

    call mpi_reduce(offset,  samleoffset, size(offset),  mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(fg,      samlefg,     size(fg),      mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)


    ! Write maxlike spectral index map and uncertainty and offset to file
    if (myid==0) then 
       where(samlekart   == 0.d0) samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       where(samleoffset == 0.d0) samleoffset = healnan
       where(samlefg     == 0.d0) samlefg = healnan
       outfile = trim(outprefix) // '_spectral_index_from_maxlike.fits'
       call write_map(samlekart, ordering, trim(outfile))
       write(*,*) '* Spectral index map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_offset.fits'
       call write_map(samleoffset, ordering, trim(outfile))
       write(*,*) '* Offset map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_uncert.fits'
       call write_map(samleuncert, ordering, trim(outfile))
       write(*,*) '* Spectral index uncert written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_fg.fits'
       call write_map(samlefg, ordering, trim(outfile))
       write(*,*) '* foreground written to file = ', trim(outfile)
       ! Write histogram to file
       outfile = trim(outprefix) // '_histogram.txt'
       open(unit, file=trim(outfile))
       do i = 0, areanpix-1
          if (healok(samleuncert(i,1))) write(unit, *) (areamap(i,1)-(-2.7))/samleuncert(i,1)
       end do
       do j = 1, nmaps
          do i = 0, areanpix-1
             if (healok(samleuncert(i,j))) then 
                offset(i,j) = (samlekart(i,j)-(-2.7))/samleuncert(i,j)
             else
                offset(i,j) = healnan
             end if
          end do
       end do
       outfile = trim(outprefix) // '_delta.fits'
       call write_map(offset, ordering, trim(outfile))
       write(*,*) '* Delta spectral index written to file = ', trim(outfile)
    end if

    ! Clean up
    if (allocated(kart1)) deallocate(kart1)
    if (allocated(kart2)) deallocate(kart2)
    if (allocated(var1)) deallocate(var1)
    if (allocated(var2)) deallocate(var2)
    if (allocated(locmask)) deallocate(locmask)
    if (allocated(areamap)) deallocate(areamap)
    if (allocated(uncert)) deallocate(uncert)
    if (allocated(offset)) deallocate(offset)
    if (allocated(fg)) deallocate(fg)
    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samlekart2)) deallocate(samlekart2)
    if (allocated(samleuncert)) deallocate(samleuncert)
    if (allocated(samleuncert2)) deallocate(samleuncert2)
    if (allocated(samleoffset)) deallocate(samleoffset)
    if (allocated(samlefg)) deallocate(samlefg)
    if (allocated(multimap)) deallocate(multimap)
    if (allocated(multiuncert)) deallocate(multiuncert)


    ! Extract maps
    allocate(kart1(n))
    allocate(kart2(n))
    if (pol==1) then 
       kart1 = map1(map2heal, 1)
       kart2 = map2(map2heal, 1)
    else if (pol==2) then
       kart1 = sqrt(map1(map2heal, 2)**2 + map1(map2heal,3)**2)
       kart2 = sqrt(map2(map2heal, 2)**2 + map2(map2heal,3)**2)
    end if
 
    if (myid==0) then 
       ! Write scatter plot to file
       outfile = trim(outprefix) // '_scatter_plot.txt'
       open(unit, file=trim(outfile))
       if (pol==1) then 
          do i = 1, n
             write(unit,*) kart2(i), kart1(i)
          end do
       else if (pol==2) then
          do i = 1, n
             write(unit,*) map2(map2heal(i), 2),  map1(map2heal(i), 2)
          end do
          write(unit,*)
          do i = 1, n
             write(unit,*) map2(map2heal(i), 3),  map1(map2heal(i), 3)
          end do
       end if
       close(unit)
       write(*,*) 'Scatter plot written to file = ', trim(outfile)
    end if
    ! Starting to clean up
    deallocate(map1)
    deallocate(map2)
    
    ! Calculate naive ratio
    if (myid==0) then 
       allocate(ratio(n))
       ratio = kart1/kart2
       write(*,*) log(sum(ratio)/n)/log(freq1/freq2), '= spectral index from average ratio'
       ratio = log(ratio)/log(freq1/freq2)
       write(*,*) sum(ratio)/n, '= average spectral index'
       ! Writing spectral index map to file
       allocate(outmap(0:npix-1,nmaps))
       outmap = healnan
       outmap(map2heal,1) = ratio
       outfile = trim(outprefix) // '_spectral_index_from_ratio.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Writing ratio-based mask to file
       outmap = 0.d0
       do i = 1, n
          if (ratio(i)>-3.7d0 .and. ratio(i)<-2.3d0) then
             outmap(map2heal(i),:) = 1.d0
          end if
       end do
       outfile = trim(outprefix) // '_mask_from_ratio.fits'
       call write_map(outmap, ordering, trim(outfile))
       write(*,*) '* Mask from ratio written to file = ', trim(outfile)
       deallocate(ratio)    
       deallocate(outmap)
    end if

    ! Clean up
    deallocate(map2heal)
    deallocate(kart1)
    deallocate(kart2)
    deallocate(mask)
    deallocate(rms1)
    deallocate(rms2)

  end subroutine spectral_index

  !-----------------------------------------------------------------------------------------------
  ! subroutine scatter_sindex
  !-----------------------------------------------------------------------------------------------

  subroutine scatter_sindex(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: map1file, map2file, outprefix, outfile, prefix
    character(len=512)          :: infreq1, infreq2, inlen, inbre, inrad
    character(len=5)            :: name
    character(len=2)            :: rad
    character(len=2)            :: comp(5)
    integer(i4b)                :: ordering2, nside2, nmaps, nmaps2, npix, n_alpha
    integer(i4b)                :: i, j, n, s, r, sources, numrad, ierr
    real(dp)                    :: freq1, freq2, len, bre, radius, vec(3), freq(2), spec(5), uncert(5)
    real(dp)                    :: maxv, minv, mid1, mid2, a, b, fwhm, elli, fwhm1, elli1, psi1
    real(dp)                    :: fwhm2, psi2, elli2, vinkel, fwhmT, elliT, psiT, beta, sigma, psi
    real(sp)                    :: fracs(2), angs(2), diff
    logical(lgt)                :: chatty

    real(dp),     allocatable, dimension(:)     :: p, rotmap1, rotmap2, fwhm_psi, ell_psi, dir_psi
    real(dp),     allocatable, dimension(:,:)   :: map1, map2, nmap1, nmap2, outmap, map4
    integer(i4b), allocatable, dimension(:)     :: healvec

    !real(dp), dimension(3) :: cvec
    
    comp(1) = 'T'
    comp(2) = 'Q'
    comp(3) = 'U'
    comp(4) = 'P'
    comp(5) = 'R'

    ! Get parameters
    if (iargc() /= 9) then
       write(*,*) 'scatter takes 8 parameters'
       call give_user_info
    else 
       call getarg(2, map1file)
       call getarg(3, infreq1)
       call getarg(4, map2file)
       call getarg(5, infreq2)
       call getarg(6, inlen)
       call getarg(7, inbre)
       call getarg(8, inrad)
       call getarg(9, outprefix)
    end if

    ! Read maps
    if (myid==0) write(*,*) 'Reading from ', trim(map1file)
    call read_map(map1, ordering,  map1file, nside=nside,  nmap=nmaps)
    if (myid==0) write(*,*) 'Reading from ', trim(map2file)
    call read_map(map2, ordering2, map2file, nside=nside2, nmap=nmaps2)
    call assert(ordering == ordering2, "Ordering mismatch between map1 and map2")
    call assert(nside    == nside2,    "Nside mismatch between map1 and map2")
    call assert(nmaps    == nmaps2,    "Nmaps mismatch between map1 and map2")
    call assert(nmaps    == 3,         "Need polarisation data!")
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,',ordering,'= ordering' 
    npix  = 12*nside**2
       
    ! Read coordinates
    read(infreq1,*) freq1
    read(infreq2,*) freq2
    read(inrad,*)   radius
    freq(1) = freq1
    freq(2) = freq2
    if (myid==0) write(*,*) 'Estimating sindex from ', trim(infreq1), ' GHz map and ', trim(infreq2), ' GHz map'
    if (trim(inlen) /= 'file') then
       read(inlen,*)   len
       read(inbre,*)   bre
       sources = 1
       numrad = 1
       name = trim(outprefix)
       chatty = .true.
       if (myid==0) write(*,*) 'using area with radius ', trim(inrad), &
            & ' degrees, centered around (l =',trim(inlen), ', b =', trim(inbre), ')'
    else
       open(13, file=trim(inbre))
       read(13,*) sources
       numrad = int(radius,i4b)
       chatty = .false.
       prefix = outprefix
    end if
    outfile = trim(outprefix) // '_data.dat'
    open(33, file=trim(outfile), recl=1024)
   
    ! Loop over sources
    do s = 1, sources
       if (trim(inlen) == 'file') then 
          read(13,*) name, len, bre, radius
          outprefix = trim(prefix) // '_' // trim(name)
          radius = max(radius- nint(numrad*0.5), 0.d0)
       else
          radius = radius - 1.d0
       end if
       bre    = bre*pi/180.d0    ! in radians
       len    = len*pi/180.d0    ! in radians
!write(*,*) bre, len
!stop
       ! Loop over radia
       do r = 1, numrad
          radius = radius + 1.d0

          ! Find your designated disk of the sky
          allocate(pixlist(0:npix))
          call ang2vec(pi/2-bre, len, cvec)
          if (ordering==1) then
             call query_disc(nside, cvec, radius*pi/180.d0, pixlist, n)
          else if (ordering==2) then
             call query_disc(nside, cvec, radius*pi/180.d0, pixlist, n, nest=1)
          end if
          if (myid==root) write(*,*) name, n, '= number of pixels inside radius', real(radius,sp)
          allocate(nmap1(n, 7))
          allocate(nmap2(n, 7))
          do i = 1, 3 
             nmap1(:,i)=map1(pixlist(0:n-1),i)
             nmap2(:,i)=map2(pixlist(0:n-1),i)
          end do
          nmap1(:,4)=sqrt(nmap1(:,2)**2+nmap1(:,3)**2)                 ! P
          nmap2(:,4)=sqrt(nmap2(:,2)**2+nmap2(:,3)**2)               

          ! Calculate spectral index
          do i = 1, 4
             call scattercalc_unni(nmap1(:,i), nmap2(:,i), n, spec(i), freq, uncert_syst=uncert(i))
             if (chatty) write(*,*) comp(i), spec(i), uncert(i)
          end do
          diff = (spec(2)-spec(3))/sqrt(uncert(2)**2 + uncert(3)**2)


          ! Write scatter plots to file
!          call int2string(int(radius,i4b), rad)
!          outfile = trim(outprefix) // '_rad' // rad // '_scatter.dat'
!          open(43, file=trim(outfile))
!          do i = 1, 4
!             maxv = maxval(nmap2(:,i))
!             minv = minval(nmap2(:,i))
!             mid2 = (maxv-minv)/2 + minv
!             mid1 = (maxval(nmap1(:,i))-minval(nmap1(:,i)))/2 + minval(nmap1(:,i))
!             nmap1(:,i) = (nmap1(:,i)-mid1)/(maxv-minv)*2.d0
!             nmap2(:,i) = (nmap2(:,i)-mid2)/(maxv-minv)*2.d0
!             do j = 1, n
!                write(43,*) nmap2(j,i), nmap1(j,i)
!             end do
!             write(43,*) ''
!          end do
!          close(43)
!          if (chatty) write(*,*) '* Scatter plot written to file = ', trim(outfile)

          nmap1(:,6)=nmap1(:,4)/nmap1(:,1)                             ! polarisation fraction
          nmap2(:,6)=nmap2(:,4)/nmap2(:,1)
          where(nmap1(:,1)==0.d0) nmap1(:,6) = 0.d0
          where(nmap2(:,1)==0.d0) nmap2(:,6) = 0.d0

          nmap1(:,7)=atan2(nmap1(:,3),nmap1(:,2))/2.d0*180.d0/pi       ! polarisation angle
          nmap2(:,7)=atan2(nmap2(:,3),nmap2(:,2))/2.d0*180.d0/pi     
!          fracs(1) = sum(nmap1(:,6)*nmap1(:,4))*100.d0/sum(nmap1(:,4))
!          fracs(2) = sum(nmap2(:,6)*nmap2(:,4))*100.d0/sum(nmap2(:,4))
          fracs(1) = sum(nmap1(:,6))*100.d0/n
          fracs(2) = sum(nmap2(:,6))*100.d0/n
          call maxlike_alpha(nmap1(:,2), nmap1(:,3), angs(1), outprefix=outprefix)
          call maxlike_alpha(nmap2(:,2), nmap2(:,3), angs(2), outprefix=outprefix)
!          angs(1)  = sum(nmap1(:,6)*nmap1(:,4))/sum(nmap1(:,4))
!          angs(2)  = sum(nmap2(:,6)*nmap2(:,4))/sum(nmap2(:,4))
          if (chatty) write(*,*) 'pol %:   ', fracs
          if (chatty) write(*,*) 'pol angs:', angs

!          comp(4) = "Q'"
!          comp(5) = "U'"
!          vinkel = (angs(1)+22.5)*pi/180.d0
!          nmap1(:,4) = nmap1(:,2)*cos(2*vinkel) + nmap1(:,3)*sin(2*vinkel)
!          nmap1(:,5) = -nmap1(:,2)*sin(2*vinkel) + nmap1(:,3)*cos(2*vinkel)
!          vinkel = (angs(2)+22.5)*pi/180.d0
!          nmap2(:,4) = nmap2(:,2)*cos(2*vinkel) + nmap2(:,3)*sin(2*vinkel)
!          nmap2(:,5) = -nmap2(:,2)*sin(2*vinkel) + nmap2(:,3)*cos(2*vinkel)

          ! Calculate beam ellipticity
          allocate(p(0:5))
          allocate(mapdata(n))
          write(*,*) ''
          write(*,*) '         fwhm         ellipticity      ang_ellip          amp'
          do i = 2, 3
             mapdata = nmap1(:,i)
             p(0) = maxval(mapdata)    !amp
             p(1) = 0.d0               !psi
             p(2) = 0.d0               !x0
             p(3) = 0.d0               !y0      
             p(4) = 1.d0*pi/180.d0     !sx
             p(5) = 1.d0*pi/180.d0     !sy
             call powell(p, ellipse, ierr)
             p(1:5) = p(1:5)*180.d0/pi
             p(1) = modulo(p(1),180.d0)
!             write(*,*) comp(i), ' map1:'
!             write(*,*) p(0),'= amp', p(1),'= psi'
!             write(*,*) p(2),'= x0',  p(3),'= y0'
!             write(*,*) p(4),'= sx',  p(5),'= sy'
             fwhm = 0.5*(abs(p(4))+abs(p(5)))*sqrt(8.d0*log(2.d0))
             a    = max(abs(p(4)), abs(p(5)))
             b    = min(abs(p(4)), abs(p(5)))
             elli = (a-b)/a
!             write(*,*) fwhm,'= fwhm', elli,'= elipticity'
             write(*,*) 'map1 ', comp(i), real(fwhm, sp), real(elli, sp), real(p(1), sp), real(p(0), sp)
          end do
          write(*,*) ''
          do i = 2, 3
             mapdata = nmap2(:,i)
             p(0) = maxval(mapdata)    !amp
             p(1) = 0.d0               !psi
             p(2) = 0.d0               !x0
             p(3) = 0.d0               !y0      
             p(4) = 1.d0*pi/180.d0     !sx
             p(5) = 1.d0*pi/180.d0     !sy
             call powell(p, ellipse, ierr)
             p(1:5) = P(1:5)*180.d0/pi
             p(1) = modulo(p(1),180.d0)
             fwhm = 0.5*(abs(p(4))+abs(p(5)))*sqrt(8.d0*log(2.d0))
             a    = max(abs(p(4)), abs(p(5)))
             b    = min(abs(p(4)), abs(p(5)))
             elli = (a-b)/a
             write(*,*) 'map2 ', comp(i), real(fwhm, sp), real(elli, sp), real(p(1), sp), real(p(0), sp)
          end do

!          fwhm2 = fwhm
!          elli2 = elli
!          psi2  = modulo(p(1),180.d0)

         ! Write map
          allocate(outmap(0:npix-1, 7))
          call int2string(int(radius,i4b), rad)
          call nmap2outmap(outmap, nmap1, pixlist(0:n-1))
          outfile = trim(outprefix) // '_rad' // rad // '_nmap1.fits'
          call write_map(outmap, ordering, outfile)
          if (chatty) write(*,*) '* nmap1 written to file = ', trim(outfile)
          call nmap2outmap(outmap, nmap2, pixlist(0:n-1))
          outfile = trim(outprefix) // '_rad' // rad // '_nmap2.fits'
          call write_map(outmap, ordering, outfile)
          if (chatty) write(*,*) '* nmap2 written to file = ', trim(outfile)
!          outfile = trim(outprefix) // '_rad' // rad // '_rotmap.fits'
!          call write_map(map1, ordering, outfile)
!          if (chatty) write(*,*) '* rotated maps written to file = ', trim(outfile)

          ! Compute scatter plot as a function of polarization angle
          outfile = trim(outprefix) // '_rad' // rad // '_beta_vs_psi.dat'
          open(42, file=trim(outfile))
          outfile = trim(outprefix) // '_rad' // rad // '_scatter_vs_psi.dat'
          open(43, file=trim(outfile))
          n_alpha = 18!0
          allocate(fwhm_psi(0:n_alpha), ell_psi(0:n_alpha), dir_psi(0:n_alpha))
          allocate(rotmap1(size(nmap1,1)), rotmap2(size(nmap2,1)))
          do i = 0, n_alpha
             psi = 90.d0/n_alpha*i
             rotmap1 = nmap1(:,2) * cos(2*psi*DEG2RAD) + nmap1(:,3) * sin(2*psi*DEG2RAD)
             rotmap2 = nmap2(:,2) * cos(2*psi*DEG2RAD) + nmap2(:,3) * sin(2*psi*DEG2RAD)
             call scattercalc_unni(rotmap1, rotmap2, n, beta, freq, uncert_syst=sigma)             
             write(42,*) psi, beta, sigma
             write(*,*) psi, beta, sigma
!             write(42,*) psi, beta, spec(4)

             p(0) = maxval(rotmap1) !amp
             p(1) = 0.d0               !psi
             p(2) = 0.d0               !x0
             p(3) = 0.d0               !y0      
             p(4) = 1.d0*pi/180.d0     !sx
             p(5) = 1.d0*pi/180.d0     !sy
             mapdata = rotmap1
             call powell(p, ellipse, ierr)
             p = p*180.d0/pi
             fwhm = 0.5*(abs(p(4))+abs(p(5)))*sqrt(8.d0*log(2.d0))
             a    = max(abs(p(4)), abs(p(5)))
             b    = min(abs(p(4)), abs(p(5)))
             elli = (a-b)/a
             fwhm_psi(i) = fwhm
             ell_psi(i) = elli
             dir_psi(i)  = modulo(p(1),180.d0)

             maxv = maxval(rotmap1)
             minv = minval(rotmap1)
             mid1 = (maxv-minv)/2 + minv
             mid2 = (maxval(rotmap2)-minval(rotmap2))/2 + minval(rotmap2)
             rotmap1 = (rotmap1-mid1)/(maxv-minv)*2.d0
             rotmap2 = (rotmap2-mid2)/(maxv-minv)*2.d0
             do j = 1, n
                write(43,*) rotmap1(j), rotmap2(j)
             end do
             write(43,*) ''
          end do
          close(42)
          close(43)
!          if (chatty) write(*,*) '* per angle written to file = ', trim(outfile)
!stop


          ! Write beam vs psi data to file
          outfile = trim(outprefix) // '_rad' // rad // '_fwhm_vs_psi.dat'
          open(42, file=trim(outfile), recl=1024)
          write(42,*) real(fwhm_psi,sp)
          close(42)
          outfile = trim(outprefix) // '_rad' // rad // '_ell_vs_psi.dat'
          open(42, file=trim(outfile), recl=1024)
          write(42,*) real(ell_psi,sp)
          close(42)
          outfile = trim(outprefix) // '_rad' // rad // '_dir_vs_psi.dat'
          open(42, file=trim(outfile), recl=1024)
          write(42,*) real(dir_psi,sp)
          close(42)
          deallocate(fwhm_psi, ell_psi, dir_psi)
          deallocate(rotmap1, rotmap2)

          ! Write scatter plots to file
          outfile = trim(outprefix) // '_rad' // rad // '_scatter.dat'
          open(42, file=trim(outfile))
!          do i = 1, 5
          do i = 2, 2
             maxv = maxval(nmap2(:,i))
             minv = minval(nmap2(:,i))
             mid2 = (maxv-minv)/2 + minv
             mid1 = (maxval(nmap1(:,i))-minval(nmap1(:,i)))/2 + minval(nmap1(:,i))
             nmap1(:,i) = (nmap1(:,i)-mid1)/(maxv-minv)*2.d0
             nmap2(:,i) = (nmap2(:,i)-mid2)/(maxv-minv)*2.d0
             do j = 1, n
                write(42,*) nmap2(j,i), nmap1(j,i)
             end do
             write(42,*) ''
          end do
          close(42)
          if (chatty) write(*,*) '* Scatter plot written to file = ', trim(outfile)

          ! Write to master file
          write(33, fmt='(a10,i4,24f10.4)') name, nint(radius), fracs, angs, real(spec,sp), &
               & real(uncert,sp), diff, fwhm1, elli1, psi1, &
               & fwhm2, elli2, psi2, fwhmT, elliT, psiT
      
          ! Clean up
          deallocate(nmap1, nmap2, outmap, pixlist, mapdata, p)
       
       end do
       write(33, *) ''
    end do
    close(33)

  end subroutine scatter_sindex

  !---------------------------------------------------------------------
  ! Likelihood function for powell search
  !---------------------------------------------------------------------

  function ellipse(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: ellipse

    integer(i4b)                                 :: i, pix

    ! some priors to take away degeneracies
!    if (p(5)<0.d0 .or. p(6)<0.d0 .or. p(2)<0.d0 .or. p(2)>pi .or. p(5)<p(6)) then
    if (p(5)<0.d0 .or. p(6)<0.d0 .or. p(5)<p(6)) then
       ellipse = 1d30
       return
    end if

    ellipse = 0.d0
    do i = 1, size(mapdata)
       pix = pixlist(i-1)
       ellipse = ellipse + (mapdata(i) - ellipse_like(pix, p))**2
    end do

  end function ellipse

  !---------------------------------------------------------------------
  ! Ellipse likelihood
  !---------------------------------------------------------------------

  function ellipse_like(pix, p)
    implicit none
    integer(i4b),            intent(in)  :: pix
    real(dp), dimension(0:), intent(in)  :: p
    real(dp)                             :: ellipse_like

    real(dp)  :: amp, psi, x0, y0, sx, sy, x, y, vec(3), dx, dy, a, b, c

    amp = p(0)
    psi = modulo(p(1),pi)
    x0  = p(2)
    y0  = p(3)
    sx  = p(4)
    sy  = p(5)

    if (ordering == 1) then
       call pix2vec_ring(nside, pix, vec)
    else if (ordering == 2) then
       call pix2vec_nest(nside, pix, vec)
    end if
    call project2xy(vec, cvec, x, y)
    x = x * pi/180.d0
    y = y * pi/180.d0

    dx = x - x0
    dy = y - y0
    a  = cos(psi)**2/(2*sx**2) + sin(psi)**2/(2*sy**2)
    b  = - 2*sin(2*psi)/(4*sx**2) + 2*sin(2*psi)/(4*sy**2)
    c  = sin(psi)**2/(2*sx**2) + cos(psi)**2/(2*sy**2)
    ellipse_like = amp*exp(-(a*dx**2 + b*dx*dy + c*dy**2))


  end function ellipse_like
  
  !---------------------------------------------------------------------
  ! project vector to local plane coordinates
  !---------------------------------------------------------------------

  subroutine project2xy(mpixvec, centervec, x, y)
    implicit none

    real(dp)               :: sigma, radius, pixdist, northpixang, psi_true, x, y
    real(dp), dimension(3) :: centervec,  northvec, centerpixvec, centernorthvec, somevec, mpixvec, testvec
   
    psi_true = 0.d0

    northvec=[0,0,1]
    somevec=[centervec(2),-centervec(1),0.d0] 

!write(*,*) somevec, 'somevec'
!write(*,*) mpixvec, 'mpixvec'
!write(*,*) centervec, 'centervec'
 
    call crossproduct(somevec, centervec, centernorthvec)
!write(*,*) centernorthvec, 'centernortvec'

    call crossproduct(centervec, mpixvec, somevec)
    call crossproduct(somevec, centervec, centerpixvec)
    call angdist(centernorthvec, centerpixvec, northpixang)
    call angdist(centervec, mpixvec, pixdist) 
    call crossproduct(centernorthvec, centerpixvec, testvec)
    if (sum(testvec*centervec) > 0) northpixang=-northpixang
    !pixdist   = pixdist*180.d0/pi    ! distance from centervec to given pix in degrees
    !pixdist_x = pixdist*cos(psi_true-northpixang)
    !pixdist_y = pixdist*sin(psi_true-northpixang)
    x = 2.d0*sin(pixdist/2.d0)*cos(psi_true-northpixang)*180.d0/pi
    y = 2.d0*sin(pixdist/2.d0)*sin(psi_true-northpixang)*180.d0/pi
     
  end subroutine project2xy

  !---------------------------------------------------------------------
  ! Crossproduct of two 3-vectors  
  !---------------------------------------------------------------------

  subroutine crossproduct(vector1, vector2, crossvector)
    implicit none

    real(dp), dimension(3), intent(in)  :: vector1, vector2
    real(dp), dimension(3), intent(out) :: crossvector

    crossvector=[vector1(2)*vector2(3)-vector1(3)*vector2(2), vector1(3)*vector2(1)-vector1(1)*vector2(3), &
                                                            & vector1(1)*vector2(2)-vector1(2)*vector2(1) ]

  end subroutine crossproduct

  !---------------------------------------------------------------------
  ! Find maxlike polarisation angle
  !----------------------------------------------------------------------

  subroutine maxlike_alpha(qmap, umap, alpha, uncert, outprefix)
    implicit none

    real(dp),     dimension(:), intent(in)   :: qmap, umap
    real(sp),                   intent(out)  :: alpha
    real(sp),         optional, intent(out)  :: uncert
    character(len=*), optional, intent(in)   :: outprefix

    real(dp),     allocatable,  dimension(:) :: grid
    real(dp)                                 :: a, da, sumq1, sumq2, beta, b
    integer(i4b)                             :: i, j, n, npix, gstep, pos(1)
    character(len=512)                       :: outfile

    npix = size(qmap)
    gstep = 10001
    a = -90.d0
    da =180.d0/real(gstep-1, dp)
    allocate(grid(gstep))
    grid = 0.d0
    do i = 1, gstep
       a = a*pi/180.d0
       do n = 1, npix
          grid(i) = grid(i) + ( umap(n)*cos(2*a) - qmap(n)*sin(2*a) )**2
       end do
       a = a*180.d0/pi
       a = a + da
    end do
    pos = minloc(grid)
    alpha = -90.d0 +real(pos(1)-1, dp)*da
    write(*,*) 'min:', pos(1), alpha, grid(pos(1))
    outfile =  trim(outprefix) // '_maxlike.dat'
    open(17, file = outfile)
    do i = 1, gstep
       write(17, *) -90+(i-1)*da, grid(i)
    end do
    close(17)
    write(*,*) 'Polarisation angle likelihood written to file = ', trim(outfile)

    sumq1 = 0.d0
    sumq2 = 0.d0
    if (alpha < 0.d0) then
       beta = alpha + 90.d0
    else
       beta = alpha - 90.d0
    end if
    a = alpha*pi/180.d0
    b = beta*pi/180.d0
    do n = 1, npix
       sumq1 = sumq1 + qmap(n)*cos(2*a) + umap(n)*sin(2*a)
       sumq2 = sumq2 + qmap(n)*cos(2*b) + umap(n)*sin(2*b)
    end do
!    write(*,*) alpha, sumq1
!    write(*,*) beta, sumq2
    if (sumq2>sumq1) alpha = beta

  end subroutine maxlike_alpha

  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scattercalc_unni(map1, map2, n, spectral, freq, uncert_tot, uncert_stat, uncert_syst)
    ! least squares method when there is noise in both y and x
    ! called the effective variance method.
    implicit none

    real(dp),          dimension(:), intent(in)  :: map1, map2, freq
    integer(i4b),                    intent(in)  :: n
    real(dp),                        intent(out) :: spectral
    real(dp),          optional,     intent(out) :: uncert_stat, uncert_syst, uncert_tot

    integer(i4b)                                 :: i, j
    real(dp)                                     :: healnan=-1.6375d30 
    real(dp)                                     :: m, dm, Vx, Vy, Cxy, dV, nn, aa, b, dtheta
    real(dp)                                     :: mu1, mu2, theta, dp2, sigma_syst, sigma_stat, m_syst(2)


    !y=map2, x=map1   y = mx + b
    nn = n*1.d0
    Vx  = sum(map1*map1)/nn - sum(map1)/nn*sum(map1)/nn
    Vy  = sum(map2*map2)/nn - sum(map2)/nn*sum(map2)/nn 
    Cxy = sum(map1*map2)/nn - sum(map1)/nn*sum(map2)/nn
    dV  = Vx-Vy
    aa  = (Vy - Vx) / ( 2*Cxy )

!    m = ( (Vy-Vx) + sqrt(4*Cxy*Cxy + dV*dV) ) / (2.d0 * Cxy)
!    m = Cxy / Vx  !the method without noise in both y and x

    m = aa + (Cxy/abs(Cxy)) * sqrt(1+aa*aa)
    b = sum(map2)/nn - m * sum(map1)/nn
    
    spectral = log( m ) / log(freq(2)/freq(1))

    if (present(uncert_stat) .or. present(uncert_tot)) then
       dtheta = sqrt( (1.d0/(nn-2.d0)) * (Vx+Vy)/ (dV*dV + 4.d0*Cxy*Cxy) ) 
       sigma_stat = sqrt((1+m**2)**2 * dtheta**2)
       if (present(uncert_stat)) uncert_stat = sqrt(sigma_stat**2 / m**2 / log(freq(2)/freq(1))**2 )
    end if

    if (present(uncert_syst) .or. present(uncert_tot)) then
       ! Divide points into two components, above and below the best-fit line
    
       do j = 1, 2
          nn = 0; Vx = 0; Vy = 0; Cxy = 0; mu1 = 0; mu2 = 0
          do i = 1, n
             if ((j == 1 .and. map2(i) > b + m*map1(i)) .or. (j == 2 .and. map2(i) < b + m*map1(i))) then
                nn  = nn  + 1
                mu1 = mu1 + map1(i)
                mu2 = mu2 + map2(i)
                Vx  = Vx  + map1(i)*map1(i)
                Vy  = Vy  + map2(i)*map2(i)
                Cxy = Cxy + map1(i)*map2(i)
             end if
          end do
          mu1 = mu1 / nn
          mu2 = mu2 / nn
          Vx  = Vx / nn - mu1**2
          Vy  = Vy / nn - mu2**2
          Cxy = Cxy / nn - mu1*mu2
          dV  = Vx-Vy
          aa  = (Vy - Vx) / ( 2*Cxy )
          write(*,*) 'HKE: Commented out line below to make gfortran compile...'
          stop
          !m_syst(j) = aa + (Cxy/abs(Cxy)) * sqrt(1+aa*aa)
       end do
       sigma_syst = 0.5 * abs(m_syst(1)-m_syst(2))
       if (present(uncert_syst)) uncert_syst = sqrt(sigma_syst**2 / m**2 / log(freq(2)/freq(1))**2 )
    end if

    if (present(uncert_tot)) uncert_tot = sqrt((sigma_stat**2 + sigma_syst**2) / m**2 / log(freq(2)/freq(1))**2 )

  end subroutine scattercalc_unni

  !---------------------------------------------------------------------
  ! Find spectral index from real plot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine maxlikeplot(map, sig1, sig2, mask, pol, index, freq1, freq2, spectral, offset, uncert, fg)

    real(dp),        dimension(:,:,:), intent(in)  :: map
    real(dp),        dimension(:),     intent(in)  :: sig1, sig2
    integer(i4b),    dimension(:),     intent(in)  :: mask
    integer(i4b),                      intent(in)  :: pol, index
    real(dp),        dimension(:),     intent(out) :: spectral, offset, uncert
    real(dp),        dimension(:,:),   intent(out) :: fg
    real(dp),                          intent(in)  :: freq1, freq2

    integer(i4b)                                 :: i, j, k, n, p, v, unit, num   
    integer(i4b)                                 :: nfreq, ngrid   
    real(dp)                                     :: healnan=-1.6375d30, tall, nevn, beta
    real(dp)                                     :: betamin, betamax, offmin, offmax, freq0
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:,:) :: polmap

!    if (index/=4677 .and. index/=4701 .and. index/=4697 .and. index/=4344 .and. index/=4333) then
!    if (index/=4701 .and. index/=4524) then
!    if (index/=4865) then
!       spectral = healnan
!       return
!    end if

    num = size(map(:,1,1))
    ! Initialize output
    spectral   = healnan
    offset     = healnan
    uncert     = healnan

    ! Find num of existing pixels for given area
    n = 0
    do i = 1, num
       if ( healok(map(i,pol,1)) .and. healok(map(i,pol,2)) ) then
          n = n+1
       end if
    end do
    ! Return if too few pixels
!    if (n<num/2) return
    if (n<1) return
    write(*,*) n,'=n', num,'=num' , index,'=hpindex'

    ! else proceed
    if (pol==1) then
       call maxlike(map(:,1,1), map(:,1,2), sig1(:), sig2(:), freq1, freq2, &
            & spectral(1), offset(1), uncert(1), n, fg=fg(:,1))!,hpindex=index)
    else if (pol==2) then
       do k=2,3
          write(*,*) k
          call maxlike(map(:,k,1), map(:,k,2), sig1(1+num*(k-2):num*(k-1)), &
               & sig2(1+num*(k-2):num*(k-1)), freq1, freq2, spectral(k), offset(k), &
               & uncert(k), n, fg=fg(:,k), polbeta=spectral(1), poloff=offset(1))
       end do
       allocate(polmap(pol*num,2))
       polmap(1:num,:)       = map(:,2,:)
       polmap(num+1:2*num,:) = map(:,3,:)
       call maxlike(polmap(:,1), polmap(:,2), sig1(:), sig2(:), freq1, freq2, spectral(1), &
            & offset(1), uncert(1), n*pol)!, hpindex=index)
       deallocate(polmap)
    end if

  end subroutine maxlikeplot

  !---------------------------------------------------------------------
  ! Find spectral index from maximum likelihood
  !----------------------------------------------------------------------

  subroutine maxlike(map1, map2, sig1, sig2, freq1, freq2, spectral, offset, uncert, n, fg, polbeta, poloff, hpindex)

    real(dp),               dimension(:), intent(in)  :: map1, map2, sig1, sig2
    real(dp),                             intent(out) :: spectral, offset, uncert
    real(dp),                             intent(in)  :: freq1, freq2
    integer(i4b),                         intent(in)  :: n
    integer(i4b), optional,               intent(in)  :: hpindex
    real(dp),     optional, dimension(:), intent(out) :: fg
    real(dp),     optional,               intent(in)  :: polbeta, poloff

    integer(i4b)                                 :: i, j, p, v, unit, num, int, k
    integer(i4b)                                 :: nfreq, ngrid   
    real(dp)                                     :: healnan=-1.6375d30, tall, beta, finebeta, fineoffset
    real(dp)                                     :: betamin, betamax, offmin, offmax, freq0
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    integer(i4b),    allocatable, dimension(:)   :: in2red
    real(dp),        allocatable, dimension(:)   :: pos, m, freq, A, a2t, prob, nevn
    real(dp),        allocatable, dimension(:,:) :: grid, y, sigsq
    logical(lgt)  :: chatty

    chatty = .false.

    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if
    num = size(map1)
    if (chatty) write(*,*) num,'=num', n ! size(fg),'= size(fg)', n

    ! Put data into place
    nfreq=2
    allocate(freq(nfreq))
    freq(1) = freq1
    freq(2) = freq2
    freq0 = 23.d0
    allocate(a2t(nfreq))
    do v = 1, nfreq
       a2t(v)=ant2thermo(freq(v))
    end do
    allocate(m(nfreq))
    m(1)=0.d0 
    allocate(y(n,nfreq))
    allocate(sigsq(n,nfreq))
    allocate(in2red(n))
    j=0
    do i = 1, num
       if ( healok(map1(i)) .and. healok(map2(i))) then
          j = j+1
          y(j,1) = map1(i)
          y(j,2) = map2(i)
          sigsq(j,1) = sig1(i)**2
          sigsq(j,2) = sig2(i)**2
          in2red(j) = i
       end if
    end do

    ! Calculate course grid = chi_square
    ngrid   = 100
    betamin = -10.d0
    betamax = 10.d0
    offmin  = -100.d0
    offmax  = 100.d0
    allocate(grid(ngrid, ngrid))
    allocate(A(n))
    allocate(nevn(n))
    allocate(pos(nfreq))

    do k = 1, 2

       ! If you want an additional loop, loop from 1 to 3, uncomment this and set the one below for k=3.
!!$       if (k == 2) then ! Calculate fine grid = chi_square
!!$          finebeta = 1.d0
!!$          fineoffset = 10.d0
!!$          beta = (betamax-betamin)/(ngrid-1)*(pos(2)-1) + betamin
!!$          m(2) = (offmax-offmin)/(ngrid-1)*(pos(1)-1) + offmin
!!$          betamin = beta -finebeta
!!$          betamax = beta +finebeta
!!$          offmin  = m(2) -fineoffset
!!$          offmax  = m(2) +fineoffset
!!$          if (chatty) write(*,*) betamin, 'betamin, betamax', betamax
!!$          if (chatty) write(*,*) offmin, 'offmin, offmax', offmax
!!$       end if

       if (k == 2) then ! Calculate grid 3 based on confidence levels on chi_square
            ! For a Gaussian distribution, the 1, 2 and 3 sigma (68%, 95% and 99.7%) confidence regions 
            ! correspond to where chi_square = -2 lnL increases by 2.3, 6.17 and 11.8 from its minimum value.
          do j = pos(2),ngrid      ! loop over beta to find where grid = 3sigma
             if ( grid(pos(1),j) >= 11.8 + minval(grid) ) exit
          end do
          do i = pos(1),ngrid      ! loop over m to find where grid = 3 sigma 
             if ( grid(i,pos(2)) >= 11.8 + minval(grid) ) exit
          end do
          finebeta = (betamax-betamin)/(ngrid-1)*(j - pos(2))
          fineoffset = (offmax-offmin)/(ngrid-1)*(i - pos(1))
          beta = (betamax-betamin)/(ngrid-1)*(pos(2)-1) + betamin
          m(2) = (offmax-offmin)/(ngrid-1)*(pos(1)-1) + offmin
          betamin = beta -finebeta
          betamax = beta +finebeta
          offmin  = m(2) -fineoffset
          offmax  = m(2) +fineoffset
          if (chatty) write(*,*) grid(pos(1),j), 'grid(pos(1),j), grid(i,pos(2))', grid(i,pos(2))
          if (chatty) write(*,*) j- pos(2), 'number grid beta, number grid m', i - pos(1)
          if (chatty) write(*,*) betamin, 'betamin, betamax', betamax
          if (chatty) write(*,*) offmin, 'offmin, offmax', offmax       
       end if

       grid= 0.d0
       do j = 1, ngrid
          beta = (betamax-betamin)/(ngrid-1)*(j-1) + betamin
          do i = 1, ngrid
             m(2) = (offmax-offmin)/(ngrid-1)*(i-1) + offmin
             A = 0.d0
             nevn = 0.d0
             do v = 1, nfreq
                A(:) = A(:) + (y(:,v)-m(v))*a2t(v)/sigsq(:,v)*(freq(v)/freq0)**beta
                nevn(:) = nevn(:) + ((a2t(v)**2)*(freq(v)/freq0)**(2.d0*beta))/sigsq(:,v)
             end do
             A = A/nevn           
             do v = 1, nfreq
                do p = 1, n ,16
                   grid(i,j) = grid(i,j) +((y(p,v)-A(p)*a2t(v)*(freq(v)/freq0)**beta-m(v))**2)/sigsq(p,v)
                end do
             end do
          end do
       end do
       pos=minloc(grid)
       if (chatty) write(*,*) 'positionon grid ',k,': ', pos(1), pos(2)
       if (pos(1)==1 .or. pos(1)==ngrid .or. pos(2)==1 .or. pos(2)==ngrid) then
          spectral = healnan
          offset   = healnan
          uncert   = healnan
          exit
       end if
    end do


    ! Check that we didn't not hit the grid
    if (pos(1)==1 .or. pos(1)==ngrid .or. pos(2)==1 .or. pos(2)==ngrid) then
       uncert   = healnan
       spectral = healnan
       offset   = healnan
       if (present(fg)) fg = healnan
    else
       ! Find spectral index and offset
       beta = (betamax-betamin)/(ngrid-1)*(pos(2)-1) + betamin
       spectral = beta
       m(2) = (offmax-offmin)/(ngrid-1)*(pos(1)-1) + offmin
       offset = m(2)
       if (chatty) write(*,*) 'Spectral index', spectral
       if (chatty) write(*,*) 'Offset', offset
       ! Calculate uncertainty
       allocate(prob(ngrid))
       tall = minval(grid)
       grid = exp(-(grid-tall)/2.d0)
       prob=sum(grid,1)
       prob=prob/sum(prob)
       prob = -2.d0*log(prob)
       if (chatty) write(*,*) 'prob', prob(pos(2)), minval(prob), pos(2)
       do i = pos(2),ngrid
          if (prob(i) >= minval(prob) + 1.d0) exit
       end do
       int=i
       do i = pos(2),1, -1
          if (prob(i) >= minval(prob) + 1.d0) exit
       end do
       if (int > ngrid .or. i < 1) then
          uncert = healnan
       else
          uncert = (betamax-betamin)/(ngrid-1)*(int-i)/2.d0
       end if
       if (uncert==0.d0) uncert=healnan
       if (chatty) write(*,*) 'delta beta', uncert, int,i
       if (chatty) write(*,*)
       deallocate(prob)
       ! Calculate foreground
       if (present(fg)) then
          if (present(polbeta) .and. present(poloff)) then
             beta = polbeta    ! obs changing to pol(=Q+U) values
             m(2) = poloff     ! obs changing to pol(=Q+U) values
          end if
          A = 0.d0
          nevn = 0.d0
          do v = 1, nfreq
             A(:) = A(:) + (y(:,v)-m(v))*a2t(v)/sigsq(:,v)*(freq(v)/freq0)**beta
             nevn(:) = nevn(:) + ((a2t(v)**2)*(freq(v)/freq0)**(2.d0*beta))/sigsq(:,v)
          end do
          A = A/nevn           
          A = A*ant2thermo(freq0)
          fg(in2red)=A(1:size(fg))
          if (chatty) write(*,*) 'fg', fg(in2red(1)), fg(1)
       end if
    end if

    if (present(hpindex)) then
       ! Write contour  plot to file
       call int2string(hpindex, filnum)
       filename='contour_'//filnum//'.txt'
       open(42, file=trim(filename))
       write(42,*) ngrid, ngrid
       do j = 1, ngrid
          beta = (betamax-betamin)/(ngrid-1)*(j-1) + betamin
          do i = 1, ngrid
             m(2) = (offmax-offmin)/(ngrid-1)*(i-1) + offmin
             write(42,*) m(2), beta, grid(i,j)
          end do
       end do
       close(42)
       write(*,*) 'written to file = ', trim(filename)
    end if

    
    ! Clean up
    deallocate(A)
    deallocate(nevn)
    deallocate(pos)
    deallocate(grid)
    deallocate(y)
    deallocate(sigsq)
    deallocate(m)
    deallocate(freq)
    deallocate(a2t)
 
  end subroutine maxlike

  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scatterplot(map1, map2, mask, pol, index, spectral, freq1, freq2, offset, outprefix)

    real(dp),          dimension(:), intent(in)  :: map1, map2
    integer(i4b),      dimension(:), intent(in)  :: mask
    integer(i4b),                    intent(in)  :: pol, index
    real(dp),          dimension(:), intent(out) :: spectral, offset
    real(dp),                        intent(in)  :: freq1, freq2
    character(len=*)                             :: outprefix 

    integer(i4b)                                 :: lmax, n, l, i, j, k, unit, num   
    real(dp)                                     :: healnan=-1.6375d30, tall
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:)   :: y, beta, vector
    real(dp),        allocatable, dimension(:,:) :: x, matrix2, h
    real(dp),                     dimension(3)   :: slope
    integer(i4b),    allocatable, dimension(:)   :: polmask

!write(*,*) 'scatterplot', index

!   if (index/=4677 .and. index/=4701 .and. index/=4697 .and. index/=4344) then
!   if (index/=191) then
!       spectral = -100.d0!healnan
!       offset = -100.d0 !healnan
!       return
!    end if

    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if

    num = size(map1)/pol

    ! Find num of existing pixels for given area
    n = 0
    do i = 1, num
!       if ( map1(i)/=healnot .and. map2(i)/=healnot .and. map1(i)>200.d0) then
       if ( healok(map1(i)) .and. healok(map2(i))) then
          n = n+1
       end if
    end do

    ! Return if too few pixels
    if (n<num/2) then
       spectral = healnan
       return
    end if

    ! Calculate scale factor for Q and U separately
    if (pol==2) then
       call scattercalc(map1(1:num), map2(1:num), mask, n, index, spectral(2), freq1, freq2)
       call scattercalc(map1(num+1:2*num), map2(num+1:2*num), mask, n, index, spectral(3), freq1, freq2)
       allocate(polmask(num*pol))
       polmask(1:num)       = mask
       polmask(num+1:2*num) = mask
       ! Calculate scale factor for T or P separately
       n=n*pol
       call scattercalc(map1, map2, polmask, n, index, spectral(1), freq1, freq2)
       deallocate(polmask)    
    else if (pol==1) then
       call scattercalc(map1, map2, mask, n, index, spectral(1), freq1, freq2, offset(1), outprefix)
    end if

return

   ! Write scatter plot to file
    call int2string(index, filnum)
    filename='scatter_'//filnum//'.txt'
    open(42, file=trim(filename))
    do i = 1, size(map1)/2
       if ( healok(map1(i)) .and. healok(map2(i))) write(42,*) map2(i), map1(i)
    end do
!    write(42,*)
    do i = size(map1)/2 +1, size(map1)
       if ( healok(map1(i)) .and. healok(map2(i))) write(42,*) map2(i), map1(i)
    end do
!    write(42,*)
!    write(42,*) 0.d0, beta(2)
!    write(42,*) -100.d0, -100*beta(1) +beta(2)
!    write(42,*) 150.d0, 150*beta(1) +beta(2)
!    close(42)
    write(*,*) 'written to file = ', trim(filename)

return

     ! Write scatter plot to file
    call int2string(index, filnum)
    filename='scatter_'//filnum//'.txt'
    open(42, file=trim(filename))
    do i = 1, n
       write(42,*) y(i), x(i,1)
    end do
    close(42)
    write(*,*) 'written to file = ', trim(filename)

  end subroutine scatterplot

  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scattercalc(map1, map2, mask, n, index, spectral, freq1, freq2, offset, outprefix)

    real(dp),          dimension(:), intent(in)  :: map1, map2
    integer(i4b),      dimension(:), intent(in)  :: mask
    integer(i4b),                    intent(in)  :: n, index
    real(dp),                        intent(out) :: spectral
    real(dp),                        intent(in)  :: freq1, freq2
    real(dp),          optional,     intent(out) :: offset
    character(len=*),  optional                  :: outprefix 

    integer(i4b)                                 :: i, j, p, k, numpairs
    real(dp)                                     :: healnan=-1.6375d30, tall, abugmed, bbugmed, amp
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:)   :: y, beta, vector, ymap, xmap, a, b, abug
    real(dp),        allocatable, dimension(:,:) :: x, matrix2, h
    real(dp)                                     :: slope(3), vec(2), amed, bmed, xmin, xmax, xlimit
   
    allocate(ymap(n))
    allocate(xmap(n))
    xmap=0.d0
    ymap=0.d0
    p=0
    amp=3.d0
    do i = 1, size(map1)
       !       if ( map1(i)/=healnot .and. map2(i)/=healnot .and. map1(i)>200.d0) then
!       if ( healok(map1(i)) .and. healok(map2(i)) .and. mask(i)==1.d0 .and. map2(i)<amp .and. -amp<map2(i)) then
       if ( healok(map1(i)) .and. healok(map2(i)) .and. mask(i)==1.d0) then
          p = p+1
          ymap(p) = map1(i)
          xmap(p) = map2(i)
       end if
    end do
!write(*,*) 'using',p,'of',n,'points'
    allocate(y(p))
    allocate(beta(2))
    allocate(x(p,2))
    allocate(matrix2(2,2))
    matrix2(1,1) = p*1.d0
    matrix2(1,2) = -sum(xmap)
    matrix2(2,1) = matrix2(1,2)
    matrix2(2,2) = sum(xmap*xmap)
    matrix2 = matrix2/(p*sum(xmap*xmap)-(sum(xmap))**2)
    vec(1) = sum(xmap*ymap)
    vec(2) = sum(ymap)
    beta = matmul(matrix2, vec) 

    numpairs = p*(p-1)/2
    allocate(a(numpairs))
    allocate(abug(numpairs))
    allocate(b(p))
   k = 0
    do i = 1, p
       do j = i+1, p	
 	     k = k+1
	     a(k) = (ymap(i)-ymap(j))/(xmap(i)-xmap(j))
!	  abug(k) = abs(ymap(i)-ymap(j))/abs(xmap(i)-xmap(j))
!	  if (xmap(i)>xmap(j)) then
!	     a(k) = (ymap(i)-ymap(j))/(xmap(i)-xmap(j))
!	  else
!	     a(k) = (ymap(j)-ymap(i))/(xmap(j)-xmap(i))
!	  end if   
       end do
    end do   
    if (k /= numpairs) then
       write(*,*) 'wrong again!'
       stop
    end if   
    amed = median(a)
!    abugmed = median(abug)
    b = ymap - amed*xmap
    bmed = median(b)
    spectral = amed
    offset = bmed
!    b = ymap - abugmed*xmap
!    bbugmed = median(b)

   ! Write scatter plot to file
    xmin = minval(xmap)
    xmax = maxval(xmap)
    call int2string(index, filnum)
    filename=trim(outprefix)//'_scatter_'//filnum//'.dat'
    open(42, file=trim(filename))
    do i = 1, n
       write(42,*) xmap(i), ymap(i)
    end do
    write(42,*)
    write(42,*) 0.d0, bmed
    write(42,*) xmin, xmin*amed + bmed
    write(42,*) xmax, xmax*amed + bmed
!    write(42,*)
!    write(42,*) 0.d0, bbugmed
!    write(42,*) xmin, xmin*abugmed + bbugmed
!    write(42,*) xmax, xmax*abugmed + bbugmed
    write(42,*)
    write(42,*) 0.d0, beta(2)
    write(42,*) xmin, xmin*beta(1) + beta(2)
    write(42,*) xmax, xmax*beta(1) + beta(2)
    close(42)
!    write(*,*) 'written to file = ', trim(filename)


    write(*,*) int(p,i2b), int(index,i2b), real(amed,sp), real(bmed,sp), real(beta(1),sp), real(beta(2),sp)

!if (abs(beta(1)-amed)>0.03d0) then
!    spectral = healnan
!    offset = healnan
!    return
!end if



   ! Clean up
    deallocate(ymap)
    deallocate(xmap)
    deallocate(a)
    deallocate(b)
    deallocate(y)
    deallocate(x)
    deallocate(beta)
    deallocate(matrix2)
return




spectral = beta(1)  
!    spectral = log(beta(1))/log(freq1/freq2)
    if (present(offset)) then
    offset = beta(2)
    end if
    slope(1)    = beta(1)
    write(*,*) beta(1), beta(2), index, n
!    write(*,*) beta(1), spectral, index, n
 !   write(*,*)
!    if (present(uncert)) then
!       ! Uncertainty
!       allocate(h(n,n))
!       allocate(vector(n))
!       h = matmul(x, matmul(matrix2, transpose(x)))     !!!!!!!!!!!!!!!!
!       h = get_identity(n) - h
!       vector = matmul(h,y)                              !!!!!!!!!!!!!!!!!!1
!       tall =sum( y*vector)
!       matrix2 = matrix2*tall/(n-2) !sjekk n-2
!       uncert = matrix2(1,1)
!       deallocate(h)
!       deallocate(vector)
!    end if

   ! Clean up
    deallocate(y)
    deallocate(x)
    deallocate(beta)
    deallocate(matrix2)


  end subroutine scattercalc

  !---------------------------------------------------------------------
  ! Estimate leakage map1 vs map2
  !----------------------------------------------------------------------

  subroutine estimate_leakage(map1, map2, nside, radius, leakage, uncert)

    real(dp),          dimension(:), intent(in)  :: map1, map2
    integer(i4b),                    intent(in)  :: nside
    real(dp),                        intent(in)  :: radius
    real(dp),                        intent(out) :: leakage, uncert

    integer(i4b)                                 :: i, j, p, n, num, nlist, pix   
    real(dp)                                     :: healnan=-1.6375d30, teller, nevner, radradius

    real(dp),        allocatable, dimension(:)   :: quiet, wmap
    real(dp),                     dimension(3)   :: gc_vec
    integer(i4b),    allocatable, dimension(:)   :: listpix

    ! Check that input maps are of same size
    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if

    ! Find pixels in a 1deg radius around gc center
    gc_vec = [1,0,0]
    radradius = radius*pi/180.d0   ! radius in radians
    allocate(listpix(0:nside**2))
    call query_disc(nside, gc_vec, radradius, listpix, nlist, nest=1) 
    
    ! Find num of existing pixels for given area
    n = 0
    do p=0,nlist-1
       pix=listpix(p)
       if ( healok(map1(pix)) .and. healok(map2(pix))) then
          n = n+1
       end if
    end do
    ! Return if too few pixels
    if (n==0) then
       leakage = healnan
       return
    end if
    write(*,*) n, '= number of pixels'

    ! Calculate leakage parameter
    allocate(quiet(n))
    allocate(wmap(n))
    j=0
    do p=0,nlist-1
       pix=listpix(p)
       if ( healok(map1(pix)) .and. healok(map2(pix))) then
          j = j+1
          wmap(j)  = map1(pix)
          quiet(j) = map2(pix)
       end if
    end do

    teller = 0.d0
    nevner = 0.d0
    do i = 1,n
       teller = teller + quiet(i)*wmap(i)
       nevner = nevner + wmap(i)*wmap(i)
    end do
    leakage = teller/nevner

    uncert = healnan

   ! Clean up
    deallocate(quiet)
    deallocate(wmap)
    deallocate(listpix)

  end subroutine 

  !---------------------------------------------------------------------
  ! Unconvolve map by beamfile, and then convolve with gaussian
  !----------------------------------------------------------------------

  subroutine beam_convolve(map, ordering, nside, nmaps, beamfile, fwhm, onlypol)

    real(dp),          dimension(0:,:), intent(inout) :: map
    character(len=256),                 intent(in)    :: beamfile
    real(dp),                           intent(in)    :: fwhm
    integer(i4b),                       intent(in)    :: nside, nmaps, ordering
    logical(lgt),      optional,        intent(in)    :: onlypol

    complex(dpc),     allocatable, dimension(:,:,:)  :: alm
    real(dp),         allocatable, dimension(:,:)    :: dw8
    real(dp),                      dimension(2)      :: z
    real(dp),         allocatable, dimension(:, :) :: beam
    integer(i4b)      :: lmax, n, l, i, j, m, unit   
    real(dp)          :: tall, sigma
    real(dp),     allocatable, dimension(:)     :: power
    logical(lgt) :: useonlypol

    useonlypol = .false.
    if (present(onlypol)) then
       if (onlypol) useonlypol = .true.
    end if

    ! Read beam file
    if (myid==0) write(*,*) 'Reading from ', trim(beamfile)
    lmax = min(3*nside, 649)
    allocate(beam(0:lmax, 1:nmaps))
    call read_beam(beamfile, beam)
    
    ! Convert map from nest to ringed
    if (ordering ==2) then                       ! map2alm requieres ringed
       call convert_nest2ring(nside ,map)
    end if

!!$    ! output beams
!!$    sigma=fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0))
!!$    open(13, file='totbeam_quiet.txt')
!!$    do l = 0, lmax
!!$       write(13,*) l, 1.d0/beam(l, 2)* exp(-0.5d0*l*(l+1)*sigma**2)
!!$    end do
!!$    close(13)
!!$
!!$    open(13, file='beam_quiet_Q.txt')
!!$    do l = 0, lmax
!!$       write(13,*) l, beam(l, 2)
!!$    end do
!!$    close(13)

    ! Finding alm's for given skymap
    allocate(dw8(1:2*nside, 1:nmaps))
    allocate(alm(1:nmaps, 0:lmax, 0:lmax))
    dw8 = 1.d0
    z = 0.d0
    if (nmaps==1) then
       call map2alm(nside, lmax, lmax, map(:,1), alm, z, dw8)
    else
       call map2alm(nside, lmax, lmax, map, alm, z, dw8)
    end if
    deallocate(dw8)

!    ! and power spectrum
!    allocate(power(0:lmax))
!    power = 0.d0
!    do i = 0, lmax
!       do m = 0, i
!          if (m == 0) then
!             power(i) = power(i) + abs(alm(2, i, m))**2             
!          else
!             power(i) = power(i) + 2*abs(alm(2, i, m))**2
!          end if
!          !write(*,*) i,m,alm(npol,i,m)       
!       end do
!       power(i) = power(i)/(2.d0*i+1.d0)*i*(i+1.d0)/2.d0*pi
!    end do
!    open(13, file='power_before_quiet.txt')
!    do l = 0, lmax
!       write(13,*) l, power(l)
!    end do
!    close(13)

    ! Modifying beam
    sigma=fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0))
    do i = 1, nmaps
       do l = 0, lmax
          alm(i, l, :) = alm(i, l, :)/beam(l, i) * exp(-0.5d0*l*(l+1.d0)*sigma**2)
       end do
    end do
    deallocate(beam)

!    ! and power spectrum
!    power = 0.d0
!    do i = 0, lmax
!       do m = 0, i
!          if (m == 0) then
!             power(i) = power(i) + abs(alm(2, i, m))**2             
!          else
!             power(i) = power(i) + 2*abs(alm(2, i, m))**2
!          end if
!    !write(*,*) i,m,alm(npol,i,m)       
!       end do
!       power(i) = power(i)/(2.d0*i+1.d0)*i*(i+1.d0)/2.d0*pi
!    end do
!    open(13, file='power_after_quiet.txt')
!    do l = 0, lmax
!       write(13,*) l, power(l)
!    end do
!    close(13)

    ! Back to map again
    if (nmaps==1) then
       call alm2map(nside, lmax, lmax, alm, map(:,1))
    else
       call alm2map(nside, lmax, lmax, alm, map)
    end if

    deallocate(alm)

    ! Convert map from ringed to nest
    call convert_ring2nest(nside ,map)

    ! put pol=sqrt(Q^2+U^2) in map(:,1)
    if (useonlypol) map(:,1) = sqrt(map(:,2)**2 + map(:,3)**2)
    
  end subroutine beam_convolve

  !-----------------------------------------------------------------------------------------------
  ! subroutine read_mapfile
  !-----------------------------------------------------------------------------------------------

  subroutine read_mapfile(infile, npix, nside, nmaps, ordering, map, nest, check_nside, check_nmaps, check_ordering)
    implicit none

    character(len=*),           intent(in)    :: infile
    real(dp),   dimension(:,:), pointer       :: map
    integer(i4b),               intent(out)   :: npix, nmaps, nside, ordering
    logical(lgt),     optional, intent(in)    :: nest
    integer(i4b),     optional, intent(in)    :: check_nmaps, check_nside, check_ordering

    integer(i4b)                   :: root=0
    real(dp)                       :: nullval
    logical(lgt)                   :: anynull, mismatch
   
    ! Read size of map
    if (myid==root) write(*,*) 'Reading from file ', trim(infile)
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,' ,ordering,'= ordering' 
    
    ! Checking that reference map and this map is of same size
    mismatch = .false.
    if (present(check_nside)) then    
       if (nside /= check_nside) mismatch = .true.
    end if
    if (present(check_nmaps)) then    
       if (nmaps /= check_nmaps) mismatch = .true.
    end if
    if (present(check_ordering)) then    
       if (ordering /= check_ordering) mismatch = .true.
    end if
    if (mismatch) then
       if (myid==root) write(*,*) trim(infile),'is not of same size as last map/reference map. Quiting'
       call mpi_finalize(ierr)
       stop
    end if
    
    ! Read map  and convert to nest if necassary
    allocate(map(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull)
    if (present(nest) .and. nest) then
       if (ordering == 1) then                       ! we want nested
          write(*,*) 'Converting from ring to nest'
          call convert_ring2nest(nside ,map)
          ordering = 2
       end if
    end if

  end subroutine read_mapfile

  subroutine process_ASCII_table(unit)
    implicit none

    character(len=256) :: infile, outfile, intext
    integer(i4b) :: i, n, firstcol, lastcol, unit
    real(dp), allocatable, dimension(:,:) :: data

    unit = 100

    ! Get parameters
    if (iargc() /= 5) then
       write(*,*) 'process_ASCII_table takes 4 parameters'
       call give_user_info
    else 
       call getarg(2, infile)
       call getarg(3, intext)
       read(intext,*) firstcol
       call getarg(4, intext)
       read(intext,*) lastcol
       call getarg(5, outfile)
    end if

    allocate(data(1000,lastcol))
    open(unit, file=trim(infile),recl=1024)
    n = 0
    do while (.true.)
       n = n+1
       read(unit,*,end=99) intext, data(n,:)
!       write(*,*) n, real(data(n,:),sp)
       write(*,*) n, trim(intext)
    end do
99  close(unit)
    
    open(unit,file=trim(outfile))
    do i = firstcol, lastcol
       write(unit,*) i, data(1,i), sqrt(variance(data(2:n-1,i)))
    end do
    close(unit)

    deallocate(data)

  end subroutine process_ASCII_table

  !-----------------------------------------------------------------------------------------------
  ! subroutine output mdmaps
  !-----------------------------------------------------------------------------------------------

  subroutine output_mdmaps(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: infile, string, outprefix, outfile, buffer(5)
    character(len=4)            :: nsidetxt
    integer(i4b)                :: i, numband, nside, npix

    real(dp),          allocatable, dimension(:,:)   :: mdval
    real(dp),          allocatable, dimension(:)     :: dmap
    character(len=64), allocatable, dimension(:)     :: name

    ! Get parameters
    if (iargc() /= 5) then
       write(*,*) 'mdmaps takes 4 parameters'
       call give_user_info
    else 
       call getarg(2, string)
       read(string,*) numband
       call getarg(3, infile)
       call getarg(4, string)
       read(string,*) nside
       call getarg(5, outprefix)
    end if

    ! Read md values from txt file
    if (myid==root) write(*,*) 'Reading from file ', trim(infile)
    allocate(mdval(1:numband,1:4))
    allocate(name(1:numband))
    open(13, file=trim(infile))
    do i = 1, numband
       read(13,*) buffer
       name(i) = buffer(1)
       read(buffer(2:5),*) mdval(i,:) 
       write(*,*) trim(name(i)), real(mdval(i,:),sp)
    end do
    close(13)

!    mdval(:,2:4) = 0.d0

    ! Make and output maps
    ordering = 1
    call int2string(nside, nsidetxt)
    npix = 12*nside**2
    allocate(dmap(0:npix-1))
    do i = 1, numband
       dmap=0.d0
       call add_dipole(nside, dmap, ordering, 2, mdval(i,:))
       outfile = trim(outprefix) //'_'//trim(name(i))//'_n'//nsidetxt//'.fits'
       call write_map(dmap, ordering, trim(outfile))
       write(*,*) '* Map written to file = ', trim(outfile)
    end do
    deallocate(mdval,name, dmap)

  end subroutine output_mdmaps

  !---------------------------------------------------------------------------------
  ! subroutine beta
  !---------------------------------------------------------------------------------

  subroutine make_weighted_map(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: outfile, weightmapfile, map1file, val_in
    integer(i4b)                :: i, j, numfiles, npix, ordering, nside, nmaps
    integer(i4b)                :: npix1, ordering1, nside1, nmaps1, comp
    real(dp)                    :: nullval, val2, minval, topval, lowval, highval, weight
    logical(lgt)                :: anynull

    real(dp),           allocatable, dimension(:,:) :: weightmap, map1, outmap

    ! Get parameters
    if (iargc() /= 7) then
       write(*,*) 'beta takes 6 parameters: weightmap, map1, val2, outfile lowcut highcut'
       write(*,*) ''
       call give_user_info
    else 
       call getarg(2, weightmapfile)
       call getarg(3, map1file)
       call getarg(4, val_in)
       read(val_in,*) val2
       call getarg(5, outfile)
       call getarg(6, val_in)
       read(val_in,*) lowval
       call getarg(7, val_in)
       read(val_in,*) highval
    end if

    ! Read maps
    npix = getsize_fits(weightmapfile, nmaps=nmaps, ordering=ordering, nside=nside)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,' ,ordering,'= ordering'
    npix1 = getsize_fits(map1file, nmaps=nmaps1, ordering=ordering1, nside=nside1)
    if (myid==root) write(*,*) nside1, '= nside,', nmaps1, '= nmaps,' ,ordering1,'= ordering'

    ! Check maps
    if (nside1 /= nside) then
       if (myid==root) write(*,*) 'Different nside. Quiting'
       stop
    end if
 
    ! Actually read maps
    allocate(weightmap(0:npix-1,nmaps))
    call read_bintab(weightmapfile, weightmap, npix, nmaps, nullval, anynull)
    allocate(map1(0:npix-1,nmaps1))
    call read_bintab(map1file, map1, npix1, nmaps1, nullval, anynull)
 
    ! convert to same ordering if necessary
    if (ordering1 /= ordering) then
       if (myid==root) write(*,*) 'Different ordering. Converting map'
       if (ordering1 == 1) then
          call convert_nest2ring(nside, weightmap)
       else if (ordering1 == 2) then
          call convert_ring2nest(nside, weightmap)
       end if
       ordering=ordering1
    end if

    ! Prepare syntehzized map
    allocate(outmap(0:npix-1,nmaps))
    j=1
    do i = 0, npix-1
       if (weightmap(i,j) > highval) then
          outmap(i,j) = map1(i,j)
       else if (weightmap(i,j) < lowval) then
          outmap(i,j) = val2
       else
          weight = (weightmap(i,j)-lowval)/(highval-lowval)
          if (weight>1.d0 .or. weight<0.d0) write(*,*) weight
          outmap(i,j) = weight * map1(i,j) + (1-weight) * val2
       end if
    end do
    write(*,*) lowval, '= lowval', highval, '= highval'

    ! Write map
    call write_map(outmap, ordering, trim(outfile))
    write(*,*) '* Map written to file = ', trim(outfile)
    
    ! Clean up
    deallocate(weightmap, map1, outmap)

  end subroutine make_weighted_map

  !---------------------------------------------------------------------------------
  ! subroutine gain
  !---------------------------------------------------------------------------------

  subroutine estimate_gain_from_map(unit)
    implicit none
    
    integer(i4b),    intent(in) :: unit
    
    character(len=512)          :: mapfile, outprefix, val_in
    character(len=4)            :: ftxt
    character(len=2)            :: sbtxt
    integer(i4b)                :: i, k, l, p, s, nfreq, nsb, nx, ny
    integer(i4b)                :: numsamp, burn, numrej, numpar
    real(dp)                    :: gain, mingain, maxgain, minbeam, maxbeam
    real(dp)                    :: r, xjup, yjup, x, y
    real(dp)                    :: gain0, fwhm0, xpos0, ypos0, mono
    real(dp)                    :: abstemp, obstemp, fwhm, nu0, dist, dist_fid
    real(dp)                    :: inittemp, initfwhm, initdist, tsys0, tsys
    real(dp)                    :: f_A, f_D, sigma, omega_b, omega_ref, ratio
    logical(lgt)                :: second_time
    real(dp),     dimension(5)  :: nu, T_p, T_p2
    real(dp),     dimension(5)  :: stepsize, parmean, parstdd
    integer(i4b), dimension(4)  :: pos

    real(dp),         allocatable, dimension(:,:,:) :: ipar, fpar
    real(dp),         allocatable, dimension(:,:)   :: params
    real(dp),         allocatable, dimension(:)     :: chisq
    character(len=4), allocatable, dimension(:)     :: parname

    type(map_type)  :: inmap

    call initialize_random_seeds(MPI_COMM_WORLD, 357888, rng_handle)

    ! Get parameters
    if (iargc() /= 5) then
       write(*,*) 'gain takes 4 parameters: map, prefix, numsamp, burn'
       write(*,*) ''
       call give_user_info
    else 
       call getarg(2, mapfile)
       call getarg(3, outprefix)
       call getarg(4, val_in)
       read(val_in,*) numsamp
       call getarg(5, val_in)
       read(val_in,*) burn
    end if

    if (numsamp-burn <= 0) then
       write(*,*) 'No samples left after burn-in. Exiting'
       call mpi_finalize(ierr)
       stop
    end if

    tsys0 = 35. !Assumed Tsys for initial gain calib

    ! Find temperature (K) of Jupiter as function of frequency nu from WMAP
    nu        = [22.85d0, 33.11d0, 40.92d0, 60.41d0, 93.25d0]
    T_p       = [136.2d0, 147.2d0, 154.4d0, 165.d0,  173.d0]
    call spline(nu, T_p, 1.d30, 1.d30, T_p2)

    ! Reference values
    omega_ref = 2.481d-8
    dist_fid = 5.2d0
    f_A = 1.d0

    ! Calculate temperature at 30 GHz for 4 arcmin fwhm in uK
    initdist = 4.4d0 !obs
    initfwhm = 4.d0
    inittemp = splint(nu, T_p, T_p2, 30.d0) * 1d6
    sigma    = initfwhm / 60.d0 / sqrt(8.d0*log(2.d0)) * pi/180.d0
    omega_b  = 2.d0*pi*sigma**2
    f_d = (initdist/dist_fid)**2
    inittemp = inittemp * omega_ref/omega_b * f_A / f_d

    ! Read map information
    call read_map_h5(mapfile, inmap)
    !write(*,*) shape(inmap%m), 'shape inmap'
    nx = size(inmap%m,1)
    ny = size(inmap%m,2)
    nfreq = size(inmap%m,3)
    nsb = size(inmap%m,4)
    if (myid==0) write(*,*) nx, ny, nfreq, nsb, 'nx, ny, nfreq, nsb'
    !write(*,*) minval(inmap%x), maxval(inmap%x), 'x range'
    !write(*,*) minval(inmap%y), maxval(inmap%y), 'y range'
    !write(*,*) numsamp, '= numsamp'
    pos = maxloc(inmap%m(:,:,:,:))

    ! Set up parameter structure and initialize
    numpar = 5
    allocate(parname(numpar))
    parname(1) = 'gain'
    parname(2) = 'fwhm'
    parname(3) = 'xpos'
    parname(4) = 'ypos'
    parname(5) = 'mono'
    allocate(params(0:numsamp,numpar))  
    allocate(ipar(0:nfreq,numpar,2))  
    allocate(fpar(0:nfreq,numpar,2))  
    allocate(chisq(0:numsamp))

    ! Loop over all frequency channels
    do l = 1, nsb
       call int2string(l, sbtxt)
       ipar = 0.d0
!       do k = 20, 20
       do k = 1+myid, nfreq, numprocs
          call int2string(k, ftxt)
          write(*,*) myid, 'processing freq:', ftxt,'      sb:',sbtxt
          ! init values
          params(0,1) = inittemp/maxval(inmap%m)    !tsys0
          params(0,2) = initfwhm
          params(0,3) = inmap%x(pos(1))
          params(0,4) = inmap%y(pos(2))
          params(0,5) = 0.d0
          ! init step size
          stepsize(1) = 0.1d0 !gain
          stepsize(2) = 0.01d0 !beam arcmin
          stepsize(3) = 0.0001d0 !ra degrees
          stepsize(4) = 0.0001d0 !dec degrees
          stepsize(5) = 1000.d0  ! offset uK
          do p = 1, numpar
             if (myid==0) write(*,*) 'init ', parname(p), ' =', params(0,p), '+-', stepsize(p)
          end do

          ! Extract observed parameters
          dist = 4.4d0
          f_d = (dist/dist_fid)**2
          !nu0 = inmap%f(k,l)
          nu0 = 27.953125d0 ! GHz
          nu0 = 30.d0 ! GHz
          if (myid==0) write(*,*) nu0, '= nu0' !, inmap%f(k,l) 
          !Find intrinsic Jupiter temperature in uK
          abstemp = splint(nu, T_p, T_p2, nu0) * 1d6

          ! Loop twice, first one to initialize start points and step size
          second_time = .false.
          do while (.true.)
             numrej = 0
             chisq  = 0.d0
             ! Calculate chisq for given gain, beam size and position
             call calculate_chisq(inmap, k, l, abstemp, omega_ref, f_A, f_D, tsys0, params(0,:), chisq(0))
             ! Loop over MCMC samples
             do s = 1, numsamp
                !if (modulo(s,10 /= 0)) write(*,*) 'Processing sample', s, 'of', numsamp
                ! Draw random parameters
                do p = 1, numpar
                   params(s,p) = params(s-1,p) + rand_gauss(rng_handle) * stepsize(p)
                end do
                ! Calculate chisq for given parameters
                call calculate_chisq(inmap, k,l, abstemp, omega_ref, f_A, f_D, tsys0, params(s,:), chisq(s))
                ! Reject if necessary poor sample
                if (chisq(s) == 0.d0) then
                   ratio = 0.d0
                else if (chisq(s) < chisq(s-1)) then
                   ratio = 1.d0
                else
                   ratio = exp(-0.5*(chisq(s)-chisq(s-1)))
                end if
                if (rand_uni(rng_handle) > ratio) then
                   params(s,:) = params(s-1,:)
                   chisq(s)    = chisq(s-1)
                   numrej = numrej + 1
                end if
             end do
             
             ! Collect and write results
             write(*,*) myid, k, '  Accept rate =', real(numsamp-numrej,dp)/real(numsamp,dp)
             do p = 1, numpar
                parmean(p) = mean(params(burn:numsamp,p))
                parstdd(p) = sqrt(variance(params(burn:numsamp,p)))
                ! Check if we have data
                if (parstdd(p) /= 0.d0) then
                   if (myid==0) write(*,*) parname(p), ' =', parmean(p), '+-', parstdd(p)
                   ! write chains to file
                   if (second_time) then
                      open(14, file=trim(outprefix)//'_'//parname(p)//'_sb'//sbtxt//'_f'//ftxt//'.dat')
                      do i = 1, numsamp
                         write(14,*) i, params(i,p)
                      end do
                      close(14)
                      ipar(k,p,1) = parmean(p)
                      ipar(k,p,2) = parstdd(p)
                   else
                      stepsize(p) = parstdd(p) * 0.3
                      params(0,p) = parmean(p)
                   end if
                else
                   write(*,*) 'No data for freq:',ftxt,' sb:',sbtxt  
                end if
             end do
            
             if (second_time) then
                exit
             end if
             second_time = .true.             
          end do !while loop

       end do !freq loop
       call mpi_reduce(ipar, fpar, size(ipar), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       if (myid==0) then 
          do p = 1, numpar
             open(15, file=trim(outprefix)//'_all_'//parname(p)//'_sb'//sbtxt//'.dat')
             do i = 1, nfreq
                write(15,*) i, fpar(i,p,:) 
             end do
             close(15)
          end do
          open(15, file=trim(outprefix)//'_all_tsys_sb'//sbtxt//'.dat')
          do i = 1, nfreq
             write(15,*) i, tsys0/fpar(i,1,1) 
          end do
          close(15)
       end if
    end do !sb loop

    deallocate(params, chisq, parname, ipar, fpar)
  end subroutine estimate_gain_from_map


  !---------------------------------------------------------------------------------
  ! subroutine calculate_chisq
  !---------------------------------------------------------------------------------

  subroutine calculate_chisq(inmap, freq, sb, abstemp, omega_ref, f_A, f_D, tsys0, params, chisq)
    implicit none
   
    type(map_type)                        :: inmap 
    integer(i4b),             intent(in)  :: freq, sb
    real(dp),                 intent(in)  :: abstemp, omega_ref, f_A, f_D, tsys0
    real(dp), dimension(4),   intent(in)  :: params
    real(dp),                 intent(out) :: chisq

    
    integer(i4b)                :: i, j, nx, ny
    real(dp)                    :: r, xjup, yjup, x, y
    real(dp)                    :: tsys, fwhm, mono, gain
    real(dp)                    :: obstemp, sigma, omega_b, sigma_n
    integer(i4b), dimension(2)  :: pos

    gain = params(1)
    fwhm = params(2)
    xjup = params(3)
    yjup = params(4)
    mono = params(5)

    ! Find observed temp at given distance and beam
    sigma    = fwhm / 60.d0 / sqrt(8.d0*log(2.d0)) * pi/180.d0
    omega_b  = 2.d0*pi*sigma**2
    obstemp = abstemp * omega_ref/omega_b * f_A / f_d
    sigma_n = 10000 ! OBS should be read from file
    !write(*,*) obstemp, abstemp, 'observed and intrinsic T_jup'

    ! Make theoretical map of Jupiter for given beam, gain and position
    ! and find chisq between this and observation for each pixel.
    nx = size(inmap%m,1)
    ny = size(inmap%m,2)

    !open(16,file='inmap.dat')
    !open(17,file='residual.dat')
    do j = 1, ny
       do i = 1, nx
          if (inmap%m(i,j,freq,sb) /= 0.) then  
             x = inmap%x(i)
             y = inmap%y(j)
             r = sqrt((yjup-y)**2 + ((xjup-x)/cos(yjup*pi/180.d0))**2) * pi/180.d0
             if (r*180.d0/pi < 3.d0*fwhm/60.) then !degrees
               !write(16,*) i, j, inmap%m(i,j,freq,sb)
               !write(17,*) i, j, (exp(-r**2/sigma**2/2.d0)*obstemp*tsys0/tsys+mono-inmap%m(i,j,freq,sb))
                !chisq  = chisq + (exp(-r**2/sigma**2/2.d0)*obstemp*tsys0/tsys+mono-inmap%m(i,j,freq,sb))**2 / sigma_n**2
                chisq  = chisq + (exp(-r**2/sigma**2/2.d0)*obstemp*gain+mono-inmap%m(i,j,freq,sb))**2 / sigma_n**2
             end if
          end if
       end do
    end do
    !close(16)
    !close(17)
    !write(*,*) 'chisq =', chisq

  end subroutine calculate_chisq

  !----------------------------------------------------------------------------------
  ! subroutine give_user_info
  !----------------------------------------------------------------------------------

  subroutine give_user_info
    implicit none

    if (myid == root) then
       write(*,*) "Usage: 'maptool sindex parameters', where parameters are:" 
       write(*,*) '(maps, beams amd mask are fits-files. freqs, fwhm and nside are numbers)' 
       write(*,*) 'map1 map2 beam1 beam2 freq1 freq2 fwhm mask rms1 rms2 nside outprefix (uncertlimit)' 
       write(*,*) ''
       write(*,*) 'scatter map1 freq1 map2 freq2 len bre radius outprefix'
       write(*,*) ''
       write(*,*) 'clusim clusters fullsim cmb beams noise cov nside outprefix at ak cmbseed noiseseed'
       write(*,*) ''
       write(*,*) 'transfer inmap outmap component outprefix'
       write(*,*) ''
       write(*,*) 'merge numfiles infiles outfile'
       write(*,*) ''
       write(*,*) 'extract numcomps infile comps outfile'
       write(*,*) ''
       write(*,*) 'mdmaps numband mdtxtfile nside outprefix'
       write(*,*) ''
       write(*,*) 'sims powspec beam pixwin mask nside numside outdir numsim (lmax)'
       write(*,*) ''
       write(*,*) 'beta weightmap, map1, val2, outfile lowcut highcut'
       write(*,*) ''
    end if
    call mpi_finalize(ierr)
    stop

  end subroutine give_user_info



end program maptool
