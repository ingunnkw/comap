program tester
   use quiet_ces_mod
   use quiet_utils
   use quiet_hdf_mod
   use quiet_fft_mod
   use ziggurat
   use rngmod
   use quiet_filter_mod
   use quiet_fft_mod
   use quiet_lx_mod
   use quiet_gain_mod
   use quiet_pointing_mod
   use quiet_target_mod
   use powell_mod
   use quiet_gain_mod
   use quiet_pointing_mod2
   use quiet_covfile_mod
   use quiet_mapfile_mod
   use quiet_pixspace_mod
   use quiet_angcorr_mod
   use quiet_mpi_mod
   use quiet_ephem_mod
   use quiet_noise_estimation_mod
   use quiet_healpix_mod
   use alm_tools
   implicit none

   real(dp)              :: powell_dk_real(4)
   integer(i4b)          :: powell_dk_int(3)
   character(len=512)    :: cmd
   real(dp), allocatable :: powell_data(:,:)

   type obs_type
      integer(i4b) :: ces, di, horn
      real(dp)     :: mjd, p0(3), p(3), tamp, amp, frac, dist, fwhm, foo(6), ddk
   end type

   ! test

   call getarg(1, cmd)
   select case(cmd)
      !case("indexsmooth"); call indexsmooth
      !case("sample_pol_multipole"); call sample_pol_multipole
      !case("print_fg_spectra"); call print_fg_spectra
      case("chisq_test"); call chisq_test
!!$      case("rdhdf");    call rdhdf
!!$      case("rdfits");   call rdfits
!!$      case("mkhdf");    call mkhdf
!!$      case("mkfits");   call mkfits
!!$      case("tcorr");    call tcorr
!!$      case("bincum");   call bincum
!!$      case("actest");   call actest
!!$      case("mktemp");   call mktemp
!!$      case("cmbfill");  call cmb_fill_brute
!!$      case("cov2map");  call cov2map
!!$      case("metro");    call metroprior
!!$      case("tsim");     call tsim
!!$      case("wnoise");   call wnoise
!!$      case("tfit");     call tfit
!!$      case("covtest");  call covtest4
!!$      case("covdump");  call covdump
!!$      case("fit_gauss");call fit_gauss
!!$      case("suncomb");  call suncomb
!!$      case("raul");     call raul
!!$      case("rstsigma"); call rstsigma
!!$      case("scorr");    call scorr
!!$      case("bleh");     call bleh
!!$      case("radconv");  call radconv
!!$      case("gainres");  call gainres
!!$      case("getgains"); call getgains2
!!$      case("eigtime");  call eigtime
!!$      case("decon");    call deconvolution
!!$      case("foo");      call foo
!!$      case("est");      call noise_est
!!$      case("mpi");      call mpitest
!!$      case("mpi2");     call anothermpitest
!!$      case("point");    call pointtest4
!!$      case("solve");    call covtest2
!!$      case("test");     call covtest3
!!$      case("r2c");      call r2c

contains

  subroutine chisq_test
    implicit none

    type(planck_rng) :: handle
    integer(i4b)     :: n, i, j, l, m
    real(dp)         :: cl, cl_theory, dl, dl_theory
    real(dp), dimension(100) :: s

!!$    n         = 1000000
!!$    l         = 2
!!$    dl_theory = 1100.d0
!!$    cl_theory = dl_theory / (l*(l+1)/2.d0/pi)
!!$    call rand_init(handle, 3741)
!!$    open(58,file='samples.dat')
!!$    do i = 1, n
!!$       cl = 0.d0
!!$       do m = 1, 2*l+1
!!$          cl = cl + cl_theory * rand_gauss(handle)**2
!!$       end do
!!$       cl = cl / (2*l+1)
!!$       dl = cl * l*(l+1)/2.d0/pi
!!$       write(58,*) i, dl
!!$    end do
!!$    close(58)

    n = 100000
    m = 100
    call rand_init(handle, 3741)
    open(58,file='samples.dat')
    do i = 1, n
       do j = 1, m
          s(j) = rand_gauss(handle)
       end do
       write(58,*) i, mean(s(1:78)*sqrt(78.d0)), mean(s(1:93)*sqrt(93.d0))
    end do
    close(58)

    
  end subroutine chisq_test
  
  subroutine sample_pol_multipole
    implicit none

    integer(i4b), parameter :: nmaps = 3
    integer(i4b) :: i, j, k, n, l, ierr
    real(dp)     :: C(nmaps,nmaps), sqrtC(nmaps,nmaps), eta(nmaps,1), cl(nmaps,nmaps)
    type(planck_rng) :: handle

    l = 2
    n = 1000000
    call rand_init(handle, 9841)

    C      = 0.d0
    C(1,1) = 1000.d0
    if (nmaps == 3) then
       C(1,2) = 1.d0
       C(2,1) = C(1,2)
       C(2,2) = 0.1d0
       C(3,3) = 0.01d0
    end if
    call cholesky_decompose(C, sqrtC)

    open(58,file='samples.dat')
    do i = 1, n
       cl = 0.d0
       do k = 1, 2*l+1
          do j = 1, nmaps
             eta(j,1) = rand_gauss(handle)
          end do
          eta = matmul(sqrtC, eta)
          cl = cl + matmul(eta,transpose(eta))
       end do
       cl = cl / (2*l+1)
       if (nmaps == 1) then
          write(58,fmt='(i10,4f16.8)') i, cl(1,1)
       else
          write(58,fmt='(i10,4f16.8)') i, cl(1,1), cl(1,2), cl(2,2), cl(3,3)
       end if
    end do
    close(58)

  end subroutine sample_pol_multipole

  subroutine indexsmooth
    implicit none
    
    integer(i4b)     :: numband, n, i, j, k, m, numsim
    real(dp)         :: A, A_rms, beta, beta_rms, fwhm, sigma, dx, xmin, xmax, A0, beta0, r, beam, w
    real(dp)         :: A_ref, beta_ref
    type(planck_rng) :: handle
    real(dp), allocatable, dimension(:)     :: x, nu_log, nu, f_smooth
    real(dp), allocatable, dimension(:,:)   :: bias
    real(dp), allocatable, dimension(:,:,:) :: f
    
    numband = 6
    fwhm    = 30d0
    xmin    = -90.d0
    xmax    =  90.d0
    n       = 11
    dx      = (xmax-xmin) / (n-1)
    sigma   = fwhm / sqrt(8.d0*log(2.d0))
    numsim  = 1000

    A        = 100.d0
    A_rms    = 10.d0
    beta     = -3.d0
    beta_rms = 0.2d0

    call rand_init(handle, 481141)

    allocate(x(n), f(n,n,numband), nu(numband), f_smooth(numband), bias(numband,numsim))
    nu     = [30.d0, 44.d0, 70.d0, 100.d0, 143.d0, 217.d0]
    do i = 1, n
       x(i) = xmin + (i-1)*dx
    end do

    open(58,file='smoothspec.dat')    
    do m = 1, numsim

       ! Set up high-resolution field
       A_ref    = 0.d0
       beta_ref = 0.d0
       w        = 0.d0
       do i = 1, n
          do j = 1, n
             r    = sqrt(x(i)**2 + x(j)**2)
             beam = exp(-0.5*(r/sigma)**2)
             w    = w + beam

             A0    = A    + A_rms    * rand_gauss(handle)
             beta0 = beta + beta_rms * rand_gauss(handle)

             A_ref    = A_ref    + A0    * beam
             beta_ref = beta_ref + beta0 * beam
             do k = 1, numband
                f(i,j,k) = A0 * (nu(k) / nu(1))**beta0
             end do
          end do
       end do
       A_ref    = A_ref / w
       beta_ref = beta_ref / w
       
       ! Do the smoothing
       f_smooth = 0.d0
       w        = 0.d0
       do i = 1, n
          do j = 1, n
             r    = sqrt(x(i)**2 + x(j)**2)
             beam = exp(-0.5*(r/sigma)**2)
             w    = w + beam
             do k = 1, numband
                f_smooth(k) = f_smooth(k) + f(i,j,k) * beam
             end do
          end do
       end do
       f_smooth = f_smooth / w
       
       ! Fit a power law
       beta0    = (numband * sum(log(nu)*log(f_smooth)) - sum(log(nu))*sum(log(f_smooth))) / &
            & (numband * sum(log(nu)**2) - (sum(log(nu)))**2)
       A0       = (sum(log(f_smooth)) - beta0*sum(log(nu))) / numband
       write(*,*) m, nu(1)**beta0 * exp(A0), beta0

       do k = 1, numband
          bias(k,m) = f_smooth(k) / (A_ref * (nu(k)/nu(1))**beta_ref)
          !write(58,*) nu(k), f_smooth(k)
          !write(58,*) nu(k), f_smooth(k) / (A * (nu(k)/nu(1))**beta)
       end do
       !write(58,*) 

    end do
    close(58)
       
    open(58,file='bias.dat')
    do k = 1, numband
       write(58,*) nu(k), mean(bias(k,:)), sqrt(variance(bias(k,:)))
    end do
    close(58)

    ! Output spectra

!!$    do i = 1, n
!!$       do j = 1, n
!!$       end do
!!$    end do

!!$    do k = 1, numband
!!$       write(58,*) nu(k), f_smooth(k)
!!$    end do
!!$    close(58)

  end subroutine indexsmooth

!!$   subroutine mkfits
!!$     implicit none
!!$     character(len=512)             :: fname, arg
!!$     character(len=16)              :: type
!!$     character(len=20), allocatable :: ttype(:), tform(:), tunit(:)
!!$     real(dp),          allocatable :: data(:,:)
!!$     integer(i4b)                   :: i, j, n, m, status, unit
!!$     type(planck_rng)               :: rng
!!$     type = "normal"
!!$     call rand_init(rng, 1)
!!$     call getarg(2, arg); read(arg,*) n
!!$     call getarg(3, arg); read(arg,*) m
!!$     call getarg(4, fname)
!!$     if(iargc() > 4) call getarg(5, type)
!!$     allocate(data(n,m))
!!$     do i = 1, m
!!$        do j = 1, n
!!$           data(j,i) = rand_gauss(rng)
!!$        end do
!!$        write(*,*) i
!!$     end do
!!$     status = 0
!!$     unit = getlun()
!!$     call ftinit(unit,"!"//trim(fname),1,status)
!!$     if(type == "primary") then
!!$        call ftiimg(unit,-64,2,[n,m], status)
!!$        call ftpprd(unit,0,1,n*m,data, status)
!!$        call ftpdat(unit,status) ! format (yyyy-mm-dd)
!!$     else
!!$        call ftphps(unit,32,0,0,status)
!!$        call ftcrhd(unit,status)
!!$        if (status > 0) call printerror(status)
!!$        allocate(ttype(m), tform(m), tunit(m))
!!$        ttype = 'unknown'; tform = '1D'; tunit = ''
!!$        call ftphbn(unit, n, m, ttype, tform, tunit, '', 0, status)
!!$        if (status > 0) call printerror(status)
!!$        deallocate(ttype, tform, tunit)
!!$        if(type == "fast") then
!!$           call ftptbb(unit, 1, 1, n*m*8, data, status)
!!$           if (status > 0) call printerror(status)
!!$        else
!!$           do i = 1, m
!!$              call ftpcld(unit, i, 1, 1, n, data(:,i), status)
!!$              if (status > 0) call printerror(status)
!!$           end do
!!$        end if
!!$     end if
!!$     call ftclos(unit, status)
!!$     if (status > 0) call printerror(status)
!!$   end subroutine
!!$
!!$   subroutine mkhdf
!!$     implicit none
!!$     character(len=512)             :: fname, arg
!!$     real(dp),          allocatable :: data(:,:)
!!$     integer(i4b)                   :: i, j, n, m
!!$     type(planck_rng)               :: rng
!!$     type(hdf_file)                 :: hfile
!!$     call rand_init(rng, 1)
!!$     call getarg(2, arg); read(arg,*) n
!!$     call getarg(3, arg); read(arg,*) m
!!$     call getarg(4, fname)
!!$     allocate(data(n,m))
!!$     do i = 1, m
!!$        do j = 1, n
!!$           data(j,i) = rand_gauss(rng)
!!$        end do
!!$        write(*,*) i
!!$     end do
!!$     call open_hdf_file("test.hdf", hfile, "w")
!!$     call write_hdf(hfile, "data", data)
!!$     call close_hdf_file(hfile)
!!$   end subroutine
!!$
!!$   subroutine rdfits
!!$     implicit none
!!$     character(len=512)             :: fname
!!$     character(len=20)              :: comment, arg
!!$     character(len=8)               :: type
!!$     integer(i4b)                   :: status, n, m, i, unit
!!$     logical(lgt)                   :: anynull
!!$     real(dp),          allocatable :: data(:,:)
!!$     type = "normal"
!!$     unit   = getlun()
!!$     status = 0
!!$     call getarg(2, fname)
!!$     if(iargc() > 2) call getarg(3, type)
!!$     if(type == "primary") then
!!$        call ftnopn(unit,trim(fname),0,status)
!!$        if (status > 0) call printerror(status)
!!$        call ftgkyj(unit,"NAXIS1",  n,  comment,status)
!!$        if (status > 0) call printerror(status)
!!$        call ftgkyj(unit,"NAXIS2",  m,  comment,status)
!!$        if (status > 0) call printerror(status)
!!$        allocate(data(n,m))
!!$        call ftgpvd(unit, 0, 1, n*m, 0d0, data, anynull, status)
!!$     else
!!$        call ftnopn(unit,trim(fname) // "[1]",0,status)
!!$        if (status > 0) call printerror(status)
!!$        call ftgkyj(unit,"NAXIS2",  n,  comment,status)
!!$        if (status > 0) call printerror(status)
!!$        call ftgkyj(unit,"TFIELDS", m,  comment,status)
!!$        if (status > 0) call printerror(status)
!!$        allocate(data(n,m))
!!$        if(type == "fast") then
!!$           call ftgtbb(unit, 1, 1, n*m*8, data, status)
!!$        else
!!$           do i = 1, m
!!$              call ftgcvd(unit, i, 1, 1, n, 0d0, data(:,i), anynull, status)
!!$              if (status > 0) call printerror(status)
!!$           end do
!!$        end if
!!$     end if
!!$     call ftclos(unit, status)
!!$     if (status > 0) call printerror(status)
!!$   end subroutine
!!$
!!$   subroutine rdhdf
!!$     implicit none
!!$     character(len=512)             :: fname
!!$     real(dp),          allocatable :: data(:,:)
!!$     type(hdf_file)                 :: hfile
!!$     call getarg(2, fname)
!!$     call open_hdf_file(fname, hfile, "r")
!!$     call read_alloc_hdf(hfile, "data", data)
!!$     call close_hdf_file(hfile)
!!$   end subroutine
!!$
!!$   subroutine tcorr
!!$     implicit none
!!$     character(len=32)         :: arg, arg2, typ
!!$     integer(i4b)              :: n, m, i, ai, nc
!!$     real(dp)                  :: srate, f, alpha, fk, asign, w, f1, f2, o, s
!!$     real(dp),   allocatable   :: tc(:), flt(:)
!!$     srate = 25
!!$     call getarg(2, arg); read(arg,*) n
!!$     m = n/2+1
!!$     allocate(flt(m))
!!$     flt = 1
!!$     ai = 2
!!$     do while(ai < iargc())
!!$        ai = ai+1
!!$        call getarg(ai, arg)
!!$        select case(arg)
!!$           case("oof")
!!$              ai = ai+1; call getarg(ai, arg2); read(arg2,*) fk
!!$              ai = ai+1; call getarg(ai, arg2); read(arg2,*) alpha
!!$              asign = 1; if(alpha < 0) asign = -1
!!$              do i = 1, m
!!$                 f = ind2freq(i, srate, m)
!!$                 flt(i) = flt(i)*(1+(f/fk)**alpha)**-1
!!$              end do
!!$           case("cos")
!!$              ai = ai+1; call getarg(ai, arg2); read(arg2,*) fk
!!$              ai = ai+1; call getarg(ai, arg2); read(arg2,*) w
!!$              f1 = fk-abs(w/2); f2 = fk+abs(w/2)
!!$              s = -1
!!$              if(w < 0) s = 1
!!$              do i = 1, m
!!$                 f = ind2freq(i, srate, m)
!!$                 if(f < f1) then
!!$                    flt(i) = flt(i)*(1+s)/2
!!$                 else if(f < f2) then
!!$                    flt(i) = flt(i)*(1+s*cos(pi*(f-f1)/(f2-f1)))/2
!!$                 else
!!$                    flt(i) = flt(i)*(1-s)/2
!!$                 end if
!!$              end do
!!$           case("sharp")
!!$              ai = ai+1; call getarg(ai, arg2); read(arg2,*) fk
!!$              ai = ai+1; call getarg(ai, arg2); read(arg2,*) s
!!$              i = freq2ind(fk, srate, m)
!!$              if(s > 0) then
!!$                 flt(i:) = 0
!!$              else
!!$                 flt(:i) = 0
!!$              end if
!!$        end select
!!$     end do
!!$     call get_inv_noise_filter_real(n, flt, nc, tc)
!!$     do i = 0, n-1
!!$        write(*,'(2e15.7)') i/srate, tc(i)
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine bincum
!!$     implicit none
!!$     character(len=512)        :: lfile, ifile, ofile, arg, opre
!!$     integer(i4b)              :: unit, fi, nside, order, ncomp, nmap, ext(3)
!!$     integer(i4b)              :: npix, nfile, i, j, k, c, m, status
!!$     integer(i4b), allocatable :: pixlist(:,:), comps(:)
!!$     real(dp),     allocatable, dimension(:,:,:) :: rhs, rhst, div, divt, bin, obin
!!$     type(hdf_file)            :: hfile
!!$     call getarg(2, lfile)
!!$     call getarg(3, opre)
!!$
!!$     m = 1 ! Which map/sim to use
!!$
!!$     unit = getlun()
!!$     ! Read the first file to get the dimensions. All files
!!$     ! will have the same pixels.
!!$     open(unit, file=lfile, action="read", status="old")
!!$     read(unit,'(a)') ifile
!!$     call open_hdf_file(ifile, hfile, "r")
!!$     call read_hdf(hfile, "nside", nside)
!!$     call read_hdf(hfile, "order", order)
!!$     call read_alloc_hdf(hfile, "comps", comps)
!!$     call read_alloc_hdf(hfile, "pixlist", pixlist)
!!$     call get_size_hdf(hfile, "rhs", ext)
!!$     ncomp = size(comps)
!!$     nmap  = ext(3)
!!$     npix  = ext(1)
!!$     allocate(rhs (npix,ncomp,nmap), div (npix,ncomp,ncomp))
!!$     allocate(rhst(npix,ncomp,nmap), divt(npix,ncomp,ncomp))
!!$     allocate(bin(npix,ncomp,nmap))
!!$     allocate(obin(size(pixlist,2),3,nmap))
!!$     rhst = 0; divt = 0
!!$     call close_hdf_file(hfile)
!!$     rewind(unit)
!!$
!!$     ! Then loop through all the files, accumulating data
!!$     fi = 0
!!$     do
!!$        read(unit,'(a)', end=1) ifile
!!$        fi = fi + 1
!!$        call dmem("Reading " // trim(ifile))
!!$        call open_hdf_file(ifile, hfile, "r")
!!$        call read_hdf(hfile, "rhs", rhs)
!!$        call read_hdf(hfile, "div", div)
!!$        call close_hdf_file(hfile)
!!$        rhst = rhst + rhs
!!$        divt = divt + div
!!$        ! Solve the current map
!!$        do m = 1, nmap
!!$           do i = 1, npix
!!$              if(divt(i,1,1) == 0) then
!!$                 bin(i,:,m) = hpx_dbadval
!!$              else
!!$                 call solve_system_eiglim(divt(i,:,:), rhst(i,:,m), bin(i,:,m), 1d-12, status)
!!$                 if(status /= 0) bin(i,:,m) = hpx_dbadval
!!$              end if
!!$           end do
!!$        end do
!!$        ! Flatten bin into the proper output format
!!$        obin = hpx_dbadval
!!$        do c = 1, ncomp
!!$           obin(:,comps(c),:) = bin(pixlist(1,:),c,:)
!!$        end do
!!$        ! And write
!!$        ofile = trim(opre) // trim(itoa(fi,3)) // ".hdf"
!!$        call write_map(obin, pixlist(2,:), nside, order, ofile)
!!$     end do
!!$     1 continue
!!$     deallocate(rhs, div, rhst, divt, bin, obin, pixlist, comps)
!!$   end subroutine
!!$
!!$   subroutine actest
!!$     implicit none
!!$     character(len=512)        :: ifile, arg
!!$     integer(i4b)              :: i, l, k, m, n, lmax, unit, nc
!!$     real(dp)                  :: bsize, C(3,3), C2(3,3), v(3,2), foo(6), rmax,r
!!$     real(dp),   allocatable   :: pslin(:,:), psmat(:,:,:), beam(:,:), corr(:)
!!$     type(angcorr_type)        :: cfun
!!$     call getarg(2, ifile)
!!$     call getarg(3, arg); read(arg, *) bsize
!!$     m = 10000
!!$     rmax = pi
!!$
!!$     unit = getlun()
!!$     open(unit,file=ifile,action="read",status="old")
!!$     lmax = 0
!!$     i    = 0
!!$     do
!!$        if(i == 0) then
!!$           read(unit,'(a)',end=1) arg
!!$           nc = num_tokens(arg, " ")-1
!!$        else
!!$           read(unit,*,end=1) l
!!$           if(l > lmax) lmax = l
!!$        end if
!!$        i=i+1
!!$     end do
!!$     1 rewind(unit)
!!$     allocate(pslin(0:lmax,nc),psmat(0:lmax,3,3),beam(0:lmax,3))
!!$     pslin = 0
!!$     do
!!$        read(unit,*,end=2) l, foo(:nc)
!!$        pslin(l,:) = foo(:nc)*2*pi/(l*(l+1))
!!$     end do
!!$     2 close(unit)
!!$     bsize = bsize*pi/180/60/sqrt(8*log(2d0))
!!$     do l = 0, lmax
!!$        beam(l,:) = exp(-0.5*l*(l+1)*bsize**2)
!!$     end do
!!$     call build_powspec_matrix(pslin, beam, psmat)
!!$     call angcorr_init(cfun, psmat, lmin=2)
!!$     v(:,1) = [0,0,1]
!!$     do i = 0, m-1
!!$        r = rmax*i/m
!!$        call ang2vec(r,pi/2,v(:,2))
!!$        !call angcorr_get(cfun, v(:,1), v(:,2), C)
!!$        call angcorr_get(cfun, cos(r), C)
!!$        call angcorr_get_raw(cfun, cos(r), C2)
!!$        write(*,'(f12.6,4e15.7)') r, C(1,1), C2(1,1), C(1,1)-C2(1,1), (C(1,1)-C2(1,1))/C2(1,1)
!!$     end do
!!$     call angcorr_free(cfun)
!!$   end subroutine
!!$
!!$   subroutine build_powspec_matrix(powspec, beam, ps)
!!$     implicit none
!!$     real(dp)     :: powspec(:,:), beam(:,:), ps(:,:,:)
!!$     integer(i4b) :: i, j, k
!!$     do i = 1, size(powspec,1)
!!$        call convert_lin2mat2(powspec(i,:), ps(i,:,:))
!!$     end do
!!$     do j = 1, size(ps,3)
!!$       do i = j, size(ps,2)
!!$          ps(:,i,j) = ps(:,i,j) * beam(:,i) * beam(:,j)
!!$          ps(:,j,i) = ps(:,i,j)
!!$       end do
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine mktemp
!!$     implicit none
!!$     character(len=512)        :: ifile, arg, ofile, cfile
!!$     integer(i4b)              :: nside, order, npix, n, ncomp, oncomp, i, j
!!$     real(dp)                  :: rad, sigma
!!$     real(dp),     allocatable :: imap(:,:), omap(:,:), ocov(:,:)
!!$     integer(i4b), allocatable :: pixels(:)
!!$     type(planck_rng)          :: rng
!!$     call getarg(2,ifile)
!!$     call getarg(3,arg);  read(arg,*) rad
!!$     call getarg(4,arg);  read(arg,*) sigma
!!$     call getarg(5,arg);  read(arg,*) oncomp
!!$     call getarg(6,ofile)
!!$     call getarg(7,cfile)
!!$     call rand_init(rng,1)
!!$     call read_map(imap, order, ifile)
!!$     npix  = size(imap,1)
!!$     ncomp = size(imap,2)
!!$     nside = npix2nside(npix)
!!$     allocate(pixels(npix))
!!$     call query_disc(nside, [1d0,0d0,0d0], rad, pixels, n, nest=1)
!!$     allocate(omap(n,ncomp),ocov(n*oncomp,n*oncomp))
!!$     omap(:,:) = imap(pixels(:n),:)*1000
!!$     do j = 1, ncomp
!!$        do i = 1, n
!!$           omap(i,j) = omap(i,j) + rand_gauss(rng)*sigma
!!$        end do
!!$     end do
!!$     ocov = 0
!!$     do i = 1, n*oncomp
!!$        ocov(i,i) = sigma**-2
!!$     end do
!!$     call write_map(omap, pixels(:n), nside, order, ofile)
!!$     call write_covmat(ocov, pixels(:n), nside, order, oncomp, cfile)
!!$   end subroutine
!!$
!!$   subroutine cmb_fill_brute
!!$     implicit none
!!$     character(len=512)        :: arg, psfile, mapfile, maskfile, ofile
!!$     integer(i4b)              :: nside, order, nside2, order2, nsim, i, j, k, n
!!$     real(dp),     allocatable :: map(:,:), mask(:,:),ps(:,:), mu(:), chol(:,:), omap(:,:)
!!$     real(dp),     allocatable :: noise(:)
!!$     integer(i4b), allocatable :: pixels(:), mpixels(:), inds(:)
!!$     type(planck_rng)          :: rng
!!$     type(hdf_file)            :: hfile
!!$     call getarg(2,mapfile)
!!$     call getarg(3,maskfile)
!!$     call getarg(4,psfile)
!!$     call getarg(5,arg); read(arg,*) nsim
!!$     call getarg(6,ofile)
!!$
!!$     call dmem("start")
!!$     call rand_init(rng, 1)
!!$
!!$     call read_map(map,  pixels,  nside,  order,  mapfile)
!!$     call read_map(mask, order2, maskfile)
!!$     call read_powspec_ascii(psfile, ps)
!!$     call dmem("read")
!!$     nside2 = npix2nside(size(mask,1))
!!$
!!$     call assert(nside2 == nside, "Inconsistent nside between map and mask!")
!!$     call assert(order2 == order, "Inconsistent order between map and mask!")
!!$
!!$     ! Convert mask into inds. Inds is the list of pixels to be replaced.
!!$     ! These are given by pixels present and nonzero in the mask.
!!$     call wherei(mask(pixels,1) /= 0 .and. healok(mask(pixels,1)), inds)
!!$
!!$     ! Ok, we are ready to make the simulations
!!$     n = size(inds)
!!$     allocate(mu(n), chol(n,n), omap(size(map,1),nsim), noise(n))
!!$     call dmem("mask")
!!$     call cmb_fill_brute_params(map(:,1), pixels, nside, order, inds, ps(:,1), mu, chol)
!!$call open_hdf_file("foo.hdf", hfile, "w")
!!$call write_hdf(hfile, "mu", mu)
!!$call write_hdf(hfile, "chol", chol)
!!$call write_hdf(hfile, "ps", ps(:,1))
!!$call close_hdf_file(hfile)
!!$
!!$     call dmem("prepare fill")
!!$     do i = 1, nsim
!!$        do j = 1, n
!!$           noise(j) = rand_gauss(rng)
!!$        end do
!!$        omap(:,i) = map(:,1)
!!$        omap(inds,i) = mu + matmul(chol,noise)
!!$     end do
!!$     call dmem("simulate")
!!$     call write_map(omap, pixels, nside, order, ofile)
!!$     call dmem("write")
!!$   end subroutine
!!$
!!$   subroutine cmb_fill_brute_params(map, pixels, nside, order, inds, ps, mu, chol)
!!$     implicit none
!!$     real(dp)              :: map(:), ps(0:), mu(:), chol(:,:)
!!$     integer(i4b)          :: pixels(:), inds(:), nside, order, n, m, i, j, k, lmax, ns
!!$     real(dp)              :: dist, v(3,2), dr, r
!!$     real(dp), allocatable :: ips(:), iC11(:,:), iC12(:,:), xy(:,:), vecs(:,:)
!!$     integer(i4b), allocatable :: uinds(:), tpix(:)
!!$     type(spline_type)     :: s
!!$     m   = size(inds)
!!$     n   = size(pixels)-m
!!$     lmax= ubound(ps,1)
!!$     allocate(ips(0:lmax), iC11(m,m), iC12(m,n), tpix(size(pixels)))
!!$ps = ps + 1e-3
!!$     ips = 0; where(ps /= 0) ips = 1/ps*(4*pi/(12*nside**2))**2
!!$
!!$     ! Build the set of unmasked pixels
!!$     tpix = pixels
!!$     tpix(inds) = -1
!!$     call wherei(tpix >= 0, uinds)
!!$     call assert(size(uinds) == n, "Duplicate pixels in cmb_fill_brute?")
!!$
!!$     ! Computing the correlation function is expensive, so we must
!!$     ! spline it.
!!$     dr = pi/lmax
!!$     ns = ceiling(pi/dr)*5
!!$     allocate(xy(ns,2))
!!$     do i = 1, ns
!!$        xy(i,1) = (i-1)*pi/(ns-1)
!!$     end do
!!$     call cl2corr(ips, xy(:,2), xy(:,1))
!!$     call spline(s, xy(:,1), xy(:,2), regular=.true.)
!!$     call dmem("build spline")
!!$
!!$!do i = 1, ns
!!$!r = (i-1)*dr/100
!!$!xy(1,1) = r
!!$!call cl2corr(ips, xy(1:1,2), xy(1:1,1))
!!$!write(*,'(f15.10,2e15.7)') r, splint(s, r), xy(1,2)
!!$!end do
!!$!stop
!!$
!!$     allocate(vecs(3,size(pixels)))
!!$     do i = 1, size(pixels)
!!$        vecs(:,i) = pix2vec(nside, order, pixels(i))
!!$     end do
!!$     call dmem("build vecs")
!!$
!!$     ! First build iC11 and iC12 from ips O(nm)
!!$     do i = 1, m
!!$        v(:,1) = vecs(:,inds(i))
!!$        do j = i, m
!!$           v(:,2) = vecs(:,inds(j))
!!$           call angdist(v(:,1), v(:,2), dist)
!!$           iC11(i,j) = splint(s, dist)
!!$           iC11(j,i) = iC11(i,j)
!!$        end do
!!$        do j = 1, n
!!$           v(:,2) = vecs(:,uinds(j))
!!$           call angdist(v(:,1), v(:,2), dist)
!!$           iC12(i,j) = splint(s, dist)
!!$        end do
!!$        call dmem("build " // trim(itoa(i)))
!!$     end do
!!$     call free_spline(s)
!!$
!!$     call eigen_pow(iC11, -0.5d0, chol)
!!$     mu = -matmul(matmul(chol,chol), matmul(iC12, map(uinds)))
!!$   end subroutine
!!$
!!$   ! Read a covariance matrix, and produce
!!$   ! maps for the diagonal and all vertical slices.
!!$   subroutine cov2map
!!$     implicit none
!!$     character(len=512)         :: arg, ifile, ofile
!!$     integer(i4b),  allocatable :: pixels(:)
!!$     real(dp),      allocatable :: cov(:,:,:,:), map(:,:,:)
!!$     integer(i4b)               :: nside, order, n, nc, ns, i, j1, j2, j, k, k2
!!$     call getarg(2, ifile)
!!$     call getarg(3, arg); read(arg,*) ns
!!$     call getarg(4, ofile)
!!$     call read_covmat(cov, pixels, nside, order, ifile)
!!$     n  = size(cov,1)
!!$     nc = size(cov,2)
!!$     allocate(map(n,nc**2,ns))
!!$     do i = 1, ns
!!$        do j1 = 1, nc
!!$           do j2 = 1, nc
!!$              j = (j1-1)*nc+j2
!!$              do k = 1, n
!!$                 k2 = i-1
!!$                 if(k2 == 0) k2 = k
!!$                 map(k,j,i) = cov(k,j1,k2,j2)
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$     call write_map(map, pixels, nside, order, ofile)
!!$   end subroutine
!!$
!!$   subroutine metroprior
!!$     implicit none
!!$     real(dp) :: p, p2, s
!!$     type(planck_rng) :: rng
!!$     call rand_init(rng, 1)
!!$     p = 0; s = 1
!!$     do
!!$        p2 = p + rand_gauss(rng)*s
!!$        if(abs(p2) > 10) cycle
!!$        if(rand_uni(rng) > 0.5 .and. abs(p2) < 10) p = p2
!!$        write(*,'(f15.10)') p
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine tsim
!!$     implicit none
!!$     character(len=512)             :: fname
!!$     real(dp),        allocatable   :: map(:,:,:,:)
!!$     integer(i4b),    allocatable   :: pixels(:)
!!$     integer(i4b)                   :: n, nside, order, ncomp, i, j, k
!!$     real(dp)                       :: A, B, S
!!$     type(planck_rng)               :: rng
!!$     call rand_init(rng, 5545)
!!$     S = 1
!!$     A = 1.3333333
!!$     B = 3
!!$     n = 12
!!$     nside = 1
!!$     order = 1
!!$     allocate(pixels(n), map(n,3,2,2))
!!$     pixels = irange(n)-1
!!$     ! Noise properties
!!$     map(:n/2,:,2,1) = 10
!!$     map(n/2+1:,:,2,1) = 40
!!$     map(:n/3,:,2,2) = 5
!!$     map(n/3+1:,:,2,2) = 100
!!$     ! Generate signal
!!$     do i = 1, n
!!$        do j = 1, size(map,2)
!!$           map(i,j,1,2) = rand_gauss(rng)*S
!!$        end do
!!$     end do
!!$     map(:,2,1,2) = map(:,2,1,2) * 3 ! Make Q and U different
!!$     ! Set up signal relation
!!$     map(:,:,1,1) = A*map(:,:,1,2)+B
!!$     ! Add noise
!!$     do i = 1, n
!!$        do j = 1, size(map,2)
!!$           do k = 1, 2
!!$              map(i,j,1,k) = map(i,j,1,k) + rand_gauss(rng)*map(i,j,2,k)
!!$           end do
!!$        end do
!!$     end do
!!$
!!$     call write_map(map(:,:,1,1), pixels, nside, order, "sim_map1.hdf")
!!$     call write_map(map(:,:,1,2), pixels, nside, order, "sim_map2.hdf")
!!$     call write_map(map(:,:,2,1), pixels, nside, order, "sim_rms1.hdf")
!!$     call write_map(map(:,:,2,2), pixels, nside, order, "sim_rms2.hdf")
!!$
!!$   end subroutine
!!$
!!$   subroutine tfit
!!$     implicit none
!!$     character(len=512)             :: arg, fname
!!$     real(dp),        allocatable   :: imap(:,:), map(:,:), var(:,:), dmap(:), dvar(:)
!!$     real(dp),        allocatable   :: res(:,:), res2(:,:,:)
!!$     integer(i4b),    allocatable   :: pixels(:)
!!$     integer(i4b)                   :: nside, npix, order, ncomp, i, j, k, n, m
!!$     integer(i4b)                   :: r(2)
!!$     real(dp)                       :: a, b, L, fudge
!!$
!!$     r = [2, 2]
!!$     call getarg(2, fname)
!!$     call read_map(imap, pixels, nside, order, fname)
!!$     npix = size(pixels)
!!$     n    = (r(2)-r(1)+1)*npix
!!$     allocate(map(n,2),var(n,2))
!!$     map(:,1) = reshape(imap(:,r(1):r(2)),[n])
!!$     call getarg(3,fname)
!!$     deallocate(imap, pixels)
!!$     call read_map(imap, pixels, nside, order, fname)
!!$     var(:,1) = reshape(imap(:,r(1):r(2)),[n])**2
!!$     call getarg(4,fname)
!!$     deallocate(imap, pixels)
!!$     call read_map(imap, pixels, nside, order, fname)
!!$     map(:,2) = reshape(imap(:,r(1):r(2)),[n])
!!$     call getarg(5,fname)
!!$     deallocate(imap, pixels)
!!$     call read_map(imap, pixels, nside, order, fname)
!!$     var(:,2) = reshape(imap(:,r(1):r(2)),[n])**2
!!$     deallocate(imap)
!!$
!!$     allocate(dmap(n), dvar(n))
!!$
!!$     ! We will fit the model map1 = A*map2 + B, with noise in both maps.
!!$     ! This requires a nonlinear search. Will try a grid first.
!!$     m = 10000
!!$     allocate(res(m,3))
!!$!     do i = 1, m
!!$!        a    = real(i-1,dp)*2/m
!!$!        dvar = var(:,1)+a**2*var(:,2)
!!$!        b    = sum((map(:,1)-a*map(:,2))/dvar)/sum(1/dvar)
!!$!        dmap = map(:,1) - a*map(:,2) - b
!!$!        L    = 0.5*sum(dmap**2/dvar) + 0.5*sum(log(2*pi*dvar))
!!$!write(*,'(5e15.7)') a, b, 0.5*sum(dmap**2/dvar), 0.5*log(sum(2*pi*dvar)), 0.5*sum(log(2*pi*dvar))
!!$!        res(i,:) = [ a, b, L ]
!!$!     end do
!!$m = 400
!!$allocate(res2(m,m,3))
!!$     do i = 1, m
!!$        do j = 1, m
!!$           a    = real(i+5,dp)*2/m
!!$           b    = real(j-m/2,dp)*40/m
!!$           dvar = var(:,1)+a**2*var(:,2)
!!$           dmap = map(:,1)-a*map(:,2) - b
!!$           L    = 0.5*sum(dmap**2/dvar)
!!$           res2(i,j,:) = [ a, b, L ]
!!$        end do
!!$     end do
!!$     res2(:,:,3) = res2(:,:,3) - minval(res2(:,:,3))
!!$     do i = 1, m
!!$        do j = 1, m
!!$           write(*,'(3e15.7)') res2(i,j,:)
!!$        end do
!!$        write(*,*)
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine wnoise
!!$     implicit none
!!$     character(len=512)             :: arg, fname
!!$     real(dp),        allocatable   :: map(:,:), omap(:,:)
!!$     integer(i4b),    allocatable   :: pixels(:)
!!$     integer(i4b)                   :: nside, npix, order, ncomp, i, j, k, n, m
!!$     real(dp)                       :: C(2,2), s
!!$     call getarg(2, fname)
!!$     call read_map(map, pixels, nside, order, fname)
!!$     allocate(omap(size(map,1),3))
!!$     s    = 6.594**2
!!$     omap = 0
!!$     do i = 1, size(map,1)
!!$        C(1,1) = map(i,2)
!!$        C(1,2) = map(i,3)
!!$        C(2,1) = map(i,3)
!!$        C(2,2) = map(i,4)
!!$        call invert_matrix(C)
!!$        omap(i,2) = s*C(1,1)*1e6
!!$        omap(i,3) = s*C(2,2)*1e6
!!$     end do
!!$     call getarg(3, fname)
!!$     call write_map(omap, pixels, nside, order, fname)
!!$   end subroutine
!!$
!!$   subroutine covtest4
!!$     implicit none
!!$     character(len=512)             :: arg
!!$     real(dp),        allocatable   :: cov(:,:)
!!$     integer(i4b),    allocatable   :: pixels(:)
!!$     integer(i4b)                   :: nside, npix, order, ncomp, i, j, k, n
!!$     nside = 2
!!$     npix  = 12*nside**2
!!$     order = nest
!!$     ncomp = 2
!!$     n     = ncomp*npix
!!$     allocate(pixels(npix))
!!$     pixels = irange(npix)-1
!!$     allocate(cov(n,n))
!!$     cov = 0
!!$     cov = cov+0.5
!!$     do i = 1, n
!!$        cov(i,i) = 1
!!$     end do
!!$     call write_covmat(cov, pixels, nside, order, ncomp,  "foo.unf")
!!$   end subroutine
!!$
!!$   subroutine covdump
!!$     implicit none
!!$     character(len=512)             :: arg
!!$     real(dp),        allocatable   :: cov(:,:)
!!$     integer(i4b),    allocatable   :: pixels(:)
!!$     integer(i4b)                   :: nside, npix, order, ncomp, i, j, k
!!$     call getarg(2, arg)
!!$     call read_covmat(cov, pixels, nside, order, ncomp, arg)
!!$     call dump_matrix(cov, fmt="(e15.7)")
!!$     do i = 1, size(cov, 1)
!!$        write(*,'(e15.7)') cov(i,i)
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine fit_gauss
!!$     implicit none
!!$     character(len=512)             :: arg
!!$     real(dp),        allocatable   :: map(:,:), data(:,:), omap(:,:)
!!$     integer(i4b),    allocatable   :: pixels(:), opix(:)
!!$     integer(i4b)                   :: i, j, k, n, nside, order, m, unit
!!$     real(dp)                       :: pos(2), mp(2), r, rmax, w, cov(2,2), vars(2)
!!$     real(dp)                       :: icov(2,2), p(7), v
!!$     call getarg(2, arg)
!!$     call read_map(map, pixels, nside, order, arg)
!!$     n = size(pixels)
!!$     rmax = 1000
!!$     allocate(data(n,3),opix(n))
!!$     m = 0
!!$     do i = 1, n
!!$        call pix2ang(nside, order, pixels(i), pos(2), pos(1))
!!$        if(pos(1) > pi) pos(1) = pos(1)-2*pi
!!$        pos(2) = pos(2)-pi/2
!!$        pos = pos*180*60/pi
!!$        r = sum(pos**2)**0.5
!!$        if(r > rmax) cycle
!!$        m = m+1
!!$        data(m,:) = [ pos(1), pos(2), map(i,1) ]
!!$        opix(m)   = pixels(i)
!!$     end do
!!$
!!$     p = [ maxval(data(:m,3)), data(maxloc(data(:m,3)),1:2), 1d0/6, 1d0/6, 0d0, 0d0 ]
!!$     call fit_helper(data(:m,:), p)
!!$
!!$     icov = reshape(p([4,6,6,5]),[2,2])
!!$     cov = icov; call invert_matrix(cov)
!!$     call get_eigenvalues(cov, vars)
!!$     write(*,'(4f10.5)') p(2:3), sqrt(vars)*sqrt(8*log(2.0))
!!$
!!$     open(unit,file="bar.txt")
!!$     do i = 1, m
!!$        write(unit,'(4f15.5)') data(i,:), evalgauss(data(i,:),p)
!!$     end do
!!$     close(unit)
!!$
!!$     allocate(omap(m,3))
!!$     do i = 1, m
!!$        pos = data(i,1:2)
!!$        v   = evalgauss(pos,p)
!!$        omap(i,:) = [ data(i,3), v, data(i,3)-v ]
!!$     end do
!!$
!!$     call write_map(omap, opix(:m), nside, order, "fit.hdf")
!!$   end subroutine
!!$
!!$   subroutine fit_helper(data, p)
!!$     implicit none
!!$     real(dp)           :: data(:,:), p(:)
!!$     integer(i4b)       :: err
!!$     if(allocated(powell_data)) deallocate(powell_data)
!!$     allocate(powell_data(size(data,1),size(data,2)))
!!$     powell_data = data
!!$     call powell(p, gausslik, err)
!!$   end subroutine
!!$
!!$   function evalgauss(pos, p) result(res)
!!$     implicit none
!!$     real(dp)       :: pos(:), p(:), res, C(2,2)
!!$     C = reshape(p([4,6,6,5]),[2,2])
!!$     res = p(1)*exp(-0.5*dot_product(pos-p(2:3),matmul(C,pos-p(2:3))))+p(7)
!!$   end function
!!$
!!$   function gausslik(p) result(res)
!!$     use healpix_types
!!$     implicit none
!!$     real(dp), dimension(:), intent(in), optional :: p
!!$     real(dp)     :: res
!!$     integer(i4b) :: i
!!$     res = 0
!!$     do i = 1, size(powell_data,1)
!!$        res = res + (powell_data(i,3)-evalgauss(powell_data(i,1:2),p))**2
!!$     end do
!!$     !write(*,'(7e10.2)') p, res
!!$   end function
!!$
!!$   subroutine suncomb
!!$     implicit none
!!$     character(len=512)             :: arg
!!$     integer(i4b)                   :: nside, order, npix, n, i, j, k
!!$     real(dp),          allocatable :: imap(:,:), irms(:,:), tsum(:,:), tw(:,:)
!!$     real(dp),          allocatable :: isig(:,:), iw(:,:)
!!$     n = (iargc()-2)/2
!!$     do i = 1, n
!!$        write(*,*) i
!!$        call getarg(i+1, arg)
!!$        call read_map(imap, order, arg)
!!$        call getarg(i+1+n, arg)
!!$        call read_map(irms, order, arg)
!!$        if(.not. allocated(tsum)) then
!!$           npix = size(imap,1)
!!$           allocate(tsum(0:npix-1,3), tw(0:npix-1,3), isig(0:npix-1,3), iw(0:npix-1,3))
!!$           tsum = 0; tw   = 0
!!$           isig = 0; iw   = 0
!!$        end if
!!$        where(healok(irms(:,2)) .and. healok(irms(:,3)))
!!$           isig(:,2) = imap(:,2)
!!$           isig(:,3) = imap(:,3)
!!$           iw(:,2)   = 1/irms(:,2)**2
!!$           iw(:,3)   = 1/Irms(:,3)**2
!!$        elsewhere
!!$           isig(:,2) = 0
!!$           isig(:,3) = 0
!!$           iw(:,2)   = 0
!!$           iw(:,3)   = 0
!!$        end where
!!$        tsum = tsum + isig*iw
!!$        tw   = tw   + iw
!!$        deallocate(imap, irms)
!!$     end do
!!$     where(tw > 0) tsum = tsum/tw
!!$     tsum(:,1) = sqrt(sum(tsum(:,2:3)**2,2))
!!$     call getarg(iargc(), arg)
!!$     call write_map(tsum, order, arg)
!!$   end subroutine
!!$
!!$   subroutine raul
!!$     implicit none
!!$     character(len=512)        :: parfile, ifile
!!$     integer(i4b)              :: unit, i, j, k, ces, mod, rdi, di, ok, ndi
!!$     real(dp)                  :: gain, sigma0, w, wsum, mjd
!!$     real(dp), allocatable     :: wtot(:,:)
!!$     call getarg(2, parfile)
!!$     call getarg(3, ifile)
!!$     call init_detector_mod(parfile)
!!$     call initialize_ces_mod(parfile)
!!$     call initialize_gain_mod(parfile)
!!$     ndi = size(quiet_diodes)
!!$     allocate(wtot(ndi,2))
!!$     wtot = 0
!!$     unit = getlun()
!!$     open(unit,file=ifile,action="read",status="old")
!!$     do
!!$     read(unit,*,end=1) ces, mod, rdi, sigma0
!!$        mjd  = ces_db%ceses(ces)%mjd(1)
!!$        di   = 4*mod+rdi+1
!!$        gain = get_gain(mjd, di)
!!$        w    = gain/sigma0**2
!!$        wtot(di,:) = wtot(di,:) + [w,1d0]
!!$     end do
!!$     1 close(unit)
!!$     wsum = sum(wtot)
!!$     do di = 1, ndi
!!$        if(wtot(di,2) == 0) cycle
!!$        write(*,'(2i3,f8.4,6i)') (di-1)/4, modulo(di-1,4), wtot(di,1)/wsum, nint(wtot(di,2))
!!$     end do
!!$   end subroutine
!!$
!!$   ! Given a set of mjds, look up the raster scans corresponding to
!!$   ! them and output their sigma0 level
!!$   subroutine rstsigma
!!$     implicit none
!!$     character(len=512)        :: parfile, ifile, ofile, l2file
!!$     integer(i4b)              :: i, j, unit, di, n
!!$     real(dp)                  :: mjds(2), mjd, pmjd, sigma0
!!$     type(Lx_struct)           :: data
!!$     call getarg(2, parfile)
!!$     call getarg(3, ifile)
!!$     call initialize_ces_mod(parfile)
!!$     pmjd = 0
!!$     unit = getlun()
!!$     open(unit,file=ifile,action="read",status="old")
!!$     do
!!$     read(unit,*,end=1) di, mjd
!!$        di = di+1
!!$        sigma0 = 0
!!$        n      = 0
!!$        do j = 1, size(ces_db%ceses)
!!$           if(abs(mean(ces_db%ceses(j)%mjd)-mjd)>0.01) cycle
!!$           l2file = ces_db%ceses(j)%l2file
!!$           if(mjd /= pmjd) then
!!$              call read_l2_file(l2file, data)
!!$              pmjd = mjd
!!$           end if
!!$           sigma0 = sigma0 + get_sigma0(real(data%tod(:,di),dp))
!!$           n = n+1
!!$        end do
!!$        sigma0 = sigma0/n
!!$        write(*,'(i4,f15.7,e15.7,i4)') di-1, mjd, sigma0, n
!!$     end do
!!$     1 close(unit)
!!$   end subroutine
!!$
!!$   function get_sigma0(tod) result(sigma0)
!!$     implicit none
!!$     real(dp)            :: tod(:), sigma0
!!$     integer(i4b)        :: n, m
!!$     complex(dpc), allocatable :: ft(:)
!!$     real(dp),     allocatable :: ps(:)
!!$     n = size(tod)
!!$     m = n/2+1
!!$     allocate(ft(m),ps(m))
!!$     call fft(tod, ft, 1)
!!$     call extract_powspec(ft, ps)
!!$     sigma0 = sqrt(mean(ps(m/2:)))
!!$     deallocate(ft,ps)
!!$   end function
!!$
!!$   subroutine scorr
!!$     implicit none
!!$     character(len=512)        :: ifile, ofile, cfile
!!$     integer(i4b)              :: unit, cid, i, b, j, d1, d2, nmax, nper, nb, dimax
!!$     integer(i4b)              :: npde
!!$     real(dp)                  :: sigma, siglim, msigma
!!$     type(hdf_file)            :: hfile
!!$     real(dp),     allocatable :: sigmas(:,:,:), bins(:,:), cov(:,:,:), corr(:,:)
!!$     real(dp),     allocatable :: means(:,:), pde(:,:,:)
!!$     nmax   = 10000
!!$     dimax  = 1000
!!$     nper   = 60
!!$     npde   = 1000
!!$     siglim = 1.5
!!$     msigma = 1e-4
!!$     call getarg(2, ifile)
!!$     call getarg(3, ofile)
!!$     allocate(sigmas(dimax,nmax,2))
!!$     sigmas = 0
!!$     nmax = 0; dimax = 0
!!$     unit = getlun()
!!$     open(unit,file=ifile,status="old",action="read")
!!$     do
!!$     read(unit,*,end=1) cid, d1, sigma
!!$        d1 = d1+1
!!$        sigmas(d1,cid,:) = [sigma,1d0]
!!$        nmax  = max(nmax,cid)
!!$        dimax = max(dimax,d1)
!!$     end do
!!$     1 close(unit)
!!$     ! Filter out outliers by making bins in cid
!!$     nb = (nmax+nper-1)/nper
!!$     allocate(bins(nb,3), means(dimax,2))
!!$     allocate(pde(npde,nb,dimax))
!!$     pde   = 0
!!$     means = 0
!!$     do d1 = 1, dimax
!!$        bins = 0
!!$        do i = 1, nmax
!!$           if(sigmas(d1,i,2)==0) cycle
!!$           b = (i-1)/nper+1
!!$           j = min(floor(sigmas(d1,i,1)*npde/msigma)+1,npde)
!!$           pde(j,b,d1) = pde(j,b,d1) + 1
!!$           bins(b,1) = bins(b,1) + sigmas(d1,i,1)
!!$           bins(b,2) = bins(b,2) + sigmas(d1,i,1)**2
!!$           bins(b,3) = bins(b,3) + 1
!!$        end do
!!$        do b = 1, nb
!!$           pde(:,b,d1) = pde(:,b,d1)/sum(pde(:,b,d1))
!!$        end do
!!$        where(bins(:,3)>0)
!!$           bins(:,1) = bins(:,1)/bins(:,3)
!!$           bins(:,2) = sqrt(bins(:,2)/bins(:,3)-bins(:,1)**2)
!!$        end where
!!$        ! Remove all outliers
!!$        do i = 1, nmax
!!$           if(sigmas(d1,i,2)==0) cycle
!!$           b = (i-1)/nper+1
!!$           if(abs(sigmas(d1,i,1)-bins(b,1))/bins(b,2) > siglim) &
!!$             sigmas(d1,i,:) = 0
!!$        end do
!!$     end do
!!$     ! Output the cleaned ones
!!$     if(iargc() > 3) then
!!$        call getarg(4, cfile)
!!$        open(unit,file=cfile)
!!$        do i = 1, nmax
!!$           do d1 = 1, dimax
!!$              if(sigmas(d1,i,2) == 0) cycle
!!$              write(unit,'(i5,i4,e15.7)') i, d1-1, sigmas(d1,i,1)
!!$           end do
!!$        end do
!!$        close(unit)
!!$     end if
!!$
!!$     do d1 = 1, dimax
!!$        means(d1,:) = sum(sigmas(d1,:,:),1)
!!$     end do
!!$     means(:,1) = means(:,1)/means(:,2)
!!$
!!$     ! Find the covariance for what remains
!!$     allocate(cov(dimax,dimax,2), corr(dimax,dimax))
!!$     cov = 0
!!$     do i = 1, nmax
!!$        do d1 = 1, dimax
!!$           if(sigmas(d1,i,2) == 0) cycle
!!$           do d2 = 1, dimax
!!$              if(sigmas(d2,i,2) == 0) cycle
!!$              cov(d1,d2,1) = cov(d1,d2,1) + (sigmas(d1,i,1)-means(d1,1))*(sigmas(d2,i,1)-means(d2,1))
!!$              cov(d1,d2,2) = cov(d1,d2,2) + 1
!!$           end do
!!$        end do
!!$     end do
!!$     cov(:,:,1) = cov(:,:,1)/cov(:,:,2)
!!$     do d1 = 1, dimax
!!$        do d2 = 1, dimax
!!$           corr(d1,d2) = cov(d1,d2,1)/sqrt(cov(d1,d1,1)*cov(d2,d2,1))
!!$        end do
!!$     end do
!!$     call open_hdf_file(ofile, hfile, "w")
!!$     call write_hdf(hfile, "cov", cov(:,:,1))
!!$     call write_hdf(hfile, "corr", corr)
!!$     call write_hdf(hfile, "pde",  pde)
!!$     call close_hdf_file(hfile)
!!$
!!$
!!$   end subroutine
!!$
!!$   subroutine bleh
!!$     implicit none
!!$     character(len=512)    :: ifile, ofile
!!$     real(dp), allocatable :: map(:,:,:)
!!$     integer(i4b), allocatable :: pixels(:)
!!$     integer(i4b)          :: nside, order, i, j
!!$     call getarg(2, ifile)
!!$     call getarg(3, ofile)
!!$     call read_map(map, pixels, nside, order, ifile)
!!$     map(:,1,:) = sqrt(sum(map(:,2:3,:)**2,2))
!!$     call write_map(map, pixels, nside, order, ofile)
!!$   end subroutine
!!$
!!$   subroutine radconv
!!$     implicit none
!!$     character(len=512)    :: ifile
!!$     real(dp), allocatable :: data(:,:), smooth(:)
!!$     real(dp)              :: tmp(3)
!!$     integer(i4b)          :: unit, n, i, j
!!$     call getarg(2, ifile)
!!$     unit = getlun()
!!$     open(unit,file=ifile,status="old",action="read")
!!$     n = 0
!!$     do
!!$        read(unit,*,end=1) tmp
!!$        n = n+1
!!$     end do
!!$     1 rewind(unit)
!!$     allocate(data(n,3),smooth(n))
!!$     do i = 1, n
!!$        read(unit,*) data(i,:)
!!$     end do
!!$     close(unit)
!!$     smooth = 0
!!$     call convolute_radial_3d(data(2:,2),data(2:,3),smooth(2:))
!!$     do i = 1, n
!!$        write(*,'(4e15.7)',err=2) data(i,:), smooth(i)
!!$     end do
!!$     2 deallocate(data, smooth)
!!$   end subroutine
!!$
!!$   subroutine dst(a,b,dir)
!!$     implicit none
!!$     real(dp)     :: a(:),b(:)
!!$     integer(i4b) :: dir, type
!!$     integer*8    :: plan
!!$     type = fftw_rodft00; if(dir < 0) type = fftw_rodft00
!!$     call dfftw_plan_r2r_1d(plan, size(a), a, b, type, &
!!$          & fftw_estimate + fftw_unaligned)
!!$     call dfftw_execute_r2r(plan, a, b)
!!$     call dfftw_destroy_plan(plan)
!!$     b = b/sqrt(real(2*size(b)+2,dp))
!!$   end subroutine
!!$
!!$   ! Pass in arrays with the 0th element omitted.
!!$   ! This is also omitted in the output, and will
!!$   ! have to be extrapolated if you need it.
!!$   subroutine convolute_radial_3d(a,b,c)
!!$     implicit none
!!$     real(dp)                  :: a(:), b(:), c(:)
!!$     real(dp),     allocatable :: x(:), fa(:), fb(:)
!!$     integer(i4b)              :: r,k,n,m
!!$     n = size(a)
!!$     allocate(fa(n),fb(n),x(n))
!!$     do r = 1, n; x(r) = a(r)*r; end do
!!$     call dst(x,fa,1)
!!$     do r = 1, n; x(r) = b(r)*r; end do
!!$     call dst(x,fb,1)
!!$     do k = 1, n; x(k) = fa(k)*fb(k)/k; end do
!!$     call dst(x,c,-1)
!!$     c = c * (2*pi)**1.5
!!$     do r = 1, n; c(r) = c(r)/r; end do
!!$   end subroutine
!!$
!!$   subroutine gainres
!!$     implicit none
!!$     character(len=512)        :: parfile, dfile
!!$     integer(i4b)              :: unit, di
!!$     real(dp)                  :: t, obs, dobs, gain
!!$     call getarg(2, parfile)
!!$     call getarg(3, dfile)
!!$     call initialize_gain_mod(parfile)
!!$     unit = getlun()
!!$     open(unit, file=dfile, status="old", action="read")
!!$     do
!!$        read(unit,*,end=1) di, t, obs, dobs
!!$        gain = get_gain(t, di+1)
!!$        write(*,'(i5,f15.7,4e15.7)') di, t, obs, dobs, gain, (obs-gain)/dobs
!!$     end do
!!$     1 close(unit)
!!$   end subroutine
!!$
!!$   subroutine getgains2
!!$     implicit none
!!$     character(len=512)        :: parfile, arg
!!$     integer(i4b)              :: unit, di
!!$     real(dp)                  :: t, gain, i, d, n, tr(2)
!!$     call getarg(2, parfile)
!!$     call getarg(3, arg); read(arg,*) n
!!$     tr = [55050d0, 55550d0]
!!$     call initialize_gain_mod(parfile)
!!$     do i = 1, n
!!$        t = tr(1)+(tr(2)-tr(1))*(i-1)/(n-1)
!!$        write(*,'(f15.7)',advance="no") t
!!$        do di = 1, size(quiet_diodes)
!!$           gain = get_gain(t,di)
!!$           write(*,'(e15.7)',advance="no") gain
!!$        end do
!!$        write(*,*)
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine getgains
!!$     implicit none
!!$     character(len=512)        :: parfile, tfile
!!$     integer(i4b)              :: unit, di
!!$     real(dp)                  :: t, gain
!!$     call getarg(2, parfile)
!!$     call getarg(3, tfile)
!!$     call initialize_gain_mod(parfile)
!!$     unit = getlun()
!!$     open(unit, file=tfile, status="old", action="read")
!!$     do
!!$        read(unit,*,end=1) di, t
!!$        write(*,'(e15.7)') get_gain(t, di+1)
!!$     end do
!!$     1 close(unit)
!!$   end subroutine
!!$
!!$   subroutine eigtime
!!$     implicit none
!!$     character(len=512)        :: arg
!!$     integer(i4b)              :: i, j, k, m, n, maxexp, status
!!$     real(dp)                  :: t1, t2
!!$     real(dp), allocatable     :: A(:,:), V(:,:), E(:)
!!$     type(planck_rng)          :: rng
!!$     call getarg(2, arg); read(arg,*) maxexp
!!$     call rand_init(rng, 1)
!!$     do i = 1, maxexp
!!$        n = 2**i
!!$        allocate(A(n,n), V(n,n), E(n))
!!$        do j = 1, n
!!$           do k = j, n
!!$              A(k,j) = rand_gauss(rng)
!!$              A(j,k) = A(k,j)
!!$           end do
!!$           A(j,j) = A(j,j) + 10*n
!!$        end do
!!$        m = 1000000/n/n/n+1
!!$        call wall_time(t1)
!!$        do j = 1, m
!!$           call eigen_decomp(A, E, V, status)
!!$        end do
!!$        call wall_time(t2)
!!$        write(*,'(i8,e15.7,i5)') n, (t2-t1)/m, status
!!$        deallocate(A,V,E)
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine deconvolution
!!$     implicit none
!!$     character(len=512)        :: fname, arg
!!$     real(dp),     allocatable :: map(:,:), ps(:,:)
!!$     complex(dpc), allocatable :: alm(:,:,:)
!!$     integer(i4b)              :: order, nside, lmax, i, j, k, l, m, unit
!!$     real(dp)                  :: b1, b2, sigma
!!$     type(planck_rng)          :: rng
!!$     call getarg(2, fname)
!!$     call getarg(3, arg);  read(arg,*) b1; b1 = b1*pi/180/sqrt(8*log(2d0))
!!$     call getarg(4, arg);  read(arg,*) b2; b2 = b2*pi/180/sqrt(8*log(2d0))
!!$     call getarg(5, arg);  read(arg,*) sigma
!!$     call rand_init(rng,1)
!!$     call read_map(map, order, fname)
!!$     nside = npix2nside(size(map,1))
!!$     lmax = 2*nside
!!$     allocate(alm(1:1, 0:lmax, 0:lmax))
!!$     call map2alm(nside, lmax, lmax, map(:,1), alm, [-1d0,1d0], get_hpix_ringweights(nside))
!!$     allocate(ps(0:lmax,1))
!!$     call alm2cl(lmax, lmax, alm, ps)
!!$     call dump_matrix(ps, "ps_orig.txt", fmt="(e15.7)")
!!$     do l = 0, lmax
!!$        alm(:,l,:) = alm(:,l,:) * exp(-0.5*l*(l+1)*b1**2)
!!$        do m = 0, lmax
!!$           alm(:,l,m) = alm(:,l,m) + rand_gauss(rng)*sigma*exp(-(l/470d0)**4)
!!$        end do
!!$     end do
!!$     call alm2cl(lmax, lmax, alm, ps)
!!$     call dump_matrix(ps, "ps_b1.txt", fmt="(e15.7)")
!!$     call alm2map(nside, lmax, lmax, alm, map(:,1))
!!$     call write_map(map, order, "map_b1.hdf")
!!$     do l = 0, lmax
!!$        alm(:,l,:) = alm(:,l,:) * exp(+0.5*l*(l+1)*b2**2)
!!$     end do
!!$     call alm2cl(lmax, lmax, alm, ps)
!!$     call dump_matrix(ps, "ps_b2.txt", fmt="(e15.7)")
!!$     call alm2map(nside, lmax, lmax, alm, map(:,1))
!!$     call write_map(map, order, "map_b2.hdf")
!!$   end subroutine
!!$
!!$   subroutine foo
!!$     implicit none
!!$     complex(dpc), allocatable :: alm(:,:,:)
!!$     real(dp),     allocatable :: map(:,:)
!!$     integer(i4b)              :: order, nside, lmax
!!$     character(len=512)        :: fname
!!$     call dmem("A")
!!$     call getarg(2, fname)
!!$     call dmem("B")
!!$     call read_map(map, order, fname)
!!$     call dmem("C")
!!$     nside = npix2nside(size(map,1))
!!$     lmax  = 2*nside
!!$     call dmem("D")
!!$     allocate(alm(1:3, 0:lmax, 0:lmax))
!!$     call map2alm(nside, lmax, lmax, map, alm, [-1d0,1d0], get_hpix_ringweights(nside))
!!$     call dmem("E")
!!$   end subroutine
!!$
!!$   subroutine noise_est
!!$     implicit none
!!$     character(len=512)    :: parfile, arg, ofile
!!$     real(dp)              :: ipar(3), opar(3), srate, chi, pstep
!!$     type(planck_rng)      :: rng
!!$     integer(i4b)          :: n, myid, nproc, err, step, i, j, r(2), np
!!$     type(hdf_file)        :: file
!!$     real(dp), allocatable :: ps(:),ps0(:),ps1(:), mychisq(:), chisq(:), px(:)
!!$     integer(i4b), allocatable :: pbin(:,:), mypbin(:,:)
!!$     n     = 100000
!!$     ipar  = [ 1d0, 0.2d0, -2d0 ]
!!$     srate = 25
!!$     pstep = 0.001
!!$     np    = 10000
!!$     call getarg(2, parfile)
!!$     call getarg(3, arg); read(arg,*) step
!!$     call getarg(4, ofile)
!!$     call initialize_noise_estimation_mod(parfile)
!!$
!!$     call mpi_init(err)
!!$     call mpi_comm_rank(MPI_COMM_WORLD, myid, err)
!!$     call mpi_comm_size(MPI_COMM_WORLD, nproc,err)
!!$     call rand_init(rng, myid+7718)
!!$
!!$     allocate(ps(n),ps0(n),ps1(n),mychisq(n),chisq(n),pbin(np,3),mypbin(np,3),px(np))
!!$     px = (irange(np)-1)*pstep
!!$     mychisq = 0
!!$     mypbin  = 0
!!$     i = 0
!!$     do
!!$        i=i+1
!!$        call makespec(ipar(1), ipar(2), ipar(3), srate, ps0)
!!$        do
!!$           call sim_noise(rng, ps0, ps)
!!$           call fit_1overf_profile(srate, 0.1d0, 0.02d0, opar(1), opar(3), opar(2), tod_ps=ps, apply_scanmask=.true.)
!!$           if(opar(3) > -8) exit
!!$        end do
!!$        call makespec(opar(1), opar(2), opar(3), srate, ps1)
!!$        mychisq = mychisq + ps/ps1*2
!!$        j=floor( opar(1)/pstep)+1; mypbin(j,1) = mypbin(j,1)+1
!!$        j=floor( opar(2)/pstep)+1; mypbin(j,2) = mypbin(j,2)+1
!!$        j=floor(-opar(3)/pstep)+1; mypbin(j,3) = mypbin(j,3)+1
!!$        if(modulo(i,step) == 0) then
!!$           call mpi_allreduce(mychisq, chisq, size(chisq), mpi_double_precision, MPI_SUM, mpi_comm_world, err)
!!$           call mpi_allreduce(mypbin, pbin, size(pbin), mpi_integer, MPI_SUM, mpi_comm_world, err)
!!$           if(myid == 0) then
!!$              write(*,*) i
!!$              call open_hdf_file(ofile, file, "w")
!!$              call write_hdf(file, "chisq", chisq(2:))
!!$              call write_hdf(file, "count", 2*i*nproc)
!!$              call write_hdf(file, "sigma", (chisq(2:)-2*i*nproc)/sqrt(2d0*2*i*nproc))
!!$              call write_hdf(file, "pbin",  pbin)
!!$              call write_hdf(file, "xbin",  px)
!!$              call close_hdf_file(file)
!!$           end if
!!$        end if
!!$     end do
!!$     call mpi_finalize(err)
!!$   end subroutine
!!$
!!$!   subroutine noise_est2
!!$!     implicit none
!!$!     character(len=512)    :: parfile, arg, ofile
!!$!     real(dp)              :: ipar(3), opar(3), srate, chi, pstep
!!$!     type(planck_rng)      :: rng
!!$!     integer(i4b)          :: n, myid, nproc, err, step, i, j, r(2), np
!!$!     type(hdf_file)        :: file
!!$!     real(dp), allocatable :: ps(:),ps0(:),ps1(:), y(:,:)
!!$!     n     = 100000
!!$!     ipar  = [ 1d0, 0.1d0, -2d0 ]
!!$!     srate = 25
!!$!     call getarg(2, parfile)
!!$!     call getarg(3, arg); read(arg,*) step
!!$!     call getarg(4, ofile)
!!$!     call initialize_noise_estimation_mod(parfile)
!!$!
!!$!     allocate(ps(n),ps0(n),ps1(n),y(step,2))
!!$!     call makespec(ipar(1), ipar(2), ipar(3), srate, ps0)
!!$!     call debug_fit(srate, 0.1d0, 0.02d0, ipar(1), ipar(2), ipar(3), ps0, y)
!!$!     call open_hdf_file(ofile, file, "w")
!!$!     call write_hdf(file, "lik", y)
!!$!     call close_hdf_file(file)
!!$!   end subroutine
!!$
!!$
!!$
!!$   subroutine makespec(sigma0, fknee, alpha, srate, spec)
!!$     implicit none
!!$     real(dp)     :: sigma0, fknee, alpha, spec(:), x, srate
!!$     integer(i4b) :: i
!!$     spec(1) = 0
!!$     do i = 2, size(spec)
!!$       x = ind2freq(i, srate, size(spec))
!!$       spec(i) = sigma0**2*(1+(x/fknee)**alpha)
!!$     end do
!!$   end subroutine
!!$  subroutine sim_noise(rng, ps0, ps)
!!$    implicit none
!!$    real(dp)                  :: sigma0, fknee, alpha, ps0(:), ps(:)
!!$    type(planck_rng)          :: rng
!!$    integer(i4b)              :: i, m
!!$    complex(dpc), allocatable :: ft(:)
!!$    m = size(ps)
!!$    allocate(ft(m))
!!$    ft(1) = 0
!!$    do i = 2, m-1
!!$       ft(i) = cmplx(rand_gauss(rng),rand_gauss(rng))/2**0.5
!!$    end do
!!$    ft(m) = rand_gauss(rng)
!!$    do i = 2, m
!!$       ft(i) = ft(i) * ps0(i)**0.5
!!$    end do
!!$    call extract_powspec(ft, ps)
!!$    deallocate(ft)
!!$  end subroutine
!!$   !subroutine noise_est
!!$   !  implicit none
!!$   !  real(dp)              :: par1(3), par(3), ext(2), x, sigma1, da, df
!!$   !  real(dp)              :: s, sa, sf, dS2f, dS2a, s0, s0f, s0a
!!$   !  integer(i4b)          :: n, i, m
!!$   !  real(dp), allocatable :: a(:), af(:), aa(:)
!!$   !  par = [ 1d0, 0.1d0, -1.5d0 ]
!!$   !  n   = 1000000
!!$   !  ext = [ 0.2d0, 0.9d0 ]
!!$   !  da  = 0.00001; df = 0.00001
!!$   !  allocate(a(n),af(n),aa(n))
!!$   !  s0  = par(1)
!!$   !  s   = s0/fix_sigma(1d0,par(2),par(3),ext)
!!$   !  s0f = fix_sigma(s,par(2)+df,par(3),ext)
!!$   !  s0a = fix_sigma(s,par(2),par(3)+da,ext)
!!$   !  call makespec(s0, par(2),   par(3),    a)
!!$   !  call makespec(s0f,par(2)+df,par(3),    af)
!!$   !  call makespec(s0a,par(2),par(3)+da,    aa)
!!$   !  call deriv_sigma2(s, par(2), par(3), ext, dS2f, dS2a)
!!$   !  write(*,*) (s0f**2-s0**2)/df, (s0a**2-s0**2)/da
!!$   !  write(*,*) dS2f, dS2a
!!$   !end subroutine
!!$   !function getsig(spec, ext) result(sig)
!!$   !  implicit none
!!$   !  real(dp) :: spec(:), ext(:), sig, x
!!$   !  integer(i4b) :: i, m
!!$   !  m = 0
!!$   !  sig = 0
!!$   !  do i = 1, size(spec)
!!$   !    x = i*1d0/size(spec)
!!$   !    if(x >= ext(1) .and. x <= ext(2)) then
!!$   !       sig = sig + spec(i)
!!$   !       m   = m+1
!!$   !    end if
!!$   !  end do
!!$   !  sig = sqrt(sig/m)
!!$   !end function
!!$   !function fix_sigma(sigma1, fknee, alpha, srange) result(sigma0)
!!$   !  implicit none
!!$   !  real(dp)     :: sigma1, sigma0, fknee,alpha, srange(2), a, b
!!$   !  a = srange(1)
!!$   !  b = srange(2)
!!$   !  sigma0 = sigma1/sqrt(1 + (b**(alpha+1)-a**(alpha+1))/(fknee**alpha*(alpha+1)*(b-a)))
!!$   !end function
!!$   !subroutine deriv_sigma2(sigma1, fknee, alpha, srange, dS2df, dS2da)
!!$   !  implicit none
!!$   !  real(dp)     :: sigma1, sigma0, fknee,alpha, srange(2), a, b, dS2df, dS2da
!!$   !  real(dp)     :: d, da, a1, fa
!!$   !  a = srange(1)
!!$   !  b = srange(2)
!!$   !  d = b-a; a1 = alpha+1; da = b**a1 - a**a1; fa = fknee**(-alpha)
!!$   !  sigma0 = sigma1/sqrt(1+da/(d*a1)*fa)
!!$   !  dS2df = alpha/a1 * sigma0**4/sigma1**2*da/d*fknee**(-alpha-1)
!!$   !  dS2da = -sigma0**4/sigma1**2/d*(1/a1*fa*(log(b/fknee)*b**a1-log(a/fknee)*a**a1)-&
!!$   !   & 1/a1**2*(b**a1-a**a1)/fknee**alpha)
!!$   !end subroutine
!!$
!!$   subroutine mpitest
!!$     implicit none
!!$     integer(i4b) :: myid, nproc, err, i, j, k, n, m, nwork, foo(10)
!!$     real(dp)     :: x, y
!!$     call mpi_init(err)
!!$     call mpi_comm_rank(mpi_comm_world, myid,  err)
!!$     call mpi_comm_size(mpi_comm_world, nproc, err)
!!$     m = 100
!!$     nwork = nproc-1
!!$     n = m*nwork
!!$     if(myid == 0) then
!!$        ! Jeg er sjefen! Hahahaha!!
!!$        do i = 1, m
!!$           do j = 1, nwork
!!$              x = real((i-1)*nwork+j,dp)/n
!!$              call mpi_send(x, 1, mpi_double_precision, j, 0, mpi_comm_world, err)
!!$           end do
!!$           do j = 1, nwork
!!$              x = real((i-1)*nwork+j,dp)/n
!!$              call mpi_recv(y, 1, mpi_double_precision, j, 0, mpi_comm_world, foo, err)
!!$              write(*,'(i4,2f9.4)') j, x, y
!!$           end do
!!$        end do
!!$        do j = 1, nwork
!!$           call mpi_send(nan, 1, mpi_double_precision, j, 0, mpi_comm_world, err)
!!$        end do
!!$     else
!!$        ! Å nei, jeg er en slave
!!$        do
!!$           call mpi_recv(x, 1, mpi_double_precision, 0, 0, mpi_comm_world, foo, err)
!!$           if(x /= x) exit
!!$           call fsleep(1d0)
!!$           y = exp(x)
!!$           call mpi_send(y, 1, mpi_double_precision, 0, 0, mpi_comm_world, err)
!!$        end do
!!$     end if
!!$     call mpi_finalize(err)
!!$   end subroutine
!!$
!!$   ! Why does mpi reduce cause crashes all of a sudden?
!!$   subroutine anothermpitest
!!$     implicit none
!!$     integer(i4b)             :: err, myid, nproc, array(100), total(100), i
!!$     call mpi_init(err) ! Starting up mpi
!!$     call mpi_comm_rank(mpi_comm_world, myid, err) ! assigning an id to each process
!!$     call mpi_comm_size(mpi_comm_world, nproc, err) ! counting the number of processes?
!!$
!!$     do i=1,100
!!$        array(i) = i
!!$     end do
!!$     call mpi_reduce(array, total, size(array), mpi_integer, mpi_sum, 0, mpi_comm_world, err)
!!$
!!$     if (myid==0) then
!!$        do i=1,100
!!$           write(*,*) total(i)
!!$        end do
!!$     end if
!!$     call mpi_finalize(err)
!!$
!!$   end subroutine anothermpitest
!!$
!!$
!!$
!!$
!!$   ! Produce akito alpha, beta from theta, phi
!!$   function akito(cin) result(cout)
!!$     real(dp) :: cin(2), cout(2)
!!$     cout(1) = -cin(1)*cos(cin(2)+pi/2)
!!$     cout(2) =  cin(1)*sin(cin(2)+pi/2)
!!$   end function
!!$
!!$   ! Produce theta, phi from akito alpha, beta
!!$   function akifrom(cin) result(cout)
!!$     real(dp) :: cin(2), cout(2)
!!$     cout(1) = sqrt(sum(cin**2))
!!$     cout(2) = atan2(cin(2),-cin(1))-pi/2
!!$     if(cout(2) < 0) cout(2) = cout(2) + 2*pi
!!$   end function
!!$
!!$   subroutine pointtest4
!!$     implicit none
!!$     character(len=512)     :: arg, params, afile, line
!!$     integer(i4b)           :: i, j, k, nmod, unit, mod
!!$     real(dp)               :: x, y, theta, phi, ab(2), coeff, a, b, c, d, az, el, dk
!!$     real(dp)               :: del, ab2(2), p(2), p2(2), dev, coaz
!!$     real(dp), allocatable  :: ahorns(:,:)
!!$     type(quiet_mount)      :: mount
!!$     call getarg(2, params)
!!$     call getarg(3, afile)
!!$     call init_detector_mod(params)
!!$     call initialize_pmac_mod(params, "")
!!$     call initialize_quiet_pointing_mod(params)
!!$     mount = mount_none
!!$     mount%akito_el  = [-0.4283d0,-1.374d-2,-1.928d-4,-5.151d-5, 54.0*DEG2RAD]
!!$     mount%akito_dk  = [ 0.1264d0,-141*DEG2RAD ]
!!$     mount%akito_dir = 28.2*DEG2RAD
!!$     call set_mount_override(.true., mount_model=mount)
!!$     nmod = size(quiet_horns)
!!$
!!$     el =  50*DEG2RAD
!!$     dk =  70*DEG2RAD
!!$     ! Apply correction to each horn
!!$     del   = (el - mount%akito_el(5))*RAD2DEG
!!$     coeff = sum(mount%akito_el(:4)*[1d0,del,del**2,del**3]) + &
!!$      & mount%akito_dk(1)*sin(dk-mount%akito_dk(2))
!!$     do mod = 0, nmod-1
!!$        p     = [quiet_horns(mod)%theta, quiet_horns(mod)%phi]
!!$        coaz  = p(1)*sin(p(2)+dk)
!!$        dev   = coaz*coeff*RAD2DEG
!!$
!!$        ab    = akito(p)
!!$        ab2   = ab - dev/60*DEG2RAD*[cos(mount%akito_dir),-sin(mount%akito_dir)]
!!$        p2    = akifrom(ab2)
!!$        write(*,'(i3,8f10.5)') mod, p*RAD2DEG, p2*RAD2DEG, ab*RAD2DEG, ab2*RAD2DEG
!!$     end do
!!$
!!$   end subroutine
!!$
!!$   ! Turn per-pixel rhs into effective inverse cov
!!$   subroutine r2c
!!$     implicit none
!!$     character(len=512)             :: ifile, ofile, arg, pfile
!!$     integer(i4b)                   :: ncomp, i, j, k, m, n, nside, order, o
!!$     type(hdf_file)                 :: hfile
!!$     real(dp),          allocatable :: imaps(:,:,:), oicov(:,:,:,:)
!!$     integer(i4b),      allocatable :: pixels(:), ppixels(:), map2mask(:)
!!$     call getarg(2, ifile)
!!$     call getarg(3, pfile)
!!$     call getarg(4, ofile)
!!$     call read_map(imaps, pixels, nside, order, ifile)
!!$     call open_hdf_file(pfile, hfile, "r")
!!$     call read_alloc_hdf(hfile, "pixels", ppixels)
!!$     call close_hdf_file(hfile)
!!$     ! We now have the pixel set for the input patch, corresponding to
!!$     ! the last index of the rhss. We want to remap this to the current pixel
!!$     ! set.
!!$     write(*,*) "A", ppixels
!!$     write(*,*) "B", pixels
!!$
!!$     allocate(map2mask(0:12*nside**2-1))
!!$     map2mask = 0
!!$     do i = 1, size(ppixels)
!!$        map2mask(ppixels(i)) = i
!!$     end do
!!$     n = size(pixels)
!!$     allocate(oicov(n,2,n,2))
!!$     oicov = 0
!!$     do i = 1, n
!!$        j = map2mask(pixels(i))
!!$        if(j <= 0) cycle
!!$        write(*,'(9i5)') i, j, shape(oicov), shape(imaps)
!!$        oicov(:,:,i,:) = imaps(:,2:3,[2*j-1,2*j])
!!$     end do
!!$     call write_covmat(oicov, pixels, nside, order, ofile)
!!$   end subroutine
!!$
!!$   subroutine pixmaps
!!$     implicit none
!!$     character(len=512)             :: ifile, ofile, arg
!!$     integer(i4b)                   :: ncomp, i, j, k, m, n, nside, order, o
!!$     real(dp),          allocatable :: imaps(:,:,:), omaps(:,:,:)
!!$     integer(i4b),      allocatable :: pixels(:)
!!$     call getarg(1, ifile)
!!$     call getarg(2, arg); read(arg,*) ncomp
!!$     call getarg(3, ofile)
!!$     call read_map(imaps, pixels, nside, order, ifile)
!!$     n = size(pixels)
!!$     allocate(omaps(n,3,n*ncomp))
!!$     o = 0; if(ncomp == 2) o = 1
!!$     omaps = 0
!!$     k = 0
!!$     do i = 1, n
!!$        do j = 1, ncomp
!!$           k = k+1
!!$           omaps(i,j+o,k) = 1
!!$        end do
!!$     end do
!!$     call write_map(omaps, pixels, nside, order, ofile)
!!$   end subroutine
!!$
!!$   ! Mini-scalapost: Read equation set, solve it and output covariance
!!$   ! matrix.
!!$   subroutine covtest3
!!$     implicit none
!!$     character(len=512) :: eqnfile, prefix
!!$     integer(i4b)       :: unit, i, j, k, m, n, order, ncomp, nmap, npix, nside, nsim, o
!!$     type(hdf_file)     :: hfile
!!$     real(dp),      allocatable :: icov1(:,:), icov2(:,:), omaps(:,:,:), joint_icov(:,:)
!!$     real(dp),      allocatable :: rhss(:,:), rhss2(:,:), maps(:,:), cov1(:,:), cov2(:,:)
!!$     real(dp),      allocatable :: eigs(:,:), cov(:,:)
!!$     integer(i4b),  allocatable :: map2mask(:), inds(:), pixels(:), tpix(:)
!!$     call getarg(2, eqnfile)
!!$     call getarg(3, prefix)
!!$     unit = getlun()
!!$     open(unit,file=eqnfile,status="old",action="read",form="unformatted")
!!$     read(unit) n
!!$     read(unit) order
!!$     read(unit) ncomp
!!$     allocate(joint_icov(n+1,n))
!!$     do i = 1, n
!!$        read(unit) joint_icov(:,i)
!!$     end do
!!$     read(unit) nmap
!!$     allocate(rhss(n,nmap))
!!$     do i = 1, nmap
!!$        read(unit) rhss(:,i)
!!$     end do
!!$     read(unit) npix
!!$     allocate(map2mask(0:npix-1))
!!$     read(unit) map2mask
!!$     nside = npix2nside(npix)
!!$     close(unit)
!!$
!!$     call open_hdf_file("foo.hdf", hfile, "w")
!!$     call write_hdf(hfile, "mat", joint_icov)
!!$     call close_hdf_file(hfile)
!!$
!!$     call dump_matrix(real(reshape(get_diag(joint_icov(:n,:n)),[n,1]),dp),"foo.txt")
!!$
!!$     ! Get the interesting indices
!!$     call wherei(get_diag(joint_icov(:n,:n)) /= 0, inds)
!!$     m = size(inds)
!!$     call dump_matrix(reshape(real(inds,dp),[m,1]), "inds.txt")
!!$
!!$     write(*,*) "n:", n, "m:", m
!!$
!!$     allocate(pixels(m/ncomp), tpix(n), icov1(m,m), icov2(m,m), rhss2(m,nmap-1))
!!$     icov1 = 0
!!$     icov2 = 0
!!$     do i = 1, m
!!$        do j = 1, i
!!$           icov1(i,j) = joint_icov(inds(j),inds(i))
!!$           icov1(j,i) = icov1(i,j)
!!$           icov2(i,j) = joint_icov(inds(i)+1,inds(j))
!!$           icov2(j,i) = icov2(i,j)
!!$        end do
!!$     end do
!!$
!!$     rhss2 = rhss(inds,2:)
!!$
!!$     j = 0
!!$     do i = 0, npix-1
!!$        if(map2mask(i) <= 0) cycle
!!$        j = j+1
!!$        tpix(j) = i
!!$     end do
!!$     pixels = tpix(inds)
!!$     deallocate(joint_icov, rhss, map2mask, tpix)
!!$
!!$     ! We can now solve for the map. We also want the matrices themselves.
!!$     nsim = nmap-1
!!$     allocate(maps(m,nsim), cov1(m,m), cov2(m,m), omaps(m/ncomp,3,nsim))
!!$     call eigen_pow(icov1, -1d0, cov1)
!!$     call eigen_pow(icov2, -1d0, cov2)
!!$
!!$     ! True cov
!!$     allocate(cov(m,m))
!!$     cov = matmul(cov1,matmul(icov2,cov1))
!!$
!!$     allocate(eigs(m,3))
!!$     call get_eigenvalues(icov1, eigs(:,1))
!!$     call get_eigenvalues(icov2, eigs(:,2))
!!$     call get_eigenvalues(cov,   eigs(:,3))
!!$     call dump_matrix(eigs, "ieigs.txt")
!!$
!!$     call dgemm('N','N',m,nsim,m,1d0,cov1,m,rhss2,m,0d0,maps,m)
!!$
!!$     o = 1; if(ncomp == 2) o = 2
!!$     omaps = 0; omaps(:,o:,:) = reshape(maps, [m/ncomp,ncomp,nsim])
!!$     call write_map(omaps,pixels,nside,order,trim(prefix)//"_maps.hdf")
!!$     omaps = 0; omaps(:,o:,:) = reshape(rhss2, [m/ncomp,ncomp,nsim])
!!$     call write_map(omaps,pixels,nside,order,trim(prefix)//"_rhss.hdf")
!!$     call write_covmat(icov1, pixels, nside, order, ncomp, trim(prefix) // "_icov1.hdf")
!!$     call write_covmat(icov2, pixels, nside, order, ncomp, trim(prefix) // "_icov2.hdf")
!!$     call write_covmat(cov1,  pixels, nside, order, ncomp, trim(prefix) // "_cov1.hdf")
!!$     call write_covmat(cov2,  pixels, nside, order, ncomp, trim(prefix) // "_cov2.hdf")
!!$     call write_covmat(cov,  pixels, nside, order, ncomp, trim(prefix) // "_cov.hdf")
!!$     deallocate(rhss2, maps, icov1, icov2, cov1, cov2, omaps, cov)
!!$   end subroutine
!!$
!!$   subroutine covtest2
!!$     implicit none
!!$     character(len=512)        :: mapfile, covfile, ocovfile
!!$     integer(i4b)              :: nside, order, ncomp, i, j, n, o, nsim, bs, bi, bj, nb
!!$     real(dp),     allocatable :: imaps(:,:,:), maps(:,:), cov(:,:), ecov(:,:), ecov2(:,:)
!!$     real(dp),     allocatable :: invcov(:,:), invecov(:,:), smaps(:,:), root(:,:)
!!$     real(dp),     allocatable :: bchisq(:,:), bicov(:,:), eigs(:,:), eeigs(:,:)
!!$     real(dp),     allocatable :: fulleigs(:), fulleeigs(:)
!!$     integer(i4b), allocatable :: pixels(:), pixels2(:)
!!$     type(hdf_file)            :: hfile
!!$     type(zig_rng)             :: rng
!!$     call getarg(2, mapfile)
!!$     call getarg(3, covfile)
!!$     call getarg(4, ocovfile)
!!$     call read_map(imaps, pixels, nside, order, mapfile)
!!$     call read_covmat(cov, pixels2, nside, order, ncomp, covfile)
!!$     n = size(cov,1)
!!$     nsim = size(imaps,3)
!!$     o = 1; if(ncomp == 2) o = 2
!!$     allocate(maps(n,nsim))
!!$     maps = reshape(imaps(:,o:,:),[n,nsim])
!!$     deallocate(imaps)
!!$
!!$     ! Generer realisasjoner selv
!!$     allocate(smaps(n,nsim), root(n,n))
!!$     call eigen_pow(cov, 0.5d0, root)
!!$     call zig_init(rng, 1)
!!$     do i = 1, nsim
!!$        do j = 1, n; smaps(j,i) = zig_gauss(rng); end do
!!$        smaps(:,i) = matmul(root, smaps(:,i))
!!$     end do
!!$
!!$     allocate(ecov(n,n), ecov2(n,n))
!!$     ecov = 0
!!$     ecov2= 0
!!$     do i = 1, nsim
!!$        if(mod(i-1,100)==0) write(*,*) i-1
!!$        ecov = ecov + matmul(reshape(maps(:,i),[n,1]),reshape(maps(:,i),[1,n]))
!!$        ecov2= ecov2+ matmul(reshape(smaps(:,i),[n,1]),reshape(smaps(:,i),[1,n]))
!!$     end do
!!$     ecov = ecov/nsim
!!$     ecov2= ecov2/nsim
!!$
!!$     ! chisquares
!!$     allocate(invcov(n,n), invecov(n,n))
!!$     call eigen_pow(cov,  -1d0, invcov)
!!$!     do i = 1, nsim
!!$!     write(*,'(i4,2f9.4)') i, (dot_product(maps(:,i),matmul(invcov,maps(:,i)))-n)/sqrt(2.0*n), (dot_product(smaps(:,i),matmul(invcov,smaps(:,i)))-n)/sqrt(2.0*n)
!!$!     end do
!!$
!!$     write(*,*) "cov"
!!$     call dump_matrix(cov(:5,:5))
!!$     write(*,*) "ecov"
!!$     call dump_matrix(ecov(:5,:5))
!!$     write(*,*) "ecov2"
!!$     call dump_matrix(ecov2(:5,:5))
!!$     write(*,*) "diff"
!!$     call dump_matrix(cov(:5,:5)-ecov(:5,:5))
!!$     write(*,*) "ratio"
!!$     call dump_matrix(cov(:5,:5)/ecov(:5,:5))
!!$
!!$     ! Eigenvalues
!!$
!!$
!!$     ! Blockwise chisquares
!!$     bs = 32
!!$     nb = n/bs
!!$     allocate(bchisq(nb,nsim), bicov(bs,bs))
!!$     allocate(eigs(bs,nb),eeigs(bs,nb))
!!$     do bi = 1, nb
!!$        call eigen_pow(cov((bi-1)*bs+1:bi*bs,(bi-1)*bs+1:bi*bs), -1d0, bicov)
!!$        do i = 1, nsim
!!$           bchisq(bi,i) = dot_product(maps((bi-1)*bs+1:bi*bs,i),matmul(bicov,maps((bi-1)*bs+1:bi*bs,i)))
!!$        end do
!!$        call get_eigenvalues(cov((bi-1)*bs+1:bi*bs,(bi-1)*bs+1:bi*bs), eigs(:,bi))
!!$        call get_eigenvalues(ecov((bi-1)*bs+1:bi*bs,(bi-1)*bs+1:bi*bs), eeigs(:,bi))
!!$     end do
!!$     allocate(fulleigs(n), fulleeigs(n))
!!$     call get_eigenvalues(cov, fulleigs)
!!$     call get_eigenvalues(ecov, fulleeigs)
!!$
!!$     call open_hdf_file("bchisq.hdf", hfile, "w")
!!$     call write_hdf(hfile, "chisq", bchisq)
!!$     call write_hdf(hfile, "eigs", eigs)
!!$     call write_hdf(hfile, "eeigs", eeigs)
!!$     call write_hdf(hfile, "feigs", fulleigs)
!!$     call write_hdf(hfile, "feeigs", fulleeigs)
!!$     call close_hdf_file(hfile)
!!$     do i = 1, n
!!$        write(*,'(f6.1)',advance="no") (dot_product(maps(i,:),maps(i,:))/cov(i,i) - nsim)/sqrt(2d0*nsim)
!!$        if(mod(i,12) == 0) write(*,*)
!!$     end do
!!$     write(*,*)
!!$
!!$     call write_covmat(ecov, pixels2, nside, order, ncomp, ocovfile)
!!$   end subroutine
!!$
!!$   subroutine smoothtest
!!$     implicit none
!!$     real(dp)     :: var(4), map(4), scov(4,4), smap(4), w(2,2), scov2(4,4), scov3(4,4)
!!$     real(dp)     :: tcov(2,2,2,2), tmap(2,2), ecov(4,4), rcov(4,4), bcov(2,2,2)
!!$     integer(i4b) :: nside, order, pixels(2), i, j, n
!!$     type(pixinfo):: pinfo, pinfo2
!!$     type(zig_rng):: rng
!!$     nside  = 1
!!$     order  = ring
!!$     pixels = [0,1]
!!$     var    = [1,1,2,2]
!!$     map    = [1,0,0,1]
!!$     w      = transpose(reshape([&
!!$               & 0.55, 0.45, &
!!$               & 0.45, 0.55 ], [2,2]))
!!$     rcov   = transpose(reshape([&
!!$               & 1, 0, 0, 0, &
!!$               & 0, 1, 0, 0, &
!!$               & 0, 0, 2, 0, &
!!$               & 0, 0, 0, 2 ], [4,4]))
!!$     bcov   = reshape([1,1,0,0,0,0,2,2],[2,2,2])
!!$
!!$     call init_pixinfo(pinfo, nside, order, pi, pixels)
!!$     call smooth_map  (reshape(map,[2,2]), pinfo, w, tmap)
!!$     smap = reshape(tmap,[4])
!!$     call smooth_noise(reshape(var,[2,2]), pinfo, w, tcov)
!!$     scov = reshape(tcov,[4,4])
!!$     call smooth_noise(reshape(rcov,[2,2,2,2]), pinfo, w, tcov)
!!$     scov2 = reshape(tcov,[4,4])
!!$     call smooth_noise(bcov, pinfo, w, tcov)
!!$     scov3 = reshape(tcov,[4,4])
!!$
!!$     write(*,*) "map:"
!!$     call dump_matrix(smap)
!!$     write(*,*) "cov:"
!!$     call dump_matrix(scov)
!!$     write(*,*) "cov2:"
!!$     call dump_matrix(scov2)
!!$     write(*,*) "cov3:"
!!$     call dump_matrix(scov3)
!!$
!!$     n    = 1000000
!!$     call zig_init(rng, 1)
!!$     ecov = 0
!!$     do i  = 1, n
!!$        if(mod(i-1,100000)==0) write(*,*) i-1
!!$        do j = 1, size(var)
!!$           map(j) = zig_gauss(rng)*sqrt(var(j))
!!$        end do
!!$        call smooth_map  (reshape(map,[2,2]), pinfo, w, tmap)
!!$        smap = reshape(tmap,[4])
!!$        ecov = ecov + matmul(reshape(smap,[4,1]),reshape(smap,[1,4]))
!!$     end do
!!$     ecov = ecov/n
!!$     write(*,*) "ecov:"
!!$     call dump_matrix(ecov)
!!$
!!$
!!$     write(*,*) "pos (6,8)"
!!$     call init_pixinfo(pinfo2, 8, order, pi, irange(12*8**2)-1)
!!$     call dump_matrix(pinfo2%rot(:,:,7,6))
!!$
!!$
!!$   end subroutine
!!$
!!$   subroutine covtest
!!$     implicit none
!!$     character(len=512)        :: mapfile, covfile
!!$     real(dp),     allocatable :: imap(:,:), cov(:,:), map(:), eigvals(:), eigvecs(:,:)
!!$     real(dp),     allocatable :: emap(:)
!!$     integer(i4b), allocatable :: pixels(:), pixels2(:), n
!!$     integer(i4b)              :: nside, ncomp, order, status
!!$     call getarg(1, mapfile)
!!$     call getarg(2, covfile)
!!$     call read_map   (imap, pixels,  nside, order, mapfile)
!!$     call read_covmat(cov,  pixels2, nside, order, ncomp, covfile)
!!$     n = size(cov,1)
!!$     allocate(map(n), emap(n))
!!$     map = reshape(imap(:,2:3),[n])
!!$     allocate(eigvals(n), eigvecs(n,n))
!!$     call eigen_decomp(cov, eigvals, eigvecs, status)
!!$
!!$     emap = matmul(transpose(eigvecs), map)
!!$     call dump_matrix(reshape(eigvals,[n,1]), "eigvals.txt")
!!$     call dump_matrix(reshape(emap,   [n,1]), "emap.txt")
!!$   end subroutine
!!$
!!$   subroutine sdss
!!$     implicit none
!!$     character(len=512) :: ifile, ofile
!!$     integer(i4b)       :: unit, n, i
!!$     type(hdf_file)     :: hfile
!!$     real(dp)           :: dat(7)
!!$     real(dp), allocatable :: time(:), pos(:,:)
!!$     real(sp), allocatable :: mag(:,:)
!!$     call getarg(1, ifile)
!!$     call getarg(2, ofile)
!!$     unit = getlun()
!!$     open(unit,file=ifile,action="read",status="old")
!!$     n = 0
!!$     do
!!$        read(unit,*,end=1) dat
!!$        if(mod(n,1000000)==0) write(*,*) n
!!$        n = n+1
!!$     end do
!!$1    allocate(time(n), pos(n,2), mag(n,4))
!!$     rewind(unit)
!!$     do i = 1, n
!!$        if(mod(i-1,1000000)==0) write(*,*) i-1
!!$        read(unit,*) pos(i,:), mag(i,:), time(i)
!!$     end do
!!$     close(Unit)
!!$     call open_hdf_file(ofile, hfile, "w")
!!$     call write_hdf(hfile, "time", time)
!!$     call write_hdf(hfile, "pos",  pos)
!!$     call write_hdf(hfile, "mag",  mag)
!!$     call close_hdf_file(hfile)
!!$   end subroutine
!!$
!!$   subroutine targtest
!!$     implicit none
!!$     character(len=512)   :: knifestr, parfile, cmd, cfile, oname
!!$     integer(i4b)         :: cnum, nces, ndi
!!$     type(acceptlist)     :: alist
!!$     type(quiet_target)   :: target
!!$     type(quiet_ces_info) :: ces
!!$     type(hdf_file)       :: hfile
!!$     type(quiet_target), allocatable :: targets(:)
!!$     type(swiss_knife),  allocatable :: knife_defs(:)
!!$     integer(i4b),       allocatable :: knife_res(:,:,:)
!!$     real(sp),           allocatable :: cstat(:,:)
!!$     real(sp),           allocatable :: dstat(:,:,:)
!!$     call getarg(1, parfile)
!!$     call getarg(2, knifestr)
!!$     call initialize_ces_mod(parfile)
!!$     call initialize_target_mod(parfile)
!!$     call initialize_accept_list("*", alist)
!!$     call init_target(target, alist)
!!$     
!!$     call getarg(3, cmd)
!!$     call getarg(4, cfile)
!!$     nces = get_num_ces()
!!$     ndi  = size(quiet_diodes)
!!$     allocate(cstat(nces,STAT_NUM))
!!$     allocate(dstat(ndi,nces,NUM_DIODE_STATS))
!!$     select case(cmd)
!!$        case("build")
!!$           do cnum = 1, get_num_ces()
!!$              write(*,*) cnum
!!$              call get_ces_info(cnum, ces)
!!$              call open_hdf_file(ces%l3file, hfile, "r")
!!$              call read_hdf(hfile, "stats",       cstat(cnum,:))
!!$              call read_hdf(hfile, "diode_stats", dstat(:,cnum,:))
!!$              call close_hdf_file(hfile)
!!$              call free_ces_info(ces)
!!$           end do
!!$           call open_hdf_file(cfile, hfile, "w")
!!$           call write_hdf(hfile, "stats", cstat)
!!$           call write_hdf(hfile, "diode_stats", dstat)
!!$           call close_hdf_file(hfile)
!!$        case("use")
!!$           call getarg(5, oname)
!!$           call open_hdf_file(cfile, hfile, "r")
!!$           call read_hdf(hfile, "stats",       cstat)
!!$           call read_hdf(hfile, "diode_stats", dstat)
!!$           call close_hdf_file(hfile)
!!$           call jackknife(target, knifestr, cstat, dstat, targets, knife_defs, knife_res)
!!$           call print_knives(oname, knife_defs, knife_res, cid_list)
!!$        case default
!!$           call assert(.false., "Unknown command " // trim(cmd))
!!$     end select
!!$   end subroutine
!!$
!!$  subroutine print_knives(prefix, defs, groups, cids)
!!$    implicit none
!!$    type(swiss_knife), intent(in) :: defs(:)
!!$    integer(i4b),      intent(in) :: groups(:,:,:)
!!$    integer(i4b),      intent(in) :: cids(:)
!!$    integer(i4b)                  :: i, j, k, n, unit
!!$    character(len=512)            :: fname, prefix
!!$    n = size(defs)
!!$    unit = getlun()
!!$    do i = 1, n
!!$       if(prefix=="-") then
!!$          unit = stdout
!!$       else
!!$          fname = trim(prefix) // "/knife" // trim(itoa(i,3)) // "_" // trim(defs(i)%full_name) // ".txt"
!!$          call mkdirs(fname,.true.)
!!$          open(unit,file=fname)
!!$       end if
!!$       do j = 1, size(groups,2)
!!$          if(.not. any(groups(:,j,i) /= 0)) cycle
!!$          write(unit,'(i5)',advance="no") cids(j)
!!$          do k = 1, size(groups,1)
!!$             write(unit,'(i4)',advance="no") groups(k,j,i)
!!$          end do
!!$          write(unit,*)
!!$       end do
!!$       if(prefix /= "-") close(unit)
!!$    end do
!!$  end subroutine
!!$
!!$!   subroutine fullval
!!$!     implicit none
!!$!     character(len=512) :: arg, ifile
!!$!     integer(i4b)       :: unit, cmax, dimax, i, j, h, mod, mod2, di, di2
!!$!     type(obs_type)     :: obs, hobs(2), pobs(2)
!!$!     integer(i4b),   allocatable :: diobs(:)
!!$!     type(obs_type), allocatable :: obslist(:,:,:), foolist(:,:,:)
!!$!     call getarg(1, ifile)
!!$!     unit = getlun()
!!$!     cmax = 0; dimax = 0
!!$!     open(unit,file=ifile,action="read",status="old")
!!$!     do
!!$!        read(unit,*,end=1,err=4) obs
!!$!        cmax  = max(cmax,  obs%ces)
!!$!        dimax = max(dimax, obs%di)
!!$!        4 continue
!!$!     end do
!!$!1    rewind(unit)
!!$!     allocate(diobs(dimax), obslist(cmax,2,dimax), foolist(cmax,2,dimax))
!!$!     obslist%ces = 0
!!$!     foolist%ces = 0
!!$!     diobs = 0
!!$!     do
!!$!        read(unit,*,end=2,err=5) obs
!!$!        if(obs%horn == 1) diobs(obs%di) = diobs(obs%di) + 1
!!$!        obslist(diobs(obs%di),obs%horn, obs%di) = obs
!!$!        foolist(obs%ces,      obs%horn, obs%di) = obs
!!$!        5 continue
!!$!     end do
!!$!2    close(unit)
!!$!
!!$!     do di = 1, dimax
!!$!        do j = 1, diobs(di)
!!$!           hobs = obslist(j,:,di)
!!$!           ! Reject if horn ratio is not -1
!!$!           if(abs(hobs(2)%frac+1) > 0.25) cycle
!!$!           ! Reject if gain is too high
!!$!           if(any(abs(hobs%amp) > 40)) cycle
!!$!           if(any(abs(hobs%amp) == 0)) cycle
!!$!           ! Reject if positions are too far away from center
!!$!           if(polangdist([pi/2,0]-hobs(1)%p([2,1])*DEG2RAD, &
!!$!            & [pi/2,0]-hobs(1)%p0([2,1])*DEG2RAD) > 0.6*DEG2RAD) cycle
!!$!           if(polangdist([pi/2,0]-hobs(2)%p([2,1])*DEG2RAD, &
!!$!            & [pi/2,0]-hobs(2)%p0([2,1])*DEG2RAD) > 0.6*DEG2RAD) cycle
!!$!           mod  = (di-1)/4
!!$!           if(modulo(mod,2) == 0) then; mod2 = mod-1; else; mod2 = mod+1; end if
!!$!           di2  = mod2*4+1
!!$!           pobs = foolist(hobs(1)%ces,:,di2)
!!$!           if(pobs(1)%ces /= hobs(1)%ces) cycle
!!$!
!!$!           ! Reject if horns don't agree
!!$!           if(polangdist([pi/2,0]-hobs(1)%p([2,1])*DEG2RAD, &
!!$!            & [pi/2,0]-pobs(1)%p([2,1])*DEG2RAD) > 7*DEG2RAD/60) cycle
!!$!
!!$!           ! And output the current entry
!!$!           3 format(i5,i4,i3,f14.7,6f10.4,3e15.7,2f10.4,6e15.7,f10.4)
!!$!           do h = 1, 2
!!$!              write(*,3) hobs(h)%ces, hobs(h)%di, hobs(h)%horn, hobs(h)%mjd, &
!!$!               & hobs(h)%p0, hobs(h)%p, hobs(h)%tamp, hobs(h)%amp, hobs(h)%frac, &
!!$!               & hobs(h)%dist, hobs(h)%foo, hobs(h)%ddk
!!$!           end do
!!$!        end do
!!$!     end do
!!$!   end subroutine
!!$
!!$   subroutine noisestuff
!!$     implicit none
!!$     character(len=512) :: arg, ifile, ofile
!!$     integer(i4b)       :: nsim, nside, order, sim, i, j
!!$     type(zig_rng)      :: rng
!!$     integer(i4b), allocatable :: pixels(:)
!!$     real(dp),     allocatable :: imap(:,:), omap(:,:,:)
!!$     call getarg(1, ifile)
!!$     call getarg(2, arg); read(arg,*) nsim
!!$     call getarg(3, ofile)
!!$     call read_map(imap, pixels, nside, order, ifile)
!!$     allocate(omap(size(pixels),size(imap,2),nsim))
!!$     call zig_init(rng, 1)
!!$     do sim = 1, nsim
!!$        do j = 1, size(imap,2)
!!$           do i = 1, size(pixels)
!!$              omap(i,j,sim) = zig_gauss(rng) * imap(i,j)
!!$           end do
!!$        end do
!!$     end do
!!$     call write_map(omap, pixels, nside, order, ofile)
!!$     deallocate(imap, omap, pixels)
!!$   end subroutine
!!$
!!$   subroutine noisestuff2
!!$     implicit none
!!$     character(len=512) :: arg, ifile, ofile
!!$     integer(i4b)       :: nsim, nside, order, i, j, n, m, sim, status
!!$     type(zig_rng)      :: rng
!!$     integer(i4b), allocatable :: pixels(:)
!!$     real(dp),     allocatable :: icov(:,:,:,:), omap(:,:,:), eigvals(:), eigvecs(:,:)
!!$     real(dp),     allocatable :: cov(:,:)
!!$     real(dp),     allocatable :: vec(:)
!!$     call getarg(1, ifile)
!!$     call getarg(2, arg); read(arg,*) nsim
!!$     call getarg(3, ofile)
!!$     call read_covmat(icov, pixels, nside, order, ifile)
!!$     call zig_init(rng, 1)
!!$     n = size(pixels)
!!$     m = size(icov,2)
!!$     allocate(omap(n,m,nsim),vec(n*m), cov(n*m,n*m))
!!$     allocate(eigvals(n*m),eigvecs(n*m,n*m))
!!$     cov = reshape(icov,[n*m,n*m])
!!$     call eigen_decomp(cov, eigvals, eigvecs, status)
!!$     call assert(status == 0, "Status is: " // trim(itoa(status)))
!!$     where(eigvals < maxval(eigvals)/1d14) eigvals = 0
!!$     do sim = 1, nsim
!!$        do i = 1, n*m
!!$           vec(i) = zig_gauss(rng)*sqrt(eigvals(i))
!!$        end do
!!$        omap(:,:,sim) = reshape(matmul(eigvecs, vec),[n,m])
!!$     end do
!!$     call write_map(omap, pixels, nside, order, ofile)
!!$     deallocate(icov, omap, pixels, vec, eigvals, eigvecs, cov)
!!$   end subroutine
!!$
!!$   subroutine pointtest3
!!$     implicit none
!!$     character(len=512) :: parfile, arg
!!$
!!$     real(dp) :: ddk, daz, del  ! encoder offset
!!$     real(dp) :: phi_c, theta_c ! collimation offset
!!$     real(dp) :: kf             ! flex
!!$     real(dp) :: etilt          ! elevation tilt
!!$     real(dp) :: omegaf, thetaf ! azimuth fixed tilt
!!$     real(dp) :: X, Y           ! tilt meter readings
!!$     real(dp) :: omegat, thetat ! azimuth dynamic tilt
!!$     real(dp) :: enc(3), mjd, phi, theta, psi, p1(3), p2(3), p3(3), t1, t2
!!$     integer(i4b) :: di, mod, i
!!$     real(dp), dimension(3,3) :: mdi, mhorn, mcol, mdk, mel, metilt, &
!!$      & maz, matilt1, matilt2, mflex, C, B, B1, B2, B3, mat
!!$     real(dp) :: params(mpar_max)
!!$     type(rotinfo) :: info
!!$
!!$     call getarg(1, parfile)
!!$     call initialize_quiet_pointing_mod(parfile, .false., .true.)
!!$     call init_quiet_pointing_mod2(parfile)
!!$
!!$     ! Our test mount model
!!$     params = 0
!!$     params(mpar_denc)  = [-44d0,10d0,100d0]
!!$     params(mpar_flex)  =  0.1
!!$     params(mpar_etilt) = -0.1
!!$     params(mpar_atilt) = [270.0,0.1]
!!$     params(mpar_col)   = [0.0,0.1]
!!$     params = params*DEG2RAD
!!$
!!$     omegat = 0
!!$     thetat = 0
!!$
!!$     daz     = params(mpar_denc(1))
!!$     del     = params(mpar_denc(2))
!!$     ddk     = params(mpar_denc(3))
!!$     etilt   = params(mpar_etilt)
!!$     kf      = params(mpar_flex)
!!$     omegaf  = params(mpar_atilt(1))
!!$     thetaf  = params(mpar_atilt(2))
!!$     phi_c   = params(mpar_col(1))
!!$     theta_c = params(mpar_col(2))
!!$
!!$     ! And test parameters
!!$     di     = 1; mod = quiet_diodes(di)%horn
!!$     enc    = [0d0,45d0,120d0]*DEG2RAD
!!$     mjd    = 55500
!!$
!!$     ! Use auto model
!!$     call init_rotinfo(info)
!!$     info%enc = enc
!!$     info%mjd = mjd
!!$     info%mod = mod
!!$     info%di  = di
!!$
!!$     call set_mount_simple(info, params)
!!$     call wall_time(t1)
!!$     do i = 1, 1!100
!!$        info%enc(1) = pi*i/400
!!$        mat = calc_rot(info, [rdef_bore])
!!$        do di = 1, 1!size(quiet_diodes)
!!$           info%di = di
!!$           info%mod = quiet_diodes(di)%horn
!!$           mat = calc_rot(info, [rdef_det,-rdef_bore,rdef_hor2equ])
!!$           !mat = calc_rot(info, [rdef_det,-rdef_bore,rdef_atm,rdef_hor2app,rdef_app2equ])
!!$           call convert_euler_matrix_to_angles_zyz(mat, phi, theta, psi)
!!$           p1 = [phi,theta,psi]
!!$        end do
!!$     end do
!!$     call wall_time(t2)
!!$     write(*,'(a8,4f15.8)') "A", p1*RAD2DEG, t2-t1
!!$
!!$     ! The corresponding with the existing model
!!$     call set_mount_override(.true., [daz,del,ddk,kf,omegaf,thetaf,etilt,theta_c,phi_c,0d0])
!!$     call wall_time(t1)
!!$     do i = 1, 1!100
!!$        enc(1) = pi*i/400
!!$        do di = 1, 1!size(quiet_diodes)
!!$           mod = quiet_diodes(di)%horn
!!$           call coord_convert(coord_tele, enc(1), enc(2), enc(3), &
!!$            & coord_equ, p2(1), p2(2), p2(3), mod=mod, diode=modulo(di-1,4), mjd=mjd)
!!$        end do
!!$     end do
!!$     call wall_time(t2)
!!$     write(*,'(a8,4f15.8)') "B", p2*RAD2DEG, t2-t1
!!$
!!$     write(*,'(a8,3f15.8)') "diff", (p1-p2)*RAD2DEG*60
!!$   end subroutine
!!$
!!$   !subroutine pointtest2
!!$   !  implicit none
!!$   !  character(len=512) :: parfile, arg
!!$   !  integer(i4b), parameter :: npar = 8
!!$   !  real(dp)     :: mjd, p(3), q(3), params(npar)
!!$   !  real(dp), dimension(3,3) :: boresight, mcorr, coll, det, hmat, tmat
!!$   !  integer(i4b) :: mod
!!$   !  call getarg(1, parfile)
!!$   !  call initialize_quiet_pointing_mod(parfile, .false.)
!!$   !  ! Set up a mount model with only collimation
!!$   !  params = 0
!!$   !  call getarg(2, arg); read(arg,*) params(6)
!!$   !  call getarg(3, arg); read(arg,*) params(7)
!!$
!!$   !  mjd = 55097.97
!!$   !  p = [ 0.5909713E+00,  0.2278542E+01, 0.5235549E+00 ]
!!$
!!$   !  call set_mount_override(.true., params)
!!$   !  coll = rot_collimate(mjd)
!!$   !  do mod = 0, 90
!!$   !     det  = rot_module2boresight(mod)
!!$   !     hmat = angles2rot(p(1),p(2),p(3))
!!$   !     tmat = matmul(hmat,matmul(transpose(det),transpose(coll)))
!!$   !     !hmat = matmul(tmat,matmul(coll,det))
!!$   !     hmat = matmul(tmat,det)
!!$   !     ! Ah, det kommuterer ikke selv med bare collimation. Hvis
!!$   !     ! vi går fram og tilbake, med collimation ene veien men ikke
!!$   !     ! andre, får vi:
!!$   !     ! hmat * det' * coll' * det /= hmat * coll'
!!$   !     call rot2angles(hmat, q(1), q(2), q(3))
!!$   !     write(*,'(i4,3f13.7)') mod, q-p
!!$   !  end do
!!$   !end subroutine
!!$
!!$!   subroutine leak
!!$!     implicit none
!!$!     character(len=512) :: parfile
!!$!     real(dp)           :: foo
!!$!     
!!$!     call getarg(1, parfile)
!!$!     call initialize_target_mod(parfile)
!!$!     write(*,*) "A", get_leak_dis()
!!$!     write(*,*) "B", get_leak_mods()
!!$!   end subroutine
!!$
!!$   subroutine test
!!$     use quiet_utils
!!$     use alm_tools
!!$     use quiet_healpix_mod
!!$     implicit none
!!$     integer(i4b)  :: nside, order, lmax, ncomp, l, i
!!$     real(dp),     dimension(:,:),   allocatable :: map, cls
!!$     complex(dpc), dimension(:,:,:), allocatable :: alms
!!$     nside = 256
!!$     order = 1
!!$     ncomp = 1
!!$     lmax  = 3*nside
!!$     allocate(map(0:12*nside**2-1,ncomp))
!!$     allocate(alms(ncomp,0:lmax,0:lmax))
!!$     allocate(cls(0:lmax,ncomp*(ncomp+1)/2))
!!$     map = 1/sqrt(4*pi)
!!$     if(ncomp == 1) then
!!$        call map2alm(nside,lmax,lmax,map(:,1),alms,[-1d0,1d0],get_hpix_ringweights(nside))
!!$     elseif(ncomp == 3) then
!!$        call map2alm(nside,lmax,lmax,map,alms,[-1d0,1d0],get_hpix_ringweights(nside))
!!$     end if
!!$     call alm2cl(lmax,lmax,alms,cls)
!!$     call dump_matrix(cls)
!!$   end subroutine
!!$
!!$   !subroutine pointtest
!!$   !  implicit none
!!$   !  character(len=512) :: parfile
!!$   !  real(dp)     :: mjd, p(3), boresight(3,3), mcorr(3,3), coll(3,3), det(3,3)
!!$   !  real(dp)     :: mat(3,3),mat2(3,3), boremount(3,3), bore2(3,3)
!!$   !  integer(i4b) :: mod
!!$   !  call getarg(1, parfile)
!!$   !  call initialize_quiet_pointing_mod(parfile)
!!$   !  mjd = 55550
!!$   !  p   = [45,45,45]*DTOR
!!$   !  mod = 2
!!$   !  call swap_coordinate_convention(p(1),p(2),p(3),coord_tele)
!!$   !  boresight = angles2rot(p(1),p(2),p(3))
!!$   !  mcorr     = rot_mount_correction(mjd, p(1), p(2), p(3), mod)
!!$   !  coll      = rot_collimate(mjd)
!!$   !  det       = rot_module2boresight(mod)
!!$
!!$   !  mat = matmul(mcorr,matmul(boresight,matmul(coll,det)))
!!$   !  mat2 = angles2hor(mjd, p(1), p(2), p(3), mod)
!!$
!!$   !  boremount = matmul(mat,transpose(matmul(coll,det)))
!!$   !  bore2     = matmul(rot_mount_invert(mjd, mod, boremount), boremount)
!!$
!!$   !  call rot2angles(bore2, p(1), p(2), p(3))
!!$   !  call swap_coordinate_convention(p(2),p(2),p(3),coord_tele)
!!$   !  write(*,'(3f12.7)') p*RTOD
!!$
!!$   !  call hor2angles(mat2, mjd, p(1), p(2), p(3), mod)
!!$   !  call swap_coordinate_convention(p(1),p(2),p(3),coord_tele)
!!$   !  write(*,'(3f12.7)') p*RTOD
!!$   !end subroutine
!!$
!!$   subroutine splinetest
!!$    implicit none
!!$    real(dp),       allocatable :: x(:), ox(:), variance(:), y(:), y2(:), oy(:)
!!$    integer(i4b)                :: d1, d2, ndi, n, i, m
!!$    n   = 10000
!!$    m   = 25
!!$    allocate(x(m),variance(m),y(m),y2(m),ox(n),oy(n))
!!$    x  = 10*real(irange(m),dp)/m-5
!!$    ox = 10*real(irange(n),dp)/n-5
!!$    y  = x**0-2*x+1
!!$    call spline(x, y, 1d30, 1d30, y2)
!!$    open(40,file="foo.txt")
!!$    do i = 1, m
!!$       write(40,'(3e15.7)') x(i), y(i), y2(i)
!!$    end do
!!$    close(40)
!!$    do i = 1, n
!!$       oy(i) = splint_plain(x, y, y2, ox(i))
!!$    end do
!!$    open(40,file="bar.txt")
!!$    do i = 1, n
!!$       write(40,'(2e15.7)') ox(i), oy(i)
!!$    end do
!!$    close(40)
!!$   end subroutine
!!$
!!$   subroutine hdftest
!!$     implicit none
!!$     character(len=512) :: fname
!!$     type(hdf_file)     :: file
!!$     integer(i4b)       :: i, j, n
!!$     integer(i4b), allocatable :: arr(:,:)
!!$     call getarg(1, fname)
!!$     n = 8
!!$     allocate(arr(n,n))
!!$     do j = 0, n-1
!!$     do i = 0, n-1
!!$        arr(i+1,j+1) = i+j*10
!!$     end do
!!$     end do
!!$     call open_hdf_file(fname, file, "w")
!!$     call write_hdf(file, "foo", arr)
!!$
!!$     call create_hdf_set(file, "bar", shape(arr), H5T_STD_I32LE)
!!$     do j = 0, 1
!!$     do i = 0, 1
!!$        call write_hdf(file, "bar", slice([n/2*i+1,n/2*i+n/2],[n/2*j+1,n/2*j+n/2]), arr(n/2*i+1:n/2*i+n/2,n/2*j+1:n/2*j+n/2))
!!$     end do
!!$     end do
!!$
!!$     call close_hdf_file(file)
!!$   end subroutine
!!$
!!$   subroutine gaintest
!!$     implicit none
!!$     character(len=512)    :: parfile, arg, ofname
!!$     integer(i4b)          :: di, n, i, j, ndi, nmax, cnum, m
!!$     real(dp)              :: mjd(2), dt
!!$     type(hdf_file)        :: hfile
!!$     type(quiet_ces_info)  :: ces
!!$     real(dp), allocatable :: gains(:,:), times(:)
!!$     call getarg(1, parfile)
!!$     call getarg(2, arg); read(arg,*) mjd(1)
!!$     call getarg(3, arg); read(arg,*) mjd(2)
!!$     call getarg(4, arg); read(arg,*) dt
!!$     call getarg(5, ofname)
!!$     call initialize_gain_mod(parfile)
!!$     call initialize_ces_mod(parfile)
!!$     nmax = (mjd(2)-mjd(1))/dt
!!$     ndi = size(quiet_diodes)
!!$     allocate(times(nmax))
!!$     ! Loop through each CES get gain samples
!!$     n = 0
!!$     do j = 1, size(cid_sort)
!!$        cnum = lookup_ces(cid_sort(j))
!!$        call get_ces_info(cnum, ces)
!!$        m = (ces%mjd(2)-ces%mjd(1))/dt
!!$        do i = 1, m
!!$           times(n+i) = ces%mjd(1) + dt*(i-1)
!!$        end do
!!$        n = n+m
!!$        call free_ces_info(ces)
!!$     end do
!!$     ! Then get all these gains
!!$     allocate(gains(n,ndi))
!!$     do di = 1, ndi
!!$        call get_gains(times(1:n), di, gains(:,di))
!!$     end do
!!$     call open_hdf_file(ofname, hfile, "w")
!!$     call write_hdf(hfile, "time", times(1:n))
!!$     call write_hdf(hfile, "gains", gains)
!!$     call close_hdf_file(hfile)
!!$     deallocate(gains, times)
!!$   end subroutine
!!$
!!$   subroutine corrtest2
!!$     implicit none
!!$     character(len=512) :: str, fname
!!$     complex(spc),allocatable :: ft(:,:)
!!$     real(sp)           :: corr(2,2)
!!$     integer(i4b)       :: n, m, i, j
!!$     type(lx_struct)    :: data
!!$     call getarg(1, fname)
!!$     call read_l3_file(fname, data)
!!$
!!$     n = 1000
!!$     allocate(ft(size(data%tod,1)/2+1,size(data%tod,2)))
!!$     call fft_multi(data%tod(:,1:2), ft(:,1:2), 1)
!!$     m = size(ft,1)
!!$     ft(1:n-1,:) = 0
!!$     call fft_multi(data%tod(:,1:2), ft(:,1:2), -1)
!!$     write(*,*) sum(data%tod(:,1)*data%tod(:,2))/size(data%tod,1)/sqrt(sum(data%tod(:,1)**2)/size(data%tod,1)*sum(data%tod(:,2)**2)/size(data%tod,1))
!!$
!!$
!!$     write(*,*) "rr", sum(real(ft(n:,1))*real(ft(n:,2))/sqrt(real(ft(n:,1))**2*real(ft(n:,2))**2))/size(ft(n:,1))
!!$     write(*,*) "ri", sum(real(ft(n:,1))*aimag(ft(n:,2))/sqrt(real(ft(n:,1))**2*aimag(ft(n:,2))**2))/size(ft(n:,1))
!!$     write(*,*) "ir", sum(aimag(ft(n:,1))*real(ft(n:,2))/sqrt(aimag(ft(n:,1))**2*real(ft(n:,2))**2))/size(ft(n:,1))
!!$     write(*,*) "ri1", sum(real(ft(n:,1))*aimag(ft(n:,1))/sqrt(real(ft(n:,1))**2*aimag(ft(n:,1))**2))/size(ft(n:,1))
!!$     write(*,*) "ir1", sum(aimag(ft(n:,1))*real(ft(n:,1))/sqrt(aimag(ft(n:,1))**2*real(ft(n:,1))**2))/size(ft(n:,1))
!!$     write(*,*) "ii", sum(aimag(ft(n:,1))*aimag(ft(n:,2))/sqrt(aimag(ft(n:,1))**2*aimag(ft(n:,2))**2))/size(ft(n:,1))
!!$
!!$     call measure_fft_corr(ft(n:m-1,1:2), corr)
!!$     write(*,*) "m"
!!$     call dump_matrix(real(corr,dp))
!!$
!!$     open(40,file="bar.txt")
!!$     do i = n, size(ft,1)-1
!!$        write(40,'(i6,16e15.7)') i, real(sqrt(ft(i,1:2)*conjg(ft(i,1:2)))), atan2(aimag(ft(i,1:2)),real(ft(i,1:2))), real(ft(i,1:2)), aimag(ft(i,1:2)), real(sqrt(ft(i+1,1:2)*conjg(ft(i+1,1:2)))), atan2(aimag(ft(i+1,1:2)),real(ft(i+1,1:2))), real(ft(i+1,1:2)), aimag(ft(i+1,1:2))
!!$
!!$     end do
!!$     close(40)
!!$
!!$     call sim_fft_corr(corr, ft(n:m-1,1:2))
!!$     call measure_fft_corr(ft(n:m-1,1:2), corr)
!!$     write(*,*) "corr again"
!!$     call dump_matrix(real(corr,dp))
!!$
!!$     !write(*,*) "real", sum(real(ft(n:,1)*conjg(ft(n:,2)))/sqrt(real(ft(n:,1)*conjg(ft(n:,1))*ft(n:,2)*conjg(ft(n:,2)))))/size(ft(n:,1))
!!$     !write(*,*) "imag", sum(aimag(ft(n:,1)*conjg(ft(n:,2)))/sqrt(real(ft(n:,1)*conjg(ft(n:,1))*ft(n:,2)*conjg(ft(n:,2)))))/size(ft(n:,1))
!!$     call fft_multi(data%tod(:,1:2), ft(:,1:2), -1)
!!$
!!$     write(*,*) sum(data%tod(:,1)*data%tod(:,2))/size(data%tod,1)/sqrt(sum(data%tod(:,1)**2)/size(data%tod,1)*sum(data%tod(:,2)**2)/size(data%tod,1))
!!$   end subroutine
!!$
!!$   subroutine measure_fft_corr(ft, res)
!!$     implicit none
!!$     complex(spc) :: ft(:,:)
!!$     real(sp)     :: res(:,:), t
!!$     integer(i4b) :: i, j
!!$     do i = 1, size(ft,2)
!!$        do j = i, size(ft,2)
!!$           res(i,j) = sum(real(ft(:,i)*conjg(ft(:,j))))/size(ft,1)
!!$           res(j,i) = res(i,j)
!!$        end do
!!$        t = res(i,i)**(-0.5)
!!$        res(i,:) = res(i,:) * t
!!$        res(:,i) = res(:,i) * t
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine sim_fft_corr(corr, ft)
!!$     implicit none
!!$     complex(spc) :: ft(:,:)
!!$     real(sp)     :: corr(:,:)
!!$     real(dp)     :: L(size(ft,2),size(ft,2))
!!$     integer(i4b) :: i, j
!!$     type(zig_rng):: rng
!!$     call zig_init(rng,1)
!!$     i = 0
!!$     call cholesky_decompose(real(corr,dp), L, i)
!!$     do i = 1, size(ft,1)
!!$        do j = 1, size(ft,2)
!!$           ft(i,j) = cmplx(zig_gauss(rng),zig_gauss(rng))/sqrt(2d0)
!!$        end do
!!$        ft(i,:) = matmul(L, ft(i,:))
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine zigtest
!!$     implicit none
!!$     type(zig_rng)            :: rng
!!$     integer(i4b), parameter  :: n = 100000, m = n/2+1
!!$     integer(i4b)             :: i
!!$     real(dp)                 :: arr(n)
!!$     complex(dpc)             :: ft(m)
!!$     call zig_init(rng, 1)
!!$     do i = 1, n
!!$        arr(i) = zig_gauss(rng)
!!$     end do
!!$     call fft(arr, ft, 1)
!!$     ft = ft*conjg(ft)
!!$     call fft(arr, ft, -1)
!!$     do i = 1, n
!!$        write(*,*) i, arr(i)
!!$     end do
!!$   end subroutine
!!$
!!$   subroutine corrtest
!!$     implicit none
!!$     character(len=512) :: arg, l3file
!!$     type(hdf_file)     :: file
!!$     integer(i4b)       :: ndi, d1, d2, i, n, ng, g, j, k, m, nbin
!!$     real(dp)           :: nu, v
!!$     integer(i4b), allocatable :: ext(:)
!!$     real(dp),     allocatable :: corr(:,:,:), eig(:)
!!$     call getarg(1, l3file)
!!$     call open_hdf_file(l3file, file, "r")
!!$     call get_size_hdf(file, "corr", ext)
!!$     allocate(corr(ext(1),ext(2),ext(3)))
!!$     call read_hdf(file, "corr", corr)
!!$     ndi = size(corr,2)
!!$     nbin= size(corr,1)
!!$     allocate(eig(ndi))
!!$     do i = 1, nbin
!!$        call get_eigenvalues(corr(i,:,:),eig)
!!$        write(*,*) maxval(eig), minval(eig)
!!$     end do
!!$     deallocate(corr,eig)
!!$   end subroutine
!!$
!!$   subroutine bintest
!!$     implicit none
!!$     character(len=512)        :: arg
!!$     integer(i4b)              :: n, m, i
!!$     real(dp)                  :: b
!!$     integer(i4b), allocatable :: bins(:,:)
!!$     call getarg(1,arg); read(arg,*) n
!!$     call getarg(2,arg); read(arg,*) m
!!$     allocate(bins(2,m))
!!$     call make_exp_bins(n, bins)
!!$     do i = 1, m
!!$        write(*,*) i, bins(:,i)
!!$     end do
!!$   end subroutine
!!$
!!$   !subroutine nantest
!!$   !  implicit none
!!$   !  real(dp) :: a, b, c, d
!!$   !  write(*,*) "A"
!!$   !  a = 0.0/0.0
!!$   !  write(*,*) "B", a
!!$   !  b = 1 + snan
!!$   !  write(*,*) "C", b
!!$   !  call fe_enable(fe_nan)
!!$   !  a = 0.0/0.0
!!$   !  write(*,*) "D", a
!!$   !  b = 1 + snan
!!$   !  write(*,*) "E", b
!!$   !end subroutine
!!$
!!$  ! Like coord convert, but the sys2 psi angle is known, while the sys1 psi angle isn't
!!$  subroutine cc_invdk(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, mod)
!!$    implicit none
!!$    integer(i4b) :: sys1, sys2, mod, err
!!$    real(dp)     :: phi1, phi2, theta1, theta2, psi1, psi2, mjd, p(1), foo
!!$    p = 0
!!$    powell_dk_real = [ phi1, theta1, psi2, mjd ]
!!$    powell_dk_int  = [ sys1, sys2, mod ]
!!$    call powell(p, powell_calc_dk, err)
!!$    psi1 = p(1)
!!$    call coord_convert(sys1, phi1, theta1, psi1, sys2, phi2, theta2, foo, mjd, mod)
!!$  end subroutine
!!$
!!$  function powell_calc_dk(p) result(res)
!!$    use healpix_types
!!$    implicit none
!!$    real(dp), dimension(:), intent(in), optional :: p
!!$    real(dp)                                     :: foo(2), psi2, res
!!$    call coord_convert(powell_dk_int(1), powell_dk_real(1), powell_dk_real(2), &
!!$     & p(1), powell_dk_int(2), foo(1), foo(2), psi2, powell_dk_real(4), powell_dk_int(3))
!!$    res = ang_diff(psi2, powell_dk_real(3))**2
!!$  end function
!!$
!!$  subroutine pmatang(mat)
!!$    implicit none
!!$    real(dp) :: mat(3,3), p(3)
!!$    call convert_euler_matrix_to_angles_zyz(mat, p(1),p(2),p(3))
!!$    write(*,'(3f15.7)') p*RAD2DEG
!!$  end subroutine
!!$
!!$  function get_observed_object_temperature(object, fwhm, mjd, nu0)
!!$    implicit none
!!$
!!$    character(len=*), intent(in) :: object
!!$    real(dp),         intent(in) :: mjd, fwhm, nu0
!!$    real(dp)                     :: get_observed_object_temperature
!!$
!!$    integer(i4b)       :: id
!!$    real(dp)           :: omega_b, omega_p, d_fid, omega_ref, f_A, f_d, d, sigma
!!$    real(dp)           :: D_w, R_pole, R_equ
!!$    real(dp), dimension(2) :: pos
!!$    real(dp), dimension(2) :: r_jup
!!$    real(dp), dimension(5) :: nu, T_p, T_p2
!!$
!!$    id       = name2eph(trim(object))
!!$    d        = ephem_dist(id, mjd)
!!$    sigma    = fwhm / 60.d0 / sqrt(8.d0*log(2.d0)) * DTOR
!!$    omega_b  = 2.d0*pi*sigma**2
!!$
!!$    if (trim(object) == 'jupiter') then
!!$       nu        = [22.85d0, 33.11d0, 40.92d0, 60.41d0, 93.25d0]
!!$       T_p       = [136.2d0, 147.2d0, 154.4d0, 165.d0,  173.d0]
!!$       call spline(nu, T_p, 1.d30, 1.d30, T_p2)
!!$       omega_ref = 2.481d-8
!!$       d_fid     = 5.2d0
!!$       D_w       = 0.d0
!!$       R_pole    = 66854.d0 ! in km; Weiland et al. 2011
!!$       R_equ     = 71492.d0 ! in km
!!$    else
!!$       get_observed_object_temperature = 0.d0
!!$       return
!!$    end if
!!$    
!!$    f_d = (d/d_fid)**2
!!$    f_A = 1.d0                 ! D_w = 0.d0 for now; sub-0.1% error
!!$
!!$!    get_observed_object_temperature = T_p * omega_ref/omega_b * f_A / f_d
!!$    get_observed_object_temperature = splint(nu, T_p, T_p2, nu0) * omega_ref/omega_b * f_A / f_d
!!$  end function get_observed_object_temperature

end program
