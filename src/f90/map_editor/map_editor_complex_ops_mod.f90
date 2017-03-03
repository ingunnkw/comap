module map_editor_complex_ops_mod
  use healpix_types
  use pix_tools
  use alm_tools
  use fitstools
  use rngmod
  use map_editor_utils
  use math_tools
  use quiet_utils
  use powell_mod
  implicit none

  real(dp), allocatable, dimension(:), private :: d, nu
  real(dp), parameter, private :: k_B      = 1.3806503d-23
  real(dp), parameter, private :: h        = 1.0545726691251021d-34 * 2.d0*pi !6.626068d-34

contains

  subroutine fit_line
    implicit none

    character(len=512) :: filename, string
    integer(i4b) :: i, j, k, n, setsize, nset
    real(dp) :: sigma_x, sigma_y
    real(dp), allocatable, dimension(:,:) :: data, data_red
    real(dp), allocatable, dimension(:)   :: a

    real(dp) :: C, b, D, V1, V2, sigma_a
    type(planck_rng) :: handle
    
    if (iargc() /= 5) then
       write(*,*) 'Usage: map_editor fit_line_to_ASCII_data [filename] [numpt] [sample size] [num sets]'
       stop
    end if
    
    call getarg(2, filename)
    call getarg(3, string)
    read(string,*) n
    call getarg(4, string)
    read(string,*) setsize
    call getarg(5, string)
    read(string,*) nset

    call rand_init(handle, 182941)

    allocate(data(n,2), data_red(setsize,2), a(nset))
    open(58,file=trim(filename))
    do i = 1, n
       read(58,*) data(i,1), data(i,2)
    end do
    close(58)

    do i = 1, nset
       do j = 1, setsize
          k = max(min(int(rand_uni(handle)*n)+1,n),1)
          data_red(j,:) = data(k,:)
       end do
       C  = mean(data_red(:,1)*data_red(:,2)) - mean(data_red(:,1))*mean(data_red(:,2))
       V2 = mean(data_red(:,2)**2)-mean(data_red(:,2))**2
       V1 = mean(data_red(:,1)**2)-mean(data_red(:,1))**2
       D  = (V1-V2)/(2.d0*C)
       a(i)  = D + C/abs(C) * sqrt(1-D**2)
    end do

    write(*,*) 'Slope = ', mean(a), '+/-', sqrt(variance(a))

    
  end subroutine fit_line
  
  subroutine print_scaled_gal_avg
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x, alpha, alpha_min, alpha_max, dalpha, dl, db, theta, phi
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr, n, comp, m, b
    integer(i4b) :: listpix(0:1000000), nlist
    real(dp), allocatable, dimension(:,:) :: maps, par
    real(dp), allocatable, dimension(:,:) :: map1, map2, diff
    real(dp), allocatable, dimension(:)   :: tot, numpt
    character(len=512) :: infile1, infile2, outprefix, temp
    character(len=10)  :: alpha_text
    character(len=80), dimension(180)         :: header

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,outprefix)

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps), diff(0:npix-1,nmaps))
    call read_bintab(infile1, map1, npix, nmaps, nullval, anynull)
    call read_bintab(infile2, map2, npix, nmaps, nullval, anynull)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map1(:,i))
          call convert_nest2ring(nside, map2(:,i))
       end do
    end if
    
    comp      = 1
    alpha_min = 0.88d0
    alpha_max = 1.12d0
    n         = 9
    dalpha    = (alpha_max-alpha_min)/(n-1)
    dl        = 1.d0*pi/180.d0 ! In degrees
    m         = 360.d0 / dl
    db        = 1.5d0*pi/180.d0
    
    allocate(tot(-180:179), numpt(-180:179))
    open(58,file=trim(outprefix)//'_mean.dat')
    do i = 1, n
       alpha = alpha_min + (i-1) * dalpha
       diff = 0.d0
       where (map1 /= -1.6375d30 .and. map2 /= -1.6375d30)
          diff  = map1 - alpha*map2
       end where

       write(alpha_text,fmt='(f4.2)') alpha
       call write_minimal_header(header, 'MAP', nside=nside, order=1, polar=.true.)
       call write_result_map(trim(outprefix)//'_'//trim(alpha_text)//'.fits',   &
            & nside, 1, header, diff)

       call query_strip(nside, 0.5d0*pi-db, 0.5d0*pi+db, listpix, nlist)

       tot   = 0.d0
       numpt = 0.d0
       do j = 0, nlist-1
          if (diff(listpix(j),comp) /= 0.d0) then
             call pix2ang_ring(nside, listpix(j), theta, phi)
             b = int((phi+0.5d0*dl)/dl)
             if (b > 179) b = b-360
             tot(b)   = tot(b)   + diff(listpix(j),comp)
             numpt(b) = numpt(b) + 1
          end if
       end do

       open(59,file=trim(outprefix)//'_'//trim(alpha_text)//'.dat')
       do b = -180, 179
          if (numpt(b) > 0) then
             write(59,*) b, tot(b)/numpt(b)
          end if
       end do
       close(59)
       
    end do
    close(58)
    

  end subroutine print_scaled_gal_avg
  
  subroutine fit_ideal_dust
    implicit none

    logical(lgt) :: anynull
    real(dp)     :: nullval, T0, beta0, x
    integer(i4b) :: numband, nside, npix, ordering, nmaps, i, j, k, ierr
    real(dp), allocatable, dimension(:,:) :: maps, par
    real(dp), allocatable, dimension(:,:) :: cmb, fg
    character(len=512) :: infile1, infile2, outprefix, temp
    character(len=80), dimension(180)         :: header

    call getarg(2,infile1)
    call getarg(3,infile2)
    call getarg(4,outprefix)
    numband = iargc()-4
    T0      = 18.d0
    beta0   = 1.5d0

    ! Read input map
    i = getsize_fits(infile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(maps(0:npix-1,numband), par(0:npix-1,5), nu(numband), d(numband), cmb(0:npix-1,1), fg(0:npix-1,1))
    call read_bintab(infile1, fg, npix, nmaps, nullval, anynull)
    call read_bintab(infile2, cmb, npix, nmaps, nullval, anynull)

    x = h / (k_B*T0)
    do i = 1, numband
       call getarg(4+i,temp)
       read(temp,*) nu(i)
       nu(i) = nu(i) * 1d9
       maps(:,i) = fg(:,1) * &
            & (exp(x*nu(1))-1.d0) / (exp(x*nu(i))-1.d0) * (nu(i)/nu(1))**(beta0+1.d0) + & 
            & cmb(:,1) / ant2thermo(nu(i)/1d9)
    end do

    maps(:,1) = 1.00d0 * (maps(:,1) + 0.)
!    maps(:,4) = 1.01d0 * (maps(:,4) + 0.)

    ! Perform the fit
    do i = 0, npix-1
       if (mod(i,1000) == 0) write(*,*) i, npix
       d = maps(i,:)
       par(i,1) = fg(i,1)
       par(i,2) = beta0
       par(i,3) = T0
       par(i,4) = cmb(i,1)
       call powell(par(i,1:4), chisq_ideal_dust, ierr)       
       par(i,5) = chisq_ideal_dust(par(i,1:3))
    end do

    ! Output result files
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=.false.)
    call write_result_map(trim(outprefix)//'_amp.fits',   nside, ordering, header, par(:,1:1))
    call write_result_map(trim(outprefix)//'_beta.fits',  nside, ordering, header, par(:,2:2))
    call write_result_map(trim(outprefix)//'_T.fits',     nside, ordering, header, par(:,3:3))
    call write_result_map(trim(outprefix)//'_cmb.fits',     nside, ordering, header, par(:,4:4))
    call write_result_map(trim(outprefix)//'_chisq.fits', nside, ordering, header, par(:,5:5))

  end subroutine fit_ideal_dust

  function chisq_ideal_dust(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                           :: chisq_ideal_dust

    integer(i4b) :: i
    real(dp)     :: x, scale

    chisq_ideal_dust = 0.d0
    x                = h / (k_B*p(3))
    do i = 1, size(d)
       scale = 1.d0; if (i == 1) scale = 1.01d0
       chisq_ideal_dust = chisq_ideal_dust + &
            & (d(i) - &
            &   (scale*p(1)*(exp(x*nu(1))-1.d0) / (exp(x*nu(i))-1.d0) * (nu(i)/nu(1))**(p(2)+1.d0) + &
            &    p(4)/ ant2thermo(nu(i)/1d9)))**2
    end do

  end function chisq_ideal_dust



  subroutine compute_index_map
    implicit none

    character(len=256) :: mapname1, mapname2, outfile, tmp
    integer(i4b)       :: i, nside, nside2, ordering, ordering2, nmaps, nmaps2, npix, npix2
    real(dp)           :: nu1, nu2, nullval
    logical(lgt)       :: anynull
    real(dp), allocatable, dimension(:,:) :: map1, map2, ind
    character(len=80), dimension(180)         :: header

    if (iargc() < 6) then
       write(*,*) 'Usage: map_editor compute_spectral_index_map [map1] [nu1] [map2] [nu2] [index map] '
       stop
    end if

    call getarg(2,mapname1)
    call getarg(3,tmp)
    read(tmp,*) nu1
    call getarg(4,mapname2)
    call getarg(5,tmp)
    read(tmp,*) nu2
    call getarg(6,outfile)

    ! Read input map
    i = getsize_fits(mapname1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps), map2(0:npix-1,nmaps))
    call read_bintab(mapname1, map1, npix, nmaps, nullval, anynull)

    i = getsize_fits(mapname2, nside=nside2, ordering=ordering2, nmaps=nmaps2)
    npix2 = nside2npix(nside2)
    call read_bintab(mapname2, map2, npix2, nmaps2, nullval, anynull)
    if (nside /= nside2) then
       write(*,*) 'Error -- incompatible Nside. Exiting'
       stop
    end if
    write(*,*) ordering, ordering2
    if (ordering2 /= ordering) then
       write(*,*) 'a'
       do i = 1, nmaps2
          if (ordering2 == 1) then
             call convert_ring2nest(nside, map2(:,i))
          else
             call convert_nest2ring(nside, map2(:,i))
          end if
       end do
    end if

    nmaps = min(nmaps, nmaps2)
    allocate(ind(0:npix-1,nmaps))
    ind = -1.6375d30
    do i = 1, nmaps
       if (nu1 > 0) then
          where (map1(:,i)*map2(:,i) > 0.d0)
             ind(:,i) = log((map1(:,i)/ant2thermo(nu1)) / (map2(:,i)/ant2thermo(nu2))) / log(nu1/nu2)
          end where
       else
          where (map1(:,i)*map2(:,i) > 0.d0)
             ind(:,i) = log(map1(:,i) / map2(:,i)) / log(abs(nu1/nu2))
          end where
       end if
    end do

    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=(nmaps==3))
    call write_result_map(outfile, nside, ordering, header, ind)

  end subroutine compute_index_map

  subroutine make_ptsrc_map
    implicit none

    integer(i4b)       :: i, n, nside, npix, nmaps, unit, col, nlist, nsrc
    character(len=128) :: partext, convention, infile, outfile, beamfile
    character(len=1024) :: line
    real(dp)           :: vec(3), cvec(3), r, sigma, dOmega, Omega, fwhm, norm, unit_conv
    real(dp)           :: u(3), v(3), len_u, len_v, theta, flux, z, sgn
    integer(i4b), allocatable, dimension(:)   :: listpix
    real(dp),     allocatable, dimension(:,:) :: map
    real(dp),     allocatable, dimension(:)   :: data, src
    real(dp),     allocatable, dimension(:)   :: b_r, b_prof, b_prof2
    character(len=80), dimension(180)         :: header

    if (iargc() < 8) then
       write(*,*) 'Usage: map_editor make_ptsrc_map [ptsrc catalog] {WMAP,Planck} [flux column] '
       write(*,*) '            [nside] [fwhm] [Kcmb to MJy/sr] [outfile] (radial beam profile for WMAP)'
       stop
    end if
    
    call getarg(2,infile)
    call getarg(3,convention)
    call getarg(4,partext)
    read(partext,*) col
    call getarg(5,partext)
    read(partext,*) nside
    call getarg(6,partext)
    read(partext,*) fwhm
    call getarg(7,partext)
    read(partext,*) unit_conv
    call getarg(8,outfile)
    if (trim(convention) == 'WMAP') then
       call getarg(9,beamfile)
    end if
    npix   = 12*nside**2
    nmaps  = 3
    unit   = 58
    fwhm   = fwhm * pi/180.d0/60.
    sigma  = fwhm / sqrt(8.d0*log(2.d0))
    dOmega = 4.d0*pi / npix

    allocate(map(0:npix-1,nmaps), listpix(0:npix-1), src(0:npix-1))
    if (trim(convention) == 'WMAP') then
       allocate(data(col+2))
       open(unit,file=trim(beamfile))
       read(unit,*) n
       allocate(b_r(n), b_prof(n), b_prof2(n))
       do i = 1, n
          read(unit,*) b_r(i), b_prof(i)
       end do
       close(unit)
       b_r = b_r * pi/180.d0
!!$       open(58,file='test.dat')
!!$       do i = 1, n
!!$          write(58,*) b_r(i)*180.d0/pi, b_prof(i)
!!$       end do
       !call spline(b_r, b_prof, 0.d0, 0.d0, b_prof2)
       call smooth_spline('uniform', 1.d-12, b_r, b_prof, 0.d0, 0.d0, b_prof2)
!!$       write(58,*)
!!$       do i = 1, n
!!$          write(58,*) b_r(i)*180.d0/pi, b_prof(i)
!!$       end do
!!$       write(58,*)
!!$       do i = 0, 50000
!!$          write(58,*) i*1d-4, splint(b_r, b_prof, b_prof2, i*1.d-4*pi/180.d0)
!!$       end do
!!$       close(58)
!!$       stop
    else
       allocate(data(7))
    end if
    map = 0.d0
    open(unit,file=trim(infile))
    nsrc = 0
    do while (.true.)
       if (trim(convention) == 'WMAP') then
          read(unit,*,end=91) data
          flux = data(2+col)
          call ang2vec(0.5d0*pi-data(2)*pi/180.d0, data(1)*pi/180.d0, cvec)
          call query_disc(nside, cvec, 5*fwhm, listpix, nlist)       
          src = 0.d0
          do i = 0, nlist-1
             call pix2vec_ring(nside, listpix(i), vec)
             r      = acos(sum(vec*cvec))
             src(i) = max(splint(b_r, b_prof, b_prof2, r),0.d0)
          end do
          ! Listed values are peak flux density
          norm = max(data(2+col),0.d0) / (sum(src(0:nlist-1))*dOmega)
       else if (trim(convention) == 'Planck') then
          read(unit,*,end=91) data
          flux      = data(4) * 1.d-3 * 1.d-6
          data(1:2) = data(1:2) * pi/180.d0
          data(5:6) = data(5:6) * pi/180.d0/60.d0 / sqrt(8.d0*log(2.d0))
          data(7)   = data(7)   * pi/180.d0
          nsrc = nsrc+1
          if (mod(nsrc,100) == 0) write(*,*) 'Processing source no. ', nsrc
          call ang2vec(0.5d0*pi-data(2), data(1), cvec)
          call query_disc(nside, cvec, 10*fwhm, listpix, nlist)       
          src = 0.d0
          do i = 0, nlist-1
             call pix2vec_ring(nside, listpix(i), vec)
             r         = acos(sum(vec*cvec))
             if (data(5) /= data(6) .and. .false.) then
                z         = sum(cvec*vec)
                sgn    = 1.d0
                if (cvec(1)*vec(2)-cvec(2)*vec(1) < 0.d0) sgn = -1.d0
                u         = cvec(3) * cvec
                u(3)      = u(3) - 1.d0
                v         = vec - z * cvec
                len_u     = sqrt(sum(u*u))
                len_v     = sqrt(sum(v*v))
                theta     = sgn*acos(max(min((z * cvec(3) - vec(3)) / (len_u*len_v), 1.d0), -1.d0))
                theta     = theta + data(7)
             else
                theta = 0.d0
             end if
             src(i)    = exp(-0.5d0 * ((r*cos(theta)/data(5))**2 + (r*sin(theta)/data(6))**2))
          end do
       else
          write(*,*) 'Unsupported flux convention: ', trim(convention)
          stop
       end if
       ! Listed values are peak flux density
       norm           = max(flux,0.d0) / (sum(src(0:nlist-1))*dOmega)
       src(0:nlist-1) = norm / unit_conv * src(0:nlist-1)
       map(listpix(0:nlist-1),1) = map(listpix(0:nlist-1),1) + src(0:nlist-1)
    end do
91  close(unit)

    call write_minimal_header(header, 'MAP', nside=nside, order=1, polar=.true.)
    call write_result_map(outfile, nside, 1, header, map)
    deallocate(map)

  end subroutine make_ptsrc_map

  subroutine convert_beam(fwhm, lmax, nmaps, outfile, infile)
    implicit none

    integer(i4b),     intent(in) :: lmax, nmaps
    real(dp),         intent(in) :: fwhm
    character(len=*), intent(in) :: outfile
    character(len=*), intent(in), optional :: infile

    integer(i4b) :: nlheader, i, l
    real(dp),     allocatable, dimension(:)     :: b_in
    real(dp),     allocatable, dimension(:,:)   :: beam
    character(len=80), dimension(80) :: header

    allocate(beam(0:lmax, nmaps), b_in(nmaps))
    if (present(infile)) then
       if (infile(len(trim(infile))-3:len(trim(infile))) == 'fits') then
          call fits2cl(infile, beam, lmax, nmaps, header)
       else
          open(58,file=trim(outfile))
          beam = 0.d0
          do while (.true.)
             read(58,*,end=94) i, b_in
             beam(i,:) = b_in
          end do
94        close(58)          
       end if
    else
       call generate_beam(fwhm, lmax, beam)
    end if

    if (outfile(len(trim(outfile))-3:len(trim(outfile))) == 'fits') then
       header   = ''
       nlheader = 1
       call write_asctab(beam, lmax, nmaps, header, nlheader, outfile)
    else
       open(58,file=trim(outfile))
       do l = 0, lmax
          write(58,*) l, real(beam(l,:),sp)
       end do
       close(58)
    end if
    deallocate(beam)

  end subroutine convert_beam

  subroutine output_pointsource
    implicit none

    character(len=100) :: string_int
    integer(i4b)       :: nside, nmaps, npix, i, j, nlist
    real(dp)           :: lat, lon, fwhm_Q, ell_Q, ell_dir_Q, amp_Q, psi, beta, nu, nu0, aq, bq, cq
    real(dp)           :: fwhm_U, ell_U, ell_dir_U, au, bu, cu, amp_U, sigma_U
    real(dp)           :: vec(3), cvec(3), x, y, f, sigma_Q, sx_Q, sy_Q, sx_U, sy_U, cos2psi, sin2psi, ratio
    character(len=128) :: filename
    character(len=80), dimension(180) :: header

    real(dp),     allocatable, dimension(:,:) :: map
    integer(i4b), allocatable, dimension(:)   :: listpix

    call getarg(2,string_int)
    read(string_int,*) lon
    call getarg(3,string_int)
    read(string_int,*) lat
    call getarg(4,string_int)
    read(string_int,*) fwhm_Q
    call getarg(5,string_int)
    read(string_int,*) ell_Q     
    call getarg(6,string_int)
    read(string_int,*) ell_dir_Q
    call getarg(7,string_int)
    read(string_int,*) amp_Q
    call getarg(8,string_int)
    read(string_int,*) fwhm_U
    call getarg(9,string_int)
    read(string_int,*) ell_U     
    call getarg(10,string_int)
    read(string_int,*) ell_dir_U
    call getarg(11,string_int)
    read(string_int,*) amp_U
!    call getarg(12,string_int)
!    read(string_int,*) psi
!    call getarg(13,string_int)
!    read(string_int,*) beta
!    call getarg(14,string_int)
!    read(string_int,*) nu
!    call getarg(15,string_int)
!    read(string_int,*) nu0
    call getarg(12,string_int)
    read(string_int,*) nside
    call getarg(13,filename)

    lat       = lat * DEG2RAD
    lon       = lon * DEG2RAD
    ell_dir_Q = ell_dir_Q * DEG2RAD
    ell_dir_U = ell_dir_U * DEG2RAD
    psi       = psi * DEG2RAD
    sigma_Q   = fwhm_Q/sqrt(8.d0*log(2.d0))
    sigma_U   = fwhm_U/sqrt(8.d0*log(2.d0))
    sx_Q      = 2*sigma_Q / (2.d0-ell_Q)
    sx_U      = 2*sigma_U / (2.d0-ell_U)
    sy_Q      = 2*sigma_Q*(1.d0-ell_Q) / (2.d0-ell_Q)
    sy_U      = 2*sigma_U*(1.d0-ell_U) / (2.d0-ell_U)

    ! Allocate data structures
    npix  = 12*nside**2
    nmaps = 3
    allocate(map(0:npix-1,nmaps), listpix(0:npix-1))
    call ang2vec(0.5d0*pi-lat, lon, cvec)
    call query_disc(nside, cvec, 5*max(sigma_Q,sigma_U)*DEG2RAD, listpix, nlist, nest=1)

!    f     = (nu/nu0)**beta
!    ratio = (2.d0*pi*(0.88/sqrt(8*log(2.d0))*DEG2RAD)**2) / (2.d0*pi*(sigma*DEG2RAD)**2)
    
    map = 0.d0
    open(58,file='dist.dat')
    do i = 0, nlist-1
       call pix2vec_nest(nside, listpix(i), vec)
       call project2xy(vec, cvec, x, y)

       aq  = cos(ell_dir_Q)**2/(2*sx_Q**2) + sin(ell_dir_Q)**2/(2*sy_Q**2)
       bq  = - sin(2*ell_dir_Q)/(4*sx_Q**2) + sin(2*ell_dir_Q)/(4*sy_Q**2)
       cq  = sin(ell_dir_Q)**2/(2*sx_Q**2) + cos(ell_dir_Q)**2/(2*sy_Q**2)
       au  = cos(ell_dir_U)**2/(2*sx_U**2) + sin(ell_dir_U)**2/(2*sy_U**2)
       bu  = - sin(2*ell_dir_U)/(4*sx_U**2) + sin(2*ell_dir_U)/(4*sy_U**2)
       cu  = sin(ell_dir_U)**2/(2*sx_U**2) + cos(ell_dir_U)**2/(2*sy_U**2)

       map(listpix(i),2) = amp_Q * exp(-(aq*x**2 + 2*bq*x*y + cq*y**2))
       map(listpix(i),3) = amp_U * exp(-(au*x**2 + 2*bu*x*y + cu*y**2))

!       map(listpix(i),1) = amp_T * ratio * f * exp(-(a*x**2 + 2*b*x*y + c*y**2))
!       map(listpix(i),2) = amp_P * ratio * f * exp(-(a*x**2 + 2*b*x*y + c*y**2)) * cos(2.d0*psi)
!       map(listpix(i),3) = amp_P * ratio * f * exp(-(a*x**2 + 2*b*x*y + c*y**2)) * sin(2.d0*psi)

       !map(listpix(i),1) = amp_T * f * exp(-0.5d0 * (x**2 + y**2)/sigma**2)
       !map(listpix(i),2) = amp_P * f * exp(-0.5d0 * (x**2 + y**2)/sigma**2) * cos(2.d0*psi)
       !map(listpix(i),3) = amp_P * f * exp(-0.5d0 * (x**2 + y**2)/sigma**2) * sin(2.d0*psi)
       write(58,*) sqrt(x**2+y**2), map(listpix(i),2)
    end do
    close(58)
!    write(*,*) 'amp_T = ', amp_T
!    write(*,*) 'amp_P = ', amp_P
!    write(*,*) 'f     = ', f
!    write(*,*) maxval(map(:,1)), maxval(map(:,2:3))

    call write_minimal_header(header, 'MAP', nside=nside, ordering='NESTED', polar=nmaps==3)
    call write_result_map(filename, nside, 2, header, map)
    deallocate(map, listpix)

  end subroutine output_pointsource



  subroutine project2xy(mpixvec, centervec, x, y)
    implicit none

    real(dp)               :: sigma, radius, pixdist, northpixang, psi_true, x, y
    real(dp), dimension(3) :: centervec,  northvec, centerpixvec, centernorthvec, somevec, mpixvec, testvec
   
    psi_true = 0.d0
    northvec=[0,0,1]
    somevec=[centervec(2),-centervec(1),0.d0] 
    call crossproduct(somevec, centervec, centernorthvec)
    call crossproduct(centervec, mpixvec, somevec)
    call crossproduct(somevec, centervec, centerpixvec)
    call angdist(centernorthvec, centerpixvec, northpixang)
    call angdist(centervec, mpixvec, pixdist) 
    call crossproduct(centernorthvec, centerpixvec, testvec)
    if (sum(testvec*centervec) > 0) northpixang=-northpixang
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

  ! Returns cos(2*psi) and sin(2*psi) of the rotation induced by
  ! parallel transport from vec1 to vec2.
  subroutine qu_transport_rot_v1(vec1, vec2, cos2psi, sin2psi)
    implicit none

    real(dp), dimension(3), intent(in)  :: vec1, vec2
    real(dp),               intent(out) :: cos2psi, sin2psi

    integer(i4b) :: i, j, nfield
    real(dp) :: len_u, len_v, z, cos_theta, sgn, c1, s1, c2, s2
    real(dp), dimension(3) :: u, v

    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       cos2psi = 1; sin2psi = 0
       return
    end if

    sgn    = 1.d0
    if (vec1(1)*vec2(2)-vec1(2)*vec2(1) < 0.d0) sgn = -1.d0

    ! Rotation from vec1 to vec 2
    u         = vec1(3) * vec1 
    u(3)      = u(3) - 1.d0
    v         = vec2 - z * vec1
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    ! Local angle from vertical
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    ! Double it and calculate cos and sin
    c1 = 2*cos_theta**2-1
    s1 = 2*sgn*sqrt(1-cos_theta**2)*cos_theta

    ! Rotation from vec2 to vec 1; sgn is opposite from 1->2
    u          = vec2(3) * vec2
    u(3)       = u(3) - 1.d0
    v          = vec1 - z * vec2
    len_u      = sqrt(sum(u*u))
    len_v      = sqrt(sum(v*v))
    cos_theta  = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    c2 =  2*cos_theta**2-1
    s2 = -2*sgn*sqrt(1-cos_theta**2)*cos_theta

    cos2psi = c1*c2+s1*s2
    sin2psi = c1*s2-s1*c2
  end subroutine


  subroutine shift_columns(filename, ncol, nside, ordering, header, map)
    implicit none

    character(len=*),                         intent(in)  :: filename
    integer(i4b),                             intent(in)  :: ncol
    integer(i4b),                             intent(out) :: nside, ordering
    real(dp),         pointer, dimension(:,:)             :: map
    character(len=80),         dimension(180)             :: header

    integer(i4b) :: j, nmaps, npix
    real(dp)     :: nullval
    logical(lgt) :: anynull

    j    = getsize_fits(filename, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = 12*nside**2

    allocate(map(0:npix-1,nmaps+ncol))
    map = 0.d0
    call read_bintab(filename, map(:,ncol+1:ncol+nmaps), npix, nmaps, nullval, anynull, header=header)

  end subroutine shift_columns

  subroutine compute_weighted_sum(infoname, nside, ordering, nmaps, outmap, header)
    implicit none

    character(len=*),                         intent(in)  :: infoname
    integer(i4b),                             intent(out) :: nside, ordering, nmaps
    real(dp),         pointer, dimension(:,:)             :: outmap
    character(len=80),         dimension(180), intent(out) :: header

    integer(i4b)       :: i, j, unit, npix, nummaps, nside_in, ordering_in, nmaps_in
    real(dp)           :: temp_r, weight, nullval, missval = -1.6375e30
    logical(lgt)       :: anynull
    character(len=128) :: filename

    real(dp), allocatable, dimension(:,:)   :: map

    unit = 24

    open(unit, file=trim(infoname))
    read(unit,*) nside
    read(unit,*) nummaps

    if (nummaps < 1) then
       write(*,*) 'Error: The number of input maps must be larger than 0.'
       stop
    end if

    ! Get parameters from the first map file, and check consistency
    read(unit,*) filename
    j = getsize_fits(filename, nside=nside, ordering=ordering, nmaps=nmaps)
    do i = 2, nummaps
       read(unit,*) filename
       j = getsize_fits(filename, nside=nside_in, ordering=ordering_in, nmaps=nmaps_in)
       if (nside_in /= nside) then
          write(*,*) 'Error: Nsides of input maps are not consistent'
          stop
       end if
       if (nmaps_in /= nmaps) then
          write(*,*) 'Error: Nmaps of input maps are not consistent'
          nmaps = 1
!          stop
       end if
    end do

    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    allocate(outmap(0:npix-1,nmaps))

    close(unit)
    open(unit, file=trim(infoname))
    read(unit,*) nside
    read(unit,*) nummaps

    outmap = 0.
    do i = 1, nummaps
       read(unit,*) filename, weight

       j = getsize_fits(filename, ordering=ordering_in)
       if (i == 1) then
          call read_bintab(filename, map, npix, nmaps, nullval, anynull, header=header)
       else
          call read_bintab(filename, map, npix, nmaps, nullval, anynull)
       end if

       if (ordering_in /= ordering) then
          if (ordering == 1) then
             do j = 1, nmaps
                call convert_nest2ring(nside, map(:,j))
             end do
          else
             do j = 1, nmaps
                call convert_ring2nest(nside, map(:,j))
             end do
          end if
       end if

       where (outmap /= missval .and. map /= missval) 
          outmap = outmap + weight * map
       elsewhere
          outmap = missval
       end where

    end do
    close(unit)

    deallocate(map)

  end subroutine compute_weighted_sum

  subroutine maskcount(maskfile)
    implicit none

    character(len=*),                   intent(in)    :: maskfile

    integer(i4b) :: npix, nside, ordering, nmaps, i
    real(dp)     :: nullval
    logical(lgt) :: anynull
    real(dp), allocatable, dimension(:,:) :: mask

    npix = getsize_fits(maskfile, nside=nside, ordering=ordering, nmaps=nmaps)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    do i=1,nmaps
       write(*,*) int(i,i2b), 'Unmasked =', count(mask(:,i) > 0.5d0)/real(npix,dp)*100.d0, &
            & 'Masked =', count(mask(:,i)<=0.5d0)/real(npix,dp)*100.d0
    end do
    write(*,*) 'Total unmasked =', count(mask>0.5d0)/real(size(mask),dp)*100.d0, &
         & 'Masked =', count(mask<=0.5d0)/real(size(mask),dp)*100.d0
    write(*,*)
    do i=1,nmaps
       write(*,*) int(i,i2b), 'Number unmasked =', count(mask(:,i)>0.5d0), &
            & ' Number masked =', count(mask(:,i)<=0.5d0)
    end do
    write(*,*) int(i,i2b), 'Total number unmasked =', count(mask>0.5d0), &
         & ', Total number masked =', count(mask<=0.5d0)

  end subroutine maskcount

  subroutine badcount(mapfile)
    implicit none

    character(len=*),                   intent(in)    :: mapfile

    integer(i4b) :: npix, nside, ordering, nmaps, i
    real(dp)     :: nullval, missval = -1.6375e30
    logical(lgt) :: anynull
    real(dp), allocatable, dimension(:,:) :: mask

    npix = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(mapfile, mask, npix, nmaps, nullval, anynull)
    write(*,*)
    do i=1,nmaps
       write(*,*) int(i,i2b), 'Bad pixels =', count(mask(:,i)==missval), ' Good pixels =', count(mask(:,i)/=missval)
    end do
    write(*,*) int(i,i2b), 'Total number bad =', count(mask==missval), ', Total number good =', count(mask/=missval)

  end subroutine badcount

  subroutine apply_mask(maskfile, nside, ordering, map, fact)
    implicit none

    integer(i4b),                       intent(in)    :: nside, ordering
    character(len=*),                   intent(in)    :: maskfile
    real(dp),         dimension(0:,1:), intent(inout) :: map
    real(dp),                           intent(in), optional :: fact

    integer(i4b) :: i, j, nside_mask, ordering_mask, npix, n, nmaps, nmaps_mask
    real(dp)     :: nullval
    logical(lgt) :: anynull
    real(dp)     :: vec(3), radius
    integer(i4b), allocatable, dimension(:) :: listpix
    real(dp), allocatable, dimension(:,:) :: mask

    nmaps = size(map,2)

    ! Read rms file and check consistency
    npix = getsize_fits(maskfile, nside=nside_mask, ordering=ordering_mask, nmaps=nmaps_mask)

    if (nside_mask /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if
    if (nmaps_mask /= nmaps) then
       write(*,*) 'Error: Different Nmaps of the two files. Putting nmaps=1'
       write(*,*) 'Mask will be applied to all layers of map.'
       nmaps=1
    end if

    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    
    do i = 1, nmaps
       if (ordering_mask /= ordering) then
          if (ordering_mask == 1) call convert_ring2nest(nside, mask(:,i))
          if (ordering_mask == 2) call convert_nest2ring(nside, mask(:,i))
       end if
    end do

    if (present(fact)) then
       allocate(listpix(0:npix-1))
       radius = fact * pi/180.d0/60.
       do i = 0, npix-1
          if (mod(i,1000) == 0) write(*,*) i, npix
          if (ordering == 1) then
             call pix2vec_ring(nside, i, vec)
          else 
             call pix2vec_nest(nside, i, vec)            
          end if
          call query_disc(nside, vec, radius, listpix, n, nest=ordering-1) !query_disc uses ring by default
          do j = 1,nmaps
             if (any(mask(listpix(0:n-1),j) == 1.d0) .and. any(mask(listpix(0:n-1),j) == 0.d0)) then
                ! Border point
                if (nmaps > 1) then
                   map(listpix(0:n-1),j) = -1.6375d30
                else
                   map(listpix(0:n-1),:) = -1.6375d30
                end if
             end if
          end do
       end do
       deallocate(listpix)
    else
       do j = 1, nmaps
          do i = 0, npix-1
             if (mask(i,j) < 0.5d0) then
                if (nmaps > 1) then
                   map(i,j) = -1.6375e30
                else
                   map(i,:) = -1.6375e30
                end if
             end if
          end do
       end do
    end if

    ! TMR hack to mask gc
!    call ang2vec(pi/2.,0.d0,vec)
!    write(*,*) vec
!    allocate(listpix(0:npix-1))
!    call query_disc(nside, vec, 0.8d0*DEG2RAD, listpix, n, nest=1)
!    write(*,*) listpix(0:n-1)
!    write(*,*) 'Masking out gc pixels'
!    map(listpix(0:n-1),:) = -1.6375e30

  end subroutine apply_mask


  subroutine output_mean_and_stddev(nside, ordering, nmaps, header, mu, rms)
    implicit none

    integer(i4b),                               intent(out)          :: nside, nmaps, ordering
    real(dp),           pointer, dimension(:,:)                      :: mu, rms
    character(len=80),           dimension(180)                      :: header

    integer(i4b)       :: npix
    integer(i4b)       :: i, j, l, m, numfiles
    real(dp)           :: nullval
    real(dp)           :: sigma_sq, nullval_dp
    logical(lgt)       :: anynull
    character(len=128) :: winfile_in, winfile_out
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    real(dp),     allocatable, dimension(:,:)   :: map
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam
    character(len=256), allocatable, dimension(:) :: filenames
    real(dp),                  dimension(2)     :: zbounds = 0.d0

    numfiles = iargc()-3
    allocate(filenames(numfiles))
    do i = 1, numfiles
       call getarg(i+3, filenames(i))
    end do

    ! Read input map
    i = getsize_fits(filenames(1), nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), mu(0:npix-1,nmaps), rms(0:npix-1,nmaps))
    
    ! Compute mu
    mu = 0.d0
    do i = 1, numfiles
       if (i == 1) then
          call read_bintab(filenames(i), map, npix, nmaps, nullval, anynull, header=header)
       else
          call read_bintab(filenames(i), map, npix, nmaps, nullval, anynull)          
       end if

       !write(*,*) i, mean(map(:,1)), sqrt(variance(map(:,1)))
       write(*,*) 'a', i, numfiles
       !call input_map(filenames(i), map, npix, nmaps)
       mu = mu + map
    end do
    mu = mu / numfiles

    ! Compute standard deviation
    rms = 0.d0
    do i = 1, numfiles
       call read_bintab(filenames(i), map, npix, nmaps, nullval, anynull)
       write(*,*) 'b', i, numfiles
!       call input_map(filenames(i), map, npix, nmaps)
       rms = rms + (map-mu)**2
    end do
    rms = sqrt(rms / (numfiles-1))

    deallocate(map)

  end subroutine output_mean_and_stddev

  subroutine smooth_map(infile, r_fill, lmin, lmax, nside, ordering, nmaps, map, header, &
       & output_EB, beamfile_in, beamfile_out, fwhm_in, fwhm_out)
    implicit none

    character(len=128),                         intent(in)           :: infile
    real(dp),                                   intent(inout)        :: r_fill
    integer(i4b),                               intent(in)           :: nside, lmin, lmax
    real(dp),                                   intent(in), optional :: fwhm_in, fwhm_out
    character(len=128),                         intent(in), optional :: beamfile_in, beamfile_out
    integer(i4b),                               intent(out)          :: nmaps, ordering
    logical(lgt),                               intent(in)           :: output_EB
    real(dp),           pointer, dimension(:,:)                      :: map
    character(len=80),           dimension(180)                      :: header

    integer(i4b)       :: nside_in, npix, npix_in
    integer(i4b)       :: i, j, k, l, m, r, nlist
    real(dp)           :: nullval, tot
    real(dp)           :: sigma_sq, nullval_dp, vec(3), theta, phi, rf
    logical(lgt)       :: anynull
    character(len=128) :: winfile_in, winfile_out
    character(len=4)   :: ntext

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map_in, map_buffer
    real(dp),     pointer,     dimension(:,:)   :: pixwin_in, pixwin_out
    real(dp),     pointer,     dimension(:,:)   :: weights
    real(dp),     allocatable, dimension(:,:)   :: beam_in, beam_out
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    integer(i4b), allocatable, dimension(:)     :: listpix

    ! Read input map
    i = getsize_fits(infile, nside=nside_in, ordering=ordering, nmaps=nmaps)
    if (nmaps /= 1 .and. nmaps /= 3) then
       if (nmaps < 3) then
          nmaps = 1
       else
          nmaps = 3
       end if
    end if

    npix    = nside2npix(nside)
    npix_in = nside2npix(nside_in)
    allocate(map_in(0:npix_in-1,nmaps), map_buffer(0:npix_in-1,nmaps))
    allocate(map(0:npix-1,nmaps))
    call read_bintab(infile, map_in, npix_in, nmaps, nullval, anynull, header=header)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside_in, map_in(0:npix_in-1,i))
       end do
    end if

    ! Read pixel windows
    call read_pixwin(nside_in, nmaps, pixwin_in)
    call read_pixwin(nside, nmaps, pixwin_out)
    call read_ringweights(nside_in, nmaps, weights)

    ! Create or read beams
    allocate(beam_in(0:lmax, nmaps))
    allocate(beam_out(0:lmax, nmaps))

    if (present(fwhm_in)) then
       call generate_beam(abs(fwhm_in), lmax, beam_in)
       if (fwhm_in < 0.) beam_in = 1.d0 / beam_in
    else
       call generate_beam(0.d0, lmax, beam_in, beamfile_in)
       !if (nmaps == 3) then
       !   beam_in(:,2) = beam_in(:,1)
       !   beam_in(:,3) = beam_in(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    if (present(fwhm_out)) then
       call generate_beam(abs(fwhm_out), lmax, beam_out)
    else
       call generate_beam(0.d0, lmax, beam_out, beamfile_out)
       !if (nmaps == 3) then
       !   beam_out(:,2) = beam_out(:,1)
       !   beam_out(:,3) = beam_out(:,1)
       !   write(*,*) 'Warning: Setting polarized beams equal to temperature beam'
       !end if
    end if

    allocate(listpix(0:npix_in-1))
    map_buffer = map_in
    do j = 1, nmaps
       do i = 0, npix_in-1
          if (abs(map_in(i,j)) > 1e30) then
             if (r_fill >= 0.d0) then
                call pix2vec_ring(nside_in, i, vec)
                do r = 1, r_fill
                   rf = (r+0.2d0) * sqrt(4.d0*pi/npix_in)
                   call query_disc(nside_in, vec, rf, listpix, nlist)
                   tot = 0.d0
                   m   = 0.d0
                   do k = 0, nlist-1
                      if (abs(map_buffer(listpix(k),j)) < 1.d30) then
                         tot = tot + map_buffer(listpix(k),j)
                         m   = m   + 1
                      end if
                   end do
                   if (m > 0) then
                      map_in(i,j) = tot / m
                      exit
                   end if
                end do

                if (r > r_fill) then
                   call pix2ang_ring(nside_in, i, theta, phi)
                   write(*,*) 'Error: Filter did not return valid pixel; increase r_fill'
                   write(*,*) '     p         = ', i
                   write(*,*) '     theta     = ', 90.d0-theta*180/pi
                   write(*,*) '     phi       = ', phi*180/pi
                   write(*,*) '     radius    = ', r_fill
                   write(*,*) '     neighbors = ', listpix(0:nlist-1)
                   write(*,*) '     vals      = ', map_buffer(listpix(0:nlist-1),j)
                   stop
                end if
             else
                map_in(i,j) = 0.d0
             end if
          end if
       end do
    end do
    deallocate(listpix, map_buffer)

    ! Compute the spherical harmonics transform
    allocate(alms(nmaps, 0:lmax, 0:lmax))

    alms = cmplx(0.,0.)
    if (nmaps == 1) then
       call map2alm(nside_in, lmax, lmax, map_in(:,1), alms, zbounds, weights)
    else
       call map2alm(nside_in, lmax, lmax, map_in, alms, zbounds, weights)
    end if

    ! Deconvolve old beam and pixel window, and convolve new beam and pixel window
    do i = 1, nmaps
       alms(i,0:lmin-1,:) = cmplx(0.,0.)
       do l = lmin, lmax
          if (abs(pixwin_in(l,i)*beam_in(l,i)) < 0.00000001d0) then
             if (l > 1 .or. i == 1) write(*,*) 'Multipole l = ', l, ' set to zero'
             alms(i,l,0:l) = cmplx(0.,0.)
          else
             alms(i,l,0:l) = alms(i,l,0:l) * &
                  & (pixwin_out(l,i)*beam_out(l,i)) / (pixwin_in(l,i)*beam_in(l,i))
          end if
       end do
    end do

    ! Compute the inverse spherical harmonics
    if (nmaps == 1) then
       call alm2map(nside, lmax, lmax, alms, map(:,1))
    else
       if (output_EB) then
          do i = 1, nmaps
             call alm2map(nside, lmax, lmax, alms(i:i,:,:), map(:,i))
          end do
       else
          call alm2map(nside, lmax, lmax, alms, map)       
       end if
    end if

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, map(:,i))
       end do
    end if

    deallocate(weights)
    deallocate(pixwin_in)
    deallocate(pixwin_out)
    deallocate(beam_in)
    deallocate(beam_out)
    deallocate(map_in)
    deallocate(alms)

  end subroutine smooth_map

  subroutine print_maximum(mapfile, component, fwhm)
    implicit none

    character(len=*),                          intent(in)           :: mapfile
    integer(i4b),                              intent(in)           :: component
    real(dp),                                  intent(in)           :: fwhm

    integer(i4b) :: npix, pix, nside, nmaps, ordering, i, nlist, nest
    real(dp)     :: nullval, max_val
    real(dp)     :: theta, phi, t1, t2, integral, radius
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)   :: map
    real(dp),                 dimension(3)     :: vector0
    integer(i4b), allocatable, dimension(:) :: listpix

    call cpu_time(t1)

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull)

    if (nmaps < component) then
       write(*,*) 'Requested component is not present in map file'
       stop
    end if

    pix    = 0
    max_val = -1.d30
    do i = 0, npix-1
       if (map(i,component) /= -1.6375e30) then
          if (map(i,component) > max_val) then
             pix    = i
             max_val = map(i,component)
          end if
       end if
    end do

    goto 1010

    ! Search for best-fit pixel
    radius = fwhm * pi/180.d0/60 * 4.d0
    allocate(listpix(0:npix-1))
    if (ordering == 1) then
       call pix2vec_ring(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist)
    else
       call pix2vec_nest(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist, nest=1)
    end if

    integral = convolve_with_gaussian(nside, ordering, map(:,component), pix, fwhm)

    pix     = 0
    max_val = -1.d30
    do i = 0, nlist-1
       integral = convolve_with_gaussian(nside, ordering, map(:,component), listpix(i), fwhm)
       if (integral > max_val) then
          pix     = listpix(i)
          max_val = integral
!          write(*,*) pix, real(map(pix,component),sp), max_val
       end if
    end do

    call cpu_time(t2)
!    write(*,*) t2-t1

1010 if (ordering == 1) then
       call pix2ang_ring(nside, pix, theta, phi)
    else
       call pix2ang_nest(nside, pix, theta, phi)
    end if

    write(*,*) pix, phi*180./pi, 90.-180./pi*theta

    deallocate(map)

  end subroutine print_maximum

  subroutine print_stats(mapfile, maskfile)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile

    integer(i4b) :: npix, pix, nside, nmaps, ordering, i, j, k, nlist, nest, order_mask
    real(dp)     :: nullval, max_val, med
    real(dp)     :: theta, phi, t1, t2, integral, radius, mu, sigma, my_min, my_max, scale
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    character(len=128) :: par
    real(dp),    allocatable, dimension(:,:)   :: map, mask
    real(dp),    allocatable, dimension(:)     :: vals

    call cpu_time(t1)

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=order_mask)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (order_mask /= ordering) then
       if (order_mask == 1) then
          do i = 1, nmaps
             call convert_ring2nest(nside, mask(:,i))
          end do
       else
          do i = 1, nmaps
             call convert_nest2ring(nside, mask(:,i))
          end do
       end if
    end if

    where (map < -1.637d30)
       mask = 0.d0
    end where

    allocate(vals(npix))
    do i = 1, nmaps
    !do i = 2, 2
       
       k = 0
       do j = 0, npix-1
          if (mask(j,i) > 0.5d0) then
             k = k+1
             vals(k) = map(j,i)
          end if
       end do

       mu    = sum(vals(1:k))/k
       sigma = sqrt(sum((vals(1:k)-mu)**2)/(k-1))
       if (iargc()==4) then
          med   = median(abs(vals(1:k)))
          call getarg(4, par)
          read(par,*) scale
          if (scale > 0.d0) then
             write(*,fmt='(f6.2)') sigma*scale
          else if (scale < 0) then
             write(*,fmt='(f6.3)') med*abs(scale)
          else if (scale == 0.d0) then
             write(*,fmt='(a,f8.3,a,f8.3,a)') '$', mu, "\pm", sigma, '$'
          end if
       else
          write(*,*) 
          write(*,*) '   Pol  = ', i
          write(*,*) '   Mean = ', mu
          write(*,*) '   RMS  = ', sigma
          write(*,*) '   Min  = ', minval(map(:,i), mask(:,i) > 0.5d0)
          write(*,*) '   Max  = ', maxval(map(:,i), mask(:,i) > 0.5d0)
       end if
    end do
    
    deallocate(map, mask)

  end subroutine print_stats

  subroutine print_stats_col(mapfile, maskfile)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile

    integer(i4b) :: npix, pix, nside, nmaps, ordering, i, j, k, nlist, nest, order_mask, col, extno
    real(dp)     :: nullval, max_val, med
    real(dp)     :: theta, phi, t1, t2, integral, radius, mu, sigma, my_min, my_max, scale
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    character(len=128) :: par
    real(dp),    allocatable, dimension(:,:)   :: map, mask
    real(dp),    allocatable, dimension(:)     :: vals

    call cpu_time(t1)

    call getarg(5,par)
    read(par,*) extno

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps, extno=extno)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, extno=extno)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=order_mask)
    allocate(mask(0:npix-1,1))
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)
    if (order_mask /= ordering) then
       if (order_mask == 1) then
          call convert_ring2nest(nside, mask(:,1))
       else
          call convert_nest2ring(nside, mask(:,1))
       end if
    end if

    where (map < -1.637d30)
       mask = 0.d0
    end where

    call getarg(4,par)
    read(par,*) col
    call getarg(6,par)
    read(par,*) scale

    allocate(vals(npix))
    k = 0
    do j = 0, npix-1
       if (mask(j,1) > 0.5d0) then
          k = k+1
          if (trim(mapfile) == 'data/products/COM_CompMap_DustPol-commander_1024_R2.00.fits' .or. &
               & trim(mapfile) == 'data/products/COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits') then
             vals(k) = sqrt(map(j,col)**2 + map(j,col+1)**2)
          else
             vals(k) = map(j,col)
          end if
       end if
    end do
    
    mu    = sum(vals(1:k))/k
    sigma = sqrt(sum((vals(1:k)-mu)**2)/(k-1))
    write(*,fmt='(a,f8.2,a,f8.2,a)') '$', mu*scale, "\pm", sigma*scale, '$'

    deallocate(map, mask)

  end subroutine print_stats_col

  function convolve_with_gaussian(nside, ordering, map, pix, fwhm)
    implicit none

    real(dp),     dimension(0:), intent(in) :: map
    integer(i4b),                intent(in) :: pix, nside, ordering
    real(dp),                    intent(in) :: fwhm
    real(dp)                                :: convolve_with_gaussian
    
    real(dp) :: s, radius, sigma
    integer(i4b) :: i, nlist, npix, nest, n
    integer(i4b), allocatable, dimension(:) :: listpix
    real(dp), dimension(3) :: vector0, vector

    npix   = size(map)
    sigma  = fwhm * pi/180.d0/60.d0 / sqrt(8.d0*log(2.d0))
    radius = 2.d0 * fwhm * pi/180.d0/60.d0
    nest   = ordering-1

    allocate(listpix(0:npix-1))
    
    if (ordering == 1) then
       call pix2vec_ring(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist)
    else
       call pix2vec_nest(nside, pix, vector0)
       call query_disc(nside, vector0, radius, listpix, nlist,nest=1)
    end if
    

    s = 0.d0
    n = 0
    do i = 0, nlist-1
       if (map(listpix(i)) /= -1.6375e30) then
          if (ordering == 1) then
             call pix2vec_ring(nside, listpix(i), vector)
          else
             call pix2vec_nest(nside, listpix(i), vector)
          end if
          if (listpix(i) /= pix) then
             radius = acos(sum(vector * vector0))
          else
             radius = 0.d0
          end if
          s = s + map(listpix(i)) * exp(-0.5d0*(radius/sigma)**2)
          n = n + 1
       end if
    end do

    if (n == nlist) then
       convolve_with_gaussian = s / real(n,dp)
    else
       convolve_with_gaussian = -1.d30
    end if

    deallocate(listpix)

  end function convolve_with_gaussian

  subroutine add_gaussian_noise(mapfile, rmsfile, seed, nside, ordering, nmaps, map, header, sigma_0)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, rmsfile
    integer(i4b),                              intent(inout)        :: seed
    integer(i4b),                              intent(out)          :: nside, ordering, nmaps
    real(dp),         pointer, dimension(:,:)                       :: map   
    character(len=80),         dimension(180)                       :: header
    real(dp),                                  intent(in), optional :: sigma_0

    integer(i4b) :: npix, nside_rms, ordering_rms, nmaps_rms, i, j
    real(dp)     :: nullval
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)   :: rmsmap

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)

    ! Read rms file and check consistency
    i = getsize_fits(rmsfile, nside=nside_rms, ordering=ordering_rms, nmaps=nmaps_rms)

    if (nside_rms /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    if (nmaps_rms /= nmaps) then
       write(*,*) 'Error: Different Nmaps of the two files.'
       stop
    end if

    allocate(rmsmap(0:npix-1,nmaps))
    call read_bintab(rmsfile, rmsmap, npix, nmaps, nullval, anynull)

    if (ordering_rms /= ordering) then
       if (ordering_rms == 1) then
          do i = 1, nmaps
             call convert_ring2nest(nside, rmsmap(:,i))
          end do
       else
          do i = 1, nmaps
             call convert_nest2ring(nside, rmsmap(:,i))
          end do
       end if
    end if

    call rand_init(rng_handle, seed)

    ! Add Gaussian noise to the input map
    if (present(sigma_0)) then
       ! Interpret RMS map as Nobs, not RMS
       do j = 1, nmaps
          do i = 0, npix-1
             if (map(i,j) /= -1.6375e30 .and. rmsmap(i,j) >= 0.) then
                map(i,j) = map(i,j) + sigma_0 / sqrt(rmsmap(i,j)) * rand_gauss(rng_handle)
             else 
                map(i,j) = -1.6375e30
             end if
          end do
       end do
    else
       do j = 1, nmaps
          do i = 0, npix-1
             if (map(i,j) /= -1.6375e30 .and. rmsmap(i,j) >= 0.) then
                map(i,j) = map(i,j) + rmsmap(i,j) * rand_gauss(rng_handle)
             else 
                map(i,j) = -1.6375e30
             end if
          end do
       end do
    end if

    deallocate(rmsmap)

  end subroutine add_gaussian_noise


  subroutine add_gaussian_noise_sqrt_N(mapfile, rmsfile, map2mask_file, seed, nside, ordering, nmaps, map, header, sigma_0)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, rmsfile, map2mask_file
    integer(i4b),                              intent(inout)        :: seed
    integer(i4b),                              intent(out)          :: nside, ordering, nmaps
    real(dp),         pointer, dimension(:,:)                       :: map   
    real(dp),                                  intent(in), optional :: sigma_0
    character(len=80),         dimension(180)                       :: header

    integer(i4b) :: npix, nside_rms, ordering_rms, nmaps_rms, i, j, k, polarization
    real(dp)     :: nullval
    integer(i8b) :: n
    logical(lgt) :: anynull, inv
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)   :: rmsmap, map2mask

    real(dp), allocatable, dimension(:,:) :: sqrt_N
    real(dp), allocatable, dimension(:)   :: eta

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)

    ! Read map2mask file
    i = getsize_fits(map2mask_file, nside=nside_rms, ordering=ordering_rms, nmaps=nmaps_rms)

    if (nside_rms /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    if (nmaps_rms /= nmaps) then
       write(*,*) 'Error: Different Nmaps of the two files.'
       stop
    end if

    allocate(map2mask(0:npix-1,nmaps))
    call read_bintab(map2mask_file, map2mask, npix, nmaps, nullval, anynull)

    ! Read covariance file
    call read_covmatrix(58, rmsfile, ordering_rms, polarization, sqrt_N, inv, n)

    if (ordering_rms /= ordering) then
       if (ordering_rms == 1) then
          do i = 1, nmaps
             call convert_ring2nest(nside, rmsmap(:,i))
          end do
       else
          do i = 1, nmaps
             call convert_nest2ring(nside, rmsmap(:,i))
          end do
       end if
    end if

    call rand_init(rng_handle, seed)

    ! Draw a random gaussian vector
    allocate(eta(n))
    do i = 1, n
       eta(i) = rand_gauss(rng_handle)
    end do
    eta = matmul(sqrt_N, eta)

    ! Add noise to the input map
    k = 1
    do j = 1, nmaps
       do i = 0, npix-1
          if (nint(map2mask(i,j)) > -1) then
             map(i,j) = map(i,j) + eta(k)
             k        = k+1
          else
             map(i,j) = -1.6375d30
          end if
       end do
    end do

    deallocate(map2mask)
    deallocate(eta)
    deallocate(sqrt_N)

  end subroutine add_gaussian_noise_sqrt_N


  subroutine subtract_mono_dipole(mapfile, maskfile, nside, ordering, map, header, md)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile
    integer(i4b),                              intent(out)          :: nside, ordering
    real(dp),         pointer, dimension(:,:)                       :: map   
    character(len=80),         dimension(180)                       :: header
    real(dp),                  dimension(4),   intent(in), optional :: md

    integer(i4b) :: npix, nside_mask, ordering_mask, nmaps, i, j, k, degree
    real(dp)     :: nullval, tot
    real(dp)     :: fmissval 
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)       :: mask, harmonics, harmonics2
    real(dp),                 dimension(0:3)       :: multipoles, b
    real(dp),                 dimension(3)         :: vector
    real(dp),                 dimension(0:3,0:3)   :: A
    real(dp),                 dimension(2)         :: zbounds = 0.d0

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if
    ordering = 1

    if (present(md)) then
       write(*,*) 'md input = ', real(md,sp)
       do i = 0, npix-1
          if (ordering == 1) then
             call pix2vec_ring(nside, i, vector)
          else
             call pix2vec_nest(nside, i, vector)          
          end if
          map(i,1) = map(i,1) - md(1)
          do j = 1, 3
             map(i,1) = map(i,1) - md(j+1) * vector(j)
          end do

!          map(i,1) = md(1)
!          do j = 1, 3
!             map(i,1) = map(i,1) + md(j+1) * vector(j)
!          end do
       end do
       return
    end if

    ! Read rms file and check consistency
    i = getsize_fits(maskfile, nside=nside_mask, ordering=ordering_mask)

    if (nside_mask /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    allocate(mask(0:npix-1,1))
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)

    if (ordering_mask == 2) then
       call convert_nest2ring(nside, mask(:,1))
    end if

    allocate(harmonics(0:npix-1,0:3))
    allocate(harmonics2(0:npix-1,0:3))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2vec_ring(nside, i, vector)
       else
          call pix2vec_nest(nside, i, vector)          
       end if
          
       if (mask(i,1) == 1.) then
          harmonics(i,0) = 1.d0
          harmonics(i,1) = vector(1)
          harmonics(i,2) = vector(2)
          harmonics(i,3) = vector(3)
       else
          harmonics(i,:) = 0.d0
       end if

       harmonics2(i,0) = 1.d0
       harmonics2(i,1) = vector(1)
       harmonics2(i,2) = vector(2)
       harmonics2(i,3) = vector(3)

    end do

    A = 0.d0
    b = 0.d0
    do j = 0, 3
       do k = 0, 3
          A(j,k) = sum(harmonics(:,j) * harmonics(:,k))
       end do
       b(j) = sum(map(:,1) * harmonics(:,j))
    end do


    call solve_system_real(A, multipoles, b)

    do i = 0, npix-1
!       if (mask(i,1) == 1) then
          do j = 0, 3
             map(i,1) = map(i,1) - multipoles(j) * harmonics2(i,j)
          end do
!       else
!          map(i,:) = -1.6375e30
!       end if
    end do

    write(*,*) 'Coefficients = ', real(multipoles,sp)

!    degree = 2
!    call remove_dipole(nside, map(:,1), ordering, degree, multipoles, zbounds, fmissval, mask(:,1))

!    write(*,*) 'multipoles = ', multipoles

    deallocate(mask)

  end subroutine subtract_mono_dipole


  subroutine subtract_mono_dipole_highl(mapfile, maskfile, lmax, lcut, nside, ordering, map, header, md)
    implicit none

    character(len=*),                          intent(in)           :: mapfile, maskfile
    integer(i4b),                              intent(in)           :: lmax, lcut
    integer(i4b),                              intent(out)          :: nside, ordering
    real(dp),         pointer, dimension(:,:)                       :: map   
    real(dp),                  dimension(4),   intent(in), optional :: md
    character(len=80),         dimension(180)                       :: header

    integer(i4b) :: npix, nside_mask, ordering_mask, nmaps, i, j, k, l, m, degree, numcomp
    real(dp)     :: nullval, tot
    real(dp)     :: nullval_sp
    real(dp)     :: fmissval 
    logical(lgt) :: anynull
    type(planck_rng) :: rng_handle
    real(dp),    allocatable, dimension(:,:)       :: mask, harmonics
    complex(dpc), allocatable, dimension(:,:,:)    :: alms
    real(dp),                 dimension(3)         :: vector
    real(dp), allocatable,    dimension(:)         :: b, multipoles
    real(dp), allocatable,    dimension(:,:)       :: A
    real(dp),                 dimension(2)         :: zbounds = 0.d0

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if
    ordering = 1

    ! Read rms file and check consistency
    i = getsize_fits(maskfile, nside=nside_mask, ordering=ordering_mask)

    if (nside_mask /= nside) then
       write(*,*) 'Error: Different Nsides of the two files.'
       stop
    end if

    allocate(mask(0:npix-1,1))
    call read_bintab(maskfile, mask, npix, 1, nullval, anynull)

    if (ordering_mask == 2) then
       call convert_nest2ring(nside, mask(:,1))
    end if

    numcomp = (lmax+1)**2
    allocate(harmonics(0:npix-1,numcomp), A(numcomp,numcomp), b(numcomp), multipoles(numcomp), alms(1,0:lmax,0:lmax))
    do i = 0, npix-1
       if (ordering == 1) then
          call pix2vec_ring(nside, i, vector)
       else
          call pix2vec_nest(nside, i, vector)          
       end if
       harmonics(i,1) = 1.d0
       harmonics(i,2) = vector(1)
       harmonics(i,3) = vector(2)
       harmonics(i,4) = vector(3)
    end do

    i = 5
    do l = 2, lmax
       do m = 0, l
          write(*,*) l, m
          alms = cmplx(0.d0, 0.d0)
          alms(1,l,m) = 1.d0
          call alm2map(nside, lmax, lmax, alms, harmonics(:,i))
          i = i+1
          if (m > 0) then
             alms(1,l,m) = cmplx(0.d0, 1.d0)
             call alm2map(nside, lmax, lmax, alms, harmonics(:,i))
             i = i+1
          end if
       end do
    end do

    A = 0.d0
    b = 0.d0
    do j = 1, numcomp
       write(*,*) numcomp, j
       do k = j, numcomp
          A(j,k) = sum(harmonics(:,j) * mask(:,1) * harmonics(:,k))
          A(k,j) = A(j,k)
       end do
       b(j) = sum(map(:,1) * mask(:,1) * harmonics(:,j))
    end do


    call solve_linear_system(A, multipoles, b)

    do i = 0, npix-1
       if (map(i,1) > -1.637d30) then
          do j = 1, (lcut+1)**2
             map(i,1) = map(i,1) - multipoles(j) * harmonics(i,j)
          end do
       end if
    end do

    write(*,*) 'Coefficients = ', real(multipoles(1:4),sp)

!    degree = 2
!    call remove_dipole(nside, map(:,1), ordering, degree, multipoles, zbounds, fmissval, mask(:,1))

!    write(*,*) 'multipoles = ', multipoles

    deallocate(mask)

  end subroutine subtract_mono_dipole_highl

  subroutine extract_multipole_range(mapfile, lmin, lmax, nside, ordering, map, header)
    implicit none

    character(len=*),                          intent(in)           :: mapfile
    integer(i4b),                              intent(in)           :: lmin, lmax
    integer(i4b),                              intent(out)          :: nside, ordering
    real(dp),         pointer, dimension(:,:)                       :: map   
    character(len=80),         dimension(180)                       :: header

    integer(i4b) :: npix, nmaps, i, j, k
    real(dp)     :: nullval, tot
    real(dp)     :: fmissval 
    logical(lgt) :: anynull
    logical(lgt), allocatable, dimension(:,:)       :: mask
    real(dp),     allocatable, dimension(:,:)       :: map_lowl
    real(dp),                  dimension(2)         :: zbounds = 0.d0
    complex(dpc), allocatable, dimension(:,:,:)     :: alms
    real(dp),     pointer,     dimension(:,:)      :: weights

    ! Read input map
    i = getsize_fits(mapfile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(:,i))
       end do
    end if
    ordering = 1

    ! Set undefined values to zero
    allocate(mask(0:npix,nmaps))
    mask = .true.
    do j = 1, nmaps
       do i = 0, npix-1
          if (map(i,j) == -1.6375e30) then
             mask(i,j) = .false.
             map(i,j)  = 0.d0
          end if
       end do
    end do

    ! Read Healpix ring weights
    call read_ringweights(nside, nmaps, weights)

    ! Compute spherical harmonic expansion
    allocate(alms(nmaps,0:lmax,0:lmax))
    if (nmaps == 1) then
       call map2alm(nside, lmax, lmax, map(:,1), alms, zbounds, weights)
    else
       call map2alm(nside, lmax, lmax, map, alms, zbounds, weights)
    end if

    alms(:,0:lmin-1,0:lmin-1) = cmplx(0.d0, 0.d0)

    ! Compute low-l map
    allocate(map_lowl(0:npix-1,nmaps))
    if (nmaps == 1) then
       call alm2map(nside, lmax, lmax, alms, map(:,1))
    else
       call alm2map(nside, lmax, lmax, alms, map)       
    end if

    ! Apply mask
    do j = 1, nmaps
       do i = 0, npix-1
          if (.not. mask(i,j)) then
             map(i,j)  = -1.6375d30
          end if
       end do
    end do
    

  end subroutine extract_multipole_range



  subroutine print_map_to_ascii(mapfile_in, maskfile, mapfile_out)
    implicit none

    character(len=128), intent(in) :: mapfile_in, maskfile, mapfile_out

    integer(i4b) :: npix, nside, ordering, ordering2, nmaps, i, j, c
    real(dp)     :: nullval, mu, sigma
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull

    real(dp), allocatable, dimension(:,:) :: map, mask

    ! Read input map
    i = getsize_fits(mapfile_in, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile_in, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do i = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))             
          end if
       end do
    end if

    open(48,file=trim(mapfile_out))
    mu = 0.d0
    c  = 1
    do i = 0, npix-1
       if (mask(i,c) /= 0.) then
          write(48,*) map(i,:)
          mu = mu + map(i,1)
       end if
    end do
    close(48)

    mu = mu / sum(mask(:,c))
    sigma = 0.d0
    do i = 0, npix-1
       if (mask(i,c) /= 0.d0) then
          sigma = sigma + (map(i,1)-mu)**2
       end if
    end do
    sigma = sqrt(sigma/(sum(mask(:,c))-1))

    write(*,*) 'Average outside mask = ', mu, ' +/- ', sigma

    deallocate(map)

  end subroutine print_map_to_ascii

  subroutine print_two_maps_to_ascii(mapfile1, mapfile2, maskfile, mapfile_out)
    implicit none

    character(len=128), intent(in) :: mapfile1, mapfile2, maskfile, mapfile_out

    integer(i4b) :: npix, nside, ordering, nmaps, i, j, ordering2
    real(dp)     :: nullval
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull

    real(dp), allocatable, dimension(:,:) :: map1, map2, mask

    ! Read input map
    i = getsize_fits(mapfile1, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map1(0:npix-1,nmaps))
    call read_bintab(mapfile1, map1, npix, nmaps, nullval, anynull)

    ! Read input map
    i = getsize_fits(mapfile2, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map2(0:npix-1,nmaps))
    call read_bintab(mapfile2, map2, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do j = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, map2(:,j))
          else
             call convert_nest2ring(nside, map2(:,j))
          end if
       end do
    end if

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do j = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,j))
          else
             call convert_nest2ring(nside, mask(:,j))
          end if
       end do
    end if

    open(48,file=trim(mapfile_out))
    do j = 1, 1 !nmaps
       do i = 0, npix-1
          if (mask(i,j) /= 0. .and. mask(i,j) /= 0.d0 .and. map1(i,j) /= -1.6375e30 .and. &
               & map2(i,j) /= -1.6375e30 .and. map1(i,j) /= -1.6375d30 .and. &
               & map2(i,j) /= -1.6375d30) write(48,*) map1(i,j), map2(i,j)
       end do
       write(48,*)
    end do
    close(48)

    deallocate(map1, map2)

  end subroutine print_two_maps_to_ascii

  subroutine print_isolatitude(mapfile_in, maskfile, outfile)
    implicit none

    character(len=128), intent(in) :: mapfile_in, maskfile, outfile

    integer(i4b) :: npix, nside, ordering, ordering2, nmaps, i, j, n, nlist, nring
    real(dp)     :: nullval, mu, sigma, dtheta
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull
    integer(i4b), allocatable, dimension(:) :: pixlist

    real(dp), allocatable, dimension(:,:) :: map, mask

    ! Read input map
    i = getsize_fits(mapfile_in, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps))
    call read_bintab(mapfile_in, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do i = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))             
          end if
       end do
    end if

    nring  = 35
    dtheta = pi/nring
    nring = nint(pi / dtheta) - 1
    allocate(pixlist(0:npix-1))
    open(58,file=trim(outfile))
    do i = 1, nring
       call query_strip(nside, (i-1)*dtheta, i*dtheta, pixlist, nlist)
       
       mu = 0.d0
       n  = 0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             mu = mu + map(pixlist(j),1)
             n  = n+1
          end if
       end do
       if (n > 0) then
          mu = mu / n
       else
          mu = 0.d0
       end if

       sigma = 0.d0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             sigma = sigma + (map(pixlist(j),1)-mu)**2
          end if
       end do
       sigma = sqrt(sigma/(n-1))

       write(58,*) (i-0.5d0)*dtheta*180.d0/pi, mu, sigma

    end do
    close(58)

  end subroutine print_isolatitude


  subroutine print_isolatitude_var(clfile, maskfile, outfile)
    implicit none

    character(len=128), intent(in) :: clfile, maskfile, outfile

    integer(i4b) :: npix, nside, ordering, ordering2, nmaps, i, j, n, nlist, nring
    real(dp)     :: nullval, mu, sigma, dtheta
    real(dp)     :: fmissval = 0.
    logical(lgt) :: anynull
    integer(i4b), allocatable, dimension(:) :: pixlist

    real(dp), allocatable, dimension(:,:) :: map, mask

!!$    ! Read input map
!!$    i = getsize_fits(mapfile_in, nside=nside, ordering=ordering, nmaps=nmaps)
!!$    npix = nside2npix(nside)
!!$    allocate(map(0:npix-1,nmaps))
!!$    call read_bintab(mapfile_in, map, npix, nmaps, nullval, anynull)

    ! Read input mask
    i = getsize_fits(maskfile, nside=nside, ordering=ordering2, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(mask(0:npix-1,nmaps))
    call read_bintab(maskfile, mask, npix, nmaps, nullval, anynull)
    if (ordering2 /= ordering) then
       do i = 1, nmaps
          if (ordering2 == 1) then
             call convert_ring2nest(nside, mask(:,i))
          else
             call convert_nest2ring(nside, mask(:,i))             
          end if
       end do
    end if

    nring  = 35
    dtheta = pi/nring
    nring = nint(pi / dtheta) - 1
    allocate(pixlist(0:npix-1))
    open(58,file=trim(outfile))
    do i = 1, nring
       call query_strip(nside, (i-1)*dtheta, i*dtheta, pixlist, nlist)
       
       mu = 0.d0
       n  = 0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             mu = mu + map(pixlist(j),1)
             n  = n+1
          end if
       end do
       if (n > 0) then
          mu = mu / n
       else
          mu = 0.d0
       end if

       sigma = 0.d0
       do j = 0, nlist-1
          if (mask(pixlist(j),1) > 0.5d0) then
             sigma = sigma + (map(pixlist(j),1)-mu)**2
          end if
       end do
       sigma = sqrt(sigma/(n-1))

       write(58,*) (i-0.5d0)*dtheta*180.d0/pi, mu, sigma

    end do
    close(58)

  end subroutine print_isolatitude_var


  subroutine summarize_detector_angles(filelist, s_max, nside, prefix)
    implicit none

    integer(i4b),     intent(in) :: s_max, nside
    character(len=*), intent(in) :: filelist, prefix

    integer(i4b)       :: unit1, unit2, npix, pix, s
    character(len=1)   :: s_string
    character(len=256) :: filename
    real(dp), allocatable, dimension(:,:,:) :: map
    real(dp), allocatable, dimension(:,:)   :: nhits
    real(dp),              dimension(6)     :: ang
    character(len=80), dimension(180)       :: header

    unit1 = 21
    unit2 = 22
    npix  = 12*nside**2

    allocate(map(0:npix-1,2,0:s_max))
    allocate(nhits(0:npix-1,1))

    map   = 0.d0
    nhits = 0.d0
    open(unit1,file=trim(filelist))
    do while (.true.)
       read(unit1,'(a)',end=990) filename
       open(unit2,file=trim(filename))
       write(*,*) 'Processing ', trim(filename)
       do while (.true.)
          read(unit2,*,end=991) ang
          ang = ang * pi / 180.d0

          ! Add first set
          call ang2pix_ring(nside, ang(2), ang(1), pix)
          nhits(pix,1) = nhits(pix,1) + 1.d0
          do s = 0, s_max
             map(pix,1,s) = map(pix,1,s) + cos(s*ang(3))
             map(pix,2,s) = map(pix,2,s) + sin(s*ang(3))
          end do

          ! Add second set
          call ang2pix_ring(nside, ang(5), ang(4), pix)
          nhits(pix,1) = nhits(pix,1) + 1.d0
          do s = 0, s_max
             map(pix,1,s) = map(pix,1,s) + cos(s*ang(6))
             map(pix,2,s) = map(pix,2,s) + sin(s*ang(6))
          end do
       end do
991    close(unit2)
    end do
990 close(unit1)

    filename = trim(prefix) // '_nhits.fits'
    call write_minimal_header(header, 'MAP', nside=nside, ordering='RING', polar=.false.)
    call write_result_map(filename, nside, 1, header, nhits)

    do s = 0, s_max
       write(*,*) 'c', s
       map(:,1,s) = map(:,1,s) / nhits(:,1)
       map(:,2,s) = map(:,2,s) / nhits(:,1)

       call int2string(s, s_string)
       filename = trim(prefix) // '_s' // s_string // '.fits'
       call write_minimal_header(header, 'MAP', nside=nside, ordering='RING')
       call write_result_map(filename, nside, 1, header, map(:,:,s))
    end do

    

  end subroutine summarize_detector_angles

  subroutine ud_grade_map_editor(filename_in, nside_out, filename_out, double_precision)
    implicit none

    character(len=*), intent(in) :: filename_in, filename_out
    integer(i4b),     intent(in) :: nside_out
    logical(lgt),     intent(in) :: double_precision

    real(dp)     :: nullval
    logical(lgt) :: anynull
    integer(i4b) :: i, j, nside_in, ordering, nmaps, npix_in, npix_out
    real(dp), allocatable, dimension(:,:) :: map_in, map_out
    character(len=80), dimension(180)     :: header

    if (iargc() /= 4) then
       write(*,*) 'Usage: map_editor ud_grade [input map] [nside_out] [output map]'
       stop
    end if

    i = getsize_fits(filename_in, nside=nside_in, ordering=ordering, nmaps=nmaps)
    npix_in  = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)

    allocate(map_in(0:npix_in-1,nmaps))
    allocate(map_out(0:npix_out-1,nmaps))
    call read_bintab(filename_in, map_in, npix_in, nmaps, nullval, anynull, header=header)
    where (map_in < -1.6d30) 
       map_in = -1.6375d30
    end where

    if (ordering == 1) then
       call udgrade_ring(map_in, nside_in, map_out, nside_out)
    else
       call udgrade_nest(map_in, nside_in, map_out, nside_out)
    end if

    call write_result_map(filename_out, nside_out, ordering, header, map_out, double_precision)

    if (count(map_in==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_in(:,i)==-1.6375d30), ' missing pixels out of ', npix_in
       end do
    end if
    if (count(map_out==-1.6375d30) > 0) then
       do i=1,nmaps
          write(*,*) count(map_out(:,i)==-1.6375d30), ' missing pixels out of ', npix_out
       end do
    end if
    deallocate(map_in)
    deallocate(map_out)

  end subroutine ud_grade_map_editor


  subroutine read_covmatrix(unit, filename, ordering, polarization, matrix, inv, n)

    integer(i4b),                          intent(in)  :: unit
    integer(i8b),                          intent(out) :: n
    integer(i4b),                          intent(out) :: ordering, polarization
    character(len=*),                      intent(in)  :: filename
    logical(lgt),                          intent(out) :: inv
    real(dp), allocatable, dimension(:,:), intent(out) :: matrix

    integer(i4b) :: n_in, i

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    open(unit, file=trim(filename), form='unformatted')
    read(unit) n_in
    n = n_in
    write(*,*) n, '= n'
    read(unit) ordering
    write(*,*) ordering, '= ordering'
    read(unit) polarization
    write(*,*) polarization, '= polarisation'
    allocate(matrix(n,n))
    do i = 1, n
       read(unit) matrix(:,i)
    end do
    read(unit) inv
    close(unit)

  end subroutine read_covmatrix


  subroutine make_co_region_map(mapname, nside, ordering, reg)
    implicit none

    character(len=*),                         intent(in)  :: mapname
    integer(i4b),                             intent(out) :: nside, ordering
    real(dp),         pointer, dimension(:,:)             :: reg

    integer(i4b) :: npix, nmaps, i, j, k, lmax, nreg, p(1), listpix(0:100000), nlist, r
    real(dp)     :: nullval, tot, fmissval, threshold, radius, vec(3), threshold0, dist, dist0
    logical(lgt) :: anynull
    real(dp),     allocatable, dimension(:,:)       :: map, map_seed, seeds

    threshold0 = 1.d0 ! K km/s
    radius     = 10.d0 * pi/180.d0 ! Distance between seeds in radians

    ! Read input map
    i = getsize_fits(mapname, nside=nside, ordering=ordering)
    npix = nside2npix(nside)
    nmaps = 1
    allocate(map(0:npix-1,nmaps), reg(0:npix-1,nmaps), map_seed(0:npix-1,nmaps), seeds(10000,3))
    call read_bintab(mapname, map, npix, nmaps, nullval, anynull)
    if (ordering == 2) call convert_nest2ring(nside, map(:,1))
    ordering = 1

    nreg = 0
    reg     = -1.6375d30

    ! Start by setting all low-amplitude pixels to region 1
    nreg = nreg+1
    where (map < threshold0) 
       reg = nreg
       map = -1.d30
    end where

    ! Set up a grid of seed points
    map_seed = map
    do while (any (map_seed >= 0.d0))
       p = maxloc(map_seed(:,1))
       call pix2vec_ring(nside, p(1), vec)
       call query_disc(nside, vec, 2.d0*radius, listpix, nlist)       

       nreg = nreg+1
       seeds(nreg,:) = vec
       reg(p,1)  = nreg
       map_seed(listpix(0:nlist-1),1) = -1.d30
       write(*,*) 'p = ', p, ', nreg = ', nreg
    end do

    ! Put remaining pixels in closest seed
    do i = 0, npix-1
       if (reg(i,1) > 0.d0) cycle
       call pix2vec_ring(nside, i, vec)

       dist0 = 1.d30
       do j = 2, nreg
          call angdist(vec, seeds(j,:), dist)
          if (dist < dist0) then
             reg(i,1) = j
             dist0    = dist
          end if
       end do
    end do

    deallocate(map, map_seed)

  end subroutine make_co_region_map



  ! Parallel transport polarization angles from a given map to centers
  ! of a lower resolution map at nside_out. 
  ! Note: Maps must be in nested format
  subroutine qu_transport_map(nside_out, map)
    implicit none

    integer(i4b),                   intent(in)    :: nside_out
    real(dp),     dimension(0:,1:), intent(inout) :: map

    integer(i4b) :: i, j, q
    real(dp)     :: cos2psi, sin2psi, m(2)
    integer(i4b), save :: nside, npix, nmaps
    real(dp), allocatable, dimension(:,:,:), save :: vec

    ! Precompute pixel vectors
    if (.not. allocated(vec)) then
       npix  = size(map,1)
       nmaps = size(map,2)
       nside = npix2nside(npix)
       q     = (nside/nside_out)**2
       if (nmaps /= 3) then
          write(*,*) 'partrans -- nmaps must be equal to 3'
          stop
       end if
       allocate(vec(0:npix-1,3,2))
       do i = 0, npix-1
          call pix2vec_nest(nside_out, i/q, vec(i,:,1)) ! Low-res pixel centers
          call pix2vec_nest(nside,     i,   vec(i,:,2)) ! High-res pixel centers
       end do
    end if

    ! Perform the parallel transport
    do i = 0, npix-1
       call qu_transport_rot(vec(i,:,1), vec(i,:,2), cos2psi, sin2psi)
       m        = map(i,2:3)
       map(i,2) =  cos2psi*m(1) - sin2psi*m(2)
       map(i,3) =  sin2psi*m(1) + cos2psi*m(2)
    end do

  end subroutine qu_transport_map

  subroutine merge_maps(double_precision)
    implicit none

    logical(lgt), intent(in) :: double_precision

    integer(i4b)       :: npix
    integer(i4b)       :: i, j, l, m
    real(dp)           :: nullval
    real(dp)           :: sigma_sq, nullval_dp, fwhm_trans, threshold, cl
    logical(lgt)       :: anynull
    character(len=128) :: beamfile1, beamfile2, infile1, infile2, outfile, partext, prefix, maskfile
    integer(i4b)       :: ltrans_start, ltrans_stop, ordering, nmaps
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside1, nside2, npix1, npix2

    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     pointer,     dimension(:,:)   :: weights1, weights2
    real(dp),     allocatable, dimension(:,:)   :: beam1, beam2, map1, map2, map, map_out, w_l
    character(len=80),         dimension(1:180) :: header
    real(dp),                  dimension(2)     :: zbounds = 0.d0
    real(dp),     pointer,     dimension(:,:)   :: pixwin1, pixwin2

    if (iargc() /= 10) then
       write(*,*) 'Usage: map_editor merge_maps [map_lowl] [beam_lowl] [map_highl] [beam_highl]'
       write(*,*) '                    [l_trans_start] [l_trans_stop] [threshold] [outmap] [prefix]'
       stop
    end if

    call getarg(2,infile1)
    call getarg(3,beamfile1)
    call getarg(4,infile2)
    call getarg(5,beamfile2)
    call getarg(6,partext)
    read(partext,*) ltrans_start
    call getarg(7,partext)
    read(partext,*) ltrans_stop
    call getarg(8,partext)
    read(partext,*) threshold
    call getarg(9,outfile)
    call getarg(10,prefix)

    ! Read input data set 1
    i = getsize_fits(infile1, nside=nside1, ordering=ordering, nmaps=nmaps)
    npix1 = nside2npix(nside1)
    allocate(map1(0:npix1-1,nmaps))
    call read_bintab(infile1, map1, npix1, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside1, map1(0:npix-1,i))
       end do
    end if
    call read_ringweights(nside1, nmaps, weights1)
    allocate(beam1(0:ltrans_stop, nmaps))
    call generate_beam(0.d0, ltrans_stop, beam1, beamfile1)
    call read_pixwin(nside1, nmaps, pixwin1)

    ! Read input data set 2
    i = getsize_fits(infile2, nside=nside2, ordering=ordering, nmaps=nmaps)
    npix2 = nside2npix(nside2)
    allocate(map2(0:npix2-1,nmaps), map_out(0:npix2-1,nmaps))
    call read_bintab(infile2, map2, npix2, nmaps, nullval, anynull)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside2, map2(0:npix2-1,i))
       end do
    end if
    call read_ringweights(nside2, nmaps, weights2)
    allocate(beam2(0:ltrans_stop, nmaps))
    call generate_beam(0.d0, ltrans_stop, beam2, beamfile2)
    call read_pixwin(nside2, nmaps, pixwin2)

    ! Threshold maps
    if (.false. .and. threshold > 0.d0) then
       map1 = max(min(map1,threshold),-threshold)
       map2 = max(min(map2,threshold),-threshold)
    end if

    ! Compute transition beam, and effective weight function
    allocate(w_l(0:ltrans_stop,nmaps))
    do i = 1, nmaps
       do l = 0, ltrans_stop
          ! Gaussian smoothing
          w_l(l,i) = 1.d0 
          if (l > ltrans_start) then
             ! Cosine apodization
             w_l(l,i) = w_l(l,i) * 0.5d0*(1.d0-cos(pi*real(ltrans_stop-l,dp)/real(ltrans_stop-ltrans_start,dp)))          
          end if
       end do
    end do

    ! Compute the spherical harmonics transform
    allocate(alms(nmaps, 0:ltrans_stop, 0:ltrans_stop))
    allocate(map(0:npix1-1,nmaps))
    map = map1
    !if (threshold > 0) map = max(min(map,threshold),-threshold)
    if (nmaps == 1) then
       call map2alm(nside1, ltrans_stop, ltrans_stop, map(:,1), alms, zbounds, weights1)
    else
       call map2alm(nside1, ltrans_stop, ltrans_stop, map, alms, zbounds, weights1)
    end if
    deallocate(map)

    ! Deconvolve old beam, convolve with new
    do i = 1, nmaps
       do l = 0, ltrans_stop
          if (beam1(l,i)*pixwin1(l,i) > 0.d0) then
             alms(i,l,0:l) = alms(i,l,0:l) / (beam1(l,i)*pixwin1(l,i)) * (beam2(l,i)*pixwin2(l,i))
          end if
       end do
    end do

    open(58,file=trim(prefix)//'_cls1.dat')
    do l = 0, ltrans_stop
       cl = abs(alms(1,l,0))**2
       do m = 1, l
          cl = cl + 2.d0*abs(alms(1,l,m))**2
       end do
       write(58,*) l, cl/real(2*l+1,dp) * l*(l+1)/2/pi
    end do
    close(58)

    ! Cosine apodization
    open(58,file=trim(prefix)//'_apod.dat')
    do i = 1, nmaps
       do l = 0, ltrans_stop
          alms(i,l,0:l) = alms(i,l,0:l) * w_l(l,i)
          if (i == 1) write(58,*) l, w_l(l,i), 1.d0-w_l(l,i)
       end do
    end do
    close(58)

    ! Compute the inverse spherical harmonics
    if (nmaps == 1) then
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map_out(:,1))
    else
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map_out)
    end if
    deallocate(map1)
    allocate(map1(0:npix2-1,nmaps))
    map1 = map_out

    call write_result_map(trim(prefix)//'_comp1.fits', nside2, ordering, header, map_out, double_precision)


    ! Add second map
    map_out = map_out + map2
    !map_out = map2


    ! Compute the spherical harmonics transform
    allocate(map(0:npix2-1,nmaps))
    map = map2
    !if (threshold > 0) map = max(min(map,threshold),-threshold)
    if (nmaps == 1) then
       call map2alm(nside2, ltrans_stop, ltrans_stop, map(:,1), alms, zbounds, weights2)
    else
       call map2alm(nside2, ltrans_stop, ltrans_stop, map, alms, zbounds, weights2)
    end if
    deallocate(map)

    open(58,file=trim(prefix)//'_cls2.dat')
    do l = 0, ltrans_stop
       cl = abs(alms(1,l,0))**2
       do m = 1, l
          cl = cl + 2.d0*abs(alms(1,l,m))**2
       end do
       write(58,*) l, cl/real(2*l+1,dp) * l*(l+1)/2/pi
    end do
    close(58)

    ! Cosine apodization
    do i = 1, nmaps
       do l = 0, ltrans_stop
          alms(i,l,0:l) = alms(i,l,0:l) * w_l(l,i)
       end do
    end do

    ! Compute the inverse spherical harmonics
    if (nmaps == 1) then
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map2(:,1))
    else
       call alm2map(nside2, ltrans_stop, ltrans_stop, alms, map2)
    end if

    ! Subtract large scales to avoid double counting
    map_out = map_out - map2

    call write_result_map(trim(prefix)//'_comp2.fits', nside2, ordering, header, map_out-map1, double_precision)

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside2, map_out(:,i))
       end do
    end if

    call write_result_map(outfile, nside2, ordering, header, map_out, double_precision)

    if (nmaps == 1) then
       call map2alm(nside2, ltrans_stop, ltrans_stop, map_out(:,1), alms, zbounds, weights2)
    else
       call map2alm(nside2, ltrans_stop, ltrans_stop, map_out, alms, zbounds, weights2)
    end if

    open(58,file=trim(prefix)//'_cls.dat')
    do l = 0, ltrans_stop
       cl = abs(alms(1,l,0))**2
       do m = 1, l
          cl = cl + 2.d0*abs(alms(1,l,m))**2
       end do
       write(58,*) l, cl/real(2*l+1,dp) * l*(l+1)/2/pi
    end do
    close(58)

  end subroutine merge_maps

  subroutine median_filter_holes
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax
    real(dp)           :: nullval, mean1, mean2
    real(dp)           :: nullval_dp, radius, vec(3), threshold
    logical(lgt)       :: anynull
    character(len=128) :: infile, outfile, partext
    integer(i4b)       :: ordering, nmaps
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix

    real(dp),     pointer,     dimension(:,:)   :: weights
    complex(dpc), allocatable, dimension(:,:,:) :: alms
    real(dp),     allocatable, dimension(:,:)   :: map, outmap, mask, mask2, mask3
    character(len=80),         dimension(1:180) :: header
    integer(i4b), allocatable, dimension(:)     :: listpix1, listpix2

    if (iargc() /= 5) then
       write(*,*) 'Usage: map_editor median_filter_source_holes [map_in] [radius arcmin] [threshold] [map_out]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) radius ! In arcmin
    call getarg(4,partext)
    read(partext,*) threshold ! In arcmin
    call getarg(5,outfile)
    radius = radius/60.d0 * pi/180.d0

    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps), mask(0:npix-1,nmaps), mask2(0:npix-1,nmaps), mask3(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if

    ! Loop through all pixels, searching for holes
    allocate(listpix1(0:1000000), listpix2(0:1000000))
    do j = 1, nmaps
       do i = 0, npix-1
          if (mod(i,1000000) == 0) write(*,*) j, i, npix
          call pix2vec_ring(nside, i, vec)
          call query_disc(nside, vec, radius, listpix1, n1)
          mean1 = mean(map(listpix1(0:n1-1),j))
          if (map(i,j) < 1.d-4 * mean1 .and. mean1 > threshold) then
             mask(i,j)   = 1.d0 
          else
             !outmap(i,j) = map(i,j)
             mask(i,j)   = 0.d0 
          end if
       end do
    end do

    ! Expand relevant region by four radii
    mask2 = 0.d0
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,j) == 1.d0) then
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, 4.d0*radius, listpix1, n1)             
             mask2(listpix1(0:n1-1),j) = 1.d0
             call query_disc(nside, vec, 2.d0*radius, listpix1, n1)             
             mask3(listpix1(0:n1-1),j) = 1.d0
          end if
       end do
    end do

    ! Median filter selected pixels over 2 radii
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask2(i,j) == 1.d0) then
             call pix2vec_ring(nside, i, vec)
             call query_disc(nside, vec, 2.d0*radius, listpix1, n1)
             outmap(i,j) = median(map(listpix1(0:n1-1),j))
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do

    ! Apodize mask
    lmax = 2*nside
    mask = mask3
    allocate(alms(1,0:lmax,0:lmax))
    call read_ringweights(nside, nmaps, weights)    
    call map2alm(nside, lmax, lmax, mask(:,1), alms, [0.d0,0.d0], weights)
    do l = 0, lmax
       alms(1,l,0:l) = alms(1,l,0:l) * &
            & exp(-0.5d0*l*(l+1.d0)*(30.d0*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
    end do
    call alm2map(nside, lmax, lmax, alms, mask(:,1))

    if (ordering == 2) then
       do i = 1, nmaps
          call convert_ring2nest(nside, outmap(0:npix-1,i))
       end do
       call convert_ring2nest(nside, mask(0:npix-1,1))
    end if

    ! Weight pixels
    do j = 1, nmaps
       do i = 0, npix-1
          if (mask(i,1) > 1.d-2) then
             outmap(i,j) = (1.d0-mask(i,1))*map(i,j) + mask(i,1) * outmap(i,j)
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do

    ! Output result map
    call write_minimal_header(header, 'MAP', nside=nside, order=ordering, polar=nmaps==3)
    call write_result_map('apomask_'//trim(outfile), nside, ordering, header, mask)
    call write_result_map('mask_'//trim(outfile), nside, ordering, header, mask2)
    call write_result_map(outfile, nside, ordering, header, outmap)

  end subroutine median_filter_holes



  subroutine median_filter_misspix
    implicit none

    integer(i4b)       :: i, j, k, l, m, n1, n2, lmax, nummiss, r, nlist
    real(dp)           :: nullval, mean1, mean2, theta, phi
    real(dp)           :: nullval_dp, radius, vec(3), threshold
    logical(lgt)       :: anynull, anymiss, dpcfix
    character(len=128) :: infile, outfile, partext
    integer(i4b)       :: ordering, nmaps, nref
    character(len=4)   :: ntext
    type(planck_rng)   :: rng_handle

    integer(i4b) :: nside, npix, r_fill
    real(dp)     :: rf, tot

    real(dp),     allocatable, dimension(:,:)   :: map, outmap
    character(len=80),         dimension(1:180) :: header
    real(dp),     allocatable, dimension(:)     :: buffer
    integer(i4b), allocatable, dimension(:)     :: listpix

    if (iargc() < 4) then
       write(*,*) 'Usage: map_editor median_filter_misspix [map_in] [radius in pixels] [map_out] [[dpc]]'
       stop
    end if

    call getarg(2,infile)
    call getarg(3,partext)
    read(partext,*) r_fill ! In pixels
    call getarg(4,outfile)
    dpcfix = .false.
    if (iargc() > 4) then
       call getarg(5,partext)
       read(partext,*) dpcfix
    end if

    ! Read input data set 1
    i = getsize_fits(infile, nside=nside, ordering=ordering, nmaps=nmaps)
    npix = nside2npix(nside)
    allocate(map(0:npix-1,nmaps), outmap(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull, header=header)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, map(0:npix-1,i))
       end do
    end if
    write(*,*) "LFI dpc =", dpcfix, "nmaps =", nmaps

    ! Loop through all pixels, searching for holes
    allocate(listpix(0:1000000), buffer(1000000))
    anymiss = .false.
    do j = 1, nmaps
       nref = j
       if (dpcfix) then
          if (nmaps==3 .and. j==3) then
             nref = 1
          else if (nmaps==10 .and. j==5) then
             nref = 1
          else if (nmaps==10 .and. j==8) then
             nref = 2
          else if (nmaps==10 .and. j==10) then
             nref = 3
          end if
       end if
       do i = 0, npix-1
          if (abs(map(i,nref)) > 1e30) then
             anymiss = .true.
             call pix2vec_ring(nside, i, vec)
             do r = 1, r_fill
                rf = (r+0.2d0) * sqrt(4.d0*pi/npix)
                call query_disc(nside, vec, rf, listpix, nlist)
                tot = 0.d0
                m   = 0.d0
                do k = 0, nlist-1
                   if (abs(map(listpix(k),nref)) < 1.d30) then
                      tot = tot + map(listpix(k),j)
                      m   = m   + 1
                      buffer(m) = map(listpix(k),j)
                   end if
                end do
                if (m > 4) then
                   !outmap(i,j) = tot / m
                   outmap(i,j) = median(buffer(1:m))
                   exit
                end if
             end do

             if (r > r_fill) then
                call pix2ang_ring(nside, i, theta, phi)
                write(*,*) 'Error: Filter did not return valid pixel; increase r_fill'
!                write(*,*) '     p         = ', i
!                write(*,*) '     theta     = ', 90.d0-theta*180/pi
!                write(*,*) '     phi       = ', phi*180/pi
!                write(*,*) '     radius    = ', r_fill
!                write(*,*) '     neighbors = ', listpix(0:nlist-1)
!                write(*,*) '     vals      = ', map(listpix(0:nlist-1),j)
                stop
             end if
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do

!    if (anymiss) then
       do i = 1, nmaps
          write(*,*) i, count(abs(map(:,i)) > 1e30), ' missing pixels out of ', npix
          if (ordering==2) call convert_ring2nest(nside, outmap(0:npix-1,i))
       end do
       ! Output result map
       call write_result_map(outfile, nside, ordering, header, outmap)
!    end if

  end subroutine median_filter_misspix

end module map_editor_complex_ops_mod
