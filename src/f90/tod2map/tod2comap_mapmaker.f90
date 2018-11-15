module tod2comap_mapmaker
  !use comap_lx_mod
  use comap_map_mod
  use comap_acceptlist_mod
  use comap_patch_mod
  !use quiet_fft_mod
  use tod2comap_utils
  implicit none



contains


  subroutine initialize_mapmaker(map, tod, parfile, pinfo)
    implicit none
    type(tod_type), dimension(:), intent(in)    :: tod
    type(map_type),               intent(inout) :: map
    type(patch_info),             intent(in)    :: pinfo
    character(len=*)                            :: parfile

    integer(i4b) :: i, j, k, l, p, q, fs, st, ierr
    real(dp)     :: x_min, x_max, y_min, y_max, pad, temp, mean_dec
    real(8), parameter :: PI = 4*atan(1.d0)

    ! Set up map grid
    if (allocated(map%x)) return
    call get_parameter(0, parfile, 'MAP_NAME', par_string=map%name)
    call get_parameter(0, parfile, 'NUMFREQ', par_int=map%nfreq)
    call get_parameter(0, parfile, 'NUM_SIDEBAND', par_int=map%nsb)

    
    !fs = 1!200
    !st = tod%nsamp!-200 ! tod%nsamp
    pad = 0.5d0 ! degrees
    !write(*,*) pinfo%resolution
    map%dthetay = pinfo%resolution ! degrees (arcmin/60), resolution
    map%mean_el = 0.d0!mean(tod(:)%mean_el)
    map%mean_az = 0.d0!mean(tod(:)%mean_az)
    if (pinfo%fixed) then 
       mean_dec = pinfo%pos(2)
    else 
       mean_dec = 0.d0
       do i = 1, size(tod)
          if (allocated(tod(i)%point)) mean_dec = mean_dec + mean(tod(i)%point(2,:,1))
       end do
       call mpi_allreduce(MPI_IN_PLACE, mean_dec, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
       mean_dec = mean_dec/size(tod)
    end if
    !map%dthetax = map%dthetay/abs(cos(map%mean_el*PI/180.d0))
    map%dthetax = map%dthetay/abs(cos(mean_dec*PI/180.d0))

    x_min = 500.d0; x_max = -500.d0
    y_min = 500.d0; y_max = -500.d0

    if (pinfo%fixed) then
       x_min = pinfo%pos(1) - pinfo%obj_rad !- 30.d0*pad
       x_max = pinfo%pos(1) + pinfo%obj_rad !+ 40.d0*pad
       y_min = pinfo%pos(2) - pinfo%obj_rad !- 10.d0*pad
       y_max = pinfo%pos(2) + pinfo%obj_rad !+ 30.d0*pad
    else
       do i = 1, size(tod)
          if (.not. allocated(tod(i)%point)) cycle
          do j = 1, tod(i)%ndet
             if (all(tod(i)%point(1:2,:,j) == 0.d0)) cycle
             temp = minval(tod(i)%point(1,:,j)) - pad
             if (temp .le. x_min) x_min = temp
             temp = maxval(tod(i)%point(1,:,j)) + pad
             if (temp .ge. x_max) x_max = temp
             temp = minval(tod(i)%point(2,:,j)) - pad
             if (temp .le. y_min) y_min = temp
             temp = maxval(tod(i)%point(2,:,j)) + pad
             if (temp .ge. y_max) y_max = temp
          end do
       end do
       call mpi_allreduce(MPI_IN_PLACE, x_min, 1, mpi_double_precision, mpi_min, mpi_comm_world, ierr)
       call mpi_allreduce(MPI_IN_PLACE, x_max, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
       call mpi_allreduce(MPI_IN_PLACE, y_min, 1, mpi_double_precision, mpi_min, mpi_comm_world, ierr)
       call mpi_allreduce(MPI_IN_PLACE, y_max, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
    end if

!!$    write(*,*) x_min, x_max
!!$    write(*,*) y_min, y_max
!!$    call mpi_finalize(i)
!!$    stop

!!$    do k = 1, size(tod)
!!$       write(*,*) k, minval(tod(k)%point(1,:,:)), maxval(tod(k)%point(1,:,:))
!!$       write(*,*) k, minval(tod(k)%point(2,:,:)), maxval(tod(k)%point(2,:,:))
!!$    end do
!!$    call mpi_finalize(i)
!!$    stop
    

    !x_min = data%point_lim(1) - pad; x_max = data%point_lim(2) - pad
    !y_min = data%point_lim(3) - pad; y_max = data%point_lim(4) - pad
    !write(*,*) map%dthetax, map%dthetay
    map%n_x = int((x_max-x_min)/map%dthetax); map%n_y = int((y_max-y_min)/map%dthetay)
    !write(*,*) map%n_x, map%n_y

    allocate(map%x(map%n_x), map%y(map%n_y))
    do i = 1, map%n_x
       map%x(i) = x_min + (i-0.5d0)*map%dthetax
    end do
    do i = 1, map%n_y
       map%y(i) = y_min + (i-0.5d0)*map%dthetay
    end do
 
    ! Set up map structures
    allocate(map%m(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%dsum(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%nhit(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%div(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%rms(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%freq(map%nfreq, map%nsb))
    map%dsum = 0.d0
    map%nhit = 0.d0
    map%div  = 0.d0
    map%m    = 0.d0
    map%rms  = 0.d0
    map%freq = 0.d0
    
  end subroutine initialize_mapmaker


  subroutine time2pix(tod, map)
    implicit none
    type(tod_type), dimension(:), intent(inout) :: tod
    type(map_type),               intent(inout) :: map

    integer(i4b) :: scan, i, j, p, q, sb, freq
    real(dp)     :: x_min, x_max, y_min, y_max, t1, t2
    real(dp), allocatable, dimension(:,:) :: sigma0

    x_min = map%x(1); x_max = map%x(map%n_x)
    y_min = map%y(1); y_max = map%y(map%n_y)
    allocate(sigma0(tod(1)%nfreq,tod(1)%nsb))
    sigma0 = 0.d0

    
    !call wall_time(t1)
    !!$OMP PARALLEL PRIVATE(scan,j,i,p,q,sb,freq)
    !!$OMP DO SCHEDULE(guided)
    do scan = 1, size(tod)
       do j = 1, tod(scan)%ndet
          if (.not. is_alive(j)) cycle
          do i = 1, tod(scan)%nsamp
             if (tod(scan)%point(1,i,j) < x_min) write(*,*) "Expand your grid! (x_min)", scan, i,j,tod(scan)%point(1,i,j), x_min
             if (tod(scan)%point(1,i,j) > x_max) write(*,*) "Expand your grid! (x_max)", scan, i,j,tod(scan)%point(1,i,j), x_max
             if (tod(scan)%point(2,i,j) < y_min) write(*,*) "Expand your grid! (y_min)"
             if (tod(scan)%point(2,i,j) > y_max) write(*,*) "Expand your grid! (y_max)"
             p = min(max(nint((tod(scan)%point(1,i,j)-x_min)/map%dthetax),1),map%n_x)
             q = min(max(nint((tod(scan)%point(2,i,j)-y_min)/map%dthetay),1),map%n_y)
             tod(scan)%pixels(i,j) = (q-1)*map%n_x + p
             do sb = 1, tod(scan)%nsb
                do freq = 1, tod(scan)%nfreq
                   if (tod(scan)%freqmask(freq,sb,j) == 0) cycle
                   !!$OMP ATOMIC
                   map%nhit(p,q,freq,sb) = map%nhit(p,q,freq,sb) + 1.d0
                   sigma0(freq,sb) = sigma0(freq,sb) + tod(scan)%sigma0(freq,sb,j)
                end do
             end do
          end do
       end do
    end do
    !!$OMP END DO
    !!$OMP END PARALLEL
    !call wall_time(t2)
    !write(*,*) 'Wall time time2pix = ', t2-t1
    sigma0 = sigma0/size(tod)
    do sb = 1, tod(1)%nsb
       do freq = 1, tod(1)%nfreq
          where(map%nhit(:,:,freq,sb) > 0)
             map%rms(:,:,freq,sb) = map%rms(:,:,freq,sb) + map%nhit(:,:,freq,sb)/sigma0(freq,sb)**2
          elsewhere
             map%rms(:,:,freq,sb) = 0.d0
          end where
       end do
    end do
    !call mpi_finalize(i)
    !stop
    deallocate(sigma0)


  end subroutine time2pix


  
  !!!!!!!!!!!!!!!!!
  ! Naive Binning !
  !!!!!!!!!!!!!!!!!

  subroutine binning(map_tot, map_scan, tod, alist, scan)
    implicit none
    type(tod_type),   intent(in)    :: tod
    type(map_type),   intent(inout) :: map_tot, map_scan
    type(acceptlist), intent(in)    :: alist

    integer(i4b) :: det, sb, freq, ndet, nsb, nfreq
    integer(i4b) :: i, j, k, l, p, q, fs, st, scan, pix
    real(dp)     :: x_min, x_max, y_min, y_max
    real(dp), allocatable, dimension(:) :: dsum, div

    x_min = map_tot%x(1); x_max = map_tot%x(map_tot%n_x)
    y_min = map_tot%y(1); y_max = map_tot%y(map_tot%n_y)

    ndet = size(tod%d,4)
    nsb = map_tot%nsb 
    nfreq = map_tot%nfreq

    !write(*,*) 'Beginning coadding'

    !fs = 200 ! starting point
    !st = tod(scan)%nsamp - 200 ! ending point
    do i = 1, tod%nsamp
       do det = 1, ndet
          if (.not. is_alive(det)) cycle
          p = min(max(nint((tod%point(1,i,det)-x_min)/map_tot%dthetax),1),map_tot%n_x)
          q = min(max(nint((tod%point(2,i,det)-y_min)/map_tot%dthetay),1),map_tot%n_y)
          do sb = 1, nsb
             do freq = 1, nfreq
                !if (tod%fknee(freq,sb,det) > 0.2d0) cycle
                if (tod%freqmask(freq,sb,det) == 0) cycle
                !write(*,*) i, det, sb, freq
                !if (any(alist%ascans(scan)%adet_sb(det,sb)%rejected == freq)) cycle
                !write(*,*) tod%rms(i,freq,sb,det)
                map_scan%dsum(p,q,freq,sb) = map_scan%dsum(p,q,freq,sb) + 1.d0 / tod%rms(i,freq,sb,det)**2 * tod%d(i,freq,sb,det)
                map_scan%div(p,q,freq,sb) = map_scan%div(p,q,freq,sb) + 1.d0 / tod%rms(i,freq,sb,det)**2
                !end if
             end do
          end do
       end do
    end do
    map_tot%dsum = map_tot%dsum + map_scan%dsum
    map_tot%div  = map_tot%div  + map_scan%div

    !write(*,*) 'Ending coadding'

  end subroutine binning


  subroutine finalize_scan_binning(map)
    implicit none
    type(map_type), intent(inout) :: map
    
    integer(i4b) :: p, q

    do p = 1, map%n_x
       do q = 1, map%n_y
          map%dsum(p,q,1,1) = sum(map%dsum(p,q,:,:))
          map%div(p,q,1,1)  = sum(map%div(p,q,:,:))
       end do
    end do

    where(map%div > 0)
       map%m   = map%dsum / map%div
       map%rms = 1.d0 / sqrt(map%div)!sqrt(map%rms) ! overwrite
    elsewhere
       map%m   = 0.d0
       map%rms = 0.d0
    end where

  end subroutine finalize_scan_binning



  subroutine finalize_binning(map)
    implicit none
    type(map_type), intent(inout) :: map

    write(*,*) count(map%div > 0)
    where(map%div > 0)
       !write(*,*) map%dsum, map%div
       map%m   = map%dsum / map%div
       map%rms = 1.d0 / sqrt(map%div) ! overwrite
    elsewhere
       map%m   = 0.d0
       map%rms = 0.d0
    end where

  end subroutine finalize_binning





  !!!!!!!!!!!!!!
  ! PCG Solver !
  !!!!!!!!!!!!!!

  ! Computes b in Ax = b for the PCG mapmaker
  ! b = P^T N^{-1} d, d = g P s
   subroutine get_rhs(rhs, map, tod, det, sb, freq)
     implicit none
     type(map_type),  intent(inout) :: map
     integer(i4b),    intent(in)    :: det, sb, freq
     type(tod_type), dimension(:), intent(in)    :: tod
     real(dp),       dimension(:), intent(inout) :: rhs

     real(dp),     allocatable, dimension(:) :: dt, Nt ! intermediate stages
     complex(dpc), allocatable, dimension(:) :: dv, Nv
     real(dp)     :: nu, x_min, x_max, y_min, y_max, n_inv
     real(dp)     :: sigma0, alpha, fknee, t1, t2
     integer(i4b) :: i, j, k, l, m, n, p, q, scan, buffer, long_samp
     !call wall_time(t1)
     !Nv = 0.d0
     rhs = 0.d0
     buffer = 0!50
     do scan = 1, size(tod)
        long_samp = 2*(tod(scan)%nsamp - 2*buffer)
        n = long_samp/2 + 1!tod(scan)%nsamp/2 + 1
        allocate(Nv(0:n-1), Nt(long_samp))!(tod(scan)%nsamp))
        allocate(dv(0:n-1), dt(long_samp))!(tod(scan)%nsamp))
        
        ! P_T F_inv N_inv F g d
        sigma0 = tod(scan)%sigma0(freq,sb,det)
        alpha  = tod(scan)%alpha(freq,sb,det)
        fknee  = tod(scan)%fknee(freq,sb,det)
        !write (*,*) sigma0, alpha, fknee
        
        !!$OMP PARALLEL PRIVATE(i)
        !!$OMP DO SCHEDULE(guided)
        do i = 1, tod(scan)%nsamp
           dt(i) = tod(scan)%d(i,freq,sb,det) 
           dt(long_samp-i+1) = tod(scan)%d(i,freq,sb,det) 
        end do
        !!$OMP END DO
        !!$OMP END PARALLEL
                
        call fft(dt, dv, 1)
        Nv(0) = 0.d0
        !open(69,file='invN.dat')
        !!$OMP PARALLEL PRIVATE(j,nu,n_inv)
        !!$OMP DO SCHEDULE(guided)
        do j = 1, n-1
           nu = ind2freq(j+1, tod(scan)%samprate, n)
           n_inv = 1.d0/(sigma0**2 * (1.d0 + (nu/fknee)**alpha))
           !if (nu < 0.5d0) then
           !n_inv = 1.d0/(1.d0 + (nu/tod(scan)%fknee(freq,sb,det))**tod(scan)%alpha(freq,sb,det))
           !n_inv = 1.d0/(1.d0 + (nu/fknee)**alpha)
           !else
           !   n_inv = 1.
           !end if
           !n_inv = 1.d0/(tod(scan)%sigma0(freq,sb,det)**2)
           Nv(j) = dv(j)*n_inv
           !write(69,*) nu, n_inv
        end do
        !!$OMP END DO
        !!$OMP END PARALLEL
        !close(69)
     
        call fft(Nt, Nv, -1)

        Nt = Nt - mean(Nt)

        !open(53,file='nt.dat',recl=1024)
        !do i = 1, tod(scan)%nsamp
        !   write(53,*) i, Nt(i)
        !end do
        !close(53)

        !open(58+scan,file='tod'//trim(itoa(scan))//'.dat')
        !do i = 1, tod(scan)%nsamp
        !   write(58+scan,*) i, tod(scan)%d(i,freq,sb,det), Nt(i)
        !end do
        !close(58+scan)
        !call mpi_finalize(i)
        !stop

        do i = 1 + buffer, tod(scan)%nsamp - buffer
           !rhs(tod(scan)%pixels(i,det)) = rhs(tod(scan)%pixels(i,det)) + tod(scan)%g(1,freq,sb,det)*Nt(i)
           rhs(tod(scan)%pixels(i,det)) = rhs(tod(scan)%pixels(i,det)) + Nt(i)
        end do

        !open(75,file='rhs.dat',recl=1024)
        !do i = 1, size(rhs)
        !   write(75,*) i, rhs(i)
        !end do
        !close(75)
        !call mpi_finalize(i)
        !stop

        deallocate(Nv, Nt)
        deallocate(dt, dv)
     

     end do

     !call wall_time(t2)
     !write(*,*) 'wall time get_rhs = ', t2-t1
    
     !call mpi_finalize(i)
     !stop


   end subroutine get_rhs


   ! Computes Ax where A = P^T N^{-1} P
   subroutine get_lhs(Ax, x, map, tod, det, sb, freq)
     implicit none
     real(dp),     dimension(:), intent(inout) :: Ax
     real(dp),     dimension(:), intent(in)    :: x
     type(tod_type),  dimension(:), intent(in) :: tod
     type(map_type),  intent(in) :: map
     integer(i4b),    intent(in) :: det, sb, freq

     integer(i4b) :: i, j, k, l, n, scan, buffer, long_samp!, p, q
     real(dp)     :: nu, sigma0, alpha, fknee, t1, t2!, x_min, x_max, y_min, y_max
     real(dp),     allocatable, dimension(:) :: xt, Nt ! intermediate stages
     complex(dpc), allocatable, dimension(:) :: xv, Nv
     

     buffer = 0!50
     Ax = 0.d0
     do scan = 1, size(tod)
        long_samp = 2*tod(scan)%nsamp
        sigma0 = tod(scan)%sigma0(freq,sb,det)
        alpha  = tod(scan)%alpha(freq,sb,det)
        fknee  = tod(scan)%fknee(freq,sb,det)

        ! long_samp or tod(scan)%nsamp
        n = long_samp/2 + 1
        allocate(xt(long_samp), xv(0:n-1))
        allocate(Nt(long_samp), Nv(0:n-1))
        xt = 0.d0
        !$OMP PARALLEL PRIVATE(i)
        !$OMP DO SCHEDULE(guided)
        do i = 1 + buffer, tod(scan)%nsamp - buffer
           !xt(i) = tod(scan)%g(1,freq,sb,det)*x(tod(scan)%pixels(i,det))
           !xt(long_samp-i+1) = tod(scan)%g(1,freq,sb,det)*x(tod(scan)%pixels(i,det))
           xt(i) = x(tod(scan)%pixels(i,det))
           xt(long_samp-i+1) = x(tod(scan)%pixels(i,det))
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        call fft(xt, xv, 1)
        Nv(0) = 0.d0!xv(1)

        !$OMP PARALLEL PRIVATE(j, nu)
        !$OMP DO SCHEDULE(guided)
        do j = 1, n-1
           nu = ind2freq(j+1, tod(scan)%samprate, n)
           Nv(j) = xv(j)/(sigma0**2 * (1.d0 + (nu/fknee)**alpha))
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        call fft(Nt, Nv, -1)
        
        Nt = Nt - mean(Nt)

        do i = 1 + buffer, tod(scan)%nsamp - buffer
           !Ax(tod(scan)%pixels(i,det)) = Ax(tod(scan)%pixels(i,det)) + tod(scan)%g(1,freq,sb,det)*Nt(i)
           Ax(tod(scan)%pixels(i,det)) = Ax(tod(scan)%pixels(i,det)) + Nt(i)
        end do

        deallocate(xt, xv)
        deallocate(Nt, Nv)
     end do

   end subroutine get_lhs


  subroutine pcg_mapmaker(tod, map, alist, det, sb, freq, parfile)
    implicit none
    type(tod_type),   dimension(:), intent(in) :: tod
    type(map_type),   intent(inout) :: map
    type(acceptlist), intent(in)    :: alist
    integer(i4b),     intent(in)    :: det, sb, freq
    character(len=*)                :: parfile

    real(dp),     allocatable, dimension(:) :: r, q, s, d, b, mp, Ax, Ad, debug, M_inv
    real(dp)     :: delta_new, delta_old, delta_0,  epsilon2, alpha, beta, t1, t2
    integer(i4b) :: i, j, k, imax, npix, nx, ny

    !if (.not. allocated(map%m)) then
    !   call initialize_mapmaker(map, tod, myid)
    !end if

    !if (sb==1 .and. freq==1) write(*,*) "Enter PCG"

    call get_parameter(0, parfile, 'CG_LIM', par_dp=epsilon2)

    npix = map%n_x*map%n_y
    allocate(mp(npix), b(npix), r(npix), s(npix), d(npix), q(npix))
    allocate(Ax(npix), Ad(npix), M_inv(npix))
    
    write(*,*) "getting rhs" 
    call get_rhs(b, map, tod, det, sb, freq)
    mp = 0.d0
    
    ! Preconditioner matrix
    !call precondition_matrix(tod, map, M_inv, det, sb, freq)
    M_inv = 1.d0
    Ax = 0.d0
    !call get_lhs(Ax, mp, map, tod, det, sb, freq)
    r = b - Ax
    d = M_inv*r
    delta_new = sum(r*d)
    delta_0 = delta_new
    !write(*,*) delta_new
    imax = 1.d3
    
    write(*,*) delta_new, epsilon2 * delta_0

    write(*,*) "starting iteration"
    i = 0
    do while (i < imax .and. delta_new > epsilon2 * delta_0)
       !write(*,*) i, delta_new, epsilon2 * delta_0
       !d = 0.d0
       !d(4000) = 1.d0
       call get_lhs(Ad, d, map, tod, det, sb, freq)
       !write(*,*) Ad(4000), sum(abs(Ad))
       !d = 0.d0
       !d(4000) = 1.d0
       !d = d*M_inv
       !write(*,*) d(4000), sum(abs(d))
       !write(*,*)
       !call mpi_finalize(i)
       !stop
       

       q = Ad
       alpha = delta_new/(sum(d*q))
       ! Update map
       mp = mp + alpha*d
       r = r - alpha*q
       s =  M_inv*r
       delta_old = delta_new
       delta_new = sum(r*s)
       beta = delta_new/delta_old
       d = s + beta*d
       i = i + 1
    end do

    
    ! Convert to map
    do j = 1, map%n_x
       do k = 1, map%n_y
          map%m(j,k,freq,sb) = map%m(j,k,freq,sb) + mp((k-1)*map%n_x + j)
       end do
    end do
    
    ! allocate(debug(tod(1)%nsamp))
    ! ! Debugging
    ! call wall_time(t1)
    ! !$OMP PARALLEL PRIVATE(i)
    ! !$OMP DO SCHEDULE(guided)
    ! do i = 1, tod(1)%nsamp
    !    debug(i) = mp(tod(1)%pixels(i,det))
    ! end do
    ! !$OMP END DO
    ! !$OMP END PARALLEL
    ! call wall_time(t2)
    ! write(*,*) 'wall time debug = ', t2-t1
    ! open(75,file='map_time.dat',recl=1024)
    ! do i = 1, tod(1)%nsamp
    !    write(75,*) i, debug(i), tod(1)%pixels(i,det), tod(1)%point_tel(1,i,det), tod(1)%point_tel(2,i,det)
    ! end do
    ! close(75)
    ! call mpi_finalize(i)
    ! stop

    deallocate(mp, Ax, Ad)
    deallocate(b, r, s, d, q)

    !write(*,*) "Exit PCG"

  end subroutine pcg_mapmaker


  subroutine precondition_matrix(tod,map,M_inv,det,sb,freq)
    implicit none
    type(tod_type), dimension(:), intent(in)    :: tod
    type(map_type),               intent(in)    :: map
    real(dp),       dimension(:), intent(inout) :: M_inv
    integer(i4b),                 intent(in)    :: det, sb, freq

    real(dp), allocatable, dimension(:) :: nobs
    integer(i4b) :: i, j, k, scan, tot_scans, npix
    real(dp)     :: sigma0, gain

    npix = map%n_x*map%n_y
    tot_scans = size(tod)

    allocate(nobs(npix))
    sigma0 = 0.d0
    gain = 0.d0
    do scan = 1, tot_scans
       sigma0 = sigma0 + tod(scan)%sigma0(freq,sb,det)**2
       !gain = gain + tod(scan)%g(1,freq,sb,det)
       gain = gain + 1.d0
    end do
    sigma0 = sigma0 / tot_scans
    gain = gain / tot_scans

    nobs = 0.d0
    do j = 1, map%n_x
       do k = 1, map%n_y
          i = (k-1)*map%n_x + j
          nobs(i) = map%nhit(j,k,freq,sb)
       end do
    end do
    !write(*,*) nobs(4000), sigma0
    where(nobs > 0)
       M_inv = 1.d0 / (gain**2 * nobs / sigma0)
    elsewhere
       M_inv = 0.d0
    end where
    !call mpi_finalize(i)
    !stop

  end subroutine precondition_matrix








  !!!!!!!!!!!!!!!!!!
  ! Gibbs sampling !
  !!!!!!!!!!!!!!!!!!

  subroutine gibbs_map(map)
    implicit none
    type(map_type), intent(inout) :: map


  end subroutine gibbs_map








end module
