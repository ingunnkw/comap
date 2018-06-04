module tod2comap_mapmaker
  !use comap_lx_mod
  use comap_map_mod
  use comap_acceptlist_mod
  !use quiet_fft_mod
  use tod2comap_utils
  implicit none



contains


  subroutine initialize_mapmaker(map, tod)
    implicit none
    type(tod_type),   intent(in)    :: tod
    type(map_type),   intent(inout) :: map

    integer(i4b) :: i, j, k, l, p, q, fs, st
    real(dp)     :: x_min, x_max, y_min, y_max, pad
    real(8), parameter :: PI = 4*atan(1.d0)

    ! Set up map grid
    if (.not. allocated(map%x)) then
       fs = 200
       st = tod%nsamp-200 ! tod%nsamp
       pad = 0.5d0 ! degrees
       map%nfreq = tod%nfreq
       map%nsb   = tod%nsb
       map%dthetay = 0.5d0/60.d0 ! degrees (arcmin/60), resolution
       map%dthetax = map%dthetay/abs(cos(tod%mean_el*PI/180.d0))

       x_min = minval(tod%point(1,fs:st)) - pad; x_max =  maxval(tod%point(1,fs:st)) + pad
       y_min = minval(tod%point(2,fs:st)) - pad; y_max =  maxval(tod%point(2,fs:st)) + pad
       !x_min = data%point_lim(1) - pad; x_max = data%point_lim(2) - pad
       !y_min = data%point_lim(3) - pad; y_max = data%point_lim(4) - pad
       map%n_x = (x_max-x_min)/map%dthetax+1; map%n_y = (y_max-y_min)/map%dthetay+1
       allocate(map%x(map%n_x), map%y(map%n_y))
       do i = 1, map%n_x
          map%x(i) = x_min + (i-1)*map%dthetax
       end do
       do i = 1, map%n_y
          map%y(i) = y_min + (i-1)*map%dthetay
       end do
    end if

    ! Set up map structures
    if (.not. allocated(map%dsum)) then
       allocate(map%m(map%n_x, map%n_y, map%nfreq, map%nsb), &
            & map%dsum(map%n_x, map%n_y, map%nfreq, map%nsb), &
            & map%nhit(map%n_x, map%n_y, map%nfreq, map%nsb), &
            & map%div(map%n_x, map%n_y, map%nfreq, map%nsb), &
            & map%rms(map%n_x, map%n_y, map%nfreq, map%nsb))
       map%dsum = 0.d0
       map%nhit = 0.d0
       map%div  = 0.d0
       map%m    = 0.d0
    end if

  end subroutine initialize_mapmaker


  subroutine time2pix(tod, map)
    implicit none
    type(tod_type), intent(inout) :: tod
    type(map_type), intent(inout) :: map

    integer(i4b) :: i, p, q
    real(dp)     :: x_min, x_max, y_min, y_max

    x_min = map%x(1); x_max = map%x(map%n_x)
    y_min = map%y(1); y_max = map%y(map%n_y)

    do i = 1, tod%nsamp
       p = min(max(int((tod%point(1,i)-x_min)/map%dthetax),1),map%n_x)
       q = min(max(int((tod%point(2,i)-y_min)/map%dthetay),1),map%n_y)
       tod%pixels(i) = (q-1)*map%n_x + p
       map%nhit(p,q,:,:) = map%nhit(p,q,:,:) + 1.d0
    end do



  end subroutine time2pix


  
  !!!!!!!!!!!!!!!!!
  ! Naive Binning !
  !!!!!!!!!!!!!!!!!

  subroutine binning(map, tod, alist)
    implicit none
    type(tod_type),   intent(in)    :: tod
    type(map_type),   intent(inout) :: map
    type(acceptlist), intent(in)    :: alist

    integer(i4b) :: i, j, k, l, p, q, fs, st
    real(dp)     :: x_min, x_max, y_min, y_max

    fs = 200 ! starting point
    st = tod%nsamp - 200 ! ending point

    x_min = map%x(1); x_max = map%x(map%n_x)
    y_min = map%y(1); y_max = map%y(map%n_y)

    write(*,*) 'Beginning coadding'

    k = 1 ! detector number

    do i = fs, st!tod%nsamp
       p = min(max(int((tod%point(1,i)-x_min)/map%dthetax),1),map%n_x)
       q = min(max(int((tod%point(2,i)-y_min)/map%dthetay),1),map%n_y)
       do l = 1, tod%nsb
          !do j = 6, 6
          do j = 1, tod%nfreq
             if (alist%status(j,k) == 0) then
                if (tod%g(1,j,l,k) .ne. 0.d0) then

                   map%dsum(p,q,j,l) = map%dsum(p,q,j,l) + tod%g(1,j,l,k)    / tod%rms(i,j,l,k)**2 * tod%d(i,j,l,k)
                   map%div(p,q,j,l)  = map%div(p,q,j,l)  + tod%g(1,j,l,k)**2 / tod%rms(i,j,l,k)**2
                   map%nhit(p,q,j,l) = map%nhit(p,q,j,l) + 1.d0

                end if
             end if
          end do
       end do
    end do

    write(*,*) 'Ending coadding'

  end subroutine binning

  subroutine finalize_binning(map)
    implicit none
    type(map_type), intent(inout) :: map

    where(map%nhit > 0)
       map%m   = map%dsum / map%div
       map%rms = 1.d0 / sqrt(map%div)
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
     real(dp)     :: sigma0, alpha, fknee
     integer(i4b) :: i, j, k, l, m, n, p, q, scan, buffer, long_samp

     !Nv = 0.d0
     rhs = 0.d0
     buffer = 0!50
     do scan = 1, size(tod)
        long_samp = 2*tod(scan)%nsamp
        n = long_samp/2 + 1!tod(scan)%nsamp/2 + 1
        allocate(Nv(0:n-1), Nt(long_samp))!(tod(scan)%nsamp))
        allocate(dv(0:n-1), dt(long_samp))!(tod(scan)%nsamp))
        
        ! P_T F_inv N_inv F g d
        sigma0 = tod(scan)%sigma0(freq,sb,det)
        alpha  = -1.5d0!tod(scan)%alpha(freq,sb,det)
        fknee  = 1.d0!tod(scan)%fknee(freq,sb,det)
        !write(*,*) tod(scan)%samprate
        !write (*,*) sigma0, alpha, fknee
        !call mpi_finalize(j)
        !stop
        do i = 1, tod(scan)%nsamp
           dt(i) = tod(scan)%d(i,freq,sb,det) 
           dt(long_samp-i+1) = tod(scan)%d(i,freq,sb,det) 
        end do
        call fft(dt, dv, 1)
        Nv(0) = 0.d0
        !open(69,file='invN.dat')
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
        !close(69)
        call fft(Nt, Nv, -1)

        !open(58,file='tod.dat')
        !do i = 1, tod(scan)%nsamp
        !   write(58,*) i, tod(scan)%d(i,freq,sb,det), Nt(i)
        !end do
        !close(58)
        !call mpi_finalize(i)
        !stop


        do i = 1 + buffer, tod(scan)%nsamp - buffer
           rhs(tod(scan)%pixels(i)) = rhs(tod(scan)%pixels(i)) + tod(scan)%g(1,freq,sb,det)*Nt(i)
        end do
        !rhs(tod(scan)%pixels) = rhs(tod(scan)%pixels) + tod(scan)%g(1,freq,sb,det)*Nt ! does this work?
        
        deallocate(Nv, Nt)
        deallocate(dt, dv)
     end do


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
     real(dp)     :: nu, sigma0, alpha, fknee!, x_min, x_max, y_min, y_max
     real(dp),     allocatable, dimension(:) :: xt, Nt ! intermediate stages
     complex(dpc), allocatable, dimension(:) :: xv, Nv
     

     buffer = 0!50
     Ax = 0.d0
     do scan = 1, size(tod)
        long_samp = 2*tod(scan)%nsamp
        sigma0 = tod(scan)%sigma0(freq,sb,det)
        alpha  = -1.5d0!tod(scan)%alpha(freq,sb,det)
        fknee  = 1.d0!tod(scan)%fknee(freq,sb,det)

        n = long_samp/2 + 1!tod(scan)%nsamp/2 + 1
        allocate(xt(long_samp), xv(0:n-1))!(xt(tod(scan)%nsamp), xv(0:n-1))
        allocate(Nt(long_samp), Nv(0:n-1))!(Nt(tod(scan)%nsamp), Nv(0:n-1))
        xt = 0.d0
        do i = 1 + buffer, tod(scan)%nsamp - buffer
           xt(i) = xt(i) + tod(scan)%g(1,freq,sb,det)*x(tod(scan)%pixels(i))
           xt(long_samp-i+1) = xt(long_samp-i+1) + tod(scan)%g(1,freq,sb,det)*x(tod(scan)%pixels(i))
        end do
        call fft(xt, xv, 1)
        Nv(0) = 0.d0!xv(1)
        do j = 1, n-1
           nu = ind2freq(j+1, tod(scan)%samprate, n)
           Nv(j) = xv(j)/(sigma0**2 * (1.d0 + (nu/fknee)**alpha))
        end do
        call fft(Nt, Nv, -1)
        do i = 1 + buffer, tod(scan)%nsamp - buffer
           Ax(tod(scan)%pixels(i)) = Ax(tod(scan)%pixels(i)) + tod(scan)%g(1,freq,sb,det)*Nt(i)
        end do

        deallocate(xt, xv)
        deallocate(Nt, Nv)
     end do
  
   end subroutine get_lhs


  subroutine pcg_mapmaker(tod, map, alist, det, sb, freq)
    implicit none
    type(tod_type),   dimension(:), intent(in) :: tod
    type(map_type),   intent(inout) :: map
    type(acceptlist), intent(in)    :: alist
    integer(i4b),     intent(in)    :: det, sb, freq

    real(dp),     allocatable, dimension(:) :: r, q, s, d, b, mp, Ax, Ad
    real(dp)     :: delta_new, delta_old, delta_0,  epsilon, alpha, beta
    real(dp)     :: M_inv
    integer(i4b) :: i, j, k, imax, npix, nx, ny

    !if (.not. allocated(map%m)) then
    !   call initialize_mapmaker(map, tod, myid)
    !end if

    !if (sb==1 .and. freq==1) write(*,*) "Enter PCG"
    
    npix = map%n_x*map%n_y
    allocate(mp(npix), b(npix), r(npix), s(npix), d(npix), q(npix))
    allocate(Ax(npix), Ad(npix))

       
    call get_rhs(b, map, tod, det, sb, freq)
    mp = 0.d0

    ! Preconditioner matrix
    !call precondition_matrix()
    M_inv = 1.d0
    Ax = 0.d0
    !call get_lhs(Ax, mp, map, tod, det, sb, freq)
    r = b - Ax
    d = M_inv*r
    delta_new = sum(r*d)
    delta_0 = delta_new
    
    epsilon = 1.d-5
    imax = 2.d2
    
    i = 0
    do while (i < imax .and. delta_new > epsilon**2 * delta_0)
       call get_lhs(Ad, d, map, tod, det, sb, freq)
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
          
    deallocate(mp, Ax, Ad)
    deallocate(b, r, s, d, q)

    !write(*,*) "Exit PCG"

  end subroutine pcg_mapmaker


  subroutine precondition_matrix()
    implicit none

  end subroutine precondition_matrix








  !!!!!!!!!!!!!!!!!!
  ! Gibbs sampling !
  !!!!!!!!!!!!!!!!!!

  subroutine gibbs_map(map)
    implicit none
    type(map_type), intent(inout) :: map


  end subroutine gibbs_map








end module
