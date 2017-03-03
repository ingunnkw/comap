program maptest
  use quiet_utils
  use quiet_filter_mod
  use tod2map_utils
  use tod2map_mapmaker
  use math_tools
  implicit none
  type(common_info)    :: info
  type(quiet_assembly) :: assembly
  type(lx_struct)      :: l3file
  type(mapinfo)        :: data
  type(mapdata)        :: out
  integer, allocatable :: map2mask(:)
  logical              :: ityp(itype_num)
  integer              :: i, j, k, l, npad, status
  real*8               :: A(5,2), mats(10,10,4), imats(10,10,4), m(10), eig(10), asym
  real*8, dimension(5,2,5,2,4) :: rhs
  real*8, allocatable  :: tod_pad(:,:)
  npad = 200*0

  info%myid     = 0
  info%nproc    = 1
  info%cid      = 1
  info%nside    = 32
  info%coord    = coord_gal
  info%myid_cov = 0
  info%nproc_cov= 1
  allocate(info%comps(2))
  info%comps    = [2,3]
  info%ncomp    = size(info%comps)
  call make_test_detector
  rhs = 0
  do i = 1, 5
     do j = 1, 2
        call get_test_l3(l3file,[i,j])
        call get_test_assembly(assembly, info, l3file)
        call get_map2mask(assembly, map2mask)
        ityp = .false.
        ityp([itype_cov,itype_rhs]+1) = .true.
        call setup_mapdata(out, count(map2mask>0), 1, info, ityp)
        call initialize_map_debug
        call init_map_maker_bench(info)
        call setup_mapinfo(assembly, info, map2mask, data, out)

        allocate(tod_pad(size(assembly%tod,1)+2*npad,size(assembly%tod,2)))
        ! Build rhs using normal method, but with padding
        tod_pad = 0
        tod_pad(npad+1:size(tod_pad,1)-npad,:) = assembly%tod
        call apply_filter_fft_matrix(assembly%N_corr_fft, tod_pad, tod_pad)
        data%N_d = tod_pad(npad+1:size(tod_pad,1)-npad,:)
        deallocate(tod_pad)
        call build_rhs(assembly, info, data, data%N_d, rhs(:,:,i,j,1))
        ! Build rhs brute force
        call build_rhs_debug(assembly, info, data, assembly%tod, rhs(:,:,i,j,2))
        ! Build rhs with brute force filtering but the normal method
        call apply_filter_real_matrix(assembly%N_corr_real, assembly%tod, data%N_d)
        call build_rhs(assembly, info, data, data%N_d, rhs(:,:,i,j,3))
        call build_cov(assembly, info, data, out)

        call dump_matrix(reshape(out%cov,[size(out%cov(:,:,1,1)),size(out%cov(1,1,:,:))]),"foocov.txt")

        write(*,*) "rhs"
        call dump_matrix(real(rhs(:,:,i,j,1),dp),fmt="(e15.7)")
        write(*,*) "rhs2"
        call dump_matrix(real(rhs(:,:,i,j,2),dp),fmt="(e15.7)")
        write(*,*) "rhs3"
        call dump_matrix(real(rhs(:,:,i,j,3),dp),fmt="(e15.7)")
        write(*,*) "A", shape(out%cov)
        do k = 1, out%n
           if(k > i) then
              A(k,:) = out%cov(i,j,k,:)
           else
              A(k,:) = out%cov(k,:,i,j)
           end if
        end do
        call dump_matrix(A,fmt="(e15.7)")
        rhs(:,:,i,j,4) = A

        call free_lx_struct(l3file)
        call deallocate_quiet_assembly(assembly)
        call free_mapdata(out)
        call free_mapinfo(data)
        deallocate(map2mask)
     end do
  end do

  ! Build our full matrices
  mats = reshape(rhs, shape(mats))

  do i = 1, 4
     call dump_matrix(mats(:,:,i),"foo"//trim(itoa(i))//".txt")
  end do

  ! Use every matric to solve every rhs
  do i = 1, size(mats,3)
     call get_eigenvalues(mats(:,:,i), eig)
     asym = get_asym(mats(:,:,i))
     write(*,'(a,e15.7)')  "A", asym
     write(*,'(a,10f8.4)') "E", eig/eig(size(eig))
     do j = 1, size(mats,3)
        do k = 1, size(rhs,2)
           do l = 1, size(rhs,1)
              call solve_linear_system(mats(:,:,i), m, mats(:,(k-1)*size(rhs,1)+l,j))
              write(*,'(4i2,10f8.4)') i, j, k, l, m
           end do
        end do
        write(*,*)
     end do
  end do

contains

  subroutine make_test_detector
    implicit none
    integer(i4b), parameter :: nmod = 1, ndi = 1
    allocate(quiet_horns(0:0), quiet_diodes(2))
    allocate(quiet_horns(0)%groups(1))
    allocate(quiet_horns(0)%diodes(2))
    allocate(quiet_horns(0)%amps(1))
    quiet_horns(0)%id    = 0
    quiet_horns(0)%board = 0
    quiet_horns(0)%slot  = 0
    quiet_horns(0)%n     = 1
    quiet_horns(0)%ndi   = 1
    quiet_horns(0)%phi   = 0
    quiet_horns(0)%theta = 0
    quiet_horns(0)%fwhm  = 60
    quiet_horns(0)%groups= [0]
    quiet_horns(0)%diodes= [1]
    quiet_horns(0)%amps  = [1]
    quiet_diodes%id      = [1,2]
    quiet_diodes%horn    = 0
    quiet_diodes%sub     = [0,1]
    quiet_diodes%ok      = .true.
    quiet_diodes%psi     = [0.0,0.5]
    quiet_diodes%freq    = 100
    quiet_diodes(1)%stokes  = [0,1,0]
    quiet_diodes(2)%stokes  = [0,1,0]
  end subroutine

  subroutine get_test_l3(l3, pos)
    implicit none
    type(lx_struct)         :: l3
    integer(i4b)            :: pos(2)
    integer(i4b), parameter :: n = 1000, nside = 32, nmod = 1, ndi = 2, ncorr = 1000
    integer(i4b), parameter :: pix_raw(5) = [3550,3423,3424,3296,3295]
    integer(i4b), parameter :: pixorder(5) = [5,4,2,3,1]
    integer(i4b)            :: QU_raw(n,2)
    integer(i4b)            :: ipix(n)
    integer(i4b)            :: pix(n)
    real(dp)                :: psi(n)
    real(dp)                :: Q(n), U(n)
    real(dp)                :: f, theta, phi, c
    integer(i4b)            :: i, j
    QU_raw = 0
    QU_raw(pixorder(pos(1)),pos(2)) = 1
    do i = 1, n
    !   j = modulo(i-1,10)+1
    !   if(j <= 5) then
    !      ipix(i) = j
    !   else
    !      ipix(i) = 11-j
    !   end if
        psi(i) = modulo(i*pi/2,2*pi)
    end do
    ipix       = 5
    ipix(:n/2) = 4
    ipix(:n/3) = 3
    ipix(:n/4) = 2
    ipix(:n/5) = 1
!    psi = 0.5


    Q = QU_raw(ipix,1)
    U = QU_raw(ipix,2)
    pix = pix_raw(ipix)

    l3%decimation = 1
    l3%samprate   = 100.0
    l3%nside      = nside
    l3%coord_sys  = coord_gal
    l3%scanfreq   = 1
    allocate(l3%time(n), l3%tod(n,ndi), l3%point(3,n,nmod), l3%orig_point(3,n))
    l3%time = irange(n)/24.0/60/60/l3%samprate
    l3%point(3,:,1) = psi
    l3%orig_point   = 0
    do i = 1, n
       call pix2ang(nside, nest, pix(i), theta, phi)
       l3%point(1:2,i,1) = [phi,theta]
    end do

    allocate(l3%time_gain(n), l3%gain(n,ndi))
    l3%time_gain = l3%time
    l3%gain      = 1
    allocate(l3%sigma0(ndi),l3%alpha(ndi),l3%fknee(ndi),l3%corr(ncorr,ndi,ndi),l3%corr_freqs(ncorr))
    l3%sigma0 = [1,20]
    l3%fknee  = [1,2]
    l3%alpha  = [-2,-2]
    do i = 1, ncorr
       l3%corr_freqs(i) = i*l3%samprate/2/ncorr
       c = 0.7!0.9/(1+(l3%corr_freqs(i)/8)**-10)
       l3%corr(i,:,:) = reshape([1d0,c,c,1d0],[2,2])
    end do

    ! The rest isn't necessary
    ! Given this, we can compute what our input tod should be.
    do i = 1, ndi
       f = 1d-9/ant2thermo(quiet_diodes(i)%freq)*quiet_diodes(i)%stokes(2)
       l3%tod(:,i) = (cos(2*(psi+quiet_diodes(i)%psi))*Q+sin(2*(psi+quiet_diodes(i)%psi))*U)*l3%gain(:,i)*f
    end do

    allocate(l3%filter_par(ndi,NUM_FILTER_PAR))
    l3%filter_par = 0
  end subroutine

  subroutine get_test_assembly(a, info, l3)
    implicit none
    type(quiet_assembly) :: a
    type(common_info)    :: info
    type(lx_struct)      :: l3
    type(filter_params)  :: opts
    integer(i4b)         :: n, i, j, k, l, nc, m
    integer(i4b), parameter   :: ndi = 2, nmod = 1
    real(dp),     allocatable :: filter(:), tmp1(:), tmp2(:), nfilter(:,:), fcov(:,:)
    n = size(l3%time)
    m = n + npad*2
    a%num_diodes = ndi
    allocate(a%diodes(ndi,2),a%diodes_abs(ndi))
    allocate(filter(m/2+1))
    a%diodes(1,:) = [0,0]
    a%diodes(2,:) = [0,1]
    a%diodes_abs  = [1,2]

    call get_assembly_data(a, info, l3)
    opts%apply_lowpass      = .true.
    opts%apply_highpass     = .false.
    opts%ignore_oneoverf    = .true.
    opts%ignore_diode_corr  = .true.
    opts%use_precomp_filter = .false.
    opts%nu_highpass        = 5d0
    opts%alpha_highpass     = -5d0
    opts%nu_lowpass         = 1d0
    opts%alpha_lowpass      = -5

    call initialize_filter_assembly(a, opts)
    write(*,*) "corrlen", a%noise_corr_length
    return

    call get_inv_noise_filter_fft(.true., m, 100d0, 1d0, 5d0, -5d0, filter)
    !filter = 1

    allocate(a%N_corr_fft(ndi,ndi,m/2+1),a%N_corr_F_sq_fft(ndi,ndi,m/2+1))
    allocate(nfilter(m/2+1,ndi),fcov(ndi,ndi))
    do i = 1, ndi
       call get_inv_noise_filter_fft(.false., m, a%samprate, assembly%sigma0(i), &
        & assembly%fknee(i), assembly%alpha(i), nfilter(:,i))
    end do
    !nfilter=1

    do i = 1, size(filter)
       j = min((i-1)*size(a%corr,1)/size(filter)+1,size(a%corr,1))
       fcov = a%corr(j,:,:)
       do j = 1, ndi
          do k = 1, ndi
             fcov(k,j) = fcov(k,j)*nfilter(i,k)*nfilter(i,j)
          end do
       end do
       call eigen_pow(fcov, -1d0, a%N_corr_fft(:,:,i))
       a%N_corr_fft     (:,:,i) = a%N_corr_fft(:,:,i)*filter(i)
       a%N_corr_F_sq_fft(:,:,i) = a%N_corr_fft(:,:,i)*filter(i)
    end do
    deallocate(filter, nfilter, fcov)

!call dump_matrix(transpose(reshape(a%N_corr_fft,[ndi*ndi,m/2+1])),"ff.txt",fmt="(e25.16)")

    nc = 0
    do i = 1, ndi
       do j = 1, ndi
          call get_inv_noise_filter_real(m, a%N_corr_fft(i,j,:),      k, tmp1)
          call get_inv_noise_filter_real(m, a%N_corr_F_sq_fft(i,j,:), l, tmp2)
          nc = max(nc,max(k,l))
write(*,*) "A", i, j, k, l, nc
          deallocate(tmp1, tmp2)
       end do
    end do
    allocate(a%N_corr_real(-nc:nc,ndi,ndi), a%N_corr_F_sq_real(-nc:nc,ndi,ndi))
    a%N_corr_real = 0
    a%N_corr_F_sq_real = 0
    do i = 1, ndi
       do j = 1, ndi
          call get_inv_noise_filter_real(m, a%N_corr_fft(i,j,:),      k, tmp1)
          call get_inv_noise_filter_real(m, a%N_corr_F_sq_fft(i,j,:), l, tmp2)
          a%N_corr_real(-k:k,i,j) = tmp1(-k:k)
          a%N_corr_F_sq_real(-l:l,i,j) = tmp2(-l:l)
          deallocate(tmp1, tmp2)
       end do
    end do
!call dump_matrix(reshape(a%N_corr_real,[2*k+1,ndi*ndi]),"fr.txt",fmt="(e25.16)")

!a%N_corr_real(-1:1,1,1) = [0.5,1,0.5]
!a%N_corr_F_sq_real(-1:1,1,1) = [0.5,1,0.5]
    a%noise_corr_length = nc
write(*,*) nc
  end subroutine

  subroutine get_map2mask(a, m2m)
    implicit none
    type(quiet_assembly) :: a
    integer, allocatable :: m2m(:), hits(:)
    integer              :: i, j, k
    allocate(m2m(0:12*a%nside**2-1), hits(0:12*a%nside**2-1))
    hits = 0
    do k = 1, size(a%diode_info)
       do i = 1, size(a%diode_info(k)%pix,2)
          do j = 1, size(a%diode_info(k)%pix,1)
             hits(a%diode_info(k)%pix(j,i)) = hits(a%diode_info(k)%pix(j,i)) + 1
          end do
       end do
    end do
    j = 0
    m2m = 0
    do i = 0, size(hits)-1
       if(hits(i) == 0) cycle
       j = j+1
       m2m(i) = j
    end do
    deallocate(hits)
  end subroutine

  function get_asym(mat) result(maxasym)
    implicit none
    real(dp)     :: mat(:,:), maxasym, asym
    integer(i4b) :: i, j, n
    n = size(mat,1)
    maxasym = 0
    do i = 1, n
       do j = i+1, n
          asym = abs(mat(j,i)-mat(i,j))/sqrt(abs(mat(i,i)*mat(j,j)))
          write(*,'(2i4,5e15.7)') i, j, asym, mat(i,i), mat(j,j), mat(i,j), mat(j,i)
          maxasym = max(asym,maxasym)
       end do
    end do
  end function

end program
