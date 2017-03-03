module map2cl_utils
  use quiet_utils
  use quiet_fileutils
  use quiet_mpi_mod
  use sort_utils
  implicit none

  type common_info
     integer(i4b)       :: myid, myid_group, myid_inter
     integer(i4b)       :: numprocs, numprocs_group, numprocs_inter, comm_group, comm_inter, comm_world, root
  end type common_info

  type powspec
     integer(i4b) :: npar, nbin, num_highl_junk_bins, num_lowl_junk_bins, num_data_set, first_bin_print, last_bin_print
     logical(lgt),              dimension(6)   :: enable_spec
     integer(i4b), allocatable, dimension(:)   :: ind2bin, ind2spec, ind2set
     integer(i4b), allocatable, dimension(:,:) :: bins
     real(dp),     allocatable, dimension(:)   :: coeff
  end type powspec

contains
  
  subroutine allocate_powspec(x, enable_spec, bins, first_bin_print, last_bin_print, &
       & num_lowl_junk_bins, num_highl_junk_bins, num_data_set)
    implicit none

    type(powspec),                   intent(out) :: x
    logical(lgt),  dimension(6),     intent(in)  :: enable_spec
    integer(i4b),  dimension(1:,1:), intent(in)  :: bins
    integer(i4b),                    intent(in)  :: first_bin_print, last_bin_print
    integer(i4b),                    intent(in)  :: num_lowl_junk_bins, num_highl_junk_bins
    integer(i4b),                    intent(in)  :: num_data_set

    integer(i4b) :: npar

    ! Clean up if needed
    call deallocate_powspec(x)

    ! Initialize structures
    x%nbin = size(bins(:,1))
    npar   = get_npar(bins, enable_spec, num_data_set, num_lowl_junk_bins, num_highl_junk_bins)
    x%npar = npar
    
    allocate(x%coeff(npar))
    allocate(x%bins(x%nbin,2))
    allocate(x%ind2bin(npar))
    allocate(x%ind2spec(npar))
    allocate(x%ind2set(npar))

    x%coeff               = 0.d0
    x%bins                = bins
    x%enable_spec         = enable_spec
    x%num_highl_junk_bins = num_highl_junk_bins
    x%num_lowl_junk_bins  = num_lowl_junk_bins
    x%num_data_set        = num_data_set
    x%first_bin_print     = first_bin_print
    x%last_bin_print      = last_bin_print
    call initialize_index_mapping(bins, enable_spec, num_data_set, &
         & num_lowl_junk_bins, num_highl_junk_bins, x%ind2bin, x%ind2spec, x%ind2set)

  end subroutine allocate_powspec


  subroutine copy_powspec(x, y)
    implicit none

    type(powspec),                   intent(in)  :: x
    type(powspec),                   intent(out) :: y

    ! Clean up if needed
    call deallocate_powspec(y)

    allocate(y%coeff(x%npar))
    allocate(y%bins(x%nbin,2))
    allocate(y%ind2bin(x%npar))
    allocate(y%ind2spec(x%npar))
    allocate(y%ind2set(x%npar))
    y%nbin                = x%nbin
    y%npar                = x%npar
    y%coeff               = x%coeff
    y%bins                = x%bins
    y%enable_spec         = x%enable_spec
    y%num_highl_junk_bins = x%num_highl_junk_bins
    y%num_lowl_junk_bins  = x%num_lowl_junk_bins
    y%first_bin_print     = x%first_bin_print
    y%last_bin_print      = x%last_bin_print
    y%num_data_set        = x%num_data_set
    y%ind2bin             = x%ind2bin
    y%ind2spec            = x%ind2spec
    y%ind2set             = x%ind2set

  end subroutine copy_powspec

  subroutine deallocate_powspec(x)
    implicit none

    type(powspec) :: x
    
    if (allocated(x%coeff))     deallocate(x%coeff)
    if (allocated(x%bins))      deallocate(x%bins)
    if (allocated(x%ind2bin))   deallocate(x%ind2bin)
    if (allocated(x%ind2spec))  deallocate(x%ind2spec)
    if (allocated(x%ind2set))   deallocate(x%ind2set)

  end subroutine deallocate_powspec


  subroutine read_bins(filename, bins, lmax)
    implicit none

    character(len=*),                          intent(in)  :: filename
    integer(i4b),     pointer, dimension(:,:)              :: bins
    integer(i4b),                              intent(out) :: lmax

    integer(i4b) :: unit, i, j, k, b_tot, b1, b2, n
    character(len=128) :: t1

    unit = getlun()
    
    open(unit, file=trim(filename))
    read(unit,*) t1
    read(unit,*) b_tot, b1, b2
    n = b2-b1+1
    allocate(bins(n,2))
    do k = 1, b_tot
       read(unit,*) i, j
       if (k < b1) then
          cycle
       else
          bins(k-b1+1,1) = i
          bins(k-b1+1,2) = j
       end if
    end do

    lmax = bins(n,2)

  end subroutine read_bins

  subroutine read_powspec_map2cl(filename, cls)
    implicit none

    character(len=128),                   intent(in)  :: filename
    real(dp),           dimension(0:,1:), intent(out) :: cls

    integer(i4b) :: l, lmax, nmaps
    real(dp)     :: sigma_sq
    real(dp),          allocatable, dimension(:,:)   :: cls_in
    character(len=80),              dimension(1:180) :: header

    lmax  = size(cls(:,1))-1
    nmaps = size(cls(0,:))

    allocate(cls_in(0:lmax,6))

    call fits2cl(filename, cls_in, lmax, 6, header)

    cls(:,1) = cls_in(:,1) ! TT
    cls(:,2) = cls_in(:,4) ! TE
    cls(:,3) = cls_in(:,5) ! TB
    cls(:,4) = cls_in(:,2) ! EE
    cls(:,5) = cls_in(:,6) ! EB
    cls(:,6) = cls_in(:,3) ! BB

    deallocate(cls_in)
    
  end subroutine read_powspec_map2cl


  subroutine compute_rotation_angle(vec1, vec2, t2a_1, t2a_2)
    implicit none

    real(dp), dimension(3),   intent(in)  :: vec1, vec2
    real(dp), dimension(:,:), intent(out) :: t2a_1, t2a_2

    integer(i4b) :: i, j, nfield
    real(dp) :: len_u, len_v, z, cos_theta, sgn
    real(dp), dimension(3) :: u, v
    
    nfield = size(t2a_1(:,1))

    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       do i = 1, nfield
          do j = 1, nfield
             if (i == j) then
                t2a_1(i,j) = 1.d0
                t2a_2(i,j) = 1.d0
             else
                t2a_1(i,j) = 0.d0
                t2a_2(i,j) = 0.d0
             end if
          end do
       end do
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
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    if (nfield == 1) then
       t2a_1(1,1) = 1.d0
    else if (nfield == 2) then
       t2a_1(1,1)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_1(2,2)  = t2a_1(1,1)             ! UU
       t2a_1(2,1)  = sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_1(1,2)  = -t2a_1(2,1)                                    ! QU
    else if (nfield == 3) then
       t2a_1(1,1)   = 1.d0 ! TT
       t2a_1(1,2)   = 0.d0 ! TQ
       t2a_1(1,3)   = 0.d0 ! TU
       t2a_1(2,1)   = 0.d0 ! QT
       t2a_1(3,1)   = 0.d0 ! UT
       t2a_1(2,2)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_1(3,3)  = t2a_1(2,2)             ! UU
       t2a_1(3,2)  = sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_1(2,3)  = -t2a_1(3,2)                                    ! QU
    end if

    ! Rotation from vec2 to vec 1; sgn is opposite from 1->2
    u         = vec2(3) * vec2
    u(3)      = u(3) - 1.d0
    v         = vec1 - z * vec2
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    cos_theta = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    if (nfield == 1) then
       t2a_2(1,1) = 1.d0
    else if (nfield == 2) then
       t2a_2(1,1)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_2(2,2)  = t2a_2(1,1)             ! UU
       t2a_2(2,1)  = -sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_2(1,2)  = -t2a_2(2,1)                                    ! QU
    else if (nfield == 3) then
       t2a_2(1,1)   = 1.d0 ! TT
       t2a_2(1,2)   = 0.d0 ! TQ
       t2a_2(1,3)   = 0.d0 ! TU
       t2a_2(2,1)   = 0.d0 ! QT
       t2a_2(3,1)   = 0.d0 ! UT
       t2a_2(2,2)  = 2.d0*cos_theta**2 -1.d0   ! QQ
       t2a_2(3,3)  = t2a_2(2,2)             ! UU
       t2a_2(3,2)  = -sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta  ! UQ
       t2a_2(2,3)  = -t2a_2(3,2)                                    ! QU
    end if

  end subroutine compute_rotation_angle

!!$  subroutine output_rotation_angles(vec, p0, filename)
!!$    implicit none
!!$
!!$    real(dp),         dimension(:,:), intent(in) :: vec
!!$    integer(i4b),                     intent(in) :: p0
!!$    character(len=*),                 intent(in) :: filename
!!$
!!$    integer(i4b) :: i, npix
!!$    real(dp), allocatable, dimension(:,:) :: map
!!$    real(dp),              dimension(2,2) :: t2a
!!$
!!$    npix = size(vec(1,:))
!!$
!!$    allocate(map(0:npix-1,4))
!!$    map = 0.d0
!!$    do i = 1, npix
!!$       call compute_rotation_angle(vec(:,p0), vec(:,i), t2a)
!!$       map(i-1,1) = 0.5d0*acos(t2a(1,1))
!!$       map(i-1,2) = 0.5d0*asin(t2a(2,1))
!!$       map(i-1,3) = 0.5d0*acos(t2a(1,2))
!!$       map(i-1,4) = 0.5d0*asin(t2a(2,2))
!!$    end do
!!$    call write_map(filename, map, 2)
!!$
!!$    deallocate(map)
!!$
!!$  end subroutine output_rotation_angles


  subroutine output_powspec(filename, append, x, sigma)
    implicit none

    character(len=*),                   intent(in) :: filename
    logical(lgt),                       intent(in) :: append
    type(powspec),                      intent(in) :: x, sigma

    integer(i4b) :: i, j, k, l, m, n, b, nbin, unit, counter, nspec
    character(len=*), parameter  :: fmt1 = '(f10.3,i6,i6)'
    character(len=*), parameter  :: fmt2 = '(f12.4,f12.4)'

    nbin  = size(x%bins(:,1))
    unit  = getlun()
    nspec = count(x%enable_spec)
    m     = size(x%coeff) / nspec
    n     = m - x%num_lowl_junk_bins - x%num_highl_junk_bins

    if (append) then
       open(unit,file=trim(filename), recl=4096, access='append', status='old')       
       write(unit,*)
    else
       open(unit,file=trim(filename), recl=4096)
       write(unit,*) '#  l_bin    l_1   l_2       TT        sigma_TT       TE       sigma_TE      TB       sigma_TB        EE       sigma_EE      EB       sigma_EB       BB       sigma_BB'
    end if

    ! Do all CMB bins
    do i = 1, m
       b = x%ind2bin(i)
       if (b < x%first_bin_print .or. b > x%last_bin_print) cycle
       if (x%ind2set(i) /= 0) write(unit,fmt='(a, i8)',advance='no') '#$', x%ind2set(i)
       write(unit,fmt1,advance='no') 0.5d0*sum(x%bins(b,:)), x%bins(b,:)
       counter = 0
       do j = 1, 6
          if (x%enable_spec(j)) then
             write(unit,fmt2,advance='no') x%coeff(i+counter*m), sigma%coeff(i+counter*m)
             counter = counter+1
          else
             write(unit,fmt2,advance='no') 0.d0, 0.d0
          end if
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine output_powspec

  subroutine compute_legendre_polynomial(z, pls)
    implicit none

    real(dp),                      intent(in)  :: z
    real(dp), dimension(0:,1:,1:), intent(out) :: pls

    integer(i4b) :: i, s_i1, s_i2, s1, s2, l, lmax
    real(dp)     :: rho
    real(dp), dimension(-1:1,-1:1) :: pl_m1, pl_00, pl_p1
    
    lmax = size(pls(:,1,1))-1

    pl_m1 = 0.d0; pl_00 = 0.d0; pl_p1 = 0.d0
    pls   = 0.d0

    ! Initialize recursions
    pl_m1(0,  0) = 1.d0
    pl_00(0,  0) = z
    pl_p1(0,  0) = 0.5d0 * (3.d0*z**2 - 1.d0)
    
    pl_p1( 1,  0) = 0.25d0 * sqrt(6.d0) * (1.d0+z) * (1.d0-z)
    pl_p1( 0,  1) = pl_p1(1, 0)
    pl_p1(-1,  0) = -pl_p1(1, 0)
    pl_p1( 0, -1) = -pl_p1(1, 0)
    
    pl_p1( 1,  1) = 0.25d0 * (1.d0+z)**2 
    pl_p1(-1, -1) = pl_p1(1, 1)
    pl_p1( 1, -1) = 0.25d0 * (1.d0-z)**2 
    pl_p1(-1,  1) = pl_p1(1,-1)
    
    pls(0,:,:) = pl_m1
    pls(1,:,:) = pl_00
    pls(2,:,:) = pl_p1
    
    pl_m1 = pl_00
    pl_00 = pl_p1
                
    ! Do the recursions
    do l = 2, lmax-1
       do s_i1 = -1, 1
          do s_i2 = -1, 1
             
             s1 = 2*s_i1
             s2 = 2*s_i2
             
             rho = sqrt(real(l**2 - s1**2,dp) * real(l**2 - s2**2,dp)) / real(l,dp)
             
             pl_p1(s_i1,s_i2) = real(2*l+1,dp) * (z - real(s1*s2,dp)/real(l*(l+1),dp)) * &
                  & pl_00(s_i1, s_i2) - rho * pl_m1(s_i1, s_i2)
             
             rho = sqrt(real((l+1)**2 - s1**2,dp) * real((l+1)**2 - s2**2,dp)) / &
                  & real(l+1,dp)
             pl_p1(s_i1,s_i2) = pl_p1(s_i1,s_i2) / rho
             
          end do
       end do
          
       pls(l+1,:,:) = pl_p1
                   
       pl_m1 = pl_00
       pl_00 = pl_p1
       
    end do

  end subroutine compute_legendre_polynomial

  subroutine compute_legendre_polynomial_QQ_UU(z, pls)
    implicit none

    real(dp),                      intent(in)  :: z
    real(dp), dimension(0:,1:,1:), intent(out) :: pls

    integer(i4b) :: i, s_i1, s_i2, s1, s2, l, lmax
    real(dp)     :: rho
    real(dp), dimension(-1:1,-1:1) :: pl_m1, pl_00, pl_p1
    
    lmax = size(pls(:,1,1))-1

    pl_m1 = 0.d0; pl_00 = 0.d0; pl_p1 = 0.d0
    pls   = 0.d0

    ! Initialize recursions
    pl_p1( 1,  1) = 0.25d0 * (1.d0+z)**2 
    pl_p1( 1, -1) = 0.25d0 * (1.d0-z)**2 
    
    pls(0,:,:) = pl_m1
    pls(1,:,:) = pl_00
    pls(2,:,:) = pl_p1
    
    pl_m1 = pl_00
    pl_00 = pl_p1
                
    ! Do the recursions
    do l = 2, lmax-1
       rho = sqrt(real(l**2 - 4,dp) * real(l**2 - 4,dp)) / real(l,dp)
       pl_p1(1, 1) = real(2*l+1,dp) * (z - 4.d0/real(l*(l+1),dp)) * &
            & pl_00(1, 1) - rho * pl_m1(1, 1)
       pl_p1(1,-1) = real(2*l+1,dp) * (z + 4.d0/real(l*(l+1),dp)) * &
            & pl_00(1,-1) - rho * pl_m1(1,-1)
       
       rho = sqrt(real((l+1)**2 - 4,dp) * real((l+1)**2 - 4,dp)) / &
            & real(l+1,dp)
       pl_p1(1, 1) = pl_p1(1, 1) / rho       
       pl_p1(1,-1) = pl_p1(1,-1) / rho       

       pls(l+1,:,:) = pl_p1
                   
       pl_m1 = pl_00
       pl_00 = pl_p1
       
    end do

  end subroutine compute_legendre_polynomial_QQ_UU


  subroutine output_mean_powspec(filename, x)
    implicit none

    character(len=*),                 intent(in) :: filename
    type(powspec),    dimension(:),   intent(in) :: x

    integer(i4b) :: i, j, npar, numsim
    type(powspec) :: mu, sigma

    numsim = size(x)
    npar   = x(1)%npar
    
    call copy_powspec(x(1), mu)
    call copy_powspec(x(1), sigma)
    mu%coeff    = 0.d0
    sigma%coeff = 0.d0
    
    ! Compute mean
    do i = 1, numsim
       mu%coeff = mu%coeff + x(i)%coeff
    end do
    mu%coeff = mu%coeff / real(numsim,dp)

    ! Compute standard deviation
    do i = 1, numsim
       sigma%coeff = sigma%coeff + (x(i)%coeff-mu%coeff)**2
    end do
    sigma%coeff = sqrt(sigma%coeff / real(numsim-1,dp))

    ! Output to file
    call output_powspec(filename, .false., mu, sigma)

    call deallocate_powspec(mu)
    call deallocate_powspec(sigma)

  end subroutine output_mean_powspec

  subroutine get_max_dist(vec, cos_theta)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: vec
    real(dp),                 intent(out) :: cos_theta

    integer(i4b) :: i, j, n

    n = size(vec(1,:))
    cos_theta = 1.d30
    do i = 1, n
       do j = i+1, n
          cos_theta = min(cos_theta, sum(vec(:,i)*vec(:,j)))
       end do
    end do

  end subroutine get_max_dist

  subroutine output_statistics(filename, lnL, chisq, np)
    implicit none

    character(len=*),                intent(in) :: filename
    real(dp),         dimension(0:), intent(in) :: lnL, chisq
    integer(i4b),                    intent(in) :: np

    integer(i4b) :: i, n, unit

    unit = getlun()
    n    = size(lnL)-1
    open(unit,file=trim(filename), recl=512)
    do i = 0, n
       write(unit,*) i, real(lnL(i),sp), real(chisq(i),sp), real((chisq(i)-np)/sqrt(2.d0*np),sp)
    end do
    close(unit)

  end subroutine output_statistics


!!$  subroutine output_p_values(filename, x)
!!$    implicit none
!!$
!!$    character(len=*),                   intent(in) :: filename
!!$    real(dp),         dimension(1:,0:), intent(in) :: x
!!$
!!$    integer(i4b) :: i, j, numsim, npar, unit
!!$    real(dp), allocatable, dimension(:,:) :: p
!!$    real(dp), allocatable, dimension(:)   :: d, p_theory, p_sim, p_sort
!!$
!!$    unit   = getlun()
!!$    npar   = size(x,1)
!!$    numsim = size(x,2)-1
!!$    allocate(p(npar,0:numsim), d(0:numsim), p_theory(0:numsim), p_sim(0:numsim), p_sort(npar))
!!$
!!$    ! Compute PTE-values for each bin
!!$    do i = 0, numsim
!!$       do j = 1, npar
!!$          p(j,i) = real(count(x(j,:) > x(j,i)),dp) / real(numsim,dp)
!!$       end do
!!$    end do
!!$
!!$    ! Compute D (= max distance between cumulative distribution and P(x)=x) values for each realization
!!$    D = 0.d0
!!$    do i = 0, numsim
!!$       p_sort = p(:,i)
!!$       call QuickSort_real(p_sort)
!!$       do j = 1, npar
!!$          D(i) = max(D(i), abs(real(j,dp)/real(npar,dp) - p_sort(j)))
!!$          D(i) = max(D(i), abs(real(j-1,dp)/real(npar,dp) - p_sort(j)))
!!$       end do
!!$    end do
!!$
!!$    ! Compute final p-values for each sim
!!$    do i = 0, numsim
!!$       P_sim(i) = real(count(D > D(i)),sp) / real(numsim,dp)
!!$    end do
!!$
!!$    ! Write results to file
!!$    open(unit,file=trim(filename), recl=512)
!!$    do i = 0, numsim
!!$       write(unit,fmt='(a,i5,a,f5.3,a,f5.3)') '# Realization = ', i, ' -- D = ', D(i), ', P_sim = ', P_sim(i)
!!$       do j = 1, npar
!!$          write(unit,*) i, j, real(p(j,i),sp)
!!$       end do
!!$       write(unit,*)
!!$    end do
!!$    close(unit)    
!!$
!!$    deallocate(p, d, p_theory, p_sim, p_sort)
!!$
!!$  end subroutine output_p_values
    
  subroutine input_powspec(cls_fid, cls, realization, filename)
    implicit none

    character(len=*),                  intent(in), optional  :: filename
    integer(i4b),                      intent(in), optional  :: realization
    real(dp),     dimension(0:,1:),    intent(in)            :: cls_fid
    real(dp),     dimension(0:,1:,1:), intent(out)           :: cls

    integer(i4b)        :: unit, lmax, i, l, l1, l2, num_data_set
    real(dp)            :: l_center, EE, sEE, BB, sBB
    real(dp), dimension(2,6) :: spec
    logical(lgt)        :: exist
    character(len=4)    :: real_text
    character(len=1024) :: clfile
    character(len=1024) :: line

    unit         = getlun()
    lmax         = size(cls(:,1,1))-1
    num_data_set = size(cls(0,1,:))
    cls          = 0.d0

    if (present(filename)) then
       clfile = filename
    else if (present(realization)) then
       call int2string(realization, real_text)
       clfile = 'cls_real' // real_text // '.dat'
    end if
    inquire(file=trim(clfile), exist=exist)
    if (exist) then
       write(*,*) 'Reading power spectrum from ', trim(clfile)
       open(unit,file=trim(clfile))
       do while (.true.)
          read(unit,'(a)',end=100) line
          line = trim(adjustl(line))
          if (line(1:2) == '#$') then
             read(line(3:),*) i, l_center, l1, l2, spec
             do l = l1, l2
                if (l <= lmax) cls(l,:,i) = spec(1,:) / (real(l*(l+1),dp) / (2.d0*pi))
             end do
          else if (line(1:1) == '#') then
             cycle
          else
             read(line,*) l_center, l1, l2, spec
             do i = 1, num_data_set
                do l = l1, l2
                   if (l <= lmax) cls(l,:,i) = spec(1,:) / (real(l*(l+1),dp) / (2.d0*pi))
                end do
             end do
          end if
       end do
100    close(unit)
    else
       do i = 1, num_data_set
          cls(:,:,i) = cls_fid
       end do
    end if
    
  end subroutine input_powspec


  subroutine input_powspec_simple(filename, cls, lmin_tab, lmax_tab)
    implicit none

    character(len=*),               intent(in)  :: filename
    real(dp),     dimension(0:,1:), intent(inout) :: cls
    integer(i4b),                   intent(in)  :: lmin_tab, lmax_tab

    integer(i4b)        :: unit, i, l
    real(dp)            :: TT, TE, EE, BB

    unit         = getlun()

    open(unit,file=trim(filename))
    if (lmin_tab > 2) then
       do l = 2, lmin_tab-1
          read(unit,*) i, TT, EE, BB, TE
       end do
    end if
    do l = lmin_tab, lmax_tab
       read(unit,*) i, TT, EE, BB, TE
       cls(l,1) = TT / (real(l*(l+1),dp) / (2.d0*pi))
       cls(l,2) = TE / (real(l*(l+1),dp) / (2.d0*pi))
       cls(l,4) = EE / (real(l*(l+1),dp) / (2.d0*pi))
       cls(l,6) = BB / (real(l*(l+1),dp) / (2.d0*pi))
    end do
    close(unit)
    
  end subroutine input_powspec_simple


  function get_npar(bins, enable_spec, num_data_set, num_lowl_junk_bins, num_highl_junk_bins)
    implicit none

    integer(i4b), dimension(1:,1:), intent(in)  :: bins
    logical(lgt), dimension(6),     intent(in)  :: enable_spec
    integer(i4b),                   intent(in)  :: num_data_set
    integer(i4b),                   intent(in)  :: num_lowl_junk_bins, num_highl_junk_bins
    integer(i4b)                                :: get_npar
    
    integer(i4b) :: nspec, num_cmb_bins, num_junk_bins

    nspec         = count(enable_spec)
    num_junk_bins = num_lowl_junk_bins + num_highl_junk_bins
    num_cmb_bins  = size(bins(:,1)) - num_junk_bins

    get_npar = nspec * (num_cmb_bins + num_data_set * num_junk_bins)

  end function get_npar


  subroutine initialize_index_mapping(bins, enable_spec, num_data_set, &
         & num_lowl_junk_bins, num_highl_junk_bins, ind2bin, ind2spec, ind2set)
    implicit none

    integer(i4b), dimension(1:,1:), intent(in)  :: bins
    logical(lgt), dimension(6),     intent(in)  :: enable_spec
    integer(i4b),                   intent(in)  :: num_data_set
    integer(i4b),                   intent(in)  :: num_lowl_junk_bins, num_highl_junk_bins
    integer(i4b), dimension(1:),    intent(out) :: ind2bin, ind2spec, ind2set

    integer(i4b) :: i, j, k, nbin, nspec, m, n, c, npar
    
    nbin  = size(bins(:,1))
    nspec = count(enable_spec)
    npar  = get_npar(bins, enable_spec, num_data_set, num_lowl_junk_bins, num_highl_junk_bins)
    m     = npar / nspec
    n     = m - num_data_set * (num_lowl_junk_bins + num_highl_junk_bins)

    ! Do all CMB bins
    c = 0
    do k = 1, 6
       if (enable_spec(k)) then
          do i = 1, m
             c = c+1
             ! Set spectrum, bin and data set
             ind2spec(c) = k
             if (i <= n) then
                ind2set(c) = 0
                ind2bin(c) = i + num_lowl_junk_bins
             else if (i <= n+num_data_set*num_lowl_junk_bins) then
                ind2set(c) = (i-n-1) / num_lowl_junk_bins + 1
                ind2bin(c) = i - n - num_lowl_junk_bins * (ind2set(c)-1)
             else
                ind2set(c) = (i-n-num_lowl_junk_bins*num_data_set-1) / num_highl_junk_bins + 1
                ind2bin(c) = i-n-num_lowl_junk_bins*num_data_set + n + num_lowl_junk_bins - &
                     & num_highl_junk_bins * (ind2set(c)-1)
             end if
          end do
       end if
    end do

  end subroutine initialize_index_mapping


  ! Cls is (0:lmax, nspec=6, nset). For CMB,
  ! cls is the same for all nsets.
  subroutine convert_x2cls(x, cls)
    implicit none

    type(powspec),                     intent(in)  :: x
    real(dp),     dimension(0:,1:,1:), intent(out) :: cls

    integer(i4b) :: i, j, k, l, nbin, counter
    real(dp)     :: x0

    cls     = 0.d0
    counter = 1
    do j = 1, 6
       if (x%enable_spec(j)) then
          
          ! Do all CMB bins
          do i = x%num_lowl_junk_bins+1, x%nbin-x%num_highl_junk_bins
             x0 = x%coeff(counter)
             do l = x%bins(i,1), x%bins(i,2)
                cls(l,j,:) = x0 / (l*(l+1)/(2.d0*pi))
             end do
             counter = counter+1
          end do

          ! Do all low-l junk bins
          do k = 1, x%num_data_set
             do i = 1, x%num_lowl_junk_bins
                x0 = x%coeff(counter)
                do l = x%bins(i,1), x%bins(i,2)
                   cls(l,j,k) = x0 / (l*(l+1)/(2.d0*pi))
                end do
                counter = counter+1
             end do
          end do

          ! Do all high-l junk bins
          do k = 1, x%num_data_set
             do i = x%nbin-x%num_highl_junk_bins+1, x%nbin
                x0 = x%coeff(counter)
                do l = x%bins(i,1), x%bins(i,2)
                   cls(l,j,k) = x0 / (l*(l+1)/(2.d0*pi))
                end do
                counter = counter+1
             end do
          end do
       end if
    end do

  end subroutine convert_x2cls


  subroutine convert_cls2x(cls, x)
    implicit none

    real(dp),     dimension(0:,1:,1:), intent(in)     :: cls
    type(powspec),                     intent(inout)  :: x

    integer(i4b) :: i, j, k, l, nbin, counter, x0

    counter = 1
    do j = 1, 6
       if (x%enable_spec(j)) then
          
          ! Do all CMB bins
          do i = x%num_lowl_junk_bins+1, x%nbin-x%num_highl_junk_bins
             l = x%bins(i,1)
             x%coeff(counter) = cls(l,j,1) * (l*(l+1)/(2.d0*pi))
             counter = counter+1
          end do

          ! Do all low-l junk bins
          do k = 1, x%num_data_set
             do i = 1, x%num_lowl_junk_bins
                l  = x%bins(i,1)
                x%coeff(counter) = cls(l,j,k) * (l*(l+1)/(2.d0*pi))
                counter = counter+1
             end do
          end do

          ! Do all high-l junk bins
          do k = 1, x%num_data_set
             do i = x%nbin-x%num_highl_junk_bins+1, x%nbin
                l  = x%bins(i,1)
                x%coeff(counter) = cls(l,j,k) * (l*(l+1)/(2.d0*pi))
                counter = counter+1
             end do
          end do
       end if
    end do

  end subroutine convert_cls2x

  subroutine compute_rotation_angle_old(vec1, vec2, t2a)
    implicit none

    real(dp), dimension(3),   intent(in)  :: vec1, vec2
    real(dp), dimension(2,2), intent(out) :: t2a

    real(dp) :: len_u, len_v, z, cos_theta, sgn
    real(dp), dimension(3) :: u, v
    
    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       t2a(1,1) = 1.d0
       t2a(2,1) = 0.d0
       t2a(1,2) = 1.d0
       t2a(2,2) = 0.d0
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
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    t2a(1,1)  = 2.d0*cos_theta**2 -1.d0
    t2a(2,1)  = -sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta

    ! Rotation from vec2 to vec 1
    u         = vec2(3) * vec2
    u(3)      = u(3) - 1.d0
    v         = vec1 - z * vec2
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    cos_theta = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    t2a(1,2)  = 2.d0*cos_theta**2 -1.d0
    t2a(2,2)  = sgn * 2.d0 * sqrt(1.d0-cos_theta**2) * cos_theta

  end subroutine compute_rotation_angle_old



  subroutine output_2d_prob(filename, q, n, lnL)
    implicit none

    character(len=*),                 intent(in) :: filename
    real(dp),         dimension(:),   intent(in) :: q, n
    real(dp),         dimension(:,:), intent(in) :: lnL

    integer(i4b) :: i, j, unit

    unit = getlun()

    open(unit,file=trim(filename), recl=1024)
    write(unit,*) size(q), size(n)
    do i = 1, size(q)
       do j = 1, size(n)
          write(unit,*) q(i), n(j), lnL(i,j)
       end do
    end do
    close(unit)

  end subroutine output_2d_prob

  subroutine output_1d_prob(filename, x, P)
    implicit none

    character(len=*),                 intent(in) :: filename
    real(dp),         dimension(:),   intent(in) :: x, P

    integer(i4b) :: i, unit

    unit = getlun()

    open(unit,file=trim(filename), recl=1024)
    do i = 1, size(x)
       write(unit,*) x(i), P(i)
    end do
    close(unit)
    
  end subroutine output_1d_prob



end module map2cl_utils
