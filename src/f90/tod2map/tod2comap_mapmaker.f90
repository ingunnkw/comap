module tod2comap_mapmaker
  !use comap_lx_mod
  use comap_map_mod
  use comap_acceptlist_mod
  use comap_patch_mod
  use comap_ephem_mod
  !use quiet_fft_mod
  use tod2comap_utils
  use comap_jackknife_mod
  implicit none



contains


  subroutine initialize_mapmaker(map, parfile, pinfo, jk_info)
    implicit none
    !type(tod_type), dimension(:), intent(in)    :: tod
    type(map_type),   intent(inout) :: map
    type(patch_info), intent(in)    :: pinfo
    type(jk_type),    intent(in)    :: jk_info
    character(len=*)                :: parfile

    integer(i4b) :: i, j, k, l, p, q, fs, st, ierr, njkfeed, n
    real(dp)     :: x_min, x_max, y_min, y_max, pad, temp, mean_dec, d1, d2
    character(len=15)   :: coord_system
    real(8), parameter :: PI = 4*atan(1.d0)

    ! Set up map grid
    if (allocated(map%x)) return
    call get_parameter(0, parfile, 'NUMFREQ', par_int=map%nfreq)
    call get_parameter(0, parfile, 'NUM_SIDEBAND', par_int=map%nsb)
    call get_parameter(0, parfile, 'NUM_DET', par_int=map%ndet_tot)
    call get_parameter(0, parfile, 'COORDINATE_SYSTEM', par_string=coord_system)
    call get_parameter(0, parfile, 'N_NOISE_SIMULATIONS', par_int=map%nsim)
    
    !fs = 1!200
    !st = tod%nsamp!-200 ! tod%nsamp
    pad = 0.5d0 ! degrees
    !write(*,*) pinfo%resolution
    map%dthetay = pinfo%resolution ! degrees (arcmin/60), resolution
    map%center = pinfo%pos
    map%nside = 4096/(pinfo%resolution*60)
    map%mean_el = 0.d0!mean(tod(:)%mean_el)
    map%mean_az = 0.d0!mean(tod(:)%mean_az)
    mean_dec = 0.d0
    if (pinfo%fixed) mean_dec = pinfo%pos(2)
    map%dthetax = map%dthetay/abs(cos(mean_dec*PI/180.d0))

    x_min = 500.d0; x_max = -500.d0
    y_min = 500.d0; y_max = -500.d0


    if (coord_system == 'horizontal') then
       x_min = 0.d0; x_max = 360.d0
       y_min = 25.d0; y_max = 80.d0
    else
       if (pinfo%fixed) then
          x_min = pinfo%pos(1) - pinfo%obj_rad /abs(cos(mean_dec*PI/180.d0))
          x_max = pinfo%pos(1) + pinfo%obj_rad /abs(cos(mean_dec*PI/180.d0))
          y_min = pinfo%pos(2) - pinfo%obj_rad 
          y_max = pinfo%pos(2) + pinfo%obj_rad 
       else
          x_min = 0.d0 - pinfo%obj_rad 
          x_max = 0.d0 + pinfo%obj_rad 
          y_min = 0.d0 - pinfo%obj_rad 
          y_max = 0.d0 + pinfo%obj_rad 
       end if
    end if

    !write(*,*) 'max min done'

!!$    write(*,*) x_min, x_max
!!$    write(*,*) y_min, y_max
!!$    call mpi_finalize(i)
!!$    stop
    
    !write(*,*) map%dthetax, map%dthetay
    !map%n_x = int((x_max-x_min)/map%dthetax); map%n_y = int((y_max-y_min)/map%dthetay)
    map%n_x = int(2.d0 * pinfo%obj_rad / pinfo%resolution); map%n_y = int(2.d0 * pinfo%obj_rad / pinfo%resolution)
    !write(*,*) map%n_x, map%n_y

    allocate(map%x(map%n_x), map%y(map%n_y))
    do i = 1, map%n_x
       map%x(i) = x_min + (i-0.5d0)*map%dthetax
    end do
    do i = 1, map%n_y
       map%y(i) = y_min + (i-0.5d0)*map%dthetay
    end do

    !write(*,*) 'grid made'
 
    ! Set up map structures
    allocate(map%m(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot), &
         & map%dsum(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot), &
         & map%nhit(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot), &
         & map%div(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot), &
         & map%rms(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot), &
         & map%dsum_co(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%nhit_co(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%div_co(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%rms_co(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%m_co(map%n_x, map%n_y, map%nfreq, map%nsb), &
         & map%freq(map%nfreq, map%nsb))

    map%dsum    = 0.0
    map%nhit    = 0
    map%div     = 0.0
    map%m       = 0.0
    map%rms     = 0.0
    map%freq    = 0.d0    
    map%dsum_co = 0.0
    map%nhit_co = 0
    map%div_co  = 0.0
    map%m_co    = 0.0
    map%rms_co  = 0.0

    ! Simulated data
    if (map%nsim > 0) then
       allocate(map%m_sim(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, map%nsim), &
            & map%dsum_sim(map%n_x, map%n_y, map%nfreq, map%nsb,map%ndet_tot, map%nsim), &
            & map%div_sim(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, map%nsim), &
            & map%rms_sim(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, map%nsim))
       map%dsum_sim = 0.0
       map%div_sim  = 0.0 
       map%m_sim    = 0.0
       map%rms_sim  = 0.0
    end if

    ! Jackknives
    map%njk = jk_info%njk
    if (map%njk > 0) then
       njkfeed = sum(jk_info%feedmap); allocate(map%jk_feed(njkfeed))
       n = 1
       do i = 1, map%njk
          if (jk_info%feedmap(i) .eq. 1) then
             map%jk_feed(n) = i
             n = n + 1
          end if
       end do
       map%jk_def = jk_info%jk_name
       allocate(map%m_jk(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, 2*jk_info%n_feed), &
            & map%rms_jk(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, 2*jk_info%n_feed), &
            & map%dsum_jk(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, 2*jk_info%n_feed), &
            & map%div_jk(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, 2*jk_info%n_feed), &
            & map%nhit_jk(map%n_x, map%n_y, map%nfreq, map%nsb, map%ndet_tot, 2*jk_info%n_feed), &
            & map%m_jkco(map%n_x, map%n_y, map%nfreq, map%nsb, 2*jk_info%n_coadd), &
            & map%rms_jkco(map%n_x, map%n_y, map%nfreq, map%nsb, 2*jk_info%n_coadd), &
            & map%dsum_jkco(map%n_x, map%n_y, map%nfreq, map%nsb, 2*jk_info%n_coadd), &
            & map%div_jkco(map%n_x, map%n_y, map%nfreq, map%nsb, 2*jk_info%n_coadd), &
            & map%nhit_jkco(map%n_x, map%n_y, map%nfreq, map%nsb, 2*jk_info%n_coadd))
       map%m_jk    = 0.0; map%m_jkco    = 0.0
       map%rms_jk  = 0.0; map%rms_jkco  = 0.0
       map%nhit_jk = 0;   map%nhit_jkco = 0
       map%dsum_jk = 0.0; map%dsum_jkco = 0.0
       map%div_jk  = 0.0; map%div_jkco  = 0.0
    end if

    ! Frequency
    d1 = 1.d0/64.d0; d2 = 2.d0/64.d0
    do i = 1, map%nsb
       temp = 26.d0 + (i-1)*2.d0 + d1
       do j = 1, map%nfreq
          map%freq(j,i) = temp + (j-1)*d2
       end do
    end do

  end subroutine initialize_mapmaker

  
  subroutine time2pix(tod, map, parfile, pinfo, jk_list)
    ! Calculates which timestep correspond to which pixelnumber, and hitmaps
    ! NOT USED ANYMORE FOR BINNING!!
    implicit none
    type(tod_type),   intent(inout) :: tod
    type(map_type),   intent(inout) :: map
    character(len=*)                :: parfile
    type(patch_info), intent(in)    :: pinfo
    integer(i4b), dimension(:,:), intent(in) :: jk_list

    integer(i4b) :: scan, i, j, p, q, sb, freq, det
    real(dp)     :: x_min, x_max, y_min, y_max, t1, t2
    character(len=15) :: coord_system, object
    real(dp), allocatable, dimension(:,:) :: sigma0
    real(dp), dimension(3) :: pos

    call get_parameter(0, parfile, 'COORDINATE_SYSTEM', par_string=coord_system)
    call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)
    call initialize_comap_ephem_mod(parfile)

    x_min = map%x(1); x_max = map%x(map%n_x)
    y_min = map%y(1); y_max = map%y(map%n_y)
    allocate(sigma0(tod%nfreq,tod%nsb))
    sigma0 = 0.d0

    map%ndet = tod%ndet
    if (.not. allocated(map%feeds)) allocate(map%feeds(map%ndet))
    map%feeds = tod%feeds
     
    !call wall_time(t1)
    !!$OMP PARALLEL PRIVATE(scan,j,i,p,q,sb,freq)
    !!$OMP DO SCHEDULE(guided)

    ! Get position of a non-fixed object (assuming ~no movement during one subscan)
    pos = 0.d0
    !write(*,*) pinfo%fixed
    if (.not. pinfo%fixed) pos = get_obj_info(object,tod%t(1))
    !write(*,*) pos

    do j = 1, tod%ndet
       det = tod%feeds(j)
       if (.not. is_alive(det)) cycle
       !write(*,*) det
       !sigma0(:,:) = sigma0(:,:) + 1.0 / tod%sigma0(:,:,det)**2
       !write(*,*) j, det
       do i = 1, tod%nsamp
          if (tod%test(i) .eq. 1) cycle
          !if (tod%point_tel(2,i,j) .lt. 35.d0) cycle
          !if (tod%point_tel(2,i,j) .gt. 65.d0) cycle
          !if (tod%point_tel(1,i,j) .lt. 40.d0) cycle
          if (trim(coord_system) .eq. 'horizontal') then
             p = nint((tod%point_tel(1,i,j)-x_min)/map%dthetax) + 1
             q = nint((tod%point_tel(2,i,j)-y_min)/map%dthetay) + 1
          else
             !if (i .eq. 1) write(*,*) det, tod%point(1,i,det)
             p = nint((tod%point(1,i,j)-pos(1)-x_min)/map%dthetax) + 1
             q = nint((tod%point(2,i,j)-pos(2)-y_min)/map%dthetay) + 1
             !p = min(max(nint((tod%point(1,i,det)-x_min)/map%dthetax),1),map%n_x)
             !q = min(max(nint((tod%point(2,i,det)-y_min)/map%dthetay),1),map%n_y)
          end if
          if ((p .ge. 1) .and. (p .le. map%n_x) .and. (q .ge. 1) .and. (q .le. map%n_y)) then!((1 <= p <= map%n_x) .and. (1 <= q <= map%n_y)) then   
             tod%pixel(i,det) = (q-1)*map%n_x + p
             do sb = 1, tod%nsb
                if (jk_list(sb,det) == 0) cycle
                do freq = 1, tod%nfreq
                   if (tod%freqmask(freq,sb,j) == 0) cycle
                   !!$OMP ATOMIC
                   ! Flipping sb 1 and 3
                   if ((sb .eq. 1) .or. (sb .eq. 3)) then
                      map%nhit(p,q,tod%nfreq-freq+1,sb,det) = map%nhit(p,q,tod%nfreq-freq+1,sb,det) + 1
                      map%nhit_co(p,q,tod%nfreq-freq+1,sb)  = map%nhit_co(p,q,tod%nfreq-freq+1,sb)  + 1                      
                   else
                      map%nhit(p,q,freq,sb,det) = map%nhit(p,q,freq,sb,det) + 1
                      map%nhit_co(p,q,freq,sb)  = map%nhit_co(p,q,freq,sb)  + 1
                   end if
                 end do
             end do
          else
             tod%pixel(i,det) = -200
          end if
       end do
    end do
    !!$OMP END DO
    !!$OMP END PARALLEL
    !call wall_time(t2)
    !write(*,*) 'Wall time time2pix = ', t2-t1
    
    !do sb = 1, tod%nsb
    !   do freq = 1, tod%nfreq
    !      if (sigma0(freq,sb) == 0) cycle
    !      where(map%nhit(:,:,freq,sb,:) > 0)
    !         map%rms(:,:,freq,sb,:) = map%rms(:,:,freq,sb,:) + map%nhit(:,:,freq,sb,:)/sigma0(freq,sb)**2 
    !      elsewhere
    !         map%rms(:,:,freq,sb,:) = 0.d0
    !      end where
    !   end do
    !end do

    !call mpi_finalize(i)
    !stop
    deallocate(sigma0)

  end subroutine time2pix


!  subroutine time2pix_sim(tod, map, parfile, pinfo)
!    implicit none
!    type(tod_type),   intent(inout) :: tod
!    type(map_type),   intent(inout) :: map
!    character(len=*)                :: parfile
!    type(patch_info), intent(in)    :: pinfo

!    integer(i4b) :: scan, i, j, p, q, sb, freq, det
!    real(dp)     :: x_min, x_max, y_min, y_max, t1, t2
!    character(len=15) :: coord_system, object
!    real(sp), allocatable, dimension(:,:) :: sigma0
!    real(dp), dimension(3) :: pos

!    call get_parameter(0, parfile, 'COORDINATE_SYSTEM', par_string=coord_system)
!    call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)
!    call initialize_comap_ephem_mod(parfile)

!    x_min = map%x(1); x_max = map%x(map%n_x)
!    y_min = map%y(1); y_max = map%y(map%n_y)

!    map%ndet = tod%ndet
!    if (.not. allocated(map%feeds)) allocate(map%feeds(map%ndet))
!    map%feeds = tod%feeds
    
!    ! Get position of a non-fixed object (assuming ~no movement during one subscan)
!    pos = 0.d0
!    !write(*,*) pinfo%fixed
!    if (.not. pinfo%fixed) pos = get_obj_info(object,tod%t(1))
!    !write(*,*) pos

!    do j = 1, tod%ndet
!       det = tod%feeds(j)
!       if (.not. is_alive(det)) cycle
!       !write(*,*) j, det
!       do i = 1, tod%nsamp
!          if (tod%test(i) .eq. 1) cycle
!          !if (tod%point_tel(2,i,j) .lt. 35.d0) cycle
!          !if (tod%point_tel(2,i,j) .gt. 65.d0) cycle
          !if (tod%point_tel(1,i,j) .lt. 40.d0) cycle
!          if (trim(coord_system) .eq. 'horizontal') then
!             p = nint((tod%point_tel(1,i,j)-x_min)/map%dthetax)
!             q = nint((tod%point_tel(2,i,j)-y_min)/map%dthetay)
!          else
!             !if (i .eq. 1) write(*,*) det, tod%point(1,i,det)
!             p = nint((tod%point(1,i,j)-pos(1)-x_min)/map%dthetax)
!             q = nint((tod%point(2,i,j)-pos(2)-y_min)/map%dthetay)
!             !p = min(max(nint((tod%point(1,i,det)-x_min)/map%dthetax),1),map%n_x)
!             !q = min(max(nint((tod%point(2,i,det)-y_min)/map%dthetay),1),map%n_y)
!          end if
!          if ((p .ge. 1) .and. (p .le. map%n_x) .and. (q .ge. 1) .and. (q .le. map%n_y)) then!((1 <= p <= map%n_x) .and. (1 <= q <= map%n_y)) then   
!             tod%pixel(i,det) = (q-1)*map%n_x + p
!             do sb = 1, tod%nsb
!                do freq = 1, tod%nfreq
!                   if (tod%freqmask(freq,sb,j) == 0) cycle
!                   ! Flipping sb 1 and 3
!                   if ((sb .eq. 1) .or. (sb .eq. 3)) then
!                      map%nhit_co(p,q,tod%nfreq-freq+1,sb)  = map%nhit_co(p,q,tod%nfreq-freq+1,sb)  + 1
!                   else
!                      map%nhit_co(p,q,freq,sb)  = map%nhit_co(p,q,freq,sb)  + 1
!                   end if
!                end do
!             end do
!          else
!             tod%pixel(i,det) = -200
!          end if
!       end do
!    end do


!  end subroutine time2pix_sim



  subroutine rotation_matrix(ra_center, dec_center, ra, dec)!x, y, z)
    implicit none
    real(dp), intent(in)    :: ra_center, dec_center
    real(dp), intent(inout) :: ra, dec!x, y, z
    !real(dp) :: ra, dec
    real(dp) :: cosd, sind, cosa, sina
    real(dp) :: x_temp, y_temp, z_temp, x, y, z
    real(8), parameter :: PI = 4*atan(1.d0)

    ra = ra_center * PI/180.d0
    dec = (90.d0 - dec_center) * PI/180.d0

    cosd = cos(dec); sind = sin(dec)
    cosa = cos(ra); sina = sin(ra)

    x_temp = cosd*cosa*x + cosd*sina*y - sind*z
    y_temp = -sina*x + cosa*y
    z_temp = sind*cosa*x + sind*sina*y + cosd*z

    x = x_temp; y = y_temp; z = z_temp

  end subroutine rotation_matrix



  subroutine gnomonic(ra, dec, ra_center, dec_center)
    implicit none
    real(dp) :: ra, dec, ra_center, dec_center
    real(dp) :: ra_c, dec_c, x, y, cosc, cosd
    real(8), parameter :: PI = 4*atan(1.d0)

    dec = (90.d0 - dec) * PI/180.d0; ra  = ra * PI/180.d0
    dec_c = (90.d0 - dec_center) * PI/180.d0
    ra_c = ra_center * PI/180.d0

    cosd = cos(dec)

    cosc = sin(dec_c)*sin(dec) + cos(dec_c)*cosd*cos(ra-ra_c)

    x = (cosd * sin(ra-ra_c)) / cosc
    y = (cos(dec_c)*sin(dec) - sin(dec_c)*cosd*cos(ra-ra_c)) / cosc

   
  end subroutine gnomonic



  
  !!!!!!!!!!!!!!!!!
  ! Naive Binning !
  !!!!!!!!!!!!!!!!!

  subroutine binning(map, tod, scan, parfile, pinfo, jk_list, jk_split)!(map_tot, map_scan, tod, alist, scan, parfile)
    implicit none
    type(tod_type),   intent(in)    :: tod
    type(patch_info), intent(in)    :: pinfo
    type(map_type),   intent(inout) :: map!_tot, map_scan
    !type(acceptlist), intent(in)    :: alist
    character(len=*)                :: parfile
    integer(i4b), dimension(:,:),   intent(in) :: jk_list
    integer(i4b), dimension(:,:,:), intent(in) :: jk_split

    integer(i4b) :: det, sb, freq, sim, ndet, nsb, nfreq, nf, nc, freq_new
    integer(i4b) :: i, j, k, l, p, q, fs, st, scan, pix, jk, split
    real(dp)     :: x_min, x_max, y_min, y_max
    character(len=15) :: coord_system, object
    real(dp), allocatable, dimension(:) :: dsum, div
    real(dp), dimension(3) :: pos

    call get_parameter(0, parfile, 'COORDINATE_SYSTEM', par_string=coord_system)
    call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)
    call initialize_comap_ephem_mod(parfile)

    x_min = map%x(1); x_max = map%x(map%n_x)
    y_min = map%y(1); y_max = map%y(map%n_y)

    map%ndet = tod%ndet ! Has this been established earlier???
    ndet = map%ndet
    nsb = map%nsb 
    nfreq = map%nfreq

    pos = 0.d0
    if (.not. pinfo%fixed) pos = get_obj_info(object, tod%t(1))
    !! Where does ephem_mod fit in ??????????

    !write(*,*) 'Beginning coadding'

    !fs = 200 ! starting point
    !st = tod(scan)%nsamp - 200 ! ending point

    do i = 1, tod%nsamp
       do j = 1, ndet
          det = tod%feeds(j)
          if (.not. is_alive(det)) cycle
          !if (tod%pixel(i,det) .lt. 0) cycle
          !if (tod%point_tel(2,i,det) > 30.d0) cycle
          if (trim(coord_system) .eq. 'horizontal') then
             p = nint((tod%point_tel(1,i,j)-x_min)/map%dthetax) + 1
             q = nint((tod%point_tel(2,i,j)-y_min)/map%dthetay) + 1
          else
             p = nint((tod%point(1,i,j)-pos(1)-x_min)/map%dthetax) + 1
             q = nint((tod%point(2,i,j)-pos(2)-y_min)/map%dthetay) + 1
          end if
          if ((p .le. 1) .or. (p .ge. map%n_x) .or. (q .le. 1) .or. (q .ge. map%n_y)) cycle ! discard points outside of the set grid
          do sb = 1, nsb
             if (jk_list(sb,det) == 0) cycle
             do freq = 1, nfreq
                !if (tod%fknee(freq,sb,det) > 0.5d0) cycle
                if (tod%freqmask(freq,sb,j) == 0) cycle
                if (tod%rms(freq,sb,j) == 0.d0) cycle
                !write(*,*) i, det, sb, freq
                !if (any(alist%ascans(scan)%adet_sb(det,sb)%rejected == freq)) cycle
                !write(*,*) tod%rms(i,freq,sb,det)
                

                ! Flipping sb 1 and 3
                freq_new = freq
                if ((sb .eq. 1) .or. (sb .eq. 3)) freq_new = nfreq-freq+1
                map%nhit(p,q,freq_new,sb,det) = map%nhit(p,q,freq_new,sb,det) + 1.d0
                map%nhit_co(p,q,freq_new,sb)  = map%nhit_co(p,q,freq_new,sb)  + 1.d0
                map%dsum(p,q,freq_new,sb,det) = map%dsum(p,q,freq_new,sb,det) + 1.d0 / tod%rms(freq,sb,j)**2 * tod%d(i,freq,sb,j)
                map%div(p,q,freq_new,sb,det)  = map%div(p,q,freq_new,sb,det)  + 1.d0 / tod%rms(freq,sb,j)**2
                map%dsum_co(p,q,freq_new,sb)  = map%dsum_co(p,q,freq_new,sb)  + 1.d0 / tod%rms(freq,sb,j)**2 * tod%d(i,freq,sb,j)
                map%div_co(p,q,freq_new,sb)   = map%div_co(p,q,freq_new,sb)   + 1.d0 / tod%rms(freq,sb,j)**2
                
                ! Jackknives
                nf = 1; nc = 1
                do jk = 1, map%njk
                   if (any(map%jk_feed == jk)) then
                      !if (jk_split(i,sb,det) == 0) then
                      split = 2*nf - 1 + jk_split(jk,sb,det)
                      map%nhit_jk(p,q,freq_new,sb,det,split) = map%nhit_jk(p,q,freq_new,sb,det,split) + 1
                      map%dsum_jk(p,q,freq_new,sb,det,split) = map%dsum_jk(p,q,freq_new,sb,det,split) + 1.0 / tod%rms(freq,sb,j)**2 * tod%d(i,freq,sb,j)
                      map%div_jk(p,q,freq_new,sb,det,split)  = map%div_jk(p,q,freq_new,sb,det,split)  + 1.0 / tod%rms(freq,sb,j)**2
                      nf = nf+1
                   else
                      split = 2*nc - 1 + jk_split(jk,sb,det)
                      map%nhit_jkco(p,q,freq_new,sb,split) = map%nhit_jkco(p,q,freq_new,sb,split) + 1
                      map%dsum_jkco(p,q,freq_new,sb,split) = map%dsum_jkco(p,q,freq_new,sb,split) + 1.0 / tod%rms(freq,sb,j)**2 * tod%d(i,freq,sb,j)
                      map%div_jkco(p,q,freq_new,sb,split)  = map%div_jkco(p,q,freq_new,sb,split)  + 1.0 / tod%rms(freq,sb,j)**2
                      nc = nc+1
                   end if
                end do

                   ! Simulations in here

             end do
          end do
       end do
    end do
 
    !map_tot%dsum    = map_tot%dsum + map_scan%dsum
    !map_tot%div     = map_tot%div  + map_scan%div
    !map_tot%dsum_co = map_tot%dsum_co + map_scan%dsum_co
    !map_tot%div_co  = map_tot%div_co  + map_scan%div_co

    !write(*,*) 'Ending coadding'

  end subroutine binning


!  subroutine binning_sim(map_tot, map_scan, tod, scan, parfile, pinfo)!(map_tot, map_scan, tod, alist, scan, parfile)
!    implicit none
!    type(tod_type),   intent(in)    :: tod
!    type(patch_info), intent(in)    :: pinfo
!    type(map_type),   intent(inout) :: map_tot, map_scan
!    !type(acceptlist), intent(in)    :: alist
!    character(len=*)                :: parfile

!    integer(i4b) :: det, sb, freq, sim, ndet, nsb, nfreq
!    integer(i4b) :: i, j, k, l, p, q, fs, st, scan, pix
!    real(dp)     :: x_min, x_max, y_min, y_max
!    character(len=15) :: coord_system, object
!    real(dp), allocatable, dimension(:) :: dsum, div
!    real(dp), dimension(3) :: pos

!    call get_parameter(0, parfile, 'COORDINATE_SYSTEM', par_string=coord_system)
!    call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)

!    x_min = map_tot%x(1); x_max = map_tot%x(map_tot%n_x)
!    y_min = map_tot%y(1); y_max = map_tot%y(map_tot%n_y)

!    ndet = map_tot%ndet
!    nsb = map_tot%nsb 
!    nfreq = map_tot%nfreq
    
!    pos = 0.d0
!    if (.not. pinfo%fixed) pos = get_obj_info(object, tod%t(1))

!    do i = 1, tod%nsamp
!       do j = 1, ndet
!          det = tod%feeds(j)
!          if (.not. is_alive(det)) cycle
!          if (tod%pixel(i,det) .lt. 0) cycle
!          !if (tod%point_tel(2,i,det) > 30.d0) cycle
!          if (trim(coord_system) .eq. 'horizontal') then
!             p = nint((tod%point_tel(1,i,j)-x_min)/map_tot%dthetax)
!             q = nint((tod%point_tel(2,i,j)-y_min)/map_tot%dthetay)
!          else
!             p = nint((tod%point(1,i,j)-pos(1)-x_min)/map_tot%dthetax)
!             q = nint((tod%point(2,i,j)-pos(2)-y_min)/map_tot%dthetay)
!             !p = min(max(nint((tod%point(1,i,det)-x_min)/map_tot%dthetax),1),map_tot%n_x)
!             !q = min(max(nint((tod%point(2,i,det)-y_min)/map_tot%dthetay),1),map_tot%n_y)
!          end if
!          do sb = 1, nsb
!             do freq = 1, nfreq
!                !if (tod%fknee(freq,sb,det) > 0.5d0) cycle
!                if (tod%freqmask(freq,sb,j) == 0) cycle
!                if (tod%rms_sim(freq,sb,j,1) == 0.d0) cycle
!                !map_scan%dsum_sim(p,q,freq,sb,det) = map_scan%dsum_sim(p,q,freq,sb,det) + 1.d0 / tod%rms_sim(freq,sb,j)**2 * tod%d_sim(i,freq,sb,j) 
!                !map_scan%div_sim(p,q,freq,sb,det) = map_scan%div_sim(p,q,freq,sb,det) + 1.d0 / tod%rms_sim(freq,sb,j)**2 
!                ! Flipping sb 1 and 3
!                if ((sb .eq. 1) .or. (sb .eq. 3)) then
!                   map_scan%dsum_sim(p,q,nfreq-freq+1,sb,1)  = map_scan%dsum_sim(p,q,nfreq-freq+1,sb,1) + 1.d0 / tod%rms_sim(freq,sb,j)**2 * tod%d_sim(i,freq,sb,j)
!                   map_scan%div_sim(p,q,nfreq-freq+1,sb,1)   = map_scan%div_sim(p,q,nfreq-freq+1,sb,1) + 1.d0 / tod%rms_sim(freq,sb,j)**2
!                else
!                   map_scan%dsum_sim(p,q,freq,sb,1)  = map_scan%dsum_sim(p,q,freq,sb,1) + 1.d0 / tod%rms_sim(freq,sb,j)**2 * tod%d_sim(i,freq,sb,j)
!                   map_scan%div_sim(p,q,freq,sb,1)   = map_scan%div_sim(p,q,freq,sb,1) + 1.d0 / tod%rms_sim(freq,sb,j)**2
!                end if
!             end do
!          end do
!       end do
!    end do

!    ! Simulated data
!    map_tot%dsum_sim = map_tot%dsum_sim + map_scan%dsum_sim
!    map_tot%div_sim = map_tot%div_sim + map_scan%div_sim

!    !write(*,*) 'Ending coadding'

!  end subroutine binning_sim



  subroutine finalize_binning(map)
    implicit none
    type(map_type), intent(inout) :: map
    
    integer(i4b) :: det, sb, freq, p, q

    where(map%div > 0)
       map%m   = map%dsum / map%div
       map%rms = 1.0 / sqrt(map%div)
    elsewhere
       map%m   = 0.0
       map%rms = 0.0
    end where

    ! Co-added feeds
    where(map%div_co > 0)
       map%m_co   = map%dsum_co / map%div_co
       map%rms_co = 1.0 / sqrt(map%div_co)
    elsewhere
       map%m_co   = 0.0
       map%rms_co = 0.0
    end where 

    ! Simulated data
    where(map%div_sim > 0)
       map%m_sim   = map%dsum_sim / map%div_sim 
       map%rms_sim = 1.0 / sqrt(map%div_sim) 
    elsewhere
       map%m_sim   = 0.0 
       map%rms_sim = 0.0 
    end where

    ! Jackknives

    where(map%div_jk > 0)
       map%m_jk   = map%dsum_jk / map%div_jk
       map%rms_jk = 1.0 / sqrt(map%div_jk)
    elsewhere
       map%m_jk   = 0.0
       map%rms_jk = 0.0
    end where 

    where(map%div_jkco > 0)
       map%m_jkco   = map%dsum_jkco / map%div_jkco
       map%rms_jkco = 1.0 / sqrt(map%div_jkco)
    elsewhere
       map%m_jkco   = 0.0
       map%rms_jkco = 0.0
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
           !rhs(tod(scan)%pixel(i,det)) = rhs(tod(scan)%pixel(i,det)) + tod(scan)%g(1,freq,sb,det)*Nt(i)
           rhs(tod(scan)%pixel(i,det)) = rhs(tod(scan)%pixel(i,det)) + Nt(i)
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
           !xt(i) = tod(scan)%g(1,freq,sb,det)*x(tod(scan)%pixel(i,det))
           !xt(long_samp-i+1) = tod(scan)%g(1,freq,sb,det)*x(tod(scan)%pixel(i,det))
           xt(i) = x(tod(scan)%pixel(i,det))
           xt(long_samp-i+1) = x(tod(scan)%pixel(i,det))
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
           !Ax(tod(scan)%pixel(i,det)) = Ax(tod(scan)%pixel(i,det)) + tod(scan)%g(1,freq,sb,det)*Nt(i)
           Ax(tod(scan)%pixel(i,det)) = Ax(tod(scan)%pixel(i,det)) + Nt(i)
        end do

        deallocate(xt, xv)
        deallocate(Nt, Nv)
     end do

   end subroutine get_lhs


  subroutine pcg_mapmaker(tod, map, det, sb, freq, parfile)!(tod, map, alist, det, sb, freq, parfile)
    implicit none
    type(tod_type),   dimension(:), intent(in) :: tod
    type(map_type),   intent(inout) :: map
    !type(acceptlist), intent(in)    :: alist
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
          map%m(j,k,freq,sb,det) = map%m(j,k,freq,sb,det) + mp((k-1)*map%n_x + j)
       end do
    end do
    
    ! allocate(debug(tod(1)%nsamp))
    ! ! Debugging
    ! call wall_time(t1)
    ! !$OMP PARALLEL PRIVATE(i)
    ! !$OMP DO SCHEDULE(guided)
    ! do i = 1, tod(1)%nsamp
    !    debug(i) = mp(tod(1)%pixel(i,det))
    ! end do
    ! !$OMP END DO
    ! !$OMP END PARALLEL
    ! call wall_time(t2)
    ! write(*,*) 'wall time debug = ', t2-t1
    ! open(75,file='map_time.dat',recl=1024)
    ! do i = 1, tod(1)%nsamp
    !    write(75,*) i, debug(i), tod(1)%pixel(i,det), tod(1)%point_tel(1,i,det), tod(1)%point_tel(2,i,det)
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
          nobs(i) = map%nhit(j,k,freq,sb,det)
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



!  subroutine l12hitmap(l1file, map, parfile, pinfo)
!    implicit none
!    type(map_type), intent(inout)   :: map
!    character(len=*)                :: parfile, l1file
!    type(patch_info), intent(in)    :: pinfo
!
!
!    type(lx_struct) :: data
!    real(dp) :: x_min, x_max, y_min, y_max, pos(3)
!    integer(i4b):: i, j, k, l, p, q, sb, freq, ndet, nsb, nfreq, nsamp
!    character(len=512) :: object
!
!
!    call get_parameter(0, parfile, 'TARGET_NAME', par_string=object)
!    call initialize_comap_ephem_mod(parfile)
!
!    ! Read the l1 file
!    write(*,*) 'reading file'
!    call read_l1_file(l1file, data, 0)
!
!    nsamp = size(data%time)
!    nfreq = size(data%nu,1)
!    nsb   = size(data%tod,3)
!    ndet  = size(data%tod,4)
!
!
!    ! Get pointing info
!    write(*,*) 'getting pointing info'
!
!    if (.not. pinfo%fixed) pos = get_obj_info(object,data%time(1))
!
!    x_min = map%x(1); x_max = map%x(map%n_x)
!    y_min = map%y(1); y_max = map%y(map%n_y)
!
!
!    do j = 1, ndet
!       if (.not. is_alive(j)) cycle
!       do i = 1, nsamp
!          p = nint((data%point_cel(1,i,j)-pos(1)-x_min)/map%dthetax)
!          q = nint((data%point_cel(2,i,j)-pos(2)-y_min)/map%dthetay)
!          if ((p .ge. 1) .and. (p .le. map%n_x) .and. (q .ge. 1) .and. (q .le. map%n_y)) then
!             !tod%pixel(i,j) = (q-1)*map%n_x + p
!             do sb = 1, nsb
!                do freq = 1, nfreq
!                   !if (tod%freqmask(freq,sb,j) == 0) cycle
!                   map%nhit(p,q,freq,sb,j) = map%nhit(p,q,freq,sb,j) + 1.d0
!                end do
!             end do
!          !else
!          !   tod%pixel(i,j) = -200
!          end if
!       end do
!    end do
!
!    call free_lx_struct(data)
!
!  end subroutine l12hitmap




end module
