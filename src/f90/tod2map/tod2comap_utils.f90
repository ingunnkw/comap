module tod2comap_utils
  use comap_lx_mod
  implicit none

  type tod_type
     real(dp)     :: samprate, Tsys
     integer(i4b) :: nsamp, ndet, nfreq, nsb
     real(dp)     :: fmin, fmax, df, mean_el
     real(dp), allocatable, dimension(:)       :: t                    ! (time) 
     real(dp), allocatable, dimension(:,:,:,:) :: d, d_long, d_raw, g, rms ! (time, freq,  sb, det)
     real(dp), allocatable, dimension(:,:,:)   :: sigma0, fknee, alpha,f ! (freq, sb, det)
     real(dp), allocatable, dimension(:,:)     :: pixels            ! (sb, freq) or (time, det)
     real(dp), allocatable, dimension(:,:,:)   :: point, point_tel     ! (det, 3, time)
  end type tod_type

contains

  subroutine get_tod(l3file, tod)
    implicit none
    character(len=*), intent(in)    :: l3file
    type(tod_type),   intent(inout) :: tod

    integer(i4b) :: i, j, k, l
    real(dp)     :: nu_cut
    type(lx_struct) :: data

    nu_cut = 0.1d0

    ! Read data
    call read_l3_file(l3file, data)
    call free_tod_type(tod)

    tod%samprate = 50.d0 !data%samprate
    tod%nsamp = size(data%time)
    tod%nfreq = size(data%nu,1)
    tod%nsb   = size(data%tod,3)
    tod%ndet  = size(data%tod,4)

!    write(*,*) tod%nsamp, tod%nfreq, tod%nsb, tod%ndet

    allocate( tod%t(tod%nsamp), tod%f(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%point(3,tod%nsamp,tod%ndet), tod%pixels(tod%nsamp, tod%ndet), &
         & tod%point_tel(3,tod%nsamp, tod%ndet), &
         & tod%d_raw(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%d(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%g(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%rms(tod%nsamp, tod%nfreq, tod%nsb, tod%ndet), &
         & tod%sigma0(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%fknee(tod%nfreq, tod%nsb, tod%ndet), &
         & tod%alpha(tod%nfreq, tod%nsb, tod%ndet))

    tod%t = data%time; tod%f = data%nu
    tod%point = data%point!_cel ! call make_angles_safe(tod%point(1,:),maxang)
    tod%point_tel = data%point_tel
    tod%g     = data%gain
    tod%sigma0 = data%sigma0
    tod%fknee  = data%fknee
    tod%alpha  = data%alpha
    tod%mean_el = mean(data%point_tel(2,:,1))
    !write(*,*) shape(data%point), tod%nsamp

    do k = 1, tod%ndet
       do l = 1, tod%nsb
          do j = 1, tod%nfreq
             !do j = 6, 6
             do i = 1, tod%nsamp
                tod%d_raw(i,j,l,k) = data%tod(i,j,l,k)
             end do

             ! Apply high pass filter
             tod%d(:,j,l,k) = tod%d_raw(:,j,l,k)
             !tod%d(:,j,l,k) = tod%d(:,j,l,k) - mean(tod%d(:,j,l,k))
             !call hp_filter(nu_cut, tod%d(:,j,l,k),tod%samprate)

             ! Estimate RMS
             tod%rms(:,j,l,k) = sqrt(variance(tod%d(:,j,l,k)))
          end do
       end do
    end do

    call free_lx_struct(data)

  end subroutine get_tod



  subroutine hp_filter(nu_cut, tod_d, samp_rate)
    implicit none   
    real(dp),                    intent(in)    :: nu_cut, samp_rate
    real(dp),     dimension(1:), intent(inout) :: tod_d
    
    integer(i4b) :: n, ind_cut
    complex(dpc), allocatable, dimension(:) :: tod_fft

    n = size(tod_d)/2 + 1
    allocate(tod_fft(n))
    ind_cut = freq2ind(nu_cut, samp_rate, n) ! n or nsamps??????
    call fft(tod_d, tod_fft, 1)
    tod_fft(1:ind_cut) = 0.d0
    call fft(tod_d, tod_fft, -1)
    
    deallocate(tod_fft)

  end subroutine hp_filter


  subroutine free_tod_type(tod)
    implicit none
    type(tod_type), intent(inout) :: tod 

    if (allocated(tod%t))      deallocate(tod%t)
    if (allocated(tod%f))      deallocate(tod%f)
    if (allocated(tod%d))      deallocate(tod%d)
    if (allocated(tod%d_raw))  deallocate(tod%d_raw)
    if (allocated(tod%g))      deallocate(tod%g)
    if (allocated(tod%rms))    deallocate(tod%rms)
    if (allocated(tod%point))  deallocate(tod%point)
    if (allocated(tod%point_tel))  deallocate(tod%point_tel)
    if (allocated(tod%sigma0)) deallocate(tod%sigma0)
    if (allocated(tod%fknee))  deallocate(tod%fknee)
    if (allocated(tod%alpha))  deallocate(tod%alpha)
    if (allocated(tod%pixels)) deallocate(tod%pixels)

  end subroutine free_tod_type

end module
