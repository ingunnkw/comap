module quiet_stat_mod
  use healpix_types
  use quiet_utils
  use math_tools
  use quiet_fft_mod
  implicit none

contains


  subroutine compute_typeB_chisq(samprate, sigma0, demod, tp, chisq)
     implicit none

     real(dp),                intent(in)  :: samprate, sigma0
     real(sp), dimension(1:), intent(in)  :: demod, tp
     real(dp),                intent(out) :: chisq

     real(dp)     :: samprate0, s0
     integer(i4b) :: deci, numsamp, i, j, k, status
     real(dp), allocatable, dimension(:)   :: d, t
     real(dp),              dimension(0:1) :: a

     if (sigma0 == 0 .or. tp(1) == tp(2)) then
        chisq = 0.d0
        return
     end if
     
     samprate0  = 1.d0 ! Samprate of downgraded data in Hz
     deci       = int(samprate/samprate0)
     numsamp    = size(demod) / deci
     s0         = sigma0 / sqrt(real(deci,dp))

     ! Generate downsampled timestreams
     allocate(d(numsamp), t(numsamp))
     do i = 1, numsamp
        d(i) = sum(demod((i-1)*deci+1:i*deci)) / real(deci,dp)
        t(i) = sum(tp((i-1)*deci+1:i*deci)) / real(deci,dp)
     end do

!     ! Fit a straight line
     call fit_polynomial(t, d, a, status)
     if(status /= 0) then
        chisq = infinity
        return
     end if

     ! Compute chi-square
     chisq = 0.d0
     do i = 1, numsamp
        chisq = chisq + (d(i) - (a(0)+a(1)*t(i)))**2 / s0**2
     end do
     chisq = (chisq-numsamp) / real(2*numsamp,dp)
     deallocate(d, t)
   end subroutine compute_typeB_chisq

   subroutine compute_weather_stat(samprate, period, tp, stat)
     implicit none

     real(dp),                intent(in)  :: samprate, period
     real(sp), dimension(1:), intent(in)  :: tp
     real(dp),                intent(out) :: stat

     integer(i4b) :: i, j, n, nrange, numsamp, deci
     real(dp)     :: mu, samprate0
     real(dp), allocatable, dimension(:) :: rms
     real(dp), allocatable, dimension(:) :: t
     samprate0  = 1.d0 ! Samprate of downgraded data in Hz
     deci       = int(samprate/samprate0)
     numsamp    = size(tp) / deci

     ! Generate downsampled timestream
     allocate(t(numsamp))
     do i = 1, numsamp
        t(i) = sum(tp((i-1)*deci+1:i*deci)) / real(deci,dp)
     end do

     n      = int(samprate0*period)
     nrange = numsamp / n
     allocate(rms(nrange))
     do i = 1, nrange
        mu     = sum(t((i-1)*n+1:i*n)) / real(n,dp)
        rms(i) = sqrt(sum((t((i-1)*n+1:i*n)-mu)**2) / real(n-1,dp))
     end do
     mu   = sum(rms) / real(nrange,dp)
     stat = sqrt(sum((rms-mu)**2) / real(nrange-1,dp))
     deallocate(rms, t)
     
   end subroutine compute_weather_stat

   subroutine jump_finder(smoothover_in, diode_in, tp, relative_jump_height)

     implicit none

     integer(i4b),            intent(in)  :: smoothover_in
     integer(i4b),            intent(in)  :: diode_in
     real(sp), dimension(1:), intent(in)  :: tp
     real(dp),                intent(out) :: relative_jump_height

     integer(i4b) :: i, numsamp, smoothover, half_smoothover, my_diode, my_position, my_position_2
     real(dp) :: dsmoothover
     real(dp) :: average_tp, dnumsamp, rms, jump_height
     real(dp), allocatable, dimension(:) :: tpsmooth

     !SKN! Disabling jump finder until:
     !SKN!  1. The hardcoded numbers are fixed
     !SKN!  2. It checks the length of the input array to make sure it doesn't overflow
     !SKN!  3. It has been tested in valgrind with very short data sets.
     relative_jump_height = 0
     return

     !set the out_counter value to zero to begin with, and see if it gets flagged
     rms=0
     jump_height=0
     relative_jump_height=-1
     smoothover = smoothover_in
     numsamp    = size(tp)
     my_diode = diode_in
     average_tp = sum(tp((1):(1000)))/1000.0
     dsmoothover = real(smoothover,dp)

     !get rms from first 1000 points
     do i = 1, 1000
        rms = rms + ((tp(i) - average_tp) * (tp(i) - average_tp))/1000.0
     end do
     !write(*,*) 'here D'
     
     rms = sqrt(rms)
     
     !if the period to smooth over is less than 5, set it equal to 5;
     !if the period to smooth over is even, add one to it such that it is an odd number;
     if (smoothover.lt.5) then 
            smoothover=5
     endif
     if (mod(smoothover,2).eq.0) then 
            smoothover=smoothover+1
     endif 
     half_smoothover = (smoothover-1)/2
     !write(*,*) 'here E'

     ! Generate smoothed total power data, smoother over a period of "smoothover"
     allocate(tpsmooth(numsamp))

     do i = 1, numsamp
             tpsmooth(i)=0
     end do

     do i = (half_smoothover + 1), (numsamp - (half_smoothover + 1))
        tpsmooth(i) = sum(tp((i-half_smoothover):(i+half_smoothover))) / dsmoothover
     end do
     !write(*,*) 'here F'
     
     do i = (half_smoothover + 1 + 1000), (numsamp - (half_smoothover + 1 + 1000))
         if (abs(tpsmooth(i+1000) - tpsmooth(i - 1000)).gt.jump_height) then
             jump_height=abs(tpsmooth(i+1000) - tpsmooth(i-1000))
             relative_jump_height=jump_height/rms
             my_position=i+1000
             my_position_2=i-1000
         end if
         if (abs(tpsmooth(i+500) - tpsmooth(i - 500)).gt.jump_height) then
             jump_height=abs(tpsmooth(i+500) - tpsmooth(i-500))
             relative_jump_height=jump_height/rms
             my_position=i+500
             my_position_2=i-500
         end if
         if (abs(tpsmooth(i+250) - tpsmooth(i - 250)).gt.jump_height) then
             jump_height=abs(tpsmooth(i+250) - tpsmooth(i-250))
             relative_jump_height=jump_height/rms
             my_position=i+250
             my_position_2=i-250
         end if
         if (abs(tpsmooth(i+20) - tpsmooth(i - 20)).gt.jump_height) then
             jump_height=abs(tpsmooth(i+20) - tpsmooth(i-20))
             relative_jump_height=jump_height/rms
             my_position=i+20
             my_position_2=i-20
         end if
     end do
     !write(*,*) 'here H'

     deallocate(tpsmooth)

   !end jump_finder     
   end subroutine jump_finder


   subroutine tp_rms_finder(tp_rms_sample_size, diode_in, tp, relative_max_rms)
     implicit none

     integer(i4b),            intent(in)  :: tp_rms_sample_size
     integer(i4b),            intent(in)  :: diode_in
     real(sp), dimension(1:), intent(in)  :: tp
     real(dp),                intent(out) :: relative_max_rms

     integer(i4b) :: i, j, numsamp, my_diode
     integer(i4b) :: number_bins
     real(dp) :: dtp_rms_sample_size
     real(dp) :: average_tp, initial_rms
     real(dp), allocatable, dimension(:) :: average_tp_section, average_rms_section, consecutive_rms_violations_value
     integer(i4b), allocatable, dimension(:) :: consecutive_rms_violations

     numsamp = size(tp)
     my_diode = diode_in
     dtp_rms_sample_size=real(tp_rms_sample_size,dp)
     number_bins=(numsamp-mod(numsamp,tp_rms_sample_size))/tp_rms_sample_size
     relative_max_rms=0


     !get rms from first 1000 points
     average_tp = sum(tp((1):(1000)))/1000.0
     do i = 1, 1000
        initial_rms = initial_rms + ((tp(i) - average_tp) * (tp(i) - average_tp))/1000.0
     end do
     
     initial_rms = sqrt(initial_rms)
     

     allocate(average_tp_section(number_bins))
     allocate(average_rms_section(number_bins))
     allocate(consecutive_rms_violations(number_bins))
     allocate(consecutive_rms_violations_value(number_bins))

     do i = 1, number_bins
         average_tp_section(i)=0
         average_rms_section(i)=0
         consecutive_rms_violations(i)=0
         consecutive_rms_violations_value(i)=0
     end do

     do i = 1, number_bins
        average_tp_section(i)= sum(tp((i*(tp_rms_sample_size-1) + 1):(i*(tp_rms_sample_size-1) + &
             & tp_rms_sample_size))) / dtp_rms_sample_size
     end do

     do i = 1, number_bins
         do j = 1, tp_rms_sample_size
            average_rms_section(i) = average_rms_section(i) + ((tp(i*(tp_rms_sample_size-1) + j) &
                 & - average_tp_section(i)) * (tp(i*(tp_rms_sample_size-1) + j) - average_tp_section(i)))
         end do
     end do


     do i = 1, number_bins
             average_rms_section(i) = sqrt(average_rms_section(i))
     end do


     do i = 1, number_bins
         if (consecutive_rms_violations(i).gt.(3*initial_rms)) then
             consecutive_rms_violations(i)=consecutive_rms_violations(i-1)+1
             consecutive_rms_violations_value(i)=max((average_rms_section(i)/initial_rms),&
                  & (average_rms_section(i-1)/initial_rms))
         end if
     end do

     !!!!!Test Case
     if (my_diode.eq.17 .or. my_diode.eq.18 .or. my_diode.eq.19 .or. my_diode.eq.20 .or. &
          & my_diode.eq.261 .or. my_diode.eq.262 .or. my_diode.eq.263 .or. my_diode.eq.264) then
         do i = 1, number_bins
            write(*,*) 'GAT', my_diode, ' ', average_tp_section(i), " ", average_rms_section(i), &
                 & " ", consecutive_rms_violations(i), " ", consecutive_rms_violations_value(i)
         end do
     end if 

     deallocate(average_tp_section, average_rms_section, consecutive_rms_violations, consecutive_rms_violations_value)

   end subroutine tp_rms_finder

   subroutine compute_single_fft_chisq(samprate, nu_low, nu_high, sigma0, powspec, chisq)
     implicit none

     real(dp),                intent(in)  :: samprate, sigma0, nu_low, nu_high
     real(sp), dimension(1:), intent(in)  :: powspec
     real(dp),                intent(out) :: chisq

     real(dp)     :: dnu, w
     integer(i4b) :: numsamp, i, j, k, ind1, ind2

     if (sigma0 == 0 .or. all(powspec == 0.d0)) then
        chisq = 0.d0
        return
     end if
     
     numsamp = size(powspec)
     dnu     = ind2freq(2, samprate, numsamp)
     ind1    = nint(nu_low / dnu)
     ind2    = nint(nu_high / dnu)

     ! Compute chi-square
     chisq = 0.d0
     w     = 0.d0
     do i = ind1, ind2
       chisq = chisq + powspec(i) / sigma0**2
       w     = w     + 1
    end do
    chisq = (chisq - w) / sqrt(w)

  end subroutine compute_single_fft_chisq

end module quiet_stat_mod
