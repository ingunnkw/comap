module comap_lx_mod
  use healpix_types
  use comap_defs
  use quiet_hdf_mod
  use quiet_fft_mod
  use quiet_utils
  implicit none

  include "mpif.h"

  type lx_struct
     !! The level1 and level2 part, present in all files.
     real(dp)                                        :: mjd_start
     integer(i4b)                                    :: decimation_time, decimation_nu
     integer(i4b)                                    :: scanmode
     real(dp)                                        :: samprate
     real(dp),     allocatable, dimension(:)         :: time
     real(dp),     allocatable, dimension(:)         :: nu          ! (freq)
     real(dp),     allocatable, dimension(:,:)       :: nu_l1       ! (freq, sideband)
     real(sp),     allocatable, dimension(:,:,:)     :: tod         ! (time, freq, detector)
     real(sp),     allocatable, dimension(:,:,:,:)   :: tod_l1      ! (time, freq, sideband, detector)
     real(sp),     allocatable, dimension(:,:)       :: orig_point  ! Hor; (az/el/dk, time)
     integer(i4b), allocatable, dimension(:)         :: scanmode_l1 ! Scanning status
     integer(i4b), allocatable, dimension(:)         :: flag        ! Status flag per time sample

     !! The level3 part, which is only present in level3-files
     integer(i4b)                                   :: coord_sys
     real(dp)                                       :: scanfreq(2), pixsize, point_lim(4)
     real(sp),     allocatable, dimension(:,:)      :: point     ! Gal; (phi/theta/psi,time,mod)

     real(dp),     allocatable, dimension(:)        :: time_gain            ! (time)
     real(sp),     allocatable, dimension(:,:,:)    :: gain                 ! (time, freq, detector)
     real(dp),     allocatable, dimension(:,:)      :: sigma0, alpha, fknee ! (freq, detector)
     real(dp),     allocatable, dimension(:,:,:,:)  :: corr                 ! (freq, freq, detector, detector)
     real(dp)                                       :: stats(ST_NUM)
     real(dp),     allocatable, dimension(:,:,:)    :: det_stats            ! (freq, detector, stat)
     real(dp),     allocatable, dimension(:,:,:)    :: filter_par           ! (freq, detector, param)

  end type lx_struct

contains

  subroutine read_l1_file(filename, data, only_point)
    implicit none
    character(len=*), intent(in)           :: filename
    logical(lgt),     intent(in), optional :: only_point
    type(lx_struct)                        :: data
    type(hdf_file)                         :: file
    integer(i4b)                           :: nsamp, nfreq, ndet, npoint, nsb, ext(7)
    logical(lgt)                           :: all
    all = .true.; if (present(only_point)) all = .not. only_point
    call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "tod", ext)
    nsamp = ext(4); nfreq = ext(3) ; nsb = ext(2); ndet = ext(1)
    call get_size_hdf(file, "point", ext)
    npoint = ext(2)
             allocate(data%time(nsamp))
             allocate(data%orig_point(npoint,nsamp))
    if (all) allocate(data%nu_l1(nsamp,nsb))
    if (all) allocate(data%tod_l1(nsamp,nfreq,nsb,ndet))
    if (all) allocate(data%flag(nsamp))
    if (all) allocate(data%scanmode_l1(nsamp))
    call read_hdf(file, "start_mjd",            data%decimation_time)
    call read_hdf(file, "decimation_time",      data%decimation_time)
    call read_hdf(file, "decimation_nu",        data%decimation_nu)
    call read_hdf(file, "samprate",             data%samprate)
    call read_hdf(file, "time",                 data%time)
    call read_hdf(file, "orig_point",           data%orig_point)
    if (all) call read_hdf(file, "nu",          data%nu_l1)
    if (all) call read_hdf(file, "tod",         data%tod_l1)
    if (all) call read_hdf(file, "flag",        data%flag)
    if (all) call read_hdf(file, "scanmode_l1", data%flag)
    call close_hdf_file(file)
  end subroutine read_l1_file


  subroutine read_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: nsamp, nfreq, ndet, npoint, ext(7)
    call free_lx_struct(data)
    call open_hdf_file(filename, file, "r")
    call get_size_hdf(file, "tod", ext)
    nsamp = ext(1); nfreq = ext(2); ndet = ext(3)
    allocate(data%time(nsamp), data%tod(nsamp,nfreq,ndet), data%flag(nsamp))
    call get_size_hdf(file, "orig_point", ext)
    npoint = ext(1); nsamp = ext(2)
    allocate(data%orig_point(npoint,nsamp), data%nu(nfreq))
    call read_hdf(file, "decimation_time",  data%decimation_time)
    call read_hdf(file, "decimation_nu",    data%decimation_nu)
    call read_hdf(file, "samprate",         data%samprate)
    call read_hdf(file, "time",             data%time)
    call read_hdf(file, "nu",               data%nu)
    call read_hdf(file, "tod",              data%tod)
    call read_hdf(file, "orig_point",       data%orig_point)
    call read_hdf(file, "flag",             data%flag)
    call close_hdf_file(file)
  end subroutine read_l2_file


  ! ! Where should this sub logically be?
  ! subroutine decimate(time, time_full, tod, tod_full, point, point_full, dec)
  !   implicit none
  !   integer(i4b)                                      :: n, nfreq, ndet, npt, dec, i, j, k, l, ind
  !   real(dp), dimension(:),      intent(inout)        :: time, time_full
  !   real(dp), dimension(:,:,:),  intent(inout)        :: tod, tod_full
  !   real(dp), dimension(:,:,0:), intent(inout)        :: point, point_full

  !   if (dec>0) then
  !      n     = size(tod,1)
  !      nfreq = size(tod,2)
  !      ndet  = size(tod,3)
  !      npt   = size(point,1)

  !      ! Averaging over every (dec) elements of the time dimension of time, tod and pointing arrays
  !      ind = 1
  !      do i=1,n
  !         time(i) = mean(time_full(ind:ind+dec-1))
  !         do j=1,ndet
  !            do l=1, nfreq
  !               tod(i,l,j) = mean(tod_full(ind:ind+dec-1,l,j))
  !            end do
  !            do k=1,npt
  !               ! do I need to make the angles safe?
  !               point(k,i,j-1) = mean(point_full(k,ind:ind+dec-1,j-1))
  !            end do
  !         end do
  !         ind = ind + dec
  !      end do

  !   else
  !      time = time_full
  !      tod = tod_full
  !      point = point_full
  !   end if
  ! end subroutine decimate


  subroutine read_l3_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b) :: npoint, nsamp, nfreq, ndet, nmod, ext(7)
    call free_lx_struct(data)
    call read_l2_file(filename, data)
    call open_hdf_file(filename, file, "r")
    call read_hdf(file, "scan_freq",  data%scanfreq)    
    !call read_hdf(file, "pixsize", data%pixsize)
    call read_hdf(file, "point_lim", data%point_lim)
    ! Read pointing
    call get_size_hdf(file, "point", ext)
    npoint = ext(1); nsamp = ext(2)!; mod = ext(3)
    allocate(data%point(npoint,nsamp))
    call read_hdf(file, "point",     data%point)
    !call read_hdf(file, "coord_sys", data%coord_sys)  
    call read_hdf(file, "coord_sys", data%coord_sys)
    ! Read gain
    call get_size_hdf(file, "time_gain", ext)
    nsamp = ext(1)
    allocate(data%time_gain(nsamp))
    call read_hdf(file, "time_gain", data%time_gain)
    call get_size_hdf(file, "gain", ext)
    nsamp = ext(1); nfreq = ext(2); ndet = ext(3)
    allocate(data%gain(nsamp,nfreq,ndet))
    call read_hdf(file, "gain", data%gain)
    ! Read noise parameters
    !call get_size_hdf(file, "corr", ext)
    call get_size_hdf(file, "sigma", ext)
    nfreq = ext(1); ndet = ext(2)
    allocate(data%sigma0(nfreq,ndet),data%alpha(nfreq,ndet),data%fknee(nfreq,ndet))
    !allocate(data%corr(nfreq,nfreq,ndet,ndet))
    !allocate(data%det_stats(nfreq,ndet,NUM_DET_STATS), data%filter_par(nfreq,ndet,NUM_FILTR_PAR))
    allocate(data%det_stats(ndet,nfreq,1), data%filter_par(4,nfreq,ndet))
    call read_hdf(file, "sigma", data%sigma0)
    call read_hdf(file, "alpha",  data%alpha)
    call read_hdf(file, "fknee",  data%fknee)
    !call read_hdf(file, "corr",   data%corr)
    ! Read stats
    !call read_hdf(file, "stats",       data%stats)
    call read_hdf(file, "stats", data%det_stats)
    ! Read filter parameters
    call read_hdf(file, "filter_par",  data%filter_par)
    call close_hdf_file(file)
  end subroutine read_l3_file

  subroutine free_lx_struct(data)
    implicit none
    type(lx_struct) :: data
    if(allocated(data%time))        deallocate(data%time)
    if(allocated(data%tod))         deallocate(data%tod)
    if(allocated(data%tod_l1))      deallocate(data%tod_l1)
    if(allocated(data%orig_point))  deallocate(data%orig_point)

    if(allocated(data%point))       deallocate(data%point)
    if(allocated(data%time_gain))   deallocate(data%time_gain)
    if(allocated(data%gain))        deallocate(data%gain)
    if(allocated(data%sigma0))      deallocate(data%sigma0)
    if(allocated(data%alpha))       deallocate(data%alpha)
    if(allocated(data%fknee))       deallocate(data%fknee)
    if(allocated(data%corr))        deallocate(data%corr)
    if(allocated(data%det_stats))   deallocate(data%det_stats)
    if(allocated(data%filter_par))  deallocate(data%filter_par)
  end subroutine

  subroutine write_l2_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "time",              data%time)
    call write_hdf(file, "decimation_time",   data%decimation_time)
    call write_hdf(file, "decimation_nu",     data%decimation_nu)
    call write_hdf(file, "samprate",          data%samprate)
    call write_hdf(file, "tod",               data%tod)
    call write_hdf(file, "orig_point",        data%orig_point)
    call close_hdf_file(file)
  end subroutine

  subroutine write_l3_file(filename, data)
    implicit none
    character(len=*), intent(in) :: filename
    type(lx_struct)              :: data
    type(hdf_file)               :: file
    integer(i4b)                 :: i

    call open_hdf_file(filename, file, "w")
    call write_hdf(file, "time",              data%time)
    call write_hdf(file, "decimation_time",   data%decimation_time)
    call write_hdf(file, "decimation_nu",     data%decimation_nu)
    call write_hdf(file, "scanfreq",          data%scanfreq)
    call write_hdf(file, "samprate",          data%samprate)
    call write_hdf(file, "coord_sys",         data%coord_sys)
    call write_hdf(file, "pixsize",           data%pixsize)
    call write_hdf(file, "point_lim",         data%point_lim)
    call write_hdf(file, "tod",               data%tod)
    call write_hdf(file, "orig_point",        data%orig_point)
    call write_hdf(file, "point",             data%point)
    call write_hdf(file, "sigma0",            data%sigma0)
    call write_hdf(file, "alpha",             data%alpha)
    call write_hdf(file, "fknee",             data%fknee)
    call write_hdf(file, "corr",              data%corr)
    call write_hdf(file, "time_gain",         data%time_gain)
    call write_hdf(file, "gain",              data%gain)
    call write_hdf(file, "stats",             data%stats)
    call write_hdf(file, "det_stats",         data%det_stats)
    call write_hdf(file, "filter_par",        data%filter_par)
    call close_hdf_file(file)
  end subroutine

  ! Concatenate a set of lx_structs. We assume that:
  ! The structs in in are compatible, and that they are
  ! either level2 or level3 (not something in between).
  ! This means that the level2-only parts always will be allocated,
  ! and that if one level3 part exists, they all do
  subroutine cat_lx_structs(in, out)
    implicit none
    type(lx_struct), intent(in)    :: in(:)
    type(lx_struct), intent(inout) :: out
    integer(i4b)                   :: i, j, jg, n, m, mg, ndet, nmod, npix, ng, nf
    logical(lgt)                   :: l3
    real(dp)                       :: f

    call free_lx_struct(out)

    ndet = size(in(1)%tod,3)
    l3   = allocated(in(1)%point)
    n    = 0
    ng   = 0
    do i = 1, size(in)
      n = n + size(in(i)%time)
      if(l3) ng = ng + size(in(i)%time_gain)
    end do
    allocate(out%time(n), out%tod(n,nf,ndet), out%orig_point(3,n))
    if(l3) then
       allocate(out%point(3,n))!,size(in(1)%point)))
       allocate(out%time_gain(ng), out%gain(ng,nf,ndet))
    end if
    j  = 0
    jg = 0
    mg = 0
    do i = 1, size(in)
       m  = size(in(i)%time)
       out%time       (  1+j:m+j  )   = in(i)%time
       out%tod        (  1+j:m+j,:,:) = in(i)%tod
       out%orig_point (:,1+j:m+j  )   = in(i)%orig_point
       if(l3) then
          mg = size(in(i)%time_gain)
          out%point(:,1+j:m+j)          = in(i)%point
          out%time_gain(1+jg:mg+jg)     = in(i)%time_gain
          out%gain(1+jg:mg+jg,:,:)      = in(i)%gain
       end if
       j  = j+m
       jg = jg + mg
    end do

    ! These have fixed length
    out%decimation_time = in(1)%decimation_time
    out%decimation_nu   = in(1)%decimation_nu
    out%samprate        = in(1)%samprate
    out%coord_sys       = in(1)%coord_sys
    out%scanfreq        = in(1)%scanfreq
    out%pixsize         = in(1)%pixsize
    out%point_lim       = in(1)%point_lim

    ! These could actually differ between the files. We make
    ! a weighted average
    if(l3) then
       allocate(out%sigma0(nf,ndet), out%alpha(nf,ndet), out%fknee(nf,ndet), out%corr(nf,nf,ndet,ndet))
       allocate(out%det_stats(ndet,nf,size(in(1)%det_stats,2)))
       allocate(out%filter_par(ndet,nf,size(in(1)%filter_par,2)))
       out%sigma0 = 0; out%alpha = 0; out%fknee = 0; out%corr = 0
       out%stats  = 0; out%det_stats = 0; out%filter_par = 0
       do i = 1, size(in)
          f = real(size(in(i)%time),dp)/n
          out%sigma0 = out%sigma0 + in(i)%sigma0 * f
          out%alpha  = out%alpha  + in(i)%alpha  * f
          out%fknee  = out%fknee  + in(i)%fknee  * f
          out%corr   = out%corr   + in(i)%corr   * f
          out%stats  = out%stats  + in(i)%stats  * f
          out%det_stats  = out%det_stats  + in(i)%det_stats*f
          out%filter_par = out%filter_par + in(i)%filter_par*f
       end do
    end if

  end subroutine


end module comap_lx_mod

