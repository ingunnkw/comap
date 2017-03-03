! Convert a 1 Hz xy tilt file to a 1/ces omega,theta tilt file.
! This will be small enough to be easily stored in pointing mod.
! We don't really want too many files of this kind, as they kept
! adding up to a large size in total the previous time, but
! automatic reading of level3-files in pointing mod was
! more problematic than this. Note! This introduces a new
! file that must be regenerated when the ces definitions change.
! Another alternative would be to store a greatly downsampled version,
! but I am worried of the edges of the ceses being degraded that way.
! One could also average within ceses, but then actually store everything
! per 5 min or something, making sure that one prefers in-ces samples.
! But for now: ces-based.
program tilt2hdf
  use quiet_utils
  use quiet_hdf_mod
  use quiet_ces_mod
  implicit none
  character(len=512)    :: ifname, ofname, parfile
  type(hdf_file)        :: hfile
  integer(i4b)          :: cnum, r(2), nces, margin
  type(quiet_ces_info)  :: ces
  real(dp)              :: xy(2), x, y, cycx, mjd0, step
  real(dp), allocatable :: idata(:,:), odata(:,:,:)
  call getarg(1, parfile)
  call getarg(2, ifname)
  call getarg(3, ofname)
  call initialize_ces_mod(parfile)
  margin = 10
  nces   = get_num_ces()
  allocate(odata(2,nces,2))
  odata  = 0
  call open_hdf_file(ifname, hfile, "r")
  call read_alloc_hdf(hfile, "tilt", idata)
  call read_hdf(hfile, "t0",   mjd0)
  call read_hdf(hfile, "step", step)
  call close_hdf_file(hfile)
  do cnum = 1, nces
     call get_ces_info(cnum, ces)
     r  = int((ces%mjd-mjd0)/step) + 1 + [margin,-margin]
     r  = max(1,min(size(idata,1),r))
     xy = sum(idata(r(1):r(2),:),1)/(r(2)-r(1)+1)*DEG2RAD/60/60
     if(r(2)-r(1) == 0) xy = 0
     x  = xy(1); y = xy(2)
     ! convert from x,y to omega,theta
     cycx = cos(y*cos(x))
     odata(:,ces%cid,1) = [ -pi/2 - atan2(sin(y*cos(x)),cycx*sin(x)), acos(cos(x)*cycx) ]
     odata(:,ces%cid,2) = xy
     call free_ces_info(ces)
  end do
  call open_hdf_file(ofname, hfile, "w")
  call write_hdf(hfile, "tilt", odata(:,:,1))
  call write_hdf(hfile, "xy",   odata(:,:,2))
  call close_hdf_file(hfile)
  deallocate(idata, odata)
end program
