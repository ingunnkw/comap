! Convert a tilt file to even-spaced hdf format by using splines
program tilt2hdf
  use quiet_utils
  use quiet_hdf_mod
  use spline_1d_mod
  implicit none
  character(len=512)    :: ifname, ofname
  integer(i4b)          :: unit, i, n, nout
  real(dp)              :: foo(4), daysec
  type(spline_type)     :: interpol
  type(hdf_file)        :: hfile
  real(dp), allocatable :: idata(:,:), odata(:,:)
  call getarg(1, ifname)
  call getarg(2, ofname)
  ! Read in the data
  unit = getlun()
  open(unit,file=ifname,action="read",status="old")
  n = 0
  do
     read(unit,*,end=1) foo
     n = n+1
  end do
1 allocate(idata(n,4))
  rewind(unit)
  do i = 1, n
     read(unit,*) idata(i,:)
  end do
  close(unit)

  ! Create regularly sampled time axis
  daysec = 24*60*60
  idata(:,1) = idata(:,1)*daysec + idata(:,2)/1e3
  nout   = int(idata(n,1)-idata(1,1)+1)
  allocate(odata(nout,3))
  do i = 1, nout
     odata(i,1) = idata(1,1) + i-1
  end do

  call spline      (interpol, idata(:,1), idata(:,3))
  call splint_multi(interpol, odata(:,1), odata(:,2))
  call spline      (interpol, idata(:,1), idata(:,4))
  call splint_multi(interpol, odata(:,1), odata(:,3))
  deallocate(idata)
  ! And output as hdf file. Since t is regularly sampled,
  ! we only need to store t0 and interval
  call open_hdf_file(ofname, hfile, "w")
  call write_hdf(hfile, "t0", odata(1,1)/daysec)
  call write_hdf(hfile, "step", 1d0/daysec)
  call write_hdf(hfile, "tilt", real(odata(:,2:3),sp))
  call close_hdf_file(hfile)

end program
