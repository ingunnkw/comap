program bin_test
  use comap_Lx_mod
  implicit none

  character(len=520) :: filename, name2
  type(lx_struct)    :: data
  real(dp)           :: cut, sum
  integer(i4b)       :: i,j, lim, num

  filename = "/mn/stornext/d5/comap/testdata/data_20h_1046_lvl1.h5"
  !write(*,*) trim(filename)
  call read_l1_file(trim(filename), data)
  data%time(:) = (data%time(:)-data%time(1))*24*3600

  cut = 1.2 ! s
  !lim = 500 ! number of bins
  do i = 1, 1000
     if (data%time(i) .ge. cut) then
        num = i ! datapoints in one bin
        exit
     end if
  end do
  lim = size(data%time)/num

  open(42, file='smooth.txt', status='replace', recl=1024)
  do i = 1, num
     sum = 0.d0
     do j = 1, lim
        sum = sum + data%tod_l1(i + (j-1)*num,100,1,1)
     end do
     write(42,*) data%time(i), sum/lim
  end do
  close(42)

end program bin_test
