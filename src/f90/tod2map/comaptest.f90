program comaptest
  use comap_Lx_mod
  implicit none

  character(len=512) :: filename
  type(lx_struct)    :: data

  filename = "/mn/stornext/d5/comap/testdata/data_20h_1046_lvl1.h5"
  call read_l1_file(filename,data)






  call free_lx_struct(data)
end program
