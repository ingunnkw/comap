! l1_mjd: Given a list of level1-files relative to the current cirectory,
! extract the start and end times, and output mjd1 mjd2 file.fits.
! Any corrupt files are simply skipped.
program l1_mjd
  use l1_read_mod
  implicit none

  character(len=512)                            :: prefix, listfile
  character(len=512), dimension(:), allocatable :: filelist
  type(point_struct)                            :: point
  integer(i4b)                                  :: i, unit, status
  call getarg(1, prefix)
  call getarg(2, listfile)
  call read_filelist(listfile, filelist)

  unit = getlun()
  do i = 1, size(filelist)
     call deallocate_point_struct(point)
     call l1_read(unit, trim(prefix) // "/" // trim(filelist(i)), pointing=point, &
      & modules=[0], selector=makesel(only=[point_time]), &
      & status_code=status)
     if(status /= 0) cycle
     write(*,fmt="(2f17.10,a)") point%time([1,size(point%time)])/24/60/60, " " // trim(filelist(i))
  end do
contains
  subroutine read_filelist(listfile, list)
    implicit none
    character(len=*)   :: listfile
    character(len=512) :: line
    character(len=*), allocatable, dimension(:)  :: list
    integer(i4b)     :: i, j, n, unit
    if(allocated(list)) deallocate(list)
    unit = getlun()
    open(unit,file=listfile,status="old",action="read")
    n = 0
    do
       read(unit,*,end=1)
       n = n+1
    end do
1   continue
    rewind(unit)
    allocate(list(n))
    do i = 1, n
       read(unit,fmt="(a)") list(i)
    end do
    close(unit)
  end subroutine
end program l1_mjd
