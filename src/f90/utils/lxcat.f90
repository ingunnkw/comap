program lxcat
  use quiet_utils
  use quiet_lx_mod
  implicit none
  character(len=512)              :: arg, ofile
  character(len=512), allocatable :: ifiles(:)
  type(lx_struct),    allocatable :: idata(:)
  type(lx_struct)                 :: odata
  integer(i4b)                    :: i, j, n, type

  ! Get the inputs. Modify this section if options are
  ! to be supported
  type = 3
  do j = 1, iargc()
     call getarg(j, arg)
     if(arg(1:1) /= '-') exit
     select case(arg)
        case("-2"); type = 2
        case("-3"); type = 3
        case("-h"); call help; stop
        case default
           write(stderr,"(a)") "Unrecognized argument: " // trim(arg)
           stop
     end select
  end do
  if(iargc() - j <= 0) then
     call help; stop
  end if

  n = iargc()-j
  allocate(ifiles(n),idata(n))
  do i = 1, n
     call getarg(j+i-1, ifiles(i))
     if(type == 2) then; call read_l2_file(ifiles(i), idata(i))
     else;               call read_l3_file(ifiles(i), idata(i)); end if
  end do
  call getarg(iargc(), ofile)

  call cat_lx_structs(idata, odata)

  if(type == 2) then; call write_l2_file(ofile, odata)
  else;               call write_l3_file(ofile, odata); end if

  do i = 1, n
     call free_lx_struct(idata(i))
  end do
  call free_lx_struct(odata)
  deallocate(idata, ifiles)

contains

  subroutine help
    implicit none
    write(stderr,'(a)') "lxcat: Concatenate level3 (or level2 with the -2 flag) files"
    write(stderr,'(a)') "into one."
    write(stderr,'(a)') "  lxcat [options] ofile ifiles"
    write(stderr,'(a)') "Options:"
    write(stderr,'(a)') " -2: Operate on level2-files."
    write(stderr,'(a)') " -3: Operate on level3-files (default)."
  end subroutine

end program
