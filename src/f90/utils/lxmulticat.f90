! Read in a list of old -> new mappings. It has almost the same
! form as a runlist.
! lxmulticat [options] runlist indir outdir
program lxmulticat
  use quiet_utils
  use quiet_lx_mod
  use quiet_fileutils
  use fake_mpi_mod
  implicit none
  character(len=512)              :: arg, idir, odir, ofile, runlist, objname, ifile
  type(lx_struct),    allocatable :: idata(:)
  type(lx_struct)                 :: odata
  integer(i4b)                    :: i, j, n, type, unit, nobj, obj, nsub, sub, iid, oid
  integer(i4b)                    :: verbosity, myid, nproc, err, k
  real(dp)                        :: from, to, az, el, dk, lon, lat

  type      = 3
  verbosity = 0
  do j = 1, iargc()
     call getarg(j, arg)
     if(arg(1:1) /= '-') exit
     select case(arg)
        case("-2"); type = 2
        case("-3"); type = 3
        case("-v"); verbosity = verbosity + 1
        case("-q"); verbosity = verbosity - 1
        case("-h"); call help; stop
        case default
           write(stderr,"(a)") "Unrecognized argument: " // trim(arg)
           stop
     end select
  end do
  if(iargc()-j+1 /= 3) then
     call help; stop
  end if

  call getarg(j+0, runlist)
  call getarg(j+1, idir)
  call getarg(j+2, odir)

  call mpi_init(err)
  call mpi_comm_rank(mpi_comm_world, myid,  err)
  call mpi_comm_size(mpi_comm_world, nproc, err)
  call dset(level=verbosity, id=myid)

  ! We do the reading of the runlist in parallel with the
  ! actual operation. That way, we don't need to create any
  ! data structures for it.
  k = -1
  unit = getlun()
  open(unit,file=runlist,action="read",status="old")
  read(unit,*) nobj
  do obj = 1, nobj
     read(unit,*) objname, nsub
     do sub = 1, nsub
        k = k+1
        read(unit,*) oid, from, to, n, az, el, dk, lon, lat
        if(modulo(k-1,nproc) /= myid) then
           do i = 1, n; read(unit,*); end do
           cycle
        end if
        call dmem("Processing " // trim(itoa(oid)),0)
        allocate(idata(n))
        do i = 1, n
           read(unit,*) iid
           ifile = getname(idir, objname, iid)
           call dmem("Reading " // trim(ifile),1)
           if(type == 2) then; call read_l2_file(ifile, idata(i))
           else;               call read_l3_file(ifile, idata(i)); end if
        end do
        call dmem("Concatenating",2)
        call cat_lx_structs(idata, odata)
        ofile = getname(odir, objname, oid)
        call mkdirs(ofile, .true.)
        call dmem("Writing " // trim(ofile),1)
        if(type == 2) then; call write_l2_file(ofile, odata)
        else;               call write_l3_file(ofile, odata); end if
        do i = 1, n
           call free_lx_struct(idata(i))
        end do
        deallocate(idata)
        call free_lx_struct(odata)
     end do
  end do
  close(unit)

contains

  function getname(dir, obj, id) result(name)
    implicit none
    character(len=*)   :: dir, obj
    character(len=512) :: name
    integer(i4b)       :: id
    name = trim(dir) // "/" // trim(obj) // "/" // &
     & trim(obj) // "_" // trim(itoa(id)) // ".hdf"
  end function

  subroutine help
    implicit none
    write(stderr,'(a)') "lxmulticat: Concatenate many levelX-files into fewer."
    write(stderr,'(a)') "  lxmulticat [options] runlist idir odir."
    write(stderr,'(a)') "Options:"
    write(stderr,'(a)') " -2: Operate on level2-files."
    write(stderr,'(a)') " -3: Operate on level3-files (default)."
  end subroutine

end program
