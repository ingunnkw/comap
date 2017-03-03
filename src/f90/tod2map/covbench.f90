program covbench
  use tod2map_utils
  use tod2map_mapmaker
  implicit none
  type(quiet_assembly) :: assembly
  type(common_info)    :: info
  type(mapinfo)        :: data
  type(mapdata)        :: out
  character(len=512)   :: fname, which
  integer(i4b)         :: bench_tot
  call getarg(1, fname)
  call getarg(2, which)
  call read_cov_indata(assembly, info, data, out, fname)
  allocate(out%cov(out%n+1,info%ncomp,out%col_start:out%col_stop,info%ncomp))
  call setup_bench(info)
  assembly%numsamp = 10000
  call bench_start(info%bench, bench_tot)
  select case(which)
     case("1"); call build_cov(assembly, info, data, out)
     case("2"); call build_cov_test(assembly, info, data, out)
     case("3"); call build_cov_test2(assembly, info, data, out)
     case("4"); call build_cov_pol(assembly, info, data, out)
  end select
  call bench_stop(info%bench, bench_tot)
  call output_benchmarks(info, "/dev/stdout")

contains

  subroutine setup_bench(info)
    implicit none
    type(common_info) :: info
    bench_cov1   = bench_init(info%bench, "cov_rbuf")
    bench_cov2   = bench_init(info%bench, "cov_fill1")
    bench_cov3   = bench_init(info%bench, "cov_fill2")
    bench_cov4   = bench_init(info%bench, "cov_sample")
    bench_tot    = bench_init(info%bench, "total")
  end subroutine

  subroutine output_benchmarks(info, benchfile)
    implicit none
    type(common_info), intent(in) :: info
    character(len=*),  intent(in) :: benchfile
    integer(i4b)                  :: i, j, n, ierr, unit, nsamp(info%bench%n)
    real(dp)                      :: dt(info%bench%n)
    n = info%bench%n
    nsamp = sum(info%bench%nsamp(:,1:n),1)
    dt    = sum(info%bench%dt   (:,1:n),1)
    unit = getlun()
    open(unit,file=benchfile)
    do i = 1, n
       write(unit,'(a16,3e15.7)') info%bench%name(i), dt(i)/nsamp(i), real(nsamp(i),dp), dt(i)
    end do
    close(unit)
  end subroutine
 
end program
