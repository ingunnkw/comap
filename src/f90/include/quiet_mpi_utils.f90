module quiet_mpi_utils
  use quiet_utils
  use quiet_system_mod
  use quiet_mpi_mod
  implicit none

contains

  subroutine print_host_mapping
    implicit none
    integer(i4b)       :: myid, nproc, err, i, pid, stat(32)
    character(len=512) :: hostname, str
    call mpi_comm_rank(mpi_comm_world, myid,  err)
    call mpi_comm_size(mpi_comm_world, nproc, err)
    call get_hostname(hostname)
    pid = get_pid()
    1 format(a8,' ',a32,' 'a10)
    2 format(i8,' ',a32,' ',i10)
    if(myid == 0) then
       write(*,1) 'Id','Host','Pid'
       write(*,2) myid, trim(hostname), pid
       do i = 1, nproc-1
          call mpi_recv(hostname, len(hostname), mpi_character, i, 0, mpi_comm_world, stat, err)
          call mpi_recv(pid, 1, mpi_integer, i, 0, mpi_comm_world, stat, err)
          write(*,2) i, trim(hostname), pid
       end do
    else
       call mpi_send(hostname, len(hostname), mpi_character, 0, 0, mpi_comm_world, err)
       call mpi_send(pid, 1, mpi_integer, 0, 0, mpi_comm_world, err)
    end if
  end subroutine

end module
