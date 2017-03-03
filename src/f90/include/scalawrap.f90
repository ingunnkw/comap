module scalawrap
  use healpix_types
  use quiet_mpi_mod
  implicit none

  type scinfo
     integer(i4b) :: context, comm, myid, nproc, gridsize(2), gridpos(2)
     logical(lgt) :: active
  end type scinfo

  type scalamat
     real(dp), dimension(:,:), allocatable :: data
     integer(i4b) :: desc(9), rglob, cglob, rloc, cloc, rblock, cblock, rorig, corig
     ! rmap, cmap map from local indices to global indices.
     ! rown, cown map from global indices to processors
     ! rpos, cpos map from global indices to local indices
     integer(i4b), dimension(:), allocatable :: rmap, cmap, rown, cown, rpos, cpos
     type(scinfo) :: info
  end type scalamat

contains

  subroutine ownership_basic(output, length, procnum, bsize, start)
    implicit none
    integer :: length, procnum, bsize, i, j, origin, output(:)
    integer, optional :: start
    origin = 0
    if(present(start)) origin = start
    do i = 0, length-1
      output(i+1) = mod(i/bsize+origin,procnum)
    end do
  end subroutine ownership_basic

  function localpos(i, nproc, bsize) result(j)
    implicit none
    integer :: i, j, nproc, bsize
    integer :: bz
    bz = bsize
    j = mod(i-1,bz) + (i-1)/(bz*nproc)*bz + 1
  end function

  subroutine full_ownership(output, desc, dir)
    implicit none
    integer :: output(:), desc(9), dir, d, procn(2), proc(2)
    d = dir-1
    call blacs_gridinfo(desc(2), procn(1), procn(2), proc(1), proc(2))
    call ownership_basic(output, desc(3+d), procn(d+1), desc(5+d), desc(7+d))
  end subroutine full_ownership

  subroutine extract_matches(full,into,whom)
    implicit none
    integer(i4b) :: full(:), whom, i, j
    integer(i4b), dimension(:), allocatable :: into
    if(.not. allocated(into)) allocate(into(count(full==whom)))
    j = 1
    do i = 1, size(full)
       if(full(i) /= whom) cycle
       into(j)=i
       j=j+1
    end do
  end subroutine

  ! Given a blacs descriptor of a distributed matric, output
  ! which rows/columns (dir=1/2) we own into the array output,
  ! which must have size loc_r/loc_c.
  subroutine ownership(output, desc, dir)
    implicit none
    integer :: desc(9), procn(2), proc(2), output(:), d, dir,n,i,j
    integer, dimension(:), allocatable :: owners
    d = dir-1
    call blacs_gridinfo(desc(2), procn(1), procn(2), proc(1), proc(2))
    allocate(owners(desc(3+d)))
    call full_ownership(owners, desc, dir)
    j = 1
    do i = 1, size(owners)
      if(owners(i) .eq. proc(dir)) then
        output(j) = i
        j = j+1
      end if
    end do
    if(j .eq. 1 .and. size(output) >=1) output(1) = 0
    deallocate(owners)
  end subroutine ownership

  subroutine proc2grid(procn, nprow, npcol)
    implicit none
    integer :: procn, nprow, npcol, bestnrow, bestprocn, i, n
    nprow = int(sqrt(real(procn,dp)))
    bestprocn = 0
    do i = nprow, max(1,int(sqrt(real(procn,dp))/2)), -1
       n = procn/i*i
       if(n > bestprocn) then; bestnrow = i; bestprocn = n; end if
    end do
    nprow = bestnrow
    npcol = procn/nprow
  end subroutine proc2grid

  subroutine sc_init(info, startcomm, active)
    implicit none
    type(scinfo) :: info
    integer(i4b), optional :: startcomm
    logical(lgt), optional :: active
    integer(i4b) :: comm, myid, ierr, numprocs, context,rc
    logical(lgt) :: redundant

    if(present(startcomm)) then
       comm = startcomm
    else
       comm = MPI_COMM_WORLD
       call mpi_init(ierr)
       if(ierr /= 0) call sc_handle_error("mpi_init", ierr)
    end if

    call mpi_comm_rank(comm, info%myid, ierr)
    if(ierr /= 0) call sc_handle_error("mpi_comm_rank", ierr)
    call mpi_comm_size(comm, info%nproc, ierr)
    if(ierr /= 0) call sc_handle_error("mpi_comm_size", ierr)
    call proc2grid(info%nproc, info%gridsize(1), info%gridsize(2))
    call sl_init(info%context, info%gridsize(1), info%gridsize(2))
    call blacs_gridinfo(info%context, info%gridsize(1), info%gridsize(2), info%gridpos(1), info%gridpos(2))
    redundant = any(info%gridpos .ge. info%gridsize)
    if(redundant) then; rc=1; else; rc=0; end if
    ! Redundant processes should not be included
    call mpi_comm_split(comm, rc, 0, info%comm, info)
    call mpi_comm_rank(info%comm, info%myid, ierr)
    if(ierr /= 0) call sc_handle_error("mpi_rank_size", ierr)
    call mpi_comm_size(info%comm, info%nproc,ierr)
    if(ierr /= 0) call sc_handle_error("mpi_comm_size", ierr)
    info%active = .not. redundant

    if(present(active)) then
       if(active .and. redundant) call sc_stop(info)
    end if
  end subroutine

  subroutine sc_stop(info)
    implicit none
    type(scinfo) :: info
    integer(i4b) :: ierr, i, j, rank
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_finalize(ierr)
    if(ierr /= 0) call sc_handle_error("mpi_finalize", ierr)
    stop
  end subroutine

  subroutine sc_alloc(mat, n, m, info, blocksize, owner, origin)
    implicit none
    type(scalamat) :: mat
    type(scinfo)   :: info
    integer(i4b), optional :: blocksize, owner, origin(2)
    integer(i4b) :: i,n, m, bsize, nproc(2), procid(2),ierr, numroc

    bsize = 32; mat%rorig = 0; mat%corig = 0; ierr = 0
    if(present(blocksize)) bsize = blocksize
    if(present(origin)) then; mat%rorig = origin(1); mat%corig = origin(2); end if
    if(present(owner)) bsize = 2147483647
    mat%rblock = bsize
    mat%cblock = bsize
    mat%info = info
    call blacs_gridinfo(info%context, nproc(1), nproc(2), procid(1), procid(2))
    if(present(owner)) call blacs_pcoord(info%context, owner, mat%rorig, mat%corig)
    mat%rglob = n; mat%cglob = m
    mat%rloc = numroc(n, bsize, procid(1), mat%rorig, nproc(1))
    mat%cloc = numroc(m, bsize, procid(2), mat%corig, nproc(2))
    call descinit(mat%desc, n, m, bsize, bsize, mat%rorig, mat%corig, info%context, max(1,mat%rloc), ierr)
    if(ierr /= 0) call sc_handle_error("descinit", ierr)
    allocate(mat%rown(n))
    allocate(mat%cown(m))
    call full_ownership(mat%rown, mat%desc, 1)
    call full_ownership(mat%cown, mat%desc, 2)

    allocate(mat%rmap(mat%rloc))
    allocate(mat%cmap(mat%cloc))
    call ownership(mat%rmap, mat%desc, 1)
    call ownership(mat%cmap, mat%desc, 2)

    allocate(mat%rpos(n))
    allocate(mat%cpos(m))
    do i = 1, n; mat%rpos(i) = localpos(i, nproc(1), bsize); end do
    do i = 1, m; mat%cpos(i) = localpos(i, nproc(2), bsize); end do

    allocate(mat%data(mat%rloc,mat%cloc))
  end subroutine

  subroutine sc_dealloc(mat)
    implicit none
    type(scalamat) :: mat
    if(allocated(mat%rown)) deallocate(mat%rown)
    if(allocated(mat%cown)) deallocate(mat%cown)
    if(allocated(mat%rmap)) deallocate(mat%rmap)
    if(allocated(mat%cmap)) deallocate(mat%cmap)
    if(allocated(mat%rpos)) deallocate(mat%rpos)
    if(allocated(mat%cpos)) deallocate(mat%cpos)
    if(allocated(mat%data)) deallocate(mat%data)
  end subroutine

  ! Read a columnwise fortran unformatted matrix
  subroutine sc_read(unit, mat,      broadcast, type)
    implicit none
    integer(i4b)   :: unit
    type(scalamat) :: mat
    integer(i4b), optional :: type
    logical(lgt), optional :: broadcast
    logical(lgt)   :: bcast
    bcast = .false.; if(present(broadcast)) bcast = broadcast
    if(bcast) then
        call sc_read_broadcast(unit,mat, type=type)
    else
        call sc_read_direct(unit,mat, type=type)
    end if
  end subroutine

  subroutine sc_read_broadcast(unit, mat, type)
    implicit none
    integer(i4b)     :: unit, i, j, k, n, m, ncol, t
    type(scalamat)   :: mat
    integer(i4b), optional :: type
    real(dp), dimension(:,:), allocatable :: arr
    real(sp), dimension(:),   allocatable :: sarr
    if(.not. mat%info%active) return
    t = dp; if(present(type)) t = type
    n = mat%rglob; m = mat%cglob
    ncol = 1
    allocate(arr(n,ncol),sarr(n))
    do i = 1, m, ncol
       do j = 1, min(m,i+ncol-1)-i+1
          if(mat%info%myid==0) then
             if(t == sp) then
                read(unit) sarr; arr(:,j) = sarr
             else
                read(unit) arr(:,j)
             end if
          else
             read(unit)
          end if
       end do
       call sc_set(mat, arr, 0, rect=(/1,i,n,min(m,i+ncol-1)/))
    end do
    deallocate(arr, sarr)
  end subroutine

  subroutine sc_read_direct(unit, mat, type)
    implicit none
    integer(i4b)     :: unit, i, j, n, m, t
    type(scalamat)   :: mat
    integer(i4b), optional :: type
    real(dp), dimension(:), allocatable :: arr
    real(sp), dimension(:), allocatable :: sarr
    if(.not. mat%info%active) return
    t = dp; if(present(type)) t = type
    n = mat%rglob; m = mat%cglob
    allocate(arr(n),sarr(n))
    do i = 1, m
      if(mat%cown(i) == mat%info%gridpos(2)) then
         if(t == sp) then
            read(unit) sarr; arr = sarr
         else
            read(unit) arr
         end if
         do j = 1, mat%rloc
            mat%data(j,mat%cpos(i)) = arr(mat%rmap(j))
         end do
      else
        read(unit)
      end if
    end do
    deallocate(arr, sarr)
  end subroutine

  ! Like sc_read, but reads an (n+1,n) combination of       ! AAAA
  ! two symmetrical matrices, A and B, which are expanded   ! BAAA
  ! internally.                                             ! BBAA
  subroutine sc_read_combo(unit, A, B,      broadcast)      ! BBBA
    implicit none                                           ! BBBB
    integer(i4b)   :: unit
    type(scalamat) :: A, B
    logical(lgt), optional :: broadcast
    logical(lgt)   :: bcast
    bcast = .false.; if(present(broadcast)) bcast = broadcast
    bcast = .true.
    if(bcast) then
       call sc_read_combo_broadcast(unit,A,B)
    else
       call sc_read_combo_direct(unit,A,B)
    end if
  end subroutine

  ! The boss reads by column, and broadcasts to
  ! the owner
  subroutine sc_read_combo_broadcast(unit, A, B)
    implicit none
    integer(i4b)     :: unit, i, j, k, n, m
    type(scalamat)   :: A, B
    real(dp), dimension(:), allocatable :: arr
    if(.not. A%info%active) return
    n = A%rglob; m = A%cglob
    allocate(arr(n+1))
    do i = 1, m
       if(A%info%myid==0) then
          read(unit) arr
       else
          read(unit)
       end if
       ! Now divide this into an A and B part.
       call sc_set(A, reshape(arr(1:i),    [i,    1]),   0, rect=(/1,i,i,i/))
       call sc_set(A, reshape(arr(1:i-1),  [1,    i-1]), 0, rect=(/i,1,i,i-1/))
       call sc_set(B, reshape(arr(i+1:n+1),[n+1-i,1]),   0, rect=(/i,i,n,i/))
       call sc_set(B, reshape(arr(i+2:n+1),[1,    n-i]), 0, rect=(/i,i+1,i,n/))
    end do
    deallocate(arr)
  end subroutine

  ! Everybody reads, and they keep only the parts they want.
  subroutine sc_read_combo_direct(unit, A, B)
    implicit none
    integer(i4b)     :: unit, i, j, k, n, m
    type(scalamat)   :: A, B
    real(dp), dimension(:), allocatable :: arr
    if(.not. A%info%active) return
    n = A%rglob; m = A%cglob
    allocate(arr(n+1))
    do i = 1, m
       read(unit) arr
       ! Get the part of A we own
       do j = 1, i
          if(all(A%info%gridpos == [ A%rown(j), A%cown(i) ])) &
           & A%data(A%rpos(j),A%cpos(i)) = arr(j)
          if(all(A%info%gridpos == [ A%rown(i), A%cown(j) ])) &
           & A%data(A%rpos(i),A%cpos(j)) = arr(j)
       end do
       ! Get the part of B we own
       do j = i, n
          if(all(B%info%gridpos == [ B%rown(j), B%cown(i) ])) &
           & B%data(B%rpos(j),B%cpos(i)) = arr(j+1)
          if(all(B%info%gridpos == [ B%rown(i), B%cown(j) ])) &
           & B%data(B%rpos(i),B%cpos(j)) = arr(j+1)
       end do
    end do
    deallocate(arr)
  end subroutine

  !subroutine sc_read_combo_direct(unit, A, B)
  !  implicit none
  !  integer(i4b)     :: unit, i, j, n, m
  !  type(scalamat)   :: A, B
  !  real(dp), dimension(:), allocatable :: arr
  !  if(.not. mat%info%active) return
  !  n = mat%rglob; m = mat%cglob
  !  allocate(arr(n))
  !  do i = 1, m
  !    if(mat%cown(i) == mat%info%gridpos(2)) then
  !       read(unit) arr
  !       do j = 1, mat%rloc
  !          mat%data(j,mat%cpos(i)) = arr(mat%rmap(j))
  !       end do
  !    else
  !      read(unit)
  !    end if
  !  end do
  !end subroutine

  subroutine sc_write(unit, mat)
    implicit none
    integer(i4b)     :: unit, i, j, k, n, m, ncol
    type(scalamat)   :: mat
    real(dp), dimension(:,:), allocatable :: arr
    if(.not. mat%info%active) return
    n = mat%rglob; m = mat%cglob
    ncol = 1
    allocate(arr(n,ncol))
    do i = 1, m, ncol
       call sc_get(mat, arr, 0, rect=(/1,i,n,min(m,i+ncol-1)/))
       if(mat%info%myid==0) then
          do j = 1, min(m,i+ncol-1)-i+1
             write(unit) arr(:,j)
          end do
       end if
    end do
    deallocate(arr)
  end subroutine

  ! Cooperative subroutine, accesses the global matrix transparently
  ! Must be called by all participants. Only gets the value for one
  ! at a time.
  subroutine sc_get_entry(mat, i, j, res, receiver)
    implicit none
    type(scalamat) :: mat
    integer(i4b)   :: i, j, me, owner, blacs_pnum, receiver, stat(8), ierr
    real(dp)       :: res, value
    if(.not. mat%info%active) return
    ierr = 0
    me = mat%info%myid
    owner = blacs_pnum(mat%info%context, mat%rown(i), mat%cown(j))

    if(receiver .ne. me .and. owner .ne. me) return
    if(receiver .eq. me) then
       if(owner .eq. me) then
          res = mat%data(mat%rpos(i),mat%cpos(j))
       else
          call mpi_recv(res, 1, mpi_double_precision, owner, 0, mat%info%comm, stat, ierr)
          if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
       end if
    else
       value = mat%data(mat%rpos(i),mat%cpos(j))
       call mpi_send(value,1, mpi_double_precision, receiver, 0, mat%info%comm, ierr)
       if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
    end if
  end subroutine

  subroutine sc_set_entry(mat, i, j, value, sender)
    implicit none
    type(scalamat) :: mat
    integer(i4b)   :: i, j, me, owner, blacs_pnum, sender, stat(8), ierr
    real(dp)       :: value
    if(.not. mat%info%active) return
    ierr = 0
    me = mat%info%myid
    owner = blacs_pnum(mat%info%context, mat%rown(i), mat%cown(j))

    if(sender .ne. me .and. owner .ne. me) return
    if(sender .eq. me) then
       if(owner .eq. me) then
          mat%data(mat%rpos(i),mat%cpos(j)) = value
       else
          call mpi_send(value, 1, mpi_double_precision, owner, 0, mat%info%comm, ierr)
          if(ierr /= 0) call sc_handle_error("mpi_send", ierr)
       end if
    else
       call mpi_recv(value,1, mpi_double_precision, sender, 0, mat%info%comm, stat, ierr)
       mat%data(mat%rpos(i),mat%cpos(j)) = value
       if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
    end if
  end subroutine

  ! Copy a rectangle of data from the global matrix mat to the
  ! local matrix res at receiver. This is not written very efficiently.
  subroutine sc_get(mat, res, receiver, rect, reverse)
    implicit none
    type(scalamat) :: mat
    integer(i4b)   :: prow, pcol, me,blacs_pnum, receiver, n, m, rectangle(4)
    integer(i4b)   :: i, j, k, l, trow, tcol, sender, stat(5), ierr
    integer(i4b), optional :: rect(4)
    logical(lgt), optional :: reverse
    logical(lgt)   :: backwards
    real(dp)       :: res(:,:)
    integer(i4b), dimension(:), allocatable :: rkey, ckey
    real(dp), dimension(:,:), allocatable :: buffer
    if(.not. mat%info%active) return
    me = mat%info%myid
    rectangle = (/1,1,mat%rglob,mat%cglob/)
    if(present(rect)) rectangle = rect
    backwards = .false.; if(present(reverse)) backwards = reverse

    if(me == receiver) then
       do prow = 0, mat%info%gridsize(1)-1
          n = count(mat%rown(rectangle(1):rectangle(3)) == prow)
          if(n < 1) cycle
          allocate(rkey(n))
          j = rectangle(1); do i = 1, n; do while(mat%rown(j) /= prow); j=j+1; end do; rkey(i) = j; j=j+1; end do
          do pcol = 0, mat%info%gridsize(2)-1
             m = count(mat%cown(rectangle(2):rectangle(4)) == pcol)
             if(m < 1) cycle
             allocate(ckey(m))
             j = rectangle(2); do i = 1, m; do while(mat%cown(j) /= pcol); j=j+1; end do; ckey(i) = j; j=j+1; end do
             sender = blacs_pnum(mat%info%context, prow, pcol)
             if(mat%info%myid == sender) then
                if(backwards) then
                   do l = 1, size(ckey); do k = 1, size(rkey)
                      mat%data(mat%rpos(rkey(k)),mat%cpos(ckey(l))) = res(rkey(k)-rectangle(1)+1,ckey(l)-rectangle(2)+1)
                   end do; end do
                else
                   do l = 1, size(ckey); do k = 1, size(rkey)
                      res(rkey(k)-rectangle(1)+1,ckey(l)-rectangle(2)+1) = mat%data(mat%rpos(rkey(k)),mat%cpos(ckey(l)))
                   end do; end do
                end if
             else
                allocate(buffer(n,m))
                if(backwards) then
                   do l = 1, size(ckey); do k = 1, size(rkey)
                      buffer(k,l) = res(rkey(k)-rectangle(1)+1,ckey(l)-rectangle(2)+1)
                   end do; end do
                   call mpi_send(buffer, n*m, mpi_double_precision, sender, 0, mat%info%comm, ierr)
                   if(ierr /= 0) call sc_handle_error("mpi_send", ierr)
                else
                   call mpi_recv(buffer, n*m, mpi_double_precision, sender, 0, mat%info%comm, stat, ierr)
                   if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
                   do l = 1, size(ckey); do k = 1, size(rkey)
                      res(rkey(k)-rectangle(1)+1,ckey(l)-rectangle(2)+1) = buffer(k,l)
                   end do; end do
                end if
                deallocate(buffer)
             end if
             deallocate(ckey)
          end do
          deallocate(rkey)
       end do
    else
       call blacs_pcoord(mat%info%context, receiver, trow, tcol)
       prow = mat%info%gridpos(1)
       pcol = mat%info%gridpos(2)
       n = count(mat%rown(rectangle(1):rectangle(3)) == prow)
       m = count(mat%cown(rectangle(2):rectangle(4)) == pcol)
       if(minval((/n,m/))>=1) then
          allocate(rkey(n))
          allocate(ckey(m))
          j = rectangle(1); do i = 1, n; do while(mat%rown(j) /= prow); j=j+1; end do; rkey(i) = j; j=j+1; end do
          j = rectangle(2); do i = 1, m; do while(mat%cown(j) /= pcol); j=j+1; end do; ckey(i) = j; j=j+1; end do
          allocate(buffer(n,m))
          if(backwards) then
             call mpi_recv(buffer, n*m, mpi_double_precision, receiver, 0, mat%info%comm, stat, ierr)
             if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
             do l = 1, size(ckey); do k = 1, size(rkey)
                mat%data(mat%rpos(rkey(k)),mat%cpos(ckey(l))) = buffer(k,l)
             end do; end do
          else
             do l = 1, size(ckey); do k = 1, size(rkey)
                buffer(k,l) = mat%data(mat%rpos(rkey(k)),mat%cpos(ckey(l)))
             end do; end do
             call mpi_send(buffer, n*m, mpi_double_precision, receiver, 0, mat%info%comm, ierr)
             if(ierr /= 0) call sc_handle_error("mpi_send", ierr)
          end if
          deallocate(buffer)
          deallocate(rkey,ckey)
       end if
    end if
  end subroutine

  subroutine sc_set(mat, from, sender, rect)
    implicit none
    type(scalamat) :: mat
    integer(i4b)   :: sender
    integer(i4b), optional :: rect(4)
    real(dp), dimension(:,:) :: from
    call sc_get(mat, from, sender, rect, reverse=.true.)
  end subroutine

  ! Slow, but handles slicing
  subroutine sc_get2(A, res, who, cols, rows, rect, which)
    implicit none
    type(scalamat) :: A, B
    real(dp)       :: res(:,:)
    integer(i4b)   :: who
    integer(i4b), optional :: cols(:), rows(:), rect(4), which
    call sc_alloc(B, size(res(:,1)), size(res(1,:)), A%info, owner=who)
    call sc_copy(A, B, rows, cols, rect, which)
    if(A%info%myid == who) res = B%data
    call sc_dealloc(B)
  end subroutine

  ! Copies (a part of) global matrix A into (a part of) global matrix B.
  ! With only A and B given, performs a full copy. If rows, cols or rect
  ! is given, copies that subset of A into B. If "which" is specified,
  ! the slicing applies to A if which is 1, B if it is 2. This gives
  ! A handy interface to both sc_get_slice and sc_set_slice.
  ! Rect format: UL, DR.
  subroutine sc_copy(A, B, rows, cols, rect, which)
    implicit none
    type(scalamat) :: A, B
    integer(i4b), optional :: rows(:), cols(:), rect(4), which
    integer(i4b), dimension(:), allocatable :: r, c
    integer(i4b) :: i, j, rc(4), n, which_
    ! Plain copy if not slicing
    if(.not. present(rows) .and. .not. present(cols) .and. .not. present(rect)) then
       call sc_copy_plain(A,B)
       return
    end if
    which_ = 1; if(present(which)) which_ = which
    call sc_range_helper(r, A%rglob, 1, inds=rows, rect=rect)
    call sc_range_helper(c, A%cglob, 2, inds=cols, rect=rect)
    if(which_ == 2) then
       call sc_copy_slice(B, A, r, c, reverse=.true.)
    else
       call sc_copy_slice(A, B, r, c)
    end if
    deallocate(r,c)
  end subroutine

  subroutine sc_copy_plain(A, B, uplo)
    implicit none
    type(scalamat) :: A, B
    character(len=1), optional :: uplo
    character(len=1) :: uplo_
    uplo_ = "A"; if(present(uplo)) uplo_ = uplo
    call pdlacpy(uplo_, A%rglob, A%cglob, A%data, 1, 1, A%desc, &
        & B%data, 1, 1, B%desc)
  end subroutine

  subroutine sc_copy_slice(A, B, rows, cols, reverse)
    implicit none
    type(scalamat) :: A, B
    integer(i4b)   :: rows(:), cols(:), nproc
    integer(i4b)   :: me, i, j
    logical(lgt), optional :: reverse
    logical(lgt)   :: done(0:A%info%nproc-1,0:A%info%nproc-1), scratch(0:A%info%nproc-1,0:A%info%nproc-1)
    logical(lgt)   :: send
    if(.not. A%info%active) return
    me = A%info%myid; nproc = A%info%nproc
    done = .false.; scratch = .false.
    do while(.not. all(done))
       scratch = done
       do i = 0, nproc-1
          do j = i, nproc-1
             if(.not. scratch(i,j)) then
                done(i,j) = .true.
                done(j,i) = .true.
                scratch(:,[i,j]) = .true.
                scratch([i,j],:) = .true.
                call sc_get_crumb(A, B, rows, cols, i, j, reverse)
                if(i/=j) call sc_get_crumb(A, B, rows, cols, j, i, reverse)
             end if
          end do
       end do
    end do
  end subroutine

  ! Handle the part of the rows,cols slice of A into B that involves
  ! sending from sender to receiver. Probably only makes sense when
  ! called by sc_get_slice. The reverse argument reverses all the
  ! data transfers themselves, but leaves the rest of the logic unchanged.
  ! This is used to implement sc_set_slice.
  subroutine sc_get_crumb(A, B, rows, cols, sender, receiver, reverse)
    implicit none
    type(scalamat) :: A, B
    integer(i4b)   :: rows(:), cols(:), sender, receiver
    integer(i4b), dimension(:), allocatable :: brmap, bcmap, arows, acols
    integer(i4b)   :: stat(5), ierr, me, rrow, rcol, srow, scol
    logical(lgt), optional :: reverse
    logical(lgt)   :: backwards
    integer(i8b)   :: s, i, j
    real(dp), dimension(:,:), allocatable :: buffer
    me = A%info%myid; ierr = 0
    if(me /= sender .and. me /= receiver) return
    backwards = .false.; if(present(reverse)) backwards = reverse
    call blacs_pcoord(A%info%context, receiver, rrow, rcol)
    call blacs_pcoord(A%info%context, sender,   srow, scol)

    ! Find the parts of B which belong in receiver
    call extract_matches(B%rown, brmap, rrow)
    call extract_matches(B%cown, bcmap, rcol)

    ! Find the parts of these that can be found in A
    call extract_matches(A%rown(rows(brmap)), arows, srow)
    call extract_matches(A%cown(cols(bcmap)), acols, scol)

    ! Ok, so we will copy into B%data(B%rpos(brmap(arows)),B%cpos(bcmap(acols)))
    ! from A%data(A%rpos(rows(brmap(arows))),A%cpos(cols(bcmap(acols))))
    s = int(size(arows),8)*int(size(acols),8)
    if(s > 0) then
       if(sender == receiver) then
          if(backwards) then
             do j = 1, size(acols); do i = 1, size(arows)
                A%data(A%rpos(rows(brmap(arows(i)))),A%cpos(cols(bcmap(acols(j))))) = &
                     & B%data(B%rpos(brmap(arows(i))),B%cpos(bcmap(acols(j))))
             end do; end do
             ! This logically equivalent statement caused stack overflow
!             A%data(A%rpos(rows(brmap(arows))),A%cpos(cols(bcmap(acols)))) = &
!              & B%data(B%rpos(brmap(arows)),B%cpos(bcmap(acols)))

          else
             do j = 1, size(acols); do i = 1, size(arows)
                B%data(B%rpos(brmap(arows(i))),B%cpos(bcmap(acols(j)))) = &
                  & A%data(A%rpos(rows(brmap(arows(i)))),A%cpos(cols(bcmap(acols(j)))))
             end do; end do
             ! This logically equivalent statement caused stack overflow
             !B%data(B%rpos(brmap(arows)),B%cpos(bcmap(acols))) = &
             !  & A%data(A%rpos(rows(brmap(arows))),A%cpos(cols(bcmap(acols))))
          end if
       elseif(sender == me) then
          allocate(buffer(size(arows),size(acols)))
          if(backwards) then
             call mpi_recv(buffer, size(buffer), mpi_double_precision, receiver, 0, A%info%comm, stat, ierr)
             if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
             A%data(A%rpos(rows(brmap(arows))),A%cpos(cols(bcmap(acols)))) = buffer
          else
             buffer = A%data(A%rpos(rows(brmap(arows))),A%cpos(cols(bcmap(acols))))
             call mpi_send(buffer, size(buffer), mpi_double_precision, receiver, 0, A%info%comm, ierr)
             if(ierr /= 0) call sc_handle_error("mpi_send", ierr)
          end if
          deallocate(buffer)
       else
          allocate(buffer(size(arows),size(acols)))
          if(backwards) then
             do j = 1, size(acols); do i = 1, size(arows)
                buffer(i,j) = B%data(B%rpos(brmap(arows(i))),B%cpos(bcmap(acols(j))))
             end do; end do
             !buffer = B%data(B%rpos(brmap(arows)),B%cpos(bcmap(acols)))
             call mpi_send(buffer, size(buffer), mpi_double_precision, sender, 0, A%info%comm, ierr)
             if(ierr /= 0) call sc_handle_error("mpi_send", ierr)
          else
             call mpi_recv(buffer, size(buffer), mpi_double_precision, sender, 0, A%info%comm, stat, ierr)
             if(ierr /= 0) call sc_handle_error("mpi_recv", ierr)
             do j = 1, size(acols); do i = 1, size(arows)
                B%data(B%rpos(brmap(arows(i))),B%cpos(bcmap(acols(j)))) = buffer(i,j)
             end do; end do
             !B%data(B%rpos(brmap(arows)),B%cpos(bcmap(acols))) = buffer
          end if
          deallocate(buffer)
       end if
    end if
    deallocate(brmap, bcmap, arows, acols)
  end subroutine

  subroutine sc_range_helper(res, n, dim,   inds, rect)
    implicit none
    integer(i4b) :: n, dim, rc(2), i, j
    integer(i4b), optional :: inds(:), rect(4)
    integer(i4b), allocatable, dimension(:) :: res
    rc = (/1,n/); if(present(rect)) rc = rect((/dim,2+dim/))
    ! Prepare initial inds and columns
    if(present(inds)) then
       allocate(res(count(inds >= rc(1) .and. inds <= rc(2))))
       j = 1; do i = 1, size(inds)
          if(inds(i) < rc(1) .or. inds(i) > rc(2)) cycle
          res(j) = inds(i); j=j+1
       end do
    else
       allocate(res(rc(2)-rc(1)+1))
       do i = 1, size(res); res(i) = rc(1)+i-1; end do
    end if
  end subroutine


  ! Solve equation system Ax = b, where b is the input value of x.
  subroutine sc_solve(A,x,  symmetric,oerr)
    implicit none
    type(scalamat) :: A, x
    integer(i4b), dimension(:), allocatable :: pivot
    logical(lgt), optional :: symmetric
    logical(lgt) :: symm
    integer(i4b) :: ierr, i, j, k, l
    integer(i4b), optional :: oerr
    if(.not. A%info%active) return
    symm = .false.; if(present(symmetric)) symm = symmetric
    if(symm) then
       call pdposv("L", A%rglob, x%cglob, A%data, 1, 1, A%desc, x%data, 1, 1, x%desc, ierr)
       if(ierr /= 0) call sc_handle_error("pdposv (in sc_solve)", ierr, oerr)
    else
       allocate(pivot(A%rloc+A%rblock))
       call pdgesv(A%rglob, x%cglob, A%data, 1, 1, A%desc, pivot, x%data, &
         & 1, 1, x%desc, ierr)
       if(ierr /= 0) call sc_handle_error("pdgesv (in sc_solve)", ierr, oerr)
       deallocate(pivot)
    end if
  end subroutine

  subroutine sc_invert(A, symmetric,oerr)
    implicit none
    type(scalamat) :: A
    integer(i4b), dimension(:), allocatable :: pivot, iwork
    real(dp),     dimension(:), allocatable :: work
    logical(lgt), optional :: symmetric
    logical(lgt) :: symm
    integer(i4b) :: ierr, i, j, k, l, iwsize
    integer(i4b), optional :: oerr
    real(dp)     :: wsize
    if(.not. A%info%active) return
    symm = .false.; if(present(symmetric)) symm = symmetric
    if(symm) then
       call pdpotrf("L", A%rglob, 1, 1, A%desc, ierr)
       if(ierr /= 0) call sc_handle_Error("pdpotrf (in sc_invert)", ierr, oerr)
       call pdpotri("L", A%rglob, 1, 1, A%desc, ierr)
       if(ierr /= 0) call sc_handle_Error("pdpotri (in sc_invert)", ierr, oerr)
    else
       allocate(pivot(A%rloc+A%rblock))
       call pdgetrf(A%rglob, A%cglob, A%data, 1, 1, A%desc, pivot, ierr)
       if(ierr /= 0) call sc_handle_Error("pdgetrf (in sc_invert)", ierr, oerr)
       ! Query work sizes
       call pdgetri(A%rglob, A%data, 1, 1, A%desc, pivot, wsize, -1, iwsize, -1, ierr)
       if(ierr /= 0) call sc_handle_Error("pdgetri (in sc_invert)", ierr, oerr)
       allocate(work(int(wsize)),iwork(iwsize))
       call pdgetri(A%rglob, A%data, 1, 1, A%desc, pivot, work, int(wsize), iwork, iwsize, ierr)
       if(ierr /= 0) call sc_handle_Error("pdgetri (in sc_invert[2])", ierr, oerr)
       deallocate(pivot, work, iwork)
    end if
  end subroutine

  subroutine sc_cholesky_decompose(A,uplo,status)
    implicit none
    type(scalamat) :: A
    integer(i4b), optional  :: status
    integer(i4b)  :: ierr
    character(len=1), optional :: uplo
    character(len=1) :: uplo_
    if(.not. A%info%active) return
    uplo_ = "L"; if(present(uplo)) uplo_ = uplo
    call pdpotrf(uplo_, A%rglob, A%data, 1, 1, A%desc, ierr)
    if (present(status)) then
       status = ierr
    else
       if(ierr /= 0) call sc_handle_error("pdpotrf (in sc_cholesky_decompose)", ierr)
    end if
  end subroutine

  subroutine sc_cholesky_invert_matrix(A,uplo,L)
    implicit none
    type(scalamat) :: A
    integer(i4b)  :: ierr
    character(len=1), optional :: uplo
    logical(lgt),     optional :: L
    logical(lgt)     :: L_
    character(len=1) :: uplo_
    uplo_ = "L"; if(present(uplo)) uplo_ = uplo
    L_ = .false.; if (present(L)) L_ = L
    if (.not. L_) call pdpotrf(uplo_, A%rglob, A%data, 1, 1, A%desc, ierr)
    call pdpotri(uplo_, A%rglob, A%data, 1, 1, A%desc, ierr)
    if(ierr /= 0) call sc_handle_error("pdpotrf (in sc_cholesky_invert_matrix)", ierr)
  end subroutine

  subroutine sc_cholesky_solve(A,x, uplo)
    implicit none
    type(scalamat) :: A, x
    integer(i4b)  :: ierr
    character(len=1), optional :: uplo
    character(len=1) :: uplo_
    if(.not. A%info%active) return
    uplo_ = "L"; if(present(uplo)) uplo_ = uplo
    call pdpotrs(uplo_, A%rglob, 1, A%data, 1, 1, A%desc, x%data, 1,1, x%desc,ierr)
    if(ierr /= 0) call sc_handle_error("pdpotrs (in sc_cholesky_solve)", ierr)
  end subroutine

  function sc_cholesky_logdet(A)
    implicit none
    type(scalamat) :: A
    real(dp)       :: sc_cholesky_logdet
    integer(i4b)  :: ierr, i, j
    real(dp) :: my_det, det
    my_det = 0.d0
    if (size(A%data) > 0) then
       do i = 1, A%rloc
          do j = 1, A%cloc
             if (A%rmap(i) == A%cmap(j)) then
                my_det = my_det + 2.d0 * log(A%data(i,j))
             end if
          end do
       end do
    end if
    call mpi_allreduce(my_det, det, 1, MPI_DOUBLE_PRECISION, MPI_SUM, A%info%comm, ierr)
    sc_cholesky_logdet = det
  end function sc_cholesky_logdet

  ! Turn A into the eigenvectors of itself, and fill M with the eigenvalues.
  subroutine sc_eigenvalue_decompose(A,M,uplo)
    implicit none
    type(scalamat) :: A, M, eigvect
    integer(i4b)   :: ierr, n
    real(dp), dimension(:), allocatable :: eigval, work
    real(dp)   :: worksize
    character(len=1), optional :: uplo
    character(len=1) :: uplo_
    if(.not. A%info%active) return
    n    = A%rglob
    uplo_ = "L"; if(present(uplo)) uplo_ = uplo
    allocate(eigval(n))
    call sc_alloc(eigvect, n, n, A%info, blocksize=A%rblock)
    ! First get the work size
    call pdsyev("V", uplo_, n, A%data, 1, 1, A%desc, eigval, eigvect%data, &
      & 1, 1, eigvect%desc, worksize, -1, ierr)
    if(ierr /= 0) call sc_handle_error("pdsyev (in sc_eigenvalue_decompose))", ierr)
    allocate(work(int(worksize)))
    call pdsyev("V", uplo_, n, A%data, 1, 1, A%desc, eigval, eigvect%data, &
      & 1, 1, eigvect%desc, work, int(worksize), ierr)
    if(ierr /= 0) call sc_handle_error("pdsyev (in sc_eigenvalue_decompose))", ierr)
    deallocate(work)
    call sc_copy(eigvect, A)
    call sc_allset(M, reshape(eigval,(/n,1/)))
    deallocate(eigval)
    call sc_dealloc(eigvect)
  end subroutine

  ! 'handle' errors by exiting
  subroutine sc_handle_error(routine, ierr,   oerr)
    implicit none
    integer(i4b)     :: ierr
    integer(i4b), optional :: oerr
    character(len=*) :: routine
    if(present(oerr)) then
       oerr = ierr
    else
       write(*,*) routine // " failed with error code ", ierr
       stop
    end if
  end subroutine

  !subroutine sc_transpose(A,B)  !!! funker ikke
  !  implicit none
  !  type(scalamat) :: A, C
  !  type(scalamat), optional :: B
  !  integer(i4b)   :: n, m
  !  if(.not. A%info%active) return
  !  n = A%rglob; m = A%cglob
  !  if(present(B)) then
  !     call pdtran(m, n, 1d0, A%data, 1, 1, A%desc, 0, B%data, 1, 1, B%desc)
  !  else
  !     ! This in-place version is surprisingly inefficient.
  !     ! I should be able to implement a faster one directly.
  !     call sc_alloc(C, m, n, A%info)
  !     call pdtran(m, n, 1d0, A%data, 1, 1, A%desc, &
  !       & 0, C%data, 1, 1, C%desc)
  !     if(n /= m) then
  !        call sc_dealloc(A)
  !        call sc_alloc(A, m, n, C%info)
  !     end if
  !     call sc_copy(C,A)
  !     call sc_dealloc(C)
  !  end if
  !end subroutine

  function sc_trace_mat_prod(A, B) ! Stupid implementation
    implicit none

    type(scalamat) :: A, B
    real(dp)       :: sc_trace_mat_prod

    integer(i4b) :: i, j, n, ierr
    real(dp)     :: tr
    real(dp), allocatable, dimension(:,:) :: row, col

    n  = A%rglob
    allocate(row(1,n), col(n,1))
    tr = 0.d0
    do i = 1, n
       call sc_get(A, row, 0, rect=[i,1,i,n])
       call sc_get(B, col, 0, rect=[1,i,n,i])
       if (A%info%myid == 0) then
          do j = 1, n
             tr = tr + row(1,j)*col(j,1)
          end do
       end if
    end do
    call mpi_bcast(tr, 1, MPI_DOUBLE_PRECISION, 0, A%info%comm, ierr)
    sc_trace_mat_prod = tr
    deallocate(row, col)
  end function sc_trace_mat_prod

  function sc_trace_mat_symm_prod(A, B)
    implicit none
    type(scalamat) :: A, B
    real(dp)       :: sc_trace_mat_symm_prod
    integer(i4b)   :: ierr
    real(dp)       :: my_tr, tr
    my_tr = sum(A%data*B%data)
    call mpi_allreduce(my_tr, tr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, A%info%comm, ierr)    
    sc_trace_mat_symm_prod = tr
  end function sc_trace_mat_symm_prod

  subroutine sc_rank1_update(A, alpha, x, y)
    implicit none
    type(scalamat) :: A, x, y
    real(dp)       :: alpha
    integer(i4b)   :: n, m
    if(.not. A%info%active) return
    n = A%rglob; m = A%cglob
    call pdger(n, m, alpha, x%data, 1, 1, x%desc, 1, y%data, 1, 1, y%desc, 1, A%data, 1, 1, A%desc)
  end subroutine

  function sc_dotprod(A,B) result(res)
    implicit none
    type(scalamat) :: A, B
    real(dp)       :: res
    integer(i4b)   :: n
    if(.not. A%info%active) return
    n = A%rglob*A%cglob
    call pddot(n, res, A%data, 1, 1, A%desc, 1, B%data, 1, 1, B%desc, 1)
  end function

  subroutine sc_matmul(A,B,C,    transa,transb,alpha,beta,symm)
    implicit none
    type(scalamat) :: A, B, C
    character(len=1), optional :: transa, transb, symm
    real(dp),         optional :: alpha, beta
    character(len=1) :: tra, trb
    real(dp)     :: alp, bet
    integer(i4b) :: m, n, k
    if(.not. A%info%active) return
    tra = 'N'; trb = 'N'; alp = 1d0; bet = 0d0
    if(present(transa)) tra = transa
    if(present(transb)) trb = transb
    m = A%rglob; n = B%cglob; k = A%cglob
    if(tra == 't' .or. tra == 'T') then; m = A%cglob; k = A%rglob; end if
    if(trb == 't' .or. trb == 'T') then; n = B%rglob; end if
    if(present(alpha))  alp = alpha
    if(present(beta))   bet = beta
    if (present(symm)) then
       call pdsymm(symm, 'l', m, n,  alp, A%data, 1, 1, A%desc, B%data, 1, 1, B%desc, bet, C%data, 1, 1, C%desc)
    else
       call pdgemm(tra, trb, m, n, k, alp, A%data, &
            & 1, 1, A%desc, B%data, 1, 1, B%desc, bet, C%data, 1, 1, C%desc)
    end if
  end subroutine

  ! Puts diagonal of matrix into vector diag
  subroutine sc_get_diag(A, diag)
    implicit none
    type(scalamat) :: A
    integer(i4b)   :: m, n, i, ierr
    real(dp)       :: diag(:)
    if(.not. A%info%active) return
    m = A%rglob; 
    n = A%cglob;
    if (n /= m) then
       write(*,*) 'Matrix not n x n. Quiting'
       return
    end if
     do i = 1, n
       call sc_get_entry(A, i, i, diag(i), 0)
       call mpi_bcast(diag(i), 1, MPI_DOUBLE_PRECISION, 0, A%info%comm, ierr)
    end do
  end subroutine

  ! Multiplies the matrix A by diag(b), where b is an n by 1 matrix
  subroutine sc_matmul_diag(A,b,    C,flip)
    implicit none
    type(scalamat) :: A, b
    type(scalamat), optional :: C
    logical(lgt),   optional :: flip
    logical(lgt)   :: flip_
    real(dp), dimension(:,:), allocatable :: arr, tarr
    integer(i4b)   :: i
    flip_ = .false.; if(present(flip)) flip_ = flip

    if(flip_) then
       allocate(arr(A%rloc,1))
       if(b%cglob > b%rglob) then
          allocate(tarr(1,A%rloc))
          call sc_allget(b,tarr, cols=A%rmap)
          arr = reshape(tarr, (/A%rloc,1/))
          deallocate(tarr)
       else; call sc_allget(b,arr, rows=A%rmap)
       end if
       if (present(C))then
          do i = 1, A%rloc; C%data(i,:) = A%data(i,:) * arr(i,1); end do
       else
          do i = 1, A%rloc; A%data(i,:) = A%data(i,:) * arr(i,1); end do
       end if
    else
       allocate(arr(A%cloc,1))
       if(b%cglob > b%rglob) then
          allocate(tarr(1,A%rloc))
          call sc_allget(b,reshape(arr,(/1,A%cloc/)), cols=A%cmap)
          arr = reshape(tarr, (/A%rloc,1/))
          deallocate(tarr)
       else; call sc_allget(b,arr, rows=A%cmap)
       end if
       if (present(C))then
          do i = 1, A%cloc; C%data(:,i) = A%data(:,i) * arr(i,1); end do
       else
          do i = 1, A%cloc; A%data(:,i) = A%data(:,i) * arr(i,1); end do
       end if
    end if
    deallocate(arr)
  end subroutine

  ! More powerful version of sc_get. All processes, not just one,
  ! get the data, and the slices can be different for all.
  ! The overhead is natrually bigger than that of sc_get.
  ! Method: Loop through all pairs. Each pair tells each other what
  ! parts they want, and then do the transfer. This is very similar
  ! to sc_copy.
  subroutine sc_allget(A,b,   rows,cols,rect,reverse)
    implicit none
    type(scalamat) :: A
    real(dp)       :: b(:,:)
    integer(i4b), optional   :: rows(:), cols(:), rect(4)
    integer(i4b), dimension(:), allocatable:: r, c
    integer(i4b)   :: me, i, j, nproc
    logical(lgt), optional :: reverse
    logical(lgt)   :: done(0:A%info%nproc-1,0:A%info%nproc-1), scratch(0:A%info%nproc-1,0:A%info%nproc-1)
    logical(lgt)   :: send
    if(.not. A%info%active) return
    call sc_range_helper(r, A%rglob, 1, inds=rows, rect=rect)
    call sc_range_helper(c, A%cglob, 2, inds=cols, rect=rect)
    me = A%info%myid; nproc = A%info%nproc
    ! This loop sets up communication between all pairs.
    done = .false.; scratch = .false.
    do while(.not. all(done))
       scratch = done
       do i = 0, nproc-1
          do j = i, nproc-1
             if(.not. scratch(i,j)) then
                done(i,j) = .true.
                done(j,i) = .true.
                scratch(:,[i,j]) = .true.
                scratch([i,j],:) = .true.
                call sc_allget_crumb(A, b, r, c, i, j, reverse)
                if(i/=j) call sc_allget_crumb(A, b, r, c, j, i, reverse)
             end if
          end do
       end do
    end do
    deallocate(r,c)
  end subroutine

  ! Performs no communication. Still has some overhead due to
  ! sharing implementation with sc_allget.
  subroutine sc_allset(A,b,rows,cols,rect)
    implicit none
    type(scalamat) :: A
    real(dp)       :: b(:,:)
    integer(i4b), optional   :: rows(:), cols(:), rect(4)
    call sc_allget(A,b,rows,cols,rect,reverse=.true.)
  end subroutine

  subroutine sc_allget_crumb(A, b, rows, cols, sender, receiver, reverse)
    implicit none
    type(scalamat) :: A
    real(dp)       :: b(:,:)
    integer(i4b)   :: rows(:), cols(:), sender, receiver
    integer(i4b), dimension(:), allocatable :: arows, acols, brows, bcols
    integer(i4b)   :: stat(5), ierr, me, rrow, rcol, srow, scol, sizes(2)
    logical(lgt), optional :: reverse
    logical(lgt)   :: backwards
    integer(i8b)   :: s, i, j
    real(dp), dimension(:,:), allocatable :: buffer
    me = A%info%myid; ierr = 0
    if(me /= sender .and. me /= receiver) return
    backwards = .false.; if(present(reverse)) backwards = reverse
    call blacs_pcoord(A%info%context, receiver, rrow, rcol)
    call blacs_pcoord(A%info%context, sender,   srow, scol)

    ! Reverse means that we will be writing to the distributed matrix.
    ! Simply doing the reverse of a get is wasteful here, as we will
    ! be writing the same thing over and over again. In reality, no
    ! communication is needed at all, so perform only the sender ==
    ! receiver part.
    if(backwards .and. .not. (sender == receiver)) return

    ! Special case of when we are talking to ourselves.
    if(sender == receiver) then
       if(backwards) then
          do j = 1, size(cols)
             if(A%cown(cols(j)) /= rcol) cycle
             do i = 1, size(rows)
                if(A%rown(rows(i)) /= rrow) cycle
                A%data(A%rpos(rows(i)),A%cpos(cols(j))) = b(i,j)
             end do
          end do
       else
          do j = 1, size(cols)
             if(A%cown(cols(j)) /= rcol) cycle
             do i = 1, size(rows)
                if(A%rown(rows(i)) /= rrow) cycle
                b(i,j) = A%data(A%rpos(rows(i)),A%cpos(cols(j)))
             end do
          end do
       end if
       return
    end if

    ! Sender must know what parts receiver wants
    if(me == receiver) then
       ! Calculate which parts of b sender has
       call extract_matches(A%rown(rows), brows, srow)
       call extract_matches(A%cown(cols), bcols, scol)
       allocate(arows(size(brows)),acols(size(bcols)))
       arows = rows(brows); acols = cols(bcols)
       ! Tell sender how many rows and cols we want, and which.
       sizes = (/size(brows),size(bcols)/)
       call mpi_send(sizes, 2, mpi_integer, sender, 0, A%info%comm, ierr)
       if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_send", ierr)
       call mpi_send(arows, size(arows), mpi_integer, sender, 0, A%info%comm, ierr)
       if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_send", ierr)
       call mpi_send(acols, size(acols), mpi_integer, sender, 0, A%info%comm, ierr)
       if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_send", ierr)
       ! And then receive that from the sender
       allocate(buffer(size(brows),size(bcols)))
       if(backwards) then
          do j = 1, size(bcols); do i = 1, size(brows)
             buffer(i,j) = b(brows(i),bcols(j))
          end do; end do
          call mpi_send(buffer, size(buffer), mpi_double_precision, sender, 0, A%info%comm, ierr)
          if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_send", ierr)
       else
          call mpi_recv(buffer, size(buffer), mpi_double_precision, sender, 0, A%info%comm, stat, ierr)
          if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_recv", ierr)
          do j = 1, size(bcols); do i = 1, size(brows)
             b(brows(i),bcols(j)) = buffer(i,j)
          end do; end do
       end if
       deallocate(buffer, arows, brows, acols, bcols)
    else
       ! Get what to send from the receiver
       call mpi_recv(sizes, 2, mpi_integer, receiver, 0, A%info%comm, stat, ierr)
       if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_recv", ierr)
       allocate(arows(sizes(1)),acols(sizes(2)))
       call mpi_recv(arows, size(arows), mpi_integer, receiver, 0, A%info%comm, stat, ierr)
       if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_recv", ierr)
       call mpi_recv(acols, size(acols), mpi_integer, receiver, 0, A%info%comm, stat, ierr)
       if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_recv", ierr)
       ! And send it
       allocate(buffer(size(arows),size(acols)))
       if(backwards) then
          call mpi_recv(buffer,size(buffer), mpi_double_precision, receiver, 0, A%info%comm, stat, ierr)
          if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_recv", ierr)
          do j = 1, size(acols); do i = 1, size(arows)
             A%data(A%rpos(arows(i)),A%cpos(acols(j))) = buffer(i,j)
          end do; end do
       else
          do j = 1, size(acols); do i = 1, size(arows)
             buffer(i,j) = A%data(A%rpos(arows(i)),A%cpos(acols(j)))
          end do; end do
          call mpi_send(buffer,size(buffer), mpi_double_precision, receiver, 0, A%info%comm, ierr)
          if(ierr /= 0) call sc_handle_error("sc_getall_crumb: mpi_send", ierr)
       end if
       deallocate(buffer, arows, acols)
    end if
  end subroutine

  ! Own implementation of transpose, since scalapack's version seems to be
  ! broken (we get strange crashes, at least). Transposes A into B.
  ! B must not be the same matrix as A!
  subroutine sc_transpose(A,B)
    implicit none
    type(scalamat)           :: A, B
    type(scinfo)             :: info
    integer(i4b), dimension(:), allocatable :: abcols, abrows, bacols, barows
    real(dp),     dimension(:), allocatable :: rbuf, cbuf
    integer(i4b)   :: me, other, oploc(2), i, j, k, l, nproc, ierr, stat(5)
    logical(lgt)   :: done(0:A%info%nproc-1,0:A%info%nproc-1), scratch(0:A%info%nproc-1,0:A%info%nproc-1)
    info = A%info
    if(.not. info%active) return
    if(A%rglob /= B%cglob .or. A%cglob /= B%rglob) call sc_handle_error("sc_transpose: incompatible arrays.", 1)
    me = info%myid; nproc = info%nproc
    ! This loop sets up communication between all pairs.
    ! Something like this is needed for all all-to-all subroutines, so
    ! it is irritating that it is so hard to factor out things like
    ! this in fortran.
    done = .false.; scratch = .false.
    do while(.not. all(done))
       scratch = done
       do i = 0, nproc-1
          do j = i, nproc-1
             if(.not. scratch(i,j)) then
                done(i,j) = .true.
                done(j,i) = .true.
                scratch(:,[i,j]) = .true.
                scratch([i,j],:) = .true.
                if(me /= i .and. me /= j) cycle
                other = i+j-me
                call blacs_pcoord(info%context, other, oploc(1), oploc(2))

                call extract_matches(B%cown(A%rmap), abrows, oploc(2))
                call extract_matches(B%rown(A%cmap), abcols, oploc(1))
                call extract_matches(A%cown(B%rmap), barows, oploc(2))
                call extract_matches(A%rown(B%cmap), bacols, oploc(1))
                allocate(cbuf(size(abrows)), rbuf(size(bacols)))
                if(i == j) then
                   ! Special case: communicating with self. We will not necessarily
                   ! be distributed symmetrically around the diagonal, so we still
                   ! need to calculate what to transfer. This is the only part
                   ! that should break if A and B are the same array, though it has
                   ! not been tested that the rest works in this case.
                   do k = 1, size(abrows)
                      do l = 1, size(abcols)
                         B%data(B%rpos(A%cmap(abcols(l))), B%cpos(A%rmap(abrows(k)))) = A%data(abrows(k),abcols(l))
                      end do
                   end do
                elseif(me == i) then
                   do k = 1, size(abcols)
                      cbuf = A%data(abrows,abcols(k))
                      call mpi_send(cbuf, size(cbuf), mpi_double_precision, other, 0, info%comm, ierr)
                      if(ierr /= 0) call sc_handle_error("sc_transpose: mpi_send A", ierr)
                   end do
                   do k = 1, size(barows)
                      call mpi_recv(rbuf, size(rbuf), mpi_double_precision, other, 0, A%info%comm, stat, ierr)
                      if(ierr /= 0) call sc_handle_error("sc_transpose: mpi_recv A", ierr)
                      B%data(barows(k), bacols) = rbuf
                   end do
                else
                   do k = 1, size(barows)
                      call mpi_recv(rbuf, size(rbuf), mpi_double_precision, other, 0, A%info%comm, stat, ierr)
                      if(ierr /= 0) call sc_handle_error("sc_transpose: mpi_recv B", ierr)
                      B%data(barows(k), bacols) = rbuf
                   end do
                   do k = 1, size(abcols)
                      cbuf = A%data(abrows,abcols(k))
                      call mpi_send(cbuf, size(cbuf), mpi_double_precision, other, 0, info%comm, ierr)
                      if(ierr /= 0) call sc_handle_error("sc_transpose: mpi_send B", ierr)
                   end do
                end if
                deallocate(cbuf, rbuf, abrows, abcols, barows, bacols)
             end if
          end do
       end do
    end do
  end subroutine

end module
