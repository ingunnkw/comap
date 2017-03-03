program fsmap
  use quiet_fileutils
  use quiet_pixspace_mod
  use quiet_healpix_mod
  use alm_tools
  use coord_v_convert
  use rngmod
  implicit none

  type stack_item
     real(dp), allocatable :: map(:,:,:)
     real(dp)              :: num
  end type

  character(len=512)        :: cmd, arg, ifile, ofile, ifile2, ofile2
  integer(i4b)              :: nside, nside2, order2, order, i, j, k, n, n2, err
  integer(i4b)              :: ai, nstack, nmap, ncomp, p, p2, r, ext(3), ext2(3)
  integer(i4b)              :: lmax, l
  real(dp)                  :: val, mat(3,3), alpha, beta, gamma
  complex(dpc), allocatable :: alm(:,:,:)
  integer(i4b), allocatable :: pixels(:), inds(:), pixels2(:), inds2(:)
  real(dp),     allocatable :: map(:,:,:), map2(:,:,:), omap(:,:,:)
  integer(i4b), allocatable :: mask(:)
  type(planck_rng)          :: rng
  type(degrade_info)        :: deg
  type(stack_item)          :: stack(64)

  call rand_init(rng, 1)
  call getarg(1, cmd)
  select case(cmd)
     case("sparse")
        call getarg (2, ifile)
        call getargd(3, ofile, ifile)
        call read_map(map, pixels, nside, order, ifile)
        allocate(inds(size(pixels)))
        n = 0
        do i = 1, size(pixels)
           if(all(healnot(reshape(map(i,:,:),[size(map(1,:,:))])))) cycle
           n = n+1
           inds(n) = i
        end do
        call write_map(map(inds(:n),:,:), pixels(inds(:n)), nside, order, ofile)
        deallocate(map, pixels, inds)
     case("delete")
        call getarg (2, ifile)
        call getarg (3, arg); read(arg,*) val
        call getargd(4, ofile, ifile)
        call read_map(map, pixels, nside, order, ifile)
        allocate(inds(size(pixels)))
        n = 0
        do i = 1, size(pixels)
           if(all(map(i,:,:) == val)) cycle
           n = n+1
           inds(n) = i
        end do
        call write_map(map(inds(:n),:,:), pixels(inds(:n)), nside, order, ofile)
        deallocate(map, pixels, inds)
     case("full")
        call getarg (2, ifile)
        call getargd(3, ofile, ifile)
        call read_map(map, order, ifile)
        call write_map(map, order, ofile)
        deallocate(map)
     case("fill")
        call getarg (2, ifile)
        call getarg (3, arg); read(arg,*) val
        call getargd(4, ofile, ifile)
        call read_map(map, order, ifile)
        do i = 1, size(map,3)
           do j = 1, size(map,2)
              do k = 1, size(map,1)
                 if(healnot(map(k,i,j))) map(k,i,j) = val
              end do
           end do
        end do
        call write_map(map, order, ofile)
        deallocate(map)
     !case("nside")
     !   call getarg (2, ifile)
     !   call getarg (3, arg); read(arg,*) nside2
     !   call getargd(4, ofile, ifile)
     !   call read_map(map, pixels, nside, order, ifile)
     !   if(nside2 < nside) then
     !      call prepare_degrade(deg, pixels, order, nside, nside2)
     !      allocate(omap(size(deg%opix),size(map,2),size(map,3)))
     !      do i = 1, size(map,3)
     !         call degrade_map(map(:,:,i), deg, omap(:,:,i))
     !      end do
     !      call write_map(omap, deg%opix, nside2, order, ofile)
     !      call free_degrade_info(deg)
     !   else
     !      r = (nside2/nside)**2
     !      allocate(pixels2(size(pixels)*r), omap(size(pixels)*r,size(map,2),size(map,3)))
     !      do i = 1, size(pixels)
     !         call set_ordering(nside, nest, order, pixels(i), p)
     !         p = p*r
     !         do j = 1, r
     !            call set_ordering(nside2, order, nest, p+j-1, p2)
     !            pixels2((i-1)*r+j) = p2
     !            omap((i-1)*r+j,:,:) = map(i,:,:)
     !         end do
     !      end do
     !      call write_map(omap, pixels2, nside2, order, ofile)
     !      deallocate(pixels2)
     !   end if
     !   deallocate(map, omap, pixels)
     case("order")
        call getarg (2, ifile)
        call getarg (3, arg)
        call getargd(4, ofile, ifile)
        select case(arg)
           case("ring"); order2 = ring
           case("nest"); order2 = nest
           case default; call assert(.false., "Unknown order: " // trim(arg))
        end select
        call read_map(map, pixels, nside, order, ifile)
        do i = 1, size(pixels)
           call set_ordering(nside, order2, order, pixels(i), pixels(i))
        end do
        call write_map(map, pixels, nside, order2, ofile)
        deallocate(map, pixels)
     case("intersect")
        call getarg (2, ifile)
        call getarg (3, ifile2)
        call getarg (4, ofile)
        call getargd(5, ofile2, "")
        call read_map(map,  pixels,  nside,  order,  ifile)
        call read_map(map2, pixels2, nside2, order2, ifile2)
        call assert(nside == nside2, "Inconsistent nside!")
        if(order /= order2) then
           do i = 1, size(pixels2)
              call set_ordering_pix(nside, order, order2, pixels2(i), pixels2(i))
           end do
        end if
        allocate(mask(0:12*nside**2-1))
        mask = 0
        mask(pixels)  = mask(pixels)  + 1
        mask(pixels2) = mask(pixels2) + 1
        call wherei(mask(pixels)  == 2, inds)
        call write_map(map(inds,:,:), pixels(inds), nside, order, ofile)
        deallocate(inds)
        if(ofile2 /= "") then
           call wherei(mask(pixels2) == 2, inds2)
           call write_map(map2(inds2,:,:), pixels2(inds2), nside, order, ofile2)
           deallocate(inds2)
        end if
        deallocate(map, map2, mask, pixels, pixels2)
     case("union")
        n = iargc() - 2
        call assert(n > 0, "Union needs at least one input map!")
        do i = 1, n
           call getarg(i+1, ifile)
           call read_map(map, pixels, nside2, order2, ifile)
           if(i == 1) then
              nside = nside2
              order = order2
              ncomp = size(map,2)
              nmap  = size(map,3)
              allocate(mask(0:12*nside**2-1), omap(0:12*nside**2-1,ncomp,nmap))
              mask = 0
           else
              call assert(nside2 == nside, "Inconsistent nside in map " // trim(itoa(i)))
              call assert(size(map,2) == ncomp, "Inconsistent ncomp in map " // trim(itoa(i)))
              call assert(size(map,3) == nmap, "Inconsistent nmap in map " // trim(itoa(i)))
              if(order /= order2) then
                 do j = 1, size(pixels)
                    call set_ordering_pix(nside, order, order2, pixels(j), pixels(j))
                 end do
              end if
           end if
           mask(pixels)     = 1
           omap(pixels,:,:) = map
           deallocate(pixels, map)
        end do
        call getarg(iargc(), ofile)
        call wherei(mask >= 1, inds)
        call write_map(omap(inds-1,:,:), inds-1, nside, order, ofile)
        deallocate(omap, mask, inds)
     case("map2mask")
        call getarg(2, ifile)
        call getarg(3, ofile)
        call read_map(map, pixels, nside, order, ifile)
        allocate(omap(0:12*nside**2-1,1,1))
        omap = -1
        do i = 1, size(pixels)
           omap(pixels(i),1,1) = i
        end do
        call write_map(omap, order, ofile)
        deallocate(map, pixels, omap)
     case("rot","hrot")
        call getarg(2, ifile)
        i = 3
        ! Two formats: fsmap rot ifile type ofile
        !              fsmap rot ifile alpha beta gamma ofile
        if(iargc() == 4) then
           call getarg(i, arg); i = i+1
           select case(arg)
              case("equ2gal")
                 call coordsys2euler_zyz(2000d0, 2000d0, 'Q', 'G', alpha, beta, gamma)
              case("gal2equ")
                 call coordsys2euler_zyz(2000d0, 2000d0, 'G', 'Q', alpha, beta, gamma)
              case default
                 call assert(.false.,"Unknown transformation " // trim(arg))
           end select
        else
           call getarg(i, arg); i=i+1; read(arg,*) alpha
           call getarg(i, arg); i=i+1; read(arg,*) beta
           call getarg(i, arg); i=i+1; read(arg,*) gamma
        end if
        call getarg(i, ofile)
        if(cmd == "hrot") then
           ! Rotate using spherical harmonics. Sky must be full.
           call read_map(map, order, ifile)
           do j = 1, size(map,3)
              call set_ordering(ring, order, map(:,:,j))
              ncomp = size(map,2); if(healnot(map(1,size(map,2),j))) ncomp=1
              nside = npix2nside(size(map,1)); lmax = 3*nside
              allocate(alm(ncomp,0:lmax,0:lmax))
              if(ncomp == 1) then
                 call map2alm(nside, lmax, lmax, map(:,1,j), alm, [-1d0,1d0], get_hpix_ringweights(nside))
              else
                 call map2alm(nside, lmax, lmax, map(:,:,j), alm, [-1d0,1d0], get_hpix_ringweights(nside))
              end if
              !! Smooth a bit to avoid aliasing
              !do l = 0, lmax
              !   alm(:,l,:) = alm(:,l,:)*exp(-2d0*l*(l+1)/(lmax**2))
              !end do
              call rotate_alm(lmax, alm, alpha, beta, gamma)
              if(ncomp == 1) then
                 call alm2map(nside, lmax, lmax, alm, map(:,1,1))
              else
                 call alm2map(nside, lmax, lmax, alm, map(:,:,j))
              end if
              deallocate(alm)
           end do
           call write_map(map, ring, ofile)
        end if
     case("op")
        ! Simple operations on exactly compatible maps (i.e. all
        ! dimensions equal)
        nstack = 0
        ext    = 0
        do ai = 2, iargc()
           call getarg(ai, arg)
           select case(arg)
              case("+","add");     call binop(stack, nstack, "+")
              case("-","sub");     call binop(stack, nstack, "-")
              case("*","mul");     call binop(stack, nstack, "*")
              case("/","div");     call binop(stack, nstack, "/")
              case("**","^","pow");call binop(stack, nstack, "^")
              case(">","greater"); call binop(stack, nstack, ">")
              case("<","less");    call binop(stack, nstack, "<")
              case("=","==","eq"); call binop(stack, nstack, "=")
              case("%","mod");     call binop(stack, nstack, "%")
              case("#","rand");    call unop (stack, nstack, "#")
              case("!","nod");     call unop (stack, nstack, "!")
              case("dup")
                 call assert(nstack > 0, "Empty stack for dup!")
                 nstack = nstack+1
                 stack(nstack)%num = stack(nstack-1)%num
                 if(allocated(stack(nstack-1)%map)) then
                    ext = shape(stack(nstack-1)%map)
                    allocate(stack(nstack)%map(ext(1),ext(2),ext(3)))
                    stack(nstack)%map = stack(nstack-1)%map
                 end if
              case("pop")
                 call assert(nstack > 0, "Cannot pop empty stack!")
                 call free_stack_item(stack(nstack))
                 nstack = nstack-1
              case("dump")
                 call assert(nstack > 0, "Cannot dump empty stack!")
                 if(allocated(stack(nstack)%map)) then
                    ext = shape(stack(nstack)%map)
                    do i = 1, ext(3)
                       do j = 1, ext(2)
                          do k = 1, ext(1)
                              write(*,'(4i9,e15.7)') i,j,k,pixels(k),stack(nstack)%map(k,j,i)
                          end do
                       end do
                    end do
                 else
                    write(*,'(e15.7)') stack(nstack)%num
                 end if
                 call free_stack_item(stack(nstack))
                 nstack = nstack-1
              case default
                 if(ai == iargc()) then
                    ! Last argument assumed to be output file if it doesn't
                    ! match a command.
                    call assert(nstack >= 1, "Stack is empty at end!")
                    call assert(allocated(stack(nstack)%map), "Top of stack is not a map!")
                    call write_map(stack(nstack)%map, pixels, nside, order, arg)
                 else
                    ! Push number or map on stack
                    nstack = nstack + 1
                    read(arg,*,iostat=err) val
                    if(err == 0) then
                       stack(nstack)%num = val
                    else
                       call read_map(stack(nstack)%map, pixels2, nside2, order2, arg)
                       ext2 = shape(stack(nstack)%map)
                       if(ext(1) == 0) then
                          ext   = ext2
                          nside = nside2
                          order = order2
                          allocate(pixels(ext(1)))
                          pixels = pixels2
                       else
                          call assert(all(ext==ext2),"Incompatible shape: " // trim(arg))
                          call assert(nside==nside2, "Incompatible nside: " // trim(arg))
                          call assert(order==order2, "Incompatible order: " // trim(arg))
                          call assert(all(pixels==pixels2), "Incompatible pixels: " // trim(arg))
                       end if
                       deallocate(pixels2)
                    end if
                 end if
           end select
        end do
     case default
        write(stderr,'(a)') "Unknown command " // trim(cmd)
        write(stderr,*) "Valid commands:"
        write(stderr,*) "fsmap sparse    ifile [ofile]"
        write(stderr,*) "fsmap full      ifile [ofile]"
        write(stderr,*) "fsmap intersect ifile1 ifile2 ofile1 [ofile2]"
        write(stderr,*) "fsmap union     ifile1 ifile2 ofile1 [ofile2]"
  end select

contains

  subroutine getargd(i, str, def)
    implicit none
    integer(i4b)     :: i
    character(len=*) :: str, def
    if(iargc() >= i) then
       call getarg(i, str)
    else
       str = def
    end if
  end subroutine

  subroutine free_stack_item(item)
    implicit none
    type(stack_item) :: item
    if(allocated(item%map)) deallocate(item%map)
  end subroutine

  subroutine binop(stack, n, op)
    implicit none
    type(stack_item) :: stack(:)
    character(len=*) :: op
    integer(i4b)     :: i, j, k, n, ext(3)
    call assert(n >= 2, "Too few arguments to " // trim(op))
    if(allocated(stack(n-1)%map)) then
       ext = shape(stack(n-1)%map)
       if(allocated(stack(n)%map)) then
          do i = 1, ext(3); do j = 1, ext(2); do k = 1, ext(1)
             if(healnot(stack(n-1)%map(k,j,i)) .or. healnot(stack(n)%map(k,j,i))) then
                stack(n-1)%map(k,j,i) = hpx_dbadval
             else
                stack(n-1)%map(k,j,i) = binop_raw(stack(n-1)%map(k,j,i),stack(n)%map(k,j,i),op)
             end if
          end do; end do; end do
       else
          do i = 1, ext(3); do j = 1, ext(2); do k = 1, ext(1)
             if(healnot(stack(n-1)%map(k,j,i)) .or. healnot(stack(n)%num)) then
                stack(n-1)%map(k,j,i) = hpx_dbadval
             else
                stack(n-1)%map(k,j,i) = binop_raw(stack(n-1)%map(k,j,i),stack(n)%num,op)
             end if
          end do; end do; end do
       end if
    else
       if(allocated(stack(n)%map)) then
          ext = shape(stack(n)%map)
          allocate(stack(n-1)%map(ext(1),ext(2),ext(3)))
          do i = 1, ext(3); do j = 1, ext(2); do k = 1, ext(1)
             if(healnot(stack(n-1)%num) .or. healnot(stack(n)%map(k,j,i))) then
                stack(n-1)%map(k,j,i) = hpx_dbadval
             else
                stack(n-1)%map(k,j,i) = binop_raw(stack(n-1)%num,stack(n)%map(k,j,i),op)
             end if
          end do; end do; end do
       else
           if(healnot(stack(n-1)%num) .or. healnot(stack(n)%num)) then
              stack(n-1)%num = hpx_dbadval
           else
              stack(n-1)%num = binop_raw(stack(n-1)%num, stack(n)%num, op)
           end if
       end if
     end if
     call free_stack_item(stack(n))
     n = n-1
  end subroutine

  subroutine unop(stack, n, op)
    implicit none
    type(stack_item) :: stack(:)
    character(len=*) :: op
    integer(i4b)     :: i, j, k, n, ext(3)
    call assert(n >= 1, "Too few arguments to " // trim(op))
    if(allocated(stack(n)%map)) then
       ext = shape(stack(n)%map)
       do i = 1, ext(3); do j = 1, ext(2); do k = 1, ext(1)
          if(healok(stack(n)%map(k,j,i))) stack(n)%map(k,j,i) = unop_raw(stack(n)%map(k,j,i),op)
       end do; end do; end do
    else
       if(healok(stack(n)%num)) stack(n)%num = unop_raw(stack(n)%num,op)
    end if
  end subroutine

  function binop_raw(a, b, op) result(c)
    implicit none
    real(dp)         :: a, b, c
    character(len=1) :: op
    select case(op)
       case("+");  c = a+b
       case("-");  c = a-b
       case("*");  c = a*b
       case("/");  c = a/b
       case("^");  c = a**b
       case(">");  c = ifeli(a> b,1,0)
       case("<");  c = ifeli(a< b,1,0)
       case("=");  c = ifeli(a==b,1,0)
       case("%");  c = modulo(a,b)
       case default
          call assert(.false., "Unknown operation " // trim(op))
    end select
  end function

  function unop_raw(a, op) result(b)
    implicit none
    real(dp)         :: a, b
    character(len=1) :: op
    select case(op)
       case("#");  b = rand_gauss(rng)*a
       case("!");  b = ifeli(a==0,1,0)
       case default
          call assert(.false., "Unknown operation " // trim(op))
    end select
  end function

end program
