program template_fit
  use quiet_fileutils
  use quiet_pixspace_mod
  use rngmod
  implicit none

  type template_type
     real(dp),     allocatable :: map(:)
     integer(i4b), allocatable :: pixels(:)
     integer(i4b)              :: nside, order
  end type

  type param_type
     real(dp) :: val, step, min, max, next
  end type

  integer(i4b),    parameter   :: nmax = 256
  type(template_type)          :: imaps(nmax)
  integer(i4b)                 :: iown(2,nmax), i, j, k, l, m, n, nmap, nsingle, comp, dcomp, nt
  integer(i4b)                 :: nside, order, npar, seed, nit, nok, thin
  real(dp)                     :: lik(2)
  type(param_type)             :: params(nmax)
  character(len=512)           :: arg, toks(2), fnames(2)
  character(len=64)            :: step_toks(nmax), tok
  type(planck_rng)             :: rng
  real(dp),        allocatable :: maps(:,:,:), tmap(:,:)
  integer(i4b),    allocatable :: pixels(:)

  iown    = 0
  nmap    = 0
  nsingle = 0
  dcomp   = 1
  seed    = 1
  thin    = 10
  i       = 0
  ! The parameters are [amps...,offset]. The first map has an amplitude of 1
  ! by definition. These defaults are not very good, and should be overridden
  params%val = 0; params%step = 0.1; params%min = -infinity; params%max = infinity
  params(1)%val = -1; params(1)%step = 0
  call dmem("Start")
  ! Read the maps and parse the arguments
  do while(i < iargc())
     i = i+1
     call getarg(i, arg)
     select case(arg)
        case("-seed")
           i = i+1
           call getarg(i, arg)
           read(arg,*) seed
        case("-step")
           i = i+1
           call getarg(i, arg)
           call get_tokens(arg, ",", step_toks, nt)
           do j = 1, nt
              read(step_toks(j),*) params(j+1)%step
           end do
        case("-start")
           i = i+1
           call getarg(i, arg)
           call get_tokens(arg, ",", step_toks, nt)
           do j = 1, nt
              read(step_toks(j),*) params(j+1)%val
           end do
        case("-prior")
           i = i+1
           call getarg(i, arg)
           call get_tokens(arg, ",", toks, n)
           do j = 1, n
              call get_tokens(toks(j), "&", step_toks, nt)
              do k = 1, nt
                 tok = step_toks(k)
                 select case(tok(1:1))
                    case("<"); read(tok(2:),*) params(j+1)%max
                    case(">"); read(tok(2:),*) params(j+1)%min
                    case default
                       read(tok,*) params(j+1)%val
                       params(j+1)%step = 0
                 end select
              end do
           end do
        case default
           if(arg(1:1) == "-") call help("Unrecognized option: " // trim(arg))
           call get_tokens(arg, ":", fnames, n)
           call assert(n == 1 .or. n == 2, "Too many noise specifiers: " // trim(arg))
           nsingle = nsingle + 1
           do j = 1, n
              nmap    = nmap + 1
              call get_tokens(fnames(j), "[]", toks, nt)
              call assert(nt > 0, "Error parsing: " // trim(fnames(j)))
              comp = dcomp; if(nt > 1) read(toks(2),*) comp
              call read_map(tmap, imaps(nmap)%pixels, imaps(nmap)%nside, imaps(nmap)%order, toks(1))
              call dmem("Read " // trim(toks(1)))
              allocate(imaps(nmap)%map(size(tmap,1)))
              imaps(nmap)%map = tmap(:,comp)
              deallocate(tmap)
              iown(j,nsingle) = nmap
           end do
     end select
  end do
  ! Prepare them for use
  call make_compatible(imaps(:nmap), iown(:,:nsingle), maps, pixels, nside, order)
  maps(:,:,2) = maps(:,:,2)**2
  call dmem("Prepared maps")

  ! Clean up input maps
  do i = 1, nmap
     deallocate(imaps(i)%map, imaps(i)%pixels)
  end do

  ! Ok, at this point we should have all our numbers in a consistent format.
  call rand_init(rng, seed)
  npar = nsingle+1
  lik  = infinity
  nit  = 0; nok = 0
  do
     do k = 1, thin
        nit = nit + 1
        ! Generate new suggestion
        do i = 1, npar
           params(i)%next = params(i)%val + params(i)%step*rand_gauss(rng)
        end do
        ! Accept?
        if(any((params(:npar)%next-params(:npar)%min)*(params(:npar)%next-params(:npar)%max) > 0)) then
           lik(2) = infinity
        else
           lik(2) = getlik(maps, params(:npar)%next)
        end if
        if(exp(lik(1)-lik(2)) > rand_uni(rng)) then
           lik(1) = lik(2)
           params(:npar)%val = params(:npar)%next
           nok = nok + 1
        end if
     end do
     ! Output sample
     do i = 1, npar
        write(*,'(e13.5)',advance="no") params(i)%val
     end do
     write(*,'(e17.9,f11.8)') lik(1), real(nok,dp)/nit
  end do

contains

  ! Model: map(1) = sum(map(i)*A(i)) + B
  function getlik(maps, params) result(lik)
    implicit none
    real(dp) :: maps(:,:,:), params(:), lik, dval, dvar
    lik = 0
    do i = 1, size(maps, 1)
       dval = 0; dvar = 0
       do j = 1, size(maps,2)
          dval = dval + maps(i,j,1)*params(j)
          dvar = dvar + maps(i,j,2)*params(j)**2
       end do
       dval = dval + params(size(maps,2)+1)
       lik  = lik + dval**2/dvar
    end do
    lik = 0.5*lik
  end function

  ! Reduce all maps to the lowest nside present, the interection of
  ! pixels and nested ordering
  subroutine make_compatible(imaps, iown, maps, pixels, nside, order)
    implicit none
    type(template_type) :: imaps(:)
    integer(i4b)        :: iown(:,:), nside, order
    integer(i4b)        :: i, j, k, l, n, npix, nmap, mtype(size(imaps))
    integer(i4b), allocatable :: hits(:,:), pixels(:)
    real(dp),     allocatable :: tmap(:,:), tmap2(:,:), maps(:,:,:)
    type(degrade_info)   :: deg
    nside = minval(imaps%nside)
    order = nest
    n     = size(imaps)
    npix  = 12*nside**2

    ! Classify each as a map or rms
    do i = 1, size(iown,2)
       do j = 1, 2
          mtype(iown(j,i)) = j
       end do
    end do

    ! Degrade and reorder
    allocate(hits(0:npix-1,0:n))
    hits  = 0
    do i = 1, n
       call prepare_degrade(deg, imaps(i)%pixels, imaps(i)%order, imaps(i)%nside, nside)
       allocate(tmap(size(imaps(i)%pixels),1), tmap2(size(deg%opix),1))
       tmap(:,1) = imaps(i)%map
       call degrade_map(tmap, deg, tmap2)
       if(mtype(i) == 2) tmap2 = tmap2 / (imaps(i)%nside/nside) ! downscale rms
       deallocate(imaps(i)%pixels, imaps(i)%map, tmap)
       allocate(imaps(i)%pixels(size(tmap2)), imaps(i)%map(size(tmap2)))
       imaps(i)%map = tmap2(:,1)
       do j = 1, size(tmap2)
          call set_ordering(nside, nest, imaps(i)%order, deg%opix(j), imaps(i)%pixels(j))
          k = imaps(i)%pixels(j)
          hits(k,0) = hits(k,0) + 1
          hits(k,i) = j
       end do
       call free_degrade_info(deg)
       deallocate(tmap2)
    end do
    ! At this point all are at the same nside and ordering, but may not share pixels
    ! Build the final pixel set
    call wherei(hits(:,0) == n, pixels)
    pixels = pixels-1
    m      = size(pixels)
    ! And fill in the map
    allocate(maps(m,size(iown,2),2))
    maps = 0
    ! Loop through the entries in the output
    do i = 1, size(iown,2)
       do j = 1, 2
          k = iown(j,i)
          if(k == 0) cycle
          ! And copy over the correct values
          do l = 1, m
             maps(l,i,j) = imaps(k)%map(hits(pixels(l),k))
          end do
       end do
    end do
    deallocate(hits)
  end subroutine

  subroutine help(msg)
    implicit none
    character(len=*), optional :: msg
    if(present(msg)) write(stderr,*) msg
    write(stderr,*) "template_fit template[:rmsmap] template[:rmsmap] [....]"
    stop
  end subroutine

end program
