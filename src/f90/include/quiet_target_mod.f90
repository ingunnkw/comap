module quiet_target_mod
  use quiet_defs
  use quiet_utils
  use quiet_ces_mod
  use quiet_acceptlist_mod
  implicit none

  type rlist
     real(dp), dimension(:,:), allocatable :: ranges
  end type

  type quiet_knife
     character(len=512) :: name, full_name
     integer(i4b)       :: type, nsplit
     real(dp)           :: zero, step
     type(rlist), dimension(:), allocatable :: lists
  end type

  type swiss_knife
    character(len=512) :: full_name
    type(quiet_knife), dimension(:), allocatable :: knives
  end type

  type diode_knife
     character(len=32), allocatable :: names(:)
     integer(i4b),      allocatable :: class(:,:), ns(:)
  end type

  type quiet_target
     character(len=512) :: object, split
     type(swiss_knife)  :: knife
     integer(i4b)       :: i, n ! this is # i of n subsets
     real(dp)           :: frac ! which fraction of this split this target is
     type(acceptlist)   :: alist
  end type

  type(diode_knife), private :: static_external, static_internal
  integer(i4b),      private, parameter :: round = 1, square = 2, curly = 3, sharp = 4

  ! Internally, jackknives will be handled by arrays of
  ! (ndi,nces), where 0 means not part of anything in this
  ! jackknife, 1 means part of part 1, and so on.
  !
  ! The input will be a corresponding array of (ndi,nces) of
  ! some property.

contains

  subroutine initialize_target_mod(parfile)
    implicit none
    character(len=*)   :: parfile
    character(len=512) :: static_file
    logical(lgt), save :: initialized = .false.
    logical(lgt) :: exist
    integer(i4b) :: unit
    if(initialized) return
    call initialize_module_mod(parfile)
    call get_parameter(0, parfile, "STATIC_JACKKNIFE_DEFS", par_string=static_file)
    unit = getlun()
    inquire(file=trim(static_file), exist=exist)
    if (exist .and. trim(static_file) /= '') call read_diode_knife(static_file, static_external)
    call initialize_internal_diode_knives(static_internal)
    initialized = .true.
  end subroutine

  subroutine init_target(target, alist)
    implicit none
    type(quiet_target) :: target
    type(acceptlist)   :: alist
    call free_target(target)
    call allocate_acceptlist(target%alist)
    target%alist = alist
  end subroutine

  subroutine copy_target(a, b)
    implicit none
    type(quiet_target) :: a, b
    call free_target(b)
    call allocate_acceptlist(b%alist)
    call copy_swiss(a%knife, b%knife)
    b = a
  end subroutine

  subroutine free_target(target)
    implicit none
    type(quiet_target) :: target
    call deallocate_acceptlist(target%alist)
    call free_swiss(target%knife)
  end subroutine

  ! Our formats:
  ! * name(nsplit)
  ! * name<edge,edge,edge..>
  ! * name[a:b+c,e:f,g], etc.
  ! * name{zero:step:nsplit}
  subroutine init_knife(knife, desc)
    implicit none
    type(quiet_knife) :: knife
    character(len=*)  :: desc
    character(len=1)  :: op, cl
    character(len=64) :: rtok(3)
    integer(i4b)      :: i, j, k, l, n, ntok, ntok2
    logical(lgt)      :: empty
    real(dp)          :: rdefault(2)
    character(len=32), dimension(:), allocatable :: tokens, tokens2
    call free_knife(knife)
    ! Find the type of parentheses we use. The last character
    ! in our string is a parenthesis, so just look at that.
    ! Nb: We might not have one, in which case we treat this as
    ! (2).
    n          = len_trim(desc)
    knife%type = round
    cl         = desc(n:n)
    empty      = .false.
    rdefault   = [ -infinity, infinity ]
    select case(cl)
       case(')'); op = '('; knife%type = round
       case(']'); op = '['; knife%type = square
       case('}'); op = '{'; knife%type = curly
       case('>'); op = '<'; knife%type = sharp
       case default; knife%type = round; empty = .true.
    end select
    if(empty) then
       knife%name      = desc
       knife%full_name = desc
       knife%nsplit    = 2
       allocate(knife%lists(0))
       return
    end if

    ! Now find the start of argument list
    j = index(desc, op)
    call assert(j > 1 .and. j < n, "Invalid jackknife format: " // trim(desc))
    knife%name      = desc(1:j-1)
    knife%full_name = desc

    ! And parse the horrible argument list
    select case(knife%type)
       case(round)
          read(desc(j+1:n-1),*) knife%nsplit
       case(curly)
          call get_tokens(desc(j+1:n-1),":",rtok, ntok)
          call assert(ntok == 3 .or. ntok == 2, "Invalid curly format: " // desc(j+1:n-1))
          read(rtok(1),*) knife%zero
          read(rtok(2),*) knife%step
          knife%nsplit = 0
          if(ntok >= 3) read(rtok(3),*) knife%nsplit
       case(square)
          ! This is the most complicated format. First parse the outer layer
          ntok = num_tokens(desc(j+1:n-1),",")
          allocate(tokens(ntok), knife%lists(ntok))
          call get_tokens(desc(j+1:n-1),",",tokens)
          do i = 1, size(tokens)
             ! Then get the set of ranges making up each side of the knife
             ntok = num_tokens(tokens(i),"+")
             allocate(tokens2(ntok), knife%lists(i)%ranges(2,ntok))
             call get_tokens(tokens(i),"+",tokens2)
             do k = 1, size(tokens2)
                ! Finally, parse the individual ranges. If no : is involved,
                ! it is assumed to be an empty range of the form a:a.
                l = index(tokens2(k),":")
                if(l == 0) then
                   read(tokens2(k),*) knife%lists(i)%ranges(1,k)
                   read(tokens2(k),*) knife%lists(i)%ranges(2,k)
                else
                   call get_tokens(tokens2(k),":",rtok, ntok, allow_empty=.true.)
                   do l = 1, ntok
                      if(len_trim(rtok(l)) == 0) then
                         knife%lists(i)%ranges(l,k) = rdefault(l)
                      else
                         read(rtok(l),*) knife%lists(i)%ranges(l,k)
                      end if
                   end do
                end if
             end do
             deallocate(tokens2)
          end do
          deallocate(tokens)
       case(sharp)
          ! This is a simpler version of square
          ntok = num_tokens(desc(j+1:n-1),",")
          allocate(tokens(ntok), knife%lists(ntok+1))
          do i = 1, ntok+1
             allocate(knife%lists(i)%ranges(2,1))
          end do
          call get_tokens(desc(j+1:n-1),",",tokens)
          do i = 1, ntok
             read(tokens(i),*) knife%lists(i  )%ranges(2,1)
             read(tokens(i),*) knife%lists(i+1)%ranges(1,1)
          end do
          knife%lists(1     )%ranges(1,1) = rdefault(1)
          knife%lists(ntok+1)%ranges(2,1) = rdefault(2)
          deallocate(tokens)
    end select

    ! Sharp is a special case of square, so the rest of the code does not
    ! need to know the difference
    if(knife%type == sharp) knife%type = square
    if(.not. allocated(knife%lists)) allocate(knife%lists(0))
  end subroutine

  subroutine free_knife(knife)
    implicit none
    type(quiet_knife) :: knife
    integer(i4b)      :: i
    if(allocated(knife%lists)) then
       do i = 1, size(knife%lists)
          if(allocated(knife%lists(i)%ranges)) deallocate(knife%lists(i)%ranges)
       end do
       deallocate(knife%lists)
    end if
  end subroutine

  subroutine copy_knife(a, b)
    implicit none
    type(quiet_knife), intent(in)    :: a
    type(quiet_knife), intent(inout) :: b
    integer(i4b)                     :: i
    call free_knife(b)
    allocate(b%lists(size(a%lists)))
    do i = 1, size(a%lists)
       allocate(b%lists(i)%ranges(size(a%lists(i)%ranges,1),size(a%lists(i)%ranges,2)))
    end do
    b = a
  end subroutine

  ! Set up a set of coupled knives (a Swiss army knife) by
  ! splitting on +
  subroutine init_swiss(swiss, desc)
    implicit none
    type(swiss_knife),  intent(inout) :: swiss
    character(len=*),   intent(in)    :: desc
    character(len=512), allocatable   :: toks(:)
    integer(i4b)                      :: i, ntok
    ntok = num_tokens(desc, "+", group="()[]{}<>")
    allocate(toks(ntok))
    call get_tokens(desc, "+", toks, group="()[]{}<>")
    call free_swiss(swiss)
    swiss%full_name = desc
    allocate(swiss%knives(ntok))
    do i = 1, ntok
       call init_knife(swiss%knives(i), toks(i))
    end do
    deallocate(toks)
  end subroutine

  subroutine free_swiss(swiss)
    implicit none
    type(swiss_knife) :: swiss
    integer(i4b)      :: i
    if(allocated(swiss%knives)) then
       do i = 1, size(swiss%knives)
          call free_knife(swiss%knives(i))
       end do
       deallocate(swiss%knives)
    end if
  end subroutine

  subroutine copy_swiss(a,b)
    implicit none
    type(swiss_knife), intent(in)    :: a
    type(swiss_knife), intent(inout) :: b
    integer(i4b)                     :: i
    call free_swiss(b)
    allocate(b%knives(size(a%knives)))
    do i = 1, size(a%knives)
       call copy_knife(a%knives(i), b%knives(i))
    end do
    b%full_name = a%full_name
  end subroutine

  subroutine filter_object(target, object)
    implicit none
    character(len=*)     :: object
    type(quiet_target)   :: target
    integer(i4b)         :: ntok, i, j
    type(quiet_ces_info) :: ces
    character(len=512), dimension(:), allocatable :: tokens
    ntok = num_tokens(object, ",+ ")
    allocate(tokens(ntok))
    call get_tokens(object, ",+ ", tokens)
    do i = 1, get_num_ces()
       call get_ces_info(i, ces)
       if(.not. (any(tokens == ces%object) .or. any(tokens == "all"))) &
        & call reject(target%alist, ces%cid, reason=REJECTED_ALIST)
    end do
    target%object = object
    deallocate(tokens)
  end subroutine

  subroutine jackknife(target, knifespecs, cstat, dstat,targets, knife_defs, knife_res)
    implicit none
    type(quiet_target) :: target
    character(len=*)   :: knifespecs
    real(sp)           :: cstat(:,:), dstat(:,:,:)
    integer(i4b)       :: i, j, k, l, m, n, ntok, ndi, nces, ndr, nmod
    type(quiet_target), dimension(:),     allocatable :: targets
    character(len=512), dimension(:),     allocatable :: tokens
    type(swiss_knife),  dimension(:),     allocatable :: knives
    integer(i4b),       dimension(:,:,:), allocatable :: groups
    integer(i4b),       dimension(:,:),   allocatable :: group, mask
    integer(i4b),       dimension(:),     allocatable :: lens
    type(swiss_knife),  optional,         allocatable :: knife_defs(:)
    integer(i4b),       optional,         allocatable :: knife_res(:,:,:)

    ! First get the list of individual knives
    ndr  = get_num_diodes()
    nmod = get_num_modules()
    ndi  = ndr*nmod
    nces = size(cstat, 1)
    ntok = num_tokens(knifespecs, ",", group="()[]{}<>")
    allocate(tokens(ntok), knives(ntok), groups(ndi,nces,ntok),lens(ntok),group(ndi, nces))
    allocate(mask(ndi, nces))
    call get_tokens(knifespecs, ",", group="()[]{}<>", toks=tokens)
    do i = 1, ntok
       call init_swiss(knives(i), tokens(i))
    end do
    deallocate(tokens)
    groups = 0
    mask   = 1-min(1,reshape(target%alist%status,[ndi,nces]))

    ! Then, for each knife, cut the original target into targets,
    ! and collect. Also transform into a format that is easier to work with
    ! than target lists. The first argument is 0 for rejected and 1 for accepted
    do i = 1, size(knives)
       do j = 1, size(knives(i)%knives)
          call jackknife_single(mask, knives(i)%knives(j), cstat,dstat, group)
          if(j == 1) then
             groups(:,:,i) = group
          else
             where(group == 0 .or. groups(:,:,i) == 0)
                groups(:,:,i) = 0
             elsewhere
                groups(:,:,i)=(groups(:,:,i)-1)*maxval(group)+group
             end where
          end if
       end do
       ! We may have ended up with some empty groups. For example,
       ! mjd<55000>+mjd[55500:] would have group 1 unpopulated.
       ! Could prune empty here, but that would make the knives unintuitive
       ! without some explanation. On the other hand, empty groups pose problems
       ! later. Not sure what to do.
       lens(i) = maxval(groups(:,:,i))
    end do
    deallocate(group, mask)
    n = sum(lens)

    ! And now reformat into the targets format (i.e. basically a fancy
    ! acceptlist).
    k = 0
    allocate(targets(n))
    do i = 1, size(lens)
       do j = 1, lens(i)
          k = k+1
          call copy_target(target, targets(k))
          call copy_swiss(knives(i), targets(k)%knife)
          where(reshape(groups(:,:,i),[ndr,nmod,nces]) /= j) &
            targets(k)%alist%status = REJECTED_SPLIT
          targets(k)%i    = j
          targets(k)%n    = lens(i)
          targets(k)%frac = real(count(groups(:,:,i) == j),dp)/count(groups(:,:,i) /= 0)
       end do
    end do

    ! This is a bit untidy, but output the raw knife information if asked for
    if(present(knife_defs)) then
       allocate(knife_defs(size(knives)))
       do i = 1, size(knives)
          call copy_swiss(knives(i), knife_defs(i))
       end do
    end if
    if(present(knife_res)) then
       allocate(knife_res(ndi, nces, size(knives)))
       knife_res = groups
    end if

    ! And clean up
    do i = 1, size(knives)
       call free_swiss(knives(i))
    end do
    deallocate(knives, lens, groups)
  end subroutine

  ! Ok, we here have to implement the jackknives requested.
  ! This is complicated by the many formats involved.
  ! Many of these can share code, though.
  ! Note: Not all knives support all types.
  subroutine jackknife_single(mask, knife, cstat, dstat, groups)
    implicit none
    integer(i4b),       intent(in) :: mask(:,:)
    type(quiet_knife),  intent(in) :: knife
    real(sp),           intent(in) :: cstat(:,:), dstat(:,:,:)
    integer(i4b),       intent(out):: groups(:,:)
    integer(i4b)       :: i, nsplit, nremap, status, ndi, mdi, nmod, n, nces, foo(size(groups,1))

    mdi    = get_num_diodes()
    nmod   = get_num_modules()
    ndi    = nmod*mdi
    nces   = size(groups,2)
    groups = 0

    ! If anything requires complicated treatment, put it here:
    ! if(split_foo(target, knife, stats, targets)) return

    ! The general ones, which support everything
    select case(knife%name)
       case("mjd");     call ces_split   (mask, cstat(:,STAT_MJD),                   knife, groups)
       case("az");      call ces_split   (mask, real(cstat(:,STAT_AZ) * RAD2DEG,sp), knife, groups)
       case("el");      call ces_split   (mask, real(cstat(:,STAT_EL) * RAD2DEG,sp), knife, groups)
       case("dk");      call ces_split   (mask, real(cstat(:,STAT_DK) * RAD2DEG,sp), knife, groups)
       case("lst");     call ces_split   (mask, cstat(:,STAT_LST),                   knife, groups)
       case("wind");    call ces_split   (mask, cstat(:,STAT_WIND),                  knife, groups)
       case("pwv");     call ces_split   (mask, cstat(:,STAT_PWV),                   knife, groups)
       case("tamb");    call ces_split   (mask, cstat(:,STAT_T_AMBIENT),             knife, groups)
       case("tenc");    call ces_split   (mask, cstat(:,STAT_TENC),                  knife, groups)
       case("dtenc");   call ces_split   (mask, cstat(:,STAT_TENC_CHANGE),           knife, groups)
       case("cryo");    call ces_split   (mask, cstat(:,STAT_CRYO),                  knife, groups)
       case("dcryo");   call ces_split   (mask, cstat(:,STAT_CRYO_CHANGE),           knife, groups)
       case("sigma0");  call cesmod_split(mask, dstat(:,:,DIODE_STAT_SIGMA0),        knife, groups)
       case("gain");    call cesmod_split(mask, dstat(:,:,DIODE_STAT_GAIN),          knife, groups)
       case("typeb");   call cesmod_split(mask, dstat(:,:,DIODE_STAT_TYPEB),         knife, groups)
       case("fknee");   call cesmod_split(mask, dstat(:,:,DIODE_STAT_FKNEE),         knife, groups)
       case("weather"); call cesmod_split(mask, dstat(:,:,DIODE_STAT_WEATHER1),      knife, groups)
       case("10hz");    call cesmod_split(mask, dstat(:,:,DIODE_STAT_10HZ),          knife, groups)
       case("sss");     call cesmod_split(mask, dstat(:,:,DIODE_STAT_SSS),           knife, groups)
       case("cid");     call ces_split   (mask, real(ces_db%ceses%cid,sp),           knife, groups)
       case("diode");   call diode_split (mask, real(irange(ndi),sp),                knife, groups)
       case("dtype");   call diode_split (mask, real(modulo(irange(ndi)-1,mdi)+1,sp),knife, groups)
       case("dleak");   call diode_split (mask, real(get_leak_dis(),sp),             knife, groups)
       case("leak");    call module_split(mask, real(get_leak_dis(),sp),             knife, groups)
       case("tleak");   call module_split(mask, real(get_ttleak_dis(), sp),          knife, groups)
       case default
          ! Special cases, which only support a few
          select case(knife%type)
          case(round)
             select case(knife%name)
             case("all"); groups = 1
             case("qu");  groups = spread(modulo(irange(ndi)/2,2)+1,2,nces)
             case default
                call static_split(knife%name, static_external, groups, status)
                if (status /= 0) call static_split(knife%name, static_internal, groups, status)
                call assert(status == 0,"Unrecognized knife " // trim(knife%name))
             end select
          case default
             call assert(.false., "Split type not supported for " // trim(knife%name))
          end select
    end select

    ! Enforce mask
    groups = groups*mask
  end subroutine

  ! Handle a ces-only split, with all types of inputs
  subroutine ces_split(mask, stats, knife, groups)
    implicit none
    integer(i4b),      intent(in) :: mask(:,:)
    real(sp),          intent(in) :: stats(:)
    type(quiet_knife), intent(in) :: knife
    integer(i4b),      intent(out):: groups(:,:)
    integer(i4b)                  :: gtmp(size(groups,2))
    call raw_split(real(sum(mask,1),dp), stats, knife, gtmp)
    groups = spread(gtmp,1,size(mask,1))
  end subroutine

  subroutine diode_split(mask, ids, knife, groups)
    implicit none
    integer(i4b),      intent(in) :: mask(:,:)
    real(sp),          intent(in) :: ids(:)
    type(quiet_knife), intent(in) :: knife
    integer(i4b),      intent(out):: groups(:,:)
    integer(i4b)                  :: gtmp(size(groups,1))
    call raw_split(real(sum(mask,2),dp), ids, knife, gtmp)
    groups = spread(gtmp,2,size(mask,2))
  end subroutine

  subroutine module_split(mask, stats, knife, groups)
    implicit none
    integer(i4b),      intent(in) :: mask(:,:)
    real(sp),          intent(in) :: stats(:)
    type(quiet_knife), intent(in) :: knife
    integer(i4b),      intent(out):: groups(:,:)
    integer(i4b)                  :: gtmp(size(quiet_horns)), mdi, nmod, nces
    mdi = get_num_diodes(); nmod = get_num_modules(); nces = size(mask,2)
    call raw_split(real(sum(sum(reshape(mask,[mdi,nmod,nces]),1),2),dp), &
     & sum(reshape(stats,[mdi,size(stats)/mdi]),1), knife, gtmp)
    groups = spread(reshape(spread(gtmp,1,mdi),[size(groups,1)]),2,size(groups,2))
  end subroutine

  subroutine cesdi_split(mask, stats, knife, groups)
    implicit none
    integer(i4b),      intent(in) :: mask(:,:)
    real(sp),          intent(in) :: stats(:,:)
    type(quiet_knife), intent(in) :: knife
    integer(i4b),      intent(out):: groups(:,:)
    integer(i4b)                  :: gtmp(size(groups))
    call raw_split(real(reshape(mask,[size(mask)]),dp), &
     & reshape(stats,[size(stats)]), knife, gtmp)
    groups = reshape(gtmp,[size(groups,1),size(groups,2)])
  end subroutine

  ! This one is needed because the ideal cesdi split would
  ! break up too many assemblies, making everything very slow.
  ! This may be a tad unreadable :/
  subroutine cesmod_split(mask, stats, knife, groups)
    implicit none
    integer(i4b),      intent(in) :: mask(:,:)
    real(sp),          intent(in) :: stats(:,:)
    type(quiet_knife), intent(in) :: knife
    integer(i4b),      intent(out):: groups(:,:)
    integer(i4b),      allocatable:: gtmp(:)
    integer(i4b)                  :: mdi, nces
    mdi = get_num_diodes()
    allocate(gtmp(size(groups)/mdi))
    call raw_split(real(sum(reshape(mask,[mdi,size(mask)/mdi]),1),dp), &
     & sum(reshape(stats,[mdi,size(stats)/mdi]),1)/mdi, knife, gtmp)
    !call raw_split(real(sum(reshape(mask,[mdi,size(mask)/mdi]),1),dp), &
    ! & reshape(stats,[size(stats)]), knife, gtmp)
    groups = reshape(spread(gtmp,1,mdi),[size(groups,1),size(groups,2)])
    deallocate(gtmp)
  end subroutine

  ! The engine for performing splits. Does things in one dimension. Callers must
  ! collapse and expand before calling this.
  subroutine raw_split(weights, stats, knife, groups)
    implicit none
    real(dp),          intent(in) :: weights(:)
    real(sp),          intent(in) :: stats(:)
    type(quiet_knife), intent(in) :: knife
    integer(i4b),      intent(out):: groups(:)
    real(dp)                      :: cum(size(weights))
    integer(i4b)                  :: inds(size(weights)), i, j, k, n, first
    n = size(weights)
    groups = 0
    select case(knife%type)
    case(round)
       ! Build up a cumulative list based on the weights
       inds = irange(n)
       call get_sort_inds(real(stats,dp), inds)
       cum   = weights
       do first = 1, size(cum); if(cum(inds(first)) > 0) exit; end do
       if(first > size(cum)) return ! What? Empty weights?
       do i = 2, n
          cum(inds(i)) = cum(inds(i)) + cum(inds(i-1))
       end do
       cum  = cum/cum(inds(n))
       ! Allocate approximately equal amounts of data to each
       groups = floor((cum-cum(inds(first)))*knife%nsplit)+1
    case(curly)
       if(knife%nsplit > 0) then
          inds = floor(modulo((stats-knife%zero)/knife%step,real(knife%nsplit,dp)))
       else
          inds = floor((stats-knife%zero)/knife%step)
       end if
       groups = inds-minval(inds)+1
    case(square)
       do i = 1, size(knife%lists)
          do j = 1, size(knife%lists(i)%ranges,2)
             do k = 1, n
                if(moresim(real(stats(k),dp), knife%lists(i)%ranges(1,j)) .and. &
                 & lesssim(real(stats(k),dp), knife%lists(i)%ranges(2,j))) groups(k) = i
             end do
          end do
       end do
    end select
  end subroutine

  subroutine static_split(name, knives, groups, status)
    implicit none
    character(len=*),  intent(in) :: name
    type(diode_knife), intent(in) :: knives
    integer(i4b),      intent(out):: groups(:,:), status
    integer(i4b)                  :: i, j, k, m, n
    status = 1
    do i = 1, size(knives%names)
       if(knives%names(i) == name) exit
    end do
    if(i > size(knives%names)) return
    groups = spread(knives%class(:,i),2, size(groups,2))
    status = 0
  end subroutine

  subroutine get_sort_inds(vals, inds)
    implicit none
    real(dp)     :: vals(:)
    integer(i4b) :: inds(:), i
    real(dp), dimension(:), allocatable :: copy
    allocate(copy(size(vals)))
    copy = vals
    do i = 1, size(inds)
       inds(i) = i
    end do
    call quicksort(inds, copy)
    deallocate(copy)
  end subroutine

  subroutine read_diode_knife(fname, dknife)
    implicit none
    character(len=*)   :: fname
    character(len=512) :: line
    type(diode_knife)  :: dknife
    integer(i4b)       :: unit, i, n, ndi, nknife
    ndi  = get_num_diodes()*get_num_modules()
    unit = getlun()
    call free_diode_knife(dknife)
    open(unit,file=fname,action="read",status="old")
    read(unit,'(a)') line
    nknife = num_tokens(line, " 	")
    allocate(dknife%names(nknife),dknife%class(ndi,nknife),dknife%ns(nknife))
    call get_tokens(line, " 	", dknife%names)
    do i = 1, ndi
       read(unit,*) dknife%class(i,:)
    end do
    close(unit)
    do i = 1, nknife
       dknife%ns(i) = maxval(dknife%class(:,i))
    end do
  end subroutine

  subroutine free_diode_knife(dknife)
    implicit none
    type(diode_knife) :: dknife
    if(allocated(dknife%names)) deallocate(dknife%names)
    if(allocated(dknife%class)) deallocate(dknife%class)
    if(allocated(dknife%ns))    deallocate(dknife%ns)
  end subroutine

  subroutine interleaved(target, nsplit, targets)
    implicit none

    type(quiet_target) :: target
    integer(i4b)       :: nsplit
    type(quiet_target), dimension(:), allocatable :: targets

  end subroutine interleaved

  subroutine initialize_internal_diode_knives(dknife)
    implicit none

    type(diode_knife), intent(inout) :: dknife

    integer(i4b) :: i, j, current, n, nmod, ndi, mab, pos
    real(dp)     :: phi
    real(dp),     allocatable, dimension(:) :: vals
    integer(i4b), allocatable, dimension(:) :: inds
    integer(i4b),              dimension(6) :: good_mab = [0,2,4,5,9,10]

    current = 0
    n       = 5
    nmod    = get_num_modules()
    ndi     = get_num_diodes()
    allocate(dknife%names(n), dknife%class(nmod*ndi,n), dknife%ns(n))

    ! QU split
    current = current+1
    call assert(current <= n, 'quiet_target_mod: Somebody forgot to update n')
    dknife%names(current) = 'qu'
    dknife%ns(current)    = 2
    do i = 1, nmod*ndi
       dknife%class(i,current) = mod(i/2,2)+1
    end do

    ! Alternating modules
    current = current+1
    call assert(current <= n, 'quiet_target_mod: Somebody forgot to update n')
    dknife%names(current) = 'modalt'
    dknife%ns(current)    = 2
    do i = 1, nmod*ndi
       dknife%class(i,current) = mod((i-1)/ndi,2)+1
    end do

    ! Module frequency
    current = current+1
    call assert(current <= n, 'quiet_target_mod: Somebody forgot to update n')
    dknife%names(current) = 'modfreq'
    dknife%ns(current)    = 2
    allocate(vals(nmod),inds(nmod))
    do i = 1, nmod
       vals(i) = get_module_freq(i-1)
    end do
    call get_sort_inds(vals,inds)
    do i = 1, ndi
       dknife%class(ndi*(inds(1:nmod/2)-1)+i,current)      = 1
       dknife%class(ndi*(inds(nmod/2+1:nmod)-1)+i,current) = 2
    end do
    deallocate(vals, inds)

    ! Good vs bad MABs
    current = current+1
    call assert(current <= n, 'quiet_target_mod: Somebody forgot to update n')
    dknife%names(current) = 'mab'
    dknife%ns(current)    = 2
    do i = 0, nmod-1
       call get_modmab(i, mab, pos)
       do j = 1, ndi
          if (any(mab == good_mab)) then
             dknife%class(ndi*i+j,current) = 1
          else
             dknife%class(ndi*i+j,current) = 2
          end if
       end do
    end do

    ! Module position
    current = current+1
    call assert(current <= n, 'quiet_target_mod: Somebody forgot to update n')
    dknife%names(current) = 'inout'
    dknife%ns(current)    = 2
    allocate(vals(nmod),inds(nmod))
    do i = 0, nmod-1
       call get_module_pos(i, vals(i+1), phi)
    end do
    call get_sort_inds(vals,inds)
    do i = 1, ndi
       dknife%class(ndi*(inds(1:nmod/2)-1)+i,current)      = 1
       dknife%class(ndi*(inds(nmod/2+1:nmod)-1)+i,current) = 2
    end do
    deallocate(vals, inds)

  end subroutine initialize_internal_diode_knives

  ! This was amazingly complicated for such a simple task!
  function get_leak_dis() result(res)
    implicit none
    real(dp)          :: res(size(quiet_diodes)), stokes(3)
    integer(i4b)      :: i, j(1)
    do i = 1, size(quiet_diodes)
       res(i) = 0
       stokes = abs(quiet_diodes(i)%stokes)
       j      = maxloc(stokes)
       if(stokes(j(1)) == 0) cycle
       stokes = stokes/stokes(j(1))
       stokes(j(1)) = 0
       res(i) = maxval(stokes)
    end do
  end function

  ! Finding TT leak for all diodes
  function get_ttleak_dis() result(res)
    implicit none
    real(dp)         :: res(size(quiet_diodes))
    integer(i4b)     :: i, ndi

    ndi = get_num_diodes() 
    do i=0, get_num_modules()-1
       res(1+i*ndi:ndi+i*ndi) = quiet_horns(i)%amps(size(quiet_horns(i)%amps))
    end do
  end function

  function get_leak_mods() result(res)
    implicit none
    real(dp)      :: res(size(quiet_horns)), tmp(size(quiet_diodes))
    integer(i4b)  :: i
    tmp = get_leak_dis()
    do i = 1, size(quiet_horns)
       res(i) = maxval(tmp(quiet_horns(i-1)%diodes))
    end do
  end function

end module
