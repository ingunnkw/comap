! ***********************************************************************
!
!      Module for handling files in the Quiet software package.
!
!     Written by Sigurd K. Næss and Hans Kristian Eriksen and IKW, 2009
!
! ***********************************************************************


module quiet_fileutils
  use healpix_types
  use quiet_utils
  use fitstools
  use head_fits
  use l1_read_mod
  use quiet_hdf_mod
  use quiet_mapfile_mod
  use quiet_covfile_mod
  implicit none

  type l1_mapping
     integer(i4b) :: n
     real(dp),     dimension(:), allocatable :: starts
     integer(i4b), dimension(:), allocatable :: lens
  end type l1_mapping

  interface read_table
     module procedure read_table_dp
  end interface

  type module_struct
     integer(i4b)                                :: module_number, nside_point
     real(dp),     allocatable, dimension(:,:)   :: pointing, tod, orig_point
     real(dp),     allocatable, dimension(:)     :: time
     integer(i4b)                                :: reduction
     integer(i4b), allocatable, dimension(:)     :: mask2map, nhits
     real(dp),     allocatable, dimension(:,:)   :: tp, tp_var
  end type module_struct

  interface deallocate_module_struct
     module procedure deallocate_module_struct_single, deallocate_module_struct_array
  end interface

contains

  subroutine read_beam(beamfile, beam)
    implicit none

    character(len=128),                   intent(in)  :: beamfile
    real(dp),           dimension(0:,1:), intent(out) :: beam

    integer(i4b) :: l, lmax, nmaps
    real(dp)     :: sigma_sq
    real(dp),          allocatable, dimension(:,:)   :: inbeam
    character(len=80),              dimension(1:180) :: header

    lmax  = size(beam(:,1))-1
    nmaps = size(beam(0,:))

    ! Seem to remember there is something weird going on with the WMAP beams when reading only 
    ! one component at a time. Remove this wrapper if you feel more comfortable with that...
    allocate(inbeam(0:lmax,4))

    call fits2cl(beamfile, inbeam, lmax, 4, header)

    beam(0:lmax,1:nmaps) = inbeam(0:lmax,1:nmaps)

    if (nmaps > 1) then
       if (sum(beam(:,2)) < 1.d0) beam(:,2) = beam(:,1)
       if (sum(beam(:,3)) < 1.d0) beam(:,3) = beam(:,2)
    end if

    deallocate(inbeam)

  end subroutine read_beam

  subroutine read_powspec(filename, cls)
    implicit none

    character(len=128),                   intent(in)  :: filename
    real(dp),           dimension(0:,1:), intent(out) :: cls

    integer(i4b) :: l, lmax, nmaps
    real(dp)     :: sigma_sq
    real(dp),          allocatable, dimension(:,:)   :: cls_in
    character(len=80),              dimension(1:180) :: header
    
    lmax  = size(cls(:,1))-1
    nmaps = size(cls(0,:))

    call fits2cl(filename, cls, lmax, nmaps, header)

  end subroutine read_powspec

  ! Read power spectrum from ascii file. Allocates ps to have
  ! the correct size. Rescales to physical scale by default.
  subroutine read_powspec_ascii(file, ps, scale)
    implicit none
    character(len=*)   :: file
    character(len=512) :: line
    integer(i4b)       :: i,j,n,lmax,l,unit,ncol
    real(dp)           :: v
    logical(lgt)       :: s
    logical(lgt),                         optional :: scale
    real(dp),          dimension(:,:), allocatable :: ps
    real(dp),          dimension(:),   allocatable :: nums
    s = .true.; if(present(scale)) s = scale
    unit = getlun()
    open(unit,file=file,action="read",status="old")
    lmax = -1
    ncol =  0
    do
       read(unit,'(a)',end=1) line
       if(line(1:1) == "#") cycle
       n    = num_tokens(line, " 	")
       read(line,*) l
       lmax = max(l,lmax)
       if(ncol == 0) ncol = n
       call assert(ncol == n, "Non-constant number of columns in " // trim(file))
    end do
1   rewind(unit)
    allocate(nums(ncol-1),ps(0:lmax,ncol-1))
    ps = 0
    do
       read(unit,'(a)',end=2) line
       if(line(1:1) == "#") cycle
       read(line,*) l, nums
       ps(l,:) = nums
    end do
2   close(unit)
    deallocate(nums)
    if(s) then
       do l = 1, lmax
          ps(l,:) = ps(l,:) * 2*pi/(l*(l+1))
       end do
    end if
  end subroutine




  subroutine read_pixwin(nside, nmaps, pixwin, filename)
    implicit none

    integer(i4b),                              intent(in)  :: nside, nmaps
    real(dp),         pointer, dimension(:,:)              :: pixwin
    character(len=*),                          intent(in), optional :: filename

    integer(i4b)        :: nc
    character(len=128)  :: pixwin_file
    character(len=4)    :: nside_text
    logical(lgt)        :: exist, anynull, binary
    real(dp)            :: nullval
    character(len=80),              dimension(1:180) :: header

    allocate(pixwin(0:4*nside,nmaps))
    
    if (nmaps == 3) then
       nc = 2
    else
       nc = 1
    end if

    if (present(filename)) then
       pixwin_file = filename
       binary      = .false.
    else
       call int2string(nside, nside_text)
       pixwin_file = 'pixel_window_n' // nside_text // '.fits'
       binary = .true.
    end if
    inquire(file=pixwin_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       if (binary) then
          call read_dbintab(pixwin_file, pixwin(0:4*nside,1:nc), 4*nside+1, nc, nullval, anynull)
       else
          call fits2cl(pixwin_file, pixwin(0:4*nside,1:nc), 4*nside, nc, header)          
       end if
       if (nmaps > 1) then
          if (sum(abs(pixwin(:,2))) == 0.d0) pixwin(:,2) = pixwin(:,1)
          if (nmaps == 3) pixwin(:,3) = pixwin(:,2)
       end if
    else
       pixwin = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Pixel window file ', trim(pixwin_file), ' not found. '
       write(*,*) 'Using unity weights.'
       write(*,*) ''
    end if

  end subroutine read_pixwin

  subroutine read_equ_set(unit, filename, n, ordering, polarization, invcov, rhs, npix, map2mask)
    implicit none
  
    character(len=256),                     intent(in) :: filename
    integer(i8b),                           intent(out):: npix, n
    integer(i4b),                           intent(in) :: unit
    integer(i4b),                           intent(out):: ordering, polarization

    real(dp), allocatable, dimension(:,:), intent(out) :: invcov, rhs
    integer(i4b), allocatable,dimension(:),intent(out) :: map2mask

    integer(i4b)                                       :: int_in, i, nrhs

    write(*,*) 'Reading from ', trim(filename)
    write(*,*) "Warning: read_equ_set uses an obsolete equation set format"

    open(unit, file=trim(filename), form='unformatted')
    read(unit) int_in
    read(unit) ordering
    read(unit) polarization
    n = int_in
    allocate(invcov(n,n))  
    do i = 1, n
       read(unit) invcov(:,i)
    end do
    read(unit) nrhs
    allocate(rhs(n,nrhs))
    do i = 1, nrhs
       read(unit) rhs(:,i)
    end do
    read(unit) int_in
    npix = int_in
    allocate(map2mask(0:npix-1))            
    read(unit) map2mask
    close(unit)
  end subroutine read_equ_set

  ! Hm, we have a whole forest of these routines now, differing in
  ! whether they are verbose or not, and in the number of dimensions
  ! out the matrix used. Oh well. Somebody should probably tidy this up
  ! at some point.

  subroutine write_covmatrix(unit, filename, ordering, polarization, matrix, inv)

    integer(i4b),                       intent(in) :: unit, ordering, polarization
    character(len=*),                   intent(in) :: filename
    logical(lgt),                       intent(in) :: inv
    real(dp),           dimension(:,:), intent(in) :: matrix

    integer(i4b) :: n, m, i

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    m = size(matrix(:,1))
    n = size(matrix(1,:))
    if (n/= m) then
       write(*,*) ' cov matrix not n x n matrix. quiting'
       stop
    end if

    open(unit, file=trim(filename), form='unformatted')
    write(unit) m
    write(unit) ordering
    write(unit) polarization
    do i = 1, n
       write(unit) matrix(:,i)
    end do
    write(unit) inv
    close(unit)

  end subroutine write_covmatrix

  subroutine read_covmatrix(unit, filename, ordering, polarization, matrix, inv, n)

    integer(i4b),                          intent(in)  :: unit
    integer(i8b),                          intent(out) :: n
    integer(i4b),                          intent(out) :: ordering, polarization
    character(len=*),                      intent(in)  :: filename
    logical(lgt),                          intent(out) :: inv
    real(dp), allocatable, dimension(:,:), intent(out) :: matrix

    integer(i4b) :: n_in, i

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    open(unit, file=trim(filename), form='unformatted')
    read(unit) n_in
    n = n_in
    write(*,*) n, '= n'
    read(unit) ordering
    write(*,*) ordering, '= ordering'
    read(unit) polarization
    write(*,*) polarization, '= polarisation'
    allocate(matrix(n,n))
    do i = 1, n
       read(unit) matrix(:,i)
    end do
    read(unit) inv
    close(unit)

  end subroutine read_covmatrix

  subroutine write_cov_mat(unit, filename, ordering, polarization, inv_cov, covmat, append)
    implicit none

    integer(i4b),                       intent(in) :: unit, ordering, polarization
    character(len=*),                   intent(in) :: filename
    logical(lgt),                       intent(in) :: append, inv_cov
    real(dp),         dimension(1:,1:), intent(in) :: covmat

    integer(i4b) :: i, n, m

    m = size(covmat(:,1))
    n = size(covmat(1,:))

    ! polarization = 1  => I only
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU

    if (append) then
       open(unit, file=trim(filename), form='unformatted', position='append')
       do i = 1, n
          write(unit) covmat(:,i)
       end do
       close(unit)
    else
       open(unit, file=trim(filename), form='unformatted')
       write(unit) m
       write(unit) ordering
       write(unit) polarization
       do i = 1, n
          write(unit) covmat(:,i)
       end do
       write(unit) inv_cov
       close(unit)
    end if

  end subroutine write_cov_mat

  ! Given a list of files and a mjd range, reads out and concatenates level1-data.
  subroutine l1_read_range(filelist, mjd, mods, sel, data, point, hk, status)
    character(len=*) :: filelist(:)
    real(dp)         :: mjd(2)
    integer(i4b)     :: mods(:), eff(2)
    type(l1_selector),  optional :: sel
    type(data_struct),  optional :: data
    type(point_struct), optional :: point
    type(hk_struct),    optional :: hk
    integer(i4b),       optional :: status

    integer(i4b)        :: i, j, k, m, n, nsamp, unit, lr(2), gr(2), samprate, s
    integer(i4b)        :: nbias, ncryo, nencl, nperi
    real(dp)            :: start_mjd
    type(l1_selector)   :: sel_
    type(data_struct)   :: data_
    type(point_struct)  :: point_
    type(hk_struct)     :: hks(size(filelist))

    if(present(sel)) then
       sel_ = sel
    else
       sel_ = l1sel(l1_std)
    end if
    if(.not. present(data))  sel_%on(data_all)  = .false.
    if(.not. present(point)) sel_%on(point_all) = .false.
    if(.not. present(hk))    sel_%on(hk_all)    = .false.
    sel_%on(point_time) = .true.

    samprate = 100
    s        = 24*60*60*samprate
    nsamp    = int((mjd(2)-mjd(1))*s)
    if(present(data))  call allocate_data_struct(size(mods), nsamp, data, sel_)
    if(present(point)) call allocate_point_struct(nsamp, point, sel_)

    ! Hack: Check for uninitialized stuff
    if(present(point)) point%encoder_elevation = NaN
    if(present(point)) point%encoder_azimuth   = NaN
    if(present(point)) point%encoder_deck      = NaN

    ! Housekeeping does not have the same samprate, but take little space, so
    ! we just read it all in, and figure out the ranges afterwards.

    gr = 0
    unit = getlun()
    do i = 1, size(filelist)
       call l1_read(unit, filelist(i), data_, hks(i), point_, modules=mods, selector=sel_, &
        & status_code=status)
       if(present(status)) then; if(status /= 0) then
          call deallocate_data_struct(data)
          call deallocate_point_struct(point)
          do j = 1, i-1
             call deallocate_hk_struct(hks(j))
          end do
          return
       end if; end if

       ! Argh! A few files have non-monotonous time ranges! Skip these ranges
       if(i == 1) then
          do j = size(point_%time), 2, -1
             if(point_%time(j) < point_%time(j-1)) exit
          end do
          eff = [ j, size(point_%time) ]
       else
          do j = 1, size(point_%time)-1
             if(point_%time(j+1) < point_%time(j)) exit
          end do
          eff = [ 1, j ]
       end if

       ! Map the times in the file into indices into range. Because of the possibility
       ! varying sample rates etc. we loop through.
       do j = eff(1), eff(2)
          if(point_%time(j) >= mjd(1)) exit
       end do
       lr(1) = j
       do j = lr(1), eff(2)
          if(point_%time(j) > mjd(2)) exit
       end do
       lr(2) = j-1
       n = lr(2)-lr(1)+1
       ! Do we have room for all this?
       if(gr(2) + n >  nsamp) lr(2) = lr(2) - (gr(2)+n-nsamp)
       ! Append this to the global range
       gr = [ gr(2)+1, gr(2)+1+(lr(2)-lr(1)) ]

       if(i == 1) start_mjd = point_%time(lr(1))


       ! And finally copy over the interesting parts
       if(present(point)) then
          if(allocated(point%time)) point%time(gr(1):gr(2)) = point_%time(lr(1):lr(2))
          if(allocated(point%mode)) point%mode(gr(1):gr(2)) = point_%mode(lr(1):lr(2))

          if(allocated(point%encoder_azimuth)) &
           & point%encoder_azimuth  (gr(1):gr(2)) = point_%encoder_azimuth(lr(1):lr(2))
          if(allocated(point%encoder_elevation)) &
           & point%encoder_elevation(gr(1):gr(2)) = point_%encoder_elevation(lr(1):lr(2))
          if(allocated(point%encoder_deck)) &
           & point%encoder_deck(gr(1):gr(2))      = point_%encoder_deck(lr(1):lr(2))

          if(allocated(point%command_azimuth)) &
           & point%command_azimuth  (gr(1):gr(2)) = point_%command_azimuth(lr(1):lr(2))
          if(allocated(point%command_elevation)) &
           & point%command_elevation(gr(1):gr(2)) = point_%command_elevation(lr(1):lr(2))
          if(allocated(point%command_deck)) &
           & point%command_deck(gr(1):gr(2))      = point_%command_deck(lr(1):lr(2))
       end if

       if(present(data)) then
          if(allocated(data%time))  data%time(gr(1):gr(2))  = data_%time(lr(1):lr(2))
          if(allocated(data%phase)) data%phase(gr(1):gr(2)) = data_%phase(lr(1):lr(2))
          do j = 1, size(data%RQ)
             if(allocated(data%RQ(j)%demod)) &
              & data%RQ(j)%demod(:,gr(1):gr(2))        = data_%RQ(j)%demod(:,lr(1):lr(2))
             if(allocated(data%RQ(j)%avg)) &
              & data%RQ(j)%avg(:,gr(1):gr(2))          = data_%RQ(j)%avg(:,lr(1):lr(2))
             if(allocated(data%RQ(j)%quad)) &
              & data%RQ(j)%quad(:,gr(1):gr(2))         = data_%RQ(j)%quad(:,lr(1):lr(2))

             if(allocated(data%RQ(j)%scale)) &
              & data%RQ(j)%scale(gr(1):gr(2))          = data_%RQ(j)%scale(lr(1):lr(2))
             if(allocated(data%RQ(j)%quad_scale)) &
              & data%RQ(j)%quad_scale(gr(1):gr(2))     = data_%RQ(j)%quad_scale(lr(1):lr(2))

             if(allocated(data%RQ(j)%demod_raw)) &
              & data%RQ(j)%demod_raw(:,gr(1):gr(2))    = data_%RQ(j)%demod_raw(:,lr(1):lr(2))
             if(allocated(data%RQ(j)%avg_raw)) &
              & data%RQ(j)%avg_raw(:,gr(1):gr(2))      = data_%RQ(j)%avg_raw(:,lr(1):lr(2))
             if(allocated(data%RQ(j)%quad_raw)) &
              & data%RQ(j)%quad_raw(:,gr(1):gr(2))     = data_%RQ(j)%quad_raw(:,lr(1):lr(2))
          end do
       end if
       call deallocate_data_struct(data_)
       call deallocate_point_struct(point_)
    end do
    if(present(point)) point%start_mjd = start_mjd
    if(present(data))  data%start_mjd  = start_mjd

    ! Finally deal with the housekeeping. Notice all the duplicate code here.
    ! This is here because bias, cryo and encl are named members, not array elements.
    if(present(hk)) then
       nbias = 0; ncryo = 0; nencl = 0; nperi = 0
       do i = 1, size(hks)
          if(.not. allocated(hks(i)%bias%time)) cycle ! File did not contain hk
          nbias = nbias + count(hks(i)%bias%time >= mjd(1) .and. hks(i)%bias%time < mjd(2))
          ncryo = ncryo + count(hks(i)%cryo%time >= mjd(1) .and. hks(i)%cryo%time < mjd(2))
          nencl = nencl + count(hks(i)%encl%time >= mjd(1) .and. hks(i)%encl%time < mjd(2))
          nperi = nperi + count(hks(i)%peri%time >= mjd(1) .and. hks(i)%peri%time < mjd(2))
       end do
       call allocate_hk_struct(nbias, ncryo, nencl, nperi, hk)
       nbias = 0; ncryo = 0; nencl = 0; nperi = 0
       do i = 1, size(hks)
          if(.not. allocated(hks(i)%bias%time)) cycle ! File did not contain hk
          do j = 1, hks(i)%bias%n_t
             if(hks(i)%bias%time(j) < mjd(1) .or. hks(i)%bias%time(j) >= mjd(2)) cycle
             nbias = nbias + 1
             hk%bias%time(nbias)    = hks(i)%bias%time(j)
             hk%bias%value(nbias,:) = hks(i)%bias%value(j,:)
          end do
          do j = 1, hks(i)%cryo%n_t
             if(hks(i)%cryo%time(j) < mjd(1) .or. hks(i)%cryo%time(j) >= mjd(2)) cycle
             ncryo = ncryo + 1
             hk%cryo%time(ncryo)    = hks(i)%cryo%time(j)
             hk%cryo%value(ncryo,:) = hks(i)%cryo%value(j,:)
          end do
          do j = 1, hks(i)%encl%n_t
             if(hks(i)%encl%time(j) < mjd(1) .or. hks(i)%encl%time(j) >= mjd(2)) cycle
             nencl = nencl + 1
             hk%encl%time(nencl)    = hks(i)%encl%time(j)
             hk%encl%value(nencl,:) = hks(i)%encl%value(j,:)
          end do
          do j = 1, hks(i)%peri%n_t
             if(hks(i)%peri%time(j) < mjd(1) .or. hks(i)%peri%time(j) >= mjd(2)) cycle
             nperi = nperi + 1
             hk%peri%time(nperi)    = hks(i)%peri%time(j)
             hk%peri%value(nperi,:) = hks(i)%peri%value(j,:)
          end do
       end do
       do i = 1, size(hks)
          call deallocate_hk_struct(hks(i))
       end do
    end if
  end subroutine

  subroutine read_ringweights(nside, polarization, weights)
    implicit none
    
    integer(i4b),                          intent(in)  :: nside
    logical(lgt),                          intent(in)  :: polarization
    real(dp),     pointer, dimension(:,:)              :: weights
    
    character(len=128)  :: weight_file
    character(len=5)    :: nside_text
    integer(i4b)        :: pol
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    if (polarization) then
       allocate(weights(1:2*nside,3))       
       pol = 3
    else 
       allocate(weights(1:2*nside,1))
       pol = 1
    end if

    call int2string(nside, nside_text)
    weight_file = 'weight_ring_n' // nside_text // '.fits'
    inquire(file=weight_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(weight_file, weights, 2*nside, pol, nullval, anynull)
       weights = 1.d0 + weights
    else
       weights = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Weight file ', trim(weight_file), ' not found. '
       write(*,*) 'Using unity weights in the spherical harmonic transforms.'
       write(*,*) ''
    end if
    
  end subroutine read_ringweights

  subroutine mkdir(path)
    implicit none
    character(len=*) :: path
    logical(lgt) :: exist
    integer(i4b) :: i, left
    exist = is_dir(path)
    if(exist) return
    if(mkdir_c(trim(path))) return
    do i = 1, 20
       call fsleep(1d0)
       if(is_dir(path)) return
       if(mkdir_c(trim(path))) return
    end do
    write(*,*) "Failed to create directory '" // trim(path) // "' after " // trim(itoa(i-1)) &
         & // " attempts! You're on your own now."
  end subroutine

  subroutine mkdirs(path, skiplast)
    implicit none
    character(len=*) :: path
    logical(lgt) :: skiplast
    integer(i4b) :: i,j,k,n,m
    n = len_trim(path)
    if(skiplast) then
       do while(path(n:n) /= '/')
          n = n-1
          if(n <= 1) return
       end do
       n = n-1
    end if

    j = 1
    do while(.true.)
       do while(path(j:j) == '/'); j = j+1; if(j > n) return; end do
       i = j
       do while(path(j:j) /= '/'); j = j+1; if(j > n) exit; end do
       call mkdir(path(1:j-1))
       if(j > n) return
    end do
  end subroutine

  subroutine mv(from, to)
    character(len=*) :: from, to
    logical(lgt) :: error
    error = .not. mv_c(trim(from), trim(to))
    if(error) write(*,*) "Failed to move '" // trim(from) // "' to '" // trim(to) // "'!"
  end subroutine

  subroutine rm(file, noerr)
    character(len=*) :: file
    logical(lgt) :: error, check
    logical(lgt), optional :: noerr
    check = .true.; if(present(noerr)) check = .not. noerr
    error = .not. rm_c(trim(file))
    if(check .and. error) write(*,*) "Failed to remove '" // trim(file) // "'!"
  end subroutine

  subroutine touch(file)
    character(len=*) :: file
    integer(i4b) :: unit
    unit = getlun()
    open(unit,file=file)
    close(unit)
  end subroutine

  function file_exists(filename) result(res)
    implicit none
    character(len=*) :: filename
    logical(lgt)     :: res
    inquire(file=filename, exist=res)
  end function

  ! Read the given modules from file into data. If modules isn't present,
  ! all modules will be read. If you specify modules 4, 8 and 9, for example,
  ! then data(1) will correspond to module 4, data(2) to module 8 and 3 to 9.
  ! It is an error to ask for a module that isn't present.
  subroutine L2_read(filename, modules, data, status, mask2map_only)
    implicit none

    character(len=*),                                  intent(in)  :: filename
    type(module_struct), allocatable, dimension(:),    intent(out) :: data
    integer(i4b),                            optional, intent(out) :: status
    integer(i4b),              dimension(:), optional, intent(in)  :: modules
    logical(lgt),                            optional, intent(in)  :: mask2map_only

    integer(i4b) :: i, j, nmod, nget, nsamp, mod, red, nred, unit, nside_point
    integer(i4b), allocatable, dimension(:)   :: mods, inds, npix
    real(dp),     allocatable, dimension(:)   :: time
    real(sp),     allocatable, dimension(:,:) :: orig_point, buffer

    unit = getlun()
    call deallocate_module_struct(data)

    ! Read L2 information
    open(unit, file=trim(filename), form='unformatted', action="read", err=1)
    read(unit, err=1, end=1) nmod
    read(unit, err=1, end=1) nsamp
    read(unit, err=1, end=1) red
    read(unit, err=1, end=1) nred
    read(unit, err=1, end=1) nside_point
    
    ! Read modules
    allocate(mods(nmod))
    read(unit, err=1, end=1) mods

    ! Select the modules we actually want
    allocate(inds(nmod))
    if(present(modules)) then
       do i = 1, nmod; inds(i) = index_of(modules, mods(i)); end do
    else
       do i = 1, nmod; inds(i) = i; end do
    end if
    nget = count(inds /= 0)

    ! Allocate output structure and insert module numbers
    allocate(data(nget))
    do i = 1, nmod
       j = inds(i)
       if (j > 0) then 
          data(j)%module_number = mods(i)
       end if
    end do

    ! Read mask2map info
    allocate(npix(nmod))
    read(unit, end=1, err=1) npix
    do i = 1, nmod
       j = inds(i)
       if(j > 0) then
          allocate(data(j)%mask2map(npix(i)))
          allocate(data(j)%nhits(npix(i)))
          read(unit, err=1, end=1) data(j)%mask2map
          read(unit, err=1, end=1) data(j)%nhits
          data(j)%nside_point = nside_point
       else
          read(unit, err=1, end=1)
          read(unit, err=1, end=1)
       end if
    end do
    deallocate(npix)
    if (present(mask2map_only)) then
       if (mask2map_only) then
          ! Return without reading the rest of the data
          close(unit)
          if (present(status)) status = 0
          return
       end if
    end if

    ! Read time info
    allocate(time(nsamp))
    allocate(orig_point(3,nsamp))
    read(unit, err=1, end=1) time
    read(unit, err=1, end=1) orig_point
    do j = 1, nget
       allocate(data(j)%time(nsamp))
       allocate(data(j)%orig_point(3,nsamp))
       data(j)%time       = time
       data(j)%orig_point = orig_point
    end do
    deallocate(time)
    deallocate(orig_point)

    ! Read pointing info
    allocate(buffer(3,nsamp))
    do i = 1, nmod
       j = inds(i)
       if (j > 0) then
          read(unit, err=1, end=1) buffer
          allocate(data(j)%pointing(3,nsamp))
          data(j)%pointing = buffer
       else
          read(unit,err=1,end=1)
       end if
    end do
    deallocate(buffer)

    ! Read TOD info
    allocate(buffer(0:3,nsamp))
    do i = 1, nmod
       j = inds(i)
       if (j > 0) then
          read(unit, err=1, end=1) buffer
          allocate(data(j)%tod(0:3,nsamp))
          data(j)%tod = buffer
       else
          read(unit,err=1,end=1)
       end if
    end do
    deallocate(buffer)

    ! Read total power info
    allocate(buffer(0:3,nred))
    do i = 1, nmod
       j = inds(i)
       if (j > 0) then
          read(unit, err=1, end=1) buffer
          allocate(data(j)%tp(0:3,nred))
          data(j)%tp = buffer
          read(unit, err=1, end=1) buffer
          allocate(data(j)%tp_var(0:3,nred))
          data(j)%tp_var = buffer
       else
          read(unit,err=1,end=1)
          read(unit,err=1,end=1)
       end if
    end do
    deallocate(buffer)

    ! Close file and clean up
    close(unit)

    if(present(status)) status=0
    return

1   if(present(status)) status=1
    if(allocated(mods))       deallocate(mods)
    if(allocated(inds))       deallocate(inds)
    if(allocated(time))       deallocate(time)
    if(allocated(orig_point)) deallocate(orig_point)
    if(allocated(buffer))     deallocate(buffer)
    if(allocated(npix))       deallocate(npix)
    call deallocate_module_struct(data)
    close(unit)

  end subroutine L2_read

  subroutine L2_write(filename, data)
    implicit none
    type(module_struct), dimension(1:), intent(in) :: data
    character(len=128),                 intent(in) :: filename
    integer(i4b) :: i, nmod, num_samples, nred, red, unit
    integer(i4b),        dimension(:),  allocatable:: npix

    unit        = getlun()
    nmod        = size(data)
    num_samples = size(data(1)%time)
    red         = data(1)%reduction
    nred        = size(data(1)%tp, 2)

    open(unit, file=trim(filename), form='unformatted')

    ! Output general information numbers
    write(unit) nmod
    write(unit) num_samples
    write(unit) red
    write(unit) nred
    write(unit) data(1)%nside_point
    
    ! Output module number
    write(unit) data%module_number

    ! Output mask2map and nhits
    allocate(npix(nmod))
    do i = 1, nmod; npix(i) = size(data(i)%mask2map); end do
    write(unit) npix
    deallocate(npix)
    do i = 1, nmod
       write(unit) data(i)%mask2map
       write(unit) data(i)%nhits
    end do

    ! Output time stamps
    write(unit) data(1)%time               ! Time stamps
    write(unit) real(data(1)%orig_point,sp)

    ! Output pointing
    do i = 1, nmod
       write(unit) real(data(i)%pointing,sp)        ! Pointing information
    end do
    
    ! Output TOD
    do i = 1, nmod
       write(unit) real(data(i)%tod,sp)             ! Time ordered data       
    end do

    ! Output total power info
    do i = 1, nmod
       write(unit) real(data(i)%tp,sp)              ! Reduced tp info
       write(unit) real(data(i)%tp_var,sp)          ! Reduced tp var
    end do

    close(unit)
  end subroutine L2_write


  subroutine deallocate_module_struct_single(data)
    implicit none
    type(module_struct) :: data

    if (allocated(data%pointing))   deallocate(data%pointing)
    if (allocated(data%tod))        deallocate(data%tod)
    if (allocated(data%time))       deallocate(data%time)
    if (allocated(data%tp))         deallocate(data%tp)
    if (allocated(data%tp_var))     deallocate(data%tp_var)
    if (allocated(data%orig_point)) deallocate(data%orig_point)
    if (allocated(data%mask2map))   deallocate(data%mask2map)
    if (allocated(data%nhits))      deallocate(data%nhits)
  end subroutine deallocate_module_struct_single

  subroutine deallocate_module_struct_array(data)
    implicit none
    type(module_struct), allocatable, dimension(:) :: data
    integer(i4b) :: i

    if(.not. allocated(data)) return
    do i = 1, size(data)
       call deallocate_module_struct_single(data(i))
    end do
    deallocate(data)
  end subroutine deallocate_module_struct_array

  subroutine read_table_dp(fname, table)
    implicit none
    character(len=*),               intent(in)    :: fname
    real(dp),          allocatable, intent(inout) :: table(:,:)
    character(len=10000)                          :: line
    integer(i4b)                                  :: i, j, k, n, m, unit
    if(allocated(table)) deallocate(table)
    unit = getlun()
    open(unit,file=fname,action="read",status="old")
    n = 0; m = 0
    do
       read(unit,'(a)',end=1) line
       k = len_trim(line)
       if(k == 0 .or. line(1:1) == "#") cycle
       n = n+1
       i = num_tokens(line(1:k)," 	")
       if(m == 0) then
          m = i
       else
          call assert(i==m, "Inconsistent number of columns in " // trim(fname))
       end if
    end do
1   allocate(table(n,m))
    rewind(unit)
    do i = 1, n
       read(unit,'(a)',end=1) line
       k = len_trim(line)
       if(k == 0 .or. line(1:1) == "#") cycle
       read(line(1:k),*) table(i,:)
    end do
    close(unit)
  end subroutine

!-------------------------------------------------------------------------
! subroutine read_filelist
!-------------------------------------------------------------------------

  subroutine read_filelist(unit, filename, numfiles, filelist)
    implicit none
    
    character(len=*),                          intent(in)  :: filename
    character(len=256), dimension(:), pointer              :: filelist
    integer(i4b),                              intent(out) :: numfiles

    integer(i4b)                                           :: i, unit

    unit=getlun()

    numfiles=0
    open(unit,file=trim(filename))
    do while (.true.)
       read(unit,*,end=6)
       numfiles = numfiles + 1
    end do
6   close(unit)

    allocate(filelist(numfiles)) 
    open(unit,file=trim(filename))
    do i=1,numfiles
       read(unit,fmt='(a)') filelist(i)
    end do
    close(unit)

  end subroutine read_filelist


end module quiet_fileutils
