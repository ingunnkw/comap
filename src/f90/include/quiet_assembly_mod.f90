module quiet_assembly_mod
  use healpix_types
  use quiet_utils
  use quiet_module_mod
  use quiet_acceptlist_mod
  implicit none

  ! In general, one might want to support that each module measures an arbitrary
  ! linear combination of several different horns, just like the temperature
  ! modules measures 17-18 and 18-17 in the Q band. The number of horns going
  ! into each module might be different inside one assembly. This is purely a
  ! question of pointing, so the pointing array, pixelized pointing array
  ! and a gain structure would need to have variable lengths in that case.
  ! Hmm. We'll think about that later.

  type assembly_diode
     real(dp),    allocatable, dimension(:,:,:) :: point ! ((phi,theta,psi),nsamp,ngroup)
     integer(i4b),allocatable, dimension(:,:)   :: pix   ! (ngroup,nsamp)
     real(dp),    allocatable, dimension(:)     :: gamp  ! (ngroup)
     real(dp),    allocatable, dimension(:)     :: gain  ! (nsamp)
  end type

  type quiet_assembly
     real(dp)                                      :: tquv_frac(4)
     real(dp)                                      :: samprate ! Samprate of data in this struct
     integer(i4b)                                  :: num_diodes, numsamp, decimation
     integer(i4b)                                  :: coordinate_system, nside
     integer(i4b), allocatable, dimension(:,:)     :: diodes
     integer(i4b), allocatable, dimension(:)       :: diodes_abs

     ! Describes the diode signal
     type(assembly_diode), allocatable, dimension(:) :: diode_info

     ! Pointing
     real(dp),     allocatable, dimension(:)       :: time
     real(dp),     allocatable, dimension(:,:)     :: orig_point! ((az,el,dk),nsamp)
     integer(i4b), allocatable, dimension(:)       :: nhits     ! (0:npix-1)
     real(dp)                                      :: scanfreq

     ! Data
     real(dp),     allocatable, dimension(:,:)     :: tod       ! (nsamp, ndi)

     ! Noise properties
     real(dp),     allocatable, dimension(:)       :: sigma0, fknee, alpha ! (ndi)
     real(dp),     allocatable, dimension(:,:,:)   :: corr       ! (nfft,ndi,ndi)

     ! Internally calculated stuff.
     integer(i4b) :: noise_corr_length
     real(dp),     allocatable, dimension(:)     :: fft_low_freq, fft_low_alpha
     real(dp),     allocatable, dimension(:)     :: fft_high_freq, fft_high_alpha
     integer(i4b), allocatable, dimension(:)     :: az_order
     real(dp),     allocatable, dimension(:,:)   :: inv_N_filter, inv_N_fft
     real(dp),     allocatable, dimension(:,:,:) :: N_corr_real, N_corr_F_sq_real
     real(dp),     allocatable, dimension(:,:,:) :: N_corr_fft,  N_corr_F_sq_fft
     real(dp),     allocatable, dimension(:,:)   :: N_eff, inv_N_eff
  end type quiet_assembly

  type shared_assembly
     integer(i4b), dimension(:), allocatable :: members
     type(quiet_assembly)                    :: assembly
  end type

  interface split_assembly
     module procedure split_assembly_alist, split_assembly_logical
  end interface

  integer(i4b)                                    :: num_assemblies
  type(quiet_assembly), allocatable, dimension(:) :: assemblies

  logical(lgt), private :: initialized = .false.

contains

  subroutine initialize_quiet_assembly_mod(unit, paramfile)
    implicit none

    integer(i4b),     intent(in) :: unit
    character(len=*), intent(in) :: paramfile

    integer(i4b) :: i, j, n
    integer(i4b), allocatable, dimension(:) :: d
    character(len=256) :: assembly_file
    if(initialized) return

    call initialize_module_mod(paramfile)
    call get_parameter(unit, paramfile, 'ASSEMBLY_LIST', par_string=assembly_file)

    ! Read assembly file
    open(58,file=trim(assembly_file)) 
    read(58,*) num_assemblies
    allocate(assemblies(num_assemblies))
    do i = 1, num_assemblies
       read(58,*) n
       allocate(assemblies(i)%diodes(n,2))
       allocate(d(2*n))
       backspace(58)
       read(58,*) n, d
       do j = 1, n
          assemblies(i)%diodes(j,1) = d(2*j-1)
          assemblies(i)%diodes(j,2) = d(2*j)
       end do
       assemblies(i)%num_diodes = n
       deallocate(d)
       allocate(assemblies(i)%diodes_abs(assemblies(i)%num_diodes))
       assemblies(i)%diodes_abs = assemblies(i)%diodes(:,1)*get_num_diodes() + assemblies(i)%diodes(:,2) + 1

       ! tquv_frac gives a rough sense of how temperature- and polarization oriented
       ! an assembly is. Should ideally be S/N-weighted, but bah. This will be used
       ! when considering whether to skip an assembly or not.
       assemblies(i)%tquv_frac = 0
       do j = 1, assemblies(i)%num_diodes
          assemblies(i)%tquv_frac = assemblies(i)%tquv_frac + &
           & abs(quiet_diodes(assemblies(i)%diodes_abs(j))%stokes([1,2,2,3]))
       end do
       assemblies(i)%tquv_frac = assemblies(i)%tquv_frac / assemblies(i)%num_diodes
    end do
    close(58)
  end subroutine initialize_quiet_assembly_mod

  subroutine cleanup_quiet_assembly_mod
    implicit none
    integer(i4b) :: i

    do i = 1, num_assemblies
       if (allocated(assemblies(i)%diodes)) deallocate(assemblies(i)%diodes)
    end do
    if (allocated(assemblies)) deallocate(assemblies)
  end subroutine cleanup_quiet_assembly_mod

  subroutine get_quiet_assembly(alist, cid, id, assembly)
    implicit none
    type(acceptlist)                  :: alist
    integer(i4b),         intent(in)  :: id, cid
    type(quiet_assembly), intent(out) :: assembly

    integer(i4b) :: i, m, num_accept

    m = 0
    do i = 1, assemblies(id)%num_diodes
       if (is_accepted(alist, cid, assemblies(id)%diodes(i,1), assemblies(id)%diodes(i,2))) then
          m = m + 1
       end if
    end do

    if (m == 0) then
       assembly%num_diodes = 0
       return
    end if

    allocate(assembly%diodes(m,2), assembly%diodes_abs(m))
    assembly%num_diodes = m
    m = 0
    do i = 1, assemblies(id)%num_diodes
       if (is_accepted(alist, cid, assemblies(id)%diodes(i,1), assemblies(id)%diodes(i,2))) then
          m = m + 1
          assembly%diodes(m,:)   = assemblies(id)%diodes(i,:)
          assembly%diodes_abs(m) = assemblies(id)%diodes_abs(i)
       end if
    end do
    assembly%tquv_frac = assemblies(id)%tquv_frac
  end subroutine get_quiet_assembly

  subroutine deallocate_quiet_assembly(assembly)
    implicit none
    type(quiet_assembly), intent(inout) :: assembly

    integer(i4b) :: i

    if (allocated(assembly%diodes))               deallocate(assembly%diodes)
    if (allocated(assembly%diodes_abs))           deallocate(assembly%diodes_abs)
    if (allocated(assembly%time))                 deallocate(assembly%time)
    if (allocated(assembly%orig_point))           deallocate(assembly%orig_point)

    if (allocated(assembly%tod))                  deallocate(assembly%tod)

    if (allocated(assembly%fft_low_freq))         deallocate(assembly%fft_low_freq)
    if (allocated(assembly%fft_low_alpha))        deallocate(assembly%fft_low_alpha)
    if (allocated(assembly%fft_high_freq))        deallocate(assembly%fft_high_freq)
    if (allocated(assembly%fft_high_alpha))       deallocate(assembly%fft_high_alpha)
    if (allocated(assembly%az_order))             deallocate(assembly%az_order)
    if (allocated(assembly%sigma0))               deallocate(assembly%sigma0)
    if (allocated(assembly%fknee))                deallocate(assembly%fknee)
    if (allocated(assembly%alpha))                deallocate(assembly%alpha)
    if (allocated(assembly%corr))                 deallocate(assembly%corr)
    if (allocated(assembly%N_corr_real))          deallocate(assembly%N_corr_real)
    if (allocated(assembly%N_corr_F_sq_real))     deallocate(assembly%N_corr_F_sq_real)
    if (allocated(assembly%N_corr_fft))           deallocate(assembly%N_corr_fft)
    if (allocated(assembly%N_corr_F_sq_fft))      deallocate(assembly%N_corr_F_sq_fft)
    if (allocated(assembly%inv_N_fft))            deallocate(assembly%inv_N_fft)
    if (allocated(assembly%inv_N_filter))         deallocate(assembly%inv_N_filter)
    if (allocated(assembly%N_eff))                deallocate(assembly%N_eff)
    if (allocated(assembly%inv_N_eff))            deallocate(assembly%inv_N_eff)

    call free_assembly_diodes(assembly%diode_info)
  end subroutine deallocate_quiet_assembly

  ! The mapmaking routines need consistent assemblies, but with autojackknifing, we
  ! might not get that. Therefore, inside the assembly loop of tod2map, we will
  ! use the jackknifed accepts to divide into groups with compatible assemblies,
  ! and then loop over those. This routine performes that division.

  subroutine split_assembly_alist(assembly, alists, cid, groups)
    implicit none
    type(quiet_assembly) :: assembly
    type(acceptlist)     :: alists(:)
    type(shared_assembly), dimension(:),   allocatable :: groups
    logical(lgt),          dimension(:,:), allocatable :: accepts
    integer(i4b) :: i, m, n, cid, cnum
    n    = alists(1)%nmod*alists(1)%ndi
    m    = size(alists)
    cnum = lookup_ces(cid)
    allocate(accepts(n,m))
    do i = 1, m
       accepts(:,i) = reshape(alists(i)%status(:,:,cnum),[n]) == REJECTED_NONE
    end do
    call split_assembly_logical(assembly, accepts, groups)
    deallocate(accepts)
  end subroutine

  subroutine split_assembly_logical(assembly, accepts, groups)
    implicit none
    logical(lgt)         :: accepts(:,:)
    type(quiet_assembly) :: assembly
    integer(i4b)         :: i, j, k, m, n, diodes(assembly%num_diodes), ndi, rep
    integer(i4b)         :: ng(size(accepts,2)), g(size(accepts,2),size(accepts,2))
    logical(lgt)         :: done(size(accepts,2)), empty(size(accepts,2))
    type(shared_assembly), dimension(:), allocatable :: groups

    ! Find groups. For each diode in the assembly, check which acceptlists,
    ! agree on which are accepted.
    ndi    = get_num_diodes()
    m      = 0
    ng     = 0
    diodes = assembly%diodes_abs
    empty  = all(.not. accepts(diodes,:), 1)
    done   = empty
    do i = 1, size(accepts,2)
       if(done(i)) cycle
       m      = m+1
       do j = i, size(accepts,2)
          if(.not. all(accepts(diodes,i).eqv.accepts(diodes,j))) cycle
          ng(m)      = ng(m) + 1
          g(ng(m),m) = j
          done(j)    = .true.
       end do
    end do

    ! Populate group struct
    call free_shared_assemblies(groups)
    allocate(groups(m))
    do i = 1, m
       allocate(groups(i)%members(ng(i)))
       groups(i)%members = g(1:ng(i),i)
       rep = groups(i)%members(1)
       groups(i)%assembly%num_diodes = count(accepts(diodes,rep))
       allocate(groups(i)%assembly%diodes(groups(i)%assembly%num_diodes,2))
       allocate(groups(i)%assembly%diodes_abs(groups(i)%assembly%num_diodes))
       j = 0
       do k = 1, size(diodes)
          if(.not. accepts(diodes(k),rep)) cycle
          j = j + 1
          groups(i)%assembly%diodes(j,:)   = [ (diodes(k)-1)/ndi, mod((diodes(k)-1),ndi) ]
          groups(i)%assembly%diodes_abs(j) = diodes(k)
       end do
       groups(i)%assembly%tquv_frac = assembly%tquv_frac
    end do
  end subroutine

  subroutine free_shared_assemblies(a)
    implicit none
    type(shared_assembly), dimension(:), allocatable :: a
    integer(i4b) :: i
    if(.not. allocated(a)) return
    do i = 1, size(a)
       call free_shared_assembly(a(i))
    end do
    deallocate(a)
  end subroutine

  subroutine free_shared_assembly(a)
    implicit none
    type(shared_assembly) :: a
    if(allocated(a%members)) deallocate(a%members)
    call deallocate_quiet_assembly(a%assembly)
  end subroutine

  subroutine free_assembly_diode(d)
    implicit none
    type(assembly_diode) :: d
    if(allocated(d%point)) deallocate(d%point)
    if(allocated(d%pix))   deallocate(d%pix)
    if(allocated(d%gamp))  deallocate(d%gamp)
    if(allocated(d%gain))  deallocate(d%gain)
  end subroutine

  subroutine free_assembly_diodes(d)
    implicit none
    type(assembly_diode), allocatable, dimension(:) :: d
    integer(i4b) :: i
    if(allocated(d)) then
       do i = 1, size(d)
          call free_assembly_diode(d(i))
       end do
       deallocate(d)
    end if
  end subroutine

  subroutine copy_assembly_diode(a,b)
    implicit none
    type(assembly_diode) :: a, b
    call free_assembly_diode(b)
    allocate(b%point(size(a%point,1),size(a%point,2),size(a%point,3)))
    allocate(b%pix(size(a%pix,1),size(a%pix,2)))
    allocate(b%gamp(size(a%gamp)))
    allocate(b%gain(size(a%gain)))
    b = a
  end subroutine

end module quiet_assembly_mod
