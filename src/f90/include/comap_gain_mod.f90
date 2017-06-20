module comap_gain_mod
  use quiet_utils
  use math_tools
!  use quiet_postutils
  implicit none

contains

  subroutine initialize_gain_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    !call get_parameter(0, parfile, "NUMFREQ", par_int=numfreq)
    initialized = .true.
  end subroutine

  subroutine estimate_gain(el,dat,g,sigma0)
    implicit none
    real(dp), dimension(:)                :: el, dat
    real(dp), dimension(:),   optional    :: sigma0
    real(dp),                 intent(out) :: g
    real(dp), dimension(:,:), allocatable :: templ, data, mat, hm
    real(dp), dimension(:),   allocatable :: amp
    integer(i4b)                          :: n, npol, nt, i, j 
    real(dp)                              :: chisq

    npol = 3
    nt = 1 + npol
    n = size(el)
    allocate(templ(n,nt), data(n,1), amp(nt), mat(nt,nt), hm(nt,1))
    templ(:,1) = 1/(sin(el*pi/180.))
    do j = 1, npol
       do i = 1, n
          templ(i,j+1) = (real(i,dp)-(real(n,dp)/2.))**(j-1)
       end do
    end do
    data(:,1) = dat
    mat = matmul(transpose(templ),templ)
    hm =  matmul(transpose(templ), data)
    call  solve_linear_system_dp(mat, amp, hm(:,1))
    g = amp(1)

!    if (exists(sigma0)) then
!       hm(:,1) = amp
!       data = (data -matmul(templ,hm))/sigma0
!       chisq = sum(data(:,1))
!    end if

  end subroutine estimate_gain

end module comap_gain_mod

