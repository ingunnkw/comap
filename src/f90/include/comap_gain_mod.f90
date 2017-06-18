module comap_gain_mod
  use quiet_utils
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

  subroutine estimate_gain(el,dat,g)
    implicit none
    real(dp), dimension(:)              :: el, dat
    real(dp),               intent(out) :: g
    real(dp), dimension(:), allocatable :: templ
    integer(i4b)                        :: n

    n = size(el)
    allocate(templ(n))
    templ = 1/(sin(el*pi/180.))
    !write(*,*) sum(templ*templ)

    g = sum(templ*dat)/sum(templ*templ)

  end subroutine estimate_gain

end module comap_gain_mod

