module comap_frequency_mod
  use quiet_utils
  implicit none
  
  integer(i4b)                 :: numfreq

contains

  subroutine init_frequency_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    call get_parameter(0, parfile, "NUMFREQ", par_int=numfreq)
    initialized = .true.
  end subroutine

  function get_num_freqs() result(res)
    implicit none
    integer(i4b) :: res
    res = numfreq
  end function

end module comap_frequency_mod

