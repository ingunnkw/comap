module comap_map_mod
  use healpix_types
  use quiet_utils
  implicit none

  type map_type
     integer(i4b) :: n_x, n_y, nfreq, n_k, ntheta ! 2^ntheta
     real(dp)     :: x0, y0, f0, df
     real(dp)     :: dtheta
     real(dp), allocatable, dimension(:)     :: x, y, f, k ! (n_x or n_y or nfreq or n_k)
     real(dp), allocatable, dimension(:,:,:) :: m, rms, dsum, nhit, div ! (n_x, n_y, nfreq)
  end type map_type


!contains

! Something

end module comap_map_mod
