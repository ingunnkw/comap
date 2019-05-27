module comap_ephem_mod
  use healpix_types
  use spline_1D_mod
  use quiet_hdf_mod
  implicit none

  character(len=512), dimension(:), allocatable, private :: obj_list
  real(dp), dimension(:,:,:), allocatable, private       :: obj_data, obj_data2    
  integer(i4b), private                                  :: n_objects

contains
  
  subroutine initialize_comap_ephem_mod(parfile)
    character(len=512), intent(in) :: parfile
    character(len=512)             :: ephem_file
    type(hdf_file)                 :: file
    integer(i4b)                   :: i, arr_dim(2)
    n_objects = 1
    allocate(obj_list(n_objects))
    obj_list(1) = 'jupiter'

    call get_parameter(0, parfile,  'EPHEMERIS_FILE', par_string=ephem_file)
    
    call open_hdf_file(ephem_file, file, "r")
    call get_size_hdf(file, obj_list(1), arr_dim)
    allocate(obj_data(n_objects, arr_dim(1), arr_dim(2)))
    allocate(obj_data2(n_objects, arr_dim(1) - 1, arr_dim(2)))
    do i = 1, n_objects
       call read_hdf(file, obj_list(i), obj_data(i, :, :))
       call spline(obj_data(i, 1, :), obj_data(i, 2, :), 1.d30, 1.d30, obj_data2(i, 2, :))
       call spline(obj_data(i, 1, :), obj_data(i, 3, :), 1.d30, 1.d30, obj_data2(i, 3, :))
       call spline(obj_data(i, 1, :), obj_data(i, 4, :), 1.d30, 1.d30, obj_data2(i, 4, :))
    end do
    call close_hdf_file(file)
  end subroutine initialize_comap_ephem_mod
  
  function get_obj_info(obj, mjd) result(res)
    character(len=*),    intent(in)   :: obj
    real(dp),            intent(in)   :: mjd
    real(dp)                          :: res(3), ra, dec, dist
    integer(i4b)                      :: i, j(1)
    j = FINDLOC(obj_list, VALUE=obj)
    i = j(1)
    ra   = splint(obj_data(i, 1, :),obj_data(i,2,:), obj_data2(i,2,:),mjd)
    dec  = splint(obj_data(i, 1, :),obj_data(i,3,:), obj_data2(i,3,:),mjd)
    dist = splint(obj_data(i, 1, :),obj_data(i,4,:), obj_data2(i,4,:),mjd)
    res(1) = ra
    res(2) = dec
    res(3) = dist
  end function get_obj_info

end module comap_ephem_mod
