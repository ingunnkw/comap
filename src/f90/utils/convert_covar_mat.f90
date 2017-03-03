program convert_covar_mat
  use healpix_types
  use quiet_fileutils
  implicit none

  integer(i4b)       :: i, n, ordering, unit, pol_flag
  integer(i8b)       :: m
  logical(lgt)       :: inv_cov
  character(len=128) :: infile, outfile, partext
  real(dp), allocatable, dimension(:,:) :: covar

  unit = 38

  if (iargc() /= 5) then
     write(*,*) 'Usage: convert_covar_mat [infile] [outfile] [ordering] [polarization status] [inv_covar]'
     stop
  end if

  call getarg(1,infile)
  call getarg(2,outfile)
  call getarg(3,partext)
  read(partext,*) ordering
  call getarg(4,partext)
  read(partext,*) pol_flag
  call getarg(5,partext)
  read(partext,*) inv_cov

  ! Read covariance matrix
  open(unit,file=trim(infile), form='unformatted')
  read(unit) n
  m = n
  allocate(covar(m,m))
  read(unit) covar
  close(unit)

  ! Write covariance matrix in new format
  call write_cov_mat(unit, outfile, ordering, pol_flag, inv_cov, covar, .false.)

end program convert_covar_mat
