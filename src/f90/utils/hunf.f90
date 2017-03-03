! convert between unf-files and hdf-files, since hdf-files are easier
! to work with.
program hunf
  use quiet_hdf_mod
  use quiet_fileutils
  implicit none
  character(len=512)                   :: arg, cmd, ifile, ofile
  real(dp),              allocatable   :: mat(:,:)
  integer(i4b),          allocatable   :: pixels(:)
  integer(i4b)                         :: nside, ncomp, order
  logical(lgt)                         :: inv
  call getarg(1, cmd)
  select case(cmd)
     case("cov")
        call getarg(2, ifile)
        call getarg(3, ofile)
        call read_covmat (mat, pixels, nside, order, ncomp, ifile, inv, verbose=.true.)
        call write_covmat(mat, pixels, nside, order, ncomp, ofile, inv, verbose=.true.)
     case default
        call help
  end select

contains
   subroutine help
     implicit none
     write(*,'(a)') "Convert between unf and hdf"
     write(*,'(a)') " hunf cov ifile ofile"
     stop
   end subroutine

end program
