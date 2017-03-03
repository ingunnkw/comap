program euler_prog
  use math_tools
  use quiet_utils
  implicit none
  character(len=512) :: line, arg, groups(256), str, angs(3)
  real(dp)           :: mat(3,3), tot(3,3), phi, theta, psi
  integer(i4b)       :: i, j, k, m, n, ng
  logical(lgt)       :: trans
  do
     read(*,'(a)',end=1) line
     call get_tokens(line, " 	,", groups, ng, '()[]', size(groups))
     tot = get_identity(3)
     do i = 1, ng
        str = groups(i)
        do j = 1, len_trim(str)
           if(str(j:j) == ')') exit
        end do
        call get_tokens(str(2:j-1), ",", angs, n, maxnum=size(angs))
        phi   = atof(angs(1))*DEG2RAD
        theta = atof(angs(2))*DEG2RAD
        psi   = atof(angs(3))*DEG2RAD
        call compute_euler_matrix_zyz(phi, theta, psi, mat)
        if(index(str(j+1:),"'") > 0) mat = transpose(mat)
        tot = matmul(mat, tot)
     end do
     call convert_euler_matrix_to_angles_zyz(tot, phi, theta, psi)
     write(*,'(a)') "(" // trim(ftoa(phi*RAD2DEG,4)) // "," // &
      & trim(ftoa(theta*RAD2DEG,4)) // "," // trim(ftoa(psi*RAD2DEG,4)) // ")"
  end do
1 continue
end program
