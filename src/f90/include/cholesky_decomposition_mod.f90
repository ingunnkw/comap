module cholesky_decomposition_mod
  use healpix_types 
  
  contains 

    subroutine cholesky_decomposition(A, n, G)
      implicit none   

      integer(i4b), intent(in)       :: n         ! Number of rows/cols in matrix
      real(dp),     intent(out)      :: G(n,n)    
      real(dp),     intent(inout)    :: A(n,n)    

      integer(i4b) :: i,j  


      G(:,:)=0.0d0
      do j = 1, n
         if (A(j,j) == 0.d0) cycle
         G(j,j) = sqrt( A(j,j) - dot_product(G(j,1:j-1),G(j,1:j-1)) )
         do i = j+1, n
            G(i,j)  = ( A(i,j) - dot_product(G(i,1:j-1),G(j,1:j-1)) ) / G(j,j)
         end do
      end do



    end subroutine cholesky_decomposition
  
  
end module cholesky_decomposition_mod
