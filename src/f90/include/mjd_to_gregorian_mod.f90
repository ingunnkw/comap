module mjd_to_gregorian_mod 
  use healpix_types
  
  contains 
    
    subroutine mjd_to_gregorian(mjd,year,month,day,hour,minute,second)
      
      implicit none
      real(dp), intent(in) :: mjd
      integer(i4b), intent(out) :: year, month, day, hour, minute, second
      real(dp) :: jt, jd
      integer(i4b) :: i,j,k,l,n

      jd = mjd + 2400000.5d0

      l = int(jd)+68569
      n = 4*l/146097
      l = l-(146097*n+3)/4
      i = 4000*(l+1)/1461001
      l = l-1461*i/4+31
      j = 80*l/2447
      k = l-2447*j/80
      l = j/11
      j = j+2-12*l
      i = 100*(n-49)+i+l

      year = i
      month = j
      day = k

      jt = dmod(jd,1.d0)*24.d0
      hour = int(jt)
      jt = dmod(jt,1.d0)*60.d0
      minute = int(jt)
      jt = dmod(jt,1.d0)*60.d0
      second = nint(jt)

      if (second == 60) then
         second = second-60
         minute = minute+1
      end if

    end subroutine mjd_to_gregorian

  end module mjd_to_gregorian_mod
