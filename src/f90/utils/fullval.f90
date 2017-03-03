program tester
   use quiet_utils
   use quiet_pointing_mod
   implicit none

   type obs_type
      integer(i4b) :: ces, di, horn
      real(dp)     :: mjd, p0(3), p(3), tamp, amp, frac, dist, fwhm, foo(6), ddk
   end type

   call fullval

contains

   subroutine fullval
     implicit none
     character(len=512) :: arg, ifile
     integer(i4b)       :: unit, cmax, dimax, i, j, h, mod, mod2, di, di2
     type(obs_type)     :: obs, hobs(2), pobs(2)
     integer(i4b),   allocatable :: diobs(:)
     type(obs_type), allocatable :: obslist(:,:,:), foolist(:,:,:)
     call getarg(1, ifile)
     unit = getlun()
     cmax = 0; dimax = 0
     open(unit,file=ifile,action="read",status="old")
     do
        read(unit,*,end=1,err=4) obs
        cmax  = max(cmax,  obs%ces)
        dimax = max(dimax, obs%di)
        4 continue
     end do
1    rewind(unit)
     allocate(diobs(dimax), obslist(cmax,2,dimax), foolist(cmax,2,dimax))
     obslist%ces = 0
     foolist%ces = 0
     diobs = 0
     do
        read(unit,*,end=2,err=5) obs
        if(obs%horn == 1) diobs(obs%di) = diobs(obs%di) + 1
        obslist(diobs(obs%di),obs%horn, obs%di) = obs
        foolist(obs%ces,      obs%horn, obs%di) = obs
        5 continue
     end do
2    close(unit)

     do di = 1, dimax
        do j = 1, diobs(di)
           hobs = obslist(j,:,di)
           ! Reject if horn ratio is not -1
           if(abs(hobs(2)%frac+1) > 0.3) cycle
           ! Reject if gain is too high
           if(any(abs(hobs%amp) > 40)) cycle
           if(any(abs(hobs%amp) == 0)) cycle
           ! Reject if positions are too far away from center
           if(polangdist([pi/2,0d0]-hobs(1)%p([2,1])*DEG2RAD, &
            & [pi/2,0d0]-hobs(1)%p0([2,1])*DEG2RAD) > 0.6*DEG2RAD) cycle
           if(polangdist([pi/2,0d0]-hobs(2)%p([2,1])*DEG2RAD, &
            & [pi/2,0d0]-hobs(2)%p0([2,1])*DEG2RAD) > 0.6*DEG2RAD) cycle
           mod  = (di-1)/4
           if(modulo(mod,2) == 0) then; mod2 = mod-1; else; mod2 = mod+1; end if
           di2  = mod2*4+1
           pobs = foolist(hobs(1)%ces,:,di2)
           if(pobs(1)%ces /= hobs(1)%ces) cycle

           ! Reject if horns don't agree
           if(polangdist([pi/2,0d0]-hobs(1)%p([2,1])*DEG2RAD, &
            & [pi/2,0d0]-pobs(1)%p([2,1])*DEG2RAD) > 7*DEG2RAD/60) cycle

           ! And output the current entry
           3 format(i5,i4,i3,f14.7,6f10.4,3e15.7,2f10.4,6e15.7,f10.4)
           do h = 1, 2
              write(*,3) hobs(h)%ces, hobs(h)%di, hobs(h)%horn, hobs(h)%mjd, &
               & hobs(h)%p0, hobs(h)%p, hobs(h)%tamp, hobs(h)%amp, hobs(h)%frac, &
               & hobs(h)%dist, hobs(h)%foo, hobs(h)%ddk
           end do
        end do
     end do
   end subroutine

end program
