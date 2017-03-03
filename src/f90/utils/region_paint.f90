! Executes a set of painting commands in a file, producing a
! 1-component output map.
program region_paint
  use quiet_utils
  use quiet_fileutils
  implicit none
  character(len=512) :: cmdfile, arg, ofname, cmd, subdir
  integer(i4b)     :: nside, ordering, unit, unit2, n, npix, isnest, i, teller, maxpix, s, ismask, jj
  real(dp)         :: val, vec(3), lat, lat1, lat2, lon, lon1, lon2, r, sub_val
  real(dp)         :: vertices(1:3,0:3), dlat, dlon, theta, phi, gmllon, h, v, grad, pmul
  real(dp),     dimension(:,:), allocatable :: map, sub_map, mask
  integer(i4b), dimension(:),   allocatable :: rpix, spix, pixels

  pmul = 0
  if(iargc() < 4) then
     write(*,*) "region_paint cmdfile.txt nside order outmap.fits"
     stop
  end if

  call getarg(1, cmdfile)
  call getarg(2, arg); read(arg,*) nside
  call getarg(3, arg)
  call getarg(4, ofname)
  select case(arg)
     case('ring'); ordering = 1; isnest = 0
     case('nest'); ordering = 2; isnest = 1
     case default
        write(*,*) "Unrecognized ordering " // trim(arg) // "!"
        stop
  end select
  do i = 5, iargc()
     call getarg(i, arg)
     select case(arg)
        case("add"); pmul = 1
        case default
           write(*,*) "Unrecognized option " // trim(arg) // "!"
           stop
     end select
  end do

  ! Set up a drawing surface, and begin to draw.
  npix = 12*nside**2
  allocate(map(0:npix-1,1), rpix(npix))
  map = 0

  ! Parse the command file while drawing
  unit = getlun()
  open(unit,file=cmdfile,action="read",status="old")
  do
     read(unit,*,end=1) val, cmd
     backspace(unit)
     select case(cmd)
        case("disc")
           read(unit,*) val, cmd, lat, lon, r
           call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
           call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
           map(rpix(1:n),1) = val + pmul*map(rpix(1:n),1)
       case("strip")
           read(unit,*) val, cmd, lat, lat2
           call query_strip(nside, pi/2-lat2*DEG2RAD, pi/2-lat*DEG2RAD, rpix, n, nest=isnest)
           map(rpix(1:n),1) = val + pmul*map(rpix(1:n),1)
       case("square","rect")
           read(unit,*) val, cmd, lat, lat2, lon, lon2
           call ang2vec(pi/2 - lat*DEG2RAD,  lon *DEG2RAD, vertices(:,0))
           call ang2vec(pi/2 - lat*DEG2RAD,  lon2*DEG2RAD, vertices(:,1))
           call ang2vec(pi/2 - lat2*DEG2RAD, lon *DEG2RAD, vertices(:,2))
           call ang2vec(pi/2 - lat2*DEG2RAD, lon2*DEG2RAD, vertices(:,3))
           call query_triangle(nside,vertices(:,0),vertices(:,1),vertices(:,2),rpix, n, nest=isnest)
           map(rpix(1:n),1) = val + pmul*map(rpix(1:n),1)
           call query_triangle(nside,vertices(:,1),vertices(:,2),vertices(:,3),rpix, n, nest=isnest)
           map(rpix(1:n),1) = val + pmul*map(rpix(1:n),1)
       ! Flat-sky rectangle, i.e. bounded by lines of constant lat and lon
       case("flat_rect")
           read(unit,*) val, cmd, lat1, lat2, lon1, lon2
           call query_strip(nside, pi/2-lat2*DEG2RAD, pi/2-lat1*DEG2RAD, rpix, n, nest=isnest)
           do i = 1, n
              call pix2ang(nside, ordering, rpix(i), theta, phi)
              dlon = modulo(phi*RAD2DEG-lon1,360.d0)
              if(dlon < lon2-lon1) map(rpix(i),1) = val + pmul*map(rpix(i),1)
           end do
       case("strip_square")
           read(unit,*) val, cmd, dlat, dlon, maxpix ,ismask
           ! dlat is width of borders in latitude (theta).
           ! dlon should be small, the width is actually dependent on maxpix
           ! maxpix is the threshold number of pixels in regions, i.e. width can be larger than dlon.
           ! ismask=1 if you want a mask. else 0
           if (ismask==1) call read_map   (mask, ordering, "sparsity.fits")
           allocate(sub_map(0:npix-1,1), spix(npix))
           sub_map = 0
           sub_val = val
           lat=90.d0
           lat2=90.d0
           teller=0
           s = 1
           do while(lat.GT.-90)
write(*,*)s
              lat2= lat
              lat = lat2-dlat
              if (lat.lt.-90) lat=-90.d0
             call query_strip(nside, pi/2-lat2*DEG2RAD, pi/2-lat*DEG2RAD, rpix, n, nest=isnest)
              h = 0.d0
              v = 0.d0
              gmllon=v
              do jj=1,nint(360./dlon+1)
!write(*,*)jj
                 h= v
                 v = h + dlon
                 do i=1,n
                    if(ismask==1 .and. .not. mask(rpix(i),1)>0.5) cycle  !if not in mask, drop pixel
                    call pix2ang_nest(nside,rpix(i),theta,phi) !phi=[0,2pi]
                    phi = phi*RAD2DEG
!write(*,*)phi
                    if (phi>180.d0) then
                       grad=phi-180.d0
                    else
                       grad=phi+180.d0
                    end if
                    if ((grad.GE.h).AND.(grad.LT.v)) then
!!$                       if(grad> 29.+gmllon) then    !uncomment to split regions 
!!$                          val=val+1               !that have to separate regions
!!$                          teller=0
!!$                       end if
                       gmllon=v
                       map(rpix(i),1) = val
                       sub_map(rpix(i),1) = sub_val
                       spix(i) = rpix(i)
                       teller = teller + 1
                    end if
                 end do
                 if(teller.GT.maxpix) then  !this is the criterium which decides the region size, in pixels.
                    val = val+1
                    sub_val = sub_val+1
                    teller = 0
                 end if
              end do  !loop stripe
              
              if (teller.LT.400) then  
                 !if last region is too small, we combine it with the one before.
                 do i = 1, n
                    if (map(rpix(i),1) == val) then
                       map(rpix(i),1) = val - 1
                       sub_map(rpix(i),1) = sub_val -1
                    end if
                 end do
                 sub_val = sub_val -1
              else if (teller.NE.0) then
                 val = val + 1 
              end if
              teller = 0

              !at the end of every stripe, create subdirs and populate with regionfiles:
              do i = 1, npix; spix(i) = i-1; end do
              subdir = "s" // trim(itoa(s)) // "/sparsity.fits"
              call mkdirs(trim(subdir), .true.)
              call write_map(sub_map, spix, nside, ordering, trim(subdir) )
              unit2 = getlun()
              open(unit2,file="s" // trim(itoa(s)) // "/regions.txt", action="write")
              write(unit2,*) sub_val
              close(unit2)
              spix = 0
              sub_map=0.d0
              sub_val=1
              s = s+1
          end do  !loop over stripes
           deallocate(sub_map,spix)
       case("full")
           read(unit,*) val, cmd
           map(:,1) = val + pmul*map(:,1)
       case default
          write(*,*) "Shape '" // trim(cmd) // "' is not supported!"
          stop
     end select
  end do
1 close(unit)
  do i = 1, npix; rpix(i) = i-1; end do
  call write_map(map, rpix, nside, ordering, trim(ofname))
  deallocate(map, rpix)
end program
