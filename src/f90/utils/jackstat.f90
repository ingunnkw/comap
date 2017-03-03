program jackstat
  use quiet_utils
  implicit none
  character(len=512)   :: arg, file
  character(len=10000) :: line
  integer(i4b)         :: nces, ndi, n, i, j, k, l, m, nj, unit, cid, counts(2), nok
  integer(i4b)         :: na, nb
  integer(i4b), allocatable :: cdis(:), knives(:,:,:), tmpi(:), map(:)
  real(dp),     allocatable :: jmat(:,:), corr(:,:), tmat(:,:)
  nj = iargc()
  do i = 1, nj
     call dmem(itoa(i))
     call getarg(i, file)
     unit = getlun()
     if(i == 1) then
        ! First one - measure nces and ndi
        open(unit, file=file, action="read", status="old")
        nces = 0
        do
           if(nces == 0) then
              read(unit,'(a)',end=1) line
              ndi = num_tokens(line, " ")-1
              allocate(cdis(ndi))
              read(line,*) cid, cdis
           else
              read(unit,*,end=1) cid, cdis
           end if
           nces = nces+1
        end do
        1 close(unit)
        allocate(knives(ndi,nces,nj))
     end if
     open(unit, file=file, action="read", status="old")
     do j = 1, nces
        read(unit,*) cid, knives(:,j,i)
     end do
     close(unit)
  end do
  ! Ok, we now have all the knives. Create a map of the useful ones.
  ! These should have exactly two subdivisions.
  j = 0
  allocate(tmpi(nj))
  do i = 1, nj
     if(maxval(knives(:,:,i)) /= 2) cycle
     if(.not. any(knives(:,:,i)==1)) cycle
     j = j+1
     tmpi(j) = i
  end do
  allocate(map(j)); map = tmpi(:j); deallocate(tmpi)
  nok = size(map)
  call dmem("map")
  ! Now create the jack matrix
  allocate(jmat(nok,ndi*nces))
  do i = 1, nok
     na = count(knives(:,:,map(i))==1)
     nb = count(knives(:,:,map(i))==2)
     do j = 1, nces
        do k = 1, ndi
           l = (j-1)*ndi+k
           select case(knives(k,j,map(i)))
              case(1); jmat(i,l) =  1d0/na
              case(2); jmat(i,l) = -1d0/nb
           end select
        end do
     end do
  end do
  call dmem("jack")
  ! And the jackknife correlation matrix
  allocate(corr(nok,nok), tmat(nok,nok))
  tmat = matmul(jmat,transpose(jmat))
  do i = 1, nok
     do j = 1, nok
        corr(j,i) = tmat(j,i)/sqrt(tmat(i,i)*tmat(j,j))
     end do
  end do
  deallocate(tmat)
  call dmem("corr")
  call dump_matrix(corr,fmt="(e15.7)")

end program
