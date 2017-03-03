module quiet_pixspace_mod
  use quiet_utils
  use quiet_fileutils
  use math_tools
  use spline_1d_mod
  implicit none

  ! This struct contains lots of information about the pixels. It is
  ! several times the size of a sparse map, but much smaller than a
  ! covariance matrix. It takes some time to build, and uses about 100*n*m
  ! bytes of memory, but it saves time for later operations.
  ! Explanation:
  !  n       Number of present pixels
  !  m       Max number of neighbors
  !  order   Healpix ordering
  !  nside   Healpix nside
  !  npix    Full sky number of pixels
  !  nneigh  Actual number of neighbors considered for each pixel
  !  neigh   Pixel index of neigbhors of each pixel, up to radius rmax
  !  pix     Pixel index to pixel number
  !  map     Pixel number to pixel index
  !  r       Radius from each pixel to all its neighbors
  !  rot     Rotation matrix to each pixel from all its neighbors
  !  rmax    Maximum radius considered.
  !  dir     Tangent-space direction from each pixel to all its neighbors.
  type pixinfo
     integer(i4b)                                  :: n, m, order, nside, npix  ! 20
     integer(i4b), dimension(:),       allocatable :: nneigh, pix, map          ! 12n
     integer(i4b), dimension(:,:),     allocatable :: neigh, nme                ! 8nm
     real(dp),     dimension(:,:),     allocatable :: r, vec, polar             ! 24nm
     real(dp),     dimension(:,:,:,:), allocatable :: rot                       ! 72nm
     real(dp),     dimension(:,:,:),   allocatable :: dir                       ! 16nm
     real(dp)                                      :: rmax, pixside             ! 16
  end type

  type degrade_info
     integer(i4b)                                  :: n, m
     integer(i4b), dimension(:),       allocatable :: nsub, ipix, opix
     integer(i4b), dimension(:,:),     allocatable :: submap
     real(dp),     dimension(:,:,:,:), allocatable :: rot
  end type

  ! Map smoothing. smooth_map is general, while smooth_map_temp is specialized for temperature
  interface smooth_map
     module procedure smooth_map_pol
  end interface

  ! General noise smoothing. Works for 1, 2 or 3 components.
  interface smooth_noise
     module procedure smooth_noise_diag_pol, smooth_noise_block_pol, smooth_noise_full_pol
  end interface

  ! Temperature specialization
    real(dp),     dimension(:,:,:,:), allocatable :: ocov
  interface smooth_noise_temp
     module procedure smooth_noise_diag_temp, smooth_noise_full_temp
  end interface

  ! Map degradation. degrade_map is general, degrade_map_temp is a specialization
  interface degrade_map
     module procedure degrade_map_pol
  end interface

  ! General noise degradation. Works for 1, 2 or 3 components
  interface degrade_noise
     module procedure degrade_noise_diag_pol, degrade_noise_block_pol, degrade_noise_full_pol
  end interface

  ! Temperature specialization
  interface degrade_noise_temp
     module procedure degrade_noise_diag_temp, degrade_noise_full_temp
  end interface

  real(dp), private :: beamtol = 1e-3

contains

  ! This fills in the pixinfo structure, which contains useful
  ! information for doing convolutions in pixel space. This is
  ! an O(nm) operation, but it is probably rather heavy due to
  ! pix2vec, query_disc etc, and their implicit trigonometric
  ! operations. This step is about twice as fast in ring ordering
  ! than in nest.
  subroutine init_pixinfo(pinfo, nside, order, maxrad, pixels)
    implicit none
    type(pixinfo) :: pinfo
    integer(i4b)  :: nside, order, i, j, k, m, n, npix, nneigh, isnest, nn, pixi, pixj
    real(dp)      :: maxrad, center(3), pos(3), matN1(3,3), matN2(3,3)
    real(dp)      :: vi(3), vj(3), s2p, c2p, t1, t2, effrad
    real(dp)      :: north(3), east(3), diffs(3,2), diff(3), basis(3,2), tmp
    real(dp),     dimension(:), allocatable :: rcopy
    integer(i4b), dimension(:), optional    :: pixels
    integer(i4b), dimension(:), allocatable :: neighs
    integer(i4b), dimension(:), allocatable :: pix
    call free_pixinfo(pinfo)
    ! First find m. Since query_disc seems to be expensive, we don't
    ! want to do two passes just two find the number of pixels covered.
    ! Instead, we use maxrad together with the pixel size to get an upper bound.
    npix = 12*nside**2
    if(present(pixels)) then
       allocate(pix(size(pixels)))
       pix = pixels
    else
       allocate(pix(npix))
       do i = 1, npix
          pix(i) = i-1
       end do
    end if
    pinfo%pixside= sqrt(4*pi/npix)
    ! Make sure to include at least the 8 neighbors.
    effrad = min(pi,maxrad+pinfo%pixside*sqrt(2d0))
    m = 10+int(npix*effrad**2/4)
    n = size(pix)
    pinfo%npix  = npix
    pinfo%order = order
    pinfo%nside = nside
    pinfo%n     = n
    pinfo%rmax  = maxrad
    pinfo%m     = m
    allocate(pinfo%nneigh(n), pinfo%pix(n), pinfo%neigh(m,n), pinfo%r(m,n), &
     & pinfo%map(0:npix-1), pinfo%dir(2,m,n), pinfo%nme(m,n), pinfo%vec(3,n),&
     & pinfo%polar(2,n))
    pinfo%pix = pix
    isnest = 1; if(order == 1) isnest = 0

    ! Cache some common information. vec will contain 3d unit vectors,
    ! while pol has the theta and phi coordinates. The latter are not
    ! really needed multiple times here, but are relatively cheap to get
    ! (compared to the rotation info, for example).
    do i = 1, n
       pinfo%vec(:,i) = pix2vec(nside, order, pix(i))
       call pix2ang(nside, order, pix(i), pinfo%polar(1,i), pinfo%polar(2,i))
    end do

    ! Set up the map2mask
    allocate(neighs(npix),rcopy(m))
    pinfo%map = 0
    do i = 1, n
       pinfo%map(pix(i)) = i
    end do

    ! Now get all the neighbors
    do i = 1, n
       center = pinfo%vec(:,i)
       call query_disc(nside, pinfo%vec(:,i), effrad, neighs, nneigh,nest=isnest)
       ! But some of these may be outside our area. These must be removed.
       pinfo%nneigh(i) = 0
       do j = 1, nneigh
          if(pinfo%map(neighs(j)) == 0) cycle
          pinfo%nneigh(i) = pinfo%nneigh(i) + 1
          pinfo%neigh(pinfo%nneigh(i),i) = pinfo%map(neighs(j))
       end do
    end do
    ! Build up distances and directions
    do i = 1, n
       ! Build up a local, north-pointing orthogonal tangent space
       center= pinfo%vec(:,i)
       north = [0d0,0d0,1d0]
       east  = crossprod(north, center)
       ! Project down to tangent space
       basis(:,2) = east / sqrt(sum(east**2))
       basis(:,1) = north - dot_product(north,center)*center
       basis(:,1) = basis(:,1) / sqrt(sum(basis(:,1)**2))

       ! Finally, get the distances and directions
       do j = 1, pinfo%nneigh(i)
          pos = pinfo%vec(:,pinfo%neigh(j,i))
          call angdist(center, pos, pinfo%r(j,i))
          diff = pos - center*dot_product(pos,center)
          tmp  = sqrt(sum(diff**2))
          if(tmp > 0) diff = diff / tmp
          do k = 1, 2
             pinfo%dir(k,j,i) = dot_product(diff, basis(:,k))
          end do
       end do
    end do

    ! We need to be able to go from (j,pixi) to (i,pixj)
    do pixi = 1, n
       do j = 1, pinfo%nneigh(pixi)
          pixj = pinfo%neigh(j,pixi)
          do k = 1, pinfo%nneigh(pixj)
             if(pinfo%neigh(k,pixj) == pixi) exit
          end do
          pinfo%nme(j,pixi) = k
       end do
    end do

    ! We also need to know how (Q,U) change when we parallel transport. This
    ! is probably the heaviest step, and is only needed for polarization.
    ! We start at the current pixel, euler to the north pole, then euler to
    ! the other pixel. The transpose of this will be the rotation we need.
    ! We can then extract the psi angle of this rotation. The result is
    ! rot(:,:,j,i) being the rotation from j to i.
    allocate(pinfo%rot(3,3,m,n))
    do i = 1, n
       vi = pinfo%vec(:,i)
       do j = 1, pinfo%nneigh(i)
          vj = pinfo%vec(:,pinfo%neigh(j,i))
          pinfo%rot(:,:,j,i) = transport_rotmat(vj, vi)
       end do
    end do

    deallocate(neighs, rcopy, pix)
  end subroutine

  subroutine free_pixinfo(pinfo)
    implicit none
    type(pixinfo) :: pinfo
    if(allocated(pinfo%nneigh)) deallocate(pinfo%nneigh)
    if(allocated(pinfo%pix))    deallocate(pinfo%pix)
    if(allocated(pinfo%neigh))  deallocate(pinfo%neigh)
    if(allocated(pinfo%r))      deallocate(pinfo%r)
    if(allocated(pinfo%rot))    deallocate(pinfo%rot)
    if(allocated(pinfo%map))    deallocate(pinfo%map)
    if(allocated(pinfo%dir))    deallocate(pinfo%dir)
    if(allocated(pinfo%vec))    deallocate(pinfo%vec)
    if(allocated(pinfo%polar))  deallocate(pinfo%polar)
    if(allocated(pinfo%nme))    deallocate(pinfo%nme)
  end subroutine

  ! Produce the radial profile of a beam. The beam that must be input is
  ! the beam in alms, not the beam power, which is the square of the latter.
  subroutine beam_to_radial(beam, r, f)
    implicit none
    real(dp)      :: beam(0:), r(:), f(:), plm(0:ubound(beam,1))
    integer(i4b)  :: i, j, k, m, n, l
    do i = 1, size(r)
       call comp_normalised_Plm(ubound(beam,1), 0, r(i), plm)
       f(i) = 0
       do l = lbound(beam,1), ubound(beam,1)
          f(i) = f(i) + beam(l)*plm(l)/sqrt(4*pi/(2*l+1))
       end do
    end do
  end subroutine

  function radial_to_rmax(r, f) result(rmax)
    implicit none
    real(dp)     :: r(:), f(:), rmax, fmax, lim
    integer(i4b) :: i
    fmax = maxval(abs(f))
    lim  = beamtol*fmax
    do i = size(r), 2, -1
       if(abs(f(i)) > lim) exit
    end do
    rmax = r(i)
  end function

  subroutine alloc_weights(pinfo, weights)
    implicit none
    type(pixinfo)                         :: pinfo
    real(dp), dimension(:,:), allocatable :: weights
    if(allocated(weights)) deallocate(weights)
    allocate(weights(pinfo%m,pinfo%n))
  end subroutine

  ! Set up the weighs based on a radial profile f(r).
  ! Uses a spline for interpolation. This could break
  ! down for undersampled beams.
  subroutine calc_weights_radial(pinfo, r, f, weights, normalize)
    implicit none
    type(pixinfo) :: pinfo
    real(dp)      :: f(:), r(:), weights(:,:)
    integer(i4b)  :: i, j, k, m, n
    logical(lgt), optional :: normalize
    real(dp), dimension(:), allocatable :: f2
    allocate(f2(size(r)))
    call spline(r, f, 1d30, 1d30, f2)
    do i = 1, pinfo%n
       do j = 1, pinfo%nneigh(i)
          if(pinfo%r(j,i) > r(size(r))) then
             weights(j,i) = 0
          else
             weights(j,i) = splint(r, f, f2, pinfo%r(j,i))
          end if
       end do
    end do
    if(present(normalize)) then; if(normalize) then
       do i = 1, pinfo%n
          n = pinfo%nneigh(i)
          weights(1:n,i) = weights(1:n,i) / sum(weights(1:n,i))
       end do
    end if; end if
    deallocate(f2)
  end subroutine

  subroutine smooth_map_temp(imap, pinfo, weights, omap)
    implicit none
    real(dp)       :: imap(:), weights(:,:), omap(:)
    type(pixinfo)  :: pinfo
    integer(i4b)   :: pixi, pixj, i, j
    omap = 0
    do pixi = 1, pinfo%n
       do j = 1, pinfo%nneigh(pixi)
          pixj = pinfo%neigh(j,pixi)
          omap(pixi) = omap(pixi) + weights(j,pixi)*imap(pixj)
       end do
    end do
  end subroutine

  subroutine smooth_map_pol(imap, pinfo, weights, omap)
    implicit none
    real(dp)       :: imap(:,:), weights(:,:), omap(:,:)
    type(pixinfo)  :: pinfo
    integer(i4b)   :: pixi, pixj, i, j, a, b
    a = 1; if(size(imap,2) == 2) a = 2
    b = a+size(imap,2)-1
    omap = 0
    do pixi = 1, pinfo%n
       do j = 1, pinfo%nneigh(pixi)
          pixj = pinfo%neigh(j,pixi)
          omap(pixi,:) = omap(pixi,:) + weights(j,pixi)*matmul(pinfo%rot(a:b,a:b,j,pixi),imap(pixj,:))
       end do
    end do
  end subroutine

  ! Given a pixel-to-neighbor mapping 'neighs' and a compatible
  ! weights array, computes the covariance matrix based on
  ! a variance vector. Scales as O(nm^2), where n is the
  ! number of pixels, and m is the max number of neighbors.
  subroutine smooth_noise_diag_temp(var, pinfo, weights, covar)
    implicit none
    real(dp)       :: var(:), weights(:,:), covar(:,:)
    type(pixinfo)  :: pinfo
    call smooth_noise_block_temp(var, pinfo, weights, covar)
  end subroutine

  ! It is hard to make this general, so I restrict myself to 3 cases:
  ! T, Q/U and T/Q/U, which are recognized by the number of components.
  ! It uses the same insane component ordering we use other places, i.e.
  ! (pix,comp,pix,comp), and for the variance map: (pix,comp)
  subroutine smooth_noise_diag_pol(var, pinfo, weights, covar)
    implicit none
    real(dp)       :: var(:,:), weights(:,:), covar(:,:,:,:), tmp(size(covar,2),size(covar,2))
    type(pixinfo)  :: pinfo
    integer(i4b)   :: pixi, pixj, pixk, i, j, k, k2, a, b, c, ncomp, off
    ncomp = size(covar,2)
    off   = 0; if(ncomp == 2) off = 1
    covar = 0
    do pixi = 1, pinfo%n
       do k = 1, pinfo%nneigh(pixi)
          pixk = pinfo%neigh(k,pixi)
          do j = 1, pinfo%nneigh(pixk)
             pixj = pinfo%neigh(j,pixk)
             k2   = pinfo%nme(j,pixk)
             tmp = 0
             ! I don't know why j and i had to be swapped in pinfo%rot here,
             ! but that was necessary to get the correct answer. :/
             do a = 1, ncomp
                do b = 1, ncomp
                   do c = 1, ncomp
                      tmp(a,b) = tmp(a,b) + pinfo%rot(b+off,c+off,k,pixi)*var(pixk,c)* &
                       & pinfo%rot(a+off,c+off,k2,pixj)
                   end do
                end do
             end do
             covar(pixj,:,pixi,:) = covar(pixj,:,pixi,:) + tmp*weights(k,pixi)*weights(k2,pixj)
          end do
       end do
    end do
  end subroutine

  subroutine smooth_noise_block_temp(var, pinfo, weights, covar)
    implicit none
    real(dp)       :: var(:), weights(:,:), covar(:,:)
    type(pixinfo)  :: pinfo
    integer(i4b)   :: pixi, pixj, pixk, i, j, k, k2
    covar = 0
    do pixi = 1, pinfo%n
       do k = 1, pinfo%nneigh(pixi)
          pixk = pinfo%neigh(k,pixi)
          do j = 1, pinfo%nneigh(pixk)
             pixj = pinfo%neigh(j,pixk)
             k2   = pinfo%nme(j,pixk)
             covar(pixj,pixi) = covar(pixj,pixi) + weights(k,pixi)*weights(k2,pixj)*var(pixk)
          end do
       end do
    end do
  end subroutine

  !  ocov(i,j) = W(i,k)*icov(k,l)*W(j,l)'
  ! But icov is diagonal (when considering pixels), so
  !  ocov(i,j) = W(i,k)*icov(k,k)*W(j,k)'
  subroutine smooth_noise_block_pol(imap, pinfo, weights, ocov)
    implicit none
    real(dp)      :: imap(:,:,:), ocov(:,:,:,:), weights(:,:)
    real(dp)      :: tmp(size(imap,2),size(imap,2))
    type(pixinfo) :: pinfo
    integer(i4b)  :: pixi, pixj, pixk, pixl, i, j, k, k2, l, a, b
    real(dp), dimension(:), allocatable :: tcov(:,:,:,:)
    a = 1; if(size(imap,2) == 2) a = 2
    b = a+size(imap,2)-1
    ocov = 0
    do pixi = 1, pinfo%n
       do k = 1, pinfo%nneigh(pixi)
          pixk = pinfo%neigh(k,pixi)
          do j = 1, pinfo%nneigh(pixk)
             pixj = pinfo%neigh(j,pixk)
             k2   = pinfo%nme(j,pixk)
             ! I don't know why j and i had to be swapped in pinfo%rot here,
             ! but that was necessary to get the correct answer. :/
             ocov(pixj,:,pixi,:) = ocov(pixj,:,pixi,:) + &
              & weights(k,pixi)*matmul(pinfo%rot(a:b,a:b,k2,pixj),&
              & matmul(imap(pixk,:,:),transpose(pinfo%rot(a:b,a:b,k,pixi)))*&
              & weights(k2,pixj))
          end do
       end do
    end do
  end subroutine

  ! As above, but for the case of a full input covariance matrix.
  ! Scales as O(n^2m). Sadly, it is not possible to implement
  ! this by just forwarding to the more general function below.
  subroutine smooth_noise_full_temp(icov, pinfo, weights, ocov)
    implicit none
    real(dp)      :: icov(:,:), ocov(:,:), weights(:,:)
    type(pixinfo) :: pinfo
    integer(i4b)  :: pixi, pixj, pixk, pixl, i, j, k, l
    real(dp), dimension(:), allocatable :: tcov(:,:)

    allocate(tcov(size(icov,1),size(icov,2)))
    ! tcov = W * icov: tcov(i,j) = W(i,k)*icov(k,j), where W = weights
    tcov = 0
    do pixj = 1, pinfo%n
       do pixi = 1, pinfo%n
          do k = 1, pinfo%nneigh(pixi)
             pixk = pinfo%neigh(k,pixi)
             tcov(pixi,pixj) = tcov(pixi,pixj) + weights(k,pixi)*icov(pixk,pixj)
          end do
       end do
    end do
    ! Then do ocov = tcov * W: ocov(i,j) = tcov(i,k) * W(k,j)
    ocov = 0
    do pixj = 1, pinfo%n
       do pixi = 1, pinfo%n
          do k = 1, pinfo%nneigh(pixj)
             pixk = pinfo%neigh(k,pixj)
             ocov(pixi,pixj) = ocov(pixi,pixj) + tcov(pixi,pixk)*weights(k,pixj)
          end do
       end do
    end do
    deallocate(tcov)
  end subroutine

  ! Conceptually, this is the matrix product ocov = W * icov * W^T.
  ! It is therefore much faster to do it in two steps, though this costs some memory.
  ! For a possible scala-version, it might actually be best to have weights be non-sparse.
  !  ocov(i,j) = W(i,k)*R(i,k)*icov(k.l)*R(j,l)'*W(j,l)'
  ! W is *not symmetric*, but R is.
  ! We split it into two steps:
  !  tcov(i,l) = W(i,k)*R(i,k)*icov(k,l)
  !  ocov(i,j) = tcov(i,l)*R(l,j)*W(j,l)'
  subroutine smooth_noise_full_pol(icov, pinfo, weights, ocov)
    implicit none
    real(dp)      :: icov(:,:,:,:), ocov(:,:,:,:), weights(:,:)
    type(pixinfo) :: pinfo
    integer(i4b)  :: pixi, pixj, pixk, pixl, i, j, k, l, a, b
    real(dp), dimension(:), allocatable :: tcov(:,:,:,:)
    a = 1; if(size(icov,2) == 2) a = 2
    b = a+size(icov,2)-1
    allocate(tcov(size(icov,1),size(icov,2),size(icov,3),size(icov,4)))
    ! tcov = W * icov: tcov(i,j) = W(i,k)*icov(k,j), where W = weights*rot
    tcov = 0
    do pixj = 1, pinfo%n
       do pixi = 1, pinfo%n
          do k = 1, pinfo%nneigh(pixi)
             pixk = pinfo%neigh(k,pixi)
             tcov(pixi,:,pixj,:) = tcov(pixi,:,pixj,:) + weights(k,pixi)* &
              & matmul(pinfo%rot(a:b,a:b,k,pixi),icov(pixk,:,pixj,:))
          end do
       end do
    end do
    ! Then do ocov = tcov * W': ocov(i,j) = tcov(i,k) * W(j,k)'
    ocov = 0
    do pixj = 1, pinfo%n
       do pixi = 1, pinfo%n
          do k = 1, pinfo%nneigh(pixj)
             pixk = pinfo%neigh(k,pixj)
             ocov(pixi,:,pixj,:) = ocov(pixi,:,pixj,:) + matmul(tcov(pixi,:,pixk,:), &
              & transpose(pinfo%rot(a:b,a:b,k,pixj))) * weights(k,pixj)
          end do
       end do
    end do
    deallocate(tcov)
  end subroutine

  subroutine degrade_map_temp(imap, deg, omap)
    implicit none
    real(dp), dimension(:) :: imap, omap
    type(degrade_info)     :: deg
    integer(i4b)           :: i, j
    do i = 1, size(deg%opix)
       omap(i) = sum(imap(deg%submap(1:deg%nsub(i),i)))/deg%nsub(i)
    end do
  end subroutine

  subroutine degrade_map_pol(imap, deg, omap)
    implicit none
    real(dp), dimension(:,:) :: imap, omap
    type(degrade_info)       :: deg
    integer(i4b)             :: i, j, a, b
    a = 1; if(size(imap,2) == 2) a = 2
    b = a+size(imap,2)-1
    omap = 0
    do i = 1, deg%n
       do j = 1, deg%nsub(i)
          omap(i,:) = omap(i,:) + matmul(deg%rot(a:b,a:b,j,i),imap(deg%submap(j,i),:))
       end do
       omap(i,:) = omap(i,:) / deg%nsub(i)
    end do
  end subroutine

  subroutine degrade_noise_diag_temp(imap, deg, omap)
    implicit none
    real(dp), dimension(:)   :: imap, omap
    type(degrade_info)       :: deg
    integer(i4b)             :: i
    call degrade_map_temp(imap, deg, omap)
    omap = omap / deg%nsub
  end subroutine

  ! Bah, let us just make a copy and call the block version
  subroutine degrade_noise_diag_pol(imap, deg, omap)
    implicit none
    real(dp), dimension(:,:)   :: imap
    real(dp), dimension(:,:,:) :: omap
    type(degrade_info)         :: deg
    integer(i4b)               :: i
    real(dp), dimension(:,:,:), allocatable :: tmap
    allocate(tmap(size(imap,1),size(imap,2),size(imap,2)))
    tmap = 0
    do i = 1, size(imap,2); tmap(:,i,i) = imap(:,i); end do
    call degrade_noise_block_pol(tmap, deg, omap)
    deallocate(tmap)
  end subroutine

  subroutine degrade_noise_block_temp(imap, deg, omap)
    implicit none
    real(dp), dimension(:)   :: imap, omap
    type(degrade_info)       :: deg
    integer(i4b)             :: i
    call degrade_map_temp(imap, deg, omap)
    omap = omap / deg%nsub
  end subroutine

  ! A small covmatrix per pixel. When degrading, pixels stay independent,
  ! so output format is stil a map like this. In the case of smoothing, the
  ! output would have been a full covmat. The block format is (n,ncomp,ncomp) to
  ! make it healpix-compatible
  subroutine degrade_noise_block_pol(imap, deg, omap)
    implicit none
    real(dp), dimension(:,:,:)   :: imap, omap
    type(degrade_info)           :: deg
    integer(i4b)                 :: pixi, pixk, k, a, b
    a = 1; if(size(imap,2) == 2) a = 2
    b = a+size(imap,2)-1
    omap = 0
    do pixi = 1, size(omap,1)
       do k = 1, deg%nsub(pixi)
          pixk = deg%submap(k,pixi)
          omap(pixi,:,:) = omap(pixi,:,:) + matmul(deg%rot(a:b,a:b,k,pixi), &
           & matmul(imap(pixk,:,:), transpose(deg%rot(a:b,a:b,k,pixi))))
       end do
       omap(pixi,:,:) = omap(pixi,:,:) / deg%nsub(pixi)**2
    end do
  end subroutine

  subroutine degrade_noise_full_temp(icov, deg, ocov)
    implicit none
    real(dp)           :: icov(:,:), ocov(:,:)
    type(degrade_info) :: deg
    integer(i4b)       :: pixi, pixj, pixk, pixl, i, j, k, l
    real(dp), dimension(:), allocatable :: tcov(:,:)
    allocate(tcov(size(ocov,1),size(icov,2)))
    ! tcov = R*icov
    do pixj = 1, size(tcov,2)
       do pixi = 1, size(tcov,1)
          tcov(pixi,pixj) = sum(icov(deg%submap(1:deg%nsub(pixi),pixi),pixj)) / deg%nsub(pixi)
       end do
    end do
    ! With this, ocov = tcov*R^T
    do pixj = 1, size(ocov,2)
       do pixi = 1, size(ocov,1)
          ocov(pixi,pixj) = sum(tcov(pixi,deg%submap(1:deg%nsub(pixj),pixj))) / deg%nsub(pixj)
       end do
   end do
   deallocate(tcov)
  end subroutine

  ! Almost the same as the smoothing case. The only difference is that
  ! the output and input pixels are different.
  subroutine degrade_noise_full_pol(icov, deg, ocov)
    implicit none
    real(dp)           :: icov(:,:,:,:), ocov(:,:,:,:)
    type(degrade_info) :: deg
    integer(i4b)       :: pixi, pixj, pixk, pixl, i, j, k, l, a, b
    real(dp), dimension(:), allocatable :: tcov(:,:,:,:)
    a = 1; if(size(icov,2) == 2) a = 2
    b = a+size(icov,2)-1
    allocate(tcov(size(ocov,1),size(ocov,2),size(icov,3),size(icov,4)))
    ! tcov = R*icov
    tcov = 0
    do pixj = 1, size(tcov,3)
       do pixi = 1, size(tcov,1)
          do k = 1, deg%nsub(pixi)
             pixk = deg%submap(k,pixi)
             tcov(pixi,:,pixj,:) = tcov(pixi,:,pixj,:) + matmul(deg%rot(a:b,a:b,k,pixi), icov(pixk,:,pixj,:))
          end do
          tcov(pixi,:,pixj,:) = tcov(pixi,:,pixj,:) / deg%nsub(pixi)
       end do
    end do
    ! With this, ocov = tcov*R^T
    ocov = 0
    do pixj = 1, size(ocov,3)
       do pixi = 1, size(ocov,1)
          do k = 1, deg%nsub(pixj)
             pixk = deg%submap(k,pixj)
             ocov(pixi,:,pixj,:) = ocov(pixi,:,pixj,:) + matmul(tcov(pixi,:,pixk,:), &
              & transpose(deg%rot(a:b,a:b,k,pixj)))
          end do
          ocov(pixi,:,pixj,:) = ocov(pixi,:,pixj,:) / deg%nsub(pixj)
       end do
   end do
   deallocate(tcov)
  end subroutine

  subroutine prepare_degrade(deg, inpix, order, inside, onside)
    implicit none
    integer(i4b)       :: inside, onside, inpix(:), order, i, j, k, m, n
    type(degrade_info) :: deg
    integer(i4b), dimension(:), allocatable :: map2mask
    real(dp)           :: v1(3), v2(3)
    call free_degrade_info(deg)
    call assert(onside <= inside, "prepare_degrade does nots upport upgrading!")
    call assert(order == NEST, "Ring ordering not supported!")
    m  = (inside/onside)**2
    ! How many reduced pixels do we have?
    allocate(map2mask(0:12*onside**2-1))
    map2mask = 0
    do i = 1, size(inpix)
       j = inpix(i)/m
       map2mask(j) = map2mask(j) + 1
    end do
    n = count(map2mask > 0)
    ! Set up output pixel ordering
    k = 1
    do i = 0, size(map2mask)-1
       if(map2mask(i) == 0) cycle
       map2mask(i) = k
       k = k+1
    end do
    ! We are now ready to fill in the structures in deg.
    allocate(deg%submap(m,n), deg%ipix(size(inpix)), deg%opix(n), deg%nsub(n))
    deg%ipix   = inpix
    deg%nsub   = 0
    deg%submap = 0
    do i = 1, size(inpix)
       j = inpix(i)/m                ! Reduced pixel number
       k = map2mask(j)               ! Reduced pixel index
       deg%opix(k) = j
       deg%nsub(k) = deg%nsub(k) + 1 ! Number of subpixels we contain
       deg%submap(deg%nsub(k),k) = i ! New index to old index
    end do
    ! Sadly, we also need parallel transport info
    allocate(deg%rot(3,3,m,n))
    do i = 1, n
       v1 = pix2vec(onside, order, deg%opix(i))
       do j = 1, deg%nsub(i)
         v2 = pix2vec(inside, order, deg%ipix(deg%submap(j,i)))
         deg%rot(:,:,j,i) = transport_rotmat(v2,v1)
       end do
    end do

    deg%n = n
    deg%m = m
    deallocate(map2mask)
  end subroutine

  subroutine free_degrade_info(deg)
    implicit none
    type(degrade_info) :: deg
    if(allocated(deg%ipix))   deallocate(deg%ipix)
    if(allocated(deg%opix))   deallocate(deg%opix)
    if(allocated(deg%nsub))   deallocate(deg%nsub)
    if(allocated(deg%submap)) deallocate(deg%submap)
    if(allocated(deg%rot))    deallocate(deg%rot)
  end subroutine



  ! Matrix for transporting (t,q,u) from v1 to v2
  function transport_rotmat(v1, v2) result(mat)
    implicit none
    real(dp) :: v1(3), v2(3), mat(3,3), c2p, s2p
    call qu_transport_rot(v1, v2, c2p, s2p)
    mat = reshape([&
     & 1d0, 0d0, 0d0, &
     & 0d0, c2p,-s2p, &
     & 0d0, s2p, c2p],[3,3])
  end function

  subroutine read_beam_real(fname, beam)
    implicit none
    character(len=*)                      :: fname
    real(dp), dimension(:,:), allocatable :: beam
    integer(i4b)                          :: unit, i, n
    real(dp)                              :: tmp(2)
    ! First find length
    unit = getlun()
    open(unit,file=fname)
    n = 0
    do
       read(unit,*,end=1)
       n = n+1
    end do
1   allocate(beam(n,2))
    rewind(unit)
    do i = 1, n
       read(unit,*) beam(i,:)
    end do
    close(unit)
  end subroutine

  ! Read and write real-space beams, which are what this
  ! program cares about.
  subroutine write_beam_real(fname, beam)
    implicit none
    character(len=*) :: fname
    real(dp)         :: beam(:,:)
    integer(i4b)     :: unit, i
    unit = getlun()
    open(unit,file=fname)
    do i = 1, size(beam,1)
       write(unit,fmt="(2e20.10)") beam(i,:)
    end do
  end subroutine

  ! Pixel-space derivatives.
  ! T: 1-1: dx, dy, |d|, d^2
  ! P: 2-2: dx, dy, |d|, d^2; 2-1: grad, curl

  ! How to estimate the derivative in a point:
  ! Given a matrix of tangent space offsets
  ! v(n,2) and values relative to this point
  ! y(n) for n points near enough that the
  ! derivative dominates, a ML estimate for
  ! the derivative deriv(2) is:
  ! (v*deriv-y)^T*(v*deriv - y) = 0
  ! => deriv = (v^T*v)^-1(v^T*y)
  !
  ! To avoid unnecessary smoothing, the points
  ! should be from as small as possible area.
  ! For the first derivative, we ideally only
  ! need the ULDR points. Actually, with 4
  ! points, we should be able to estimate
  ! four parameters, almost enough to go to
  ! second order. With all 8 (7) neighbors,
  ! we can estimate dx, dy, dxx, dyy, dxy.
  ! But must calculate formula for this first.

end module
