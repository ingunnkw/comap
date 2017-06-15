module comap_patch_detect_mod
  use quiet_utils
  use comap_patch_mod
  use comap_pointing_mod
  use comap_detector_mod
!  use quiet_ephem_mod
  implicit none

  !interface get_patch_pos
  !   module procedure get_patch_pos_single, get_patch_pos_multi
  !end interface

  ! General procedure: Allocate space for patch positions:
  ! allocate(patch_pos(3,nsamp,npatch))
  ! call get_patch_pos(patches, t, coords, patch_pos)
  ! And then use get_patch_hits once for each pointing
  ! stream you have.
  ! For efficiency reasons, it can be a good idea to
  ! call get_patch_pos with t sparsely sampled, for
  ! example t(::100). Make sure that nsamp in patch
  ! pos matches this (nsamp=size(t(::100))).

contains

  subroutine initialize_patch_detect_mod(parfile)
    implicit none
    character(len=*)   :: parfile
    logical(lgt), save :: initialized = .false.
    if(initialized) return
    call initialize_comap_patch_mod(parfile)
    call initialize_comap_pointing_mod(parfile)
    call initialize_detector_mod(parfile)
    initialized = .true.
  end subroutine

  function get_patch_pos_single(pinfo, t, sys) result(pos)
    implicit none
    real(dp)               :: t, pos(2), psi
    type(patch_info)       :: pinfo
    integer(i4b)           :: sys
    pos = pinfo%pos
    call coord_convert(coord_gal,pos(1),pos(2),0d0,sys,pos(1),pos(2),psi)
!!$    if(pinfo%fixed) then
!!$       pos = pinfo%pos
!!$       call coord_convert(coord_gal,pos(1),pos(2),0d0,sys,pos(1),pos(2),psi)
!!$    else
!!$       pos = ephem(pinfo%eph, t)
!!$       call coord_convert(coord_cel,pos(1),pos(2),0d0,sys,pos(1),pos(2),psi)
!!$    end if
  end function
!!$
!!$  function get_patch_dist(pinfo, t) result(dist)
!!$    implicit none
!!$    real(dp)               :: t, dist
!!$    type(patch_info)       :: pinfo
!!$    if(pinfo%fixed) then
!!$       dist = infinity
!!$    else
!!$       dist = ephem_dist(pinfo%eph, t)
!!$    end if
!!$  end function

  subroutine get_patch_pos_multi(patches, t, sys, pos, angles)
    implicit none
    type(patch_info)       :: patches(:)
    real(dp)               :: t(:), pos(:,:,:)
    real(dp)               :: p(2)
    logical(lgt), optional :: angles
    logical(lgt)           :: ang
    integer(i4b)           :: i, j, m, np, sys
    ang = .false.; if(present(angles)) ang = angles
    m  = size(t)
    np = size(patches)
    do i = 1, m
       do j = 1, np
          p = get_patch_pos_single(patches(j),t(i),sys)
          if(ang) then
             pos(:,i,j) = [ p(1), p(2), 0d0 ]
          else
             call ang2vec(p(2),p(1),pos(:,i,j))
          end if
       end do
    end do
  end subroutine

  ! For a given set of pointings, find out which objects are hit, using
  ! precomputed patch location information from get_patch_pos_multi.
  ! This does not need to have the same sampling rate as phi and theta,
  ! the nearest point is used.
  subroutine get_patch_hits(patches, patch_pos, phi, theta, rad, objs, dists)
    implicit none
    type(patch_info)   :: patches(:)
    real(dp)           :: phi(:), theta(:), patch_pos(:,:,:)
    real(dp)           :: dists2(size(patches)), v(3), rad, bval
    real(dp), optional :: dists(:)
    integer(i4b)       :: objs(:)
    integer(i4b)       :: i, j, k, l, n, m, np, sys, nh, bloc, mod
    n  = size(phi,1)
    np = size(patches)
    if(np == 0) return
    do j = 1, n
       i = (j-1)*size(patch_pos,2)/n+1
       objs(j) = 0
       call ang2vec(theta(j),phi(j),v)
       do k = 1, np
          call angdist(patch_pos(:,i,k),v, dists2(k))
       end do
       bloc = maxloc(patches%priority, 1, dists2 < patches%obj_rad + rad)
       if(bloc == 0) cycle
       if(dists2(bloc) >= patches(bloc)%obj_rad+rad) cycle
       bval = patches(bloc)%priority
       if(bval == 0) cycle
       objs(j) = minloc(dists2, 1, patches%priority == bval)
       if(present(dists)) dists(j) = dists2(objs(j))
    end do
  end subroutine

  function reduce_hits(patches, hits, dists) result(best)
    implicit none
    type(patch_info)   :: patches(:)
    integer(i4b)       :: hits(:), best, i, j
    real(dp), optional :: dists(:)
    real(dp)           :: bd, bv
    bv = 0
    bd = pi
    best = 0
    do i = 1, size(hits)
       if(hits(i) == 0) cycle
       if(patches(hits(i))%priority > bv) then
          best = hits(i)
          bv   = patches(best)%priority
       elseif(present(dists) .and. patches(hits(i))%priority == bv) then
          if(dists(i) < bd) then
             best = hits(i)
             bd = dists(i)
          end if
       end if
    end do
  end function

end module
