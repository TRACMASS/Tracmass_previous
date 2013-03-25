! Binned Random Walk I from Thygesen & Aadlandsvik 2007

Subroutine bin_rand_walk

! h=timestep
! D=eddy diffusivity
! k=width of depth cell sigma
! sigma=depth cell
! sigma-1=depth cell below
! sigma+1=depth cell above

implicit none
Real*8      D,h,U,q,p, dz
Integer     k, ntrac

!loop over particles (ntrac)
!loop over time (dtmin) [Time loop is elsewhere.]

do ntrac=1,ntractot
  D=AKt
  h=dtmin
  U=ran(iseed)
  dz=s_w-(s_w-1)
  q = D(k-1)*2 / ((dz(k-1)+dz(k))/k)
  p = D(k)*2 / ((dz(k)+dz(k+1))/k)

! Hsbl = depth of oceanic surface boundary layer
! Tcline = S-coordinate surface/bottom layer width
! s_rho = S-coordinate at RHO-points
! s_w = S-coordinate at W-points


  If (0. .le. U) and (U .lt. q*h) then
     k(ntrac) = k(ntrac) - 1  !move to cell below
  Elseif (q*h .le U) and (U .lt. (1. - p*h)) then
     k(ntrac) = k(ntrac)      !stay in same cell
  Elseif ((1. - p*h .le. U)) and (U .lt. 1.) then
     k(ntrac) = k(ntrac) + 1  !move to cell above
  Endif
end do

End Subroutine bin_rand_walk
