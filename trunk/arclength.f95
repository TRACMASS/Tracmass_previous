!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
! Computes the length of the trajectory
! Note that the calculation is not exact 
 
subroutine arclength(ia,ja,ka,dt,rr,arc)

USE mod_param
USE mod_coord
USE mod_grid
USE mod_vel
IMPLICIT none

REAL*8 rr,rg,dt,arc,uu,vv,ww
INTEGER ii,iim,jj,jjm,kk,kkm,ia,ja,ka

rg=1.-rr

! velocities
! u
ii=ia
iim=ia-1
if(iim.eq.0) iim=imt
uu=0.5*ff*( (rg*u(ii ,ja,ka,NST)+rr*u(ii ,ja,ka,1)) &
          + (rg*u(iim,ja,ka,NST)+rr*u(iim,ja,ka,1)) ) / (dy*deg*dz(ka))
! v
jj=ja
jjm=ja-1
vv=0.5*ff*(  (rg*v(ia,jj ,ka,NST)+rr*v(ia,jj ,ka,1)) &
           + (rg*v(ia,jjm,ka,NST)+rr*v(ia,jjm,ka,1)) ) /(dx*deg*cst(ja)*dz(ka))
! w
kk=ka
kkm=ka-1
ww=0.5*( w(kk ) + w(kkm) )/dxdy(ia,ja)

! arclength
arc=dsqrt(uu**2 + vv**2 + ww**2 ) * dt

return
end subroutine arclength

