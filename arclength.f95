! Compu%tes the length of the trajectory
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
  ! === u ===
  ii=ia
  iim=ia-1
  if(iim.eq.0) iim=imt
  uu=0.5*ff*( (rg*uflux(ii ,ja,ka,NST)+rr*uflux(ii ,ja,ka,1)) &
       + (rg*uflux(iim,ja,ka,NST)+rr*uflux(iim,ja,ka,1)) ) / (dy*deg*dz(ka))
  
  ! === v ===
  jj=ja
  jjm=ja-1
  vv=0.5*ff*(  (rg*vflux(ia,jj ,ka,NST)+rr*vflux(ia,jj ,ka,1)) &
       + (rg*vflux(ia,jjm,ka,NST)+rr*vflux(ia,jjm,ka,1)) ) /(dx*deg*cst(ja)*dz(ka))
  
  ! === w ===
  kk=ka
  kkm=ka-1
#ifdef  full_wflux
  ww=0.5*( wflux(ia ,ja ,kk ,1) + wflux(ia ,ja ,kkm ,1) )/dxdy(ia,ja)
#else  
  ww=0.5*( wflux(kk)            + wflux(kkm)            )/dxdy(ia,ja)
#endif
  ! arclength
  arc=dsqrt(uu**2 + vv**2 + ww**2 ) * dt
  return
end subroutine arclength

