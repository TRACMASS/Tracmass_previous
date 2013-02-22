subroutine pos_orgn(ijk,ia,ja,ka,r0,r1,ds,rr)
#ifndef timeanalyt 
  !====================================================================
  ! computes new position (r0 --> r1) of trajectory after time ds
  ! the new coordinate is still on one of the faces of box at ia,ja,ka
  !
  !  Input:
  !
  !    ijk      : considered direction (i=zonal,2=meridional,3=vertical)
  !    ia,ja,ka : original position in integers
  !    r0       : original non-dimensional position in the ijk-direction
  !                   of particle (fractions of a grid box side in the 
  !                                corresponding direction)
  !    rr       : time interpolation constant between 0 and 1 
  !    sp       : crossing time to reach the grid box wall (units=s/m3)
  !
  !  Output:
  !    
  !    r1       : the new position (coordinate)
  !====================================================================
  
  USE mod_param
  USE mod_vel
  USE mod_turb
  IMPLICIT none

  real*8 r0,r1,rr,rg,ds,uu,um,vv,vm,en
  integer ijk,ia,ja,ka,ii,im
  
  rg=1.d0-rr
  
#ifdef twodim   
  if(ijk.eq.3) then
     r1=r0
     return
  endif
#endif
  
  
  if(ijk.eq.1) then
     ii=ia
     im=ia-1
     if(im.eq.0) im=IMT
     uu=(rg*uflux(ia,ja,ka,NST)+rr*uflux(ia,ja,ka,1))*ff
     um=(rg*uflux(im,ja,ka,NST)+rr*uflux(im,ja,ka,1))*ff
#ifdef turb    
     if(r0.ne.dble(ii)) then
        uu=uu+upr(1,2)  
     else
        uu=uu+upr(1,1)  
        ! add u' from previous iterative time step if on box wall
     endif
     if(r0.ne.dble(im)) then
        um=um+upr(2,2)
     else
        um=um+upr(2,1)  
        ! add u' from previous iterative time step if on box wall
     endif
#endif
  elseif(ijk.eq.2) then
     ii=ja
     uu=(rg*vflux(ia,ja  ,ka,NST)+rr*vflux(ia,ja  ,ka,1))*ff
     um=(rg*vflux(ia,ja-1,ka,NST)+rr*vflux(ia,ja-1,ka,1))*ff
#ifdef turb    
     if(r0.ne.dble(ja  )) then
        uu=uu+upr(3,2)  
     else
        uu=uu+upr(3,1)  
        ! add u' from previous iterative time step if on box wall
     endif
     if(r0.ne.dble(ja-1)) then
        um=um+upr(4,2)
     else
        um=um+upr(4,1)  
        ! add u' from previous iterative time step if on box wall
     endif
#endif
  elseif(ijk.eq.3) then
     ii=ka
#if defined full_wflux | defined explicit_w
     uu=wflux(ia ,ja ,ka   ,1)
     um=wflux(ia ,ja ,ka-1 ,1)
#else
     uu=rg*wflux(ka  ,NST)+rr*wflux(ka  ,1)
     um=rg*wflux(ka-1,NST)+rr*wflux(ka-1,1)
#endif
#ifndef twodim   
#ifdef turb    
     if(r0.ne.dble(ka  )) then
        uu=uu+upr(5,2)  
     else
        uu=uu+upr(5,1)  
        ! add u' from previous iterative time step if on box wall
     endif
     if(r0.ne.dble(ka-1)) then
        uu=uu+upr(6,2)  
     else
        uu=uu+upr(6,1)  
        ! add u' from previous iterative time step if on box wall
     endif
#endif
#endif
  endif
  
  !
  ! note: consider in future to improve the code below for accuracy 
  ! in case of um-uu = small; also see subroutine cross
  if(um.ne.uu) then
     r1= (r0+(-dble(ii-1) + um/(uu-um))) * dexp( (uu-um)*ds ) + dble(ii-1) - um/(uu-um)
  else
     r1=r0+uu*ds
  endif
  !if(abs(um/(uu-um)).gt.1.d10) print *,'possible precision problem?',um/(uu-um),uu,um,ijk,ia,ja,ka,r0,r1,ds,rr
  
  return
#endif
end subroutine pos_orgn

