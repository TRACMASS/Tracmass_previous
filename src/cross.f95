#ifndef timeanalyt 

subroutine cross(ijk,ia,ja,ka,r0,sp,sn,rr)
  
  ! subroutine to compute time (sp,sn) when trajectory 
  ! crosses face of box (ia,ja,ka) 
  ! two crossings are considered for each direction:  
  ! east and west for longitudinal directions, etc.  
  ! 
  !  Input:
  !
  !  ijk      : considered direction (i=zonal, 2=meridional, 3=vertical)
  !  ia,ja,ka : original position in integers
  !  r0       : original non-dimensional position in the 
  !             ijk-direction of particle 
  !             (fractions of a grid box side in the 
  !              corresponding direction)
  !  rr       : time interpolation constant between 0 and 1 
  !
  !
  !  Output:
  !
  !    sp,sn   : crossing time to reach the grid box wall 
  !              (in units of s/m3)
  

USE mod_param
USE mod_vel
USE mod_turb
IMPLICIT none

real*8 r0,ba,sp,sn,uu,um,rr,rg,vv,vm
integer ijk,ia,ja,ka,ii,im

rg=1.d0-rr

#ifdef mod2 
uu=rg*u2(ii ,ja,ka,ns)+rr*u2(ii ,ja,ka,nsn)
um=rg*u2(iim,ja,ka,ns)+rr*u2(iim,ja,ka,nsn)
stop 2567 ! Kolla på gammal traj.F för OCCAM
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
  uu=uu+upr(1,1)  ! add u' from previous iterative time step if on box wall
 endif
 if(r0.ne.dble(im)) then
  um=um+upr(2,2)
 else
  um=um+upr(2,1)  ! add u' from previous iterative time step if on box wall
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
  uu=uu+upr(3,1)  ! add u' from previous iterative time step if on box wall
 endif
 if(r0.ne.dble(ja-1)) then
  um=um+upr(4,2)
 else
  um=um+upr(4,1)  ! add u' from previous iterative time step if on box wall
 endif
#endif
elseif(ijk.eq.3) then
 ii=ka
#ifdef full_wflux
 uu=wflux(ia ,ja ,ka   ,1)
 um=wflux(ia ,ja ,ka-1 ,1)
#elif defined timeanalyt || timestep
 uu=rg*wflux(ka  ,NST)+rr*wflux(ka  ,1)
 um=rg*wflux(ka-1,NST)+rr*wflux(ka-1,1)
#else
 uu=wflux(ka)
 um=wflux(ka-1)
#endif

#ifndef twodim   
#ifdef turb   
 if(r0.ne.dble(ka  )) then
  uu=uu+upr(5,2)  
 else
  uu=uu+upr(5,1)  ! add u' from previous iterative time step if on box wall
 endif
 if(r0.ne.dble(ka-1)) then
  uu=uu+upr(6,2)  
 else
  uu=uu+upr(6,1)  ! add u' from previous iterative time step if on box wall
 endif
#endif
#endif
endif

! east, north or upward crossing
if(uu.gt.0.d0 .and. r0.ne.dble(ii)) then
 if(um.ne.uu) then
  ba=(r0+dble(-ii+1)) * (uu-um) + um
  if(ba.gt.0.d0) then
  sp=-1.d0/(um-uu)*( dlog(uu) - dlog(ba) )
  else
   sp=1.d20
  endif
 else
  sp=(dble(ii)-r0)/uu
 endif
else
 sp=1.d20
endif

if(sp.le.0.d0) sp=1.d20

! west, south or downward crossing
if(um.lt.0.d0 .and. r0.ne.dble(ii-1)) then
 if(um.ne.uu)then
  ba=-((r0-dble(ii))*(uu-um)+uu) 
  if(ba.gt.0.d0) then
   sn=-1.d0/(um-uu)*( dlog(-um) - dlog(ba)  )
  else
   sn=1.d20
  endif
 else
  sn=(dble(ii-1)-r0)/uu
 endif
else
 sn=1.d20
endif

if(sn.le.0.d0) sn=1.d20

return
end subroutine cross
#endif
