!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
 
subroutine pos(ijk,ia,ja,ka,r0,r1,ds,rr)

! computes new position (r0 --> r1) of trajectory after time ds
! the new coordinate is still on one of the faces of box at ia,ja,ka
!
!  Input:
!
!    ijk      : considered direction (i=zonal, 2=meridional, 3=vertical)
!    ia,ja,ka : original position in integers
!    r0       : original non-dimensional position in the ijk-directionof particle 
!                (fractions of a grid box side in the corresponding direction)
!    rr       : time interpolation constant between 0 and 1 
!    sp       : crossing time to reach the grid box wall (in units of s/m3)
!
!  Output:
!    
!    r1       : the new position (coordinate)

USE mod_param
USE mod_vel
USE mod_turb
IMPLICIT none

real*8 r0,r1,rr,rg,ds,uu,um,vv,vm,en
integer ijk,ia,ja,ka,ii,im

rg=1.d0-rr

#ifdef mod2 
stop 2567 ! Kolla på gammal traj.F för OCCAM
#endif

if(ijk.eq.1) then
 ii=ia
 im=ia-1
 if(ii.eq.0) im=IMT
 uu=(rg*u(ia,ja,ka,NST)+rr*u(ia,ja,ka,1))*ff
 um=(rg*u(im,ja,ka,NST)+rr*u(im,ja,ka,1))*ff
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
 uu=(rg*v(ia,ja  ,ka,NST)+rr*v(ia,ja  ,ka,1))*ff
 um=(rg*v(ia,ja-1,ka,NST)+rr*v(ia,ja-1,ka,1))*ff
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
 uu=w(ka  )
 um=w(ka-1)
#ifdef turb    
 if(r0.ne.dble(ka  )) then
  uu=uu+upr(5,2)  
 else
  uu=uu+upr(5,1)  ! add u' from previous iterative time step if on box wall
 endif
#endif
endif

!
! note: consider in future to improve the code below for accuracy 
! in case of um-uu = small; also see subroutine cross
!
if(um.ne.uu) then
 r1= (r0+(-dble(ii-1) + um/(uu-um))) * dexp( (uu-um)*ds ) + dble(ii-1) - um/(uu-um)
else
 r1=r0+uu*ds
endif

return
end subroutine pos
