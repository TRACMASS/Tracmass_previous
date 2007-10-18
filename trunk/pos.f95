!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
 
subroutine pos(ijk,ia,ja,ka,r0,r1,ds,rr)

! computes new coordinate (r0 --> r1) of trajectory after time ds
! the new coordinate is still on one of the faces of box at ia,ja,ka

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

!#ifdef turb
! if(iim.eq.0) iim=IMT
! uu=(rg*u(ia  ,ja,ka,NST)+rr*u(ia  ,ja,ka,1))*ff
! um=(rg*u(ia-1,ja,ka,NST)+rr*u(ia-1,ja,ka,1))*ff
! vv=(rg*v(ia,ja  ,ka,NST)+rr*v(ia,ja  ,ka,1))*ff
! vm=(rg*v(ia,ja-1,ka,NST)+rr*v(ia,ja-1,ka,1))*ff
! en=0.25*sqrt(uu**2+um**2+vv**2+vm**2)
!#endif

if(ijk.eq.1) then
 ii=ia
 im=ia-1
 if(ii.eq.0) im=IMT
 uu=(rg*u(ia,ja,ka,NST)+rr*u(ia,ja,ka,1))*ff
 um=(rg*u(im,ja,ka,NST)+rr*u(im,ja,ka,1))*ff
#ifdef turb    
 uu=uu+upr(1)  
 um=um+upr(2)
#endif
elseif(ijk.eq.2) then
 ii=ja
 uu=(rg*v(ia,ja  ,ka,NST)+rr*v(ia,ja  ,ka,1))*ff
 um=(rg*v(ia,ja-1,ka,NST)+rr*v(ia,ja-1,ka,1))*ff
#ifdef turb    
 uu=uu+upr(3)  
 um=um+upr(4)
#endif
elseif(ijk.eq.3) then
 ii=ka
 uu=w(ka  )
 um=w(ka-1)
#ifdef turb    
! uu=uu*rand(5)
! um=um*rand(6)
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
