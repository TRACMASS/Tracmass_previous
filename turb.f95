#ifdef turb    
subroutine turbuflux(ia,ja,ka,rr)
  
  ! computes the paramterised turbulent velocities u' and v' into upr
  
  USE mod_param
  USE mod_vel
  USE mod_turb
  IMPLICIT none
  
  real*8 uv(6),rr,rg,en,localW
  integer ia,ja,ka,im,jm,n
  
  call random_number(rand)
  rand=2.*rand-1. ! === Max. amplitude of turb. 
                  ! === (varies with the same aplitude as the mean vel)
  !rand=1.*rand-0.5  ! Reduced amplitude of turb.
  
  rg=1.d0-rr
  
  im=ia-1
  if(im.eq.0) im=IMT
  jm=ja-1
  
  ! time interpolated velocities
  uv(1)=(rg*uflux(ia,ja,ka,NST)+rr*uflux(ia,ja,ka,1))*ff ! western u
  uv(2)=(rg*uflux(im,ja,ka,NST)+rr*uflux(im,ja,ka,1))*ff ! eastern u
  uv(3)=(rg*vflux(ia,ja,ka,NST)+rr*vflux(ia,ja,ka,1))*ff ! northern v
  uv(4)=(rg*vflux(ia,jm,ka,NST)+rr*vflux(ia,jm,ka,1))*ff ! southern v
  
  !upr=uv*rand
  
  !upr(:,1)=upr(:,2) ! store u' from previous time iteration step
  upr(:,2)=uv(:)*rand(:)
  upr(:,1)=upr(:,2)
  
#ifdef turb 
#ifndef twodim   
  ! === Calculates the w' at the top of the box from 
  ! === the divergence of u' and v' in order
  ! === to respect the continuity equation even 
  ! === for the turbulent velocities
  do n=1,2
     
     !Detta ser ut som en bugg   ! t2
     !  upr(5,n) = w(ka-1) - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     !Detta gör att man bara justerar vertikala hastigheten på ovansidan av boxen ! u2
     !  upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     
     ! === The vertical velocity is calculated from 
     ! === the divergence of the horizontal turbulence !v2
     ! === The turbulent is evenly spread between the 
     ! === top and bottom of box if not a bottom box

#ifdef full_wflux
    localW=wflux(ia,ja,ka-1,1)
#else
    localW=wflux(ka-1)
#endif
    
    if(localW.eq.0.d0) then
       upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
       upr(6,n) = 0.d0
    else
       upr(5,n) = - 0.5d0 * ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
       upr(6,n) = - upr(5,n)
    endif
 enddo
#endif
#endif
 
 
 return
end subroutine turbuflux
#endif
!_______________________________________________________________________
