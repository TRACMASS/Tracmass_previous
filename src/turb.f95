#ifdef turb    
subroutine turbuflux(ia,ja,ka,rr)
  
  ! computes the paramterised turbulent velocities u' and v' into upr
  
  USE mod_param
  USE mod_vel
  USE mod_turb
  IMPLICIT none
  
  REAL*8         :: uv(12),rr,rg,en,localW
  REAL*4         :: qran(12)
  INTEGER        :: ia,ja,ka,im,jm,n

!  call random_number(qran) ! generates a random number between 0 and 1
  
! random generated numbers between 0 and 1
do n=1,12
 qran(n)=rand()
enddo
	
  qran=2.*qran-1. ! === Max. amplitude of turbulence with random numbers between -1 and 1
                  ! === (varies with the same aplitude as the mean vel)
  !rand=1.*rand-0.5  ! Reduced amplitude of turb.
  
  rg=1.d0-rr
  
  im=ia-1
  if(im.eq.0) im=IMT
  jm=ja-1
  
#if defined timestep
  ! time interpolated velocities
  uv(1)=(rg*uflux(ia,ja,ka,NST)+rr*uflux(ia,ja,ka,1))*ff ! western u
  uv(2)=(rg*uflux(im,ja,ka,NST)+rr*uflux(im,ja,ka,1))*ff ! eastern u
  uv(3)=(rg*vflux(ia,ja,ka,NST)+rr*vflux(ia,ja,ka,1))*ff ! northern v
  uv(4)=(rg*vflux(ia,jm,ka,NST)+rr*vflux(ia,jm,ka,1))*ff ! southern v
#elif defined timeanalyt
  uv( 1)=uflux(ia,ja,ka,1)*ff ! western u at t-1
  uv( 2)=uflux(im,ja,ka,1)*ff ! eastern u at t-1
  uv( 3)=vflux(ia,ja,ka,1)*ff ! northern v at t-1
  uv( 4)=vflux(ia,jm,ka,1)*ff ! southern v at t-1
  uv( 7)=uflux(ia,ja,ka,2)*ff ! western u at t
  uv( 8)=uflux(im,ja,ka,2)*ff ! eastern u at t
  uv( 9)=vflux(ia,ja,ka,2)*ff ! northern v at t
  uv(10)=vflux(ia,jm,ka,2)*ff ! southern v at t
#endif  

!   upr(:,1)=upr(:,2) ! store u' from previous time iteration step
   
! multiply the time interpolated velocities by random numbers
   
! 4 different random velocities at the 4 walls
!   upr(:,2)=uv(:)*rand(:) 
! or same random velocities at the eastern and western as well as nothern and southern
  upr(1,2)=uv(1)*qran(1)
  upr(2,2)=uv(2)*qran(1)
  upr(3,2)=uv(3)*qran(3)
  upr(4,2)=uv(4)*qran(3)

  upr(:,1)=upr(:,2) ! impose same velocities for t-1 and t (this is more stabel but why? K.D��s)
  
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
     !Detta g�r att man bara justerar vertikala hastigheten p� ovansidan av boxen ! u2
     !  upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     
     ! === The vertical velocity is calculated from 
     ! === the divergence of the horizontal turbulence !v2
     ! === The turbulent is evenly spread between the 
     ! === top and bottom of box if not a bottom box

#ifdef full_wflux
    localW=wflux(ia,ja,ka-1,1)
#else
    localW=wflux(ka-1,1)
#endif
    
    if(localW.eq.0.d0) then
       upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
       upr(6,n) = 0.d0
    else
       upr(5,n) = - 0.5d0 * ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
       upr(6,n) = - upr(5,n)
    endif

#if defined timeanalyt
    if(localW.eq.0.d0) then
       upr(11,n) = - ff * ( upr(7,n) - upr(8,n) + upr(9,n) - upr(10,n) )
       upr(12,n) = 0.d0
    else
       upr(11,n) = - 0.5d0 * ff * ( upr(7,n) - upr(8,n) + upr(9,n) - upr(10,n) )
       upr(12,n) = - upr(11,n)
    endif
#endif

 enddo
#endif
#endif
 
 
 return
end subroutine turbuflux
#endif
!_______________________________________________________________________