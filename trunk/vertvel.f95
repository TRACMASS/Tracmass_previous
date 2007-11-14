!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

subroutine vertvel(rr,ia,iam,ja,ka)

! Computes the vertical velocity by integrating the continuity eq. from the bottom 


USE mod_param
USE mod_vel
USE mod_turb
#if defined sediment
USE mod_sediment
USE mod_orbital
USE mod_grid
#endif

IMPLICIT none

#if defined sediment
REAL*8 wsedtemp
REAL kin
#endif

real*8 rr,rg,uu,um,vv,vm
integer ia,iam,ja,ka,k

rg=1.d0-rr

w(0)=0.d0
do k=1,ka
#ifdef twodim
 w(k)=0.d0
#else 
 uu=rg*u(ia ,ja  ,k,NST)+rr*u(ia ,ja  ,k,1)
 um=rg*u(iam,ja  ,k,NST)+rr*u(iam,ja  ,k,1)
 vv=rg*v(ia ,ja  ,k,NST)+rr*v(ia ,ja  ,k,1)
 vm=rg*v(ia ,ja-1,k,NST)+rr*v(ia ,ja-1,k,1)
 w(k) = w(k-1) - ff * ( uu - um + vv - vm )
#endif
enddo


#ifdef turb 
#ifndef twodim   
! Calculates the w' at the top of the box from the divergence of u' and v' in order
! to respect the continuity equation even for the turbulent velocities
 do k=1,2
  uu=upr(1,k)  
  um=upr(2,k)
  vv=upr(3,k)
  vm=upr(4,k)
  upr(5,k) = w(ka-1) - ff * ( uu - um + vv - vm )
  upr(6,k) = 0.d0
 enddo
#endif
#endif


#ifdef sediment
! Godtyckligt värde på kinetiska energin där wsed inte längre påverkar, 3e6.

do k=0,km
   wsedtemp=0.
   kin=(u(ia,ja,k,1)*u(ia,ja,k,1)+v(ia,ja,k,1)*v(ia,ja,k,1))*0.5
   !if (kin.le.3000000) then   !för RCO
      !wsedtemp=wsed*(3000000-kin)/3000000
   if (kin.le.kincrit) then   !för SKB
      wsedtemp=wsed*(kincrit-kin)/kincrit
   endif
   w(k)=w(k) +  wsedtemp * dxdy(ia,ja)     ! *dx *dy *deg**2 *cst(jb)
enddo
#endif

return
end subroutine vertvel

