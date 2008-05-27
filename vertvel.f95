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
#if defined ifs
   w(k) = w(k-1) - ff * ( uu - um + vv - vm )
#else
   w(k) = w(k-1) + ff * ( uu - um + vv - vm )
#endif
#endif
enddo

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

