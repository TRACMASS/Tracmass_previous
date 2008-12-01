
subroutine vertvel(rr,ia,iam,ja,ka)

#ifndef explicit_w
  
  ! === Computes the vertical velocity by integrating ===
  ! === the continuity eq. from the bottom            ===
  
  USE mod_param
  USE mod_vel
  USE mod_turb
#ifdef ifs
  USE mod_grid
#endif
#ifdef sediment
  USE mod_sed
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
  
  wflux=0.d0
#ifdef twodim
  return
#else
  do k=1,ka
     uu=rg*uflux(ia ,ja  ,k,NST)+rr*uflux(ia ,ja  ,k,1)
     um=rg*uflux(iam,ja  ,k,NST)+rr*uflux(iam,ja  ,k,1)
     vv=rg*vflux(ia ,ja  ,k,NST)+rr*vflux(ia ,ja  ,k,1)
     vm=rg*vflux(ia ,ja-1,k,NST)+rr*vflux(ia ,ja-1,k,1)
#if defined ifs
!     wflux(k) = wflux(k-1) + ff * ( uu - um + vv - vm )
     wflux(k) = wflux(k-1) + ff * &
     ( uu - um + vv - vm - (dztb(ia,ja,k,2)-dztb(ia,ja,k,1))/21600. )
#elif full_wflux
     wflux(ia,ja,k,1)=wflux(ia,ja,k-1,1) - ff * ( uu - um + vv - vm )
#else
     wflux(k) = wflux(k-1) + ff * ( uu - um + vv - vm )
#endif
  enddo
  
#ifdef sediment
  ! === Godtyckligt vaerde paa kinetiska energin ===
  ! === daer wsed inte laengre paaverkar, 3e6.   ===
  
  kloop: do k=0,km
     wsedtemp=0.d0
     kin=(uflux(ia,ja,k,1)*uflux(ia,ja,k,1)+ &
          vflux(ia,ja,k,1)*vflux(ia,ja,k,1))*0.5d0
     !if (kin.le.3000000) then   !för RCO
     !wsedtemp=wsed*(3000000-kin)/3000000
     if (kin.le.kincrit) then   !för SKB
        wsedtemp=wsed*(kincrit-kin)/kincrit
     endif
#ifdef full_wflux
     wflux(k)=wflux(ia,ja,k,1) +  wsedtemp * dxdy(ia,ja)     ! *dx *dy *deg**2 *cst(jb)
#else
     wflux(k)=wflux(k) +  wsedtemp * dxdy(ia,ja)     ! *dx *dy *deg**2 *cst(jb)
#endif
     
  end do kloop
#endif
  
#endif
  return
#endif
end subroutine vertvel

 
