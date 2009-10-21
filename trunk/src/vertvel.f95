
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
  integer ia,iam,ja,ka,k,n
    
  rg=1.d0-rr
  wflux=0.d0
  
#ifdef twodim
  return
  
! start 3D code
#else
  kloop: do k=1,ka
     uu=rg*uflux(ia ,ja  ,k,NST)+rr*uflux(ia ,ja  ,k,1)
     um=rg*uflux(iam,ja  ,k,NST)+rr*uflux(iam,ja  ,k,1)
     vv=rg*vflux(ia ,ja  ,k,NST)+rr*vflux(ia ,ja  ,k,1)
     vm=rg*vflux(ia ,ja-1,k,NST)+rr*vflux(ia ,ja-1,k,1)
! start ifs code
#if defined ifs
#if defined timeanalyt || timestep
    do n=1,NST
     wflux(k,n) = wflux(k-1,n) - ff * &
     ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n)  &
     + (dzt(ia,ja,k,2)-dzt(ia,ja,k,1))/tseas )  ! time change of the mass the in grid box
!     print *,k,ia,ja,wflux(k,n)
    enddo
#else
  wflux(k) = wflux(k-1) - ff * &
     ( uu - um + vv - vm + (dzt(ia,ja,k,2)-dzt(ia,ja,k,1))/tseas )
!     print *,k,ia,ja,wflux(k)
#endif
#endif
!end ifs code

! start ocean code
#ifndef ifs
#ifdef  full_wflux
     wflux(ia,ja,k,1)=wflux(ia,ja,k-1,1) - ff * ( uu - um + vv - vm )
#elif defined timeanalyt || timestep
    do n=1,NST
     wflux(k,n) = wflux(k-1,n) - ff * &
     ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n) )
    enddo
#else
     wflux(k) = wflux(k-1) - ff * ( uu - um + vv - vm )
#endif
#endif
!end ocean code
  end do kloop

#endif
! end 3D code

! start sediment code
#ifdef sediment  
  ! === Godtyckligt vaerde paa kinetiska energin ===
  ! === daer wsed inte laengre paaverkar, 3e6.   ===
  
  k2loop: do k=0,km
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
#elif defined timeanalyt || timestep
    do n=1,NST
     wflux(k,n)=wflux(k,n) +  wsedtemp * dxdy(ia,ja)     ! *dx *dy *deg**2 *cst(jb)
    enddo
#else
     wflux(k)=wflux(k) +  wsedtemp * dxdy(ia,ja)     ! *dx *dy *deg**2 *cst(jb)
#endif
  end do k2loop
#endif   
! end sediment code
  
!#endif
  return
#endif
end subroutine vertvel

 