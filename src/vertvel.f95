
subroutine vertvel(ia,iam,ja,ka)
  
  ! === Computes the vertical velocity by integrating ===
  ! === the continuity eq. from the bottom            ===
  ! === for the nsm and nsp velocity time steps       ===
  USE mod_vel,              only: uflux, vflux, wflux, nsm, nsp, ff
  USE mod_time,             only: intrpr, intrpg, tseas
  USE mod_active_particles, only: upr
  USE mod_grid
#ifdef sediment
  USE mod_sed
  USE mod_orbital
#endif
  
  IMPLICIT none
    
  real*8                                     :: uu, um, vv, vm
  integer                                    :: ia, iam, ja, ka, k,n
  integer                                    :: n1, n2

  REAL                                       :: kin 
  
#if defined twodim || explicit_w
  return
#else
   
  n1=min(nsm,nsp)
  n2=max(nsm,nsp)
  kloop: do k=1,ka
#if defined zgrid3D
     do n=n1,n2
        ! time change of the mass the in grid box
        wflux(k,n) = wflux(k-1,n) - ff * &
             (  uflux(ia,ja,k,n) - uflux(iam, ja,   k, n)   & 
              + vflux(ia,ja,k,n) - vflux(ia,  ja-1, k, n)   & 
#if defined ifs
              + (dzt(ia,ja,k,n2)-dzt(ia,ja,k,n1))*dxdy(ia,ja)/tseas )
#else
              - (dzt(ia,ja,k,n2)-dzt(ia,ja,k,n1))*dxdy(ia,ja)/tseas )
#endif
              ! ska det vara plus för ifs och minus för orca???
    enddo
#else
#if defined  full_wflux
     uu = intrpg * uflux(ia ,ja  ,k,nsp) + intrpr * uflux(ia ,ja  ,k,nsm)
     um = intrpg * uflux(iam,ja  ,k,nsp) + intrpr * uflux(iam,ja  ,k,nsm)
     vv = intrpg * vflux(ia ,ja  ,k,nsp) + intrpr * vflux(ia ,ja  ,k,nsm)
     vm = intrpg * vflux(ia ,ja-1,k,nsp) + intrpr * vflux(ia ,ja-1,k,nsm)
     wflux(ia,ja,k,nsm)=wflux(ia,ja,k-1,nsm) - ff * ( uu - um + vv - vm )
#else
    do n=n1,n2
     wflux(k,n) = wflux(k-1,n) - ff * &
     ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n) )
    enddo
#endif
#endif
!end ocean code
  end do kloop


! Make sure the vertical velocity is always zero at the bottom and below
#if defined ifs
wflux(0,:) = 0.d0
#elif defined orca1 || orca12
do k=0,KM-kmt(ia,ja)
  wflux(k,:)=0.d0
enddo
#endif

#endif

#ifdef sediment  
  ! === Godtyckligt vaerde paa kinetiska energin ===
  ! === daer wsed inte laengre paaverkar, 3e6.   ===
  
  k2loop: do k=0,km
     wsedtemp=0.d0
     kin=(uflux(ia,ja,k,nsm)*uflux(ia,ja,k,nsm)+ &
          vflux(ia,ja,k,nsm)*vflux(ia,ja,k,nsm))*0.5d0
     !if (kin.le.3000000) then   !for RCO
     !wsedtemp=wsed*(3000000-kin)/3000000
     if (kin.le.kincrit) then   !for SKB
        wsedtemp=wsed*(kincrit-kin)/kincrit
     endif
#ifdef explicit_w || full_wflux
     wflux(k)=wflux(ia,ja,k,nsm) +  wsedtemp * dxdy(ia,ja) 
#else
    do n=n1,n2
     wflux(k,n)=wflux(k,n) +  wsedtemp * dxdy(ia,ja)
    enddo
#endif
  end do k2loop
#endif   
! end sediment code
  
  return

end subroutine vertvel

 
