
subroutine vertvel(rr,ia,iam,ja,ka)

  
  ! === Computes the vertical velocity by integrating ===
  ! === the continuity eq. from the bottom            ===
  
  USE mod_param
  USE mod_vel
  USE mod_turb
#if defined ifs || defined larval_fish
  USE mod_grid
#endif
#ifdef sediment
  USE mod_sed
  USE mod_orbital
  USE mod_grid
#endif
#if defined larval_fish
  USE mod_fish
#endif /*fish*/
  
  IMPLICIT none
  
#if defined sediment
  REAL*8 wsedtemp
  REAL*8 kin
#endif
  
  real*8 rr,rg,uu,um,vv,vm
  integer ia,iam,ja,ka,k,n
    
#ifndef explicit_w
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
#if defined ifs || defined roms
    do n=1,NST
     wflux(k,n) = wflux(k-1,n) - ff * &
     ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n)  &
     + (dzt(ia,ja,k,2)-dzt(ia,ja,k,1))*dxdy(ia,ja)/tseas )  ! time change of the mass the in grid box
    enddo
#endif
!end ifs code

! start ocean code
#ifndef ifs
#ifdef  full_wflux
     wflux(ia,ja,k,1)=wflux(ia,ja,k-1,1) - ff * ( uu - um + vv - vm )
#else
    do n=1,NST
     wflux(k,n) = wflux(k-1,n) - ff * &
     ( uflux(ia,ja,k,n) - uflux(iam,ja,k,n) + vflux(ia,ja,k,n) - vflux(ia,ja-1,k,n) )
    enddo
#endif
#endif
!end ocean code
  end do kloop

#ifdef ifs
wflux(0,:) = 0.d0
wflux(km,:) = 0.d0
#endif

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
#if defined full_wflux || explicit_w
     wflux(ia,ja,k,1)=wflux(ia,ja,k,1) + wsedtemp * dxdy(ia,ja)   ! *dx *dy *deg**2 
#else
    do n=1,NST
     wflux(k,n)=wflux(k,n) +  wsedtemp * dxdy(ia,ja)     ! *dx *dy *deg**2 
    enddo
#endif
  end do k2loop
#endif   
! end sediment code
  
#ifdef larval_fish
  ! === fisk!   ===

  k3loop: do k=0,km
#if defined full_wflux || defined explicit_w
     wflux(ia,ja,k,1)=wflux(ia,ja,k,1) + wfish * dxdy(ia,ja)   ! *dx *dy *deg**2
!     wflux(k)=wflux(ia,ja,k,1) +  wfish * dxdy(ia,ja)     ! *dx *dy *deg**2
#else
    do n=1,NST
     wflux(k,n)=wflux(k,n) +  wfish * dxdy(ia,ja)     ! *dx *dy *deg**2
    enddo
#endif
  end do k3loop
#endif
! end fish code
#endif /* !explicit_w */

  return
end subroutine vertvel
