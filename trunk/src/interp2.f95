
subroutine interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,ns)

  ! === interpolation of the temperature, salinity, ===
  ! === density at the position of the trajectory   ===
! from the ia,ja,ka and ib,ib,kb  
  
  USE mod_param
  USE mod_dens
  IMPLICIT none
  
  real temp,salt,dens
  
  integer ib,jb,kb,ia,ja,ka,ns
  
#ifdef mod2
  stop 2567 ! Kolla på gammal traj.F för OCCAM
#endif
  
  temp=0.5*(tem(ia,ja,ka,ns)+tem(ib,jb,kb,ns))
  salt=0.5*(sal(ia,ja,ka,ns)+sal(ib,jb,kb,ns))
  dens=0.5*(rho(ia,ja,ka,ns)+rho(ib,jb,kb,ns))
!#ifdef ifs
!  temp=temp-273.15
!  if(temp.lt.-99.98) temp=-99.99
  !salt=0.
!  dens=dens/9.81/1000. ! convert geopotential into geopotential height in km
!#endif
  
  return
end subroutine interp2


