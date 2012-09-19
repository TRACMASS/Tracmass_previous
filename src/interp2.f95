#ifdef tempsalt
subroutine interp2(i,j,k,temp,salt,dens)

  ! === NO interpolation of the temperature, salinity, and density===
  ! === just their value at the centre of the T-box interpolated in time  ===
  ! === can be called as either ia,ja,ka or ib,jb,kb
  ! === used to calculate the thermohaline stram function with -Dstream_thermohaline  
  
  USE mod_param
  USE mod_loopvars
  USE mod_dens
  IMPLICIT none
  
  real*8 temp,salt,dens
  
  integer i,j,k
  
  temp=rbg*tem(i,j,k,NST)+rb*tem(i,j,k,1)
  salt=rbg*sal(i,j,k,NST)+rb*sal(i,j,k,1)
  dens=rbg*rho(i,j,k,NST)+rb*rho(i,j,k,1)
  
  return
end subroutine interp2
#endif

