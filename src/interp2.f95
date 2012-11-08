#ifdef tempsalt
subroutine interp2(i,j,k,temp,salt,dens)

  ! === NO interpolation of the temperature, salinity, and density===
  ! === just their value at the centre of the T-box interpolated in time  ===
  ! === can be called as either ia,ja,ka or ib,jb,kb
  ! === used to calculate the thermohaline stram function with -Dstream_thermohaline  
  
  USE mod_param
  USE mod_loopvars
  USE mod_time
  USE mod_dens
  IMPLICIT none
  
  real temp,salt,dens
  
  integer i,j,k
  
  temp=rbg*tem(i,j,k,nsp)+rb*tem(i,j,k,nsm)
  salt=rbg*sal(i,j,k,nsp)+rb*sal(i,j,k,nsm)
  dens=rbg*rho(i,j,k,nsp)+rb*rho(i,j,k,nsm)
  
  return
end subroutine interp2
#endif

