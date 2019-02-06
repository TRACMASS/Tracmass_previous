#ifdef tempsalt
subroutine interp2(i,j,k,temp,salt,dens)

  ! === NO interpolation of the temperature, salinity, and density===
  ! === just their value at the centre of the T-box interpolated in time  ===
  ! === can be called as either ia,ja,ka or ib,jb,kb
  ! === used to calculate the thermohaline stram function with -Dstream_thermohaline  
  
  USE mod_grid
  USE mod_param
  USE mod_loopvars
  USE mod_time
  USE mod_vel
  USE mod_dens
  USE mod_tempsalt
  IMPLICIT none
  
  real temp,salt,dens
  
  integer i,j,k

  intrpbg=dmod(ts,1.d0) 
  intrpb =1.d0-intrpbg

  temp=intrpbg*tem(i,j,k,nsp)+intrpb*tem(i,j,k,nsm)
  salt=intrpbg*sal(i,j,k,nsp)+intrpb*sal(i,j,k,nsm)
  dens=intrpbg*rho(i,j,k,nsp)+intrpb*rho(i,j,k,nsm)

  return
end subroutine interp2
#endif

