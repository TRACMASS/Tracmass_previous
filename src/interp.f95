#ifdef tempsalt

subroutine interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,ns)

!     computes temperature, salinity, density at position of trajectory
!     by interpolating data from the center of eight nearest boxes  
!
!     This subroutine should be improved in order to include time interpolation   

USE mod_grid
USE mod_dens
USE mod_vel
IMPLICIT none

REAL*8  :: x1,y1,z1,ax,ay,az

REAL    :: tppp,tppm,tpmp,tpmm,tmpp,tmpm,tmmp,tmmm
REAL    :: sppp,sppm,spmp,spmm,smpp,smpm,smmp,smmm
REAL    :: rppp,rppm,rpmp,rpmm,rmpp,rmpm,rmmp,rmmm
REAL    :: temp,salt,dens

INTEGER :: ib,jb,kb,ip,im,jp,jm,kp,kn,ns
! determining nearest centers of boxes 
      if(x1.le.dble(ib)-dble(.5)) then
       ip=ib
       im=ib-1
#ifndef regional
       if(im.eq.0) im=imt
#endif
      else
       ip=ib+1
       im=ib
#ifndef regional
       if(ip.gt.imt) ip=1
#endif
      endif
      if(y1.le.dble(jb)-dble(.5)) then
       jp=jb
       jm=jb-1
       if(jm.eq.0) jm=1
      else
       jp=jb+1
       jm=jb
       if(jp.gt.jmt) jp=jmt
      endif

      if(z1.le.dble(kb)-dble(.5)) then
       kp=kb
       kn=kb-1
       if(kn.le.0) kn=1
      else
       kp=kb+1
       if(kp.gt.km) kp=km
       kn=kb
      endif

      ax=dble(ip)-x1
      if(ax.gt.100.d0) then
!       print *,ax,ip,im,x1,ib
       ax=ax-dble(imt)
!       stop 49678
      elseif(ax.lt.-100.d0) then
       ax=ax+dble(imt)
!       stop 49679
      endif
      ay=(dble(jp)-y1)
      az=(dble(kp)-z1)

! temperature, salinity, density calculation 
      tppp=tem(ip,jp,kp,ns)
      tppm=tem(ip,jp,kn,ns)
      tpmp=tem(ip,jm,kp,ns)
      tpmm=tem(ip,jm,kn,ns)
      tmpp=tem(im,jp,kp,ns)
      tmpm=tem(im,jp,kn,ns)
      tmmp=tem(im,jm,kp,ns)
      tmmm=tem(im,jm,kn,ns)
#ifdef roms
      if (tppp.gt.1.e30) tppp=tem(ib,jb,kb,ns)
      if (tppm.gt.1.e30) tppm=tem(ib,jb,kb,ns)
      if (tpmp.gt.1.e30) tpmp=tem(ib,jb,kb,ns)
      if (tpmm.gt.1.e30) tpmm=tem(ib,jb,kb,ns)
      if (tmpp.gt.1.e30) tmpp=tem(ib,jb,kb,ns)
      if (tmpm.gt.1.e30) tmpm=tem(ib,jb,kb,ns)
      if (tmmp.gt.1.e30) tmmp=tem(ib,jb,kb,ns)
      if (tmmm.gt.1.e30) tmmm=tem(ib,jb,kb,ns)
#else
      if (tppp.eq.0.) tppp=tem(ib,jb,kb,ns)
      if (tppm.eq.0.) tppm=tem(ib,jb,kb,ns)
      if (tpmp.eq.0.) tpmp=tem(ib,jb,kb,ns)
      if (tpmm.eq.0.) tpmm=tem(ib,jb,kb,ns)
      if (tmpp.eq.0.) tmpp=tem(ib,jb,kb,ns)
      if (tmpm.eq.0.) tmpm=tem(ib,jb,kb,ns)
      if (tmmp.eq.0.) tmmp=tem(ib,jb,kb,ns)
      if (tmmm.eq.0.) tmmm=tem(ib,jb,kb,ns)
#endif

      sppp=sal(ip,jp,kp,ns)
      sppm=sal(ip,jp,kn,ns)
      spmp=sal(ip,jm,kp,ns)
      spmm=sal(ip,jm,kn,ns)
      smpp=sal(im,jp,kp,ns)
      smpm=sal(im,jp,kn,ns)
      smmp=sal(im,jm,kp,ns)
      smmm=sal(im,jm,kn,ns)
#ifdef roms
      if (sppp.gt.1.e30) sppp=sal(ib,jb,kb,ns)
      if (sppm.gt.1.e30) sppm=sal(ib,jb,kb,ns)
      if (spmp.gt.1.e30) spmp=sal(ib,jb,kb,ns)
      if (spmm.gt.1.e30) spmm=sal(ib,jb,kb,ns)
      if (smpp.gt.1.e30) smpp=sal(ib,jb,kb,ns)
      if (smpm.gt.1.e30) smpm=sal(ib,jb,kb,ns)
      if (smmp.gt.1.e30) smmp=sal(ib,jb,kb,ns)
      if (smmm.gt.1.e30) smmm=sal(ib,jb,kb,ns)
#else
      if (sppp.eq.0.) sppp=sal(ib,jb,kb,ns)
      if (sppm.eq.0.) sppm=sal(ib,jb,kb,ns)
      if (spmp.eq.0.) spmp=sal(ib,jb,kb,ns)
      if (spmm.eq.0.) spmm=sal(ib,jb,kb,ns)
      if (smpp.eq.0.) smpp=sal(ib,jb,kb,ns)
      if (smpm.eq.0.) smpm=sal(ib,jb,kb,ns)
      if (smmp.eq.0.) smmp=sal(ib,jb,kb,ns)
      if (smmm.eq.0.) smmm=sal(ib,jb,kb,ns)
#endif

      rppp=rho(ip,jp,kp,ns)
      rppm=rho(ip,jp,kn,ns)
      rpmp=rho(ip,jm,kp,ns)
      rpmm=rho(ip,jm,kn,ns)
      rmpp=rho(im,jp,kp,ns)
      rmpm=rho(im,jp,kn,ns)
      rmmp=rho(im,jm,kp,ns)
      rmmm=rho(im,jm,kn,ns)
#ifdef roms
      if (rppp.gt.1.e30) rppp=rho(ib,jb,kb,ns)
      if (rppm.gt.1.e30) rppm=rho(ib,jb,kb,ns)
      if (rpmp.gt.1.e30) rpmp=rho(ib,jb,kb,ns)
      if (rpmm.gt.1.e30) rpmm=rho(ib,jb,kb,ns)
      if (rmpp.gt.1.e30) rmpp=rho(ib,jb,kb,ns)
      if (rmpm.gt.1.e30) rmpm=rho(ib,jb,kb,ns)
      if (rmmp.gt.1.e30) rmmp=rho(ib,jb,kb,ns)
      if (rmmm.gt.1.e30) rmmm=rho(ib,jb,kb,ns)
#else
      if (rppp.eq.0.) rppp=rho(ib,jb,kb,ns)
      if (rppm.eq.0.) rppm=rho(ib,jb,kb,ns)
      if (rpmp.eq.0.) rpmp=rho(ib,jb,kb,ns)
      if (rpmm.eq.0.) rpmm=rho(ib,jb,kb,ns)
      if (rmpp.eq.0.) rmpp=rho(ib,jb,kb,ns)
      if (rmpm.eq.0.) rmpm=rho(ib,jb,kb,ns)
      if (rmmp.eq.0.) rmmp=rho(ib,jb,kb,ns)
      if (rmmm.eq.0.) rmmm=rho(ib,jb,kb,ns)
#endif

      temp=tppp*(1.-ax)*(1.-ay)*(1.-az) &
        + tmpp*    ax *(1.-ay)*(1.-az) &
        + tpmp*(1.-ax)*    ay *(1.-az) &
        + tmmp*    ax *    ay *(1.-az) &
        + tppm*(1.-ax)*(1.-ay)*    az  &
        + tmpm*    ax *(1.-ay)*    az  &
        + tpmm*(1.-ax)*    ay *    az  &
        + tmmm*    ax *    ay *    az  
 
      salt=sppp*(1.-ax)*(1.-ay)*(1.-az) & 
        + smpp*    ax *(1.-ay)*(1.-az) &
        + spmp*(1.-ax)*    ay *(1.-az) &  
        + smmp*    ax *    ay *(1.-az) &
        + sppm*(1.-ax)*(1.-ay)*    az  &
        + smpm*    ax *(1.-ay)*    az  &
        + spmm*(1.-ax)*    ay *    az  &
        + smmm*    ax *    ay *    az  
 
      dens=rppp*(1.-ax)*(1.-ay)*(1.-az) &
        + rmpp*    ax *(1.-ay)*(1.-az) &
        + rpmp*(1.-ax)*    ay *(1.-az) &
        + rmmp*    ax *    ay *(1.-az) &
        + rppm*(1.-ax)*(1.-ay)*    az  &
        + rmpm*    ax *(1.-ay)*    az  &
        + rpmm*(1.-ax)*    ay *    az  &
        + rmmm*    ax *    ay *    az

return
end subroutine interp

#endif

