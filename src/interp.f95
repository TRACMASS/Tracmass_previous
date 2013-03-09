#ifdef tempsalt

subroutine interp(ib,jb,kb,x1,y1,z1,ns)

!     computes temperature, salinity, density at position of trajectory
!     by interpolating data from the center of eight nearest boxes
!
!     This subroutine should be improved in order to include time interpolation

USE mod_grid
USE mod_dens
USE mod_particle
USE mod_vel
USE mod_loopvars
IMPLICIT none

REAL*8  :: x1,y1,z1,ax,ay,az

REAL    :: tppp,tppm,tpmp,tpmm,tmpp,tmpm,tmmp,tmmm
REAL    :: sppp,sppm,spmp,spmm,smpp,smpm,smmp,smmm
REAL    :: rppp,rppm,rpmp,rpmm,rmpp,rmpm,rmmp,rmmm
REAL    :: zppp,zppm,zpmp,zpmm,zmpp,zmpm,zmmp,zmmm
REAL    :: latpp,latpm,latmp,latmm
REAL    :: lonpp,lonpm,lonmp,lonmm

INTEGER :: ib,jb,kb,ip,im,jp,jm,kp,kn,ns,it

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
      do it = 1,ns
        tppp=tem(ip,jp,kp,it)
        tppm=tem(ip,jp,kn,it)
        tpmp=tem(ip,jm,kp,it)
        tpmm=tem(ip,jm,kn,it)
        tmpp=tem(im,jp,kp,it)
        tmpm=tem(im,jp,kn,it)
        tmmp=tem(im,jm,kp,it)
        tmmm=tem(im,jm,kn,it)
#ifdef roms
        if (tppp.gt.1.e30) tppp=tem(ib,jb,kb,it)
        if (tppm.gt.1.e30) tppm=tem(ib,jb,kb,it)
        if (tpmp.gt.1.e30) tpmp=tem(ib,jb,kb,it)
        if (tpmm.gt.1.e30) tpmm=tem(ib,jb,kb,it)
        if (tmpp.gt.1.e30) tmpp=tem(ib,jb,kb,it)
        if (tmpm.gt.1.e30) tmpm=tem(ib,jb,kb,it)
        if (tmmp.gt.1.e30) tmmp=tem(ib,jb,kb,it)
        if (tmmm.gt.1.e30) tmmm=tem(ib,jb,kb,it)
#else
        if (tppp.eq.0.) tppp=tem(ib,jb,kb,it)
        if (tppm.eq.0.) tppm=tem(ib,jb,kb,it)
        if (tpmp.eq.0.) tpmp=tem(ib,jb,kb,it)
        if (tpmm.eq.0.) tpmm=tem(ib,jb,kb,it)
        if (tmpp.eq.0.) tmpp=tem(ib,jb,kb,it)
        if (tmpm.eq.0.) tmpm=tem(ib,jb,kb,it)
        if (tmmp.eq.0.) tmmp=tem(ib,jb,kb,it)
        if (tmmm.eq.0.) tmmm=tem(ib,jb,kb,it)
#endif

        sppp=sal(ip,jp,kp,it)
        sppm=sal(ip,jp,kn,it)
        spmp=sal(ip,jm,kp,it)
        spmm=sal(ip,jm,kn,it)
        smpp=sal(im,jp,kp,it)
        smpm=sal(im,jp,kn,it)
        smmp=sal(im,jm,kp,it)
        smmm=sal(im,jm,kn,it)
#ifdef roms
        if (sppp.gt.1.e30) sppp=sal(ib,jb,kb,it)
        if (sppm.gt.1.e30) sppm=sal(ib,jb,kb,it)
        if (spmp.gt.1.e30) spmp=sal(ib,jb,kb,it)
        if (spmm.gt.1.e30) spmm=sal(ib,jb,kb,it)
        if (smpp.gt.1.e30) smpp=sal(ib,jb,kb,it)
        if (smpm.gt.1.e30) smpm=sal(ib,jb,kb,it)
        if (smmp.gt.1.e30) smmp=sal(ib,jb,kb,it)
        if (smmm.gt.1.e30) smmm=sal(ib,jb,kb,it)
#else
        if (sppp.eq.0.) sppp=sal(ib,jb,kb,it)
        if (sppm.eq.0.) sppm=sal(ib,jb,kb,it)
        if (spmp.eq.0.) spmp=sal(ib,jb,kb,it)
        if (spmm.eq.0.) spmm=sal(ib,jb,kb,it)
        if (smpp.eq.0.) smpp=sal(ib,jb,kb,it)
        if (smpm.eq.0.) smpm=sal(ib,jb,kb,it)
        if (smmp.eq.0.) smmp=sal(ib,jb,kb,it)
        if (smmm.eq.0.) smmm=sal(ib,jb,kb,it)
#endif

        rppp=rho(ip,jp,kp,it)
        rppm=rho(ip,jp,kn,it)
        rpmp=rho(ip,jm,kp,it)
        rpmm=rho(ip,jm,kn,it)
        rmpp=rho(im,jp,kp,it)
        rmpm=rho(im,jp,kn,it)
        rmmp=rho(im,jm,kp,it)
        rmmm=rho(im,jm,kn,it)
#ifdef roms
        if (rppp.gt.1.e30) rppp=rho(ib,jb,kb,it)
        if (rppm.gt.1.e30) rppm=rho(ib,jb,kb,it)
        if (rpmp.gt.1.e30) rpmp=rho(ib,jb,kb,it)
        if (rpmm.gt.1.e30) rpmm=rho(ib,jb,kb,it)
        if (rmpp.gt.1.e30) rmpp=rho(ib,jb,kb,it)
        if (rmpm.gt.1.e30) rmpm=rho(ib,jb,kb,it)
        if (rmmp.gt.1.e30) rmmp=rho(ib,jb,kb,it)
        if (rmmm.gt.1.e30) rmmm=rho(ib,jb,kb,it)
#else
        if (rppp.eq.0.) rppp=rho(ib,jb,kb,it)
        if (rppm.eq.0.) rppm=rho(ib,jb,kb,it)
        if (rpmp.eq.0.) rpmp=rho(ib,jb,kb,it)
        if (rpmm.eq.0.) rpmm=rho(ib,jb,kb,it)
        if (rmpp.eq.0.) rmpp=rho(ib,jb,kb,it)
        if (rmpm.eq.0.) rmpm=rho(ib,jb,kb,it)
        if (rmmp.eq.0.) rmmp=rho(ib,jb,kb,it)
        if (rmmm.eq.0.) rmmm=rho(ib,jb,kb,it)
#endif

#ifdef larval_fish
        zppp=z_r(ip,jp,kp,it)
        zppm=z_r(ip,jp,kn,it)
        zpmp=z_r(ip,jm,kp,it)
        zpmm=z_r(ip,jm,kn,it)
        zmpp=z_r(im,jp,kp,it)
        zmpm=z_r(im,jp,kn,it)
        zmmp=z_r(im,jm,kp,it)
        zmmm=z_r(im,jm,kn,it)
#endif

      temp2(it)=tppp*(1.-ax)*(1.-ay)*(1.-az) &
              + tmpp*    ax *(1.-ay)*(1.-az) &
              + tpmp*(1.-ax)*    ay *(1.-az) &
              + tmmp*    ax *    ay *(1.-az) &
              + tppm*(1.-ax)*(1.-ay)*    az  &
              + tmpm*    ax *(1.-ay)*    az  &
              + tpmm*(1.-ax)*    ay *    az  &
              + tmmm*    ax *    ay *    az

      salt2(it)=sppp*(1.-ax)*(1.-ay)*(1.-az) &
              + smpp*    ax *(1.-ay)*(1.-az) &
              + spmp*(1.-ax)*    ay *(1.-az) &
              + smmp*    ax *    ay *(1.-az) &
              + sppm*(1.-ax)*(1.-ay)*    az  &
              + smpm*    ax *(1.-ay)*    az  &
              + spmm*(1.-ax)*    ay *    az  &
              + smmm*    ax *    ay *    az

      dens2(it)=rppp*(1.-ax)*(1.-ay)*(1.-az) &
              + rmpp*    ax *(1.-ay)*(1.-az) &
              + rpmp*(1.-ax)*    ay *(1.-az) &
              + rmmp*    ax *    ay *(1.-az) &
              + rppm*(1.-ax)*(1.-ay)*    az  &
              + rmpm*    ax *(1.-ay)*    az  &
              + rpmm*(1.-ax)*    ay *    az  &
              + rmmm*    ax *    ay *    az

#ifdef larval_fish
      depth2(it)=zppp*(1.-ax)*(1.-ay)*(1.-az) &
               + zmpp*    ax *(1.-ay)*(1.-az) &
               + zpmp*(1.-ax)*    ay *(1.-az) &
               + zmmp*    ax *    ay *(1.-az) &
               + zppm*(1.-ax)*(1.-ay)*    az  &
               + zmpm*    ax *(1.-ay)*    az  &
               + zpmm*(1.-ax)*    ay *    az  &
               + zmmm*    ax *    ay *    az
#endif
      end do

      if (ns == 1) then
        temp=temp2(1)
        salt=salt2(1)
        dens=dens2(1)
      else
        temp=rg*temp2(nsp)+rr*temp2(nsm)
        salt=rg*salt2(nsp)+rr*salt2(nsm)
        dens=rg*dens2(nsp)+rr*dens2(nsm)
      endif

#ifdef larval_fish
      if (ns == 1) then
        depth1=depth2(1)
      else
        depth1=rg*depth2(nsp)+rr*depth2(nsm)
      endif

      latpp=lat(ip,jp)
      latpm=lat(ip,jm)
      latmp=lat(im,jp)
      latmm=lat(im,jm)
      lonpp=lon(ip,jp)
      lonpm=lon(ip,jm)
      lonmp=lon(im,jp)
      lonmm=lon(im,jm)

      flon1 = lonpp*(1.-ax)*(1.-ay) &
            + lonmp*    ax *(1.-ay) &
            + lonpm*(1.-ax)*    ay  &
            + lonmm*    ax *    ay
      flat1 = latpp*(1.-ax)*(1.-ay) &
            + latmp*    ax *(1.-ay) &
            + latpm*(1.-ax)*    ay  &
            + latmm*    ax *    ay
#endif

return
end subroutine interp

#endif

