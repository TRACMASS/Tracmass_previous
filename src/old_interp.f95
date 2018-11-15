#ifdef tempsalt

subroutine interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,ns)

!     computes temperature, salinity, density at position of trajectory
!     by interpolating data from the center of eight nearest boxes  
!
!     This subroutine should be improved in order to include time interpolation   

USE mod_grid
USE mod_vel
USE mod_dens
USE mod_vel
USE mod_tempsalt
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
       if(im.eq.0) im=imt
      else
       ip=ib+1
       im=ib
       if(ip.gt.imt) ip=1
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
!      tppp=tem(ip,jp,kp,ns)
!      if(tppp.eq.0.) tppp=tem(ib,jb,kb,ns)
!      tppm=tem(ip,jp,kn,ns)
!      if(tppm.eq.0.) tppm=tem(ib,jb,kb,ns)
!      tpmp=tem(ip,jm,kp,ns)
!      if(tpmp.eq.0.) tpmp=tem(ib,jb,kb,ns)
!      tpmm=tem(ip,jm,kn,ns)
!      if(tpmm.eq.0.) tpmm=tem(ib,jb,kb,ns)
!      tmpp=tem(im,jp,kp,ns)
!      if(tmpp.eq.0.) tmpp=tem(ib,jb,kb,ns)
!      tmpm=tem(im,jp,kn,ns)
!      if(tmpm.eq.0.) tmpm=tem(ib,jb,kb,ns)
!      tmmp=tem(im,jm,kp,ns)
!      if(tmmp.eq.0.) tmmp=tem(ib,jb,kb,ns)
!      tmmm=tem(im,jm,kn,ns)
!      if(tmmm.eq.0.) tmmm=tem(ib,jb,kb,ns)
!
!      sppp=sal(ip,jp,kp,ns)
!      if(sppp.eq.0.) sppp=sal(ib,jb,kb,ns)
!      sppm=sal(ip,jp,kn,ns)
!      if(sppm.eq.0.) sppm=sal(ib,jb,kb,ns)
!      spmp=sal(ip,jm,kp,ns)
!      if(spmp.eq.0.) spmp=sal(ib,jb,kb,ns)
!      spmm=sal(ip,jm,kn,ns)
!      if(spmm.eq.0.) spmm=sal(ib,jb,kb,ns)
!      smpp=sal(im,jp,kp,ns)
!      if(smpp.eq.0.) smpp=sal(ib,jb,kb,ns)
!      smpm=sal(im,jp,kn,ns)
!      if(smpm.eq.0.) smpm=sal(ib,jb,kb,ns)
!      smmp=sal(im,jm,kp,ns)
!      if(smmp.eq.0.) smmp=sal(ib,jb,kb,ns)
!      smmm=sal(im,jm,kn,ns)
!      if(smmm.eq.0.) smmm=sal(ib,jb,kb,ns)
!
!      rppp=rho(ip,jp,kp,ns)
!      if(rppp.eq.0.) rppp=rho(ib,jb,kb,ns)
!      rppm=rho(ip,jp,kn,ns)
!      if(rppm.eq.0.) rppm=rho(ib,jb,kb,ns)
!      rpmp=rho(ip,jm,kp,ns)
!      if(rpmp.eq.0.) rpmp=rho(ib,jb,kb,ns)
!      rpmm=rho(ip,jm,kn,ns)
!      if(rpmm.eq.0.) rpmm=rho(ib,jb,kb,ns)
!      rmpp=rho(im,jp,kp,ns)
!      if(rmpp.eq.0.) rmpp=rho(ib,jb,kb,ns)
!      rmpm=rho(im,jp,kn,ns)
!      if(rmpm.eq.0.) rmpm=rho(ib,jb,kb,ns)
!      rmmp=rho(im,jm,kp,ns)
!      if(rmmp.eq.0.) rmmp=rho(ib,jb,kb,ns)
!      rmmm=rho(im,jm,kn,ns)
!      if(rmmm.eq.0.) rmmm=rho(ib,jb,kb,ns)

      tppp=tem(ip,jp,kp,ns)
      sppp=sal(ip,jp,kp,ns)
      rppp=rho(ip,jp,kp,ns)
      if(tppp==0. .and. sppp==0.) then
       tppp=tem(ip,jp,kn,ns)
       sppp=sal(ip,jp,kn,ns)
       rppp=rho(ip,jp,kn,ns)
      endif
      
      tppm=tem(ip,jp,kn,ns)
      sppm=sal(ip,jp,kn,ns)
      rppm=rho(ip,jp,kn,ns)
      if(tppm==0. .and. sppm==0.) then
       tppm=tem(ip,jp,kn,ns)
       sppm=sal(ip,jp,kn,ns)
       rppm=rho(ip,jp,kn,ns)
      endif
      
      tpmp=tem(ip,jm,kp,ns)
      spmp=sal(ip,jm,kp,ns)
      rpmp=rho(ip,jm,kp,ns)
      if(tpmp==0. .and. spmp==0.) then
       tpmp=tem(ip,jp,kn,ns)
       spmp=sal(ip,jp,kn,ns)
       rpmp=rho(ip,jp,kn,ns)
      endif
      
      tpmm=tem(ip,jm,kn,ns)
      spmm=sal(ip,jm,kn,ns)
      rpmm=rho(ip,jm,kn,ns)
      if(tpmm==0. .and. spmm==0.) then
       tpmm=tem(ip,jp,kn,ns)
       spmm=sal(ip,jp,kn,ns)
       rpmm=rho(ip,jp,kn,ns)
      endif
      
      tmpp=tem(im,jp,kp,ns)
      smpp=sal(im,jp,kp,ns)
      rmpp=rho(im,jp,kp,ns)
      if(tmpp==0. .and. smpp==0.) then
       tmpp=tem(ip,jp,kn,ns)
       smpp=sal(ip,jp,kn,ns)
       rmpp=rho(ip,jp,kn,ns)
      endif
      
      tmpm=tem(im,jp,kn,ns)
      smpm=sal(im,jp,kn,ns)
      rmpm=rho(im,jp,kn,ns)
      if(tmpm==0. .and. smpm==0.) then
       tmpm=tem(ip,jp,kn,ns)
       smpm=sal(ip,jp,kn,ns)
       rmpm=rho(ip,jp,kn,ns)
      endif
      
      tmmp=tem(im,jm,kp,ns)
      smmp=sal(im,jm,kp,ns)
      rmmp=rho(im,jm,kp,ns)
      if(tmmp==0. .and. smmp==0.) then
       tmmp=tem(ip,jp,kn,ns)
       smmp=sal(ip,jp,kn,ns)
       rmmp=rho(ip,jp,kn,ns)
      endif
      
      tmmm=tem(im,jm,kn,ns)
      smmm=sal(im,jm,kn,ns)
      rmmm=rho(im,jm,kn,ns)
      if(tmmm==0. .and. smmm ==0.) then
       tmmm=tem(ip,jp,kn,ns)
       smmm=sal(ip,jp,kn,ns)
       rmmm=rho(ip,jp,kn,ns)
      endif



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

