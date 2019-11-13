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
REAL    :: tp,tm,sp,sm,rp,rm

INTEGER :: ib,jb,kb,ip,im,jp,jm,kp,kn,ns
! determining nearest centers of boxes 

PRINT *, x1, y1, z1

 IF(z1.LE.DBLE(kb)-DBLE(.5)) THEN
    kp=kb
    kn=kb-1
    IF(kn.LE.0) kn=1
 ELSE
    kp=kb+1
    IF(kp.GT.km) kp=km
    kn=kb
 ENDIF

 az=(DBLE(kp)-(z1 + .5)) !Sara
 
 tp = tem(ib,jb,kp,ns)
 tm = tem(ib,jb,kn,ns)
 sp = sal(ib,jb,kp,ns)
 sm = sal(ib,jb,kn,ns)
 rp = rho(ib,jb,kp,ns)
 rm = rho(ib,jb,kn,ns)
 
 IF (kp == km) THEN
    tp = tem(ib,jb,kb,ns)
    tm = tem(ib,jb,kb,ns)
 ENDIF
 
 temp = tp*(1.-az) + tm*az
 salt = sp*(1.-az) + sm*az
 dens = rp*(1.-az) + rm*az
      


 RETURN
end subroutine interp

#endif

