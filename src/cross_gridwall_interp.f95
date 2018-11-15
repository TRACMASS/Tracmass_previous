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


 IF(x1 == INT(x1)) THEN
    !PRINT *, 'x', x1, ib
    !PRINT *, 'y', y1, jb --> y1 < jb
    !PRINT *, 'z', z1, kb --> z1 < kb

    IF (x1 == ib) THEN ! Crossing the left grid wall

       temp = 0.5*(tem(ib,jb,kb,ns)+tem(ib+1,jb,kb,ns))
       salt = 0.5*(sal(ib,jb,kb,ns)+sal(ib+1,jb,kb,ns))
       dens = 0.5*(rho(ib,jb,kb,ns)+rho(ib+1,jb,kb,ns))

    ELSEIF (x1 < ib) THEN ! Crossing the right grid wall

       temp = 0.5*(tem(ib-1,jb,kb,ns)+tem(ib,jb,kb,ns))
       salt = 0.5*(sal(ib-1,jb,kb,ns)+sal(ib,jb,kb,ns))
       dens = 0.5*(rho(ib-1,jb,kb,ns)+rho(ib,jb,kb,ns))

    ELSE
       PRINT *, 'Error: The trajectory is not in the expected gridbox!'
    ENDIF

 ELSEIF (y1 == INT(y1)) THEN
    
    IF (y1 == jb) THEN ! Crossing the left grid wall

       temp = 0.5*(tem(ib,jb,kb,ns)+tem(ib,jb+1,kb,ns))
       salt = 0.5*(sal(ib,jb,kb,ns)+sal(ib,jb+1,kb,ns))
       dens = 0.5*(rho(ib,jb,kb,ns)+rho(ib,jb+1,kb,ns))

    ELSEIF (y1 < jb) THEN ! Crossing the right grid wall

       temp = 0.5*(tem(ib,jb-1,kb,ns)+tem(ib,jb,kb,ns))
       salt = 0.5*(sal(ib,jb-1,kb,ns)+sal(ib,jb,kb,ns))
       dens = 0.5*(rho(ib,jb-1,kb,ns)+rho(ib,jb,kb,ns))

    ELSE
       PRINT *, 'Error: The trajectory is not in the expected gridbox!'
    ENDIF
   
 ELSEIF (z1 == INT(z1)) THEN

    IF (z1 == kb) THEN ! Crossing the left grid wall

       temp = 0.5*(tem(ib,jb,kb,ns)+tem(ib,jb,kb+1,ns))
       salt = 0.5*(sal(ib,jb,kb,ns)+sal(ib,jb,kb+1,ns))
       dens = 0.5*(rho(ib,jb,kb,ns)+rho(ib,jb,kb+1,ns))

    ELSEIF (z1 < kb) THEN ! Crossing the right grid wall

       temp = 0.5*(tem(ib,jb,kb-1,ns)+tem(ib,jb,kb,ns))
       salt = 0.5*(sal(ib,jb,kb-1,ns)+sal(ib,jb,kb,ns))
       dens = 0.5*(rho(ib,jb,kb-1,ns)+rho(ib,jb,kb,ns))

    ELSE
       PRINT *, 'Error: The trajectory is not in the expected gridbox!'
    ENDIF

 ENDIF

      


 RETURN
end subroutine interp

#endif

