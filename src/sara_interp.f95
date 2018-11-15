#ifdef tempsalt

subroutine interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,ns)

!     computes temperature, salinity, density at position of trajectory
!     by interpolating data from the center of the nearest boxes in depth.  
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
REAL    :: tb,tm,sb,sm,rb,rm

INTEGER :: ib,jb,kb,ip,im,jp,jm,kp,kn,ns

! ----- CROSSING X -----
 IF(x1 == INT(x1)) THEN

    IF (x1 == ib) THEN ! Crossing the right grid wall
       
       ! The box the trajectory crosses. ib,jb,kb -> ib+1,jb,kb. None of these can be land, since the trajectory can cross.
       tb = 0.5 * ( tem(ib,jb,kb,ns)   + tem(ib+1,jb,kb,ns)   )
       sb = 0.5 * ( sal(ib,jb,kb,ns)   + sal(ib+1,jb,kb,ns)   )
       rb = 0.5 * ( rho(ib,jb,kb,ns)   + rho(ib+1,jb,kb,ns)   )
       
       ! Determine position in z. 
       IF ( z1 < (kb - 0.5) ) THEN ! closer to kb-1 than kb+1
       
          IF (sal(ib,jb,kb-1,ns) == 0 .AND. tem(ib,jb,kb-1,ns) == 0) THEN ! Land in gridbox ib,jb,kb-1
             tm = tem(ib+1,jb,kb-1,ns)
             sm = sal(ib+1,jb,kb-1,ns)
             rm = rho(ib+1,jb,kb-1,ns)

          ELSEIF (sal(ib+1,jb,kb-1,ns) == 0 .AND. tem(ib+1,jb,kb-1,ns) == 0) THEN ! Land in gridbox ib+1,jb,kb-1
             tm = tem(ib,jb,kb-1,ns)
             sm = sal(ib,jb,kb-1,ns)
             rm = rho(ib,jb,kb-1,ns)

          ELSE ! No land or land in both boxes.
             tm = 0.5 * ( tem(ib,jb,kb-1,ns) + tem(ib+1,jb,kb-1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb-1,ns) + sal(ib+1,jb,kb-1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb-1,ns) + rho(ib+1,jb,kb-1,ns) )

          ENDIF
         
          ! --- Computes the total temperature of the trajektory ---
	  ! ---    by giving weigt to tb,sb,rb and tm,sm,rm      ---

          IF (sm == 0 .AND. tm == 0) THEN ! If both gridboxes are land.s
             temp = tb
             salt = sb
             dens = rb

          ELSE ! No land
             az = kb - z1 - 0.5 ! Distance to trajectory from kb.(?)

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF


       ELSE ! If the trajectory is closer to kb + 1

          IF (sal(ib,jb,kb+1,ns) == 0 .AND. tem(ib,jb,kb+1,ns) == 0) THEN ! Land in gridbox ib,jb,kb+1
             tm = tem(ib+1,jb,kb+1,ns)
             sm = sal(ib+1,jb,kb+1,ns)
             rm = rho(ib+1,jb,kb+1,ns)

          ELSEIF (sal(ib+1,jb,kb+1,ns) == 0 .AND. tem(ib+1,jb,kb+1,ns) == 0) THEN ! Land in gridbox ib+1,jb,kb+1
             tm = tem(ib,jb,kb+1,ns)
             sm = sal(ib,jb,kb+1,ns)
             rm = rho(ib,jb,kb+1,ns)

          ELSE ! No land or land in both boxes.
             tm = 0.5 * ( tem(ib,jb,kb+1,ns) + tem(ib+1,jb,kb+1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb+1,ns) + sal(ib+1,jb,kb+1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb+1,ns) + rho(ib+1,jb,kb+1,ns) ) 
          
          ENDIF


	! --- Computes the total temperature of the trajectory ---
	! ---    by giving weigt to tb,sb,rb and tm,sm,rm      ---
	
          IF (sm == 0 .AND. tm == 0) THEN ! Land in both gridboxes. 
             temp = tb
             salt = sb
             dens = rb

          ELSE	! No land or land in only one box.
             az = z1 - kb + 0.5

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF
          
       ENDIF
       
    ! Crossing the left grid wall
    ELSEIF (x1 < ib) THEN 

       ! The box the trajectory crosses, same for both z-cases
       tb = 0.5 * ( tem(ib,jb,kb,ns)   + tem(ib-1,jb,kb,ns)   )
       sb = 0.5 * ( sal(ib,jb,kb,ns)   + sal(ib-1,jb,kb,ns)   )
       rb = 0.5 * ( rho(ib,jb,kb,ns)   + rho(ib-1,jb,kb,ns)   )

       IF ( z1 < kb - 0.5 ) THEN ! Determening position of trajectory in z.
          
          IF (sal(ib,jb,kb-1,ns) == 0 .AND. tem(ib,jb,kb-1,ns) == 0) THEN ! Land
             tm = tem(ib-1,jb,kb-1,ns)
             sm = sal(ib-1,jb,kb-1,ns)
             rm = rho(ib-1,jb,kb-1,ns)
             
          ELSEIF (sal(ib-1,jb,kb-1,ns) == 0 .AND. tem(ib-1,jb,kb-1,ns) == 0) THEN ! Land
             tm = tem(ib,jb,kb-1,ns)
             sm = sal(ib,jb,kb-1,ns)
             rm = rho(ib,jb,kb-1,ns)

          ELSE ! No land or both land. 
             tm = 0.5 * ( tem(ib,jb,kb-1,ns) + tem(ib-1,jb,kb-1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb-1,ns) + sal(ib-1,jb,kb-1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb-1,ns) + rho(ib-1,jb,kb-1,ns) )

          ENDIF

          
          IF (sm == 0 .AND. tm == 0) THEN ! if both land
             temp = tb
             salt = sb
             dens = rb

          ELSE ! If no land
             az = z1 - kb + 0.5

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF

       ELSE

          IF (sal(ib,jb,kb+1,ns) == 0 .AND. tem(ib,jb,kb+1,ns) == 0) THEN ! Land
             tm = tem(ib-1,jb,kb+1,ns)
             sm = sal(ib-1,jb,kb+1,ns)
             rm = rho(ib-1,jb,kb+1,ns)
             
          ELSEIF (sal(ib-1,jb,kb+1,ns) == 0 .AND. tem(ib-1,jb,kb+1,ns) == 0) THEN ! Land
             tm = tem(ib,jb,kb+1,ns)
             sm = sal(ib,jb,kb+1,ns)
             rm = rho(ib,jb,kb+1,ns)

          ELSE ! No land or both land. 
             tm = 0.5 * ( tem(ib,jb,kb+1,ns) + tem(ib-1,jb,kb+1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb+1,ns) + sal(ib-1,jb,kb+1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb+1,ns) + rho(ib-1,jb,kb+1,ns) )

          ENDIF

          
          IF (sm == 0 .AND. tm == 0) THEN ! if both land
             temp = tb
             salt = sb
             dens = rb

          ELSE ! If no land
             az = z1 - kb + 0.5

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF
       ENDIF


    ELSE
       PRINT *, 'Error: The trajectory is not in the expected gridbox!'
    ENDIF





! ############# CROSSING Y #############
 ELSEIF (y1 == INT(y1)) THEN
 !   
    IF (y1 == jb) THEN ! Crossing the right grid wall
!
! The box the trajectory crosses, will be the same for both z-cases
       tb = 0.5 * ( tem(ib,jb,kb,ns)   + tem(ib,jb+1,kb,ns)   )
       sb = 0.5 * ( sal(ib,jb,kb,ns)   + sal(ib,jb+1,kb,ns)   )
       rb = 0.5 * ( rho(ib,jb,kb,ns)   + rho(ib,jb+1,kb,ns)   )

 ! Determine position in z. 
       IF ( z1 < (kb - 0.5) ) THEN ! closer to kb-1 than kb+1
       
          IF (sal(ib,jb,kb-1,ns) == 0 .AND. tem(ib,jb,kb-1,ns) == 0) THEN ! Land
             tm = tem(ib,jb+1,kb-1,ns)
             sm = sal(ib,jb+1,kb-1,ns)
             rm = rho(ib,jb+1,kb-1,ns)

          ELSEIF (sal(ib,jb+1,kb-1,ns) == 0 .AND. tem(ib,jb+1,kb-1,ns) == 0) THEN ! Land
             tm = tem(ib,jb,kb-1,ns)
             sm = sal(ib,jb,kb-1,ns)
             rm = rho(ib,jb,kb-1,ns)

          ELSE ! No land or both land. 
             tm = 0.5 * ( tem(ib,jb,kb-1,ns) + tem(ib,jb+1,kb-1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb-1,ns) + sal(ib,jb+1,kb-1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb-1,ns) + rho(ib,jb+1,kb-1,ns) )

          ENDIF
         
          
          IF (sm == 0 .AND. tm == 0) THEN ! If both are land.
             temp = tb
             salt = sb
             dens = rb

          ELSE ! No land
             az = kb - z1 - 0.5 ! Distance to trajectory.

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF


       ELSE ! If the trajectory is closer to kb + 1

          IF (sal(ib,jb,kb+1,ns) == 0 .AND. tem(ib,jb,kb+1,ns) == 0) THEN ! Land
             tm = tem(ib,jb+1,kb+1,ns)
             sm = sal(ib,jb+1,kb+1,ns)
             rm = rho(ib,jb+1,kb+1,ns)

          ELSEIF (sal(ib,jb+1,kb+1,ns) == 0 .AND. tem(ib,jb+1,kb+1,ns) == 0) THEN ! Land
             tm = tem(ib,jb,kb+1,ns)
             sm = sal(ib,jb,kb+1,ns)
             rm = rho(ib,jb,kb+1,ns)

          ELSE
             tm = 0.5 * ( tem(ib,jb,kb+1,ns) + tem(ib,jb+1,kb+1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb+1,ns) + sal(ib,jb+1,kb+1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb+1,ns) + rho(ib,jb+1,kb+1,ns) ) 
          
          ENDIF

          IF (sm == 0 .AND. tm == 0) THEN
             temp = tb
             salt = sb
             dens = rb

          ELSE
             az = z1 - kb + 0.5

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF

       ENDIF


! ----- CROSSING THE LEFT GRID WALL -----
    ELSEIF (y1 < jb) THEN 


! The box the trajectory crosses, same for both z-cases
       tb = 0.5 * ( tem(ib,jb,kb,ns)   + tem(ib,jb-1,kb,ns)   )
       sb = 0.5 * ( sal(ib,jb,kb,ns)   + sal(ib,jb-1,kb,ns)   )
       rb = 0.5 * ( rho(ib,jb,kb,ns)   + rho(ib,jb-1,kb,ns)   )

       IF ( z1 < kb - 0.5 ) THEN
          
          IF (sal(ib,jb,kb-1,ns) == 0 .AND. tem(ib,jb,kb-1,ns) == 0) THEN ! Land
             tm = tem(ib,jb-1,kb-1,ns)
             sm = sal(ib,jb-1,kb-1,ns)
             rm = rho(ib,jb-1,kb-1,ns)
             
          ELSEIF (sal(ib,jb-1,kb-1,ns) == 0 .AND. tem(ib,jb-1,kb-1,ns) == 0) THEN ! Land
             tm = tem(ib,jb,kb-1,ns)
             sm = sal(ib,jb,kb-1,ns)
             rm = rho(ib,jb,kb-1,ns)

          ELSE ! No land or both land. 
             tm = 0.5 * ( tem(ib,jb,kb-1,ns) + tem(ib,jb-1,kb-1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb-1,ns) + sal(ib,jb-1,kb-1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb-1,ns) + rho(ib,jb-1,kb-1,ns) )

          ENDIF

          
          IF (sm == 0 .AND. tm == 0) THEN ! if both land
             temp = tb
             salt = sb
             dens = rb

          ELSE ! If no land
             az = z1 - kb + 0.5

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF

       ELSE


          IF (sal(ib,jb,kb+1,ns) == 0 .AND. tem(ib,jb,kb+1,ns) == 0) THEN ! Land
             tm = tem(ib,jb-1,kb+1,ns)
             sm = sal(ib,jb-1,kb+1,ns)
             rm = rho(ib,jb-1,kb+1,ns)
             
          ELSEIF (sal(ib,jb-1,kb+1,ns) == 0 .AND. tem(ib,jb-1,kb+1,ns) == 0) THEN ! Land
             tm = tem(ib,jb,kb+1,ns)
             sm = sal(ib,jb,kb+1,ns)
             rm = rho(ib,jb,kb+1,ns)

          ELSE ! No land or both land. 
             tm = 0.5 * ( tem(ib,jb,kb+1,ns) + tem(ib,jb-1,kb+1,ns) )
             sm = 0.5 * ( sal(ib,jb,kb+1,ns) + sal(ib,jb-1,kb+1,ns) )
             rm = 0.5 * ( rho(ib,jb,kb+1,ns) + rho(ib,jb-1,kb+1,ns) )

          ENDIF

          
          IF (sm == 0 .AND. tm == 0) THEN ! if both land
             temp = tb
             salt = sb
             dens = rb

          ELSE ! If no land
             az = z1 - kb + 0.5

             temp = (1. - az) * tb + az * tm
             salt = (1. - az) * sb + az * sm
             dens = (1. - az) * rb + az * rm

          ENDIF


       ENDIF

    ELSE
       PRINT *, 'Error: The trajectory is not in the expected gridbox!'
    ENDIF
  


 ! ---- CROSSING Z -----
 ELSEIF (z1 == INT(z1)) THEN

    IF (z1 == kb) THEN ! Crossing the left upper grid wall

       temp = 0.5*(tem(ib,jb,kb,ns)+tem(ib,jb,kb+1,ns))
       salt = 0.5*(sal(ib,jb,kb,ns)+sal(ib,jb,kb+1,ns))
       dens = 0.5*(rho(ib,jb,kb,ns)+rho(ib,jb,kb+1,ns))

    ELSEIF (z1 < kb) THEN ! Crossing the lower grid wall

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

