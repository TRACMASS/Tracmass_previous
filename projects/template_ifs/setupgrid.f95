SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  USE mod_tempsalt

  IMPLICIT NONE
  
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of the GCM. Run once
  ! before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Area of horizontal cell-walls.
  !  dz   - Height of k-cells in 1 dim. |\
  !  dzt  - Height of k-cells i 3 dim.  |- Only one is needed
  !  kmt  - Number of k-cells from surface to seafloor.
  !
  ! The following might be needed to calculate
  ! dxdy, uflux, and vflux
  !
  !  dzu - Height of each u-gridcell.
  !  dzv - Height of each v-gridcell.
  !  dxu -
  !  dyu -
  ! -------------------------------------------------------------



  ! === Init local variables for the subroutine ===
  INTEGER                                    :: i,j,k
  REAL*8                                     :: rlatt
  
  ! -------------------------------------------------------------
  
  !!
  !! All levels are active in the IFS model
  !!
  kmt = KM
  mask = 1  

  
  !!
  !! Grid box sizes
  !!
  dx = 360./float(IMT)   
  dy = 180./float(JMT)  ! Horizontal resolution in degrees [deg]
!  print *,'dy=',dy
  stlon1=0. 
  stlat1=-90000.
  dxdeg = dx*deg        
  dydeg = dy*deg        ! Horizontal resolution in radians [rad]
  
  !! Cosine at each v-point (C-grid) 
  DO j=0,JMT
     phi(j) = -90. + dy * FLOAT(j)    ! Latitude for each v-point (C-grid)
     csu(j) = DCOS (phi(j) * radian) 
!     print *,jj,phi(jj),csu(jj)
  END DO
  
  !! Cosine at each u-point (C-grid)
  DO j=1,JMT
     rlatt = 0.5 * ( phi(j) + phi(j-1) )   ! Latitude at each u-point (C-grid)
     cst(j) = DCOS ( rlatt * radian )
  END DO
  
  !! Total horizontal area of grid box
  DO i=1,IMT
     DO j=1,JMT
        dxdy(i,j) = dx * deg * cst(j) * dy * deg
     END DO
  END DO
  
  ! -------------------------------------------------------------
  
  !!
  !! Some atmosphere constants
  !!
  
  
  !!
  !! Min/max for pressure, temperature, specific humidity
  !!
  rmin  =    0.d0 ![hPa]
  rmax  = 1100.d0 ![hPa]
  tmin  =  173.d0 ![K]
  tmax  =  323.d0 ![K]
  smin  =    0.d0 ![g/kg]
  smax  =   25.d0 ![g/kg]
#ifdef energy
  tmin  =  240.d0
  tmax  =  400.d0
  smin  =  240.d0
  smax  =  450.d0
#endif
#ifdef pottemp
  tmin  =  240.d0
  tmax  =  400.d0
  smin  =  240.d0
  smax  =  450.d0
#endif
  dr    = (rmax-rmin)/dble(MR-1) ![hPa]
  dtemp = (tmax-tmin)/dble(MR-1) ![K]
  dsalt = (smax-smin)/dble(MR-1) ![g/kg]
  

  
  !!
  !! Read A(k) and B(k) to determine hybrid coordinate levels.
  !! p(k) = A(k) + B(k) * p_surface
  !!
  ALLOCATE ( aa(0:KM), bb(0:KM) )
  
!  OPEN (12,FILE=TRIM(inDataDir)//'topo/model_60lev.txt')
!  OPEN (12,FILE='/Users/doos/data/ifs/topo/model_60lev.txt')
  OPEN (12,FILE='/Users/doos/Dropbox/data_cylinder/ifs/topo/model_60lev.txt')
99 FORMAT(10x,f12.6,4x,f10.8)
  
  DO k=0,KM
     READ (12,99) aa(k),bb(k)
  END DO
  
  CLOSE (12)
  
  !!
  !! Setting necessary dummy argument
  !!
!  zw=9.e10  ! should not be used
  DO k=1,KM
!     print *,k, aa(k),bb(k)
!     dz(jk) = zw(jk)-zw(jk-1)
  END DO
  



  ! -------------------------------------------------------------
  
END SUBROUTINE setupgrid
