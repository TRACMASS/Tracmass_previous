SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile

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
  INTEGER                                    :: ji,jj,jk
  REAL*8                                     :: rlatt
  
  ! -------------------------------------------------------------
  
  !!
  !! All levels are active in the IFS model
  !!
  kmt = KM
  
  !!
  !! Grid box sizes
  !!
  dx = 360./float(IMT)   
  dy = 180./float(JMT)  ! Horizontal resolution in degrees [deg]
  stlon1=0. 
  stlat1=-90000.
  dxdeg = dx*deg        
  dydeg = dy*deg        ! Horizontal resolution in radians [rad]
  
  !! Cosine at each v-point (C-grid) 
  DO jj=0,JMT
     phi(jj) = -90. + dy * FLOAT(jj)    ! Latitude for each v-point (C-grid)
     csu(jj) = DCOS (phi(jj) * radian) 
  END DO
  
  !! Cosine at each u-point (C-grid)
  DO jj=1,JMT
     rlatt = 0.5 * ( phi(jj) + phi(jj-1) )   ! Latitude at each u-point (C-grid)
     cst(jj) = DCOS ( rlatt * radian )
  END DO
  
  !! Total horizontal area of grid box
  DO ji=1,IMT
     DO jj=1,JMT
        dxdy(ji,jj) = dx * deg * cst(jj) * dy * deg
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
  tmin  =  150.d0
  tmax  = 1350.d0
  smin  =  150.d0
  smax  = 1400.d0
#endif
#ifdef pottemp
  tmin  =  200.d0
  tmax  = 2000.d0
  smin  =  200.d0
  smax  = 2000.d0
#endif
  dr    = (rmax-rmin)/dble(MR-1) ![hPa]
  dtemp = (tmax-tmin)/dble(MR-1) ![K]
  dsalt = (smax-smin)/dble(MR-1) ![g/kg]
  

  
  !!
  !! Read A(k) and B(k) to determine hybrid coordinate levels.
  !! p(k) = A(k) + B(k) * p_surface
  !!
  ALLOCATE ( aa(0:KM), bb(0:KM) )
  
  OPEN (12,FILE=TRIM(inDataDir)//'topo/model_60lev.txt')
99 FORMAT(10x,f12.6,4x,f10.8)
  
  DO jk=0,KM
     READ (12,99) aa(jk),bb(jk)
  END DO
  
  CLOSE (12)
  
  !!
  !! Setting necessary dummy argument
  !!
  zw=9.e10  ! should not be used
  DO jk=1,KM
     dz(jk) = zw(jk)-zw(jk-1)
  END DO


  ! -------------------------------------------------------------
  
END SUBROUTINE setupgrid
