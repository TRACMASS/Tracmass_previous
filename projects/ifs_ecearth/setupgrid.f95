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
  REAL*8, DIMENSION(JMT)                     :: rlatv
  REAL*8, DIMENSION(0:JMT)                   :: rlatu
  !REAL*8, DIMENSION(JMT)                     :: dydegv  
  
  !REAL*8, DIMENSION(0:KM),SAVE               :: aa,bb
  ! -------------------------------------------------------------
  
  !!
  !! All levels are active in the IFS model
  !!
  kmt = KM
  mask = 1  

  ! LATITUDE
  ! -------------------------------------------------------------
   OPEN(13,FILE=TRIM(topoDataDir)//'latitude160.txt') ! Update path [AITOR]
  DO j = JMT,1,-1
    READ(13,"(39x,f9.5)") phi(j)
  END DO
  CLOSE(13)
  phi(0)   = -90.
  phi(JMT) = 90.

  ! Dy/Dx length
  ! -------------------------------------------------------------
  ALLOCATE(dydegv(JMT))
  dydegv(1:JMT)  = 0.5*(phi(1:JMT)-phi(0:JMT-1))
  dydegv = dydegv*deg

  dx = 360./float(IMT)
  dxdeg = dx*deg

  ! Cosines definition
  ! -------------------------------------------------------------

  !! Cosine at each v-point (C-grid)
  rlatu = phi
  csu(:) = cos(rlatu(:)*radian)

  !! Cosine at each u-point (C-grid)
  rlatv(:) = 0.5*(phi(1:JMT)+phi(0:JMT-1))
  cst(:) = cos(rlatv(:)*radian)

  !! Total horizontal area of grid box
  DO i=1,IMT
    DO j=1,JMT
        dxdy(i,j) = dx*deg*cst(j)*dydegv(j)
    END DO
  END DO
  
  ! -------------------------------------------------------------

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
  
  ALLOCATE(aa(0:KM),bb(0:KM))
  
  OPEN (12,FILE=TRIM(topoDataDir)//'model_62lev.txt')

  99 FORMAT(10x,f12.6,4x,f10.8)
  
  DO k=0,KM
     READ (12,99) aa(k),bb(k)
  END DO
  
  CLOSE (12)

  ! -------------------------------------------------------------
  
END SUBROUTINE setupgrid
