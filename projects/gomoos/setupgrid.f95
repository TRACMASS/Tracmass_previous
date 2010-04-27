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

  IMPLICIT none
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
  INTEGER                                    :: i ,j ,k ,kk
  REAL,          ALLOCATABLE, DIMENSION(:)   :: valsz,fort
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:) :: e1v,e1t,e2u,e2t
  CHARACTER (len=200)                        :: gridfile

  allocate ( valsz(km),fort(imt*jmt*2) )
  allocate ( depth(imt,jmt), ang(imt,jmt),  mask(imt,jmt)  )
  call coordinat
   
  start1d  = [  1]
  count1d  = [ km]


  start2d  = [ subGridImin, subGridJmin, 1 ,1]
  count2d  = [ subGridImax, subGridJmax, 1 ,1]
  !Order is   t  k  i  j
  map2d    = [0, 0, 1, 2]

  start3d  = [  1 ,  1 ,subGridImin ,subGridJmin]
  count3d  = [  1 , km ,subGridImax ,subGridJmax]
  map3d    = [  3 ,  4 ,          2 ,          1]   
  gridfile = trim(inDataDir) // 'grid.cdf'
  
  ! === Read  and setup horizontal Grid===
  OPEN(41,FILE=trim(inDataDir) // 'fort.41')
  READ(41,*) fort
  dyu(1:imt,:) = reshape(fort(1:imt*jmt),(/ imt, jmt/) )
  dxv(1:imt,:) = reshape(fort(imt*jmt+1:imt*jmt*2),(/ imt, jmt/) )
  dxdy = dyu * dxv
  ang  = get2DfieldNC(trim(gridfile) , 'ang')  / 180 * pi
    
  ! === Load the bathymetry ===
  depth = int(floor(get2DfieldNC(trim(gridfile) , 'depth')))
  mask = 1
  where (depth<0) 
     depth = 0
     mask = 0
  end where
  kmt = 22 * mask
  
end SUBROUTINE setupgrid
