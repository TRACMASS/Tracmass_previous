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
  !  dzt  - Height of k-cells i 3 dim.  |--- Only one is needed
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
  allocate ( depth(imt,jmt),  mask(imt,jmt)  )
  call coordinat

  gridfile = trim(inDataDir) // '../grid_spec_cm2_ocean_tripolar.nc'
   
  start1d  = [  1]
  count1d  = [ km]

  start2d  = [ subGridImin, subGridJmin, 1 ,1]
  count2d  = [ subGridImax, subGridJmax, 1 ,1]
  !Order is   t  k  i  j
  map2d    = [0, 0, 2, 1]

  start3d  = [1, subGridImin, subGridJmin,  1]
  count3d  = [1, subGridImax, subGridJmax, km]
  !Use  t=1  i=2  j=3  k=4
  map3d    = [2, 3, 4, 1]   

  dxv = int(floor(get2DfieldNC(trim(gridfile) , 'dxu')))
  dyu = int(floor(get2DfieldNC(trim(gridfile) , 'dyu')))
  dxdy = dxv*dyu
  dzt = int(floor(get3DfieldNC(trim(gridfile) , 'dz_t')))

  ! === Load the bathymetry ===
  depth = int(floor(get2DfieldNC(trim(gridfile) , 'depth_t')))
  kmt = int(floor(get2DfieldNC(trim(gridfile) , 'kmt')))

  mask = 1
  where (kmt==0) 
     mask = 0
  end where
  
end SUBROUTINE setupgrid
