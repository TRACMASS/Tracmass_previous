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
  REAL,          ALLOCATABLE, DIMENSION(:)   :: valsz
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:) :: e1v,e1t,e2u,e2t
  CHARACTER (len=200)                        :: gridfile

  allocate ( valsz(km) )
  allocate ( depth(imt,jmt), ang(imt,jmt),  mask(imt,jmt)  )
  call coordinat
   
  !Use  t=1  k=2  i=3  j=4
  map2d    = [3, 4, 1, 1]

  start3d  = [  1, subGridImin, subGridJmin,  1]
  count3d  = [  1, subGridImax, subGridJmax, km]


  gridfile = trim(inDataDir) // 'grid.cdf'

  ! === Read  and setup horizontal Grid===
  ncTpos = 1
  dyu(1:imt,:) = get2DfieldNC(trim(gridfile) , 'h1')
  dxv(1:imt,:) = get2DfieldNC(trim(gridfile) , 'h2')
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
