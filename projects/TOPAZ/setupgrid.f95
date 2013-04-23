SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel

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
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
  CHARACTER (len=200)                        :: gridfile

  allocate ( dzu(imt,jmt,km),dzv(imt,jmt,km),dzt0surf(imt,jmt) )
  allocate ( depth(imt,jmt),  mask(imt,jmt) )
  call coordinat

  gridfile = trim(inDataDir) // '../grid_spec_cm2_ocean_tripolar.nc'
   
  start1d  = [  1]
  count1d  = [ km]

  start2d  = [1, subGridImin, subGridJmin, 1 ]
  count2d  = [1, subGridImax, subGridJmax, 1 ]
  start3d  = [1, subGridImin, subGridJmin,  1]
  count3d  = [1, subGridImax, subGridJmax, km]

  ncTpos = 1
  !Use  t=1  i=2  j=3
  map2d    = [2, 3, 1, 1]
  !Use  t=1  i=2  j=3  k=4
  map3d    = [2, 3, 4, 1]   

  dxv = int(floor(get2DfieldNC(trim(gridfile) , 'dxu')))
  dyu = int(floor(get2DfieldNC(trim(gridfile) , 'dyu')))
  dxv(1:imt,:) = cshift (dxv(1:imt,:), shift=79, dim=1)
  dyu(1:imt,:) = cshift (dyu(1:imt,:), shift=79, dim=1)
  dxdy = dxv*dyu

  ! === Set up cell heights ===
  dzt = int(floor(get3DfieldNC(trim(gridfile) , 'dz_t')))
  dzt = cshift (dzt, shift=79, dim=1)
  dzu = int(floor(get3DfieldNC(trim(gridfile) , 'dz_c')))
  dzu = cshift (dzu, shift=79, dim=1)
  do  k=1,km
     dzt(:,:,km-k+1) =  dzt(:,:,k)
     dzu(:,:,km-k+1) =  dzu(:,:,k)
  end do
  where (dzt>100000) 
     dzt = 0
  end where
 where (dzu>100000) 
     dzu = 0
  end where
  dzv = dzu
  dzt0surf = dzu(:,:,km)
  ! === Load the bathymetry ===
  depth = int(floor(get2DfieldNC(trim(gridfile) , 'depth_t')))
  kmt = int(floor(get2DfieldNC(trim(gridfile) , 'kmt')))
  depth = cshift (depth, shift=79, dim=1)
  kmt = cshift (kmt, shift=79, dim=1)

  mask = 1
  where (kmt==0) 
     mask = 0
  end where
  
end SUBROUTINE setupgrid
