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
  CHARACTER (len=200)                        :: gridfile


! === Template for setting up grids. Move the code from readfile.f95
  allocate ( mask(imt,jmt), depth(imt,jmt) )

  !Order is   t  k  i  j 
  map2d    = [3, 4, 1, 2]
  map3d    = [2, 3, 4, 1]

  gridfile =  trim(inDataDir) // "NEP6_grd_depth_min_15m.nc"

  ncTpos = 1
  print *, trim(gridfile)
  dxv(:imt,:) = get2DfieldNC(trim(gridfile), 'pm')
  dyu(:imt,:) = get2DfieldNC(trim(gridfile), 'pn')
  print *,"Read pm and pn"


!  dxv(1:imt-1,:) = dxv(2:imt,:)-dxv(1:imt-1,:)
!  dyu(:,1:jmt-1) = dyu(:,2:jmt)-dyu(:,1:jmt-1)
!  dxv(imt:imt+1,:) = dxv(imt-2:imt-1,:)
!  dyu(:,jmt) = dyu(:,jmt-1)
  dxv(:imt,:) = 1/dxv(:imt,:)
  dyu(:imt,:) = 1/dyu(:imt,:)
  dxdy = dyu(:imt,:)*dxv(:imt,:)

  
  depth = get2DfieldNC(trim(gridfile), 'h')
  mask = get2DfieldNC(trim(gridfile), 'mask_rho')
  kmt = 50 

  !where (mask(2:imt,:) == 0) 
  !   mask(1:imt-1,:) = 0
  !end where
  
  !where (mask(1:imt-1,:) == 0) 
  !   mask(2:imt,:) = 0
  !end where
  !where (mask(:, 2:jmt) == 0) 
  !   mask(:,1:jmt-1) = 0
  !end where
  !where (mask(:, 1:jmt-1) == 0) 
  !   mask(:, 2:jmt) = 0
  !end where

  !where (mask==0) kmt=0

end SUBROUTINE setupgrid
