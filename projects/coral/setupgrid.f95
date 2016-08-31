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
  !allocate ( mask(imt,jmt), depth(imt,jmt) )  !FC
  allocate ( depth(imt,jmt) )  !FC
  ALLOCATE ( z_r(imt,jmt,km,nst) )   !BJ
  ALLOCATE ( z_w(imt,jmt,0:km,nst) ) !BJ

  print*, 'imt=', imt
  print*, 'jmt=', jmt

  !Order is   t  k  i  j 
  map2d    = [3, 4, 1, 2]
  map3d    = [2, 3, 4, 1]

  gridfile =  "/Volumes/P2/ROMS/CORAL/Inputs/Grid/coral_grd.nc"

  ncTpos = 1
  print *, trim(gridfile)
  dxv(:,:) = 1./get2DfieldNC(trim(gridfile), 'pm')
  dyu(:,:) = 1./get2DfieldNC(trim(gridfile), 'pn')

  dxdy = dyu*dxv

  
  depth = get2DfieldNC(trim(gridfile), 'h')
  mask = get2DfieldNC(trim(gridfile), 'mask_rho')
  kmt = 50 

end SUBROUTINE setupgrid
