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


  ! = ECCO Grid fields
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)     :: gridDRC ,gridDRF
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridDXC ,gridDXG
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridDYC ,gridDYG
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: hFacW   ,hFacS
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridRAC
  REAL, DIMENSION(2)                        :: ttest1  ,ttest2
  ! === Init local variables for the subroutine ===
  INTEGER                                    :: i ,j ,k ,kk
  CHARACTER (len=200)                        :: gridfile


! === Template for setting up grids. Move the code from readfile.f95
  allocate ( mask(imt,jmt), depth(imt,jmt) )

  !Order is   t  k  i  j 

  start1d  =  1
  count1d  = 42
  gridfile = trim(inDataDir) // '/grid/' // 'DRC.data'
  gridDRC  = get1dfield()
  gridfile = trim(inDataDir) // '/grid/' // 'DRF.data'
  gridDRF  = get1dfield()
  
  start2d  = [   1,  1]
  count2d  = [2160,320]
  start3d  = [   1,  1, 1]
  count3d  = [2160,320,42]
  gridfile = trim(inDataDir) // '/grid/' // 'DXC.data'
  gridDXC  = get2dfield()
  gridfile = trim(inDataDir) // '/grid/' // 'RAC.data'
  gridRAC  = get2dfield()
  gridfile = trim(inDataDir) // '/grid/' // 'DXG.data'
  gridDXG  = get2dfield()
  gridfile = trim(inDataDir) // '/grid/' // 'DYC.data'
  gridDYC  = get2dfield()
  gridfile = trim(inDataDir) // '/grid/' // 'DXG.data'
  gridDYG  = get2dfield() 
  gridfile = trim(inDataDir) // '/grid/' // 'hFacW.data'
  hFacW    = get3dfield()
  gridfile = trim(inDataDir) // '/grid/' // 'hFacS.data'
  hFacW    = get3dfield()

  dxdy     = gridDXC*gridDYC
  dz       = gridDRC(km:1:-1)
  kmt      = sum(ceiling(hFacW),3)













  ncTpos = 1
  print *, trim(gridfile)
  dxv(:-2,:) = get2DfieldNC(trim(gridfile), 'x_rho')
  dyu(:-2,:) = get2DfieldNC(trim(gridfile), 'y_rho')

  dxv(1:imt-1,:) = dxv(2:imt,:)-dxv(1:imt-1,:)
  dyu(:,1:jmt-1) = dyu(:,2:jmt)-dyu(:,1:jmt-1)
  dxv(imt:imt+1,:) = dxv(imt-2:imt-1,:)
  dyu(:,jmt) = dyu(:,jmt-1)
  dxdy = dyu*dxv

  depth = get2DfieldNC(trim(gridfile), 'h')
  mask = get2DfieldNC(trim(gridfile), 'mask_rho')
  kmt = 40 

  where (mask(2:imt,:) == 0) 
     mask(1:imt-1,:) = 0
  end where
  
  where (mask(1:imt-1,:) == 0) 
     mask(2:imt,:) = 0
  end where
  where (mask(:, 2:jmt) == 0) 
     mask(:,1:jmt-1) = 0
  end where
  where (mask(:, 1:jmt-1) == 0) 
     mask(:, 2:jmt) = 0
  end where

  where (mask==0) kmt=0

end SUBROUTINE setupgrid
