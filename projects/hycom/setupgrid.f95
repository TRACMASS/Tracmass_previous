
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

  INTEGER                                    :: i ,j ,k ,kk

  REAL,          ALLOCATABLE, DIMENSION(:)   :: lat,lon,dz_inv
  REAL,          ALLOCATABLE, DIMENSION(:,:) :: dytt,dxtt
  CHARACTER (len=200)                        :: gridfile

  allocate ( lon(imt), lat(jmt), dz_inv(km) )
  allocate ( dxtt(imt,jmt), dytt(imt,jmt), depth(imt,jmt) )
  call coordinat

  map2d    = [3, 4, 1, 1]
  map3d    = [2, 3, 4, 1]
  ncTpos = 1
  !gridfile = trim(inDataDir) // 'hycom_glb_910_2013121800_t000_uv3z.nc'
  gridfile = trim(inDataDir) // 'hycom_glb_909_2013010100_t000_ts3z.nc'
  start1D  = [subGridImin]
  count1d  = [imt]
  lon =  get1DfieldNC(trim(gridfile) , 'lon')
  start1D  = [subGridJmin]
  count1d  = [jmt]
  lat =  get1DfieldNC(trim(gridfile) , 'lat')
  start1D  = [subGridKmin]
  count1d  = [km]

  do i=1,imt-1
     do j=1,jmt-1
        dxtt(i,j) = l2d( lon(i), lon(i+1), lat(j), lat(j) )
        dytt(i,j) = l2d( lon(i), lon(i), lat(j), lat(j+1) )
     end do
  end do

  do j=1,jmt-1
     dxtt(imt,j) = l2d( lon(imt),lon(1)+360,lat(j),lat(j) )
     dytt(imt,j) = l2d( lon(imt),lon(1)+360,lat(j),lat(j+1) )
  end do

  dxtt(:,jmt) = dxtt(:,jmt-1)
  dytt(:,jmt) = dytt(:,jmt-1)

  dxv(1:imt-1,:) = dxtt(1:imt-1,:)/2 + dxtt(2:imt,:)/2
  dyu(:,1:jmt-1) = dytt(:,1:jmt-1)/2 + dytt(:,2:jmt)/2
  dyu(:,jmt) = dyu(:,jmt-1)
  dxv(imt,:) = dxtt(imt,:)/2 + dxtt(1,:)/2
  dxdy = dyu * dxv                                                          

  dz_inv = get1DfieldNC(trim(gridfile), 'depth')
  if  (km > 1) then
     dz_inv(1:km-1) = dz_inv(2:km)-dz_inv(1:km-1)
     dz_inv(km) = dz_inv(km-1)
  end if
  dz = dz_inv(km:1:-1)
  uvel = get3DfieldNC(trim(gridfile), 'salinity')
  mask = 1
  where (uvel(:,:,1) < 0) mask = 0  
  kmt = km
  end SUBROUTINE setupgrid
