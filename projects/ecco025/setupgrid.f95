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

  INTEGER                                    :: i ,j ,k ,kk

  REAL,          ALLOCATABLE, DIMENSION(:)   :: lat,lon,dz_inv
  REAL,          ALLOCATABLE, DIMENSION(:,:) :: dytt,dxtt
  CHARACTER (len=200)                        :: gridfile

  allocate ( lon(imt), lat(jmt), dz_inv(km) )
  allocate ( dxtt(imt,jmt), dytt(imt,jmt), depth(imt,jmt), mask(imt,jmt) )
  call coordinat

  map2d    = [3, 4, 1, 1]
  map3d    = [2, 3, 4, 1]
  ncTpos = 1
  gridfile = trim(inDataDir) // 'UVEL.1440x720x50.20050204.nc'

  count1d  = imt
  lon =  get1DfieldNC(trim(gridfile) , 'LONGITUDE_T')
  count1d  = jmt
  lat =  get1DfieldNC(trim(gridfile) , 'LATITUDE_T')
  count1d  = km

  do i=1,imt-1
     do j=1,jmt-1
        dxtt(i,j) = l2d( lon(i),lon(i+1),lat(j),lat(j) )
        dytt(i,j) = l2d( lon(i),lon(i),lat(j),lat(j+1) )
     end do
  end do

  do j=1,jmt-1
     dxtt(imt,j) = l2d( lon(imt),lon(1)+360,lat(j),lat(j) )
     !dytt(imt,j) = l2d( lon(imt),lon(1)+360,lat(j),lat(j+1) )
  end do

  dxv(1:imt-1,:) = dxtt(1:imt-1,:)/2 + dxtt(2:imt,:)/2
  dyu(:,1:jmt-1) = dytt(:,1:jmt-1)/2 + dytt(:,2:jmt)/2
  dxv(imt,:) = dxtt(imt,:)/2 + dxtt(1,:)/2
  dxdy = dyu * dxv                                                          

  dz_inv = get1DfieldNC(trim(gridfile), 'DEPTH_T')
  dz_inv(1:km-1) = dz_inv(2:km)-dz_inv(1:km-1)
  dz_inv(km) = dz_inv(km-1)
  dz = dz_inv(km:1:-1)
 
  kmt = 50
  mask = 1

end SUBROUTINE setupgrid
