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
  REAL,          ALLOCATABLE, DIMENSION(:)   :: lat,lon
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
  REAL,          ALLOCATABLE, DIMENSION(:,:) :: dytt,dxtt
  CHARACTER (len=200)                        :: gridfile

  allocate ( lon(imt), lat(jmt) )
  allocate ( dxtt(imt,jmt), dytt(imt,jmt), depth(imt,jmt), mask(imt,jmt)  )
  call coordinat

  gridfile = trim(inDataDir) // 'oscar_vel2009.nc'
   
  start1d  = [  1]
  count1d  = [ imt]
  lon =  get1DfieldNC(trim(gridfile) , 'longitude')
  count1d  = [ jmt]
  lat =  get1DfieldNC(trim(gridfile) , 'latitude')

  
  WHERE ( lon >360 )
     lon = lon -360
  end WHERE
  WHERE ( lon >180 )
     lon = lon-360
  end WHERE
  lon = cshift(lon,481,1)

  do i=1,imt-1
     do j=1,jmt-1
        dxtt(i,j) = l2d( lon(i),lon(i+1),lat(j),lat(j) )
        dytt(i,j) = l2d( lon(i),lon(i),lat(j),lat(j+1) )
     end do
  end do

  do j=1,jmt-1
     dxtt(imt,j) = l2d( lon(imt),lon(1)+360,lat(j),lat(j) )
     dytt(imt,j) = l2d( lon(imt),lon(1)+360,lat(j),lat(j+1) )
  end do


  dxv(1:imt-1,:) = dxtt(1:imt-1,:)/2 + dxtt(2:imt,:)/2
  dyu(:,1:jmt-1) = dytt(:,1:jmt-1)/2 + dytt(:,2:jmt)/2
  dxv(imt,:) = dxtt(imt,:)/2 + dxtt(1,:)/2
  dxdy = dyu * dxv                                                          
  

  start2d  = [ subGridImin, subGridJmin, 1 ,1]
  count2d  = [ subGridImax, subGridJmax, 1 ,1]
  !Use  t=1  i=2  j=3  k=4
  !map2d    = [2, 3, 1, 1]

  kmt = 1
  dz  = 10
  ncTpos = 1
  
  contains 

    function l2d(lon1,lon2,lat1,lat2)
      real :: lon1,lon2,lat1,lat2,l2d
      real :: dlon,dlat,a,c
      dlon = lon2 - lon1
      dlat = lat2 - lat1
      a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
      c = 2 * asin(min(1.0,sqrt(a)))
      l2d = 6367 * c * 1000
    end function l2d
end SUBROUTINE setupgrid
