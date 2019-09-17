SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_precdef
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
  INTEGER                                     :: i ,j ,k, ip, jp, im, jm
  INTEGER                                     :: kk, ii
  REAL(PP)                                    :: dlon, dlat
  CHARACTER (len=200)                         :: gridFile
  
  allocate ( kmu(imt,jmt), kmv(imt,jmt) )
  
  dlon = 0.25 
  dlat = 0.25
  
  do j = 1, jmt
     do i = 1, imt
        dxv(i,j) = dlon * pi/180. * radius * cos( (-90+dlat*j)*pi/180. )
        dyu(i,j) = dlat * pi/180. * radius 
     end do
  end do
  
  dx   = dxv(imt/2, jmt/2)
  dy   = dyu(imt/2, jmt/2)
  dxdy(1:imt,1:jmt) = dxv(1:imt,1:jmt) * dyu(1:imt,1:jmt)

  dz(:) = 1. !assume 10m thick layer
  
  kmt(:,:) = 1.
  
  ! comment this out first time you run it through the entire time series
   open(unit=111,file='/Users/doos/data/aviso/topo/kmt.bin',form='unformatted')
   read(111) kmt
   close(111) 
  
   kmu=0 ; kmv=0
   do j=1,jmt
      jp=j+1
      if(jp == jmt+1) jp=jmt
      do i=1,imt
         ip=i+1
         if(ip == imt+1) ip=1
         kmu(i,j)=min(kmt(i,j), kmt(ip,j),KM)
         kmv(i,j)=min(kmt(i,j), kmt(i,jp),KM)
      enddo
   enddo

  
end SUBROUTINE setupgrid
