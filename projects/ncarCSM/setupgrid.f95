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
  CHARACTER (len=200)                        :: gridFile
  
  !Order is   t  k  i  j
  map2d    = [3, 4, 1, 1]  
  map3d    = [2, 3, 4, 1]  
  
  gridFile = trim(indatadir) // 'ccsm3t31savedd.pop.h.0002-03-01.nc'
  
  ncTpos = 1
  dz   = get1DfieldNC(trim(gridFile)  ,'dz') / 100.
  dyu  = get2DfieldNC(trim(gridFile) ,'DYU') / 100.
  dxv  = get2DfieldNC(trim(gridFile) ,'DXU') / 100.
  dxdy = dxv * dyu
  kmt  = get2DfieldNC(trim(gridFile) ,'KMT')
  do j=1,jmt
     do i=1,imt
        do k=1,km
           dzt(i,j,k) = dz(km+1-k)
        enddo
     enddo
  enddo

  where (kmt>0)
     mask = 1
  end where


end SUBROUTINE setupgrid
