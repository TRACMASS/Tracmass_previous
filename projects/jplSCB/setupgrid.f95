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


! === Template for setting up grids. Move the code from readfile.f95
  !Order is     t    k            i            j 
  map2d    = [  1 ,  4 ,          2 ,          3]
  map3d    = [2, 3, 4, 1]


!!$  CHARACTER (len=200)                        :: gridFileXY, gridFileZ
!!$  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: kmask
!!$
!!$  alloCondGrid: if ( .not. allocated (kmask) ) then
!!$     allocate ( kmask(IMT+2,JMT,KM) )
!!$  end if alloCondGrid
!!$  
!!$  start1d  = [  1]
!!$  count1d  = [ km]
!!$  !Order is     t    k            i            j
!!$  start2d  = [  1 ,  1 ,subGridImin ,subGridJmin]
!!$  count2d  = [  1 ,  1 ,subGridImax ,subGridJmax]
!!$  map2d    = [  4 ,  3 ,          1 ,          2]  
!!$  start3d  = [  1 ,  1 ,subGridImin ,subGridJmin]
!!$  count3d  = [  1 , km ,subGridImax ,subGridJmax]
!!$  map3d    = [  4 ,  3 ,          2 ,          1]  
!!$  
!!$  gridFileXY = trim(inDataDir)//'grid_cell_xy.nc'
!!$  gridFileZ  = trim(inDataDir)//'grid_cell_z.nc'
!!$  
!!$  dz   = get1DfieldNC(trim(gridFileZ)  ,'dz')  / 100.
!!$  dxv  = get2DfieldNC(trim(gridFileXY) ,'DXU') / 100.
!!$  dyu  = get2DfieldNC(trim(gridFileXY) ,'DYU') / 100.
!!$  dxdy = dxv * dyu
!!$
!!$  dzt = 0
!!$  kmask  = get3DfieldNC(trim(gridFileZ) ,'SALT')
!!$  do j=1,jmt
!!$     do i=1,imt
!!$        do k=1,km
!!$           kk=km+1-k
!!$           if(kmask(i,j,k) .le. 1000.) then
!!$              kmt(i,j)=k
!!$              dzt(i,j,k) = dz(kk)
!!$           end if
!!$        enddo
!!$     enddo
!!$  enddo

end SUBROUTINE setupgrid
