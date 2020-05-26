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
 map2d    = [3, 4, 2, 1]
!  map3d    = [2, 3, 4, 1]
!  map3d    = [4, 3, 1, 2]

  gridFile = trim(topoDataDir)//trim(coordFile)

  ncTpos = 1
  print *, trim(gridfile)
  dxv=0. ; dyu=0.
  dxv(1:IMT,1:JMT) = get2DfieldNC(trim(gridfile), dxv_name)
  dyu(1:IMT,1:JMT) = get2DfieldNC(trim(gridfile), dyu_name)
  
  do i=1,imt+2
   do j=1,jmt
    if(dxv(i,j)/=0.) dxv(i,j)=1./dxv(i,j)
    if(dyu(i,j)/=0.) dyu(i,j)=1./dyu(i,j)
   enddo
  enddo

  dxdy = dyu*dxv
  
  depth(1:IMT,1:JMT) = get2DfieldNC(trim(gridfile), 'h')
  mask(1:IMT,1:JMT) = get2DfieldNC(trim(gridfile), 'mask_rho')
  kmt = KM

! do j=JMT,1,-1
!  write( *, '(400i1)') (mask(i,j),i=1,IMT)
! enddo

  

end SUBROUTINE setupgrid
