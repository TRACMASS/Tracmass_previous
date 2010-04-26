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
  REAL,          ALLOCATABLE, DIMENSION(:)   :: valsz,fort
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:) :: e1v,e1t,e2u,e2t
  CHARACTER (len=200)                        :: gridfile

  allocate ( valsz(km),fort(imt*jmt*2) )
  allocate ( depth(imt,jmt), ang(imt,jmt),  mask(imt,jmt)  )
  call coordinat
   
  start1d  = [  1]
  count1d  = [ km]
  !Order is     t    k            i            j
  start2d  = [  1 ,  1 ,subGridImin ,subGridJmin]
  count2d  = [  1 ,  1 ,subGridImax ,subGridJmax]
  map2d    = [  1 ,  4 ,          3 ,          2]  
  start3d  = [  1 ,  1 ,subGridImin ,subGridJmin]
  count3d  = [  1 , km ,subGridImax ,subGridJmax]
  map3d    = [  3 ,  4 ,          2 ,          1]   
  gridfile = trim(inDataDir) // 'grid.cdf'
  
  ! === Read  and setup horizontal Grid===
  OPEN(41,FILE=trim(inDataDir) // 'fort.41')
  READ(41,*) fort
  dyu(1:imt,:) = reshape(fort(1:imt*jmt),(/ imt, jmt/) )
  dxv(1:imt,:) = reshape(fort(imt*jmt+1:imt*jmt*2),(/ imt, jmt/) )
  dxdy = dyu * dxv
  ang  = get2DfieldNC(trim(gridfile) , 'ang')
    
  ! === Load the bathymetry ===
  depth = int(floor(get2DfieldNC(trim(gridfile) , 'depth')))
  mask = 1
  do i=1,imt
     do j=1,jmt
        if (depth(i,j) .lt. 0) then
           depth(i,j)=0
           mask(i,j)=0
        end if
     end do
  end do
  kmt = depth
    
  !============================
  !===== Read in uv-masks =====
  !============================
  
!!$  ierr=nf90_OPEN('/data/GOM/uvmasks.cdf',nf90_NOWRITE,ncid)
!!$  if(ierr.ne.0) stop 5001
!!$  
!!$  startB(1)=1
!!$  startB(2)=1
!!$  startB(3)=1
!!$  startB(4)=1
!!$  countB(1)=imt
!!$  countB(2)=jmt
!!$  countB(3)=km
!!$  countB(4)=1
!!$  
!!$  ! === Load the masks ===
!!$  ierr=nf90_INQ_VARID(ncid,'u0mask',varid)
!!$  if(ierr.ne.0) stop 5002
!!$  ierr=NF90_GET_VAR (ncid,varid,u0mask)
!!$  if(ierr.ne.0) stop 5003
!!$  u0mask=-(u0mask-1)*9999
!!$  
!!$  ierr=nf90_INQ_VARID(ncid,'u2mask',varid)
!!$  if(ierr.ne.0) stop 5012
!!$  ierr=NF90_GET_VAR (ncid,varid,u2mask)
!!$  if(ierr.ne.0) stop 5013
!!$  u2mask=(u2mask-1)*9999
!!$  
!!$  ierr=nf90_INQ_VARID(ncid,'v0mask',varid)
!!$  if(ierr.ne.0) stop 5022
!!$  ierr=NF90_GET_VAR (ncid,varid,v0mask)
!!$  if(ierr.ne.0) stop 5023
!!$  v0mask=-(v0mask-1)*9999
!!$  
!!$  ierr=nf90_INQ_VARID(ncid,'v2mask',varid)
!!$  if(ierr.ne.0) stop 5032
!!$  ierr=NF90_GET_VAR (ncid,varid,v2mask)
!!$  if(ierr.ne.0) stop 5033
!!$  v2mask=(v2mask-1)*9999
  
end SUBROUTINE setupgrid
