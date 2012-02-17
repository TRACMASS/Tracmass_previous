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

  REAL*8,  ALLOCATABLE, DIMENSION(:,:)        :: temp2d_doub

  CHARACTER (len=200)                         :: gridFile
  logical                                    :: around

  ! === Init local variables for the subroutine ===
  INTEGER  :: i ,j ,k , kk

  start1D  = [ 1]
  count1D  = [KM]
  start2D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count2D  = [         imt,        jmt ,  1 , 1 ]
  map2D    = [          1 ,          2 ,  3 , 4 ]  
  start3D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count3D  = [         imt,        jmt , KM , 1 ]
  map3D    = [          1 ,          2 ,  3 , 4 ] 

! === Template for setting up grids. Move the code from readfile.f95
  allocate ( mask(imt,jmt), depth(imt,jmt) )

  
    ! === Open mesh file ===
  gridFile = trim(inDataDir)//'topo/roms_grd.nc'
  inquire(file=gridfile,exist=around)
  if(.not.around) stop 4556

  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3751

  allocate( temp2d_doub(IMT,JMT) )
  
! dx
  ierr=NF90_INQ_VARID(ncid,'pm',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF   
  
  dxv=0.
  do j=1,JMT
   do i=1,IMT
    if(temp2d_doub(i,j).ne.0.) dxv(i,j)=1./temp2d_doub(i,j)
   enddo
  enddo
    
    
! dy
  ierr=NF90_INQ_VARID(ncid,'pn',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF   
  dyu=0.
  do j=1,JMT
   do i=1,IMT
    if(temp2d_doub(i,j).ne.0.) dyu(i,j)=1./temp2d_doub(i,j)
   enddo
  enddo
  
  do i=1,IMT
    dxdy(i,:) = dyu(i,:) * dxv(i,:)
  enddo  
  
! total depth
  ierr=NF90_INQ_VARID(ncid,'h',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF   
  do i=1,IMT
    depth(i,:)=temp2d_doub(i,:)
  enddo
    
! mask
    ierr=NF90_INQ_VARID(ncid,'mask_rho',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF   
  do i=1,IMT
    mask(i,:)=temp2d_doub(i,:)
  enddo

  kmt = KM
  where (mask==0) kmt=0
  
  
    
end SUBROUTINE setupgrid

!  inquire(file=gridfile,exist=around)
!  if(.not.around) stop 4556
!
!
!   ierr = NF90_OPEN (gridfile,NF90_NOWRITE,ncid)
!   IF (ierr /= 0) THEN
!    PRINT*,NF90_STRERROR (ierr)
!     END IF
!     
!   ierr = NF90_INQ_VARID (ncid,'ds_01_21_N',varid)
!   ierr = NF90_GET_VAR (ncid,varid,tem2d,start2d,count2d)
!   IF (ierr /= 0) THEN
!    PRINT*,NF90_STRERROR (ierr)
!    STOP
!   END IF

