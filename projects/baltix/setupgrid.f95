SUBROUTINE setupgrid
! =============================================================
!    ===  Set up the grid ===
! =============================================================
! Subroutine for defining the grid of BaltiX. Run once
! before the loop starts.
! -------------------------------------------------------------

  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile

  IMPLICIT none

  INTEGER                                    :: i ,j ,k ,kk, ip, jp
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: temp2d_int
  
  REAL*4, ALLOCATABLE, DIMENSION(:,:)        :: e1t
  REAL*8, ALLOCATABLE, DIMENSION(:,:)        :: temp2d_doub
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:)      :: temp3d_doub
  
  CHARACTER(LEN=200)                         :: gridFile

! ================================================

  start1D  = [ 1]
  count1D  = [KM]
  start2D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count2D  = [         imt,        jmt ,  1 , 1 ]
  map2D    = [          1 ,          2 ,  3 , 4 ]  
  start3D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count3D  = [         imt,        jmt , KM , 1 ]
  map3D    = [          1 ,          2 ,  3 , 4 ]
  
! =================================================
  
  ALLOCATE( temp2d_doub(IMT,JMT), e1t(IMT,JMT) )
  
  ! === dx for T points ===
  gridFile = trim(topoDataDir)//'mesh_mask_baltix.nc'
  ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 3751
  ierr=NF90_INQ_VARID(ncid,'e1t',varid)
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
      e1t(i,:)=temp2d_doub(i,:) ! [m]
  enddo
  
  ! === dy for T points     ===
  ! === gives dxdy for grid ===
  ierr=NF90_INQ_VARID(ncid,'e2t',varid)
  if(ierr.ne.0) stop 3765
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
      dxdy(i,:) = e1t(i,:) * temp2d_doub(i,:) ! [m2]
  enddo
  
  DEALLOCATE( e1t )
  
  ! === dy for u points ===
  ierr=NF90_INQ_VARID(ncid,'e2u',varid)
  if(ierr.ne.0) stop 3764
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
      dyu(i,:)=temp2d_doub(i,:) ! [m]
  enddo
  
  ! === dx for v points ===
  ierr=NF90_INQ_VARID(ncid,'e1v',varid)
  if(ierr.ne.0) stop 3766
  ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
  if(ierr.ne.0) stop 3799
  do i=1,IMT
      dxv(i,:)=temp2d_doub(i,:) ! [m]
  enddo
  
  DEALLOCATE( temp2d_doub )
  
! =======================================
  
  ALLOCATE( temp2d_int(IMT,JMT), kmu(IMT,JMT), kmv(IMT,JMT) )
  
  ! === Bathymetry ===
  if(ierr.ne.0) stop 5751
  ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! kmt field
  if(ierr.ne.0) stop 3767
  ierr=NF90_GET_VAR(ncid,varid,temp2d_int,start2d,count2d)
  do i=1,IMT
      kmt(i,:)=temp2d_int(i,:)
  enddo
  
  kmu=0 ; kmv=0
  do j=1,jmt
      jp=j+1
      if(jp.eq.jmt+1) jp=jmt
      do i=1,imt
          ip=i+1
          if(ip.eq.IMT+1) ip=1
          kmu(i,j)=min(temp2d_int(i,j),temp2d_int(ip,j),KM)
          kmv(i,j)=min(temp2d_int(i,j),temp2d_int(i,jp),KM)
      enddo
  enddo
  
  DEALLOCATE( temp2d_int ) 
  
  ALLOCATE( temp3d_doub(IMT+2,JMT,KM) )
  
  ! === dz at T points ===
  ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
  if(ierr.ne.0) stop 3763
  ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
  do i=1,IMT
      do j=1,JMT
          if(kmt(i,j).ne.0) dztb(i,j,1)=temp3d_doub(i,j,kmt(i,j))
      enddo
  enddo
  
  DEALLOCATE( temp3d_doub )
  
  ! Done reading grid
  ierr=NF90_CLOSE(ncid)

! ==============================================

  ! === Depth coordinates ===
  zw(0)    = 0.d0
  zw(1:KM) = get1DfieldNC(trim(gridFile) ,'e3t_0')
  do k=1,KM
      kk=KM+1-k
      dz(kk)=zw(k)
      zw(k)=zw(k)+zw(k-1)
  end do
  
  ! dz is independent of x,y except at bottom
  do j=1,JMT
      do i=1,IMT
          if(kmt(i,j).eq.0) then
              dztb(i,j,1)=0.
          endif
      enddo
  enddo
  
  dztb(:,:,2)=dztb(:,:,1)  ! this and the third "time" dimension is probably unnecessary

! ===============================================

end SUBROUTINE setupgrid

