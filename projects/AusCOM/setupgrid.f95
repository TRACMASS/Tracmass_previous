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
  !  dzt  - Height of k-cells i 3 dim.  |--- Only one is needed
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
  INTEGER                                    :: i ,j ,k ,kk, im, jm, ip, jp
!  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask
  CHARACTER (len=200)                        :: gridfile
  REAL tem3d(imt,jmt,km),tem2d(imt,jmt)
  REAL geolon(IMT,JMT),geolat(IMT,JMT),aa
  INTEGER itemp2d(IMT,JMT)
  logical                                    :: around

  allocate ( dzu(imt,jmt,km),dzv(imt,jmt,km),dzt0surf(imt,jmt) )
  allocate ( depth(imt,jmt),  mask(imt,jmt), kmu(imt,jmt), kmv(imt,jmt) )
  
  
  call coordinat

gridfile = trim(inDataDir) // 'grid_auscom.nc'
print *,trim(gridfile),' tday=',tday

  map2D    = [  1, 2, 3, 4 ]  
  map3D    = [  1, 2, 3, 4 ]  
  start1d  = [  1]
  count1d  = [ KM]
  start2d  = [1  ,   1,  1, 1 ]
  count2d  = [IMT, JMT,  1, 1 ]
  start3d  = [1  ,   1,  1, 1 ]
  count3d  = [IMT, JMT, KM, 1]

  inquire(file=gridfile,exist=around)
  if(.not.around) stop 4556


   ierr = NF90_OPEN (gridfile,NF90_NOWRITE,ncid)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
     END IF
     
   ierr = NF90_INQ_VARID (ncid,'ds_01_21_N',varid)
   ierr = NF90_GET_VAR (ncid,varid,tem2d,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF
   
   dxv(1:imt,:)=tem2d
   
   
   ierr = NF90_INQ_VARID (ncid,'ds_10_12_E',varid)
   ierr = NF90_GET_VAR (ncid,varid,tem2d,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF   
   
   dyu(1:imt,:)=tem2d
   dxdy = dxv(1:imt,:)*dyu(1:imt,:)
  
   dx=dxv(IMT/2,JMT/2) ; dy=dyu(IMT/2,JMT/2) ! rough resolution for arclength
  
  ! === Set up cell heights ===
  
  gridfile = trim(inDataDir) // 'monthlydata/ocean_month_20061231.nc'
  print *,trim(gridfile)
  
   ierr = NF90_OPEN (gridfile,NF90_NOWRITE,ncid)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
     END IF
     
   ierr = NF90_INQ_VARID (ncid,'dht',varid)
   ierr = NF90_GET_VAR (ncid,varid,tem3d,start3d,count3d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF
  
  where (tem3d<-10000.) 
   tem3d = 0.
  end where

  do  k=1,km
   dzt(:,:,km-k+1) =  tem3d(:,:,k)
  end do
  dz=0.
  do i=1,IMT
   do j=1,JMT
    do k=1,KM
     dz(k)=max(dz(k),dzt(i,j,k))
    enddo
   enddo
  enddo
!  aa=0.
!  do k=KM,1,-1
!  aa=aa+dzt(IMT/2,JMT/2,k)
!  print *,k,aa,dzt(IMT/2,JMT/2,k)
!  enddo
!  stop 406
  
  do i=1,IMT
   ip=i+1
   if(ip.eq.IMT+1) ip=1
   do j=1,JMT
    jp=j+1
    if(jp.eq.JMT+1) jp=JMT ! there should be an f-fold here i think K.Döös
     do k=1,KM
      dzu(i,j,k) = min(dzt(i,j,k),dzt(ip,j,k))
      dzv(i,j,k) = min(dzt(i,j,k),dzt(i,jp,k))
!      if(k.eq.KM .and. i==1) print *,j, dzt(i,j,k), dzt(ip,j,k), dzu(i,j,k)
     enddo
    enddo
   enddo
! find the total depth
depth=0.
 do k=1,KM
  depth(:,:)=depth(:,:)+dzt(:,:,k)
 enddo

  ! === Load the bathymetry ===
   ierr = NF90_INQ_VARID (ncid,'kmt',varid)
   ierr = NF90_GET_VAR (ncid,varid,tem2d,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF  
   
  kmt=tem2d
  where (kmt<-1000) 
   kmt = 0
  end where


   ierr = NF90_INQ_VARID (ncid,'kmu',varid)
   ierr = NF90_GET_VAR (ncid,varid,tem2d,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF  

  kmu=tem2d
  where (kmu<-1000) 
   kmu = 0
  end where
  
  
   ierr = NF90_INQ_VARID (ncid,'kmv',varid)
   ierr = NF90_GET_VAR (ncid,varid,tem2d,start2d,count2d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF  

  kmv=tem2d
  where (kmv<-1000) 
   kmv = 0
  end where

  
  mask = 1
  where (kmt==0) 
     mask = 0
  end where

! Extra for extracting and testing variables which should be commented 
!   ierr = NF90_INQ_VARID (ncid,'geolon_c',varid)
!   ierr = NF90_GET_VAR (ncid,varid,geolon,start2d,count2d)
!   IF (ierr /= 0) THEN
!    PRINT*,NF90_STRERROR (ierr)
!    STOP
!   END IF  
!   
!   ierr = NF90_INQ_VARID (ncid,'geolat_c',varid)
!   ierr = NF90_GET_VAR (ncid,varid,geolat,start2d,count2d)
!   IF (ierr /= 0) THEN
!    PRINT*,NF90_STRERROR (ierr)
!    STOP
!   END IF  
!
!open(21,file='/Users/doos/data/auscom/topo/longlat',form='unformatted')
!write(21) geolon
!write(21) geolat
!write(21) kmt
!write(21) kmu
!close(21)
!
!do j=JMT,1,-1
! write(*,"(370f13.3)") (geolat(i,j),i=1,10)
! write(*,"(370i1)") (mask(i,j),i=1,IMT/2)
!enddo 
!
!stop 3495
  
end SUBROUTINE setupgrid
