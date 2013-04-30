SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_traj
  USE mod_getfile
  use mod_seed

#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1,hr
  INTEGER                                    :: yr2 ,mn2 ,dy2
  
  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  ! = Variables for converting from S to Z
  REAL,       ALLOCATABLE, DIMENSION(:)      :: sc_r,Cs_r
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL*4, ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  REAL*4, ALLOCATABLE, DIMENSION(:,:)    :: temp2d_simp
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)  :: temp3d_simp

  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( temp3d_simp(IMT-1,JMT,KM), temp2d_simp(IMT,JMT)  )
     allocate ( ssh(imt,jmt), dzt0(imt,jmt) )
     allocate ( sc_r(km), Cs_r(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ

  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_r = 0
  Cs_r = 0

  call datasetswap
  call updateClock
  if(currYear.ge.yearmax+1) currYear=yearmin

  ! === update the time counting ===
  dstamp      = '1993M01.nc'

  write (dstamp(1:4),'(i4)') currYear
  if(currMon.le.9) write (dstamp(6:7),'(i1,i1)') 0,currMon
  if(currMon.gt.9) write (dstamp(6:7),'(i2)') currMon
  gridFile  = trim(inDataDir) // 'fields/roms_avg_Y' // dstamp

print *,gridFile
! u velocity
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'u',varid) 
if(ierr.ne.0) stop 3769
  count3D  = [         imt-1,        jmt , KM , 1 ]
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF 
  do k=1,km
     kk=km+1-k
     uvel(1:imt-1,:,k)=temp3d_simp(1:imt-1,:,kk)
  enddo
  
ierr=NF90_INQ_VARID(ncid,'v',varid) 
if(ierr.ne.0) stop 3769
  count3D  = [         imt,        jmt-1 , KM , 1 ]
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
    PRINT*,NF90_STRERROR (ierr)
    STOP
   END IF      
  do k=1,km
     kk=km+1-k
     vvel(1:imt-1,1:JMT-1,k)=temp3d_simp(1:imt-1,1:JMT-1,kk)
  enddo
!  vvel(1:imt-1,1:JMT-1,:)=temp3d_simp(1:imt-1,1:JMT-1,:)
   
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'sc_r', sc_r)
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'Cs_r', Cs_r)
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'hc', hc)


!ierr=NF90_INQ_VARID(ncid,'sc_w',varid) 
!   IF (ierr /= 0) THEN
!    PRINT*,NF90_STRERROR (ierr)
!    STOP
!   END IF    
!   if(ierr.ne.0) stop 3770
!  count3D  = [         imt,        jmt-1 , KM , 1 ]
!ierr=NF90_GET_VAR(ncid,varid,sc_r,start1d,count1d)
!   IF (ierr /= 0) THEN
!    PRINT*,NF90_STRERROR (ierr)
!    STOP
!   END IF     
   
ierr=NF90_CLOSE(ncid)


!  uvel      = get3DfieldNC(trim(dataprefix) ,   'u')
!  vvel      = get3DfieldNC(trim(dataprefix) ,   'v')
!  ssh       = get2dfieldNC(trim(dataprefix) ,'zeta')
!  hs(:,:,2) = ssh
!
!  ierr = NF90_OPEN(trim(dataprefix) ,NF90_NOWRITE ,ncid)
!  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'sc_w', sc_r)
!  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'Cs_w', Cs_r)
!  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'hc', hc)

  do k=1,km
     dzt0 = (sc_r(k)-Cs_r(k))*hc + Cs_r(k) * depth
     dzt(:,:,k)= dzt0 + ssh*(1.0 + dzt0/depth)
  end do

  dzt0 = dzt(:,:,km)
  dzt(:,:,1:km-1)=dzt(:,:,2:km)-dzt(:,:,1:km-1)
  dzt(:,:,km) = ssh-dzt0

!     print *,k,dzt(100,100,:)
!stop 4967

  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:)*0.5 + dzt(2:imt,:,:)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:)*0.5 + dzt(:,2:jmt,:)*0.5

  do k=1,km
     kk=km+1-k
     uflux(1:imt-1,:,k,2)   = uvel(1:imt-1,:,k) * dzu(1:imt-1,:,k) * dyu(1:imt-1,:)
     vflux(1:imt-1,1:JMT-1,k,2)   = vvel(1:imt-1,1:JMT-1,k) * dzv(1:imt-1,1:JMT-1,k) * &
     dxv(1:imt-1,1:JMT-1)
  end do

  return

end subroutine readfields
