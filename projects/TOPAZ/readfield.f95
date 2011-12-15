SUBROUTINE readfields
  
  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile


#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none

  CHARACTER(len=200), SAVE         :: ncFile, maskfile
  INTEGER                          :: i, im, ip, j, jm, jp, k, kk
  INTEGER                          :: year, month, day, nread, ndates
  INTEGER, SAVE                    :: varid_u,varid_v,varid_ssh,varid_sst,varid_sss,varid_rho
  LOGICAL                          :: around
  CHARACTER(len=100)               :: fileName
  REAL*4                           :: du,dv
  REAL*4,  ALLOCATABLE, DIMENSION(:,:)     :: temp2d_simp
  REAL*4,  ALLOCATABLE, DIMENSION(:,:,:)   :: temp3d_simp
    
  
!  REAL temp2d(IMT,JMT)
 
alloCondUVW: IF (.NOT. ALLOCATED (temp2d_simp)) THEN
      ALLOCATE (temp3d_simp(IMT,JMT,KM), temp2d_simp(IMT,JMT))
!#ifdef tempsalt
!      ALLOCATE (tempb(KMM), saltb(KMM), rhob(KMM), depthb(KM), latb(KM))
!#endif
END IF alloCondUVW


  call datasetswap

  ! === NetCDF file and fields ===
  fileName = 'monthlydata/ocean-20071231_month.nc'
  ncFile   = trim(inDataDir)//fileName
  nread=mod(ints,2)+1
  nread=1
!  print *,'ncFile=',ints,nread,ncFile

  map2D    = [  1, 2, 3, 4 ]  
  map3D    = [  1, 2, 3, 4 ]  
  start1d  = [  1]
  count1d  = [ km]
  start2d  = [1,  1, nread, 1   ]
  count2d  = [IMT, JMT, 1, 1 ]
  start3d  = [1,  1, 1, nread]
  count3d  = [IMT, JMT, KM, 1]

  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556
  
  if(ints==1) then
   ierr = NF90_OPEN (ncFile,NF90_NOWRITE,ncid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID (ncid,'eta_u',varid_ssh)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'u',varid_u) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'v',varid_v) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'sst',varid_sst) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'sss',varid_sss) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'rho',varid_sss) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
  endif
  
  
     ! Read SSH
      
   ierr = NF90_GET_VAR (ncid,varid_ssh,temp2d_simp,start2d,count2d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF

  where (temp2d_simp<-100000.) 
   temp2d_simp = 0.
  end where
  
     DO i=1,IMT
      DO j=1,JMT
   		hs(i,j,2) = temp2d_simp(i,j)
      END DO
   END DO

   ! Read U velocity
   ierr = NF90_GET_VAR(ncid,varid_u,uvel,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ! Read V velocity
   ierr = NF90_GET_VAR(ncid,varid_v,vvel,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
#ifdef tempsalt
   ! Read SST 
   ierr = NF90_GET_VAR(ncid,varid_sst,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   where (temp3d_simp<-1000.) 
    temp3d_simp = 0.
   end where
   do k=1,KM
    tem(:,:,km-k+1,2) =  temp3d_simp(:,:,k)
   enddo
   ! Read SSS
   ierr = NF90_GET_VAR(ncid,varid_sss,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   where (temp3d_simp<-1000.) 
    temp3d_simp = 0.
   end where
   do k=1,KM
    sal(:,:,km-k+1,2) =  temp3d_simp(:,:,k)
   enddo
   ! Read rho
   ierr = NF90_GET_VAR(ncid,varid_rho,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   where (temp3d_simp<-1000.) 
    temp3d_simp = 0.
   end where
   do k=1,KM
    rho(:,:,km-k+1,2) =  temp3d_simp(:,:,k)
   enddo
#endif

!   ierr=NF90_CLOSE(ncid)
!   IF (ierr /= 0) THEN
!      PRINT*,NF90_STRERROR (ierr)
!      STOP
!   END IF

  ! === Velocities ===

  where (uvel<-1000.) 
   uvel = 0.
  end where
  where (vvel<-1000.) 
   vvel = 0.
  end where

  do i=1,IMT
   ip=i+1
   im=i-1
   if(i.eq.IMT) ip=1
   if(i.eq.1  ) im=IMT
   do j=1,JMT
    jp=j+1
    jm=j-1
    if(j.eq.JMT) jp=JMT
    if(j.eq.  1) jm=1
    do k=1,km
     kk=km-k+1
     du=dzu(i,j,k)
     dv=dzv(i,j,k)
     if(k.eq.KM) then
      du=du+0.5*(hs(i,j,2)+hs(ip,j,2))
      dv=dv+0.5*(hs(i,j,2)+hs(i,jp,2))
     endif
     uflux(i,j,k,2) = 0.5 * (uvel(i,j,kk)+uvel(im,j,kk)) * dyu(i,j) * du
     vflux(i,j,k,2) = 0.5 * (vvel(i,j,kk)+vvel(i,jm,kk)) * dxv(i,j) * dv
    enddo
   enddo
  enddo


   DEALLOCATE ( temp3d_simp, temp2d_simp )
#ifdef tempsalt
!   DEALLOCATE ( tempb, saltb, rhob, depthb, latb )
#endif


  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  

contains
  
  
  subroutine datasetswap
    hs(:,:,1)      = hs(:,:,2)
    uflux(:,:,:,1) = uflux(:,:,:,2)
    vflux(:,:,:,1) = vflux(:,:,:,2)
#ifdef explicit_w
    wflux(:,:,:,1) = wflux(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)   = tem(:,:,:,2)
    sal(:,:,:,1)   = sal(:,:,:,2)
    rho(:,:,:,1)   = rho(:,:,:,2)
#endif
  end subroutine datasetswap


end subroutine readfields
