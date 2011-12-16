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
  INTEGER                          :: i, im, ip, j, jm, jp, k, kk, ii
  INTEGER                          :: nread, ndates
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
  
! initialise  
  if(ints.eq.1) then
   uflux=0. ; vflux=0.
   imon=0 ; iyear=startYear 
  endif
  
  ! time count
  if(imon.eq.0) iyear=startYear
  imon=imon+1
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
   if(iyear.eq.2008) iyear=2006
  endif
 
  ! === NetCDF file and fields ===
  fileName = ''
  write(fileName(1:4), '(i4) ') iyear
  fileName = 'monthlydata/ocean_month_'//trim(fileName)//'1231.nc'
  ncFile   = trim(inDataDir)//fileName
  nread=mod(ints,12)+1
  nread=imon
!  print *,'ncFile=',imon,iyear,ints,nread,ncFile

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
   ierr = NF90_INQ_VARID (ncid,'eta_t',varid_ssh)
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
   ierr = NF90_GET_VAR(ncid,varid_u,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
     ! === Velocities ===

  where (temp3d_simp<-1000.) 
   temp3d_simp = 0.
  end where
  uvel(1:IMT,:,:)=temp3d_simp

   
   ! Read V velocity
   ierr = NF90_GET_VAR(ncid,varid_v,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
  where (temp3d_simp<-1000.) 
   temp3d_simp = 0.
  end where
  vvel(1:IMT,:,:)=temp3d_simp
   
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


  
!  dyu=1000000. ; dxv=1000000. 
  
!  do j=JMT,1,-1
! write(*,"(370i1)") (int(100000.*uvel(i,j,1)),i=1,IMT/2)
!enddo 
!stop 3495

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
    do k=1,KM
     kk=km-k+1
     du=dzu(i,j,k)
     dv=dzv(i,j,k)
     if(k.eq.KM) then
      du=du+0.5*(hs(i,j,2)+hs(ip,j,2))
      dv=dv+0.5*(hs(i,j,2)+hs(i,jp,2))
     endif
!     if(du.lt.0. .and. kk.le.kmu(i,j)) then
!      print *,i,j,k,dzu(i,j,k),du
!      print *,'hs',hs(i,j,2),hs(ip,j,2)
!      print *,'uv',uvel(i,j,kk),uvel(i,jm,kk)
!      write(*,"(370i1)") (kmt(ii,jp),ii=im,ip)
!      write(*,"(370i1)") (kmt(ii,j),ii=im,ip)
!      write(*,"(370i1)") (kmt(ii,jm),ii=im,ip)
!      stop 4967
!     endif
!      if(dv.lt.0. .and. kk.le.kmv(i,j)) then
!      print *,i,j,k,dzv(i,j,k),dv
!      print *,'hs',hs(i,j,2),hs(i,jp,2)
!      print *,'uv',uvel(i,j,kk),uvel(i,jm,kk)
!      write(*,"(370i1)") (kmt(ii,jp),ii=im,ip)
!      write(*,"(370i1)") (kmt(ii,j),ii=im,ip)
!      write(*,"(370i1)") (kmt(ii,jm),ii=im,ip)
!      stop 4968
!     endif
     if(kk.le.kmu(i,j)) uflux(i,j,k,2) = 0.5 * (uvel(i,j,kk)+uvel(i,jm,kk)) * dyu(i,j) * du
     if(kk.le.kmv(i,j)) vflux(i,j,k,2) = 0.5 * (vvel(i,j,kk)+vvel(im,j,kk)) * dxv(i,j) * dv
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
