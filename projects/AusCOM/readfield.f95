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

  CHARACTER(len=200), SAVE         :: ncFile, maskfile
  INTEGER                          :: i, im, ip, j, jm, jp, k, kk, ii
  INTEGER                          :: nread, ndates
  INTEGER, SAVE                    :: varid_u,varid_v,varid_ssh,varid_temp,varid_salt,varid_rho
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
  
! === Initialising fields ===
  initFieldcond: if(ints.eq.intstart) then
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
   imon=0 ; iyear=startYear 
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif
ihour=startHour
iday=startDay
imon=startMon
iyear=startYear

else
  
  ! time count
  imon=imon+nff
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
   if(iyear.eq.yearmax+1) iyear=yearmin
  elseif(imon.eq.0) then
   imon=12
   iyear=iyear-1
   if(iyear.eq.yearmin-1) iyear=yearmax
  endif
  
  endif initFieldcond

  
! === Time number ===
ntime=10000*iyear+100*imon+iday

!print *,ints,intstart,ntime

 
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
  count1d  = [ KM]
  start2d  = [1,  1, nread, 1   ]
  count2d  = [IMT, JMT, 1, 1 ]
  start3d  = [1,  1, 1, nread]
  count3d  = [IMT, JMT, KM, 1]

  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556
  
  if(ints==intstart) then
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
   ierr = NF90_INQ_VARID(ncid,'temp',varid_temp) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'salt',varid_salt) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   ierr = NF90_INQ_VARID(ncid,'rho',varid_rho) 
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
   ! Read temp 
   ierr = NF90_GET_VAR(ncid,varid_temp,temp3d_simp,start3d,count3d)
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
   ! Read salt
   ierr = NF90_GET_VAR(ncid,varid_salt,temp3d_simp,start3d,count3d)
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
   do i=1,IMT
    do j=1,JMT   
     do k=1,KM
      if(temp3d_simp(i,j,k) <- 1000.) then
       rho(i,j,km-k+1,2) = 0.
      else
       rho(i,j,km-k+1,2) = temp3d_simp(i,j,k) - 1000.
      endif
     enddo
    enddo
   enddo
#endif



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
     if(kk.le.kmu(i,j)) uflux(i,j,k,2) = 0.5 * (uvel(i,j,kk)+uvel(i,jm,kk)) * dyu(i,j) * du
     if(kk.le.kmv(i,j)) vflux(i,j,k,2) = 0.5 * (vvel(i,j,kk)+vvel(im,j,kk)) * dxv(i,j) * dv
    enddo
   enddo
  enddo


   DEALLOCATE ( temp3d_simp, temp2d_simp )
#ifdef tempsalt
!   DEALLOCATE ( tempb, saltb, rhob, depthb, latb )
#endif

end subroutine readfields
