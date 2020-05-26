SUBROUTINE readfields

  USE mod_precdef
  USE mod_calendar
  USE mod_param
  USE mod_vel
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  
  USE mod_tempsalt
  USE mod_dens
  USE mod_stat

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
  INTEGER                                    :: t ,i ,j ,k ,kk
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  ! = Variables for converting from S to Z
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_r,Cs_r
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_w,Cs_w
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL*8,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt00
  ! ===   ===   ===
  REAL :: uuu(IMT-1,JMT,KM),vvv(IMT,JMT-1,KM),ttt(IMT,JMT,KM)

   ! Variables to calculate potential density
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: rhozvec, depthzvec, latvec
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: tmpzvec, salzvec
  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt), dzt00(imt,jmt) )
     allocate ( sc_r(km), Cs_r(km) )
     allocate ( sc_w(km), Cs_w(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km,1), dzv(imt,jmt,km,1) )
  end if alloCondDZ
  
     !
   ! Allocate variables for density calculations
   !
   if (readTS .and. .not. allocated(tmpzvec)) then
      allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
   end if

  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_r = 0
  Cs_r = 0
  sc_w = 0
  Cs_w = 0


  call datasetswap
  call updateClock
  
!if(ints==intstart) then
! intpart1=startmon-1
!endif
!intpart1=intpart1+1
 
 

  ! === update the time counting ===
!  intpart1    = mod(ints,12)
!  intpart2    = floor((ints)/24.)
!  dstamp      = 'roms_avg.nc'

  !write (dstamp(11:15),'(I5.5)') & 
       !int(currJDtot) - 731365
!  write (dstamp(11:15),'(I5.5)') & 
!       int(currJDtot) - 714777
  dataprefix  = trim(physDataDir) // RunID
!  ncTpos        = mod(ints-8,12)+1
  ncTpos        = currMon
!  print *, "curr JD:", currJDtot, "read file:", dataprefix

!  print *,'tiden',ncTpos,ints,intstart, currMon


 map3D    = [  1    , 2 ,  3,   4 ]
 start3D  = [ 1,  1,  1,   ncTpos ]

 IF (ints==intstart) THEN
 ierr=NF90_CLOSE(ncid)
 ierr = NF90_OPEN (trim(dataprefix),NF90_NOWRITE,ncid)
 IF (ierr /= 0) THEN
 PRINT*,'hello file', dataprefix
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 END IF

 count3D  = [ IMT-1    , JMT, KM, 1]
 ierr=NF90_INQ_VARID(ncid, ueul_name,varid)
 IF (ierr /= 0) THEN
 PRINT*, 'No U'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 ierr=NF90_GET_VAR(ncid,varid,uuu,start3d,count3d)
 IF (ierr /= 0) THEN
 PRINT*, 'No u'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF

 count3D  = [ IMT    , JMT-1, KM, 1]
 ierr=NF90_INQ_VARID(ncid, veul_name,varid)
 IF (ierr /= 0) THEN
 PRINT*, 'No V'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 ierr=NF90_GET_VAR(ncid,varid,vvv,start3d,count3d)
 IF (ierr /= 0) THEN
 PRINT*, 'No v'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF



  ssh (1:IMT,1:JMT)       = get2DfieldNC(trim(dataprefix) , ssh_name)
  uvel(1:IMT-1,1:JMT,:)=uuu(1:IMT-1,1:JMT,:)
  vvel(1:IMT,1:JMT-1,:)=vvv(1:IMT,1:JMT-1,:)

  where (uvel > 1000)
     uvel = 0
  end where
  where (vvel > 1000)
     vvel = 0
  end where
  where (ssh > 1000)
     ssh = 0
  end where

  hs(:imt,:jmt,2) = ssh(:imt,:jmt)
  

  Cs_w = get1DfieldNC (trim(dataprefix), 'Cs_w')
  sc_w = get1DfieldNC (trim(dataprefix), 's_w')
  Cs_r = get1DfieldNC (trim(dataprefix), 'Cs_r')
  sc_r = get1DfieldNC (trim(dataprefix), 's_rho')
  hc   = getScalarNC (trim(dataprefix), 'hc')


  z_w(:,:,0,2) = depth(:,:)

  do k=1,km
     dzt00 = (hc*sc_r(k) + depth*Cs_r(k)) / (hc + depth)
     z_r(:,:,k,2) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt00(:imt,:)
     dzt00 = (hc*sc_w(k) + depth*Cs_w(k)) / (hc + depth)
     z_w(:,:,k,2) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt00(:imt,:)
     dzt(:,:,k,2) = z_w(:,:,k,2)
  end do
  dzt00 = dzt(:,:,km,2)
  dzt(:,:,1:km-1,2)=dzt(:,:,2:km,2)-dzt(:,:,1:km-1,2)
  dzt(:,:,km,2) = ssh(:imt,:) - dzt00
  dzt(:,:,:,1)=dzt(:,:,:,2)
  dzu(1:imt-1,:,:,1) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
  dzv(:,1:jmt-1,:,1) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5

  do k=1,km
     uflux(:,:,k,2)       = uvel(:imt,:,k) * dzu(:,:,k,1) * dyu(:imt,:)
     vflux(:,1:jmt,k,2)   = vvel(:imt,:,k) * dzv(:,:,k,1) * dxv(:imt,:)
  end do

! do j=JMT-1,1,-1
!   write( *, '(400i1)') (int(0.1* vflux(i,j,1,2)),i=1,IMT-1)
! enddo
!stop 4968

#ifdef tempsalt
 count3D  = [ IMT , JMT, KM, 1]
 ierr=NF90_INQ_VARID(ncid, temp_name,varid)
 IF (ierr /= 0) THEN
 PRINT*, 'No TEMP'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 ierr=NF90_GET_VAR(ncid,varid,ttt,start3d,count3d)
 IF (ierr /= 0) THEN
 PRINT*, 'No temp'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 tem(1:IMT,1:JMT,:,2)=ttt(1:IMT,1:JMT,:)
 
 ierr=NF90_INQ_VARID(ncid, salt_name,varid)
 IF (ierr /= 0) THEN
 PRINT*, 'No TEMP'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 ierr=NF90_GET_VAR(ncid,varid,ttt,start3d,count3d)
 IF (ierr /= 0) THEN
 PRINT*, 'No temp'
 PRINT*,NF90_STRERROR (ierr)
 STOP
 END IF
 sal(1:IMT,1:JMT,:,2)=ttt(1:IMT,1:JMT,:)

! ierr=NF90_INQ_VARID(ncid, 'rho',varid)
! IF (ierr /= 0) THEN
! PRINT*, 'No TEMP'
! PRINT*,NF90_STRERROR (ierr)
! STOP
! END IF
! ierr=NF90_GET_VAR(ncid,varid,ttt,start3d,count3d)
! IF (ierr /= 0) THEN
! PRINT*, 'No temp'
! PRINT*,NF90_STRERROR (ierr)
! STOP
! END IF
! rho(1:IMT,1:JMT,:,2)=ttt(1:IMT,1:JMT,:)
 
 
       ! Calculate potential density
      depthzvec = 0.
      do j=1,jmt
         latvec=-80+1./12.*float(j+subGridJmin-1)
         do i=1,IMT
            tmpzvec = tem(i,j,:,nsp)
            salzvec = sal(i,j,:,nsp)
            call statvd(tmpzvec, salzvec, rhozvec ,km ,depthzvec ,latvec)
            rho(i,j,:,nsp)=rhozvec - 1000.
         end do
      end do





#endif

  return

end subroutine readfields
