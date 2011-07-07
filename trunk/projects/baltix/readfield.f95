SUBROUTINE readfields
!!------------------------------------------------------------------------------
!!
!!
!!       SUBROUTINE: readfields
!!
!!          Needs netcdf and gunzip, since data is compressed in .nc.gz format.
!!
!!          Populates the matrices
!!             uflux          -  zonal volume flux
!!             vflux          -  meridional volume flux
!!               dzt          -  depth of each grid box
!!              dztb          -  depth of bottom grid box
!!               tem          -  temperature         (if option tempsalt)
!!               sal          -  salinity            (if option tempsalt)
!!               rho          -  density             (if option tempsalt)
!!
!!
!!       Last change: Joakim Kjellsson, 5 July 2011
!!
!!------------------------------------------------------------------------------
   USE netcdf
   USE mod_param
   USE mod_vel
   USE mod_coord
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_traj
   USE mod_getfile
   USE mod_seed
#ifdef tempsalt
   USE mod_dens
   USE mod_stat
#endif

   IMPLICIT NONE

!!------------------------------------------------------------------------------

   INTEGER                                      :: ji, jj, jk, ik, ntrac,      &
   &                                               jhour, jday, jmon, jyear,   &
   &                                               kbot,ktop, ntempusb,nread
   INTEGER, PARAMETER                           :: IMTG = 619,                 &
                                                   JMTG = 523,                 &
                                                    KMM = 84
   INTEGER, DIMENSION(1:12,1:31,0:21), SAVE     :: ntempus
   
   REAL*4                                       :: dd,dmult,uint,vint,zint
   REAL*4,  ALLOCATABLE, DIMENSION(:,:)         :: temp2d_simp
   REAL*4,  ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_simp
#ifdef tempsalt
   REAL*4,  ALLOCATABLE, DIMENSION(:)           :: tempb, saltb, rhob, &
   &                                               depthb,latb
#endif
   
   CHARACTER (len=200)                          :: fieldFile, dataprefix,      &
   &                                               zfile, rfile
   
   LOGICAL                                      :: around

!!------------------------------------------------------------------------------
!!--------------------------- ALLOCATE DATA MATRICES ---------------------------
!!------------------------------------------------------------------------------

alloCondUVW: IF (.NOT. ALLOCATED (temp2d_simp)) THEN
      ALLOCATE (temp3d_simp(IMT,JMT,KMM), temp2d_simp(IMT,JMT))
#ifdef tempsalt
      ALLOCATE (tempb(KMM), saltb(KMM), rhob(KMM), depthb(KM), latb(KM))
#endif
END IF alloCondUVW

!!------------------------------------------------------------------------------

!! Definition of the slices of data to read
   start1D  = [ 1]
   count1D  = [KM]
   start2D  = [subGridImin ,subGridJmin ,  1 , 1 ]
   count2D  = [         imt,        jmt ,  1 , 1 ]
   map2D    = [          1 ,          2 ,  3 , 4 ]  
   start3D  = [subGridImin ,subGridJmin ,  1 , 1 ]
   count3D  = [         imt,        jmt ,KMM , 1 ]
   map3D    = [          1 ,          2 ,  3 , 4 ]  
    
    
!! Swap between data sets
   hs(:,:,1)      =  hs(:,:,2)
   uflux(:,:,:,1) =  uflux(:,:,:,2)
   vflux(:,:,:,1) =  vflux(:,:,:,2)
   dzt(:,:,:,1)   =  dzt(:,:,:,2)
#ifdef tempsalt 
   tem(:,:,:,1)   =  tem(:,:,:,2)
   sal(:,:,:,1)   =  sal(:,:,:,2)
   rho(:,:,:,1)   =  rho(:,:,:,2)
#endif

!!------------------------------------------------------------------------------
!!--------------------- IF THIS IS THE FIRST TIME STEP -------------------------
!!------------------------------------------------------------------------------

  
   initFieldcond: IF (ints == intstart) THEN
      
      hs     = 0.
      uflux  = 0.
      vflux  = 0.
#ifdef tempsalt
      tem    = 0.
      sal    = 0.
      rho    = 0.
#endif
      ntempus=89 ; ntempusb=0
      
      ! Set date to start date
      ihour=startHour
      iday=startDay
      imon=startMon
      iyear=startYear
      
      ! Fill ntempus which is a matrix containing the prefix number for files
      ! ntempus is 0 at 1 Jan 03:00 and 2918 at 31 Dec 21:00 (except leap years)
      ik=-1
      DO jmon=1,12
         DO jday=1,idmax(jmon,iyear)
            DO jhour=0,21,3
               ntempus(jmon,jday,jhour) = ik
               ik=ik+1
            END DO
         END DO
      END DO      
    
   ELSE
   
      ! Update date
      ihour = ihour + ngcm
      
      IF (ihour == 24) THEN
         iday = iday+1
         ihour=0
      END IF

      IF (iday > idmax(imon,iyear)) THEN
      
         iday = iday-idmax(imon,iyear)
         imon = imon+1
      
         IF (imon == 13) THEN
            imon  = 1
            iyear = iyear+1
         END IF
      
      END IF
   
   END IF initFieldcond

!!------------------------------------------------------------------------------
!!---------------------------- DETERMINE FILE NAME -----------------------------
!!------------------------------------------------------------------------------
   
   ntime=1000000*iyear+10000*imon+100*iday+ihour

   ! === 1 Jan at 00:00 is named as 31 Dec 24:00 the previous year ===
   IF (imon == 1 .AND. iday == 1 .AND. ihour == 0) THEN
      
      IF (idmax (2,iyear) == 29) THEN
         dataprefix='2927_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
      ELSE
         dataprefix='2919_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
      END IF
   
      WRITE (dataprefix(17:20),'(i4)') iyear-1
      WRITE (dataprefix(26:29),'(i4)') iyear-1

   ELSE 
      
      IF (ntempus(imon,iday,ihour) <= 9) THEN
      
         dataprefix=   'x_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
         WRITE (dataprefix(1:1)  ,'(i1)') ntempus(imon,iday,ihour)
         WRITE (dataprefix(14:17),'(i4)') iyear
         WRITE (dataprefix(23:26),'(i4)') iyear
      
      ELSE IF (ntempus(imon,iday,ihour) <= 99) THEN
         
         dataprefix=  'xx_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
         WRITE (dataprefix(1:2),'(i2)') ntempus(imon,iday,ihour)
         WRITE (dataprefix(15:18),'(i4)') iyear
         WRITE (dataprefix(24:27),'(i4)') iyear
      
      ELSE IF (ntempus(imon,iday,ihour) <= 999) THEN
   
         dataprefix= 'xxx_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
         WRITE (dataprefix(1:3),'(i3)') ntempus(imon,iday,ihour)
         WRITE (dataprefix(16:19),'(i4)') iyear
         WRITE (dataprefix(25:28),'(i4)') iyear
      
      ELSE IF (ntempus(imon,iday,ihour) <= 9999) THEN
         
         dataprefix='xxxx_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
         WRITE (dataprefix(1:4),'(i4)') ntempus(imon,iday,ihour)
         WRITE (dataprefix(17:20),'(i4)') iyear
         WRITE (dataprefix(26:29),'(i4)') iyear
      
      ELSE
         
         PRINT*,' This date is not recognized!'
         PRINT*,' year:',iyear
         PRINT*,' month:',imon
         PRINT*,' day:',iday
         PRINT*,' hour:',ihour
         STOP
      
      END IF
   
   END IF

   fieldFile = trim(inDataDir)//trim(dataprefix)
   nread=1
   start2D  = [subGridImin ,subGridJmin ,  1 , nread ]
   start3D  = [subGridImin ,subGridJmin ,  1 , nread ]

!!------------------------------------------------------------------------------

!!------------------------------------------------------------------------------
!!------------------------ UNZIP AND READ IN DATA ------------------------------
!!------------------------------------------------------------------------------

   ! Unzip the T file
   ! T file contains temperature, salinity, and sea surface height
   
   fieldFile = trim(inDataDir)//trim(dataprefix)//'T.nc'
   INQUIRE (file=trim(fieldFile)//'.gz',exist=around)
   IF (.NOT. around) THEN
      PRINT*,'This file is missing:',fieldFile
      STOP
   END IF
   
   zfile = 'gzip -c -d '//trim(fieldFile)//'.gz > '                            &
   &       //trim(outDataDir)//'tmp/'//trim(outDataFile)
   CALL system (zfile)
   
   rfile = trim(outDataDir)//'tmp/'//trim(outDataFile)
   INQUIRE (file=trim(rfile),exist=around)
   IF (.NOT. around) THEN
      PRINT*,'This file is missing:',rfile
   END IF
   
   ! Open T file
   ierr = NF90_OPEN (trim (rfile),NF90_NOWRITE,ncid)
   
   ! Read SSH
   ierr = NF90_INQ_VARID (ncid,'sossheig',varid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ierr = NF90_GET_VAR (ncid,varid,temp2d_simp,start2d,count2d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   hs(:,:,2) = temp2d_simp(:,:)

#ifdef tempsalt 

   ! Read in temperature
   ierr=NF90_INQ_VARID(ncid,'votemper',varid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
 
   ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   DO ji=1,IMT
      DO jj=1,JMT
         DO jk=1,kmt(ji,jj)
            ik = KM+1-jk
            tem (ji,jj,ik,2) = temp3d_simp (ji,jj,jk)
         END DO
      END DO
   END DO

   ! Read in salinity
   ierr = NF90_INQ_VARID(ncid,'vosaline',varid) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   DO ji=1,IMT
      DO jj=1,JMT
         DO jk=1,kmt(ji,jj)
            ik=KM+1-jk
            sal (ji,jj,ik,2) = temp3d_simp (ji,jj,jk)
         END DO
      END DO
   END DO
   
   ierr = NF90_CLOSE (ncid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF


   ! Compute density
   depthb = 0.
   DO jj=1,JMT
      latb = -80 + 0.25 * FLOAT (jj+subGridJmin-1)
      DO ji=1,IMT
         
         DO jk=1,kmt(ji,jj)
            ik=KM+1-jk
            tempb (jk) = tem (ji,jj,ik,2)
            saltb (jk) = sal (ji,jj,ik,2)
         END DO
         
         CALL statvd(tempb, saltb, rhob ,KM ,depthb ,latb)
         
         DO jk=1,kmt(ji,jj)
            ik=KM+1-jk
            rho(ji,jj,ik,2)=rhob(jk)-1000.
         END DO
      
      END DO
   END DO

#endif     

!!------------------------------------------------------------------------------

   dmult = 1. ! amplification of the velocity amplitude by simple multiplication

   ! Open U file

   fieldFile = trim(inDataDir)//trim(dataprefix)//'U.nc'
   INQUIRE (file=trim(fieldFile)//'.gz',exist=around)
   IF (.NOT. around) THEN
      PRINT *,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
      STOP
   END IF
   
   zfile = 'gzip -c -d '//trim(fieldFile)//'.gz > '                            &
   &       //trim(outDataDir)//'tmp/'//trim(outDataFile)
   CALL system(zfile)
   rfile=trim(outDataDir)//'tmp/'//trim(outDataFile)
   INQUIRE (file=trim(rfile),exist=around)
   IF (.NOT. around) THEN
      PRINT*,'This file is missing:',rfile
      STOP
   END IF

   ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ! Read U velocity
   ierr = NF90_INQ_VARID(ncid,'vozocrtx',varid) 
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ierr = NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ierr=NF90_CLOSE(ncid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   DO ji=1,IMT-1
      DO jj=1,JMT
         DO jk=1,kmu(ji,jj)
            ik=KM+1-jk
            dd = dz(ik)
            IF (jk == kmu(ji,jj)) THEN
               dd = MIN (dztb(ji,jj,2),dztb(ji+1,jj,2))
            ENDIF
            IF (kmu(ji,jj) <= 0) THEN
               dd = 0.
            END IF
            
            ! Depth of grid box
            dd = dd * ( zw(kmt(ji,jj)) + (hs(ji,jj,2) + hs(ji+1,jj,2))/2.      &
            &         )/zw(kmt(ji,jj))
            
            uflux(ji,jj,ik,2) = temp3d_simp(ji,jj,jk) * dyu(ji,jj) * dd * dmult
         END DO
      END DO
   END DO

!-------------------------------------------------------------------------------

   ! Open V file
   fieldFile = trim(inDataDir)//trim(dataprefix)//'V.nc'
   INQUIRE (file=trim(fieldFile)//'.gz',exist=around)
   IF (.NOT. around) THEN
      PRINT*,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
      STOP 4555
   END IF
   zfile = 'gzip -c -d '//trim(fieldFile)//'.gz > '                            &
           //trim(outDataDir)//'tmp/'//trim(outDataFile)
   CALL system(zfile)
   rfile=trim(outDataDir)//'tmp/'//trim(outDataFile)
   INQUIRE (file=trim(rfile),exist=around)
   IF (.NOT. around) THEN
      PRINT*,'This file is missing:',rfile
      STOP 4556
   END IF

   ierr = NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ! Read V velocity
   ierr=NF90_INQ_VARID(ncid,'vomecrty',varid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   ierr=NF90_CLOSE(ncid)
   IF (ierr /= 0) THEN
      PRINT*,NF90_STRERROR (ierr)
      STOP
   END IF
   
   DO ji=1,IMT
      DO jj=1,JMT-1
         DO jk=1,kmv(ji,jj)
            ik = KM+1-jk
            dd = dz(ik)
            IF (jk == kmv(ji,jj)) THEN
               dd = MIN (dztb(ji,jj,2),dztb(ji,jj+1,2))
            END IF
            IF (kmv(ji,jj) <= 0) THEN
               dd = 0.
            END IF
         
            dd = dd * ( zw(kmt(ji,jj)) + (hs(ji,jj,2) + hs(ji,jj+1,2))/2.      &
            &         ) / zw(kmt(ji,jj))
            
            vflux (ji,jj,ik,2) = temp3d_simp(ji,jj,jk) * dxv(ji,jj) * dd * dmult
         END DO
      END DO
   END DO

!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------


   ! Compute z-star coordinates
   ! Layer thicknesses dz* = dz (H+ssh)/H 
   DO ji=1,IMT
      DO jj=1,JMT
         DO jk=1,KM
            ik = KM+1-jk
            IF (kmt(ji,jj) == ik) THEN
               dztb(ji,jj,1) = dztb(ji,jj,1) * ( zw(kmt(ji,jj))+hs(ji,jj,2) )  &
               &               / zw(kmt(ji,jj))
               dzt(ji,jj,jk,2) = dztb(ji,jj,1)
            ELSE IF (kmt(ji,jj) /= 0) THEN
               dzt(ji,jj,jk,2) = dz(jk)
               dzt(ji,jj,jk,2) = dzt(ji,jj,jk,2) *                             &
               &                 (zw(kmt(ji,jj)) + hs(ji,jj,2)) / zw(kmt(ji,jj))
            ELSE
               dzt(ji,jj,jk,2) = 0.
            END IF
         END DO
      END DO
   END DO

!!------------------------------------------------------------------------------

#ifdef drifter
   ! Average velocity/transport over surface drifter drogue depth 
   ! to simulate drifter trajectories
   kbot = 79 ; ktop = 80 ! k=78--80 is z=12--18m
   DO ji=1,IMT
      DO jj=1,JMT
         
         uint=0. ; vint=0. ; zint=0.
         IF (ktop == KM) THEN
            zint = hs(ji,jj,2)
         END IF
         
         DO jk=kbot,ktop
            uint = uint+uflux(ji,jj,jk,2) ! integrated transport
            vint = vint+vflux(ji,jj,jk,2)
            zint = zint+dz(jk)            ! total depth of drougued drifter
         END DO
   
         ! weighted transport for each layer
         DO jk=kbot,KM
            IF (jk /= KM) THEN
               uflux(ji,jj,jk,2) = uint*dz(jk)/zint 
               vflux(ji,jj,jk,2) = vint*dz(jk)/zint
            ELSE
               uflux(ji,jj,jk,2) = uint*(hs(ji,jj,2)+dz(jk))/zint
               vflux(ji,jj,jk,2) = vint*(hs(ji,jj,2)+dz(jk))/zint
            END IF
         END DO
      END DO
   END DO

#endif

   DEALLOCATE ( temp3d_simp, temp2d_simp )
#ifdef tempsalt
   DEALLOCATE ( tempb, saltb, rhob, depthb, latb )
#endif

   RETURN


END SUBROUTINE readfields

!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------



