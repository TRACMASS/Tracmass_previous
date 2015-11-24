SUBROUTINE readfields
!!------------------------------------------------------------------------------
!!
!!
!!   Subroutine to read the ERA-Interim data fields, and fill the matrices
!!   uflux  -  Zonal mass flux [kg/s]
!!   vflux  -  Meridional mass flux [kg/s]
!!   tem    -  Temperature [K]
!!   sal    -  Specific humidity [g/kg]
!!   rho    -  Pressure [hPa]
!!
!!   If option -Dpottemp is activated in Makefile, tem is potential temperature
!!   [K], and sal is equivalent ("moist") potential temperature [K].
!!   If option -Denergy is activated in Makefile, tem is dry static energy
!!   [kJ/kg] and sal is moist static energy [kJ/kg].
!!
!!
!!------------------------------------------------------------------------------

  USE mod_param
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_dens
  USE mod_stat
  USE mod_tempsalt

  IMPLICIT none

!!------------------------------------------------------------------------------
 
 INTEGER                                            ::  i, j, k, n, ii, kk,    &
 &                                                      im, jj, jm, l
 INTEGER, PARAMETER                                 ::  NY=145
 INTEGER*8, SAVE                                    ::  nlon(NY)

 REAL*4                                             ::  pp, tv, pc, pm, pref,  &
 &                                                      Rd, cp, Lv
 REAL*8, SAVE                                       ::  punit, eunit
 REAL*4, ALLOCATABLE, DIMENSION(:,:)                ::  txy, zxy, uxy, vxy,    &
 &                                                      qxy, pxy,              &
 &                                                      th, zh, uh, vh, qh, ph
 REAL*4, ALLOCATABLE, DIMENSION(:,:,:)              ::  zeta, td, tw, zg
  
 CHARACTER (LEN=200)                                ::  gridFile, fieldFile,   &
 &                                                      string, prefix
 
CHARACTER(LEN=200)                         :: outDataDir, outDataFile
 LOGICAL around


!!------------------------------------------------------------------------------

IF ( .NOT. ALLOCATED(txy) ) THEN
   ALLOCATE ( txy(IMT,NY), zxy(IMT,NY), uxy(IMT,NY), vxy(IMT,NY), qxy(IMT,NY), &
   &          pxy(IMT,NY), th(IMT,NY), zh(IMT,NY), uh(IMT,NY), vh(IMT,NY),     &
   &          qh(IMT,NY), ph(IMT,NY) )
END IF

!!------------------------------------------------------------------------------

!!
!! Update the time counter
!!

ihour = ihour + int(ff) * ngcm  

IF (ihour >= 24) THEN   !Forward scheme

   ihour = 0
   iday  = iday + 1
   IF (iday > idmax(imon,iyear)) THEN
      iday = 1
      imon = imon + 1
      IF (imon == 13) THEN
         imon  = 1
         iyear = iyear + 1
      END IF
   END IF

ELSE IF (ihour < 0) THEN  !Backward scheme
   
   ihour = 18
   iday  = iday - 1
   IF (iday == 0) THEN
      imon = imon - 1
      IF (imon == 0) THEN
         imon  = 12
         iyear = iyear - 1
      END IF
      iday=idmax(imon,iyear)
   END IF

END IF


ntime = 1000000 * iyear + 10000 * imon + 100 * iday + ihour


!!
!! Initialise on the first time step
!!

IF (ints == intstart) THEN
   
   hs    = 0.
   uflux = 0.
   vflux = 0.
#ifdef tempsalt
   tem   = 0.
   sal   = 0.
   rho   = 0.
#endif

   punit = 1.e-2 ! Scale factor for pressure units  
               ! 1.=Pa, 1.e-2=hPa, 1.e.-3=kPa 
   iyear = startYear
   imon  = startMon
   iday  = startDay
   ihour = startHour

END IF


!!
!! Swap data sets
!!

   uflux(:,:,:,1) = uflux(:,:,:,2)
   vflux(:,:,:,1) = vflux(:,:,:,2)
   dzt(:,:,:,1)   = dzt(:,:,:,2)
#ifdef tempsalt 
   tem(:,:,:,1)   = tem(:,:,:,2)
   sal(:,:,:,1)   = sal(:,:,:,2)
   rho(:,:,:,1)   = rho(:,:,:,2)
#endif



!!
!! Construct string to read input files
!!

prefix = '0000/uvtqzp_00000000.0000'

WRITE (prefix(1:4),'(i4)') iyear
WRITE (prefix(13:16),'(i4)') iyear

IF (imon < 10) THEN
   WRITE (prefix(18:18),'(i1)') imon
ELSE
   WRITE (prefix(17:18),'(i2)') imon
END IF

IF (iday < 10) THEN
   WRITE (prefix(20:20),'(i1)') iday
ELSE
   WRITE (prefix(19:20),'(i2)') iday
END IF

IF (ihour < 10) THEN
   WRITE (prefix(23:23),'(i1)') ihour
ELSE
   WRITE (prefix(22:23),'(i2)') ihour
END IF

fieldFile = TRIM(inDataDir)//'era/'//TRIM(prefix)//'.grb'

ntime     = 1000000 * iyear + 10000 * imon + 100 * iday + ihour

INQUIRE (FILE = fieldFile, EXIST = around)
IF (.not. around) THEN
   PRINT*,'ERROR: cannot find: ',fieldFile
   STOP 39467
END IF


!!
!! Read in data from the A-grid using wgrib
!! Set a environment variable WGRIB, e.g. export WGRIB=/home/user/bin/wgrib
!!

!string = '$WGRIB '//TRIM(fieldFile)//' -o '//TRIM(inDataDir)//                 &
!&        TRIM(outDataFile)//'.bin -d all -bin -nh -V > log.txt'

!string = '/Users/doos/Dropbox/bibliotek/wgrib '//TRIM(fieldFile)//' -o '//TRIM(inDataDir)//    &   
!                    &        TRIM(outDataFile)//'.bin -d all -bin -nh -V > log.txt'
                    
string = '/Users/doos/Dropbox/bibliotek/wgrib '//TRIM(fieldFile)//' -o '//TRIM(inDataDir)//    &   
                    &        'tempo.bin -d all -bin -nh -V > log.txt'

!&
!string = '$WGRIB '//TRIM(fieldFile)//' -o '//TRIM(inDataDir)//                 &
!string = 'wgrib '//TRIM(fieldFile)//' -o '//TRIM(inDataDir)//                 &
!string = '/Volumes/hav1/Applications/wgrib/wgrib '//TRIM(fieldFile)//' -o '//TRIM(inDataDir)// &
!print *,'a',string

CALL SYSTEM(string)


!OPEN (14,FORM='UNFORMATTED',FILE=trim(inDataDir)//trim(outDataFile)//'.bin',   &
!&     ACCESS='DIRECT',RECL=IMT*145*4,CONVERT='little_endian')

!print *,'b',trim(inDataDir)//'tempo.bin'
OPEN (14,FORM='UNFORMATTED',FILE=trim(inDataDir)//'tempo.bin',   &
&     ACCESS='DIRECT',RECL=IMT*145*4,CONVERT='little_endian')


!!
!! Read for each model level
!!

kk = 0
DO k=1,KM

   l=k
   kk=kk+1
   READ (14,rec=kk) txy ! read tempeterature in [K]
   kk=kk+1
   READ (14,rec=kk) qxy  ! read specific humidity in [Kg/Kg]
   
   IF (k == 1) THEN
      kk=kk+1
      READ (14,rec=kk) zxy ! read geopotential at the orography in [m2/s2]
      kk=kk+1
      READ (14,rec=kk) pxy ! read the ln surface pressure 
      DO j=1,NY
         ph(:,j) = EXP (pxy(:,NY+1-j))   ![Pa]
         zh(:,j) = zxy(:,NY+1-j)         ![m2/s2]
      END DO
   END IF
   
   kk=kk+1
   READ (14,rec=kk) uxy ! read the zonal velocity in [m/s]
   kk=kk+1
   READ (14,rec=kk) vxy !  read the meridional velocity in [m/s]

   ! Reverse latitude order
   DO j=1,NY
      jj=NY+1-j
      th(:,j) = txy(:,jj)
      qh(:,j) = qxy(:,jj)*1.e3 !  [g/Kg]
      uh(:,j) = uxy(:,jj)
      vh(:,j) = vxy(:,jj)
   END DO
 
   ! In ERA-Interim data, vh is not zero at NP, but the zonal mean is.
   vh(:,NY) = 0.
   
   ! A-grid -> C-grid & store in matrixes
   DO j=1,JMT
      
      jj = j+1
      jm = j
      
      DO i=1,IMT
         
         im=i-1
         IF (im == 0) THEN
            im=IMT
         END IF
         
         tem  (i,j,l,2) = 0.25*(th(i,jj)+th(im,jj)+th(i,jm)+th(im,jm)) ![K]
         sal  (i,j,l,2) = 0.25*(qh(i,jj)+qh(im,jj)+qh(i,jm)+qh(im,jm)) ![g/kg]
         
         pp = 0.25*(ph(i,jj)+ph(im,jj)+ph(i,jm)+ph(im,jm)) ![Pa]
         dzt  (i,j,l,2) = ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*pp )*punit / grav
         
         rho  (i,j,l,2) = 0.5*( aa(k)+aa(k-1) + (bb(k)+bb(k-1))*pp )*punit ![hPa]
                  
         uflux(i,j,k,2)=0.5*( uh(i,jj)+uh(i ,jm) ) * dydeg / grav * &
         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(i,jm)) )*punit
         vflux(i,j,k,2)=0.5*( vh(i,jj)+vh(im,jj) ) * dxdeg*csu(j) / grav * &
         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(im,jj)) )*punit
         
         
      END DO
   END DO

END DO


!!
!! Compute the geopotential, potential temperature, and energies
!!

ALLOCATE (zg(IMT,JMT,0:KM), td(IMT,JMT,KM), tw(IMT,JMT,KM))

zg = 0.
td = 0.
tw = 0.
DO i = 1,IMT
   im = i-1
   IF (im == 0) THEN
      im = IMT
   END IF
   zg(i,:,KM) = 0.25 * (zh(i,1:JMT)+zh(i,2:NY)+zh(im,1:JMT)+zh(im,2:NY))
END DO

pref = 100000. * punit ![hPa]
eunit = 1e-3 ![kJ/kg]

Rd = 287.05d0   ! Gas constant for dry air
Lv = 2.5d+6     ! Latent heat for condensation of water vapor [J/kg]
cp = 1004.d0    ! Specific heat for dry air


DO k = KM,1,-1
   DO j = 1,JMT
      DO i = 1,IMT
         
         im = i-1
         IF (im == 0) THEN
            im = IMT
         END IF
         
         ! Virtual temperature in layer k
         tv = ( 1.0 + 0.61 * sal(i,j,k,2)/1000. ) * tem(i,j,k,2)
         ! Pressure at current interface k
         pc = aa(k) + bb(k) * &
         &         0.25 * (ph(i,j+1)+ph(im,j+1)+ph(i,j)+ph(im,j)) ![Pa]
         ! Pressure at interface k-1
         pm = aa(k-1) + bb(k-1) * &
         &         0.25 * (ph(i,j+1)+ph(im,j+1)+ph(i,j)+ph(im,j)) ![Pa]
         IF (k-1 == 0) THEN
            pm = 10. ![Pa]
         END IF
         ! Geopotential at interface k-1
         zg(i,j,k-1) = zg(i,j,k) + Rd * tv * LOG (pc/pm) ![m2/s2]
      
      END DO
   END DO
ENDDO

#ifdef pottemp
! Potential temperature (dry)
td(:,:,:) = tem(:,:,:,2) * ( pref/rho(:,:,:,2) )**(Rd/cp) ![K]
! Potential temperature (wet)
tw(:,:,:) = td(:,:,:) * Lv/cp * sal(:,:,:,2)/1000. ![K]

tem(:,:,:,2) = td
sal(:,:,:,2) = tw
#endif

#ifdef energy
! Dry static energy
td(:,:,:) = cp * tem(:,:,:,2) + 0.5 * (zg(:,:,1:KM) + zg(:,:,0:KM-1)) ![J/kg]
! Moist static energy
tw(:,:,:) = td(:,:,:) + Lv * sal(:,:,:,2)/1000. ![J/kg]

tem(:,:,:,2) = td * eunit ![kJ/kg]
sal(:,:,:,2) = tw * eunit ![kJ/kg]
#endif

CLOSE (14) 

RETURN

END SUBROUTINE readfields


