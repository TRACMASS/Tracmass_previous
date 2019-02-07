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
  USE mod_getfile
  USE mod_tempsalt

  IMPLICIT none



!!------------------------------------------------------------------------------
 
 INTEGER                                   ::  i, j, k, n, ii, kk, m, jj, jm, im
 INTEGER, SAVE                             ::  ncidp, ncidt, ncidq     ! [AITOR]
 INTEGER                                   ::  ncTpos1                ! [AITOR]
 INTEGER, PARAMETER                        ::  NY=160 ! size of raw data: JMT+1

 REAL(DP), SAVE                            ::  punit, eunit,  pp, tv, pc, pm, pref, Rd, cp, Lv
 REAL*4, SAVE                              ::  celsius0, qmin
 !REAL*4, DIMENSION(KM)                     ::  aa, bb
 REAL*4, ALLOCATABLE, DIMENSION(:,:)       ::  pxy,  th, zh, uh, vh, qh, ph, dpgh
 REAL*4, ALLOCATABLE, DIMENSION(:,:,:)     ::  td, tw, zg, qxyz,txyz,zxyz,uxyz,vxyz
  
 CHARACTER (LEN=200)                       ::  gridFile, fieldFile, prefix
 LOGICAL around

 REAL(DP)    :: scale_factor, add_offset, vertsummam, vertsummap,emp(IMT,JMT)
 REAL*4, DIMENSION(NY)        :: rlat

 REAL*4, ALLOCATABLE, DIMENSION(:,:)   :: ishort
 REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: ishortxyz


!!------------------------------------------------------------------------------

IF ( .NOT. ALLOCATED(txyz) ) THEN
   ALLOCATE ( zxyz(IMT,NY,KM), uxyz(IMT,NY,KM), vxyz(IMT,NY,KM), qxyz(IMT,NY,KM), txyz(IMT,NY,KM), &
   &          pxy(IMT,NY), th(IMT,NY), zh(IMT,NY), uh(IMT,NY), vh(IMT,NY),     &
   &          qh(IMT,NY), ph(IMT,NY),  dpgh(IMT,NY), ishort(IMT,NY), ishortxyz(IMT,NY,KM) )
END IF

!!------------------------------------------------------------------------------


!! Update the time counter

ihour = ihour + int(ff) * ngcm  

IF (ihour >= 24) THEN   !Forward scheme

   ihour = 0
   iday  = iday + 1
   IF (iday > idmax(imon,iyear)) THEN
      iday = 1
      imon = imon + 1
      ncTpos=0
      IF (imon == 13) THEN
         imon  = 1
         iyear = iyear + 1
         ncTpos1 = 0 ! Correction for the u-v file [AITOR]
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


!! Initialise on the first time step

IF (ints == intstart) THEN
   emp=0.
   hs    = 0.
   uflux = 0.
   vflux = 0.
#ifdef tempsalt
   tem   = 0.
   sal   = 0.
   rho   = 0.
#endif

   punit = 1. ! Scale factor for pressure units  
               ! 1.=Pa, 1.e-2=hPa, 1.e.-3=kPa 
   pref = 100000. * punit ![hPa]
   eunit = 1e-3 ![kJ/kg]
   Rd = 287.05d0   ! Gas constant for dry air
   Lv = 2.5d+6     ! Latent heat for condensation of water vapor [J/kg]
   cp = 1004.d0    ! Specific heat for dry air
   celsius0=273.15

!   qmin=1.e-6 ! A minimum of water mass is required to compute water mass trajectories
!   qmin=1.e-7 ! A minimum of water mass is required to compute water mass trajectories
   qmin=1.e-10 ! A minimum of water mass is required to compute water mass trajectories

   iyear = startYear
   imon  = startMon
   iday  = startDay
   ihour = startHour

END IF

ntime = 1000000 * iyear + 10000 * imon + 100 * iday + ihour

ncTpos  = ncTpos+1
ncTpos1 = ncTpos1+1



!! Swap data sets
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

 prefix = '0000/X000000'

 WRITE (prefix(1:4),'(i4)') iyear
 WRITE (prefix(7:10),'(i4)') iyear

 IF (imon < 10) THEN
   WRITE (prefix(12:12),'(i1)') imon
 ELSE
   WRITE (prefix(11:12),'(i2)') imon
 END IF

 start1D  = [ 1]
 count1D  = [NY]
 start2D  = [  1,  1,   1, ncTpos ]
 count2D  = [imt, NY,   1,      1 ]
 map2D    = [ 1 , 2 ,   3 ,     4 ]
 start3D  = [   1,  1,  1, ncTpos ]
 count3D  = [ imt, NY, KM,      1 ]
 map3D    = [  1 , 2 ,  3,      4 ]


 ! READ orography (geopotential at the surface) [AITOR]
 IF(ncTpos==1) THEN

    ierr=NF90_CLOSE(ncid)
    fieldFile = TRIM(topoDataDir)//'Z0.nc'
    ierr=NF90_OPEN (trim(fieldFile),NF90_NOWRITE,ncid)
    IF (ierr /= 0) THEN
       PRINT*,NF90_STRERROR (ierr)
       STOP
    END IF

    ierr=NF90_INQ_VARID(ncid,'Z',varid)
    IF (ierr /= 0) THEN
        PRINT*, 'No Z'
        PRINT*,NF90_STRERROR (ierr)
        STOP
    END IF

    ierr=NF90_GET_VAR(ncid,varid,ishort,start2d,count2d)
    IF (ierr /= 0) THEN
        PRINT*, 'No Z'
        PRINT*,NF90_STRERROR (ierr)
        STOP
    END IF

    ! asking if there are the scale_factor and add_offset attributes
    ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
    IF (ierr == -43) scale_factor =1.0
    ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
    IF (ierr == -43) add_offset = 0.0

    !zxyz(:,:,1) = float(ishort)* scale_factor+ add_offset
    zxyz(:,:,1) = ishort* scale_factor+ add_offset
 END IF

 ! READ surface pressure
 ! =========================================================================
 WRITE (prefix(6:6),'(A)') 'P'
 fieldFile = TRIM(inDataDir)//TRIM(prefix)//'.nc'

 IF (ncTpos==1) THEN
    ierr=NF90_CLOSE(ncidp)
    ierr = NF90_OPEN (trim(fieldFile),NF90_NOWRITE,ncidp)
    IF (ierr /= 0) THEN
        PRINT*,'hello_p', fieldFile
        PRINT*,NF90_STRERROR (ierr)
        STOP
    END IF
 END IF

 ierr=NF90_INQ_VARID(ncidp,'LNSP',varid)

 IF (ierr /= 0) THEN
    PRINT*, 'No P', ncidp, varid
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ierr=NF90_GET_VAR(ncidp,varid, ishort,start2d,count2d)
 IF (ierr /= 0) THEN
    PRINT*, 'No P, size'
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF
   
 ! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncidp, varid,"scale_factor", scale_factor)
 IF (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncidp, varid,"add_offset", add_offset)
 IF (ierr == -43) add_offset = 0.0
   
 !pxy=float(ishort)* scale_factor+ add_offset
 pxy=ishort* scale_factor+ add_offset  


 ! READ humidity
 ! =========================================================================
 WRITE (prefix(6:6),'(A)') 'Q'
 fieldFile = TRIM(inDataDir)//TRIM(prefix)//'.nc'

 IF (ncTpos==1) THEN
    ierr=NF90_CLOSE(ncidq)
    ierr = NF90_OPEN (trim(fieldFile),NF90_NOWRITE,ncidq)
    IF (ierr /= 0) THEN
        PRINT*,'hello_q', fieldFile
        PRINT*,NF90_STRERROR (ierr)
        STOP
    END IF
 END IF

 ierr=NF90_INQ_VARID(ncidq,'Q',varid)
 IF (ierr /= 0) THEN
    PRINT*, 'No Q' 
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ierr=NF90_GET_VAR(ncidq,varid,ishortxyz,start3d,count3d)
 IF (ierr /= 0) THEN
    PRINT*, 'No Q'
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncidq, varid,"scale_factor", scale_factor)
 IF (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncidq, varid,"add_offset", add_offset)
 IF (ierr == -43) add_offset = 0.0
  
 !qxyz=float(ishortxyz)* scale_factor+ add_offset
 qxyz=ishortxyz* scale_factor+ add_offset
 
#if defined hydro
! A minimum of water mass is required to compute water mass trajectories
do i=1,IMT
 do j=1,NY
  do k=1,KM
   if(qxyz(i,j,k)<qmin) qxyz(i,j,k)=qmin
  enddo
 enddo
enddo
#endif

! READ temperature
! =========================================================================
 WRITE (prefix(6:6),'(A)') 'T'
 fieldFile = TRIM(inDataDir)//TRIM(prefix)//'.nc'

 IF (ncTpos==1) THEN
    ierr=NF90_CLOSE(ncidt)
    ierr = NF90_OPEN (trim(fieldFile),NF90_NOWRITE,ncidt)
    IF (ierr /= 0) THEN
        PRINT*,'hello_t', fieldFile
        PRINT*,NF90_STRERROR (ierr)
        STOP
    END IF
 END IF

 ierr=NF90_INQ_VARID(ncidt,'T',varid)
 IF (ierr /= 0) THEN
      PRINT*, 'No T'
      PRINT*,NF90_STRERROR (ierr)
      STOP
 END IF

 ierr=NF90_GET_VAR(ncidt,varid,ishortxyz,start3d,count3d)
 IF (ierr /= 0) THEN
      PRINT*, 'No T, size'
      PRINT*,NF90_STRERROR (ierr)
      STOP
 END IF

 ! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncidt, varid,"scale_factor", scale_factor)
 IF (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncidt, varid,"add_offset", add_offset)
 IF (ierr == -43) add_offset = 0.0
   
 !txyz=float(ishortxyz)* scale_factor+ add_offset
 txyz=ishortxyz* scale_factor+ add_offset


! READ U & V
! =========================================================================

 start3D  = [   1,  1,  1, ncTpos1 ] ! Correct start3D

 fieldFile = TRIM(inDataDir)//TRIM(prefix(1:4))//"/U_V_SHC1_"//TRIM(prefix(1:4))//"_6h_ml.nc"
 
 IF (ncTpos==1) THEN
    ierr=NF90_CLOSE(ncid)
    ierr = NF90_OPEN (trim(fieldFile),NF90_NOWRITE,ncid)
    IF (ierr /= 0) THEN
        PRINT*,'hello_uv', fieldFile
        PRINT*,NF90_STRERROR (ierr)
        STOP
    END IF
 END IF

 ! U
 ierr=NF90_INQ_VARID(ncid,'u',varid)
 IF (ierr /= 0) THEN
    PRINT*, 'No U'
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
 IF (ierr /= 0) THEN
    PRINT*, 'No u'
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 IF (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 IF (ierr == -43) add_offset = 0.0
   
 uxyz=ishortxyz* scale_factor+ add_offset
 !uxyz=float(ishortxyz)* scale_factor+ add_offset
 
 ! V
 ierr=NF90_INQ_VARID(ncid,'v',varid)
 IF (ierr /= 0) THEN
    PRINT*, 'No v'
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
 IF (ierr /= 0) THEN
    PRINT*, 'No V'
    PRINT*,NF90_STRERROR (ierr)
    STOP
 END IF

 ! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 IF (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 IF (ierr == -43) add_offset = 0.0
 
 vxyz=ishortxyz* scale_factor+ add_offset  
 !vxyz=float(ishortxyz)* scale_factor+ add_offset


 DO j=1,NY
    ph(:,j) = EXP (pxy(:,NY+1-j))   ![Pa]
    zh(:,j) =     zxyz(:,NY+1-j,1)  ![m2/s2]
 ENDDO

 DO k=1,KM
   
   ! Reverse latitude order
   DO j=1,NY
      jj=NY+1-j
      th(:,j) = txyz(:,jj,k)
      qh(:,j) = qxyz(:,jj,k)*1.e3 !  [g/Kg]
      uh(:,j) = uxyz(:,jj,k)
      vh(:,j) = vxyz(:,jj,k)
   END DO
   
   uh(:,1)=0. ; vh(:,1)=0.
   
   dpgh (:,:) = ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*ph(:,:) ) / grav

   DO j=1,JMT
      
      jj = j+1
      jm = j 

      DO i=1,IMT
         
         im=i-1
         IF (im == 0) im=IMT
         
   ! A-grid -> middle of grid cells for temperature humidity, layer thickness and pressure
         tem  (i,j,k,2) = 0.25*(th(i,jj)+th(im,jj)+th(i,jm)+th(im,jm)) ![K]
         sal  (i,j,k,2) = 0.25*(qh(i,jj)+qh(im,jj)+qh(i,jm)+qh(im,jm)) ![g/kg]
 
        !print*, sal(IMT/2,JMT/2,:,2)
        !print*, qh(IM/2-1,JMT/2),qh(IM/2,JMT/2),qh(IM/2-1,JMT/2+1),qh(IM/2,JMT/2+1)        
                
        pp = 0.25*(ph(i,jj)+ph(im,jj)+ph(i,jm)+ph(im,jm)) ![Pa]
         dzt  (i,j,k,2) = ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*pp )*punit / grav
         rho  (i,j,k,2) = 0.5*( aa(k)+aa(k-1) + (bb(k)+bb(k-1))*pp )*punit ![Pa]

#if defined hydro
         dzt  (i,j,k,2) = dzt  (i,j,k,2) * sal  (i,j,k,2)*0.001 ! water mass instead of air mass

        ! Use the A-grid as much as possible
         uflux(i,j,k,2)=0.5*( uh(i ,jj)*qh(i ,jj)*dpgh(i ,jj) + &
            &                     uh(i ,jm)*qh(i ,jm)*dpgh(i ,jm)    ) *dydeg       *0.001
         vflux(i,j,k,2)=0.5*( vh(i ,jj)*qh(i ,jj)*dpgh(i ,jj) + &
            &                     vh(im,jj)*qh(im,jj)*dpgh(im,jj)    ) *dxdeg*csu(j)*0.001
#else
         uflux(i,j,k,2)=0.5*( uh(i,jj)+uh(i ,jm) ) * dydegv(j) / grav * &
         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(i,jm)) )*punit
         vflux(i,j,k,2)=0.5*( vh(i,jj)+vh(im,jj) ) * dxdeg*csu(j) / grav * &
         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(im,jj)) )*punit
#endif
      END DO
   END DO
 END DO

 ! Impose zero velocitites at the poles
 vflux(:  ,0,:,:)=0.
 vflux(:,JMT,:,:)=0.


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

#if defined pottemp

 ! Potential temperature (dry)
 td(:,:,:) = tem(:,:,:,2) * ( pref/rho(:,:,:,2) )**(Rd/cp) ![K]
 ! Potential temperature (wet)
 tw(:,:,:) = td(:,:,:) * Lv/cp * sal(:,:,:,2)/1000. ![K]

 tem(:,:,:,2) = td
 sal(:,:,:,2) = tw

#elif defined energy
 
 ! Dry static energy
 td(:,:,:) = cp * tem(:,:,:,2) + 0.5 * (zg(:,:,1:KM) + zg(:,:,0:KM-1)) ![J/kg]
 ! Moist static energy
 !tw(:,:,:) = td(:,:,:) + Lv * sal(:,:,:,2)/1000. ![J/kg]
 ! Latent heat
 tw(:,:,:) = Lv * sal(:,:,:,2)/1000. ![J/kg]

 tem(:,:,:,2) = td * eunit ![kJ/kg]
 sal(:,:,:,2) = tw * eunit ![kJ/kg]

#else
 tem(:,:,:,2) = tem(:,:,:,2)-celsius0 ![C]
 rho(:,:,:,2) = 0.5/cp * (zg(:,:,1:KM) + zg(:,:,0:KM-1)) ![C]

#endif

!PRINT*,'Z'
!PRINT*, zg(IMT/2,JMT/2,:)
!print*, 'T'
!print*, tem(IMT/2,JMT/2,:,2)
!PRINT*, 'Q'
!print*, sal(IMT/2,JMT/2,:,2)
!PRINT*, 'RHO'
!print*, rho(IMT/2,JMT/2,:,2)
!PRINT*,'UFLUX'
!print*,uflux(IMT/2,JMT/2,:,2)
!PRINT*,'VFLUX'
!print*,vflux(IMT/2,JMT/2,:,2)
!STOP
RETURN

END SUBROUTINE readfields

