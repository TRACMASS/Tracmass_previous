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
!!   if option -Dpottemp is activated in Makefile, tem is potential temperature
!!   [K], and sal is equivalent ("moist") potential temperature [K].
!!   if option -Denergy is activated in Makefile, tem is dry static energy
!!   [kJ/kg] and sal is moist static energy [kJ/kg].
!!
!!
!!------------------------------------------------------------------------------

  USE mod_param
  USE mod_calendar
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
 INTEGER, PARAMETER                        ::  NY=241 ! size of raw data: JMT+1

 REAL(DP), SAVE                            ::  punit, eunit,  pps, tv, pc, pm, pref, Rd, cp, Lv
 REAL*4, SAVE                              ::  celsius0, qmin
 REAL*4, ALLOCATABLE, DIMENSION(:,:)       ::  pxy,  sst, th, zh, uh, vh, qh, ph, dpgh
 REAL*4, ALLOCATABLE, DIMENSION(:,:,:)     ::  td, tw, zg, qxyz,txyz,zxyz,uxyz,vxyz
  
 CHARACTER (LEN=200)                       ::  gridFile, fieldFile, fieldFileSST, prefix
 LOGICAL around
 
 REAL(DP)    :: scale_factor, add_offset, vertsummam, vertsummap,emp(IMT,JMT)
 REAL*4, DIMENSION(NY)        :: rlat

 INTEGER*4, ALLOCATABLE, DIMENSION(:,:)   :: ishort
 INTEGER*4, ALLOCATABLE, DIMENSION(:,:,:) :: ishortxyz


!!------------------------------------------------------------------------------

if ( .NOT. ALLOCATED(txyz) ) then
   ALLOCATE ( zxyz(IMT,NY,KM), uxyz(IMT,NY,KM), vxyz(IMT,NY,KM), qxyz(IMT,NY,KM), txyz(IMT,NY,KM), &
   &          pxy(IMT,NY), sst(IMT,NY), th(IMT,NY), zh(IMT,NY), uh(IMT,NY), vh(IMT,NY),     &
   &          qh(IMT,NY), ph(IMT,NY),  dpgh(IMT,NY), ishort(IMT,NY), ishortxyz(IMT,NY,KM) )
endif

!!------------------------------------------------------------------------------

!if(ints/=intstart) call update_calendar

!print *,'hej'


!! Update the time counter
ihour = ihour + nff * ngcm  
if (ihour >= 24) then   !Forward scheme
 ihour = 0
 iday  = iday + 1
 if (iday > idmax(imon,iyear)) then
  iday = 1
  imon = imon + 1
  ncTpos=0
  if (imon == 13) then
   imon  = 1
   iyear = iyear + 1
  endif
 endif
elseif (ihour < 0) then  !Backward scheme
  ihour = 18
  iday  = iday - 1
  if (iday == 0) then
   imon = imon - 1
   if (imon == 0) then
    imon  = 12
    iyear = iyear - 1
   endif
   iday=idmax(imon,iyear)
   ncTpos=4*iday
  endif
endif


!! Initialise on the first time step

if (ints == intstart) then
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
   qmin=1.e-7     ! A minimum of water mass is required to compute water mass trajectories

   iyear = startYear
   imon  = startMon
   iday  = startDay
   ihour = startHour


endif

ntime = 1000000 * iyear + 10000 * imon + 100 * iday + ihour
ncTpos=4*(iday-1)+ (ihour)/6 +1
!ncTpos=4*(currDay-1)+ (currHour)/6 +1

!! Swap data sets
!   uflux(:,:,:,1) = uflux(:,:,:,2)
!   vflux(:,:,:,1) = vflux(:,:,:,2)
!     dzt(:,:,:,1)   = dzt(:,:,:,2)
!#ifdef tempsalt 
!   tem(:,:,:,1)   = tem(:,:,:,2)
!   sal(:,:,:,1)   = sal(:,:,:,2)
!   rho(:,:,:,1)   = rho(:,:,:,2)
!#endif
   call datasetswap ! Swap between current and previous step


!!
!! Construct string to read input files
!!

prefix = '0000/'//trim(RunID)//'000000'
WRITE (prefix(1:4),'(i4)') iyear
WRITE (prefix(9:12),'(i4)') iyear
if (imon < 10) then
 WRITE (prefix(14:14),'(i1)') imon
else
 WRITE (prefix(13:14),'(i2)') imon
endif
fieldFile = TRIM(physDataDir)//TRIM(prefix)//trim(fileSuffix)

prefix = '0000/sst000000'
WRITE (prefix(1:4),'(i4)') iyear
WRITE (prefix(9:12),'(i4)') iyear
if (imon < 10) then
 WRITE (prefix(14:14),'(i1)') imon
else
 WRITE (prefix(13:14),'(i2)') imon
endif
fieldFileSST = TRIM(physDataDir)//TRIM(prefix)//trim(fileSuffix)
!print *,'fieldFile ', fieldFile

if (log_level >= 1) then 
   print *,'     Read fieldFile: ', fieldFile
end if 

start1D  = [ 1]
count1D  = [NY]
start2D  = [  1,  1,   1, ncTpos ]
count2D  = [imt, NY,   1,      1 ]
map2D    = [ 1 , 2 ,   3 ,     4 ]  
start3D  = [   1,  1,  1, ncTpos ]
count3D  = [ imt, NY, KM,      1 ]
map3D    = [  1 , 2 ,  3,      4 ] 



ierr=NF90_CLOSE(ncid)
ierr = NF90_OPEN (trim(fieldFileSST),NF90_NOWRITE,ncid)
if (ierr /= 0) then
 PRINT*,NF90_STRERROR (ierr)
 STOP
endif   


!__________________________ Read SST
   ierr=NF90_INQ_VARID(ncid,'sst',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif
   ierr=NF90_GET_VAR(ncid,varid, ishort,start2d,count2d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif
   
! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
! pxy=float(ishort)* scale_factor+ add_offset
 th=float(ishort)* scale_factor+ add_offset

! land values
do j=1,NY 
 do i=1,IMT
 if(ishort(i,j)==-32767) th(i,j)=0.
 enddo
enddo


do j=1,NY
 sst(:,j)=th(:,NY+1-j) 
enddo


!print *,sst(1,1),th(1,NY),ishort(1,NY)
!print *,sst(1,:)

!stop 48956



ierr=NF90_CLOSE(ncid)
ierr = NF90_OPEN (trim(fieldFile),NF90_NOWRITE,ncid)
if (ierr /= 0) then
 PRINT*,NF90_STRERROR (ierr)
 STOP
endif   


!__________________________ Read surface pressure
   ierr=NF90_INQ_VARID(ncid,'lnsp',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif
   ierr=NF90_GET_VAR(ncid,varid, ishort,start2d,count2d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif
   
! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
 pxy=float(ishort)* scale_factor+ add_offset
   


!!_________ Read q
   ierr=NF90_INQ_VARID(ncid,'q',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

   ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
 qxyz=float(ishortxyz)* scale_factor+ add_offset


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



!!_________ Read t
   ierr=NF90_INQ_VARID(ncid,'t',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

   ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
 txyz=float(ishortxyz)* scale_factor+ add_offset



!!_________ Read z
   ierr=NF90_INQ_VARID(ncid,'z',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

   ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
 zxyz=float(ishortxyz)* scale_factor+ add_offset


!!_________ Read u
   ierr=NF90_INQ_VARID(ncid,'u',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

   ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
 uxyz=float(ishortxyz)* scale_factor+ add_offset


!!_________ Read v
   ierr=NF90_INQ_VARID(ncid,'v',varid)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

   ierr=NF90_GET_VAR(ncid,varid,ishortxyz,start3d,count3d)
   if (ierr /= 0) then
      PRINT*,NF90_STRERROR (ierr)
      STOP
   endif

! asking if there are the scale_factor and add_offset attributes
 ierr = nf90_get_att(ncid, varid,"scale_factor", scale_factor)
 if (ierr == -43) scale_factor =1.0
 ierr = nf90_get_att(ncid, varid,"add_offset", add_offset)
 if (ierr == -43) add_offset = 0.0
   
 vxyz=float(ishortxyz)* scale_factor+ add_offset


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
   
   if(k==KM) th=sst
   
   uh(:,1)=0. ; vh(:,1)=0.
   
   dpgh (:,:) = ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*ph(:,:) ) / grav

   DO j=1,JMT
      
      jj = j+1
      jm = j
      
      DO i=1,IMT
         
         im=i-1
         if (im == 0) im=IMT
         
   ! A-grid -> middle of grid cells for temperature humidity, layer thickness and pressure
         tem  (i,j,k,2) = 0.25*(th(i,jj)+th(im,jj)+th(i,jm)+th(im,jm)) ![K]
         if(min(th(i,jj),th(im,jj),th(i,jm),th(im,jm))==0.) then
          tem  (i,j,k,2)=max(th(i,jj),th(im,jj),th(i,jm),th(im,jm))
          tem  (i,j,k,2)=min(tem(i,j,k,2),celsius0+27.)
         endif
          
          
         sal  (i,j,k,2) = 0.25*(qh(i,jj)+qh(im,jj)+qh(i,jm)+qh(im,jm)) ![g/kg]
         
         pps = 0.25*(ph(i,jj)+ph(im,jj)+ph(i,jm)+ph(im,jm)) ![Pa]
         dzt  (i,j,k,2) = ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*pps )*punit / grav
         rho  (i,j,k,2) = 0.5*( aa(k)+aa(k-1) + (bb(k)+bb(k-1))*pps )*punit ![Pa]

#if defined hydro
         dzt  (i,j,k,2) = dzt  (i,j,k,2) * sal  (i,j,k,2)*0.001 ! water mass instead of air mass

!  old version
!         uflux(i,j,k,2)=0.5*( uh(i,jj)*qh(i,jj)+uh(i ,jm)*qh(i ,jm) ) * dydeg / grav * &
!         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(i,jm)) )*punit *0.001
!         vflux(i,j,k,2)=0.5*( vh(i,jj)*qh(i,jj)+vh(im,jj)*qh(im,jj) ) * dxdeg*csu(j) / grav * &
!         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(im,jj)) )*punit *0.001

! Use the A-grid as much as possible
         uflux(i,j,k,2)=ff*0.5*( uh(i ,jj)*qh(i ,jj)*dpgh(i ,jj) + &
        &                     uh(i ,jm)*qh(i ,jm)*dpgh(i ,jm)    ) *dydeg       *0.001
         vflux(i,j,k,2)=ff*0.5*( vh(i ,jj)*qh(i ,jj)*dpgh(i ,jj) + &
        &                     vh(im,jj)*qh(im,jj)*dpgh(im,jj)    ) *dxdeg*csu(j)*0.001

! First the variables on a C-grid then compute the fluxes
!         uflux(i,j,k,2)=0.25*( uh(i ,jj)+uh(i ,jm) ) * (qh(i ,jj)*dpgh(i ,jj) + qh(i ,jm)*dpgh(i ,jm) ) *dydeg       *0.001
 !        vflux(i,j,k,2)=0.25*( vh(i ,jj)+vh(im,jj) ) * (qh(i ,jj)*dpgh(i ,jj) + qh(im,jj)*dpgh(im,jj) ) *dxdeg*csu(j)*0.001         
         
! d      
!         uflux(i,j,k,2)=0.25* ( uh(i,jj)+uh(i ,jm) ) * ( qh(i,jj)+qh(i ,jm) ) * dydeg / grav * &
!         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(i,jm)) )*punit *0.001
!         vflux(i,j,k,2)=0.25* ( vh(i,jj)+vh(im,jj) ) * ( qh(i,jj)*qh(im,jj) ) * dxdeg*csu(j) / grav * &
!         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(im,jj)) )*punit *0.001
#else
         uflux(i,j,k,2)=ff*0.5*( uh(i,jj)+uh(i ,jm) ) * dydeg / grav * &
         &     ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*0.5*(ph(i,jj)+ph(i,jm)) )*punit
         vflux(i,j,k,2)=ff*0.5*( vh(i,jj)+vh(im,jj) ) * dxdeg*csu(j) / grav * &
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
   if (im == 0) then
      im = IMT
   endif
   zg(i,:,KM) = 0.25 * (zh(i,1:JMT)+zh(i,2:NY)+zh(im,1:JMT)+zh(im,2:NY))
END DO


DO k = KM,1,-1
   DO j = 1,JMT
      DO i = 1,IMT
               
         im = i-1
         if (im == 0) then
            im = IMT
         endif
         
         ! Virtual temperature in layer k
         tv = ( 1.0 + 0.61 * sal(i,j,k,2)/1000. ) * tem(i,j,k,2)
         ! Pressure at current interface k
         pc = aa(k) + bb(k) * &
         &         0.25 * (ph(i,j+1)+ph(im,j+1)+ph(i,j)+ph(im,j)) ![Pa]
         ! Pressure at interface k-1
         pm = aa(k-1) + bb(k-1) * &
         &         0.25 * (ph(i,j+1)+ph(im,j+1)+ph(i,j)+ph(im,j)) ![Pa]
         if (k-1 == 0) then
            pm = 10. ![Pa]
         endif
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
!rho(:,:,:,2) = rho(:,:,:,2)*0.01 ![hPa]
rho(:,:,:,2) = 0.5/cp * (zg(:,:,1:KM) + zg(:,:,0:KM-1)) ![C]


#endif


tracers3D(1)%data(:,:,:,nsp) = tem(:,:,:,nsp)
tracers3D(2)%data(:,:,:,nsp) = sal(:,:,:,nsp)
tracers3D(3)%data(:,:,:,nsp) = rho(:,:,:,nsp)

!print *,'hejdu'

RETURN

END SUBROUTINE readfields


