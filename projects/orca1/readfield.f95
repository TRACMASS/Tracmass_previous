SUBROUTINE readfields

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
  use mod_seed

#ifdef tempsalt
!  USE mod_dens
  USE mod_stat
#endif
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER (len=200)                     :: zfile
  CHARACTER (len=200), SAVE               :: dataprefix,gridFileT,gridFileU,gridFileV,fieldFile

  ! = Loop variables
  INTEGER                                 :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
  INTEGER                                 :: kbot,ktop
  INTEGER, SAVE                           :: nread,npremier,ndernier
  INTEGER, SAVE                           :: ncidt,ncidu,ncidv,varidt,varidu,varidv

  ! = Variables used for getfield procedures

 ! INTEGER, PARAMETER :: IMTG=1440,JMTG=1021,KMM=75

  REAL*4, ALLOCATABLE, DIMENSION(:,:)        :: temp2d_simp
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: temp3d_simp

  REAL*4 						             :: dd,dmult,uint,vint,zint
  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob, depthb,latb
#endif

 LOGICAL around
 
 

!---------------------------------------------------------------


 alloCondUVW: if(.not. allocated (temp2d_simp)) then
   allocate ( temp3d_simp(IMT,JMT,KM), temp2d_simp(IMT,JMT)  )
#ifdef tempsalt
   allocate ( tempb(KM), saltb(KM), rhob(KM), depthb(KM), latb(KM))
#endif
   end if alloCondUVW


  

! === Initialising fields ===
initFieldcond: if(ints==intstart) then
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif

 currHour=startHour
 currDay=startDay
 currMon=startMon
 currYear=startYear
 nread=0

 if(nff==1) then
  npremier=1             ; ndernier= fieldsPerFile
  currYear=yearmin
 elseif(nff==-1) then
  npremier=fieldsPerFile ; ndernier=1
  currYear=yearmax
  currMon=12
 endif
 
 
#ifdef seasonal
  nsp=0
  if(ngcm/=365*24/min(12,NST)) stop 7777
#else
  nsp=2 ; nsm=1
#endif

     
else ! i.e. when ints/=intstart

! === Update clockworks ===

#ifdef seasonal
  currMon=currMon+nff*12/min(12,NST)
#else
  currMon=currMon+nff
#endif

 if(ngcm>=730) then ! monthly or coarser time resolution

  if(currMon > 12) then
   currMon=currMon-12
   currYear=currYear+1
   if(currYear==yearmax+1) currYear=yearmin
  elseif(currMon < 1) then
   currMon=currMon+12
   currYear=currYear-1
  endif
     
 else ! For the 5 day or more often GCM data
  
! === Update clockworks ===
  currHour= currHour+nff*ngcm
  
  if(currHour>=24) then
   currHour= currHour-24
   currDay=currDay+1
  elseif(currHour<0) then 
   currHour=currHour+24
   currDay=currDay-1
  endif
  
  if(currDay > idmax(currMon, currYear)) then 
   currDay=currDay-idmax(currMon, currYear)
   currMon=currMon+1
   if(currMon == 13) then
    currMon=1
    currYear=currYear+1
    if(currYear==yearmax+1) currYear=yearmin
   endif
  elseif(currDay <=0) then
   currMon=currMon-1
   if(currMon == 0) then
    currMon=12
    currYear=currYear-1
    if(currYear.eq.yearmin-1) currYear=yearmax
   endif
   currDay=currDay+idmax(currMon, currYear)
  endif
   
 endif

endif initFieldcond

ntime=100*currYear+currMon

#ifdef seasonal
 nsp=nsp+1
 nsm=nsp-1
 if(nsp==NST+1) then
  nsp=1
 elseif(nsp==0) then
  nsp=NST
 endif
 if(nsm==NST+1) then
  nsm=1
 elseif(nsp==0) then
  nsm=NST
 endif
 

 if(abs(ints)>NST) return
#else
  call datasetswap
#endif

nread=CurrMon

 !print *,'nsm=',nsm,' nsp=',nsp,' ints=',ints


! === Find the file for this timestep ===
start2D  = [subGridImin ,subGridJmin ,  1 , nread ]
start3D  = [subGridImin ,subGridJmin ,  1 , nread ]
  
dataprefix='ORCA1-MSP1_MM_xxxx0101_xxxx1231_grid_'
dataprefix='ORCA1-SS21_MM_xxxx0101_xxxx1231_grid_'

write(dataprefix(15:18),'(i4)') currYear
write(dataprefix(24:27),'(i4)') currYear

fieldFile = trim(inDataDir)//'fields/'//trim(dataprefix)

!print *,nread, fieldFile

if(nread==npremier) then

 gridFileT=trim(fieldFile)//'T.nc'   ! SSH + T + S
! print *, trim(gridFileT)
 inquire(file=trim(gridFileT),exist=around)
 if(.not.around) then
  zfile='gunzip -c '//trim(gridFileT)//'.gz > tmpT'
  CALL system(zfile)
  gridFileT='tmpT'
 endif
 ierr=NF90_OPEN(trim(gridFileT),NF90_NOWRITE,ncidt)
 if(ierr.ne.0) stop 5053
endif


ierr=NF90_INQ_VARID(ncidt,'sossheig',varidt)
if(ierr.ne.0) then
 print *,'file not found:',trim(gridFileT)
 stop 3768
endif



ierr=NF90_GET_VAR(ncidt,varidt,temp2d_simp,start2d,count2d)
if(ierr.ne.0) stop 3799
!ierr=NF90_CLOSE(ncid)

do j=1,JMT
 do i=1,IMT+1
  ii=i
  if(ii.eq.IMT+1) ii=1
  hs(i,j,nsp)=temp2d_simp(ii,j)
 enddo
enddo 

do i=4,IMT
 ii=IMT+4-i
 hs(i,JMT+1,nsp)=hs(ii,JMT-3,nsp)  !  north fold 
enddo

#ifdef tempsalt 
! Temperature
!gridFile = trim(fieldFile)//'T.nc'
!ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncidt,'votemper',varidt) 
if(ierr.ne.0) stop 3769
ierr=NF90_GET_VAR(ncidt,varidt,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
!ierr=NF90_CLOSE(ncid)
do i=1,IMT
 do j=1,JMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   tem(i,j,kk,nsp)=temp3d_simp(i,j,k)
  enddo
 enddo
enddo

!print *,'ntime=',ntime,nread,trim(dataprefix),tem(IMT/2,JMT/2,KM,nsp)
! Salinity
!ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!if(ierr.ne.0) stop 5752
ierr=NF90_INQ_VARID(ncidt,'vosaline',varidt) 
if(ierr.ne.0) stop 3770
ierr=NF90_GET_VAR(ncidt,varidt,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
if(nread==ndernier) then
ierr=NF90_CLOSE(ncidt)
endif
do i=1,IMT
 do j=1,JMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   sal(i,j,kk,nsp)=temp3d_simp(i,j,k)
  enddo
 enddo
enddo

depthb=0.
do j=1,JMT
 latb=-80+0.25*float(j+subGridJmin-1)
 do i=1,IMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   tempb(k)=tem(i,j,kk,nsp)
   saltb(k)=sal(i,j,kk,nsp)
  enddo
  call statvd(tempb, saltb, rhob ,KM ,depthb ,latb)
  do k=1,kmt(i,j)
   kk=KM+1-k
   rho(i,j,kk,nsp)=rhob(k)-1000.
  enddo
 enddo
enddo

#endif     

dmult=1.  ! amplification of the velocity amplitude by simple multiplication
if(nread.eq.npremier) then
 ! u velocity
 gridFileU=trim(fieldFile)//'U.nc'
 inquire(file=trim(gridFileU),exist=around)
 if(.not.around) then
  zfile='gunzip -c '//trim(gridFileU)//'.gz > tmpU'
  CALL system(zfile)
  gridFileU='tmpU'
 endif
 ierr=NF90_OPEN(trim(gridFileU),NF90_NOWRITE,ncidu)
 if(ierr.ne.0) stop 5753
endif

ierr=NF90_INQ_VARID(ncidu,'vozocrtx',varidu) 
if(ierr.ne.0) stop 3771
ierr=NF90_GET_VAR(ncidu,varidu,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
if(nread==ndernier) then
ierr=NF90_CLOSE(ncidu)
endif

do i=1,IMT
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,nsp) + hs(i+1,j,nsp))
   if(k.eq.kmu(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,1)
   uflux(i,j,kk,nsp)=temp3d_simp(i,j,k) * dyu(i,j) * dd * dmult
  enddo
 enddo
enddo


! v velocity
if(nread.eq.npremier) then
gridFileV=trim(fieldFile)//'V.nc'
inquire(file=trim(gridFileV),exist=around)
if(.not.around) then
 zfile='gunzip -c '//trim(gridFileV)//'.gz > tmpV'
 CALL system(zfile)
 gridFileV='tmpV'
endif

ierr=NF90_OPEN(trim(gridFileV),NF90_NOWRITE,ncidv)
if(ierr.ne.0) stop 5754
endif
ierr=NF90_INQ_VARID(ncidv,'vomecrty',varidv) ! kmt field
if(ierr.ne.0) stop 3770
ierr=NF90_GET_VAR(ncidv,varidv,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
if(nread==ndernier) then
ierr=NF90_CLOSE(ncidv)
endif

!  north fold 
!do i=4,IMT
! ii=IMT+4-i
! vflux(i,JMT,:,nsp)=-vflux(ii,JMT-3,:,nsp)
!enddo

do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,nsp) + hs(i,j+1,nsp))
   if(k.eq.kmv(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,2)
   vflux(i,j,kk,nsp)=temp3d_simp(i,j,k) * dxv(i,j) * dd * dmult
  enddo
 enddo
enddo



#ifdef drifter
! average velocity/transport over surface drifter drogue depth to simulate drifter trajectories
kbot=59 ; ktop=60 ! number of surface layers to integrat over
do i=1,imt
 do j=1,jmt
 
  uint=0. ; vint=0. ; zint=0.
  if(ktop.eq.KM) zint=hs(i,j,nsp)
  do k=kbot,ktop
   uint=uint+uflux(i,j,k,nsp) ! integrated transport
   vint=vint+vflux(i,j,k,nsp)
   zint=zint+dz(k)          ! total depth of drougued drifter
  enddo
  ! weighted transport for each layer
  do k=kbot,KM
   if(k.ne.KM) then
    uflux(i,j,k,nsp)=uint*dz(k)/zint 
    vflux(i,j,k,nsp)=vint*dz(k)/zint
   else
    uflux(i,j,k,nsp)=uint*(hs(i,j,nsp)+dz(k))/zint
    vflux(i,j,k,nsp)=vint*(hs(i,j,nsp)+dz(k))/zint
   endif
  enddo

 enddo
enddo

#endif


 deallocate ( temp3d_simp, temp2d_simp )
     
#ifdef tempsalt
 deallocate ( tempb, saltb, rhob, depthb, latb )
#endif

  return
  
end subroutine readfields

