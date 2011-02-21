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
  USE mod_dens
  USE mod_stat
#endif
  IMPLICIT none
  

  ! = Loop variables
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm,ii,jmm,jpp,ntrac, l
  INTEGER                                      :: kbot,ktop
  INTEGER, SAVE                                :: ntempus,ntempusb,nread

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  
    ! = Variables for filename generation
  CHARACTER (len=200)                        :: dataprefix !, dstamp

#ifdef orca1
  INTEGER, PARAMETER :: IMTG=???,JMTG=???
#elif orca025
  INTEGER, PARAMETER :: IMTG=1440,JMTG=1021,KMM=64
#elif orca025l75h6
  INTEGER, PARAMETER :: IMTG=1440,JMTG=1021,KMM=75
!  INTEGER, PARAMETER :: IMTG=1440,JMTG=1021,KMM=13
#endif

  REAL*4, ALLOCATABLE, DIMENSION(:,:)         :: ssh, temp2d_simp
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_simp
  REAL*8, ALLOCATABLE, DIMENSION(:,:)         :: temp2d_doub
  
  REAL*4,  SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: botbox !,rhom
  REAL*4,  SAVE, ALLOCATABLE, DIMENSION(:,:)   :: e1t,e2t !,rhom
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: kmu,kmv
#ifdef initxyt
  INTEGER, PARAMETER :: NTID=73
!#ifdef orca1
!  INTEGER, PARAMETER :: IJKMAX2=254
!#else
!  INTEGER, PARAMETER :: IJKMAX2=? ! for distmax=0.05 and 32 days
!  INTEGER, PARAMETER :: IJKMAX2=? ! for distmax=0.10 and 32 days
  INTEGER, PARAMETER :: IJKMAX2=7392 ! for distmax=0.25 and 32 days

!  INTEGER, PARAMETER :: IJKMAX2=169 ! for single particle and 1024 days
!  INTEGER, PARAMETER :: IJKMAX2=3305 ! for single particle and 365 days
!  INTEGER, PARAMETER :: IJKMAX2=34855 ! for single particle and 64 days
!  INTEGER, PARAMETER :: IJKMAX2=73756 ! for single particle and 32 days

!#endif
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
  REAL*4 , SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
#endif
!  REAL*4 :: temp2d_simp(IMT,JMT)
  INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: itemp
!  INTEGER itemp(IMT,JMT)
!  INTEGER :: varid, ncid !, idzrho(IMT,JMT)
  REAL*4 dd,dmult,uint,vint,zint

  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob, depthb,latb
! integer kmm
#endif

LOGICAL around

REAL*8 temp3d_doub(IMT,JMT,KM) ! to be used when botbox is written

  
  alloCondGrid: if ( .not. allocated (botbox) ) then
     allocate (  botbox(IMTG,JMTG,3) ) ! , rhom(IMT,JMT) 
     allocate ( kmu(IMT,JMT)    ,kmv(IMT,JMT) )
#ifdef initxyt
     allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
#endif
  end if alloCondGrid

  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt),temp3d_simp(IMT,JMT,KMM), temp2d_doub(IMT,JMT), temp2d_simp(IMT,JMT), itemp(IMT,JMT)   )
     allocate ( e1t(IMT,JMT) , e2t(IMT,JMT) )
#ifdef tempsalt
     allocate ( tempb(KMM), saltb(KMM), rhob(KMM), depthb(KM), latb(KM))
#endif
  end if alloCondUVW

  start1D  = [ 1]
  count1D  = [KM]
  start2D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count2D  = [         imt,        jmt ,  1 , 1 ]
  map2D    = [          1 ,          2 ,  3 , 4 ]  
  start3D  = [subGridImin ,subGridJmin ,  1 , 1 ]
  count3D  = [         imt,        jmt ,KMM , 1 ]
  map3D    = [          1 ,          2 ,  3 , 4 ]  
    
    
    ! === swap between datasets ===
  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif

  ! === initialise ===
  initFieldcond: if(ints.eq.intstart) then
     ! call coordinat
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif
     ntempus=0 ; ntempusb=0
!     rhom=0.

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     
gridFile = trim(inDataDir)//'topo/mesh_hgr.nc'
!print *,gridFile
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 3751
ierr=NF90_INQ_VARID(ncid,'e1t',varid) ! the main data fields
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e1t(i,:)=temp2d_doub(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e2u',varid) ! dyu in meters
if(ierr.ne.0) stop 3764
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  dyu(i,:)=temp2d_doub(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e2t',varid) ! the main data fields
if(ierr.ne.0) stop 3765
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  e2t(i,:)=temp2d_doub(i,:)
  dxdy(i,:) = e1t(i,:) * e2t(i,:)
enddo

ierr=NF90_INQ_VARID(ncid,'e1v',varid) ! the main data fields
if(ierr.ne.0) stop 3766
ierr=NF90_GET_VAR(ncid,varid,temp2d_doub,start2d,count2d)
if(ierr.ne.0) stop 3799
do i=1,IMT
  dxv(i,:)=temp2d_doub(i,:)
enddo

#ifdef ocra025
ierr=NF90_CLOSE(ncid)
gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
#endif
print *,gridFile


! Depth coordinates
zw(0)    = 0.d0
zw(1:km) = get1DfieldNC(trim(gridFile) ,'e3t_0')
do k=1,km
 kk=km+1-k
 dz(kk)=zw(k)
 zw(k)=zw(k)+zw(k-1)
! print *,k,zw(k),kk,dz(kk)
end do


! The bathymetry
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! kmt field
if(ierr.ne.0) stop 3767
ierr=NF90_GET_VAR(ncid,varid,itemp,start2d,count2d)

kmt=itemp
!do j=JMT,JMT-100,-1
!write(6,'(i4,1x,999i1)') j,(kmt(i,j),i=930,1230)
!enddo
!stop 3956

botbox=0.
! layer thickness at T points
ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
do i=1,IMT
 do j=1,JMT
  if(kmt(i,j).ne.0) then
   botbox(i,j,3)=temp3d_doub(i,j,kmt(i,j))
   if(botbox(i,j,3).eq.0.) then
    print *,i,j,kmt(i,j),botbox(i,j,3),temp3d_doub(i,j,kmt(i,j))
    botbox(i,j,3)=dz(km+1-kmt(i,j))
   endif
  endif
 enddo
enddo



kmu=0 ; kmv=0

do j=1,JMT
 jp=j+1
 if(jp.eq.jmt+1) jp=jmt  ! should be north fold instead
 do i=1,imt
  ip=i+1
  if(ip.eq.IMT+1) ip=1
  kmu(i,j)=min(kmt(i,j),kmt(ip,j),KMM)
  kmv(i,j)=min(kmt(i,j),kmt(i,jp),KMM)
!   if(kmt(i,j).lt.kmt(ip,j)) then ! kmu
!    kmu(i,j)=kmt(i,j)
!    botbox(i,j,1)=botbox(i,j,3)
!   elseif(kmt(i,j).gt.kmt(ip,j)) then
!    kmu(i,j)=kmt(ip,j)
!    botbox(i,j,1)=botbox(ip,j,3)
!   else
!    kmu(i,j)=kmt(i,j)
!    botbox(i,j,1)=min(botbox(i,j,3),botbox(ip,j,3))
!   endif
!   
!   if(kmt(i,j).lt.kmt(i,jp)) then ! kmu
!    kmv(i,j)=kmt(i,j)
!    botbox(i,j,2)=botbox(i,j,3)
!   elseif(kmt(i,j).gt.kmt(i,jp)) then
!    kmv(i,j)=kmt(i,jp)
!    botbox(i,j,2)=botbox(ip,j,3)
!   else
!    kmv(i,j)=kmt(i,j)
!    botbox(i,j,2)=min(botbox(i,j,3),botbox(i,jp,3))
!   endif
   
   
 enddo
enddo

!  north fold 
do i=4,IMT
 ii=IMT+4-i
 kmv(i,JMT)=kmv(ii,JMT-3)
! botbox(i,JMT,2)=botbox(ii,JMT-3,2)
enddo

! layer thickness at u points
ierr=NF90_INQ_VARID(ncid,'e3u',varid) 
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
do i=1,IMT
do j=1,JMT
if(kmu(i,j).ne.0) botbox(i,j,1)=temp3d_doub(i,j,kmu(i,j))
enddo
enddo

! layer thickness at v points
ierr=NF90_INQ_VARID(ncid,'e3v',varid)
if(ierr.ne.0) stop 3763
ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
do i=1,IMT
do j=1,JMT
if(kmv(i,j).ne.0) botbox(i,j,2)=temp3d_doub(i,j,kmv(i,j))
enddo
enddo


!open(21,file=trim(inDataDir)//'topo/botbox',form='unformatted')
!!write(21) botbox
!read(21) botbox
!close(21)
!stop 4966

!ierr=NF90_INQ_VARID(ncid,'nav_lon',varid) ! londitude mesh
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
!do i=100,1000
!ii=IMT+4-i
!print *,i,temp2d_simp(ii,JMT-2)-temp2d_simp(i,JMT)
!enddo
!open(21,file=trim(inDataDir)//'topo/long',form='unformatted')
!write(21) temp2d_simp
!close(21)
!
!ierr=NF90_INQ_VARID(ncid,'nav_lat',varid) ! londitude mesh
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
!open(21,file=trim(inDataDir)//'topo/lat',form='unformatted')
!write(21) temp2d_simp
!close(21)

!do i=100,1000
!ii=IMT+4-i
!print *,i,temp2d_simp(ii,JMT-2)-temp2d_simp(i,JMT)
!enddo
!stop 3957

ierr=NF90_CLOSE(ncid)



!open(21,file=trim(inDataDir)//'topo/kmt',form='unformatted')
!write(21)kmt
!close(21)
!stop 4985



!stop 4956

!i=500 ; j=500
!print *,'subGridImin',subGridImin,subGridJmin,imt,jmt
!print *,'i+subGridImin-1,j+subGridJmin-1',i+subGridImin-1,j+subGridJmin-1
!print *,'kmt',kmt(i,j)
!print *,botbox(i+subGridImin-1,j+subGridJmin-1,1)
!print *,botbox(i+subGridImin-1,j+subGridJmin-1,2)
!print *,botbox(i+subGridImin-1,j+subGridJmin-1,3)
!stop 3957
! Bottom box 
do j=1,JMT
 do i=1,IMT
  if(kmt(i,j).ne.0) then
   dztb(i,j,1)=botbox(i+subGridImin-1,j+subGridJmin-1,3)
   if(botbox(i+subGridImin-1,j+subGridJmin-1,3).eq.0.) then
    print *,i,j,kmt(i,j),botbox(i+subGridImin-1,j+subGridJmin-1,3)
    stop 4957
   endif
!   dztb(i,j,1)=dz(KM+1-kmt(i,j))
  else
   dztb(i,j,1)=0.
  endif
 enddo
enddo

#ifdef initxyt
! Time for individual start positions

!open(84,file=trim(inDataDir)//'topo/masktime_32_005',form='unformatted')
!open(84,file=trim(inDataDir)//'topo/masktime_32_010',form='unformatted')
if(IJKMAX2.eq.7392) open(84,file=trim(inDataDir)//'topo/masktime_32_025',form='unformatted')
!open(84,file=trim(inDataDir)//'topo/masktime_orca1_32_025',form='unformatted')
if(IJKMAX2.eq.169) open(84,file=trim(inDataDir)//'topo/masktime_orca025_1024_single',form='unformatted')
if(IJKMAX2.eq.3305) open(84,file=trim(inDataDir)//'topo/masktime_orca025_365_single',form='unformatted')
if(IJKMAX2.eq.34855) open(84,file=trim(inDataDir)//'topo/masktime_orca025_64_single',form='unformatted')
if(IJKMAX2.eq.73756) open(84,file=trim(inDataDir)//'topo/masktime_orca025_32_single',form='unformatted')
 read(84) trajinit
close(84)
j=0
do k=1,NTID
 do i=1,IJKMAX2
  if(trajinit(k,i,3).ne.0.) then
   j=j+1
#if orca025l75h6
   trajinit(k,i,3)=float(kst2)-0.5
!   print *,j,trajinit(k,i,:)
#endif
  endif
 enddo
! print *,k,j
enddo
ijkst=0
! print *,'ijkmax=',j,IJKMAX2,ijkmax
if(j.ne.IJKMAX2) then
 stop 4396
endif
#endif

ihour=startHour
iday=startDay
imon=startMon-2
if(ngcm.le.24) then
ihour=24-ngcm
iday=startDay-5
imon=startMon
endif
iyear=startYear
     
endif initFieldcond

! date update

#if defined orca1  || orca025
  iday=iday+5
#elif orca025l75h6
  if(ngcm.le.24) then
   ihour=ihour+ngcm
   if(ihour.eq.24) then
    iday=iday+1
    ihour=0
   endif
  else
   imon=imon+1
   if(imon.eq.13) then
    imon=1
    iyear=iyear+1
    if(iyear.eq.2007) iyear=2000
   endif
  endif
#endif

  if(ngcm.le.24) then
  if(iday.gt.idmax(imon,1999)) then
     iday=iday-idmax(imon,1999)
     imon=imon+1
     if(imon.eq.13) then
        imon=1
        iyear=iyear+1
        ! if kan skrivas här om man vill börja om från iyear0
     endif
  endif
  endif
!iyear=2000 ! quick and dirty fix for years
ntime=10000*iyear+100*imon+iday
! file names
#ifdef orca1
 dataprefix='xxxx/ORCA1-N202_xxxxxxxx'
 write(dataprefix(17:24),'(i8)') ntime
 write(dataprefix(1:4),'(i4)') iyear
 fieldFile = trim(inDataDir)//trim(dataprefix)/'d05'
#elif orca025
 dataprefix='xxxx/ORCA025-N112_xxxxxxxx'
 write(dataprefix(19:26),'(i8)') ntime
 write(dataprefix(1:4),'(i4)') iyear
 fieldFile = trim(inDataDir)//trim(dataprefix)//'d05'
#elif orca025l75h6
!print *,'ngcm=',ngcm
!stop 3496
 if(ngcm.eq.  3) dataprefix='ORCA025.L75-SLB2_3h_19900101_19901231_grid_'
! if(ngcm.eq.6) dataprefix='ORCA025.L75-SLB2_6h_19580101_19581231_grid_'
 if(ngcm.eq.  6) dataprefix='ORCA025.L75-SLB2_6h_20050101_20051231_grid_'
 if(ngcm.eq.730) dataprefix='ORCA025.L75-SLB0_730h_19770101_19771231_grid_'
! write(dataprefix(19:26),'(i8)') ntime
 write(dataprefix(23:26),'(i4)') iyear
 write(dataprefix(32:35),'(i4)') iyear
! print *,dataprefix
! stop 4906
! fieldFile = trim(inDataDir)//trim(dataprefix)
 fieldFile = trim(inDataDir)//'fields/'//trim(dataprefix)
! print *,fieldFile
  nread=mod(ints-1,intmax)+1
  start2D  = [subGridImin ,subGridJmin ,  1 , nread ]
  start3D  = [subGridImin ,subGridJmin ,  1 , nread ]
#else
 stop 39573
#endif

    
! Sea surface height
ierr=NF90_OPEN(trim(fieldFile)//'T.nc',NF90_NOWRITE,ncid)
ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
if(ierr.ne.0) then
print *,ints,trim(fieldFile)//'T.nc'
stop 3768
endif
ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

do j=1,JMT
 do i=1,IMT+1
  ii=i
  if(ii.eq.IMT+1) ii=1
  hs(i,j,2)=temp2d_simp(ii,j)
 enddo
enddo

do i=4,IMT
 ii=IMT+4-i
 hs(i,JMT+1,2)=hs(ii,JMT-3,2)  !  north fold 
enddo

#ifdef tempsalt 

! Temperature
gridFile = trim(fieldFile)//'T.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'votemper',varid) 
if(ierr.ne.0) stop 3769
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)
do i=1,IMT
 do j=1,JMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   tem(i,j,kk,2)=temp3d_simp(i,j,k)
  enddo
 enddo
enddo

   
! Salinity
!gridFile = trim(fieldFile)//'d05T.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'vosaline',varid) 
if(ierr.ne.0) stop 3769
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)
do i=1,IMT
 do j=1,JMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   sal(i,j,kk,2)=temp3d_simp(i,j,k)
!   print *,i,j,k,sal(i,j,kk,2)
  enddo
 enddo
enddo

!kmm=KM
depthb=0.
do j=1,JMT
 latb=-80+0.25*float(j+subGridJmin-1)
 do i=1,IMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   tempb(k)=tem(i,j,kk,2)
   saltb(k)=sal(i,j,kk,2)
  enddo
  call statvd(tempb, saltb, rhob ,KM ,depthb ,latb)
!  call statvd(tempb, saltb, rhob ,kmm ,depthb ,latb)
  do k=1,kmt(i,j)
   kk=KM+1-k
   rho(i,j,kk,2)=rhob(k)-1000.
!   rhom(i,j,k)=rhom(i,j,k)+rho(i,j,kk,2)
  enddo
 enddo
enddo

!if(ints.eq.intstart+intend-1) then
!print *,'write rho',ints,intstart,intend,intend+intstart-1,trim(inDataDir)//'rho0_1year'
!rhom=rhom/float(intend)
!if(ints.eq.intstart) open(38,file=trim(inDataDir)//'dzrho',form='unformatted')
!write(38) rhom
!close(38)
!print *,rhom
!endif

! Norske havet
!if(ints.eq.intstart) open(38,file=trim(inDataDir)//'dzrho_1996',form='unformatted')
!idzrho=0
!do i=1,IMT
! do j=1,JMT
!!  do k=KM,1,-1
!!   if(rho(i,j,k,2).ge.27.9) idzrho(i,j)=k
!!  enddo
!  do k=2,KM
!   if(rho(i,j,k-1,2).ge.27.9 .and. rho(i,j,k,2).le.27.9) idzrho(i,j)=KM+1-k
!  enddo
! enddo
!enddo
!write(38) idzrho


#endif     

dmult=1.  ! amplification of the velocity amplitude by simple multiplication

! u velocity
gridFile = trim(fieldFile)//'U.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'vozocrtx',varid) 
if(ierr.ne.0) stop 3769
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)
do i=1,IMT
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i+1,j,2))
   if(k.eq.kmu(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,1)
   uflux(i,j,kk,2)=temp3d_simp(i,j,k) * dyu(i,j) * dd * dmult
!   if(botbox(i+subGridImin-1,j+subGridJmin-1,1).eq.0.) then
!    print *,i,j,kk,k,uflux(i,j,kk,2)
!    stop 4967
!   endif
  enddo
 enddo
enddo

! v velocity
gridFile = trim(fieldFile)//'V.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! kmt field
if(ierr.ne.0) stop 3770
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)


!  north fold 
do i=4,IMT
 ii=IMT+4-i
! vflux(i,JMT,:,2)=-vflux(ii,JMT-3,:,2)
enddo

do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
!  do k=1,KM
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i,j+1,2))
   if(k.eq.kmv(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,2)
   vflux(i,j,kk,2)=temp3d_simp(i,j,k) * dxv(i,j) * dd * dmult
!   if(i.eq.IMT) print *,i,j,kk,k,vflux(i,j,kk,2)
!   if(botbox(i+subGridImin-1,j+subGridJmin-1,2).eq.0.) then
!    print *,i,j,kk,k,vflux(i,j,kk,2)
!    stop 4968
!   endif

  enddo
 enddo
enddo



#ifdef drifter
! average velocity/transport over surface drifter drogue depth to simulate drifter trajectories
kbot=65 ; ktop=66 ! number of surface layers to integrat over
do i=1,imt
 do j=1,jmt
 
  uint=0. ; vint=0. ; zint=0.
  if(ktop.eq.KM) zint=hs(i,j,2)
  do k=kbot,ktop
   uint=uint+uflux(i,j,k,2) ! integrated transport
   vint=vint+vflux(i,j,k,2)
   zint=zint+dz(k)          ! total depth of drougued drifter
  enddo
  ! weighted transport for each layer
  do k=kbot,KM
   if(k.ne.KM) then
    uflux(i,j,k,2)=uint*dz(k)/zint 
    vflux(i,j,k,2)=vint*dz(k)/zint
   else
    uflux(i,j,k,2)=uint*(hs(i,j,2)+dz(k))/zint
    vflux(i,j,k,2)=vint*(hs(i,j,2)+dz(k))/zint
   endif
  enddo

 enddo
enddo

#endif


#ifdef initxyt
! Set the initial trajectory positions
!ijkst(:,5)=ntimask(ntempus,:)
#ifdef orca025l75h6
if( mod(ints,24/ngcm*5).eq.1 .or. ints.le.2) ntempus=ntempus+1
if(ntempus.ne.ntempusb .and. ntempus.le.NTID) then
ntempusb=ntempus
!print *,'ints=',ints,' ntempus=',ntempus,' ntempusb=',ntempusb
#else
if(ints.le.NTID) then
#endif
do ntrac=1,ijkmax
 if(trajinit(ntempus,ntrac,3).ne.0.) then
  ijkst(ntrac,4)=0
  ijkst(ntrac,5)=5
  ijkst(ntrac,6)=ijkst(ntrac,6)+1
  do l=1,3
   ijkst(ntrac,l)=trajinit(ntempus,ntrac,l)+1
   trj(ntrac,l)=trajinit(ntempus,ntrac,l)
!   if(l.eq.1) print *,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
   if(trj(ntrac,l).gt.float(ijkst(ntrac,l)) .or. trj(ntrac,l).lt.float(ijkst(ntrac,l)-1)) then
    print *,l,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
    stop 3946
   endif
  enddo
 else
  ijkst(ntrac,5)=0
  ijkst(ntrac,6)=0
 endif
enddo
endif
#ifdef orca025l75h6

#endif
if( mod(ints,24/ngcm*5).ne.1 .and. ints.gt.2) then
 ijkst=0 
endif
#endif


!uflux=0.0001 ; vflux=0.0001  ! special case with no velocities

! multiply the velocity field by a constant
!uflux=1.3*uflux
!vflux=1.3*vflux




     deallocate ( ssh , temp3d_simp, temp2d_doub, temp2d_simp, itemp )
     deallocate ( e1t , e2t )
#ifdef tempsalt
     deallocate ( tempb, saltb, rhob, depthb, latb )
#endif

!print *,'uv',uflux(1157,257,61,2),uflux(1158,257,61,2),vflux(1158,257,61,2),vflux(1158,256,61,2)

  return
end subroutine readfields



