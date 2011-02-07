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
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm,ipp,jmm,jpp,ntrac, l
  INTEGER                                      :: kbot,ktop
  INTEGER, SAVE                                :: ntempus,ntempusb,nread

  CHARACTER (len=200)                          :: gridFile ,fieldFile,dataprefix,zfile,rfile

  INTEGER, PARAMETER :: IMTG=619,JMTG=523,KMM=84

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

!REAL*8 temp3d_doub(IMT+2,JMT,KM) ! to be used when botbox is written

  
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
  dzt(:,:,:,1)=dzt(:,:,:,2)
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
     ntempus=89 ; ntempusb=0
!     rhom=0.

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     
gridFile = trim(inDataDir)//'topo/mesh_mask_baltix.nc'
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
!print *,dyu
!stop 4965

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

!#ifdef ocra025
!ierr=NF90_CLOSE(ncid)
!gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
!#endif

! The bathymetry
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'mbathy',varid) ! kmt field
if(ierr.ne.0) stop 3767
ierr=NF90_GET_VAR(ncid,varid,itemp,start2d,count2d)

do i=1,IMT
  kmt(i,:)=itemp(i,:)
enddo

kmu=0 ; kmv=0

do j=1,jmt
jp=j+1
if(jp.eq.jmt+1) jp=jmt  ! should be north fold instead
 do i=1,imt
 ip=i+1
 if(ip.eq.IMT+1) ip=1
   kmu(i,j)=min(itemp(i,j),itemp(ip,j),KMM)
   kmv(i,j)=min(itemp(i,j),itemp(i,jp),KMM)
 enddo
enddo

!botbox=0.
!! layer thickness at u points
!ierr=NF90_INQ_VARID(ncid,'e3u',varid) 
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!do i=1,IMT
!do j=1,JMT
!if(kmu(i,j).ne.0) botbox(i,j,1)=temp3d_doub(i,j,kmu(i,j))
!enddo
!enddo
!
!! layer thickness at v points
!ierr=NF90_INQ_VARID(ncid,'e3v',varid)
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!do i=1,IMT
!do j=1,JMT
!if(kmv(i,j).ne.0) botbox(i,j,2)=temp3d_doub(i,j,kmv(i,j))
!enddo
!enddo
!
!! layer thickness at T points
!ierr=NF90_INQ_VARID(ncid,'e3t',varid) 
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp3d_doub,start3d,count3d)
!do i=1,IMT
!do j=1,JMT
!if(kmt(i,j).ne.0) botbox(i,j,3)=temp3d_doub(i,j,kmt(i,j))
!enddo
!enddo


open(21,file=trim(inDataDir)//'topo/botbox',form='unformatted')
!write(21) botbox
read(21) botbox
close(21)
!stop 4966

!ierr=NF90_INQ_VARID(ncid,'nav_lon',varid) ! londitude mesh
!if(ierr.ne.0) stop 3763
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
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


ierr=NF90_CLOSE(ncid)



!open(21,file=trim(inDataDir)//'topo/kmt',form='unformatted')
!write(21)kmt
!close(21)
!stop 4985


! Depth coordinates
zw(0)    = 0.d0
zw(1:km) = get1DfieldNC(trim(gridFile) ,'e3t_0')
do k=1,km
 kk=km+1-k
 dz(kk)=zw(k)
 zw(k)=zw(k)+zw(k-1)
! print *,k,zw(k),kk,dz(kk)
end do
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

ihour=24-ngcm
iday=startDay-5
imon=startMon
iyear=startYear
     
endif initFieldcond

! date update

  ntempus=ntempus+1
  ihour=ihour+ngcm
  if(ihour.eq.24) then
   iday=iday+1
   ihour=0
  endif

  if(iday.gt.idmax(imon,1999)) then
     iday=iday-idmax(imon,1999)
     imon=imon+1
     if(imon.eq.13) then
        imon=1
        iyear=iyear+1
        ! if kan skrivas här om man vill börja om från iyear0
     endif
  endif
!iyear=2000 ! quick and dirty fix for years
ntime=10000*iyear+100*imon+iday
! file names
if(ntempus.le.9) then
 dataprefix='xxxx/x_BALTIX4_3h_19840101_19841231_grid_'
 write(dataprefix(1:4),'(i4)') iyear
 write(dataprefix(6:6),'(i1)') ntempus
elseif(ntempus.le.99) then
 dataprefix='xxxx/xx_BALTIX4_3h_19840101_19841231_grid_'
 write(dataprefix(1:4),'(i4)') iyear
 write(dataprefix(6:7),'(i2)') ntempus
elseif(ntempus.le.999) then
 dataprefix='xxxx/xxx_BALTIX4_3h_19840101_19841231_grid_'
 write(dataprefix(1:4),'(i4)') iyear
 write(dataprefix(6:8),'(i3)') ntempus
else
 stop 4956
endif
! write(dataprefix(19:26),'(i8)') ntime

 fieldFile = trim(inDataDir)//trim(dataprefix)
! print *,fieldFile
!  nread=mod(ints-1,intmax)+1
  nread=1
  start2D  = [subGridImin ,subGridJmin ,  1 , nread ]
  start3D  = [subGridImin ,subGridJmin ,  1 , nread ]
!#else
! stop 39573
!#endif


fieldFile = trim(inDataDir)//trim(dataprefix)//'T.nc'
!print *,fieldFile
inquire(file=trim(fieldFile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
stop 4555
endif
zfile='gzip -c -d '//trim(fieldFile)//'.gz > '//trim(inDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(inDataDir)//'tmp/'//trim(outDataFile)
inquire(file=trim(rfile),exist=around)
if(.not.around) stop 4556



ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
if(ierr.ne.0) then
print *,ints,trim(fieldFile)//'T.nc'
stop 3768
endif
ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

do j=1,JMT
 do i=1,IMT
  hs(i,j,2)=temp2d_simp(i,j)
 enddo
enddo

do i=4,IMT+1
 hs(i,jmt+1,2) =hs(IMT+4-i,jmt-3,2)  !  north fold 
enddo

! z-star calculations of the layer thicknesses   e3w(k)* = (H - eta )/(H + eta)  *e3w(k)
! eller DZ=DZ*×(H+eta)/H 
do i=1,IMT
 do j=1,JMT
  do k=1,KM
   kk=KM+1-k
   if(kmt(i,j).eq.kk) then
    dzt(i,j,k,2)=botbox(i+subGridImin-1,j+subGridJmin-1,3)
    dzt(i,j,k,2)=dzt(i,j,k,2)*( zw(kmt(i,j))+hs(i,j,2) )/zw(kmt(i,j))
   elseif(kmt(i,j).ne.0) then
    dzt(i,j,k,2)=dz(k)
    dzt(i,j,k,2)=dzt(i,j,k,2)*( zw(kmt(i,j))+hs(i,j,2) )/zw(kmt(i,j))
   else
    dzt(i,j,k,2)=0.
   endif
  enddo
 enddo
enddo


#ifdef tempsalt 

! Read the net downward heat flux
!ierr=NF90_OPEN(trim(fieldFile)//'d05T.nc',NF90_NOWRITE,ncid)
!ierr=NF90_INQ_VARID(ncid,'sohefldo',varid) ! the main data fields
!if(ierr.ne.0) stop 3768
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
!if(ierr.ne.0) stop 3799
!ierr=NF90_CLOSE(ncid)
!!rhom=rhom+temp2d_simp
!rhom=temp2d_simp


! Temperature
!gridFile = trim(fieldFile)//'T.nc'
ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5752
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
ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5753
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


fieldFile = trim(inDataDir)//trim(dataprefix)//'U.nc'
!print *,fieldFile
inquire(file=trim(fieldFile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
stop 4555
endif
zfile='gzip -c -d '//trim(fieldFile)//'.gz > '//trim(inDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(inDataDir)//'tmp/'//trim(outDataFile)
inquire(file=trim(rfile),exist=around)
if(.not.around) stop 4556

! u velocity
ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5754
ierr=NF90_INQ_VARID(ncid,'vozocrtx',varid) 
if(ierr.ne.0) stop 3769
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)
do i=1,IMT-1
 do j=1,JMT
  do k=1,kmu(i,j)
   kk=KM+1-k
   dd = 0.5*(dzt(i,j,kk,2)+dzt(i+1,j,kk,2))
!   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i+1,j,2))
!   if(k.eq.kmu(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,1)
   uflux(i,j,kk,2)=temp3d_simp(i,j,k) * dyu(i,j) * dd * dmult
!   print *,i,j,kk,k,uflux(i,j,kk,2)
  enddo
 enddo
enddo

! v velocity
!temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05V.nc' ,'vomecrty')

fieldFile = trim(inDataDir)//trim(dataprefix)//'V.nc'
!print *,fieldFile
inquire(file=trim(fieldFile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
stop 4555
endif
zfile='gzip -c -d '//trim(fieldFile)//'.gz > '//trim(inDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(inDataDir)//'tmp/'//trim(outDataFile)
inquire(file=trim(rfile),exist=around)
if(.not.around) stop 4556


ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5755
ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! kmt field
if(ierr.ne.0) stop 3770
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)
do i=1,IMT
 do j=1,JMT-1
  do k=1,kmv(i,j)
   kk=KM+1-k
   dd = 0.5*(dzt(i,j,kk,2)+dzt(i,j+1,kk,2))
!   dd = dz(kk) 
!   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i,j+1,2))
!   if(k.eq.kmv(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,2)
   vflux(i,j,kk,2)=temp3d_simp(i,j,k) * dxv(i,j) * dd * dmult
!   if(i.eq.100) print *,i,j,kk,k,vflux(i,j,kk,2)
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



