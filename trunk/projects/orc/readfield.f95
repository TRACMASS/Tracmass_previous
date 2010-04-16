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
  INTEGER, SAVE                                :: ntempus

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  
    ! = Variables for filename generation
  CHARACTER (len=200)                        :: dataprefix !, dstamp

#ifdef orca1
  INTEGER, PARAMETER :: IMTG=???,JMTG=???
#else
  INTEGER, PARAMETER :: IMTG=1440,JMTG=1021
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
!#endif
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
  REAL*4 , SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
#endif
!  REAL*4 :: temp2d_simp(IMT,JMT)
  INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: itemp
!  INTEGER itemp(IMT,JMT)
  INTEGER :: varid, ncid !, idzrho(IMT,JMT)
  REAL*4 dd,dmult
  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob, depthb,latb
 integer kmm
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
     allocate ( ssh(imt,jmt),temp3d_simp(IMT,JMT,KM), temp2d_doub(IMT,JMT), temp2d_simp(IMT,JMT), itemp(IMT,JMT)   )
     allocate ( e1t(IMT,JMT) , e2t(IMT,JMT) )
#ifdef tempsalt
     allocate ( tempb(KM), saltb(KM), rhob(KM), depthb(KM), latb(KM))
#endif
  end if alloCondUVW

  start1D  = [ 1]
  count1D  = [km]
  start2D  = [subGridImin ,subGridJmin , 1 ,1]
  count2D  = [      imt,        jmt , 1 ,1]
  map2D    = [          1 ,          2 , 3 ,4]  
  start3D  = [subGridImin ,subGridJmin , 1 ,1]
  count3D  = [      imt,        jmt ,km ,1]
  map3D    = [          1 ,          2 , 3 ,4]  
  
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
     ntempus=0
!     rhom=0.

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     
gridFile = trim(inDataDir)//'topo/mesh_hgr.nc'

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
ierr=NF90_CLOSE(ncid)

! The bathymetry
gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
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
   kmu(i,j)=min(itemp(i,j),itemp(ip,j))
   kmv(i,j)=min(itemp(i,j),itemp(i,jp))
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
  else
   dztb(i,j,1)=0.
  endif
 enddo
enddo

#ifdef initxyt
! Time for individual start positions

!open(84,file=trim(inDataDir)//'topo/masktime_32_005',form='unformatted')
!open(84,file=trim(inDataDir)//'topo/masktime_32_010',form='unformatted')
open(84,file=trim(inDataDir)//'topo/masktime_32_025',form='unformatted')
!open(84,file=trim(inDataDir)//'topo/masktime_orca1_32_025',form='unformatted')
! read(84) ntimask
 read(84) trajinit
close(84)
j=0
do k=1,NTID
 do i=1,IJKMAX2
  if(trajinit(k,i,3).ne.0.) j=j+1
 enddo
enddo
ijkst=0
! print *,'ijkmax=',j,IJKMAX2,ijkmax
if(j.ne.IJKMAX2) then
 stop 4396
endif
#endif

iday=startDay-5
imon=startMon
iyear=startYear
     
endif initFieldcond

! date update
  iday=iday+5
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
!ntempus=ntempus+1 ! quick and dirty fix to be generalised
ntempus=ints
!print *,'ntempus=',ints,ntempus,ntime
! file names
#ifdef orca1
 dataprefix='xxxx/ORCA1-N202_xxxxxxxx'
write(dataprefix(17:24),'(i8)') ntime
write(dataprefix(1:4),'(i4)') iyear
#else
 dataprefix='xxxx/ORCA025-N112_xxxxxxxx'
write(dataprefix(19:26),'(i8)') ntime
write(dataprefix(1:4),'(i4)') iyear
#endif

  fieldFile = trim(inDataDir)//trim(dataprefix)
    
! temp, salt and ssh
!print *,trim(fieldFile)//'d05T.nc'

ierr=NF90_OPEN(trim(fieldFile)//'d05T.nc',NF90_NOWRITE,ncid)
ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
if(ierr.ne.0) stop 3768
ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)

do j=1,JMT
 do i=1,IMT+1
  ip=i
  if(ip.eq.IMT+1) ip=1
  hs(i,j,2)=temp2d_simp(ip,j)
 enddo
enddo

do i=4,IMT+1
 hs(i,jmt+1,2) =hs(IMT+4-i,jmt-3,2)  !  north fold 
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
gridFile = trim(fieldFile)//'d05T.nc'
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

kmm=KM
depthb=0.
do j=1,JMT
 latb=-80+0.25*float(j+subGridJmin-1)
 do i=1,IMT
  do k=1,kmt(i,j)
   kk=KM+1-k
   tempb(k)=tem(i,j,kk,2)
   saltb(k)=sal(i,j,kk,2)
  enddo
  call statvd(tempb, saltb, rhob ,kmm ,depthb ,latb)
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
gridFile = trim(fieldFile)//'d05U.nc'
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
  enddo
 enddo
enddo

! v velocity
!temp3d_simp = get3DfieldNC(trim(fieldFile)//'d05V.nc' ,'vomecrty')
gridFile = trim(fieldFile)//'d05V.nc'
ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 5751
ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! kmt field
if(ierr.ne.0) stop 3770
ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
if(ierr.ne.0) stop 3799
ierr=NF90_CLOSE(ncid)
do i=1,IMT
 do j=1,JMT
  do k=1,kmv(i,j)
   kk=KM+1-k
   dd = dz(kk) 
   if(k.eq.1) dd = dd + 0.5*(hs(i,j,2) + hs(i,j+1,2))
   if(k.eq.kmv(i,j)) dd = botbox(i+subGridImin-1,j+subGridJmin-1,2)
   vflux(i,j,kk,2)=temp3d_simp(i,j,k) * dxv(i,j) * dd * dmult
  enddo
 enddo
enddo


#ifdef initxyt
! Set the initial trajectory positions
!ijkst(:,5)=ntimask(ntempus,:)
j=0
do ntrac=1,ijkmax
 if(trajinit(ntempus,ntrac,3).ne.0.) then
!  j=j+ntimask(ntempus,ntrac,3)
  ijkst(ntrac,4)=0
  ijkst(ntrac,5)=5
  ijkst(ntrac,6)=1
  do l=1,3
   ijkst(ntrac,l)=trajinit(ntempus,ntrac,l)+1
   trj(ntrac,l)=trajinit(ntempus,ntrac,l)
!   if(ntrac.eq.12) print *,l,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
!   if(ntrac.eq.12) print *,ntempus,ntimask(ntempus,ntrac,l)
   if(trj(ntrac,l).gt.float(ijkst(ntrac,l)) .or. trj(ntrac,l).lt.float(ijkst(ntrac,l)-1)) then
    print *,l,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
    !trj(ntrac,l)=float(ijkst(ntrac,l))-0.5
    stop 3946
   endif
  enddo
 else
  ijkst(ntrac,5)=0
  ijkst(ntrac,6)=0
 endif
enddo
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

  return
end subroutine readfields



