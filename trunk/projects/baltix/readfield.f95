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
  
  
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm,ipp,jmm,jpp,ntrac, l
  INTEGER                                      :: kbot,ktop
  INTEGER, SAVE                                :: ntempus(1:12,1:31,0:21),ntempusb,nread

  CHARACTER (len=200)                          :: fieldFile,dataprefix,zfile,rfile

  INTEGER, PARAMETER :: IMTG=619,JMTG=523,KMM=84
  
  REAL*4,  ALLOCATABLE, DIMENSION(:,:)         :: ssh, temp2d_simp
  REAL*4,  ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_simp
#ifdef tempsalt
  REAL*4,  ALLOCATABLE, DIMENSION(:)           :: tempb, saltb, rhob, depthb,latb
#endif

  REAL*4 dd,dmult,uint,vint,zint
  LOGICAL around
#ifdef initxyt
  INTEGER, PARAMETER :: NTID=73
!  INTEGER, PARAMETER :: IJKMAX2=254
!  INTEGER, PARAMETER :: IJKMAX2=? ! for distmax=0.05 and 32 days
!  INTEGER, PARAMETER :: IJKMAX2=? ! for distmax=0.10 and 32 days
  INTEGER, PARAMETER :: IJKMAX2=7392 ! for distmax=0.25 and 32 days
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
  REAL*4 , SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
#endif  


! === ALLOCATE READ MATRICES ===

#ifdef initxyt
    alloCondGrid: if ( .not. allocated (ntimask) ) then
        allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
    end if alloCondGrid
#endif


alloCondUVW: if(.not. allocated (ssh)) then
    allocate ( ssh(imt,jmt),temp3d_simp(IMT,JMT,KMM), temp2d_simp(IMT,JMT))
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
    
    ! Fill ntempus which is a matrix containing the prefix number for the files
    ! ntempus is 0 at 1 Jan 03:00 and 2918 at 31 Dec 21:00
    k=-1
    do imon=1,12
        do iday=1,idmax(imon,1999)
            do ihour=0,21,3
                ntempus(imon,iday,ihour)=k
                k=k+1
            enddo
        enddo
    enddo

    ihour=startHour-ngcm
    iday=startDay
    imon=startMon
    iyear=startYear

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
    
endif initFieldcond

  ! === date update ===
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
  
ntime=1000000*iyear+10000*imon+100*iday+ihour
!print*,iyear,imon,iday,ihour,ntempus(imon,iday,ihour)

! === 1 Jan at 00:00 is named as 31 Dec 24:00 the previous year ===
if(imon == 1 .and. iday == 1 .and. ihour == 0) then
 dataprefix='2919_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
 write(dataprefix(17:20),'(i4)') iyear-1
 write(dataprefix(26:29),'(i4)') iyear-1
else 
 if(ntempus(imon,iday,ihour).le.9) then
  dataprefix=   'x_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
  write(dataprefix(1:1)  ,'(i1)') ntempus(imon,iday,ihour)
  write(dataprefix(14:17),'(i4)') iyear
  write(dataprefix(23:26),'(i4)') iyear
 elseif(ntempus(imon,iday,ihour).le.99) then
  dataprefix=  'xx_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
  write(dataprefix(1:2),'(i2)') ntempus(imon,iday,ihour)
  write(dataprefix(15:18),'(i4)') iyear
  write(dataprefix(24:27),'(i4)') iyear
 elseif(ntempus(imon,iday,ihour).le.999) then
  dataprefix= 'xxx_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
  write(dataprefix(1:3),'(i3)') ntempus(imon,iday,ihour)
  write(dataprefix(16:19),'(i4)') iyear
  write(dataprefix(25:28),'(i4)') iyear
 elseif(ntempus(imon,iday,ihour).le.9999) then
  dataprefix='xxxx_BALTIX4_3h_xxxx0101_xxxx1231_grid_'
  write(dataprefix(1:4),'(i4)') ntempus(imon,iday,ihour)
  write(dataprefix(17:20),'(i4)') iyear
  write(dataprefix(26:29),'(i4)') iyear
 else
  stop 4956
 endif
endif

 fieldFile = trim(inDataDir)//trim(dataprefix)
!  nread=mod(ints-1,intmax)+1
  nread=1
  start2D  = [subGridImin ,subGridJmin ,  1 , nread ]
  start3D  = [subGridImin ,subGridJmin ,  1 , nread ]
!#else
! stop 39573
!#endif

! === Unzip the T file ===
fieldFile = trim(inDataDir)//trim(dataprefix)//'T.nc'
inquire(file=trim(fieldFile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',fieldFile,ntempus(imon,iday,ihour),iyear,imon,iday,ihour
stop 4555
endif
zfile='gzip -c -d '//trim(fieldFile)//'.gz > '//trim(outDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(outDataDir)//'tmp/'//trim(outDataFile)
inquire(file=trim(rfile),exist=around)
if(.not.around) stop 4556

! === Open T file ===
ierr=NF90_OPEN(trim(rfile),NF90_NOWRITE,ncid)
ierr=NF90_INQ_VARID(ncid,'sossheig',varid) ! the main data fields
if(ierr.ne.0) then
!print *,ints,trim(fieldFile)//'T.nc'
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

#ifdef tempsalt 

! Temperature
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
  enddo
 enddo
enddo

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
  do k=1,kmt(i,j)
   kk=KM+1-k
   rho(i,j,kk,2)=rhob(k)-1000.
  enddo
 enddo
enddo


#endif     

dmult=1.  ! amplification of the velocity amplitude by simple multiplication

fieldFile = trim(inDataDir)//trim(dataprefix)//'U.nc'
inquire(file=trim(fieldFile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
stop 4555
endif
zfile='gzip -c -d '//trim(fieldFile)//'.gz > '//trim(outDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(outDataDir)//'tmp/'//trim(outDataFile)
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
   !dd = 0.5*(dzt(i,j,kk,2)+dzt(i+1,j,kk,2))
   dd = dz(kk)
   if(k == kmu(i,j)) dd = min(dztb(i,j,2),dztb(i+1,j,2))
   if(kmu(i,j) <= 0) dd = 0.
   dd = dd*( zw(kmt(i,j))+(hs(i,j,2)+hs(i+1,j,2))/2. )/zw(kmt(i,j))
   !print*,k,kk,kmu(i,j),dz(kk),dztb(i,j,2),dztb(i+1,j,2),dd
   uflux(i,j,kk,2)=temp3d_simp(i,j,k) * dyu(i,j) * dd * dmult
  enddo
 enddo
enddo

! v velocity
fieldFile = trim(inDataDir)//trim(dataprefix)//'V.nc'
inquire(file=trim(fieldFile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',fieldFile,ntempus,iyear,imon,iday,ihour
stop 4555
endif
zfile='gzip -c -d '//trim(fieldFile)//'.gz > '//trim(outDataDir)//'tmp/'//trim(outDataFile)
CALL system(zfile)
rfile=trim(outDataDir)//'tmp/'//trim(outDataFile)
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
   !dd = 0.5*(dzt(i,j,kk,2)+dzt(i,j+1,kk,2))
   dd = dz(kk)
   if(k == kmv(i,j)) dd = min(dztb(i,j,2),dztb(i,j+1,2))
   if(kmv(i,j) <= 0) dd = 0.
   dd = dd*( zw(kmt(i,j))+(hs(i,j,2)+hs(i,j+1,2))/2. )/zw(kmt(i,j))
   vflux(i,j,kk,2)=temp3d_simp(i,j,k) * dxv(i,j) * dd * dmult
  enddo
 enddo
enddo

! z-star calculations of the layer thicknesses dz* = dz (H+eta)/H 
do i=1,IMT
 do j=1,JMT
  do k=1,KM
   kk=KM+1-k
   if(kmt(i,j).eq.kk) then
    dztb(i,j,1)=dztb(i,j,1)*( zw(kmt(i,j))+hs(i,j,2) )/zw(kmt(i,j))
    dzt(i,j,k,2)=dztb(i,j,1)
   elseif(kmt(i,j).ne.0) then
    dzt(i,j,k,2)=dz(k)
    dzt(i,j,k,2)=dzt(i,j,k,2)*( zw(kmt(i,j))+hs(i,j,2) )/zw(kmt(i,j))
   else
    dzt(i,j,k,2)=0.
   endif
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


     deallocate ( ssh , temp3d_simp, temp2d_simp )
#ifdef tempsalt
     deallocate ( tempb, saltb, rhob, depthb, latb )
#endif


  return
end subroutine readfields



