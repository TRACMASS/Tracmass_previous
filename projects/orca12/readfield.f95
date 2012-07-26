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
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
  INTEGER                                      :: kbot,ktop
  INTEGER, SAVE                                :: ntempus,ntempusb,nread

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  ! = Variables for filename generation
  CHARACTER (len=200)                        :: dataprefix


  REAL*4, ALLOCATABLE, DIMENSION(:,:)         :: temp2d_simp
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_simp

#ifdef initxyt
  INTEGER, PARAMETER :: NTID=73
  INTEGER, PARAMETER :: IJKMAX2=7392 ! for distmax=0.25 and 32 days

  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
  REAL*4 , SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
#endif
  REAL*4 dd,hu,hv,uint,vint,zint
  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob, depthb,latb
#endif

 LOGICAL around

!---------------------------------------------------------------

#ifdef initxyt
  alloCondGrid: if ( .not. allocated (ntimask) ) then
     allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
  endif alloCondGrid
#endif

 alloCondUVW: if(.not. allocated (temp2d_simp)) then
   allocate ( temp3d_simp(IMT,JMT,KM), temp2d_simp(IMT,JMT)  )
#ifdef tempsalt
   allocate ( tempb(KM), saltb(KM), rhob(KM), depthb(KM), latb(KM))
#endif
 endif alloCondUVW

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

! -------------------------------------------------------------

! === Initialising fields ===
  initFieldcond: if(ints.eq.intstart) then
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
     dzt=0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif
     ntempus=0 ; ntempusb=0

#ifdef initxyt
! Time for individual start positions
if(IJKMAX2.eq.7392) open(84,file=trim(inDataDir)//'topo/masktime_32_025',form='unformatted')
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
imon=startMon
iyear=startYear
if(ngcm.le.24) then
ihour=24-ngcm
iday=startDay-5
imon=startMon
endif
     
else

! ----------------------------------------------------------------

! === Update clockworks ===
  iday=iday+5
  if(iday > idmax(imon,1999)) then
    iday=iday-idmax(imon,1999)
    imon=imon+1
    if(imon == 13) then
       imon=1
       iyear=iyear+1
    endif
  endif
#if orca025l75h6
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

endif initFieldcond

! === Time number ===
ntime=10000*iyear+100*imon+iday

! ------------------------------------------------------------

! === Find the file for this timestep ===

 dataprefix='xxxx/ORCA0083-N01_xxxxxxxx'
 write(dataprefix(19:26),'(i8)') ntime
 write(dataprefix(1:4),'(i4)') iyear
 fieldFile = trim(inDataDir)//'fields/'//trim(dataprefix)//'d05'

!print *,'fieldFile=',trim(fieldFile)
    
! Sea surface height
ierr=NF90_OPEN(trim(fieldFile)//'T.nc',NF90_NOWRITE,ncid)
if(ierr.ne.0) stop 3766
ierr=NF90_INQ_VARID(ncid,'sossheig',varid)
if(ierr.ne.0) stop 3768
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

!!------------------------------------------------------------------------------
   ! Compute the level thickness of all boxes tking account of the z-star coordinates
   ! withayer thicknesses dz* = dz (H+ssh)/H and the variable bottom box
   do i=1,IMT
      do j=1,JMT
         do k=1,KM
            kk = KM+1-k
            if (kmt(i,j) == kk) then ! for the bottom box
               dzt(i,j,k,2) = dztb(i,j,1) * (zw(kmt(i,j)) + hs(i,j,2)) / zw(kmt(i,j))
               if(dzt(i,j,k,2).eq.0.) then
                print *,i,j,kmt(i,j),dztb(i,j,1),hs(i,j,2),zw(kmt(i,j))
                stop 4967
               endif
            elseif (kmt(i,j) /= 0) then ! for the levels above the bottom box
               dzt(i,j,k,2) = dz(k)       * (zw(kmt(i,j)) + hs(i,j,2)) / zw(kmt(i,j))
               if(dzt(i,j,k,2).eq.0.) then
                print *,i,j,kmt(i,j),dz(k),hs(i,j,2),zw(kmt(i,j))
                stop 4968
               endif
            else
               dzt(i,j,k,2) = 0.
            endif
         enddo
      enddo
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
  enddo
 enddo
enddo

depthb=0.
do j=1,JMT
 latb=-80+1./12.*float(j+subGridJmin-1)
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
   	  ip=i+1
   	  if(i.eq.IMT) ip=1
      do j=1,JMT
       hu=min(zw(kmt(i,j)),zw(kmt(ip,j))) ! total depth at u point
         do kk=1,kmu(i,j)
            k=KM+1-kk
            dd = dz(k)
            if (kk == kmu(i,j)) THEN
               dd = MIN (dztb(i,j,2),dztb(ip,j,2))
            endif
            if (kmu(i,j) <= 0) THEN
               dd = 0.
            endif
            ! thickness of the wall at the u-point
            if (kmu(i,j) /= 0) then
             dd = dd * ( hu + 0.5*(hs(i,j,2) + hs(ip,j,2)) ) / hu 
             uflux(i,j,k,2) = temp3d_simp(i,j,kk) * dyu(i,j) * dd 
            endif
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


   do i=1,IMT
      do j=1,JMT-1
       hv=min(zw(kmt(i,j)),zw(kmt(i,j+1))) ! total depth at v point
         do kk=1,kmv(i,j)
            k = KM+1-kk
            dd = dz(k)
            if (kk == kmv(i,j)) THEN
               dd = MIN (dztb(i,j,2),dztb(i,j+1,2))
            endif
            if (kmv(i,j) <= 0) then
               dd = 0.
            endif
            ! thickness of the wall at the u-point
            if (kmv(i,j) /= 0) then
             dd = dd * ( hv + 0.5*(hs(i,j,2) + hs(i,j+1,2)) ) / hv 
             vflux(i,j,k,2) = temp3d_simp(i,j,kk) * dxv(i,j) * dd 
            endif
         enddo
      enddo
   enddo



!  north fold 
do i=4,IMT
 ii=IMT+4-i
! vflux(i,JMT,:,2)=-vflux(ii,JMT-3,:,2)
enddo

#ifdef drifter
! average velocity/transport over surface drifter drogue depth to simulate drifter trajectories
kbot=65 ; ktop=66 ! number of surface layers to integrat over
do i=1,imt
 do j=1,jmt
 
  uint=0. ; vint=0. ; zint=0.
  do k=kbot,ktop
   uint=uint+uflux(i,j,k,2) ! integrated transport
   vint=vint+vflux(i,j,k,2)
   zint=zint+dz(k)          ! total depth of drougued drifter
  enddo
  ! weighted transport for each layer
  do k=kbot,KM
    uflux(i,j,k,2)=uint*dz(k)/zint 
    vflux(i,j,k,2)=vint*dz(k)/zint
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


 deallocate ( temp3d_simp, temp2d_simp )
     
#ifdef tempsalt
 deallocate ( tempb, saltb, rhob, depthb, latb )
#endif

  return
end subroutine readfields



