!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
subroutine readfields

  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_stat
  
#ifdef tempsalt
  USE mod_dens
#endif

  IMPLICIT none
  
#ifdef tempsalt
  REAL*4, ALLOCATABLE, DIMENSION(:) :: tempb, saltb, rhob
  integer kmm
#endif
  
!  integer NLEN,NSNAPS,NSNAPS_2D,NSNAPS_3D,nd,l_tke,nlength
  
!  parameter (NLEN=29557,NSNAPS_2D=1,NSNAPS_3D=0,nd=43200,NSNAPS=1,nlength=2379)
!  parameter (l_tke=0)

  INTEGER, PARAMETER :: NLEN=29557,NSNAPS_2D=1,NSNAPS_3D=0,nd=43200,NSNAPS=1
  INTEGER, PARAMETER :: nlength=2379,l_tke=0
  
  INTEGER :: itt0,year,month,day,hour,minute,second,ittstart,itt,dtts
  INTEGER :: imt0,jmt0,km0,nt0,NLEN0,NSNAPS0
  INTEGER ::  i,j,k,m,kz,ii,ints2,kk,i0
  INTEGER,  SAVE, ALLOCATABLE, DIMENSION(:,:) :: kmu

  REAL*8 :: ird0,ird20,ird30,ird40,stlon,stlat,dxdeg,dydeg
  REAL*4 :: snapd,totsec,ird,ird2,ird3,ird4
  REAL*4, ALLOCATABLE, DIMENSION(:)   :: rd1d_a, rd1d_b, zdzz,dzw,dxt
  REAL*4, ALLOCATABLE, DIMENSION(:)   :: phit,yu,snap1d
  REAL*4, ALLOCATABLE, DIMENSION(:,:) :: rd2d
  
  REAL :: snap2d(imt,jmt) ! ??????????????????

  CHARACTER ofile*20,infile*78,zfile*123,rfile*59,a_exp1*3,a_exp2*2

  LOGICAL around
  
  REAL*8, SAVE :: dxa,dya

  alloCondGrid: if(.not. allocated (snap1d)) then
     allocate ( snap1d(NLEN),rd2d(IMT,JMT) )
     allocate ( rd1d_a(NSNAPS),rd1d_b(NSNAPS) )
     allocate ( zdzz(KM),dzw(0:km),dxt(imt) ) 
     allocate ( phit(jmt),yu(jmt) )
     allocate ( tempb(KM), saltb(KM), rhob(KM) )
  end if alloCondGrid

 if ( .not. allocated (kmu) ) then
     allocate ( kmu(imt,jmt) )
  end if

 
!_______________________ update the time counting ________________________________________
ihour=ihour+6
if(ihour.eq.24) then
 ihour=0
 iday=iday+1
 if(iday.gt.idmax(imon,iyear)) then
  iday=1
  imon=imon+1
  if(imon.eq.13) then
   imon=1
   iyear=iyear+1
   if(iyear.gt.yearmax) iyear=yearmin ! recycle over gcm outputdata
  endif
 endif
endif
ntime=1000000*iyear+10000*imon+100*iday+ihour

!print *,'tiden=',ntime,iyear,imon,iday,ihour,iyear0

!____________________________ initialise ___________________________
if(ints.eq.intstart) then
hs=0.
uflux=0.
vflux=0.
#ifdef tempsalt
tem=0.
sal=0.
rho=0.
#endif
iyear=startyear
imon=startMon
iday=startDay
ihour=0

endif
ntime=1000000*iyear+10000*imon+100*iday+ihour

!_______________________________________________________________________

! swap between datasets

    ! === swap between datasets ===
    hs(:,:,1)=hs(:,:,2)
    uflux(:,:,:,1)=uflux(:,:,:,2)
    vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef explicit_w
    wflux(:,:,:,1)=wflux(:,:,:,2)
#endif
#ifdef tempsalt
    tem(:,:,:,1)=tem(:,:,:,2)
    sal(:,:,:,1)=sal(:,:,:,2)
    rho(:,:,:,1)=rho(:,:,:,2)
#endif


!  print *,'hhhhhhhhhhh',ints
                     
!     === Create filenames for the snap-files to be used ===
ofile='d0000000000.snap1'
write(ofile(2:11),'(i10)') ntime

infile=trim(inDataDir)//ofile
!print *,infile
inquire(file=trim(infile)//'.gz',exist=around)
if(.not.around) then
print *,'This file is missing:',infile,ntime
stop 4555
endif
zfile='gunzip -c '//trim(infile)//'.gz > '//trim(inDataDir)//'tmp/'//trim(outDataFile)
!print *,zfile
!print *,'outDataFile=',outDataFile
!print *,'trim(outDataFile)=',trim(outDataFile)
CALL system(zfile)
rfile=trim(inDataDir)//'tmp/'//trim(outDataFile)
inquire(file=trim(rfile),exist=around)
if(.not.around) stop 4556



!     === open snap file ===
!open(unit=30,file=trim(rfile),status='old',form='unformatted',err=4000,convert='little_endian')
open(unit=30,file=trim(rfile),status='old',form='unformatted',err=4000,convert='big_endian')

! start read header:
read(30) ird
itt0 = ird
read(30) ird
km0 = ird
!print *,'km0',ird
if(km0.ne.km) stop 2844
read(30) ird
!print *,'ird',ird
nt0 = ird
read(30) ird
imt0 = ird
!print *,'imt0',ird
if(imt0.ne.imt) stop 2845
read(30) ird
!print *,'ird',ird
jmt0 = ird
if(jmt0.ne.jmt) stop 2846
read(30) ird
NLEN0 = ird
read(30) ird
NSNAPS0 = ird
read(30) ird
year = ird
read(30) ird
month = ird
read(30) ird
day = ird
read(30) ird
hour = ird
read(30) ird
minute = ird
read(30) ird
second = ird


read(30) ird0,ird20,ird30
dtts = ird0
totsec = ird20
snapd = ird30
!print *,'dtts',ird0,ird20,ird30
read(30) ird0,ird20,ird30,ird40
!dx = ird0 !model grid in cm
!dy = ird20
dxdeg = ird30
dydeg = ird40
!print *,'dxdeg',ird0,ird20,ird30,ird40
read(30,err=2000) ird0,ird20

read(30,err=2000) rd1d_a,rd1d_b
!print *,'rd1d_a,rd1d_b=',rd1d_a,rd1d_b
!do i=1,NSNAPS
! ispvar(i) = rd1d_a(i)
! isplev(i) = rd1d_b(i)
!enddo
read(30,err=2000) rd2d

if(ints.eq.intstart) then

stlon1 = ird0
stlat1 = ird20
!print *,'stlon,stlat=',stlon1,stlat1

!dya=0.005*dydeg
!dxa=0.005*dxdeg


do i=1,imt
 do j=1,jmt
  kmt(i,j) = rd2d(i,j)
 enddo
enddo

 do i=1,imt-1
  do j=1,jmt-1
   kmu(i,j)=min(kmt(i,j),kmt(i+1,j),kmt(i,j+1),kmt(i+1,j+1))
  enddo
 enddo

call coordinat

dya=0.005*dy*deg  
dxa=0.005*dx*deg

endif

!print *,'dx=',dxa,dya,dx,dy,deg


! ssh 

read(30,err=2000) ird

i0=nint(ird)
!print *,'i0',ird,i0
!print *,'snap1d',snap1d

read(30,err=2000) snap1d(1:i0)


ii=0
do j=1,jmt
 do i=1,imt
  hs(i,j,2)=0.
  if(kmt(i,j).ge.1) then
   ii=ii+1
   hs(i,j,2)=0.01*snap1d(ii) 
!   if(kmt(i,j).ne.41) print *,ii,i,j,kmt(i,j),hs(1,j,2)
  endif
 enddo
enddo

!print *,'hs',hs


! ubt
read(30,err=2000) ird
i0=nint(ird)
read(30,err=2000) snap1d(1:i0)
!ii=0
!do j=1,jmt
! do i=1,imt
!  ubt(i,j)=vmask(1)
!  if(kmu(i,j).ge.1) then
!   ii=ii+1
!   ubt(i,j)=snap1d(ii)/hr(i,j)
!  endif
! enddo
!enddo

!print *,'ubt'


! vbt
read(30,err=2000) ird
i0=nint(ird)
read(30,err=2000) snap1d(1:i0)
!ii=0
!do j=1,jmt
! do i=1,imt
!  vbt(i,j)=vmask(1)
!  if(kmu(i,j).ge.1) then
!   ii=ii+1
!   vbt(i,j)=snap1d(ii)/hr(i,j)
!  endif
! enddo
!enddo

!print *,'vbt'


! vad är detta?
do m=4,58
 read(30,err=2000) ird
 i0=nint(ird)
 read(30,err=2000) snap1d(1:i0)
! print*,'2d snapshots',m,i0,snap1d(1)
enddo


! temperature
do k=1,km
 read(30,err=2001) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2001) snap1d(1:i0)
 endif
! print*,'temperature',k,i0,snap1d(1)
#ifdef tempsalt
 kk=km+1-k
 ii=0
 do j=1,jmt
  do i=1,imt
   tem(i,j,kk,2)=0.
   if(kmt(i,j).ge.k) then
    ii=ii+1
    tem(i,j,kk,2)=snap1d(ii)
   endif
  enddo
 enddo
#endif
enddo

do k=1,km
 read(30,err=2002) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2002) snap1d(1:i0)
 endif
! print*,'salinity',k,i0,snap1d(1)
#ifdef tempsalt
 kk=km+1-k
 ii=0
 do j=1,jmt
  do i=1,imt
   sal(i,j,kk,2)=0.
   if(kmt(i,j).ge.k) then
    ii=ii+1
    sal(i,j,kk,2)=snap1d(ii)
   endif
  enddo
 enddo
#endif
enddo

! zonal velocity
do k=1,km
 kk=km+1-k
 read(30,err=2003) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2003) snap1d(1:i0)
 endif
! print*,'u3',k,i0,snap1d(1)
 ii=0
 do j=1,jmt
  do i=1,imt
   snap2d(i,j)=0.
   if(kmu(i,j).ge.k) then
    ii=ii+1
    snap2d(i,j)=snap1d(ii)
   endif
  enddo
 enddo

! u -> transports
 do j=2,jmt
  do i=1,imt-1
    if(kk.ne.km) then
     uflux(i,j,kk,2)=dya*(snap2d(i,j)+snap2d(i,j-1))*dz(kk)
    else
     uflux(i,j,kk,2)=dya*(snap2d(i,j)+snap2d(i,j-1))*( dz(kk)+0.5*(hs(i,j,2)+hs(i+1,j,2)) )
    endif
  enddo
 enddo

enddo

! meridional velocity
do k=1,km
 kk=km+1-k
 read(30,err=2004) ird
 i0=nint(ird)
 if(i0.gt.0) then
  read(30,err=2004) snap1d(1:i0)
 endif
 ii=0
 do j=1,jmt
  do i=1,imt
   snap2d(i,j)=0.
   if(kmu(i,j).ge.k) then
    ii=ii+1
    snap2d(i,j)=dxa*csu(j)*snap1d(ii)
   endif
  enddo
 enddo
! v -> transport
 do j=1,jmt-1
  do i=2,imt
   if(kk.ne.km) then
    vflux(i,j,kk,2)=(snap2d(i,j)+snap2d(i-1,j))*dz(kk)
   else
    vflux(i,j,kk,2)=(snap2d(i,j)+snap2d(i-1,j))*( dz(kk)+0.5*(hs(i,j,2)+hs(i,j+1,2)) )
   endif
  enddo
 enddo

enddo



goto 9000

2001  continue
2002  continue
2003  continue
2004  continue
    
4000  continue
2000  print*,'error 2000'
stop 2000
9000 continue

close(30) 


#ifdef tempsalt
! the density
do i=1,imt
 do j=1,jmt
  if(kmt(i,j).ne.0) then
   kmm=kmt(i,j)
   do k=1,kmm
    kk=km+1-k
    tempb(k)=tem(i,j,kk,2)
!    saltb(k)=(sal(i,j,kk,2)-35.)/1000.
    saltb(k)=sal(i,j,kk,2)
    if(saltb(k).lt.0.) saltb(k)=0.
    enddo
    call statv(tempb,saltb,rhob,kmm)
    do k=1,kmm
     kk=km+1-k
     rho(i,j,kk,2)=rhob(k)
    enddo
   endif
  enddo
 enddo
#endif


!deallocate ( snap1d, rd2d )
!deallocate ( rd1d_a, rd1d_b )
!deallocate ( zdzz,dzw,dxt )
!deallocate ( phit, yu )
!deallocate ( tempb, saltb, rhob )
!deallocate ( kmu )


!print *,'readfield slut',ints
!if(ints.gt.1) stop 49678

return
end subroutine readfields

!________________________________________________________________________________________
      
