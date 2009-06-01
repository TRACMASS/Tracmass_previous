
subroutine readfields

  USE mod_param
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_dens
  USE mod_stat
  
  
  IMPLICIT none
  
 INTEGER, PARAMETER ::  NY=145

 REAL*4, ALLOCATABLE, DIMENSION(:,:) :: txy,zxy,uxy,vxy,qxy,pxy,th,zh,uh,vh,qh,ph
 REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: zeta
  
 INTEGER :: i,j,k,n,ii,kk,im,jj,jm,l
 REAL*4 :: pp
 CHARACTER (len=200)                        :: gridFile ,fieldFile,string
 CHARACTER hour(4)*4,month(12)*2,date(31)*2,year(1989:2009)*4
 LOGICAL around
 REAL*8, SAVE :: dxdeg,dydeg,rconst,aa(0:60),bb(0:60),punit
 INTEGER*8, SAVE :: nlon(NY)

data year /'1989','1990','1991','1992','1993','1994','1995','1996','1997','1998','1999',&
                  '2000','2001','2002','2003','2004','2005','2006','2007','2008','2009'/
data month /'01','02','03','04','05','06','07','08','09','10','11','12'/
data date /'01','02','03','04','05','06','07','08','09','10',&
           '11','12','13','14','15','16','17','18','19','20',&
           '21','22','23','24','25','26','27','28','29','30','31'/
data hour /'0000','0600','1200','1800'/


!  print *,'readfield startar',ints
  
  if ( .NOT. ALLOCATED(txy) ) then
    allocate ( txy(IMT,NY),zxy(IMT,NY),uxy(IMT,NY),vxy(IMT,NY),qxy(IMT,NY),pxy(IMT,NY) )
    allocate ( th(IMT,NY),zh(IMT,NY),uh(IMT,NY),vh(IMT,NY),qh(IMT,NY),ph(IMT,NY) )
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


!____________________________ initialise ________________________________________________
if(ints.eq.intstart) then
hs=0.
uflux=0.
vflux=0.
#ifdef tempsalt
tem=0.
sal=0.
rho=0.
#endif
kmt=KM

dxdeg=dx*deg
dydeg=dy*deg
rconst=287.05
punit=1.e-2 ! Pressure units  1.=Pa, 1.e-2=hPa, 1.e.-3=kPa 
iyear=startYear
imon=startMon
iday=startDay
ihour=startHour
print *,'iyear=',iyear,imon,iday,ihour,dxdeg,dydeg

! read the vertical levales
open(12,file=trim(inDataDir)//'topo/model_60lev.txt')
99 format(10x,f12.6,4x,f10.8)
do k=0,KM
 read(12,99) aa(k),bb(k)
! print *,k,aa(k)
enddo
close(12)
!aa=100.*aa   !  [hPa] -->  [Pa]

endif

  ! === swap between datasets ===
!  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
  dztb(:,:,:,1)=dztb(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif


!____ construct format of time to read files _______________________

!print *,'inDataDir',inDataDir

!_____________________________ read ifs fields _________________________________________
!fil=trim(dirgcm)//year(nyear)//'/uvtqzp_'//year(nyear)//month(nmonth)//date(nday)//'.'//trim(hour(nhour))//'.grb'

fieldFile=trim(inDataDir)//'era/'//year(iyear)//'/uvtqzp_'//year(iyear)//month(imon)//date(iday)//'.'//trim(hour(ihour/6+1))//'.grb'

ntime=1000000*iyear+10000*imon+100*iday+ihour
!print *,'hour=',ihour,iday,imon,iyear,fieldFile
!print *,'ntime=',ntime
inquire(file = fieldFile, exist = around )
if(.not.around) then
 print *,'cannot find: ',fieldFile
 stop 39467
endif


! first from u,v,temperature and geopotential from the sprectral gird
!string='wgrib '//trim(fieldFile)//' -o '//trim(inDataDir)//'tmp.bin -d all -bin -nh -V > log.txt'
string='/Volumes/hav1/Applications/wgrib/wgribb '//trim(fieldFile)//' -o '//trim(inDataDir)//'tmp.bin -d all -bin -nh -V > log.txt'

call system(string)
! read
open(14,form='unformatted',file=trim(inDataDir)//'tmp.bin',access='direct',recl=IMT*145*4,convert='little_endian')

kk=0
do k=1,KM
! l=KM+1-k
 l=k
 kk=kk+1
 read(14,rec=kk) txy ! read tempeterature in [K]
 kk=kk+1
 read(14,rec=kk) qxy  ! read specific humidity in [Kg/Kg]
 if(k.eq.1) then
  kk=kk+1
  read(14,rec=kk) zxy ! read geopotential at the orography in [m2/s2]
  kk=kk+1
!print *,'kk=',kk
  read(14,rec=kk) pxy ! read the ln surface pressure 
  do j=1,NY
   ph(:,j)=exp(pxy(:,NY+1-j))   !  [Pa]
!   zh(:,j)=zxy(:,NY+1-j)
  enddo

 endif
 kk=kk+1
!print *,'kk=',kk
 read(14,rec=kk) uxy ! read the zonal velocity in [m/s]
 kk=kk+1
!print *,'kk=',kk
 read(14,rec=kk) vxy !  read the meridional velocity in [m/s]

! reverse the latitudenal order 
  do j=1,NY
   jj=NY+1-j
   th(:,j)=txy(:,jj)
   qh(:,j)=qxy(:,jj)*1.e3 !  [g/Kg]
   uh(:,j)=uxy(:,jj)*( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*ph(:,j) )*punit
   vh(:,j)=vxy(:,jj)*( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*ph(:,j) )*punit
!   uh(:,j)=uxy(:,jj)
!   vh(:,j)=vxy(:,jj)
  enddo

! vh(:,NY)=0.  ! sets velocity to zero at North pole ??????????????????

! A-grid -> C-grid & store in matrixes
 do j=1,JMT
  jj=j+1
  jm=j
  do i=1,IMT
   im=i-1
   if(im.eq.0) im=IMT
   tem  (i,j,l,2)=0.25*(th(i,jj)+th(im,jj)+th(i,jm)+th(im,jm))
   sal  (i,j,l,2)=0.25*(qh(i,jj)+qh(im,jj)+qh(i,jm)+qh(im,jm))
   pp=0.25*(ph(i,jj)+ph(im,jj)+ph(i,jm)+ph(im,jm))
   rho  (i,j,l,2)=0.5*( aa(k)+aa(k-1) + (bb(k)+bb(k-1))*pp )*punit
   dztb (i,j,l,2)= ( aa(k)-aa(k-1) + (bb(k)-bb(k-1))*pp )*punit*dxdy(i,j)

   uflux(i,j,k,2)=0.5*( uh(i,jj)+uh(i ,jm) ) * dydeg
   vflux(i,j,k,2)=0.5*( vh(i,jj)+vh(im,jj) ) * dxdeg*csu(jj)

!   if(k.eq.1) print *,i,j,uflux(i,j,k,2),vflux(i,j,k,2),dztb(i,j,k,2)
  enddo
 enddo
!stop 4067
enddo ! enddo k-loop


close(14) 




!deallocate ( pxy,zxy,uxy,vxy,qxy,pxy  )
!#ifdef tempsalt
!deallocate ( gaus )
!#endif
return
end subroutine readfields
