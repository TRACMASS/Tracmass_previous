subroutine readfields

  USE mod_param

  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_dens
  USE mod_stat
  
  
  IMPLICIT none
  
 INTEGER, PARAMETER ::  NY=160,NGAUS=35718

 REAL*4, ALLOCATABLE, DIMENSION(:,:) :: snap,psurf
 REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: zeta
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:) :: gaus
#endif
  
 INTEGER :: nhour,i,j,k,n,ii,kk,im,jj,jm

 REAL*8 :: x0,x1,dxx,dxm,dym

 CHARACTER string*91,hour*6

 LOGICAL around
  
 REAL*8, SAVE :: dxdeg,dydeg,rconst,dsigma,pi,radian,radius,rlatt,rlatu,rmax,smax,tmax

 INTEGER*8, SAVE :: nlon(NY)

data nlon /18 , 25, 36, 40, 45, 54, 60, 64, 72, 72, 80, 90, 96,100,108, &
           120,120,128,135,144,144,150,160,160,180,180,180,192,192,200, &
           200,216,216,216,225,225,240,240,240,256,256,256,256,288,288, &
           288,288,288,288,288,288,288,300,300,300,300,320,320,320,320, &
           320,320,320,320,320,320,320,320,320,320,320,320,320,320,320, &
           320,320,320,320,320,320,320,320,320,320,320,320,320,320,320, &
           320,320,320,320,320,320,320,320,320,320,320,320,320,320,300, &
           300,300,300,288,288,288,288,288,288,288,288,288,256,256,256, &
           256,240,240,240,225,225,216,216,216,200,200,192,192,180,180, &
           180,160,160,150,144,144,135,128,120,120,108,100, 96, 90, 80, &
            72, 72, 64, 60, 54, 45, 40, 36, 25, 18/

!  print *,'readfield startar',ints
  
  if ( .NOT. ALLOCATED(snap) ) then
     allocate ( snap(IMT,NY), zeta(IMT,NY,KM), psurf(IMT,0:JMT) )
#ifdef tempsalt
     allocate ( gaus(NGAUS) )
#endif
  end if

 
!_______________________ update the time counting ________________________________________
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
ntime=1000000*iyear+10000*imon+100*iday+ihour

!____________________________ initialise ________________________________________________
if(ints.eq.intstart) then
hs=0.
u=0.
v=0.
#ifdef tempsalt
tem=273.15
sal=0.
rho=0.
#endif

! Earth constants
radius = 6371229.d0 ! earth radius in metre
pi = 2.d0 * dasin(1.d0)
radian = pi/180.d0
deg=radius*radian ! ~ 111000 metre
grav=9.81 ! m/s2tyear=365.25d0 * 24.d0 * 3600.d0
dx = 360./float(IMT)
dxdeg=dx*deg
rconst=287.05  ! The gas constant J/(K kg)
dsigma=1.d0/real(KM) !?????????????????????????? inte sgimakoordinater utan etha

stlon1=0. 

! geopotential height in km
rmin=0.d0
rmax=200.d0
dr=(rmax-rmin)/dble(MR-1)
! temperature
tmin=-100.d0
tmax=50.d0 
dtemp=(tmax-tmin)/dble(MR-1)
! specific humidity
smin= 0.d0
smax=25.d0
dsalt=(smax-smin)/dble(MR-1)

! read in the guassian latitude coordinates
open(12,file=directory//'/topo/ifsn256lat.txt')
99 format(39x,f9.5)
do j=1,JMT-1
read(12,99) x0
phi(j)=-x0
print 99,phi(j)
enddo
phi(0  )=-90. ! add the north pole as a fictive grid point
phi(JMT)= 90. ! add the south pole as a fictive grid point

do j=1,JMT
dyt(j)=(phi(j)-phi(j-1))*deg
print *,j,phi(j),phi(j-1),dyt(j)
enddo

do j=1,JMT
 rlatt=0.5*(phi(j)+phi(j-1))
 rlatu=phi(j)
 csu(j)=dcos(rlatu*radian)
 cst(j)=dcos(rlatt*radian)
enddo

do i=1,IMT
 do j=1,JMT
  dxdy(i,j)=dx*deg*cst(j)*dyt(j)
 enddo
enddo

kmt=KM  ! all grid cells are active in an AGCM




endif
stop 3057

!____________________________ swap between datasets _____________________________________

hs(:,:,1)=hs(:,:,2)
u(:,:,:,1)=u(:,:,:,2)
v(:,:,:,1)=v(:,:,:,2)
#ifdef tempsalt 
tem(:,:,:,1)=tem(:,:,:,2)
sal(:,:,:,1)=sal(:,:,:,2)
rho(:,:,:,1)=rho(:,:,:,2)
#endif

!____ construct format of time to read files _______________________
nhour=(ints-1)*24
      open(45,file='tmp.hour')
 81   format(6i1)
 82   format(4i1,i2)
 83   format(3i1,i3)
 84   format(2i1,i4)
 85   format(1i1,i5)
 86   format(    i6)

 89   format(3a)

rewind(45)
if(nhour.ge.100000) then
write(45,86) nhour
elseif(nhour.ge.10000) then
write(45,85) 0,nhour
elseif(nhour.ge.1000) then
write(45,84) 0,0,nhour
elseif(nhour.ge.100) then
write(45,83) 0,0,0,nhour
elseif(nhour.ge.10) then
write(45,82) 0,0,0,0,nhour
elseif(nhour.lt.10) then
write(45,81) 0,0,0,0,0,nhour
endif
rewind(45)
read (45,89) hour
!print *,'hour=',hour,nhour,' day=',nhour/24

!_____________________________ read ifs fields ___________________________________________
! first from u,v,temperature and geopotential from the sprectral gird
string='cdo dv2uvl /Users/doos/data/ifs12/spectral/ICMSHPLewhz+'//hour//' tmp'
!print *,string
call system(string)
call system('cdo sp2gpl tmp tmp.gaus > log1.txt')
call system('wgrib tmp.gaus -o tmp.bin -d all -bin -nh -V > log2.txt')
! read
!open(14,form='unformatted',file='tmp.bin',access='direct',recl=IMT*NY*4)
open(14,form='unformatted',file='tmp.bin',access='direct',recl=IMT*NY*4,convert='little_endian')

! read tempeterature in K
do kk=1,km
!print *,kk
 k=KM+1-kk
 read(14,rec=kk) snap
 do i=1,IMT
  do j=1,NY
   jj=NY+1-j
   tem(i,j,k,2)=snap(i,jj)
  enddo
 enddo
enddo
! read the geopotential
do kk=1,km
!print *,kk
 k=KM+1-kk
 read(14,rec=kk+km) snap
 do i=1,IMT
  do j=1,NY
   jj=ny+1-j
   rho(i,j,k,2)=snap(i,jj)  ! rho is now the geopotential
  enddo
 enddo
enddo
! read the zonal wind & reverse the latitude order & A-grid -> C-grid
do kk=1,km
!print *,kk,kk+2*km
 k=KM+1-kk
 read(14,rec=kk+2*km) snap
 do i=1,IMT
  do j=1,JMT
   jj=ny+1-j
   if(j.eq.1) then
    u(i,j,k,2)=0.5*snap(i,jj) 
   elseif(j.eq.JMT) then
    u(i,j,k,2)=0.5*snap(i,jj-1) 
   else
    u(i,j,k,2)=0.5*(snap(i,jj-1)+snap(i,jj))   
   endif
  enddo
 enddo
enddo


! read the meridional wind & reverse the latitude order & A-grid -> C-grid
do kk=1,km
!print *,kk,kk+3*km
 k=KM+1-kk
 read(14,rec=kk+3*km) snap
 do i=1,IMT
 im=i-1
 if(im.eq.0) im=IMT
  do j=1,NY
   jj=ny+1-j
   v(i,j,k,2)=0.5*(snap(i,jj)+snap(im,jj))  
  enddo
 enddo
enddo

close(14) 

#ifdef tempsalt
! secondly read the specific humidity from the gausian grid

string='wgrib /Users/doos/data/ifs12/gaus/ICMGGewhz+'//hour//' -o tmp.bin -d all -bin -nh -V > log3.txt'
!print *,string

call system(string)

open(14,form='unformatted',file='tmp.bin',access='direct',recl=NGAUS*4,convert='little_endian')

do k=1,KM
read(14,rec=29+k) gaus

n=0
do j=1,NY
 dxx=dxdeg*dble(IMT)/dble(nlon(j))
 do ii=1,nlon(j)
  n=n+1
  x0=dxdeg*dble(ii-1)
  x1=dxx*dble(ii)
  do i=1,IMT
! units are here projected on a regular grid and converted into g/kg of specific humidity
   if(x0.lt.dxdeg*dble(i-1).and.dxdeg*dble(i-1).le.x1) sal(i,j,k,2)=gaus(n)*1.e3
  enddo
  enddo
enddo

enddo

close(14)

#endif

!___________________ Calculates the surface pressure at the orography ____________________
! air density=pressure/(R*temperature)
do i=1,IMT
 do j=1,NY
! surface pressure= p0*(1+geopotential/R*temp)
  psurf(i,j)=1.e3*(1. + rho(i,j,KM,2)/(rconst*tem(i,j,KM,2)) )
 enddo
  psurf(i,0  )=psurf(i,1 )
  psurf(i,JMT)=psurf(i,NY)
enddo

!__________________ from speed to "mass transport"  m/s -> kg/s ___________________________
! This is not in m2*Pa/s but in m2*hPa....this should probably be changed!
! u
 do i=1,IMT   
 im=i-1
 if(im.eq.0) im=IMT
  do j=1,JMT
   jm=j-1
   hs(i,j,2)=0.25*(psurf(i,j)+psurf(im,j)+psurf(i,jm)+psurf(im,jm))
   dxm=dxdeg*csu(j)*dsigma*0.5*(psurf(i,j)+psurf(im,j))/grav
   dym=dydeg       *dsigma*0.5*(psurf(i,j)+psurf(i,jm))/grav
   do k=1,km
    u(i,j,k,2)=u(i,j,k,2)*dxm
    v(i,j,k,2)=v(i,j,k,2)*dym
   enddo

  enddo
 enddo


deallocate ( snap, zeta, psurf )
#ifdef tempsalt
deallocate ( gaus )
#endif
return
end subroutine readfields

!_______________________________________________________________________
      
