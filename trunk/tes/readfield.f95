!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

SUBROUTINE readfields

IMPLICIT none
 
#include "../param.h"

common/vel/u(imt,0:JMAX,km,NST),v(imt,0:JMAX,km,NST),hs(imt,JMAX,NST),w(0:km),ff
REAL u,v,hs
REAL*8 w,ff

common /coord/dx,dy,deg,stlon1,stlat1,csu(jmt),cst(jmt),zw(0:km)
REAL*8        dx,dy,deg,stlon1,stlat1,csu,cst,zw

common /grid/dxdy(imt,jmt),dztb(imt,jmt,KD),dz(km),rmin,dr,tmin,dtemp,smin,dsalt,kmt(imt,jmt)
REAL*8 dxdy,dztb,dz,rmin,dr,tmin,dtemp,smin,dsalt
INTEGER kmt


#ifdef tempsalt
common/dens/ tem(imt,JMAX,km,NST),sal(imt,JMAX,km,NST),rho(imt,JMAX,km,NST)
REAL tem,sal,rho
#endif

#ifdef time
common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
integer    ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
#endif

common/namn/name,namep,directory
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=27) :: directory

INTEGER i,j,k,kk,im,jm

REAL :: dxyz(km),dzt(imt,jmt,km),dzu(imt,jmt,km),dzv(imt,jmt,km)

REAL :: d1,d2,pi,psi0,omtime,a,b,c,co,d3,t0,om

CHARACTER :: hour*3,month*2,day*4

LOGICAL around

#if defined sigma
REAL dsig(KM),depth(IMT,JMT)
#endif

save dxyz

if(mod(ints,365).eq.0) then
print *,'psi written for ints=',ints
call writepsi
endif

! swap between time steps
hs(:,:,1)=hs(:,:,2)
u(:,:,:,1)=u(:,:,:,2)
v(:,:,:,1)=v(:,:,:,2)
#ifdef tempsalt 
tem(:,:,:,1)=tem(:,:,:,2)
sal(:,:,:,1)=sal(:,:,:,2)
rho(:,:,:,1)=rho(:,:,:,2)
#endif

#if defined sigma
dsig=1./real(KM)
! open(86,file='/Volumes/sjo5/data/sig/depth',form='unformatted')
! read(86) depth
! close(86)
do j=1,JMT
 do i=1,IMT
  depth(i,j)=1000.-500.*real(j-1)/real(JMT-1)
!  depth(i,j)=100.
 enddo
enddo
#endif


! initialise
if(ints.eq.intstart) then
 call coordinat
 zw=0.
 hs=0.
 kmt=km
 dzt=0.
 u=0.
 v=0.
#ifdef tempsalt
 tem=0.
 sal=0.
 rho=0.
#endif
endif

! transports
! grid box volumes in m3
! grid box depths
do i=1,imt
 do j=1,jmt
  do k=1,km

#if defined sigma
   dzt(i,j,k)=dsig(k)*(depth(i,j)+hs(i,j,2))
   dztb(i,j,k)=dzt(i,j,k)
!if(i.eq.1 .and. j.eq.1) print *,k,dzt(i,j,k)
#else
   dzt(i,j,k)=dz(k)
   if(k.eq.km) dzt(i,j,k)=dzt(i,j,k)+hs(i,j,2)
#endif

  enddo
 enddo
enddo

!print *,(dztb(9,j,1),j=jmt,1,-1)
!stop 5976

do i=1,imt
 im=i-1
 if(im.eq.0) im=IMT
 do j=1,jmt
 jm=j-1
 if(jm.eq.0) jm=1
  do k=1,km

   dzu(i,j,k)=0.5*(dzt(i,j,k)+dzt(im,j,k))
   dzv(i,j,k)=0.5*(dzt(i,j,k)+dzt(i,jm,k))

  enddo
 enddo
enddo


do i=1,imt-1
! dxzv(i,:,:)=0.5*(dzt(i,:,:)+dzt(i+1,:,:))
enddo

pi = 2.*asin(1.)
psi0=100.e6/200.
t0=pi*float((imt-1)*(jmt-1))*dx*dy*deg**2/(8.*psi0)
!time=float(ns-1)/float(nstot)*365.25*3600.*24.
om=2.*pi/(365.25*24.*3600.)
!omtime=2.*pi*float(ns-1)/float(nstot)
omtime=2.*pi*float(mod(ints-1,10))/float(10)
!print *,ints,mod(ints-1,10)
!a=float(imt-3)/2.
b=float(jmt-3)/2.
!a=a/1.5
b=4.*b
c=b/5.
!c=0.
co=2.*pi/t0*float(jmt-1)/float(imt-1)

!print *,co,t0,pi,a,b,c

d1=-8.*psi0/float(jmt-1)**2
d2= 8.*psi0/float(imt-1)**2
d3=c*cos(omtime)*deg
      
do j = 1,jmt
 do i = 1,imt
  do k = 1,km
!stationary
!   u(i,j,k,2)=dy*deg*dz(k)*ff*(d1*(float(j)-0.5-float(jmt-1)/2.)/deg)
!   v(i,j,k,2)=dx*deg*dz(k)*ff*(d2*(float(i)-0.5-float(imt-1)/2.)/deg)
! time dependent
!u(i,j,k,2)=dy*deg*dz(k)*ff*(co*(0.5+float(jmt-1)/2.-float(j))*dy*deg+dx*d3)
!v(i,j,k,2)=dx*deg*dz(k)*ff*(co*(float(i)-0.5-float(imt-1)/2.-c*cos(omtime))*dy*deg)

!u(i,j,k,2)=dy*deg*dz(k)*ff*( co*( 0.5+float(jmt-1)/2.-float(j))*dy*deg -dx*deg*c*om*cos(omtime) )
!v(i,j,k,2)=dx*deg*dz(k)*ff*( co*(-0.5-float(imt-1)/2.+float(i))*dx*deg -dy*deg*c*om*cos(omtime) )

!u(i,j,k,2)=0.

u(i,j,k,2)=0.
v(i,j,k,2)=-0.5

u(i,j,k,2)=u(i,j,k,2)*dy*deg*ff*dzu(i,j,k)
v(i,j,k,2)=v(i,j,k,2)*dx*deg*ff*dzv(i,j,k)

#ifdef tempsalt
   tem(i,j,k,2)=20.*float(k)/float(km)
   sal(i,j,k,2)=30.
   rho(i,j,k,2)=(28.-20.)*float(km-k)/float(km) +20.
#endif

  enddo
 enddo
enddo


! Boundary conditions
do k = 1,km
 do j = 1,JMAX
  u(1  ,j,k,2)=0.
  u(IMT,j,k,2)=0.
 enddo
 do i = 1,imt
  v(i,0   ,k,2)=0.
  v(i,JMAX,k,2)=0.
 enddo
enddo

!stop 45956
return
end subroutine readfields


      
