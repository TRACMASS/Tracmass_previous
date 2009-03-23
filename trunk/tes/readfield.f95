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
    
 INTEGER :: i,j,k,n,ii,kk,im,jj,jm,l
 CHARACTER hour(4)*4,month(12)*2,date(31)*2,year(1989:2009)*4
 REAL*8, SAVE :: dxdeg,dydeg,psi0,t0,omtime,d0,d1,d2,d3,pi,b,c,om,co

data year /'1989','1990','1991','1992','1993','1994','1995','1996','1997','1998','1999',&
                  '2000','2001','2002','2003','2004','2005','2006','2007','2008','2009'/
data month /'01','02','03','04','05','06','07','08','09','10','11','12'/
data date /'01','02','03','04','05','06','07','08','09','10',&
           '11','12','13','14','15','16','17','18','19','20',&
           '21','22','23','24','25','26','27','28','29','30','31'/
data hour /'0000','0600','1200','1800'/




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
hs=0. ; uflux=0. ; vflux=0.
#ifdef tempsalt
tem=0. ; sal=0. ; rho=0.
#endif
kmt=KM

dxdeg=dx*deg
dydeg=dy*deg
iyear=startYear
imon=startMon
iday=startDay
ihour=startHour
print *,'iyear=',iyear,imon,iday,ihour,dxdeg,dydeg
print *,'dx=',dx,dxdeg
print *,'dy=',dy,dydeg

call coordinat
endif

  ! === swap between datasets ===
!  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
  dztb(:,:,1)=dztb(:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif
ntime=1000000*iyear+10000*imon+100*iday+ihour


!____ construct format of time to read files _______________________

pi = 2.*asin(1.)
psi0=100.e6/200.
t0=pi*float((imt-1)*(jmt-1))*dx*dy*deg**2/(8.*psi0)
!time=float(ns-1)/float(nstot)*365.25*3600.*24.
om=2.*pi/(365.25*24.*3600.)
!omtime=2.*pi*float(ns-1)/float(nstot)
omtime=2.*pi*float(mod(ints-1,10))/float(10)
!print *,ints,mod(ints-1,10)

b=float(jmt-3)/2.
b=4.*b
c=b/5.
co=2.*pi/t0*float(jmt-1)/float(imt-1)
d1=-8.*psi0/float(jmt-1)**2
d2= 8.*psi0/float(imt-1)**2
d3=float(jmt-3)*2./5.*cos(omtime)


!C-grid & store in matrixes

do k=1,KM
 do j=1,JMT
  do i=1,IMT
   im=i-1
   if(im.eq.0) im=IMT
   tem  (i,j,k,2)=20.*float(k)/float(km)
   sal  (i,j,k,2)=30.
   rho  (i,j,k,2)=(28.-20.)*float(km-k)/float(km) +20.
   dztb (i,j,2)=dz(k)*dxdy(i,j)
   uflux(i,j,k,2)=dy*deg*dz(k)*( co*( 0.5+float(jmt-1)/2.-float(j))*dy*deg -dx*deg*c*om*cos(omtime) )
   vflux(i,j,k,2)=dx*deg*dz(k)*( co*(-0.5-float(imt-1)/2.+float(i))*dx*deg -dy*deg*c*om*cos(omtime) )
  enddo
 enddo
enddo ! enddo k-loop

return
end subroutine readfields

!_______________________________________________________________________

