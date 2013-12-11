subroutine readfields

  USE mod_param
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
!  USE mod_dens
  USE mod_stat
  
  
  IMPLICIT none
    
 INTEGER :: i,j,k,n,ii,kk,im,jj,jm,l
 CHARACTER hour(4)*4,month(12)*2,date(31)*2
 REAL*8, SAVE ::omtime,cox,coy,om,co,uwe,dl
 REAL*8, SAVE ::ug,u0,gamma,gammag,fcor

data month /'01','02','03','04','05','06','07','08','09','10','11','12'/
data date /'01','02','03','04','05','06','07','08','09','10',&
           '11','12','13','14','15','16','17','18','19','20',&
           '21','22','23','24','25','26','27','28','29','30','31'/
data hour /'0000','0600','1200','1800'/




!_______________________ update the time counting ________________________________________
 ihour=ihour+1
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

iyear=startYear
imon=startMon
iday=startDay
ihour=startHour
print *,'iyear=',iyear,imon,iday,ihour,dxdeg,dydeg
!print *,'dx=',dx,dxdeg
!print *,'dy=',dy,dydeg

!call coordinat
endif

  ! === swap between datasets ===
!  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif
ntime=1000000*iyear+10000*imon+100*iday+ihour


!____ construct format of time to read files _______________________

!omtime=2.d0*pi*dble(ints)/dble(2)
omtime=dble(ints)/dble(5)+1.
!omtime=2.d0*pi*dble(ints)/dble(1000)
!print *,'omtime',ints,omtime

cox=0.5d0+0.5d0*dcos(omtime)
coy=0.5d0+0.5d0*dcos(omtime+pi)

!cox=0.d0  ! stationary


uwe=-0.4d0
dl=dble(ints)*0.01d0*pi

! Nicoletta Fabboni velocities, which have analytical solutions
ug=0.04  ; u0=0.5 ; fcor=2.*2.*pi/(24.*3600.)*cos(45.*pi/180.)
gamma=1./(2.89*24.*3600.)
gammag=1./(28.9*24.*3600.)
omtime=dble(ints-intstart)*dble(ngcm*3600)

!print *,'fcor=',fcor,omtime,omtime*fcor,cos(omtime*fcor)


!co=0.

!C-grid & store in matrixes

do k=1,KM
 do j=1,JMT
  do i=1,IMT
   im=i-1
   if(im.eq.0) im=IMT
#ifdef tempsalt 
   tem  (i,j,k,2)=20.*float(k)/float(km)
   sal  (i,j,k,2)=30.
   rho  (i,j,k,2)=(28.-20.)*float(km-k)/float(km) +20.
#endif

! time evolving oscilating field
!   uflux(i,j,k,2)=dy*deg*dz(k)*cox*( dcos( pi*dble(i-1-imt/2)/dble(imt-1) + dl) *  &
!                                     dsin(-pi*dble(j-1-jmt/2)/dble(jmt-1) )  + uwe &
!                                     + dcos(omtime)                                )
!   vflux(i,j,k,2)=dx*deg*dz(k)*coy*( dsin( pi*dble(i-1-imt/2)/dble(imt-1) + dl) *  &
!                                     dcos( pi*dble(j-1-jmt/2)/dble(jmt-1) )        &
!                                     + dsin(omtime)                                )
! Circular stationary field
!   uflux(i,j,k,2)=-cox*dy*deg*dz(k)*dble(j-1-jmt/2)/dble(jmt-1) 
!   vflux(i,j,k,2)= cox*dx*deg*dz(k)*dble(i-1-imt/2)/dble(imt-1)

! Spatially constant field which oscilates in time
!   uflux(i,j,k,2)=dy*deg*dz(k)*(dcos(omtime)-0.01  +dcos(omtime/pi))
!   vflux(i,j,k,2)=dx*deg*dz(k)*(dsin(omtime)-0.001 +dsin(omtime/pi))
   
!   uflux(i,j,k,2)=dy*deg*dz(k)*(-0.01  +dcos(omtime/pi))
!   vflux(i,j,k,2)=dx*deg*dz(k)*(-0.001 +dsin(omtime/pi))
   
!   uflux(i,j,k,2)=dy*deg*dz(k)*cox*( -0.05 +omtime/5000.+ dcos(omtime) + 0.5*dsin(omtime*2.) + 2.*dsin(omtime*10.))
!   vflux(i,j,k,2)=dx*deg*dz(k)*coy*( -0.025 + dsin(omtime) + 0.5*dcos(omtime*2.) + 0.5*dcos(omtime/2))

! Interia oscillation velocities, which have analytical solutions for the trajectories
   uflux(i,j,k,2)=dy*dz(k)*( ug*dexp(-gammag*omtime) &
                             + (u0-ug)*dexp(-gamma*omtime)*cos(fcor*omtime) )
   vflux(i,j,k,2)=dx*dz(k)*(  -(u0-ug)*dexp(-gamma*omtime)*sin(fcor*omtime) )
! Munk Gyre velocities, which have analytical solutions for the trajectories (barotropic stream lines)
!   uflux(i,j,k,2)=dy*dz(k)*( ? )
!   vflux(i,j,k,2)=dx*dz(k)*( ? )





  enddo
 enddo
enddo ! enddo k-loop

! zero at north and south boundaries 
!vflux(:,0  ,:,:)=0.
!vflux(:,JMT,:,:)=0.

do j=jmt,1,-10
!print *,(vflux(i,j,2,2),i=1,imt,10)
enddo

!stop 39567

!print *,'uv1=',uflux(95,16,5,1),vflux(95,16,5,1),ints
!print *,'uv2=',uflux(95,16,5,2),vflux(95,16,5,2),ints

!print *,dztb

return
end subroutine readfields

!_______________________________________________________________________

