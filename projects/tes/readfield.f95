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
    
 INTEGER :: i,j,k,n,ii,kk,im,jj,jm,l,isort
 CHARACTER hour(4)*4,month(12)*2,date(31)*2
 REAL*8, SAVE :: omtime,omtimed,cox,coy,om,co,uwe,dl
 REAL*8, SAVE ::ug,u0,gamma,gammag,fcor

data month /'01','02','03','04','05','06','07','08','09','10','11','12'/
data date /'01','02','03','04','05','06','07','08','09','10',&
           '11','12','13','14','15','16','17','18','19','20',&
           '21','22','23','24','25','26','27','28','29','30','31'/
data hour /'0000','0600','1200','1800'/
!_________________________________________________________________________________________
! Different possible analytical velocity fields can be chosen
! 1=time evolving oscilating field
! 2=Circular stationary field
! 3=Spatially constant field which oscilates in time
! 9=Stommel f-plane Gyre velocities
! 10=Stommel beta-plane Gyre velocities
! 5=Munk Gyre velocities
! 4=Interia oscillation velocities, which have analytical solutions for the trajectories
! 6=random generated velocities
! 7=Interia oscillation + random generated velocities
! 8=Interia oscillation, which decays and the grows

isort=4


!_______________________ update the time counting ________________________________________
 ihour=ihour+1
 if(ihour==24) then
  ihour=0
  iday=iday+1
  if(iday.gt.idmax(imon,iyear)) then
   iday=1
   imon=imon+1
   if(imon==13) then
    imon=1
    iyear=iyear+1
    if(iyear.gt.yearmax) iyear=yearmin ! recycle over gcm outputdata
   endif
  endif
 endif

!____________________________ initialise ________________________________________________
if(ints==intstart) then
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

if(isort==1) then ! 
!omtime=2.d0*pi*dble(ints)/dble(2)
!omtime=dble(ints)/dble(5)+1.
!omtime=2.d0*pi*dble(ints)/dble(1000)
!print *,'omtime',ints,omtime
omtime=dble(ints-intstart)*dble(ngcm*3600)

cox=0.5d0+0.5d0*dcos(omtime)
coy=0.5d0+0.5d0*dcos(omtime+pi)

!cox=0.d0  ! stationary


!uwe=-0.4d0
!dl=dble(ints)*0.01d0*pi

elseif(isort==2) then ! 
ug=0.04  ; u0=0.5 ; fcor=2.*2.*pi/(24.*3600.)*cos(45.*pi/180.)
omtime=dble(ints-intstart)*dble(ngcm*3600)


elseif(isort==4 .or.isort==7 .or. isort==8) then ! Interia oscillations
ug=0.04  ; u0=0.5 ; fcor=2.*2.*pi/(24.*3600.)*cos(45.*pi/180.)
gamma=1./(2.89*24.*3600.)
gammag=1./(28.9*24.*3600.)
omtime=dble(ints-intstart)*dble(ngcm*3600)
!omtimed=dble(ints-intstart)*dble(ngcm*3600)
!if(ints>=200) omtimed=dble(-ints+400-intstart)*dble(ngcm*3600)
!!if(ints>=400) omtimed=dble( ints+600-intstart)*dble(ngcm*3600)
omtimed =0.5-0.5*cos(float(ints-intstart)/200.*pi)
omtimed =omtimed*200*dble(ngcm*3600)
!print *,ints, omtime ,omtimed
elseif(isort==6) then ! random velocitites
ug=0.5
elseif(isort==9 .or. isort==10) then ! random velocitites
ug=0.001
endif

!co=0.

!C-grid & store in matrixes

do k=1,KM
 do j=1,JMT
  do i=1,IMT
   im=i-1
   if(im==0) im=IMT
#ifdef tempsalt 
   tem  (i,j,k,2)=20.*float(k)/float(km)
   sal  (i,j,k,2)=30.
   rho  (i,j,k,2)=(28.-20.)*float(km-k)/float(km) +20.
#endif

if(isort==1) then ! time evolving oscilating field
   uflux(i,j,k,2)=dy*deg*dz(k)*cox*( dcos( pi*dble(i-1-imt/2)/dble(imt-1) + dl) *  &
                                     dsin(-pi*dble(j-1-jmt/2)/dble(jmt-1) )  + uwe &
                                     + dcos(omtime)                                )
   vflux(i,j,k,2)=dx*deg*dz(k)*coy*( dsin( pi*dble(i-1-imt/2)/dble(imt-1) + dl) *  &
                                     dcos( pi*dble(j-1-jmt/2)/dble(jmt-1) )        &
                                     + dsin(omtime)                                )
elseif(isort==2) then ! Circular stationary field
!   uflux(i,j,k,2)=-cox*dy*deg*dz(k)*dble(j-1-jmt/2)/dble(jmt-1) 
!   vflux(i,j,k,2)= cox*dx*deg*dz(k)*dble(i-1-imt/2)/dble(imt-1)

   uflux(i,j,k,2)=dy*dz(k)*( ug+ (u0-ug)*cos(fcor*omtime) )
   vflux(i,j,k,2)=dx*dz(k)*(    -(u0-ug)*sin(fcor*omtime) )

elseif(isort==3) then ! Spatially constant field which oscilates in time
!   uflux(i,j,k,2)=dy*deg*dz(k)*(dcos(omtime)-0.01  +dcos(omtime/pi))
!   vflux(i,j,k,2)=dx*deg*dz(k)*(dsin(omtime)-0.001 +dsin(omtime/pi))
   
!   uflux(i,j,k,2)=dy*deg*dz(k)*(-0.01  +dcos(omtime/pi))
!   vflux(i,j,k,2)=dx*deg*dz(k)*(-0.001 +dsin(omtime/pi))
   
!   uflux(i,j,k,2)=dy*deg*dz(k)*cox*( -0.05 +omtime/5000.+ dcos(omtime) + 0.5*dsin(omtime*2.) + 2.*dsin(omtime*10.))
!   vflux(i,j,k,2)=dx*deg*dz(k)*coy*( -0.025 + dsin(omtime) + 0.5*dcos(omtime*2.) + 0.5*dcos(omtime/2))

elseif(isort==4) then ! Interia oscillation
   uflux(i,j,k,2)=dy*dz(k)*( ug*dexp(-gammag*omtime) &
                             + (u0-ug)*dexp(-gamma*omtime)*cos(fcor*omtime+pi/2.d0) )
   vflux(i,j,k,2)=dx*dz(k)*(  -(u0-ug)*dexp(-gamma*omtime)*sin(fcor*omtime+pi/2.d0) )
elseif(isort==5) then ! Munk Gyre velocities
!   uflux(i,j,k,2)=dy*dz(k)*( ? )
!   vflux(i,j,k,2)=dx*dz(k)*( ? )
elseif(isort==6) then ! 
   uflux(i,j,k,2)=dy*dz(k)*ug*(rand()-0.5+0.05)
   vflux(i,j,k,2)=dx*dz(k)*ug*(rand()-0.5+0.005)
elseif(isort==7) then ! 
!   uflux(i,j,k,2)=dy*dz(k)*( ug*dexp(-gammag*omtime) &
!                             + (u0-ug)*dexp(-gamma*omtime)*cos(fcor*omtime) + 0.29*(rand()-0.5) )
!   vflux(i,j,k,2)=dx*dz(k)*(  -(u0-ug)*dexp(-gamma*omtime)*sin(fcor*omtime) + 0.29*(rand()-0.5) )
   uflux(i,j,k,2)=dy*dz(k)*( ug*dexp(-gammag*omtimed) &
                             + (u0-ug)*dexp(-gamma*omtimed)*cos(fcor*omtime) + 0.29*(rand()-0.5) )
!   vflux(i,j,k,2)=dx*dz(k)*(  -(u0-ug)*dexp(-gamma*omtimed)*sin(fcor*omtime) + 0.29*(rand()-0.5) )
   vflux(i,j,k,2)=dx*dz(k)*( ug/10.*dexp(-gammag*omtimed) &
                              -(u0-ug)*dexp(-gamma*omtimed)*sin(fcor*omtime) + 0.29*(rand()-0.5) )
elseif(isort==8) then ! 
   uflux(i,j,k,2)=dy*dz(k)*( ug*dexp(-gammag*omtimed) &
                             + (u0-ug)*dexp(-gamma*omtimed)*cos(fcor*omtime) )
   vflux(i,j,k,2)=dx*dz(k)*(  -(u0-ug)*dexp(-gamma*omtimed)*sin(fcor*omtime) )
elseif(isort==9) then ! Stommel f-plane gyre
   uflux(i,j,k,2)=-dy*dz(k)*pi*ug*cos(pi*dble(j-1)/dble(jmt-1)) &
                &  *( 1. - (1.-dexp(-1.d0))/(dexp(1.d0)-dexp(-1.d0))*dexp( dble(i-1)/dble(IMT-1)) &
                &          -  (dexp(1.d0)-1.d0)/(dexp(1.d0)-dexp(-1.d0))*dexp(-dble(i-1)/dble(IMT-1)) )
                     
   vflux(i,j,k,2)=-dx*dz(k)*ug*sin(pi*dble(j-1)/dble(jmt-1)) &
                &        *( (1-dexp(-1.d0))/(dexp(1.d0)-dexp(-1.d0))*dexp( dble(i-1)/dble(IMT-1)) &
                &           - (dexp(1.d0)-1.d0)/(dexp(1.d0)-dexp(-1.d0))*dexp(-dble(i-1)/dble(IMT-1)) )
elseif(isort==10) then ! Stommel beta-plane gyre
  uflux(i,j,k,2)=-dy*dz(k)*pi*ug*cos(pi*dble(j-1)/dble(jmt-1)) * &
     & ( 1. - dexp( -1.d0/0.02d0*dble(i-1)/dble(IMT-1)) - dble(i-1)/dble(IMT-1) )
  vflux(i,j,k,2)=dx*dz(k)*ug*sin(pi*dble(j-1)/dble(jmt-1))  * &
                &        (  1.d0/0.02d0*dexp( -1.d0/0.02d0*dble(i-1)/dble(IMT-1)) - 1. )

endif

if(j==1 .or. j==JMT) vflux(i,j,k,2)=0.



  enddo
 enddo
enddo ! enddo k-loop

!print *,(uflux(i,1*JMT/4,KM/2,2),i=1,IMT)
!print *,(uflux(3*IMT/4,j,KM/2,2),j=1,JMT)
!stop 395

! zero at north and south boundaries 
!vflux(:,0  ,:,:)=0.
!vflux(:,JMT,:,:)=0.

!do j=jmt,1,-10
!print *,(uflux(i,j,2,2),i=1,imt,10)
!enddo

!stop 39567

!print *,'uv1=',uflux(95,16,5,1),vflux(95,16,5,1),ints
!print *,'uv2=',uflux(95,16,5,2),vflux(95,16,5,2),ints

!print *,dztb

return
end subroutine readfields

!_______________________________________________________________________

