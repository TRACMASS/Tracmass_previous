subroutine readfields

  USE mod_param
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  !  USE mod_dens
  USE mod_stat
  USE mod_tempsalt
  
  IMPLICIT none
    
  INTEGER :: i,j,k,n,ii,kk,im,jj,jm,l
  REAL*8, SAVE ::omtime,cox,coy,om,co,uwe,dl
  REAL*8, SAVE ::ug,u0,gamma,gammag,fcor
  
  call datasetswap
  call updateClock

  ntime=1000000*currYear + 10000*currMon + 100*currDay + currHour


  !____ construct format of time to read files _______________________
  !omtime = 2.d0 *pi * dble(ints)/dble(2)
  omtime = dble(ints) / dble(5) + 1.
  !omtime = 2.d0 * pi * dble(ints)/dble(1000)

  cox    = 0.5d0 + 0.5d0 * dcos(omtime)
  coy    = 0.5d0 + 0.5d0 * dcos(omtime+pi)

  uwe    = -0.4d0
  dl     = dble(ints) * 0.01d0 * pi

  ! Parameters for Nicoletta Fabboni velocities.
  ug     = 0.04 
  u0     = 0.5 
  fcor   = 2. * 2. * pi/(24.*3600.) * cos(45.*pi/180.)
  gamma  = 1./(2.89*24.*3600.)
  gammag = 1./(28.9*24.*3600.)
  omtime = dble(ints-intstart) * dble(ngcm*3600)
  
  ug=1. +0.1*dsin(omtime)
!  vg=1. +0.1*dsin(omtime)


!C-grid & store in matrixes

dzt(:,:,:,2)=10.-9.5*mod(ints,2)

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
! ------------------------------
!         uflux(i,j,k,2) = dy * deg * dz(k) * cox *  &
!              ( dcos( pi*dble(i-1-imt/2)/dble(imt-1) + dl) *  &
!              dsin(-pi*dble(j-1-jmt/2)/dble(jmt-1) ) + uwe + dcos(omtime) )
!         vflux(i,j,k,2) = dx * deg * dz(k) * coy * & 
!              ( dsin( pi*dble(i-1-imt/2)/dble(imt-1) + dl) *  &
!              dcos( pi*dble(j-1-jmt/2)/dble(jmt-1) ) + dsin(omtime) )

! Circular stationary field
! -------------------------
!        uflux(i,j,k,2)=-cox*dy*deg*dz(k)*dble(j-1-jmt/2)/dble(jmt-1) 
!        vflux(i,j,k,2)= cox*dx*deg*dz(k)*dble(i-1-imt/2)/dble(imt-1)

! Spatially constant field which oscilates in time
! ------------------------------------------------
!        uflux(i,j,k,2) = dy * deg*dz(k)*(dcos(omtime)-0.01  +dcos(omtime/pi))
!        vflux(i,j,k,2) = dx * deg*dz(k)*(dsin(omtime)-0.001 +dsin(omtime/pi))
!        uflux(i,j,k,2) = dy * deg*dz(k)*(-0.01  +dcos(omtime/pi))
!        vflux(i,j,k,2) = dx * deg* dz(k)*(-0.001 +dsin(omtime/pi))
   
!         uflux(i,j,k,2) = dy * deg * dz(k) * cox * & 
!              ( -0.05 +omtime/5000.+ dcos(omtime) + 0.5 * dsin(omtime*2.) + &
!              2. * dsin(omtime*10.))
!         vflux(i,j,k,2) = dx * deg * dz(k) *coy * & 
!              ( -0.025 + dsin(omtime) + 0.5 * dcos(omtime*2.) + & 
!              0.5*dcos(omtime/2))

! Nicoletta Fabboni velocities, which have analytical solutions
! -------------------------------------------------------------
!         uflux(i,j,k,2) = dyu(i,j) * dzt(i,j,k,2) * ( ug*dexp(-gammag*omtime) + &
!                                               (u0-ug) * dexp(-gamma*omtime) * cos(fcor*omtime+pi/2.d0) )
!         vflux(i,j,k,2) = dxv(i,j) * dzt(i,j,k,2) * ( -(u0-ug) * dexp(-gamma*omtime) * sin(fcor*omtime+pi/2.d0) )

! One single grid cell
! -------------------------------------------------------------
    if(   i== 2 .and. j==3 ) then
     uflux(i,j,k,2) = -0.3 * dyu(i,j) * dzt(i,j,k,2) 
    elseif(i==3 .and. j==3 ) then
     uflux(i,j,k,2) =  0.1 * dyu(i,j) * dzt(i,j,k,2) 
    else
     uflux(i,j,k,2) = 0.
    endif
    
    if(    i== 3 .and. j==2 ) then
     vflux(i,j,k,2) = 1. * dxv(i,j) * dzt(i,j,k,2) 
    elseif(i==3 .and. j==3 ) then
     vflux(i,j,k,2) =  0.4 * dxv(i,j) * dzt(i,j,k,2)
    else
     vflux(i,j,k,2) = 0.
    endif
    
    print *,i,j,uflux(i,j,k,2),vflux(i,j,k,2)


    !if(j==1 .or. j==JMT) vflux(i,j,k,2)=0.

      end do
   end do
end do ! enddo k-loop


!uflux(:,:,:,2)=uflux(:,:,:,2)*ug
!vflux(:,:,:,2)=vflux(:,:,:,2)/ug


! set velocity on one grid box wall to ensure no divergence
 vflux(3,3,:,:) = vflux(3,2,:,:) + uflux(2,3,:,:) - uflux(3,3,:,:)
 print *,'vflux(3,3,:,:)',vflux(3,3,1,:),dzt(1,1,1,2)

tracers3D(1)%data(:,:,:,nsp) = tem(:,:,:,nsp)
tracers3D(2)%data(:,:,:,nsp) = sal(:,:,:,nsp)
tracers3D(3)%data(:,:,:,nsp) = rho(:,:,:,nsp)

! zero at north and south boundaries 
!vflux(:,0  ,:,:)=0.
!vflux(:,JMT,:,:)=0.


return
end subroutine readfields
