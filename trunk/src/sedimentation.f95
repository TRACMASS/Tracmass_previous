#ifdef sediment
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
!
! Subroutine of a sedimentation model for the Baltic Sea by Hanna Corell.
! The routine calculates the settling velocity of the particle size modelled.
! The formulas are taken from 
! F. Shepard, Submarine Geology, Harper InternationalEdition, 1948
! p. Nielsen, Coastla bottom boundary layers ..., World Scientific, 1992
!
!      subroutine sedvel(rhos,D,visc,wsed) ! If called from main.F

subroutine sedvel(temp,dens)  ! If called from loop.F
USE mod_param
USE mod_sed
USE mod_orbital
IMPLICIT none

REAL*4   :: grav,rho,visc,xny,Dm,temp,dens,a1,a2,a3,a4

grav = 9.81 ! gravity
Dm = partdiam/1000.0 ! change diameter from mm to m
a1=0.003869
a2=0.0248
a3=0.011607
a4=0.07440

! mean density of the Baltic sea
! rho = 1002  ! If called from main.F
! temp = 10.
rho=dens+1000.0  ! If called from loop.F
if(rho.lt.900.) stop 3095
if(rho.gt.1200.) stop 3096

if (temp.lt.2.5) then
 visc=0.001787
elseif ((temp.ge.2.5).and.(temp.lt.7.5)) then 
 visc=0.001519
elseif ((temp.ge.7.5).and.(temp.lt.12.5)) then 
 visc=0.00131
elseif ((temp.ge.12.5).and.(temp.lt.17.5)) then 
 visc=0.001139
elseif ((temp.ge.17.5).and.(temp.lt.22.5)) then 
 visc=0.001002
elseif (temp.ge.22.5)then 
 visc=0.0008904
endif

xny = visc/rho

if (partdiam.gt.0.2) then
 wsed = (-3*xny+SQRT(9*xny**2+grav*partdiam**2*((rhos-rho)/rho)*(a1+a2*partdiam)))/(a3+a4*partdiam)
 if(wsed.gt.0.) stop 4967
elseif (partdiam.le.0.2 .and. partdiam.gt.0.) then
 wsed = -(1.0/18.0)*((rhos-rho)/visc)*grav*(Dm**2) 
elseif (partdiam.le.0.) then
 print*,'illegal diameter'
 stop 4956
endif

!print*,rho,rhos,temp,visc,wsed
!print*,'wsed=  ',wsed

return
end subroutine sedvel

!______________________________________________________________________________________
!_Check if water velocities are large enough for resuspention of sedimentated trajectories
subroutine resusp(res,ib,jb,kb)
USE mod_param
USE mod_sed
USE mod_orbital
USE mod_vel
USE mod_coord
USE mod_grid
IMPLICIT none


REAL botvel,uorb,surfvel,kin   !kin ska bort
INTEGER ib,jb,kb
logical res 

 
surfvel=sqrt( ((vflux(ib,jb,km,1)+vflux(ib,jb-1,km,1))/(2*dx*deg*cst(jb)*dz(km)))**2+ &
              ((uflux(ib,jb,km,1)+uflux(ib-1,jb,km,1))/(2*dy*deg*        dz(km)))**2)
!alpha=cwamp*surfvel
!uorb=orb(41-kb)*alpha
uorb=orb(km+1-kb)*cwamp*surfvel

botvel=sqrt( ((vflux(ib,jb,kb,1)+vflux(ib,jb-1,kb,1))/(2*dx*deg*cst(jb)*dz(kb)))**2+ &
             ((uflux(ib,jb,kb,1)+uflux(ib-1,jb,kb,1))/(2*dy*deg*        dz(kb)))**2)

! 549  format(f9.7 ,'  ',f9.7,' ',f9.7,'  ',f9.7,' ',i8)
!       write(548,549) surfvel,botvel,uorb,wsed,kb


if(botvel+uorb.ge.critvel) then
! print*,'resuspension, botvel+uorb= ',botvel+uorb,' critvel= ',critvel
 res=.true.
else
 res=.false.
endif

return
end subroutine resusp


!______________________________________________________________________________________
! Subroutine to sedimentation model for the Baltic sea. The routine calculates the a 
! meta parameter "orb" that is multiplied with an approximation of the wave amplitude in 
! loop.F, and added to the velocity in the bottom box to simmulate the water movements due
! to short surface waves.

subroutine orbitalv(H)
USE mod_param
USE mod_sed
USE mod_orbital
USE mod_grid
IMPLICIT none

INTEGER k,l
REAL dh,omega,pi,grav,H(km),zk(km),zk49(km),zk50(km),test

dh=0.0
pi = 2.*asin(1.)
omega=2.*pi/twave
grav=9.81
test=0.0 !


! H = vector with the depth at the lower wall in each box in meters,
! vector orb from formula for orbital velocity.

do k=1,km  ! depth --> km levels
 dh=dh+dz(k)
 H(k)=dh
 zk(k)=1./H(k)         ! first guess
 do l=1,50   ! iterativ loop for wave number k
  zk(k)=0.9*zk(k)+0.1*omega**2/(grav*tanh(zk(k)*H(k)))
  if(l.eq.49) then 
   zk49(k)=zk(k)
  endif
  if(l.eq.50)then
   zk50(k)=zk(k)
  endif
 enddo
 orb(k)=omega/sinh(zk(k)*H(k))
enddo

do l=1,km
 test=test+abs(zk49(l)-zk50(l))
!           print*,orb(l),' ',zk(l)
enddo
print*,'sum of total error iterations 49 and 50 ',test

return
end subroutine orbitalv
#endif
