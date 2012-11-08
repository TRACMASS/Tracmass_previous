#ifdef larval_fish
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
!
! Subroutine of a sedimentation model for the Baltic Sea by Hanna Corell.
! The routine calculates the settling velocity of the particle size modelled.
! The formulas are taken from 
! F. Shepard, Submarine Geology, Harper InternationalEdition, 1948
! p. Nielsen, Coastla bottom boundary layers ..., World Scientific, 1992
!
!      subroutine fishvel(rhof,D,visc,wfish) ! If called from main.F

subroutine fishvel(rhof,temp,dens)  ! If called from loop.F
USE mod_param
USE mod_fish
IMPLICIT none

REAL*8   :: grav,rho,visc,xny,Dm,rhof,temp,dens,a1,a2,a3,a4

grav = 9.81 ! gravity
Dm = fishdiam/1000.0 ! change diameter from mm to m
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

if (fishdiam.gt.0.2) then
 wfish = (-3*xny+SQRT(9*xny**2+grav*fishdiam**2*((rhof-rho)/rho)*&
     &       (a1+a2*fishdiam)))/(a3+a4*fishdiam)
elseif (fishdiam.le.0.2 .and. fishdiam.gt.0.) then
 wfish = -(1.0/18.0)*((rhof-rho)/visc)*grav*(Dm**2) 
elseif (fishdiam.le.0.) then
 print*,'illegal diameter'
 stop 4956
endif

!print*,rho,rhof,temp,visc,wfish
!print*,'wfish=  ',wfish

return
end subroutine fishvel
#endif
