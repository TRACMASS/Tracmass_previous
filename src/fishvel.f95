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
!      by Kate Hedstrom
! Modified to add egg density as a function of age group
! Age group is based on stage
! Stage is a function of age and temperature
! CMP 11/13/12

subroutine fishvel(temp,dens)  ! If called from loop.F
USE mod_param
USE mod_fish
USE mod_loopvars
IMPLICIT none

REAL*8   :: grav,rho,visc,Dm,rhof,temp,dens
REAL*8   :: age, egg_hatch, yolk_len, length
REAL*8   :: est, group, egg_sg, hatch_hrs, hatch_len
REAL*8   :: ec1,ec2,ec3,ec4,ec5,ec6,ec7,ec8,ec9,ec10
REAL*8   :: ec11,ec12,ec13,ec14,ec15,ec16,ec17,ec18,ec19,ec20
REAL*8   :: ec21,ec22,ec23
REAL*8   :: g_day, g_dt
INTEGER  :: istage
!INTEGER, intent(in)  :: ntrac

age = fish(ntrac,i_age)
egg_hatch = fish(ntrac,i_hatchtime)
yolk_len = fish(ntrac,i_hatchlength)
length = fish(ntrac,i_length)

!--------------DEVELOPMENT / GROWTH-----------------------
!EGGS
If (stage(ntrac) .eq. f_egg)then
ec1 = 2.835455d-1 !a
ec2 = -4.766954d-3 !a2
ec3 = 1.436144d-1 !T
ec4 = -4.228248d-2 !T2
ec5 = 3.565329d-2 !Ta
ec6 = -1.362398d-3 !Ta2
ec7 = 4.305655d-5 !a3
ec8 = 2.173404d-5 !Ta3
ec9 = -2.284878d-7 !a4
ec10 = -1.691007d-07 !Ta4
ec11 = 7.793980d-10 !a5
ec12 = 7.520298d-10 !Ta5
ec13 = -1.772402d-12 !a6
ec14 = -2.054619d-12 !Ta6
ec15 = 2.679053d-15 !a7
ec16 = 3.517559d-15 !Ta7
ec17 = -2.581572d-18 !a8
ec18 = -3.681552d-18 !Ta8
ec19 = 1.430693d-21 !a9
ec20 = 2.153112d-21 !Ta9
ec21 = -3.457755d-25 !a10
ec22 = -5.387357d-25 !Ta10
ec23 = -3.659724d-1 !intercept

age = age + dtmin/3600.d0
fish(ntrac,i_age) = age
est = ec1*age + ec2*age**2 + ec3*temp + ec4*temp**2 + ec5*temp*age + &
    & ec6*temp*age**2 + ec7*age**3 + ec8*temp*age**3 + ec9*age**4 + &
    & ec10*temp*age**4 + ec11*age**5 + ec12*temp*age**5 + ec13*age**6 + &
    & ec14*temp*age**6 + ec15*age**7 + ec16*temp*age**7 + ec17*age**8 + &
    & ec18*temp*age**8 + ec19*age**9 + ec20*temp*age**9 + ec21*age**10 + &
    & ec22*temp*age**9 + ec23
istage = int(est)

If (istage .le. 6)then
   group = 1.
Elseif ((istage .ge. 7).and.(istage .le. 8))then
   group = 2.
Elseif ((istage .ge. 9).and.(istage .le. 12))then
   group = 3.
Elseif ((istage .ge. 13).and.(istage .le. 15))then
   group = 4.
Elseif ((istage .ge. 16).and.(istage .le. 18))then
   group = 5.
Elseif (istage .ge. 19)then
   group = 6.
Endif

!Time to hatch
hatch_hrs = 895.97d0*exp(-0.194d0*temp) !hrs to hatch at current temp
egg_hatch = egg_hatch + (1 / hatch_hrs) * dtmin/3600.d0 !accum. time to hatch
fish(ntrac,i_hatchtime) = egg_hatch

If (egg_hatch .ge. 1)then
   stage(ntrac) = f_yolk
   If (yolk_len .lt. 4.5)then
      yolk_len=4.5
   Endif
   fish(ntrac,i_length) = yolk_len
   fish(ntrac,i_age) = 0.0                 !age is now days post hatch
   print*,hatch_hrs,egg_hatch,hatch_len,yolk_len
   stop
Endif

!Length at hatch
If (temp .le. 0.0)then
   hatch_len = 4.0
Else
   hatch_len = 4.8462*temp**0.046d0 !hatch length (mm) at current T
Endif
yolk_len = yolk_len + hatch_len * (1 / hatch_hrs)*dtmin/3600.d0 !accum. hatch length
fish(ntrac,i_hatchlength) = yolk_len

Endif !Egg if

!YOLKSAC
If (stage(ntrac) .eq. f_yolk)then
   length = 4.505 * exp(7.854d0 * (1.d0 - exp(-0.004d0 * age))) ! SL in mm
   fish(ntrac,i_length) = length

   If (length .ge. 5.5)then
      stage(ntrac) = f_pre
   Endif

Endif !Yolk if

!PREFLEXION
!If (stage(ntrac) .eq. f_yolk)then
!   length = 4.505 * exp(7.854 * (1 - exp(-0.004 * age))) ! SL in mm
!   fish(ntrac,i_length) = length

 !  If (length .ge. 10.0)then
 !     stage(ntrac) = f_post
 !  Endif

!Endif !Pre if


!POSTFLEXION
!If (stage(ntrac) .eq. f_yolk)then
 !  length = 4.505 * exp(7.854 * (1 - exp(-0.004 * age))) ! SL in mm
  ! fish(ntrac,i_length) = length

   !If (length .ge. 40)then
    !  stage(ntrac) = f_juv
   !Endif

!Endif !Post if



!-------------BUOYANCY / VERTICAL BEHAVIOR-------------
!EGGS
!If (stage(ntrac) .eq. f_egg)then

grav = 9.81 ! gravity
Dm = fishdiam/1000.0 ! change diameter from mm to m

rho=dens+1025.0  ! If called from loop.F
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

egg_sg = -9.0d-5*group**2 + 0.0007d0*group + 1.0227d0  !specific gravity
rhof = 1.d3 * egg_sg				!density

!Stokes eq determines depth
 wfish = -(1.0d0/18.0d0)*((rhof-rho)/visc)*grav*(Dm**2) !add random variance
 
!if (istage .eq. 21)then
  !print*,hatch_hrs,egg_hatch,hatch_len,yolk_len
  !stop
!endif

!Endif !Egg if

!YOLKSAC
!If (stage(ntrac) .eq. f_yolk)then

!Endif !Yolk if


!PREFLEXION
!If (stage(ntrac) .eq. f_pre)then

!Endif !Pre if

!POSTFLEXION
!If (stage(ntrac) .eq. f_post)then

!Endif !Post if

return
end subroutine fishvel
#endif
