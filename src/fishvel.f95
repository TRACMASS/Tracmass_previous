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
! CMP 11/13/12
! Modified to add egg density as a function of age group
! Age group is based on stage
! Stage is a function of age and temperature
! CMP 1/29/13
! Size at hatch is a random number from normal distr with mean and std from data
! Egg, yolksac, and preflexion will do random walk within mixed layer

subroutine fishvel        ! If called from loop.F
USE mod_param
USE mod_fish
USE mod_loopvars
USE mod_particle
USE mod_time
USE mod_traj
IMPLICIT none

 REAL*8   :: hatchJD
 REAL*8   :: age, egg_hatch, yolk_len, length
 REAL*8   :: hatch_hrs, hatch_len
 REAL*8   :: r1, r2, r, theta
 INTEGER  :: clock, seed
 !DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462

 length = fish(ntrac,i_length)
 CALL SYSTEM_CLOCK(COUNT=clock)
 seed = clock
 CALL RANDOM_SEED(seed)
 !call init_random_seed()

!--------------DEVELOPMENT / GROWTH-----------------------
!EGGS
 If (stage(ntrac) .eq. f_egg)then
   !Time to hatch
   egg_hatch = fish(ntrac,i_hatchtime)
   hatch_hrs = 895.97*exp(-0.194*temp) !hrs to hatch at current temp
   egg_hatch = egg_hatch + (1.0 / hatch_hrs) * dtmin/3600. !accum time to hatch in hrs
   fish(ntrac,i_hatchtime) = egg_hatch

   If (egg_hatch .ge. 1)then
      stage(ntrac) = f_yolk
      !Length at hatch
      !draw yolk length from a normal distribution with mean=5.125 and stdev=0.46
      CALL RANDOM_NUMBER(r1)   
      CALL RANDOM_NUMBER(r2) 
      r =(-2.0d0*log(r1))**0.5
      theta = 2.0d0*PI*r2
      yolk_len = 5.125 + (0.46**2)*r*sin(theta)
      fish(ntrac,i_hatchlength) = yolk_len
      fish(ntrac,i_jd) = currJDtot         !Use JD to calc dph
   Endif
 Endif !Egg if

!YOLKSAC
 If (stage(ntrac) .eq. f_yolk)then
   hatchJD = fish(ntrac,i_jd)
   age = currJDtot - hatchJD
   yolk_len = fish(ntrac,i_hatchlength)
   length = yolk_len * exp(7.854d0 * (1.d0 - exp(-0.004d0 * age))) ! SL in mm
   fish(ntrac,i_length) = length
   If (length .ge. 6.0)then
      stage(ntrac) = f_pre
   Endif
 Endif !Yolk if

!PREFLEXION
 If (stage(ntrac) .eq. f_pre)then
   hatchJD = fish(ntrac,i_jd)
   age = currJDtot - hatchJD
   yolk_len = fish(ntrac,i_hatchlength)
   length = yolk_len * exp(7.854d0 * (1.d0 - exp(-0.004d0 * age))) ! SL in mm
   fish(ntrac,i_length) = length
   If (length .ge. 10.0)then
      stage(ntrac) = f_post
   Endif   
 Endif !Pre if


!POSTFLEXION
 If (stage(ntrac) .eq. f_post)then
   hatchJD = fish(ntrac,i_jd)
   age = currJDtot - hatchJD
   yolk_len = fish(ntrac,i_hatchlength)
   length = yolk_len * exp(7.854d0 * (1.d0 - exp(-0.004d0 * age))) ! SL in mm
   fish(ntrac,i_length) = length
   If (length .ge. 40.0)then
      stage(ntrac) = f_juv
   Endif
 Endif !Post if



!----------------VERTICAL BEHAVIOR----------------
!EGGS
!If (stage(ntrac) .eq. f_egg)then
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

