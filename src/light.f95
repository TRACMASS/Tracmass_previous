#ifdef larval_fish
!
! Compute light now from average light over the day at that point.
!
subroutine light_spot
USE mod_param
USE mod_particle
USE mod_time
IMPLICIT none

INTEGER :: month, year
real*8, parameter :: deg2rad = pi / 180.0
real*8  :: Dangle, Hangle, LatRad
real*8  :: cff, cff1, cff2, hour, yday

!  DIURNAL_SRFLUX option: Modulate shortwave radiation SRFLX (which is
!                         read and interpolated elsewhere) by the local
!  diurnal cycle (a function of longitude, latitude and day-of-year).
!  This option is provided for cases where SRFLX computed by SET_DATA is
!  an average over >= 24 hours. For "diurnal_srflux" to work ana_srflux
!  must be undefined. If you want a strictly analytical diurnal cycle
!  enter it explicitly at the end of this subroutine or use the "albedo"
!  option.
!
!  For a review of shortwave radiation formulations check:
!
!    Niemela, S., P. Raisanen, and H. Savijarvi, 2001: Comparison of
!      surface radiative flux parameterizations, Part II, Shortwave
!      radiation, Atmos. Res., 58, 141-154.
!
!-----------------------------------------------------------------------
!
!  Assume time is in modified Julian day.  Get hour and year day.
!
      call updateClock
      hour = currHour + currMin/60.0 + currSec/3600.0
      yday = currJDyr
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44d0*COS((172.0d0-yday)*2.0d0*pi/365.25d0)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
      Hangle=(12.0d0-hour)*pi/12.0d0
!
!  Local daylight is a function of the declination (Dangle) and hour
!  angle adjusted for the local meridian (Hangle-lonr(i,j)/15.0).
!  The 15.0 factor is because the sun moves 15 degrees every hour.
!
      LatRad=flat1*deg2rad
      cff1=SIN(LatRad)*SIN(Dangle)
      cff2=COS(LatRad)*COS(Dangle)
!
!  SRFLX is reset on each time step in subroutine SET_DATA which
!  interpolates values in the forcing file to the current date.
!  This DIURNAL_SRFLUX option is provided so that SRFLX values
!  corresponding to a greater or equal daily average can be modulated
!  by the local length of day to produce a diurnal cycle with the
!  same daily average as the original data.  This approach assumes
!  the net effect of clouds is incorporated into the SRFLX data.
!
!  Normalization = (1/2*pi)*INTEGRAL{ABS(a+b*COS(t)) dt}  from 0 to 2*pi
!                = (a*ARCCOS(-a/b)+SQRT(b**2-a**2))/pi    for |a| < |b|
!
          light=MAX(0.0d0, srflx)
          IF (ABS(cff1) > ABS(cff2)) THEN
           IF (cff1*cff2.gt.0.0d0) THEN
              cff=cff1                                 ! All day case
              light=MAX(0.0d0, light/cff*        &
     &               (cff1+cff2*COS(Hangle-flon1*deg2rad)))
            ELSE
              light=0.0d0                        ! All night case
            END IF
          ELSE
            cff=(cff1*ACOS(-cff1/cff2)+SQRT((cff2+cff1)*(cff2-cff1)))/pi
            IF (cff .lt. 10.e-10) THEN
              light=0.0d0
            ELSE
              light=MAX(0.0d0, light/cff*         &
     &               (cff1+cff2*COS(Hangle-flon1*deg2rad)))
            END IF
          END IF



return
end subroutine light_spot
#endif
