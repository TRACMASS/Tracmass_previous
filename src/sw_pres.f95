
function sw_pres(DEPTH,LAT)
  
  ! SW_PRES    Pressure from depth
  !==================================================================
  ! SW_PRES   $Revision: 1.5 $  $Date: 1994/10/11 01:23:32 $
  !           Copyright (C) CSIRO, Phil Morgan 1993.
  !
  ! USAGE:  pres = sw_pres(depth,lat)
  !
  ! DESCRIPTION:
  !    Calculates pressure in dbars from depth in meters.
  !
  ! INPUT:  (all must have same dimensions)
  !   depth = depth [metres]  
  !   lat   = Latitude in decimal degress north [-90..+90]
  !           (LAT may have dimensions 1x1 or 1xn where depth(mxn) )
  !
  ! OUTPUT:
  !  pres   = Pressure    [db]
  !
  ! AUTHOR:  Phil Morgan 93-06-25  (morgan@ml.csiro.au)
  !
  ! DISCLAIMER:
  !   This software is provided "as is" without warranty of any kind.  
  !   See the file sw_copy.m for conditions of use and licence.
  !
  ! REFERENCES:
  !    Saunders, P.M. 1981
  !    "Practical conversion of Pressure to Depth"
  !    Journal of Physical Oceanography, 11, 573-574
  !
  ! CHECK VALUE: 
  !    P=7500.00 db for LAT=30 deg, depth=7321.45 meters
  !==================================================================

  IMPLICIT NONE
  REAL, INTENT(IN)                           :: depth(:),lat(:)       
  REAL, ALLOCATABLE, DIMENSION (:)           :: sw_pres ,x ,C1
  REAL                                       :: deg2rad
  REAL, PARAMETER                            :: pi=3.14159265358979
  INTEGER                                    :: km
  km = size(depth)
  allocate ( sw_pres(km) ,x(km) ,C1(km) )


  deg2rad = pi/180;
  x       = sin(abs(LAT)*deg2rad);  ! convert to radians
  C1      = 5.92E-3+x**2*5.25E-3;
  sw_pres = ((1-C1)-sqrt(((1-C1)**2)-(8.84E-6*DEPTH)))/4.42E-6;
  
end function sw_pres
