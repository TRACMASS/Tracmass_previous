function sw_seck(t,s,pres)
  
  ! SW_SECK    Secant bulk modulus (K) of sea water
  !=========================================================================
  ! SW_SECK  $Revision: 1.3 $  $Date: 1994/10/10 05:50:45 $
  !          Copyright (C) CSIRO, Phil Morgan 1992.
  !
  ! USAGE:  dens = sw_seck(S,T,P)
  !
  ! DESCRIPTION:
  !    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980. 
  !    UNESCO polynomial implementation.
  !
  ! INPUT:  (all must have same dimensions)
  !   S = salinity    [psu      (PSS-78) ]
  !   T = temperature [degree C (IPTS-68)]
  !   P = pressure    [db]
  !       (alternatively, may have dimensions 1*1 or 1*n 
  !        where n is columns in S)
  !
  ! OUTPUT:
  !   K = Secant Bulk Modulus  [bars]
  ! 
  ! AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
  !
  ! DISCLAIMER:
  !   This software is provided "as is" without warranty of any kind.  
  !   See the file sw_copy.m for conditions of use and licence.
  !
  ! REFERENCES:
  !    Fofonoff, P. and Millard, R.C. Jr
  !    Unesco 1983. Algorithms for computation of fundamental properties of 
  !    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !    Eqn.(15) p.18
  !
  !    Millero, F.J. and  Poisson, A.
  !    International one-atmosphere equation of state of seawater.
  !    Deep-Sea Res. 1981. Vol28A(6) pp625-629.
  !=========================================================================
  ! Pure water terms of the secant bulk modulus at atmos pressure.
  ! UNESCO eqn 19 p 18
  IMPLICIT NONE

  REAL*8, INTENT(IN)                       :: t(:),s(:),pres(:)
  REAL*8, ALLOCATABLE, DIMENSION (:)       :: sw_seck
  REAL*8, ALLOCATABLE, DIMENSION (:)       :: p,AW,A,BW,B,KW,k,SR
  INTEGER                                  :: km

  REAL*8, PARAMETER                        :: h3 = -5.77905E-7
  REAL*8, PARAMETER                        :: h2 = +1.16092E-4
  REAL*8, PARAMETER                        :: h1 = +1.43713E-3
  REAL*8, PARAMETER                        :: h0 = +3.239908![-0.1194975]

  REAL*8, PARAMETER                        :: k2 =  5.2787E-8
  REAL*8, PARAMETER                        :: k1 = -6.12293E-6
  REAL*8, PARAMETER                        :: k0 = +8.50935E-5![+3.47718E-5]    
  REAL*8, PARAMETER                        :: e4 = -5.155288E-5
  REAL*8, PARAMETER                        :: e3 = +1.360477E-2
  REAL*8, PARAMETER                        :: e2 = -2.327105
  REAL*8, PARAMETER                        :: e1 = +148.4206
  REAL*8, PARAMETER                        :: e0 = 19652.21![-1930.06]  

  REAL*8, PARAMETER                        :: j0 = 1.91075E-4

  REAL*8, PARAMETER                        :: i2 = -1.6078E-6
  REAL*8, PARAMETER                        :: i1 = -1.0981E-5
  REAL*8, PARAMETER                        :: i0 =  2.2838E-3

  REAL*8, PARAMETER                        :: m2 =  9.1697E-10
  REAL*8, PARAMETER                        :: m1 = +2.0816E-8
  REAL*8, PARAMETER                        :: m0 = -9.9348E-7

  REAL*8, PARAMETER                        :: f3 =  -6.1670E-5
  REAL*8, PARAMETER                        :: f2 =  +1.09987E-2
  REAL*8, PARAMETER                        :: f1 =  -0.603459
  REAL*8, PARAMETER                        :: f0 = +54.6746

  REAL*8, PARAMETER                        :: g2 = -5.3009E-4
  REAL*8, PARAMETER                        :: g1 = +1.6483E-2
  REAL*8, PARAMETER                        :: g0 = +7.944E-2

  km = size(t)
  allocate ( sw_seck(km) ,AW(km) ,BW(km) ,KW(km) ,SR(km) )
  allocate ( p(km) ,A(km) ,B(km) ,K(km) )

  p       = pres/10;    !convert from db to atmospheric pressure units

  AW      = h0 + (h1 + (h2 + h3*T)*T)*T
  BW      = k0 + (k1 + k2*T) * T
  KW      = e0 + (e1 + (e2 + (e3 + e4*T)*T)*T)*T     ! eqn 19
  SR      = sqrt(s)
  A       = AW + (i0 + (i1 + i2*T)*T + j0 * SR) * S
  B       = BW + (m0 + (m1 + m2*T)*T)*S               ! eqn 18
  K       = KW + ( f0 + (f1 + (f2 + f3*T)*T)*T &
               + ( g0 + (g1 + g2*T)*T)*SR         )*S  ! eqn 16    
  
  sw_seck = K  + (A + B*p)*p  ! eqn 15
  
end function sw_seck

