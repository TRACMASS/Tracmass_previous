
  
function sw_dens0 (t,s)
  
  IMPLICIT NONE
  REAL*8, INTENT(IN)                        :: t(:),s(:)
  REAL*8, ALLOCATABLE, DIMENSION (:)        :: sw_dens0
  
  REAL*8, ALLOCATABLE, DIMENSION (:)        :: dens_temp
  
  
  INTEGER                                   :: kmm
  INTEGER                                   :: km
  REAL*8, PARAMETER                         :: a0 = 999.842594
  REAL*8, PARAMETER                         :: a1 =   6.793952e-2
  REAL*8, PARAMETER                         :: a2 =  -9.095290e-3
  REAL*8, PARAMETER                         :: a3 =   1.001685e-4
  REAL*8, PARAMETER                         :: a4 =  -1.120083e-6
  REAL*8, PARAMETER                         :: a5 =   6.536332e-9
    
  REAL*8, PARAMETER                         :: b0 =  8.24493e-1
  REAL*8, PARAMETER                         :: b1 = -4.0899e-3
  REAL*8, PARAMETER                         :: b2 =  7.6438e-5
  REAL*8, PARAMETER                         :: b3 = -8.2467e-7
  REAL*8, PARAMETER                         :: b4 =  5.3875e-9
  
  REAL*8, PARAMETER                         :: c0 = -5.72466e-3
  REAL*8, PARAMETER                         :: c1 = +1.0227e-4
  REAL*8, PARAMETER                         :: c2 = -1.6546e-6
  
  REAL*8, PARAMETER                         :: d0 = 4.8314e-4

  km = size(s)
  allocate ( dens_temp(km) ,sw_dens0(km) )
  
  !  do k=1,kmm
  !     dens_temp(k) = a0+(a1+(a2+(a3+(a4+a5*t(k))*t(k))*t(k))*t(k))*t(k)
  !     dens_zero(k) = dens_temp(k) &
  !          + (b0 + (b1 + (b2 + (b3 + b4*t(k))*t(k))*t(k))*t(k))*s(k) &
  !          + (c0 + (c1 + c2*T(k))*T(k))*s(k)*sqrt(S(k)) + d0*s(k)**2
  !     sw_dens0(k)=dens_zero(k)
  !  enddo
    
  dens_temp = a0+(a1+(a2+(a3+(a4+a5*t)*t)*t)*t)*t
  sw_dens0  = dens_temp &
       + (b0+(b1+(b2+(b3+b4*t)*t)*t)*t)*s &
       + (c0+(c1+c2*T)*T)*s*sqrt(S) + d0*s**2
  

  !  densP0 = sw_dens0(S,T);
  !  K      = sw_seck(S,T,P);
  !  P      = P/10;  ! convert from db to atm pressure units
  !  dens   = densP0./(1-P./K);
end function sw_dens0

