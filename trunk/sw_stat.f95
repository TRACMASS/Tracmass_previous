module mod_stat

  interface
     function sw_dens0(t,s)
       REAL, INTENT(IN)                        :: t(:),s(:)       
       REAL, ALLOCATABLE, DIMENSION (:)        :: sw_dens0
     end function sw_dens0
     function sw_pres(depth,lat)
       REAL, INTENT(IN)                        :: depth(:),lat(:)       
       REAL, ALLOCATABLE, DIMENSION (:)        :: sw_pres
     end function sw_pres
     function sw_seck(t,s,pres)
       REAL, INTENT(IN)                        :: t(:),s(:),pres(:)       
       REAL, ALLOCATABLE, DIMENSION (:)        :: sw_seck
     end function sw_seck
  end interface

contains
  subroutine statvd(t, s, dens ,km ,depth ,lat)

    IMPLICIT NONE
    
    INTEGER                                     :: km
    REAL, ALLOCATABLE, DIMENSION(:)             :: t ,s ,dens0 ,pres 
    REAL, ALLOCATABLE, DIMENSION(:)             :: depth ,lat ,seck ,dens
    
    allocate ( dens0(km) ,pres(km) ,seck(km))
    
    dens0  = sw_dens0(t,s)
    pres   = sw_pres(depth,lat)
    seck   = sw_seck(t,s,pres)
    
    dens   = dens0/(1-(pres/10)/seck)
  end subroutine statvd


  
!   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  subroutine statv(t, s, rho, km)
    
    IMPLICIT NONE

!    REAL*4, ALLOCATABLE, DIMENSION (:)       :: dens_temp, dens_zero, rho
!    REAL*4, ALLOCATABLE, DIMENSION (:)       :: t, s
    INTEGER                                  :: KM,k
    
    REAL, PARAMETER                          :: a0 = 999.842594
    REAL, PARAMETER                          :: a1 =   6.793952e-2
    REAL, PARAMETER                          :: a2 =  -9.095290e-3
    REAL, PARAMETER                          :: a3 =   1.001685e-4
    REAL, PARAMETER                          :: a4 =  -1.120083e-6
    REAL, PARAMETER                          :: a5 =   6.536332e-9
    
    REAL, PARAMETER                          :: b0 =  8.24493e-1
    REAL, PARAMETER                          :: b1 = -4.0899e-3
    REAL, PARAMETER                          :: b2 =  7.6438e-5
    REAL, PARAMETER                          :: b3 = -8.2467e-7
    REAL, PARAMETER                          :: b4 =  5.3875e-9
    
    REAL, PARAMETER                          :: c0 = -5.72466e-3
    REAL, PARAMETER                          :: c1 = +1.0227e-4
    REAL, PARAMETER                          :: c2 = -1.6546e-6
    
    REAL, PARAMETER                          :: d0 = 4.8314e-4

!    allocate ( dens_temp(KM),dens_zero(KM) ,rho(KM) )
    REAL*4 t(KM),s(KM)  
    REAL*4 dens_temp(KM),dens_zero(KM) ,rho(KM) 

    do k=1,km
       dens_temp(k) = a0+(a1+(a2+(a3+(a4+a5*t(k))*t(k))*t(k))*t(k))*t(k)
       dens_zero(k) = dens_temp(k) &
            + (b0 + (b1 + (b2 + (b3 + b4*t(k))*t(k))*t(k))*t(k))*s(k) &
            + (c0 + (c1 + c2*T(k))*T(k))*s(k)*sqrt(S(k)) + d0*s(k)**2
       rho(k)=dens_zero(k)-1000.
!print *,'k=',k,rho(k),t(k),s(k)
     enddo

!print *,'statv11111111=',km,rho(km)

!  densP0 = sw_dens0(S,T);
!  K      = sw_seck(S,T,P);
!  P      = P/10;  ! convert from db to atm pressure units
!  dens   = densP0./(1-P./K);

end subroutine statv
end module mod_stat
