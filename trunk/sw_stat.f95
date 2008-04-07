
module mod_stat
contains
  subroutine statv(t, s, d, rho, kmm)
    
    REAL*4, ALLOCATABLE, DIMENSION (:) :: dens_temp, dens_zero, t, s, rho, p
    INTEGER KMM,k
    
    REAL, PARAMETER ::  a0 = 999.842594
    REAL, PARAMETER ::  a1 =   6.793952e-2
    REAL, PARAMETER ::  a2 =  -9.095290e-3
    REAL, PARAMETER ::  a3 =   1.001685e-4
    REAL, PARAMETER ::  a4 =  -1.120083e-6
    REAL, PARAMETER ::  a5 =   6.536332e-9
    
    REAL, PARAMETER :: b0 =  8.24493e-1
    REAL, PARAMETER :: b1 = -4.0899e-3
    REAL, PARAMETER :: b2 =  7.6438e-5
    REAL, PARAMETER :: b3 = -8.2467e-7
    REAL, PARAMETER :: b4 =  5.3875e-9
    
    REAL, PARAMETER :: c0 = -5.72466e-3
    REAL, PARAMETER :: c1 = +1.0227e-4
    REAL, PARAMETER :: c2 = -1.6546e-6
    
    REAL, PARAMETER :: d0 = 4.8314e-4
    
    REAL, PARAMETER :: h3 = -5.77905E-7
    REAL, PARAMETER :: h2 = +1.16092E-4
    REAL, PARAMETER :: h1 = +1.43713E-3
    REAL, PARAMETER :: h0 = +3.239908
    
    REAL, PARAMETER :: k2 =  5.2787E-8
    REAL, PARAMETER :: k1 = -6.12293E-6
    REAL, PARAMETER :: k0 =  +8.50935E-5
 
    REAL, PARAMETER :: e4 = -5.155288E-5
    REAL, PARAMETER :: e3 = +1.360477E-2
    REAL, PARAMETER :: e2 = -2.327105
    REAL, PARAMETER :: e1 = +148.4206
    REAL, PARAMETER :: e0 = 19652.21
   
    REAL, PARAMETER :: j0 = 1.91075E-4;

    REAL, PARAMETER :: i2 = -1.6078E-6;
    REAL, PARAMETER :: i1 = -1.0981E-5;
    REAL, PARAMETER :: i0 =  2.2838E-3;

    REAL, PARAMETER :: m2 =  9.1697E-10
    REAL, PARAMETER :: m1 = +2.0816E-8
    REAL, PARAMETER :: m0 = -9.9348E-7
    
    REAL, PARAMETER :: f3 =  -6.1670E-5
    REAL, PARAMETER :: f2 =  +1.09987E-2
    REAL, PARAMETER :: f1 =  -0.603459
    REAL, PARAMETER :: f0 = +54.6746
    
    REAL, PARAMETER :: g2 = -5.3009E-4
    REAL, PARAMETER :: g1 = +1.6483E-2
    REAL, PARAMETER :: g0 = +7.944E-2
    REAL :: C1, AW,BW,KW, SR, A, B, K0, K

    allocate ( p(KM), dens_temp(KM),dens_zero(KM) )

    C1=0
    P=((1-C1)-sqrt(((1-C1).^2)-(8.84E-6*z)))/4.42E-6)/10
    
    do k=1,kmm
       dens_temp(k)=a0+(a1+(a2+(a3+(a4+a5*t(k))*t(k))*t(k))*t(k))*t(k)
       dens_zero(k)=dens_temp(k)+ & 
            (b0+(b1+(b2+(b3+b4*t(k))*t(k))*t(k))*t(k))*s(k) + &
            (c0+(c1+c2*T(k))*T(k))*s(k)*sqrt(S(k))+d0*s(k)**2
       
       AW  = h0+(h1+(h2+h3*T(k))*T(k))*T(k)
       BW  = k0+(k1+k2*T(k))*T(k)
       KW  = e0+(e1+(e2+(e3+e4*T(k))*T(k))*T(k))*T(k)
       SR = sqrt(S(k))
       A=AW+(i0+(i1+i2*T(k))*T(k)+j0*SR)*S(k)
       B=BW+(m0+(m1+m2*T(k))*T(k))*S(k)
       K0=KW+(f0+(f1+(f2+f3*T(k))*T(k))*T(k)+ &
            (g0+(g1+g2*T(k))*T(k))*SR)*S(k)
       K = K0+(A+B*P(k))*P(k)
       
       rho(k)=dens_zero(k)/(1-P./K)-1000.
    enddo
    
    !  densP0 = sw_dens0(S,T);
    !  K      = sw_seck(S,T,P);
    !  P      = P/10;  ! convert from db to atm pressure units
    !  dens   = densP0./(1-P./K);
    
  end subroutine statv
end module mod_stat
