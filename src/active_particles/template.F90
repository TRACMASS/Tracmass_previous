MODULE active_velocities
  REAL                                       :: upr(12,2) = 0

CONTAINS

  subroutine active_niter
    ! Skeleton subroutine instead of active velocities
    USE mod_param
    USE mod_vel, only:  uflux, vflux, wflux, ff
    USE mod_time, only: intrpr, intrpg
    
    IMPLICIT none
    
    INTEGER                                  :: ia, ja, ka, im, jm, n
  
#if defined timestep
    ! time interpolated velocities
    uv(1) = 0 ! Added west  uvel(ia,ja,ka,nsm)
    uv(2) = 0 ! Added east  uvel(im,ja,ka,nsm)
    uv(3) = 0 ! Added north vvel(ia,ja,ka,nsm)
    uv(4) = 0 ! Added south vvel(ia,jm,ka,nsm)
#elif defined timeanalyt
    uv( 1) = 0 ! western u at t-1
    uv( 2) = 0 ! eastern u at t-1
    uv( 3) = 0 ! northern v at t-1
    uv( 4) = 0 ! southern v at t-1
    uv( 7) = 0 ! western u at t
    uv( 8) = 0 ! eastern u at t
    uv( 9) = 0 ! northern v at t
    uv(10) = 0 ! southern v at t
#endif  

    upr(:,1) = upr(:,2) ! store u' from previous time iteration step.


  upr(1,2) = 0 ! Added west uvel(ia,ja,ka,nsm)
  upr(2,2) = 0 ! Added east uvel(im,ja,ka,nsm)
  upr(3,2) = 0 ! Added north vvel(ia,ja,ka,nsm)
  upr(4,2) = 0 ! Added south vvel(ia,jm,ka,nsm)

#ifndef twodim   
  do n=1,2
     
     !Detta ser ut som en bugg   ! t2
     !  upr(5,n) = w(ka-1) - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     !Detta gör att man bara justerar vertikala hastigheten på ovansidan av boxen ! u2
     !  upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
     !  upr(6,n) = 0.d0
     
     ! === The vertical velocity is calculated from 
     ! === the divergence of the horizontal turbulence !v2
     ! === The turbulent is evenly spread between the 
     ! === top and bottom of box if not a bottom box

#ifdef full_wflux
    localW = wflux(ia,ja,ka-1,nsm)
#else
    localW = wflux(ka-1,nsm)
#endif
    
    if(localW.eq.0.d0) then
       upr(5,n) = - ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
       upr(6,n) = 0.d0
    else
       upr(5,n) = - 0.5d0 * ff * ( upr(1,n) - upr(2,n) + upr(3,n) - upr(4,n) )
       upr(6,n) = - upr(5,n)
    endif

#if defined timeanalyt
    if(localW.eq.0.d0) then
       upr(11,n) = - ff * ( upr(7,n) - upr(8,n) + upr(9,n) - upr(10,n) )
       upr(12,n) = 0.d0
    else
       upr(11,n) = - 0.5d0 * ff * ( upr(7,n) - upr(8,n) + upr(9,n) - upr(10,n) )
       upr(12,n) = - upr(11,n)
    endif
#endif
 enddo
#endif
#endif
 
 return
#endif
end subroutine actve_substep
ENDMODULE active_velocities
!_______________________________________________________________________
