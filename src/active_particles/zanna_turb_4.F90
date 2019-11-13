MODULE mod_active_particles

  USE mod_time
  USE mod_traj
  USE mod_grid
  USE mod_param
  USE mod_deformation
  USE mod_vel
  USE mod_loopvars

  IMPLICIT none

  REAL                                       :: upr(12,2) = 0.0
  REAL                                       :: kappa = -0.2
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)      :: mem_param
  REAL(DP), ALLOCATABLE, DIMENSION(:)        :: lapu_b, lapv_b, lapu_n, lapv_n, dlapu_n, dlapv_n, dt_n
  REAL(DP), DIMENSION(4)                     :: test_sum


CONTAINS

  subroutine active_init
    ALLOCATE ( mem_param(ntracmax,12) )
    ALLOCATE ( lapu_b(ntracmax), lapu_n(ntracmax), lapv_b(ntracmax), lapv_n(ntracmax), &
             & dlapu_n(ntracmax),  dlapv_n(ntracmax), dt_n(ntracmax) )
    mem_param(:,:) = 0.
    lapu_b(:) = 0.
    lapu_n(:) = 0.
    lapv_b(:) = 0.
    lapv_n(:) = 0.
    dlapu_n(:)  = 0.
    dlapv_n(:)  = 0.
    dt_n(:) = 0.
    test_sum(:) = 0.
    return
  end subroutine active_init

  subroutine active_ints(ints)
    integer, intent(in) :: ints    
    return
  end subroutine active_ints

  subroutine active_ntrac(ntrac)
     
     integer, intent(in) :: ntrac
     integer :: ibp
     real(DP) :: ddx,dlapu,dlapv,lapu1,lapu2,lapv1,lapv2,frac !! temporary variables
     real(DP) :: xx,yy,zz,zuu,zum,zvv,zvm
     integer  :: ic,im,ip,jm,jc,jp
          
    return
  end subroutine active_ntrac
  
  subroutine active_niter
     !! Calculate new turbulent fluxes
     !! for each step a particle takes
     !!
     !! RE tensor parameterisation says, after some assumptions, 
     !! du/dt = kappa * dx**2 * d/dt(laplacian(u))
     !!
     !! We then write
     !! u2 = u1 + kappa * dx**2 * dt * d/dt(lapu)
     !! where lapu = d2u/dx2 + d2u/dy2
     !!
     !! Now, the dt * d/dt do not necessarily cancel. 
     !! dt is the time it takes for the particle to 
     !! go travel in this niter iteration. 
     !! d/dt is the inverse of the time between lapu calculations
     !! which maybe between GCM steps. 
     !!
     !! Also, dt is not known yet, since we have not called cross yet.
     !!
     integer  :: ic,im,ip,jm,jc,jp,jm2,iam2 !! temporary indices
     real(dp) :: zup, zvp, zx, zy, zluim, zlui, zdlu, zlvjm, zlvj, zdlv, ddx, zuu, zum, zvv, zvm, zdt
     real(dp) :: zludu, zlvdv, zludu1, zlvdv1, zludu2, zlvdv2, &
               & zux_11, zux_12, zux_21, zux_22, &
               & zvx_11, zvx_12, zvx_21, zvx_22, &
               & zuy_11, zuy_12, zuy_21, zuy_22, &
               & zvy_11, zvy_12, zvy_21, zvy_22, zvx, zvy, zux, zuy !! temporary variables     
     real(dp) :: zfrac1, zfrac2, zfrac3, zfrac4
          
     !! calculate Laplacian of u and v 
     !! Interpolate in time and x for u, and y for v
     !! so that even if particle says in the same grid cell, d lap(u)/dt is not zero
     zx = x0 - dble(iam) ! between 0 and 1
     zy = y0 - dble(ja-1)
     
     lapu_b(ntrac) = lapu_n(ntrac)
     lapv_b(ntrac) = lapv_n(ntrac)
     lapu_n(ntrac) = (intrpg * lapu(iam,ja  ,ka,nsp) + (1.d0 - intrpg) * lapu(iam,ja  ,ka,nsm)) * (1.d0 - zx) + &
                     (intrpg * lapu(ia ,ja  ,ka,nsp) + (1.d0 - intrpg) * lapu(ia ,ja  ,ka,nsm)) * zx
     lapv_n(ntrac) = (intrpg * lapv(ia ,ja-1,ka,nsp) + (1.d0 - intrpg) * lapv(ia ,ja-1,ka,nsm)) * (1.d0 - zy) + &
                     (intrpg * lapv(ia ,ja  ,ka,nsp) + (1.d0 - intrpg) * lapv(ia ,ja  ,ka,nsm)) * zy
     
     dlapu_n(ntrac) = (lapu_n(ntrac) - lapu_b(ntrac)) !/ dt_n(ntrac)                                                                               
     dlapv_n(ntrac) = (lapv_n(ntrac) - lapv_b(ntrac)) !/ dt_n(ntrac) 
     
     !print*,intrpg,lapu(iam,ja  ,ka,nsp),(1.d0 - intrpg),lapu(iam,ja  ,ka,nsm),(1.d0 - zx),&
     !       intrpg,lapu(ia ,ja  ,ka,nsp),(1.d0 - intrpg),lapu(ia ,ja  ,ka,nsm),zx
     !print*,'lapu_b, lapu_n ',lapu_b(ntrac),lapu_n(ntrac)
     !! find grid box size                                                                                                             
     ddx    = 0.5 * (dxv(ia,ja)+dxv(ia,ja-1)) + 0.5 * (dyu(ia,ja)+dyu(iam,ja))
          
#if defined add_lapu_adv
     ip = ia+1
     iam2 = iam-1
     jp = ja+1
     jm = ja-1
     jm2 = jm-1
     if (ip > imt) ip=ip-imt
     if (iam2 < 1) iam2=iam
     if (jp > jmt) jp=ja
     if (jm < 1) jm=ja
     if (jm2 <1) jm2=jm
     !! Find du/dx
     zux_11 = (uflux(ia,ja,ka,nsm) - uflux(iam2,ja,ka,nsm)) / (dxv(ia,ja)*2.d0) ! du/dx at (iam,nsm)
     zux_12 = (uflux(ia,ja,ka,nsp) - uflux(iam2,ja,ka,nsp)) / (dxv(ia,ja)*2.d0) !          (iam,nsp)
     zux_21 = (uflux(ip,ja,ka,nsm) - uflux(iam ,ja,ka,nsm)) / (dxv(ia,ja)*2.d0) !          (ia,nsm)
     zux_22 = (uflux(ip,ja,ka,nsp) - uflux(iam ,ja,ka,nsp)) / (dxv(ia,ja)*2.d0) !          (ia,nsp)
     
     !! Find du/dy 
     zuy_11 = (uflux(ia,jp,ka,nsm) - uflux(ia,jm,ka,nsm)) / (dyu(ia,ja)*2.d0)   ! du/dy at (ia,nsm)
     zuy_12 = (uflux(ia,jp,ka,nsp) - uflux(ia,jm,ka,nsm)) / (dyu(ia,ja)*2.d0)   !          (ia,nsp)
     zuy_21 = (uflux(iam,jp,ka,nsm) - uflux(iam,jm,ka,nsm)) / (dyu(ia,ja)*2.d0) !          (iam,nsm)
     zuy_22 = (uflux(iam,jp,ka,nsp) - uflux(iam,jm,ka,nsp)) / (dyu(ia,ja)*2.d0) !          (iam,nsp)
     
     !! Find dv/dx 
     zvx_21 = (vflux(ip,ja,ka,nsm) - vflux(iam,ja,ka,nsm)) / (dxv(ia,ja)*2.d0)  ! dv/dx at (ja,nsm)
     zvx_22 = (vflux(ip,ja,ka,nsp) - vflux(iam,ja,ka,nsp)) / (dxv(ia,ja)*2.d0)  !          (ja,nsp)
     zvx_11 = (vflux(ip,jm,ka,nsm) - vflux(iam,jm,ka,nsm)) / (dxv(ia,jm)*2.d0)  !          (jm,nsm)
     zvx_12 = (vflux(ip,jm,ka,nsp) - vflux(iam,jm,ka,nsp)) / (dxv(ia,jm)*2.d0)  !          (jm,nsp)
     
     !! Find dv/dy
     zvy_21 = (vflux(ia,jp,ka,nsm) - vflux(ia,jm,ka,nsm)) / (dyu(ia,ja)*2.d0) ! dv/dy at (ja,nsm)
     zvy_22 = (vflux(ia,jp,ka,nsp) - vflux(ia,jm,ka,nsp)) / (dyu(ia,ja)*2.d0) !       at (ja,nsp)
     zvy_11 = (vflux(ia,ja,ka,nsm) - vflux(ia,jm2,ka,nsm)) / (dyu(ia,ja)*2.d0) !      at (jm,nsm)
     zvy_12 = (vflux(ia,ja,ka,nsp) - vflux(ia,jm2,ka,nsp)) / (dyu(ia,ja)*2.d0) !      at (jm,nsp)
     
     zux    = (intrpg * zux_12 + (1.d0 - intrpg) * zux_11) * (1.d0 - zx) + &
            & (intrpg * zux_22 + (1.d0 - intrpg) * zux_21) * zx
     zuy    = (intrpg * zuy_12 + (1.d0 - intrpg) * zuy_11) * (1.d0 - zx) + &
            & (intrpg * zuy_22 + (1.d0 - intrpg) * zuy_21) * zx
     zvx    = (intrpg * zvx_12 + (1.d0 - intrpg) * zvx_11) * (1.d0 - zy) + &
            & (intrpg * zvx_22 + (1.d0 - intrpg) * zvx_21) * zy
     zvy    = (intrpg * zvy_12 + (1.d0 - intrpg) * zvy_11) * (1.d0 - zy) + &
            & (intrpg * zvy_22 + (1.d0 - intrpg) * zvy_21) * zy
     
     zux = zux / (dyu(ia,ja) * dzt(ia,ja,ka,nsm))
     zuy = zuy / (dyu(ia,ja) * dzt(ia,ja,ka,nsp))
     zvx = zvx / (dxv(ia,ja) * dzt(ia,ja,ka,nsm))
     zvy = zvy / (dxv(ia,ja) * dzt(ia,ja,ka,nsp))
     
     !! Add lap(u) * nabla(u)                                                                                                                    
     !! i.e. (d2/dx2 + d2/dy2) u * du/dx + (d2/dx2 + d2/dy2) v * du/dy                                                                           
     !! and  (d2/dx2 + d2/dy2) u * dv/dx + (d2/dx2 + d2/dy2) v * dv/dy                                                                           
          
     zludu = -(lapu_n(ntrac) * zux + lapv_n(ntrac) * zuy) * dt_n(ntrac) 
     zlvdv = -(lapu_n(ntrac) * zvx + lapv_n(ntrac) * zvy) * dt_n(ntrac) 

#elif defined add_lapu_adv2
     ip = ia+1
     iam2 = iam-1
     jp = ja+1
     jm = ja-1
     jm2 = jm-1
     if (ip > imt) ip=ip-imt
     if (iam2 < 1) iam2=iam
     if (jp > jmt) jp=ja
     if (jm < 1) jm=ja
     if (jm2 <1) jm2=jm
     
     !! Find du/dx 
     zux_11 = (uflux(ia,ja,ka,nsm) - uflux(iam,ja,ka,nsm)) / dxv(ia,ja) ! du/dx at (ia,nsm)  
     zux_12 = (uflux(ia,ja,ka,nsp) - uflux(iam,ja,ka,nsp)) / dxv(ia,ja) !          (ia,nsp) 
     
     zuy_11 = (uflux(ia,jp,ka,nsm) - uflux(ia,jm,ka,nsm)) / (dyu(ia,ja)*2.d0)   ! du/dy at (ia,nsm)
     zuy_12 = (uflux(ia,jp,ka,nsp) - uflux(ia,jm,ka,nsm)) / (dyu(ia,ja)*2.d0)   !          (ia,nsp) 
     zuy_21 = (uflux(iam,jp,ka,nsm) - uflux(iam,jm,ka,nsm)) / (dyu(ia,ja)*2.d0) !          (iam,nsm) 
     zuy_22 = (uflux(iam,jp,ka,nsp) - uflux(iam,jm,ka,nsp)) / (dyu(ia,ja)*2.d0) !          (iam,nsp)                                                                               
     !! Find dv/dx 
     zvx_21 = (vflux(ip,ja,ka,nsm) - vflux(iam,ja,ka,nsm)) / (dxv(ia,ja)*2.d0)  ! dv/dx at (ja,nsm)  
     zvx_22 = (vflux(ip,ja,ka,nsp) - vflux(iam,ja,ka,nsp)) / (dxv(ia,ja)*2.d0)  !          (ja,nsp) 
     zvx_11 = (vflux(ip,jm,ka,nsm) - vflux(iam,jm,ka,nsm)) / (dxv(ia,jm)*2.d0)  !          (jm,nsm) 
     zvx_12 = (vflux(ip,jm,ka,nsp) - vflux(iam,jm,ka,nsp)) / (dxv(ia,jm)*2.d0)  !          (jm,nsp)                                                                              
     !! Find dv/dy 
     zvy_11 = (vflux(ia,ja,ka,nsm) - vflux(ia,jm,ka,nsm)) / dyu(ia,ja) ! dv/dy at (ja,nsm) 
     zvy_12 = (vflux(ia,ja,ka,nsp) - vflux(ia,jm,ka,nsp)) / dyu(ia,ja) !       at (ja,nsp)
     
     zux    = (intrpg * zux_12 + (1.d0 - intrpg) * zux_11) 
     zuy    = (intrpg * zuy_12 + (1.d0 - intrpg) * zuy_11) * (1.d0 - zx) + &
            & (intrpg * zuy_22 + (1.d0 - intrpg) * zuy_21) * zx
     zvx    = (intrpg * zvx_12 + (1.d0 - intrpg) * zvx_11) * (1.d0 - zy) + &
            & (intrpg * zvx_22 + (1.d0 - intrpg) * zvx_21) * zy
     zvy    = (intrpg * zvy_12 + (1.d0 - intrpg) * zvy_11) 
     
     zux = zux / (dyu(ia,ja) * dzt(ia,ja,ka,nsm))
     zuy = zuy / (dyu(ia,ja) * dzt(ia,ja,ka,nsp))
     zvx = zvx / (dxv(ia,ja) * dzt(ia,ja,ka,nsm))
     zvy = zvy / (dxv(ia,ja) * dzt(ia,ja,ka,nsp))
     
     !! Add lap(u) * nabla(u)                                                               
     !! i.e. (d2/dx2 + d2/dy2) u * du/dx + (d2/dx2 + d2/dy2) v * du/dy                      
     !! and  (d2/dx2 + d2/dy2) u * dv/dx + (d2/dx2 + d2/dy2) v * dv/dy 
     zludu = -(lapu_n(ntrac) * zux + lapv_n(ntrac) * zuy) * dt_n(ntrac)
     zlvdv = -(lapu_n(ntrac) * zvx + lapv_n(ntrac) * zvy) * dt_n(ntrac)
     
#else
     zludu = 0.
     zlvdv = 0.
#endif
     
     !! delta t is length of time step between new velocity fields
     !! however, it should be the time step of the particle, ds 
     !! but we havent called cross and calculated ds yet
     !zdt = dtmin/tseas !=1/iter
     
     !zdt = dt_n(ntrac)     
     zdt = 1.d0
     
     !! store turbulent velocity increments in upr
     upr(1,:) = kappa * ddx**2 * (dlapu_n(ntrac) + zludu) * zdt  !! param at ia     
     upr(2,:) = kappa * ddx**2 * (dlapu_n(ntrac) + zludu) * zdt  !! param at iam      
     upr(3,:) = kappa * ddx**2 * (dlapv_n(ntrac) + zlvdv) * zdt  !! param at ja     
     upr(4,:) = kappa * ddx**2 * (dlapv_n(ntrac) + zlvdv) * zdt  !! param at ja-1        
     upr(5,:) = 0.
     upr(6,:) = 0.
     upr(7,:) = kappa * ddx**2 * (dlapu_n(ntrac) + zludu) * zdt
     upr(8,:) = kappa * ddx**2 * (dlapu_n(ntrac) + zludu) * zdt
     upr(9,:) = kappa * ddx**2 * (dlapv_n(ntrac) + zlvdv) * zdt
     upr(10,:) = kappa * ddx**2 *(dlapv_n(ntrac) + zlvdv) * zdt
     upr(11,:) = 0.
     upr(12,:) = 0.     
     
     !! integral of all turbulent velocities                                                                           
     mem_param(ntrac,1) = mem_param(ntrac,1) + upr(1,1)
     mem_param(ntrac,2) = mem_param(ntrac,2) + upr(2,1)
     mem_param(ntrac,3) = mem_param(ntrac,3) + upr(3,1)
     mem_param(ntrac,4) = mem_param(ntrac,4) + upr(4,1)

     upr(1,:) = mem_param(ntrac,1)
     upr(2,:) = mem_param(ntrac,2)
     upr(3,:) = mem_param(ntrac,3)
     upr(4,:) = mem_param(ntrac,4)
     upr(7,:) = mem_param(ntrac,1)
     upr(8,:) = mem_param(ntrac,2)
     upr(9,:) = mem_param(ntrac,3)
     upr(10,:) = mem_param(ntrac,4)
     
     !print*,'upr ',upr(1,1),upr(2,1),upr(3,1),upr(4,1)
     
     !! Calculate fluxes on grid cell walls, 
     !! we may need it for later... 
     zuu = (intrpg*uflux(ia ,ja  ,ka,nsp)+intrpr*uflux(ia ,ja  ,ka,nsm))*ff
     zum = (intrpg*uflux(iam,ja  ,ka,nsp)+intrpr*uflux(iam,ja  ,ka,nsm))*ff
     zvv = (intrpg*vflux(ia ,ja  ,ka,nsp)+intrpr*vflux(ia ,ja  ,ka,nsm))*ff
     zvm = (intrpg*vflux(ia ,ja-1,ka,nsp)+intrpr*vflux(ia ,ja-1,ka,nsm))*ff

     !! Adjust upr so that we dont change the sign of 
     !! uflux and vflux on grid box walls
     if (upr(1,1) /= 0.d0) then 
        !! Set limiter in case to ensure uu+upr has same sign as uu 
        if (abs(upr(1,1)) > abs(zuu)) then
           zfrac1 = 0.99 * zuu/upr(1,1)
        else
           zfrac1 = 1.d0
        end if 
     else 
        zfrac1 = 0.d0
     end if
     
     if (upr(2,1) /= 0.d0) then
        if (abs(upr(2,1)) > abs(zum)) then
           zfrac2 = 0.99 * zum/upr(2,1)
        else
           zfrac2 = 1.d0
        end if
     else
        zfrac2 = 0.d0
     end if
     
     if (upr(3,1) /= 0.d0) then
        if (abs(upr(3,1)) > abs(zvv)) then
           zfrac3 = 0.99 * zvv/upr(3,1)
        else
           zfrac3 = 1.d0
        end if
     else
        zfrac3 = 0.d0
     end if
     
     if (upr(4,1) /= 0.d0) then
        if (abs(upr(4,1)) > abs(zvm)) then
           zfrac4 = 0.99 * zvm/upr(4,1)
        else
           zfrac4 = 1.d0
        end if
     else
        zfrac4 = 0.d0
     end if
     
     !! Make sure to add same fluxes on i and i-1, and j and j-1                                                                           
     !! so that divergence in grid cells remain unchanged                                                        
     if (abs(zfrac1*upr(1,1)) <= abs(zfrac2*upr(2,1))) then
        upr(1,:) = zfrac1*upr(1,:)
        upr(2,:) = upr(1,:)
     else
        upr(1,:) = zfrac2*upr(2,:)
        upr(2,:) = zfrac2*upr(2,:)
     end if
     if (abs(zfrac3*upr(3,1)) <= abs(zfrac4*upr(4,1))) then
        upr(3,:) = zfrac3*upr(3,:)
        upr(4,:) = upr(3,:)
     else
        upr(3,:) = zfrac4*upr(4,:)
        upr(4,:) = zfrac4*upr(4,:)
     end if
     upr(5,:) = 0.
     upr(6,:) = 0.
     upr(7,:) = upr(1,:)
     upr(8,:) = upr(2,:)
     upr(9,:) = upr(3,:)
     upr(10,:) = upr(4,:)
     upr(11,:) = 0.
     upr(12,:) = 0.
     
     !! set turbulent velocities to zero
     !! if original velocity is zero
     if (zuu == 0.d0 .or. zum == 0.d0) then
        upr(1,:) = 0.
        upr(2,:) = 0.
        upr(7,:) = 0.
        upr(8,:) = 0.
     end if
     if (zvv == 0.d0 .or. zvm == 0.d0) then
        upr(3,:) = 0.
        upr(4,:) = 0.
        upr(9,:) = 0.
        upr(10,:) = 0.
     end if

    return
  end subroutine active_niter
  
  subroutine active_niter_2
     
     real(dp) :: zx, zy
     
     !! calculate Laplacian of u and v                                         
     !! Interpolate in time and x for u, and y for v                     
     !! so that even if particle says in the same grid cell, d lap(u)/dt is not zero          
     zx = x1 - dble(ibm) ! between 0 and 1                                         
     zy = y1 - dble(jb-1)
     !! d/dt laplacian(u) 
     dt_n(ntrac) = ds * dxyz 
 
     return    
  end subroutine active_niter_2
  

END MODULE mod_active_particles



