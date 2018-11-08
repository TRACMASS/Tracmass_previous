SUBROUTINE loop
!!---------------------------------------------------------------------------
!!
!!       SUBROUTINE loop:
!!
!!          The main loop where new trajectory positions are 
!!          calculated, i.e. trj and nrj are updated each time step.
!!          Sets the flags for reruns.
!!
!!          Contains subroutines for computing the grid box volume, time,
!!          and writing data to files. 
!!
!!          See tracmass manual for schematic of the structure.
!!
!!
!!---------------------------------------------------------------------------          
  USE mod_param,    only: ntracmax, undef, tday 
  USE mod_loopvars, only: dse, dsw, dsmin, ds, dsu, dsd, dsn, dss, &
                          niter, lbas, scrivi, subvol
  USE mod_grid,     only: imt, jmt, km, kmt, dyu, dxv, dxdy, dxyz, dz, dzt, &
                          mask, iter, nsm, nsp, hs, calc_dxyz, nperio !joakim
  use mod_vel,      only: uflux, vflux, wflux
  USE mod_seed,     only: ff, nff, seedTime, seed
  USE mod_domain,   only: timax, jens, jenn, iene, ienw
  USE mod_vel,      only: degrade_counter, degrade_time
  USE mod_write,    only: writedata
  USE mod_pos,      only: pos
  USE mod_print,    only: print_start_loop, print_cycle_loop, print_end_loop
  USE mod_tempsalt
  ! === Selectable moules ===
  USE mod_active_particles
  USE mod_streamfunctions, only: intpsi
  USE mod_tracer
  USE mod_sed
  
  USE mod_deformation
  
  IMPLICIT none

  ! === Loop variables ===
  INTEGER                                    :: i,  j,  k, l, m, n
  ! === Variables to interpolate fields ===
  REAL                                       :: temp, salt, dens
  REAL                                       :: temp2, salt2, dens2
  ! === Error Evaluation ===
  INTEGER                                    :: errCode
  INTEGER                                    :: landError=0, boundError=0
  REAL                                       :: zz
  
  call print_start_loop
  
  dstep = 1.d0 / dble(iter)
  dtmin = dstep * tseas
    
  !==========================================================
  !===   Read in end positions from a previous run        === 
  !==========================================================
  
#ifdef rerun
 I=0 ; j=0 ; k=0 ; l=0
 !print *,'rerun with initial points from ', & 
!      trim(outDataDir)//trim(outDataFile)//'_rerun.asc'
! open(67,file=trim(outDataDir)//trim(outDataFile)//'_rerun.asc')
40 continue
 read(67,566,end=41,err=41) ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens
!  print 566, ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens

#if defined orca025 || orca1 
  if(    y1 == dble(jenn(1))) then
     trajectories(ntrac)%lbas = 1 ! Southern boundary
     i=i+1
     nout=nout+1
  elseif(y1 == dble(jens(2))) then
     trajectories(ntrac)%lbas = 2 ! Northern boundary
     j=j+1
     nout=nout+1
!  elseif(temp > tmaxe .and. salt < smine .and. tt-t0>365.) then
!     nrj(8,ntrac)=1    ! back to the warm pool
!     k=k+1
  else
     trajectories(ntrac)%lbas = 0
     l=l+1   
     nout=nout+1
  endif
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )


#elif defined orca025L75
  if( tt-t0 >365.*1000. .and. temp > tmax0 ) then
     trajectories(ntrac)%lbas = 1 ! warm end points
     i=i+1
  elseif( tt-t0 >365.*1000. .and. temp <= tmax0 ) then
     trajectories(ntrac)%lbas = 2 ! cold end points
     j=j+1
  else                 ! too short
     trajectories(ntrac)%lbas = 0
     k=k+1
     nout=nout+1
  endif
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )

#elif defined ifs
  nendLoop: do n=1,nend
     if( dble(ienw(n)) <= x1 .and. x1 <= dble(iene(n)) .and. &
         dble(jens(n)) <= y1 .and. y1 <= dble(jenn(n))  ) then  
        trajectories(ntrac)%lbas = n
        dist(n) = dist(n) + 1
        cycle lbasLoop
     endif
  enddo nendLoop
  
  if( trajectories(ntrac)%lbas == 0 ) then
     print 566,ntrac,niter,x1,y1,zz
     stop 4957
  endif
  
566 format(i8,i7,2f8.2,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#endif
 
  goto 40
41 continue
#ifdef ifs
  print *,'Lagrangian decomposition distribution in %: '
  do n=1,nend
     PRINT*,100.*float(dist(n))/float(sum(dist))
  enddo
#else
  m=i+j+k+l
  print *,'Lagrangian decomposition distribution in %: ',     &
  100.*float(i)/float(i+j+k+l),100.*float(j)/float(i+j+k+l),  &
  100.*float(k)/float(i+j+k+l),100.*float(l)/float(i+j+k+l)
#endif
  do ntrac=1,ntracmax ! eliminate the unwanted trajectories
   if (trajectories(ntrac)%lbas == 0) trajectories(ntrac)%active=.false.
  enddo
  
#else
  lbas=1 ! set to 1 if no rerun
#endif /*rerun*/
  !==========================================================
  !=== read ocean/atmosphere GCM data files               ===
  !==========================================================
  
  call fancyTimer('initialize dataset','start')
  ff=dble(nff)
  ints = intstart
  call updateclock  
  call readfields   ! initial dataset
  call active_init
  ntrac = 0
  call fancyTimer('initialize dataset','stop')

  !==========================================================
  !=== Start main time loop                               ===
  !==========================================================
  intsTimeLoop: do ints=intstart+nff, intstart+intrun, nff
     call fancyTimer('reading next datafield','start')
     tt = ints*tseas
     if (degrade_counter < 1) call readfields
     degrade_counter = degrade_counter + 1
     if (degrade_counter > degrade_time) degrade_counter = 0
     call fancyTimer('reading next datafield','stop')
     
     !=======================================================
     !=== write stream functions and "particle tracer"    ===
     !=======================================================
     if(mod(ints,intpsi) == 0) then 
      call write_streamfunctions
      call writetracer
     endif

     intspinCond: if(ints*nff <= (intstart+intspin)*nff) then
        call fancyTimer('seeding','start')
        call seed (tt,ts)
        call fancyTimer('seeding','stop')
        t0 = tt
        dt = 0.d0
     end if intspinCond

     call active_ints(ints)
     !=== Check if the output file should be switched. ===
     call writedata(99)! switch

     !=======================================================
     !=== Loop over all trajectories and calculate        ===
     !=== a new position for this time step.              ===
     !=======================================================
     
     call fancyTimer('advection','start')
     
     ntracLoop: do ntrac=1,ntractot
        !print *,ntrac, ntractot
        ! === Test if the trajectory is dead   ===
        if (trajectories(ntrac)%active .eqv. .false.) cycle ntracLoop
        
        ! === Read in the position, etc at the === 
        ! === beginning of new time step       ===
        x1     = trajectories(ntrac)%x1
        y1     = trajectories(ntrac)%y1
        z1     = trajectories(ntrac)%z1
        tt     = trajectories(ntrac)%tt
        subvol = trajectories(ntrac)%subvol
        t0     = trajectories(ntrac)%t0
        
        ib     = trajectories(ntrac)%ib
        jb     = trajectories(ntrac)%jb
        kb     = trajectories(ntrac)%kb
        niter  = trajectories(ntrac)%niter
        ts     = DBLE(trajectories(ntrac)%nts)
        tss    = 0.d0
        
        !! at start of a new time step, calc laplacian u and v
        !! and remember the old value
        !! The next four lines are laplacian of volume fluxes
        !trajectories(ntrac)%lapu1   = trajectories(ntrac)%lapu2
        !trajectories(ntrac)%lapv1   = trajectories(ntrac)%lapv2
        !trajectories(ntrac)%lapu2   = lapu(ib,jb,kb,1) 
        !trajectories(ntrac)%lapv2   = lapv(ib,jb,kb,1) 
        !! Here we translate laplacian of volume flux to velocity
        lapu1  = trajectories(ntrac)%lapu1 / (dyu(ib,jb) * dzt(ib,jb,kb,1))   
        lapu2  = trajectories(ntrac)%lapu2 / (dyu(ib,jb) * dzt(ib,jb,kb,1))
        lapv1  = trajectories(ntrac)%lapv1 / (dxv(ib,jb) * dzt(ib,jb,kb,1))      
        lapv2  = trajectories(ntrac)%lapv2 / (dxv(ib,jb) * dzt(ib,jb,kb,1))
        !dlapu = trajectories(ntrac)%dlapu
        !dlapv = trajectories(ntrac)%dlapv
        !dhdiv  = hdiv2 - hdiv1
        !dvort  = vort2 - vort1
        !trajectories(ntrac)%dlapu = dlapu
        !trajectories(ntrac)%dlapv = dlapv
        ! put new step statistics as current
        !trajectories(ntrac)%lapu1 = lapu2
        !trajectories(ntrac)%lapv1 = lapv2
        !trajectories(ntrac)%hdiv1 = hdiv2
        !trajectories(ntrac)%vort1 = vort2
        
#ifdef rerun
        lbas = trajectories(ntrac)%lbas
        if(lbas.lt.1 .or.lbas.gt.nend) then
           print *,'lbas=',lbas,'ntrac=',ntrac
           print *,'x1=',trajectories(ntrac)%x1
           print *,'y1=',trajectories(ntrac)%y1
           print *,'z1=',trajectories(ntrac)%z1
           print *,'tt=',trajectories(ntrac)%tt
           print *,'t0=',trajectories(ntrac)%t0
           print *,'subvol=',trajectories(ntrac)%subvol
           print *,'ib=',trajectories(ntrac)%ib
           print *,'jb=',trajectories(ntrac)%jb
           print *,'kb=',trajectories(ntrac)%kb
           print *,'nts=',trajectories(ntrac)%nts
           print *,'icycle=',trajectories(ntrac)%icycle
           print *,'active=',trajectories(ntrac)%active
           exit intsTimeLoop
        endif
#endif /*rerun*/

        call active_ntrac(ntrac)
          ! ===  start loop for each trajectory ===
        scrivi=.true.
        niterLoop: do           
           niter=niter+1 ! iterative step of trajectory
           ! === change velocity fields &  === 
           ! === store trajectory position ===
           if( niter /= 1 .and. tss == dble(iter) &
                .and. trajectories(ntrac)%icycle /= 1 ) then
              
              trajectories(ntrac)%x1     = x1
              trajectories(ntrac)%y1     = y1
              trajectories(ntrac)%z1     = z1
              trajectories(ntrac)%tt     = tt
              trajectories(ntrac)%subvol = subvol
              trajectories(ntrac)%ib     = ib
              trajectories(ntrac)%jb     = jb
              trajectories(ntrac)%kb     = kb
              trajectories(ntrac)%niter  = niter
              trajectories(ntrac)%nts    = idint(ts)
              trajectories(ntrac)%icycle = 1
              
              !trajectories(ntrac)%lapu2   = lapu(ib,jb,kb,2)
              !trajectories(ntrac)%lapv2   = lapv(ib,jb,kb,2)
              !trajectories(ntrac)%hdiv2   = hdiv(ib,jb,kb,2)
              !trajectories(ntrac)%vort2   = vort(ib,jb,kb,2)
              
              cycle ntracLoop
           endif
           
           trajectories(ntrac)%icycle = 0
#if defined fixedtimestep 
           intrpg = 0.d0  ! mimics Ariane's lack of linear interpolation of the velocity fields
#else
           intrpg = dmod(ts,1.d0) ! time interpolation constant between 0 and 1
#endif
           intrpr = 1.d0-intrpg
           if(intrpg.lt.0.d0 .or.intrpg.gt.1.d0) then
              print *,'intrpg=',intrpg
              exit intsTimeLoop
           endif
           
           if (nperio /= 0) then
              ! === Cyclic world ocean/atmosphere === 
              IF (ib == 1 .AND. x1 >= DBLE (IMT)) THEN
                 x1 = x1 - DBLE(IMT)
              END IF
           end if
           
           x0  = x1
           y0  = y1
           z0  = z1
           ia  = ib
           iam = ia-1
           if(iam == 0) iam = IMT
           ja  = jb
           ka  = kb
           
           call calc_dxyz(intrpr, intrpg)
           call errorCheck('dxyzError'     ,errCode)
           call errorCheck('coordBoxError' ,errCode)
           call errorCheck('infLoopError'  ,errCode)
           if (errCode.ne.0) cycle ntracLoop

           ! === write trajectory ===                       
#ifdef tracer
           if(ts == dble(int(ts, 8))) then 
              tra(ia,ja,ka)=tra(ia,ja,ka)+real(subvol)
           end if
#endif /*tracer*/
           call writedata(11)
           
           !==============================================! 
           ! calculate the 3 crossing times over the box  ! 
           ! choose the shortest time and calculate the   !
           ! new positions                                !
           !                                              !
           !-- solving the differential equations ---     !
           ! note:                                        !
           ! space variables (x,...) are dimensionless    !
           ! time variables (ds,...) are in seconds/m^3   !
           !==============================================! 
#ifdef regulardt
           dtreg=dtmin * ( dble(int(tt/tseas*dble(iter),8)) +  & 
                1.d0 - tt/tseas*dble(iter) )
           dt=dtreg
           dsmin=dt/dxyz
#else
           dsmin=dtmin/dxyz
#endif /*regulardt*/ 
           call active_niter 
           
           !call turbuflux(ia,ja,ka,dt)
           ! === calculate the vertical velocity ===
           call vertvel(ia,iam,ja,ka)
#ifdef timeanalyt
!           ss0=dble(int(ts,8))*tseas/dxyz or should ssp be in the call cross_time?
           call cross_time(1,ia,ja,ka,x0,dse,dsw) ! zonal
           call cross_time(2,ia,ja,ka,y0,dsn,dss) ! merid
           call cross_time(3,ia,ja,ka,z0,dsu,dsd) ! vert
#else
           call cross_stat(1,ia,ja,ka,x0,dse,dsw) ! zonal
           call cross_stat(2,ia,ja,ka,y0,dsn,dss) ! meridional
           call cross_stat(3,ia,ja,ka,z0,dsu,dsd) ! vertical
#endif /*timeanalyt*/

           ds = min(dse, dsw, dsn, dss, dsu, dsd, dsmin)
           call errorCheck('dsCrossError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           call calc_time
           
           if (ntrac == 365) print*,'ds',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
           if (ntrac == 365) print*,'hej1',x0,y0,z0,ia,ja,ka
           
           ! === calculate the new positions of the particle ===    
           call pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1)
           !call errorCheck('longjump', errCode)
           if (ntrac == 365) print*,'hej2',x1,y1,z1,ib,jb,kb

           if (nperio == 6) then
              ! === north fold cyclic for the ORCA grids ===
              if( y1 == dble(JMT-1) ) then ! North fold for ntrac
                 x1 = dble(IMT+2) - x1
                 ib=idint(x1)+1
                 jb=JMT-1
                 x0=x1 ; y0=y1 ; ia=ib ; ja=jb
              elseif(y1 > dble(JMT-1)) then
!                print *,'north of northfold for ntrac=',ntrac
                x1 = dble(IMT+2) - x1
                ib=idint(x1)+1
                jb=JMT-1
                y1= dble(JMT-1) -y1 + dble(JMT-1)
                x0=x1 ; y0=y1 ; ia=ib ; ja=jb
              endif

           else if (nperio == 4) then
              ! === another north fold implementation 
              if( y1 == dble(JMT-1) ) then
                 x1 = dble(IMT+3) - x1
                 y1 = dble(JMT-2)
                 ib=idint(x1)
                 jb=JMT-2
                 x0=x1 ; y0=y1 ; ia=ib ; ja=jb
              elseif(y1 > dble(JMT-1)) then
                 print *,ia,ib,x0,x1
                 print *,ja,jb,y0,y1
                 print *,ka,kb,z0,z1
                 print *,ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
                 nerror=nerror+1
                 trajectories(ntrac)%active = .false.
                 cycle ntracLoop
              endif
              
              ! === Cyclic Arctic in a global cylindrical projection ===
              if( y1 == dble(JMT-1) ) then ! North fold for ntrac
                 x1 = dble(IMT+2) - x1
                 ib=idint(x1)+1
                 jb=JMT-1
                 x0=x1 ; y0=y1 ; ia=ib ; ja=jb
              elseif(y1 > dble(JMT-1)) then
                 print *,'north of northfold for ntrac=',ntrac
                 x1 = dble(IMT+2) - x1
                 ib=idint(x1)+1
                 jb=JMT-1
                 y1= dble(JMT-1) -y1 + dble(JMT-1)
                 x0=x1 ; y0=y1 ; ia=ib ; ja=jb
              endif
           end if
           
           if (nperio /= 0) then
              ! === East-west cyclic 
              if(x1 <  0.d0    ) then
                 print*,'<0',ntrac,x1
                 x1=x1+dble(IMT)       
                 print*,ntrac,x1
              end if
              if(x1 > dble(IMT)) then
                 print*,'>imt',ntrac,x1
                 x1=x1-dble(IMT)   
                 print*,ntrac,x1
              end if
              IF (ib == 1 .AND. x1 >= DBLE (IMT)) THEN
                 x1 = x1 - DBLE(IMT)
              endif
              if(ib > IMT      ) ib=ib-IMT 
           end if
           
           ! === make sure that trajectory ===
           ! === is inside ib,jb,kb box    ===
           if(x1 /= dble(idint(x1))) ib=idint(x1)+1 
           if(y1 /= dble(idint(y1))) jb=idint(y1)+1
           if(z1 /= dble(idint(z1))) kb=idint(z1)+1 
           
           if (ja>jmt) ja = jmt - (ja - jmt)
           if (jb>jmt) jb = jmt - (jb - jmt)
           
           call active_niter_2 !!joakim edit
           
           if (ntrac == 365) print*,'hej3',x1,y1,z1,ib,jb,kb
           
           call errorCheck('boundError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           call errorCheck('landError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           call errorCheck('bottomError', errCode)
       !    if (errCode.ne.0) cycle ntracLoop
           call errorCheck('airborneError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           
           call errorCheck('corrdepthError', errCode)
!           if (errCode.ne.0) cycle ntracLoop
           call errorCheck('cornerError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           
           ! === diffusion, which adds a random position ===
           ! === position to the new trajectory          ===
#if defined diffusion     
           call diffuse(x1,y1,z1,ib,jb,kb,dt)
#endif
           ! === end trajectory if outside chosen domain === 
           nendloop: do k=1,nend
              if(ienw(k) <= x1 .and. x1 <= iene(k) .and. &
                 jens(k) <= y1 .and. y1 <= jenn(k) ) then
                 nexit(k)=nexit(k)+1
                 exit niterLoop                                
              endif
           enddo nendloop
           
          if(ia>imt .or. ib>imt .or. ja>jmt .or. jb>jmt &
               .or. ia<1 .or. ib<1 .or. ja<1 .or. jb<1) then
             print *,'Warning: Trajectory leaving model area'
             call writedata(17)
             trajectories(ntrac)%active = .false.
             exit niterLoop                                
          end if







           
           
#if defined tempsalt
           call interp (ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
           ! if (temp < tmine .or. temp > tmaxe .or. &
           ! &   salt < smine .or. salt > smaxe .or. &
           ! &   dens < rmine .or. dens > rmaxe      ) then
                if (temp > tmaxe .and. salt < smine .and.  &
               &   (tt-t0)/tday > 365.      ) then
                 nexit(NEND)=nexit(NEND)+1
                 exit niterLoop                                
               endif
#endif       
           ! === stop trajectory if the choosen time or ===
           ! === water mass properties are exceeded     ===
           if(tt-t0.gt.timax) then
              nexit(NEND)=nexit(NEND)+1
              exit niterLoop
           endif
        end do niterLoop

        nout=nout+1
        call writedata(17)
        trajectories(ntrac)%active = .false.
     end do ntracLoop
 
     call print_cycle_loop()

     IF (ntractot /= 0 .AND. ntractot - nout - nerror == 0  .AND. &
          seedTime /= 2) THEN
        EXIT intsTimeLoop
     END IF
     
  end do intsTimeLoop
  call print_end_loop
  
return
     
     



























   CONTAINS
     
     subroutine errorCheck(teststr,errCode)
       CHARACTER (len=*),intent(in)        :: teststr    
       INTEGER                             :: verbose = 1
       INTEGER                             :: strict  = 1
       INTEGER,intent(out)                 :: errCode
       REAL, save                          :: dxmax = 0, dymax = 0
       INTEGER, save                       :: dxntrac, dyntrac
       CHARACTER(79)                       :: thinline, thickline
       
       thickline = "===============================================" // &
                   "==============================================="
       thinline  = "-----------------------------------------------" // &
                   "-----------------------------------------------"
       errCode = 0
       
       select case (trim(teststr))
       case ('infLoopError')
          if(niter-trajectories(ntrac)%niter > 30000) then ! break infinite loops
             if (verbose == 2) then
                print *, thickline !========================================
                print *,'Warning: Particle in infinite loop '
                print *, thinline !-----------------------------------------
                print '(A,I7.7,A,I6.6,A,I7.7)', ' ntrac : ', ntrac,     & 
                                          ' niter : ', niter,     &
                                          '    trajectories(ntrac)%niter : ', trajectories(ntrac)%niter
                call print_grd
                call print_pos

                print '(A,F7.0,A,F7.0,A,F7.0,A,F7.0)',            &
                     ' ufl(ia) : ',(intrpbg*uflux(ia ,ja,ka,nsp) +    &
                                    intrpb*uflux(ia ,ja,ka,nsm))*ff,  &
                     ' ufl(ib) : ', (intrpbg*uflux(iam,ja,ka,nsp) +   & 
                                    intrpb*uflux(iam,ja,ka,nsm))*ff,  &
                     ' vfl(ja) : ', (intrpbg*vflux(ia,ja  ,ka,nsp) +  & 
                                    intrpb*vflux(ia,ja  ,ka,nsm))*ff, &
                     ' vfl(jb) : ', (intrpbg*vflux(ia,ja-1,ka,nsp) +  & 
                                    intrpb*vflux(ia,ja-1,ka,nsm))*ff 
                print *, thinline !-----------------------------------------
             end if
!             z1=dble(kb)-0.5d0
             trajectories(ntrac)%x1     = x1 !dble(ib)-0.5d0  !x1
             trajectories(ntrac)%y1     = y1 !dble(jb)-0.5d0 ! y1
             trajectories(ntrac)%z1     = z1 !dble(kb)-0.5d0 ! z1
             trajectories(ntrac)%tt     = tt
             trajectories(ntrac)%subvol = subvol
             trajectories(ntrac)%ib     = ib
             trajectories(ntrac)%jb     = jb
             trajectories(ntrac)%kb     = kb
             trajectories(ntrac)%niter  = niter
             trajectories(ntrac)%nts    = idint(ts)
             trajectories(ntrac)%active = .false.  ! 0=continue trajectory, 1=end trajectory
             trajectories(ntrac)%icycle = 1
             
             nloop=nloop+1             
             errCode = -48
          end if
       case ('ntracGTntracmax')
          if(ntrac > ntracmax) then
             print *, thickline !========================================
             print *,'ERROR: to many trajectories,'
             print *, thinline !-----------------------------------------
             print *,'increase ntracmax since'
             print *,'ntrac >',ntrac
             print *,'when ints=',ints,' and ' 
             print *,'intspin=',intspin
             print *,',(intspin-ints)/ints*ntrac='
             print *,(intspin-ints)/ints*ntrac
             print *, thinline !-----------------------------------------
             print *,'The run is terminated'
             print *, thickline !========================================
             errCode = -38
             stop
          endif
          
       case ('dxyzError')
          if(dxyz == 0.d0) then
             if (verbose == 1) then                 
                print *, thickline !========================================
                print *,'ERROR: dxyz is zero'
                print *, thinline !-----------------------------------------
                print *,'ntrac=',ntrac,' ints=', ints
                call print_pos
                print *,'kmt=',kmt(ib,jb)
                print *,'dz=',dz(kb)
                print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb)
                print *, thinline !-----------------------------------------
                print *,'The trajectory is killed'
                print *, thickline !========================================
             end if
             nerror=nerror+1
             errCode = -39
             if (strict==1) stop 40961
             call writedata(40)
             trajectories(ntrac)%active = .false.
          endif          

       case ('boundError')
          if(ia<1 .or. ia>imt .or. ib<1 .or. ib>imt .or.    &
             ja<1 .or. ja>jmt .or. jb<1 .or. jb>jmt .or.    &
             y0<1 .or. y0>jmt .or. y1<1 .or. y1>jmt         &
             ) then
             if (verbose == 1) then
                print *, thickline !========================================
                print *,'Warning: Trajectory leaving model area'
                print *, thinline !-----------------------------------------
                call print_pos
                call print_ds
                print *,'tt=',tt,ts
                print *,'ntrac=',ntrac
                print *, thinline !-----------------------------------------
                print *,'The trajectory is killed'
                print *, thickline !========================================
             end if
             call writedata(40)
             nerror=nerror+1
             boundError = boundError +1
             errCode = -50
             if (strict==1) stop
             call writedata(40)
             trajectories(ntrac)%active = .false.
          endif

       case ('landError')
          if(kmt(ib,jb) == 0) then
             if (verbose == 1) then
                print *, thickline !========================================
                print *,'Warning: Trajectory on land'
                print *, thinline !-----------------------------------------
                call print_pos
                call print_ds
                print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb),dxdy(ia,ja)
                print *,'hs=',hs(ia,ja,nsm),hs(ia,ja,nsp),hs(ib,jb,nsm),hs(ib,jb,nsp)
                print *,'uflux=',uflux(ia,ja,ka,nsm),uflux(ia,ja,ka,nsp),uflux(ib,jb,kb,nsm),uflux(ib,jb,kb,nsp)
                print *,'vflux=',vflux(ia,ja,ka,nsm),vflux(ia,ja,ka,nsp),vflux(ib,jb,kb,nsm),vflux(ib,jb,kb,nsp)
                print *,'kmt=',kmt(ia-1,ja),kmt(ia,ja),kmt(ia+1,ja),kmt(ia,ja-1),kmt(ia,ja+1)
                print *,'tt=',tt,ts,tt/tday,t0/tday
                print *,'ntrac=',ntrac
                print *,'niter=',niter
#ifdef turb
                print *,'upr=',upr
#endif
                print *, thinline !-----------------------------------------
                print *,'The trajectory is killed'
                print *, thickline !========================================
             end if
             nerror=nerror+1
             landError = landError +1
             errCode = -40             
             call writedata(40)
             trajectories(ntrac)%active = .false.
             !if (strict==1) stop  !!joakim edit
          endif
          case ('coordboxError')
          ! ===  Check that coordinates belongs to   ===
          ! ===  correct box. Valuable for debugging ===
          if( dble(ib-1).gt.x1 .or. dble(ib).lt.x1 )  then
             print *, thickline !========================================
             print *,'ERROR: Particle overshoot in i direction'
             print *, thinline !-----------------------------------------
             call print_pos
             print *, thinline !-----------------------------------------
             print *,'The run is terminated'
             print *, thickline !========================================
             errCode = -42
             stop
          elseif( dble(jb-1).gt.y1 .or. dble(jb).lt.y1 )  then
             print *, thickline !========================================
             print *,'ERROR: Particle overshoot in j direction'
             print *, thinline !-----------------------------------------
             call print_pos
             print *, thinline !-----------------------------------------
             print *,'The run is terminated'
             print *, thickline !========================================
             errCode = -44
             stop
          elseif((dble(kb-1).gt.z1.and.kb.ne.KM).or. & 
               dble(kb).lt.z1 ) then
             print *, thickline !========================================
             print *,'ERROR: Particle overshoot in k direction'
             print *, thinline !-----------------------------------------
             call print_grd
             call print_pos
             print *, thinline !-----------------------------------------
             print *,'The run is terminated'
             print *, thickline !========================================
             errCode = -46
             stop
          end if
       case ('bottomError')
          ! if trajectory under bottom of ocean, 
          ! then put in middle of deepest layer 
          ! (this can happen when using time dependent vertical coordinates)
           if( z1.le.dble(KM-kmt(ib,jb)) ) then
              print *,'Particle below bottom',z1,dble(KM-kmt(ib,jb))
              print *,'x1,y1',x1,y1
              print *,'ntrac=',ntrac,niter 
              nerror=nerror+1
 !             trajectories(ntrac)%iend = 1
              !stop 3957
              z1=dble(KM-kmt(ib,jb))+0.5d0
              errCode = -49
           end if
        case ('airborneError')
           ! if trajectory above sea level,
           ! then put back in the middle of shallowest layer (evaporation)
           if( z1.ge.dble(KM) ) then
              z1=dble(KM)-0.5d0
              kb=KM
              errCode = -50
           endif
        case ('corrdepthError')
           ! sets the right level for the corresponding trajectory depth
           if(z1.ne.dble(idint(z1))) then
              kb=idint(z1)+1
              if(kb == KM+1) kb=KM  ! (should perhaps be removed)
              errCode = -52
           endif
           case ('cornerError')
              ! problems if trajectory is in the exact location of a corner
           if(x1 == dble(idint(x1)) .and. y1 == dble(idint(y1))) then
              !print *,'corner problem',ntrac,x1,x0,y1,y0,ib,jb
              !print *,'ds=',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
              !stop 34957
              ! corner problems may be solved the following way 
              ! but should really not happen at all
              if(ds == dse .or. ds == dsw) then
                 if(y1.ne.y0) then
                    y1=y0 ; jb=ja
                 else
                    y1=dble(jb)-0.5d0
                 endif
              elseif(ds == dsn .or. ds == dss) then
                 if(y1.ne.y0) then
                    x1=x0 ; ib=ia 
                 else
                    x1=dble(ib)-0.5d0
                 endif
              else
                 x1=dble(ib)-0.5d0
                 y1=dble(jb)-0.5d0
              endif
              errCode = -54
           endif
        case ('dsCrossError')
           ! === Can not find any path for unknown reasons ===
           if(ds == UNDEF .or.ds == 0.d0)then 
              if (verbose == 0) then
                 print *, " "
                 print *, " "
                 print *, thickline !========================================
                 print *,'Warning: not find any path for unknown reason '
                 print *, " "
                 write (*,'(A, E9.3, A, E9.3)'), ' uflux= ', &
                      uflux(ia,ja,ka,nsm),'  vflux= ', vflux(ia,ja,ka,nsm)
                 call print_ds
                 print *,'---------------------------------------------------'
                 print *,"   ntrac = ",ntrac
                 call print_pos
                 write (*,'(A7, I10, A7, I10, A7, I10)'), & 
                      ' k_inv= ', KM+1-kmt(ia,ja), ' kmt= ', kmt(ia,ja), &
                      'lnd= ', mask(ia,ja)
                 print *, thinline !-----------------------------------------
                 print *,'The trajectory is killed'
                 print *, thickline !========================================
              end if
              nerror=nerror+1
              trajectories(ntrac)%active = .false.
              errCode = -56
           end if
        case ('longjump')
           ! === Check if the particles are too large jumps ===
           !print *, 'x0 = ', x0, 'x1 = ', x1, 'dx = ', x1-x0
           if (abs(x1-x0) > dxmax) then 
              print *,'New dxmax: ', dxmax, ' ntrac = ', ntrac
              dxntrac = ntrac
              dxmax = abs(x1-x0)
           end if
           if (abs(y1-y0) > dymax) then 
              print *,'New dymax: ', dymax, ' ntrac = ', ntrac
              dyntrac = ntrac
              dymax = abs(y1-y0)
           end if
           if (dxmax>imt) then
              print *, thickline !========================================
              print *,'dx unrealistic, ntrac = ', ntrac, ' dx = ', dxmax
              print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb),dxdy(ia,ja)
              print *, thinline !-----------------------------------------
              stop
           end if
           if (dymax>jmt) then
              print *, thickline !========================================
              print *,'dy unrealistic, ntrac = ', ntrac, ' dy = ', dymax
              print *,'ds',ds,dse,dsw,dsn,dss,dsu,dsd
              print *,'dsmin=',ds,dsmin,dtmin
              print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb),dxdy(ia,ja)
              print *, thinline !-----------------------------------------
              stop
           end if
        end select
      end subroutine errorCheck

  subroutine print_ds
    print *, '   ds = ', ds 
    print *, '  dse = ', dse
    print *, '  dsw = ', dsw
    print *, '  dsn = ', dsn
    print *, '  dss = ', dss
    print *, '  dsu = ', dsu
    print *, '  dsd = ', dsd
    print *, 'dsmin = ', dsmin
    print *, 'dtmin = ', dtmin
    print *, 'dxyz = ', dxyz
    !write (*,FMT='(A, 5E9.2)'),' ds=',ds,dse,dsw,dsn,dss
    !             write (*,FMT='(4E9.2)'), dsu,dsd,dsmin,dxyz
  end subroutine print_ds

  subroutine print_grd
    print '(A,F16.4, A,E12.3)', '  dxdy :',dxdy(ib,jb), &
         '          dxyz :  ',dxyz
    print '(A,I4,A,F7.2,A,F7.2)',    &
         '    kmt: ', kmt(ib,ja), &
#if defined zgrid3D
         '    dz(k) : ', dz(kb), '   dzt :  ', dzt(ib,jb,kb,1)

#else
         '    dz(k) : ', dz(kb), '   dzt :  ', 0.0
#endif
       end subroutine print_grd

  subroutine print_pos
    print '(A,I4,A,I4,A,I4,A,I4,A,I4,A,I4)', &
         '      ia : ', ia,    '         ib : ', ib, &
         '         ja : ', ja, '         jb : ', jb
    print '(A,F7.2,A,F7.2,A,F7.2,A,F7.2)','      x1 : ', x1, &
         '      x0 : ', x0, '      y1 : ', y1, '      y0 : ', y0 
    print '(A,I4,A,I4)',    '      ka : ', ka, '         kb : ', kb
    print '(A,F7.2,A,F7.2)','      z1 : ', z1, '      z0 : ', z0
  end subroutine print_pos

  subroutine fancyTimer(timerText ,testStr)
    IMPLICIT NONE

    CHARACTER (len=*)                          :: timerText ,testStr
    REAL ,SAVE                                 :: fullstamp1 ,fullstamp2
    REAL ,SAVE ,DIMENSION(2)                   :: timestamp1 ,timestamp2
    REAL                                       :: timeDiff
!!$    select case (trim(testStr))
!!$    case ('start')
!!$       WRITE (6, FMT="(A)", ADVANCE="NO") ,' - Begin '//trim(timerText)
!!$       call etime(timestamp1,fullstamp1)
!!$    case ('stop')
!!$       call etime(timestamp2,fullstamp2)
!!$       timeDiff=fullstamp2-fullstamp1
!!$       write (6 , FMT="(A,F6.1,A)") ', done in ' ,timeDiff ,' sec'
!!$    end select
  end subroutine fancyTimer
end subroutine loop

