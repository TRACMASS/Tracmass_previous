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
  USE mod_param, only:
  USE mod_name, only:
  USE mod_time, only:
  USE mod_loopvars
  USE mod_grid
  USE mod_buoyancy
  USE mod_seed
  USE mod_domain
  USE mod_vel
  USE mod_traj
  USE mod_write
  USE mod_pos
  USE mod_print
  ! === Selectable moules ===
  USE mod_turb
  USE mod_streamfunctions
  USE mod_tracer
  USE mod_sed

  IMPLICIT none

  ! === Loop variables ===
  INTEGER                                    :: i,  j,  k, l, m
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
 print *,'rerun with initial points from ', & 
      trim(outDataDir)//trim(outDataFile)//'_rerun.asc'
 open(67,file=trim(outDataDir)//trim(outDataFile)//'_rerun.asc')
40 continue
 read(67,566,end=41,err=41) ntrac,niter,rlon,rlat,z1,tt,t0,subvol,temp,salt,dens
!  print 566, ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens

#if defined orca025
 if(rlat == dble(jenn(1))) then
    nrj(ntrac,8)=0     ! Southern boundary
    i=i+1
    nout=nout+1
 elseif(rlon == dble(iene(3))) then
    nrj(ntrac,8)=0    ! East
     j=j+1
     nout=nout+1
  elseif(temp > tmaxe .and. salt < smine .and. tt-t0>365.) then
     nrj(ntrac,8)=1    ! back to the warm pool
     k=k+1
  else
     nrj(ntrac,8)=0 
     l=l+1   
     nout=nout+1
  endif
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#elif defined orca025L75
  if( tt-t0 >365.*1000. .and. temp > tmax0 ) then
     nrj(ntrac,8)=1     ! warm end points
     i=i+1
  elseif( tt-t0 >365.*1000. .and. temp <= tmax0 ) then
     nrj(ntrac,8)=2    ! cold end points
     j=j+1
  else                 ! too short
     nrj(ntrac,8)=0 
     k=k+1
     nout=nout+1
  endif
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )

#elif defined ifs
  if(rlat == dble(jenn(1))) then
     nrj(ntrac,8)=1    ! Southern boundary
  elseif(rlat == dble(jens(2))) then
     nrj(ntrac,8)=2    ! Northern boundary 
  else
     nrj(ntrac,8)=0
     print 566,ntrac,niter,rlon,rlat,zz
     stop 4957
  endif /*orc*/
566 format(i8,i7,2f8.2,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#endif
 
  goto 40
41 continue
  m=i+j+k+l
  print *,'Lagrangian decomposition distribution in %: ',     &
       100.*float(i)/float(i+j+k+l),100.*float(j)/float(i+j+k+l),  &
       100.*float(k)/float(i+j+k+l),100.*float(l)/float(i+j+k+l)

  do ntrac=1,ntracmax ! eliminate the unwanted trajectories
   if(nrj(ntrac,8) == 0) nrj(ntrac,6)=1 
  enddo
  
#else
  lbas=1 ! set to 1 if no rerun
#endif /*rerun*/
  !==========================================================
  !=== read ocean/atmosphere GCM data files               ===
  !==========================================================
  
  call fancyTimer('initialize dataset','start')
  ff=dble(nff)
!  tstep=dble(intstep) 
  ints=intstart
  call updateclock
  call readfields   ! initial dataset
  ntrac=0
  call fancyTimer('initialize dataset','stop')

  !==========================================================
  !==========================================================
  !=== Start main time loop                               ===
  !==========================================================
  !==========================================================
  intsTimeLoop: do ints=intstart+1,intstart+intrun
!  intsTimeLoop: do ints=intstart+nff,intstart+intrun,nff
     call fancyTimer('reading next datafield','start')
     tt = ints*tseas
     if (degrade_counter < 1) call readfields
     degrade_counter = degrade_counter + 1
     if (degrade_counter > degrade_time) degrade_counter = 0
     call fancyTimer('reading next datafield','stop')
     
     !=======================================================
     !=== write stream functions and "particle tracer"    ===
     !=======================================================
     if(mod(ints,120) == 0) then 
      call write_streamfunctions
      call writetracer
     endif

!    intspinCond: if(nff*ints <= nff*(intstart+intspin)) then
    intspinCond: if(ints <= intstart+intspin) then
        call fancyTimer('seeding','start')
        call seed (tt,ts)
        call fancyTimer('seeding','stop')
        t0 = tt
        dt = 0.d0
     end if intspinCond

     !=== Check if the output file should be switched. ===
     call writedata(99)! switch

     !=======================================================
     !=== Loop over all trajectories and calculate        ===
     !=== a new position for this time step.              ===
     !=======================================================
     
     call fancyTimer('advection','start')

     ntracLoop: do ntrac=1,ntractot  
     
        ! === Test if the trajectory is dead   ===
        if(nrj(ntrac,6) == 1) cycle ntracLoop
        
        ! === Read in the position, etc at the === 
        ! === beginning of new time step       ===
        x1     =  trj(ntrac,1)
        y1     =  trj(ntrac,2)
        z1     =  trj(ntrac,3)
        tt     =  trj(ntrac,4)
        subvol =  trj(ntrac,5)
        t0     =  trj(ntrac,7)
        
        ib     =  nrj(ntrac,1)
        jb     =  nrj(ntrac,2)
        kb     =  nrj(ntrac,3)
        niter  =  nrj(ntrac,4)
        ts     =  dble(nrj(ntrac,5))
        tss    =  0.d0
        
#ifdef rerun
        lbas=nrj(ntrac,8)
        if(lbas.lt.1 .or.lbas.gt.LBT) then
           print *,'lbas=',lbas,'ntrac=',ntrac
           print *,'trj(ntrac,:)=',trj(ntrac,:)
           print *,'nrj(ntrac,:)=',nrj(ntrac,:)
           exit intsTimeLoop
        endif
#endif /*rerun*/
#ifdef sediment
        ! === Check if water velocities are === 
        ! === large enough for resuspention ===
        ! === of sedimentated trajectories  ===   
        if( nrj(ntrac,6) == 2 ) then
           call resusp(res,ib,jb,kb)
           if(res) then
              ! === updating model time for  ===
              ! === resuspended trajectories ===
              ts=dble(ints-2)
              nrj(ntrac,5)=ints-2
              tt=tseas*dble(ints-2)
              trj(ntrac,4)=tt
              ! === resuspension to bottom layer ===
              ! === kb same as before            ===
              nrj(ntrac,3)=kb
              z1=z1+0.5d0
              ! z1=z1+0.1  !resusp l?gre i boxen
              trj(ntrac,3)=z1
              ! === change flag to put trajectory back in circulation ===
              nrj(ntrac,6)=0
              nsed=nsed-1
              nsusp=nsusp+1
           else
              cycle ntracLoop 
           endif
        endif
#endif  /*sediment*/    
          ! ===  start loop for each trajectory ===
        scrivi=.true.
        niterLoop: do        
           niter=niter+1 ! iterative step of trajectory
#ifdef sediment
           ! Find settling velocity for active gridbox ===
           call sedvel(temp,dens) 
#endif /*sediment*/
           ! === change velocity fields &  === 
           ! === store trajectory position ===
           if( niter.ne.1 .and. tss == dble(iter) &
                .and. nrj(ntrac,7).ne.1 ) then
              trj(ntrac,1)=x1
              trj(ntrac,2)=y1
              trj(ntrac,3)=z1
              trj(ntrac,4)=tt
              trj(ntrac,5)=subvol
              nrj(ntrac,1)=ib
              nrj(ntrac,2)=jb
              nrj(ntrac,3)=kb
              nrj(ntrac,4)=niter
              nrj(ntrac,5)=idint(ts)
              nrj(ntrac,7)=1
              cycle ntracLoop
           endif
           nrj(ntrac,7)=0
           intrpg = dmod(ts,1.d0) ! time interpolation constant between 0 and 1
           intrpr = 1.d0-intrpg
           if(intrpg.lt.0.d0 .or.intrpg.gt.1.d0) then
              print *,'intrpg=',intrpg
              exit intsTimeLoop
           endif
           
           ! === Cyclic world ocean/atmosphere === 
           IF (ib == 1 .AND. x1 >= DBLE (IMT)) THEN
              x1 = x1 - DBLE(IMT)
           END IF
           
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
           if(ts == dble(idint(ts))) then 
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
           dtreg=dtmin * ( dble(int(tt/tseas*dble(iter))) +  & 
                1.d0 - tt/tseas*dble(iter) )
           dt=dtreg
           dsmin=dt/dxyz
#else
           dsmin=dtmin/dxyz
#endif /*regulardt*/ 

           call turbuflux(ia,ja,ka,dt)
           ! === calculate the vertical velocity ===
           call vertvel(ia,iam,ja,ka)
#ifdef timeanalyt
           ss0=dble(idint(ts))*tseas/dxyz
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
           ! === calculate the new positions of the particle ===    
           call pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1)
           !call errorCheck('longjump', errCode)

           ! === north fold cyclic for the ORCA grids ===
#if defined orc || orca1 || orca12 
            if( y1 == dble(JMT-1) ) then ! North fold for ntrac
              x1 = dble(IMT+2) - x1
              ib=idint(x1)+1
              jb=JMT-1
              x0=x1 ; y0=y1 ; ia=ib ; ja=jb
 !             print *,'Changed to',ntrac,ib,jb,kb,x1,y1,z1,kmt(ib,jb+1),kmt(ib,jb)
           elseif(y1 > dble(JMT-1)) then
              print *,'north of northfold for ntrac=',ntrac
              print *,ia,ib,x0,x1
              print *,ja,jb,y0,y1
              print *,ka,kb,z0,z1
              print *,ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
              nerror=nerror+1
              nrj(ntrac,6)=1
              cycle ntracLoop
              !stop 4967
           endif
#elif defined orca025 || orca025L75
           if( y1 == dble(JMT-1) ) then
 !              print *,'North fold for',ntrac
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
              nrj(ntrac,6)=1
              cycle ntracLoop
           endif
#endif
           ! === Cyclic world ocean/atmosphere === 
           if(x1 <  0.d0    ) x1=x1+dble(IMT)       
           if(x1 > dble(IMT)) x1=x1-dble(IMT)   
           IF (ib == 1 .AND. x1 >= DBLE (IMT)) THEN
            x1 = x1 - DBLE(IMT)
           endif    
           if(ib > IMT      ) ib=ib-IMT 
            
           ! === make sure that trajectory ===
           ! === is inside ib,jb,kb box    ===
           if(x1 /= dble(idint(x1))) ib=idint(x1)+1 
           if(y1 /= dble(idint(y1))) jb=idint(y1)+1
           if(z1 /= dble(idint(z1))) kb=idint(z1)+1 

           call errorCheck('boundError', errCode)
           if (errCode.ne.0) cycle ntracLoop

           call errorCheck('landError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           
           call errorCheck('bottomError', errCode)
           if (errCode.ne.0) cycle ntracLoop

           call errorCheck('airborneError', errCode)
           call errorCheck('corrdepthError', errCode)
           call errorCheck('cornerError', errCode)
           ! === diffusion, which adds a random position ===
           ! === position to the new trajectory          ===
#if defined diffusion     
           call diffuse(x1,y1,z1,ib,jb,kb,dt)
#endif
           ! === end trajectory if outside chosen domain === 
           LBTloop: do k=1,LBT
              if(dble(ienw(k)) <= x1 .and. x1 <= dble(iene(k)) .and. &
                 dble(jens(k)) <= y1 .and. y1 <= dble(jenn(k))  ) then
                 nexit(k)=nexit(k)+1
                 exit niterLoop                                
              endif
           enddo LBTLOOP
           if (x1 < 0) exit niterloop

           
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
        call writedata(17) !out
        nrj(ntrac,6)=1
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
       INTEGER                             :: strict  = 0
       INTEGER,intent(out)                 :: errCode
       REAL, save                          :: dxmax = 0, dymax = 0
       INTEGER, save                       :: dxntrac, dyntrac
       CHARACTER(79)                       :: thinline, thickline
       
       thickline = "===============================================" // &
                   "==============================================="
       thinline  = "-----------------------------------------------" // &
                   "-----------------------------------------------"
       errCode = 0
       !return
       select case (trim(teststr))
       case ('ntracGTntracmax')
          if(ntrac.gt.ntracmax) then
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
             nrj(ntrac,6)=1
          endif          

       case ('boundError')
          if(ia>imt .or. ib>imt .or. ja>jmt .or. jb>jmt &
               .or. ia<1 .or. ib<1 .or. ja<1 .or. jb<1) then
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
             nrj(ntrac,6)=1
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
             nrj(ntrac,6)=1
             if (strict==1) stop 
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
       case ('infLoopError')
          if(niter-nrj(ntrac,4).gt.30000) then ! break infinite loops
             if (verbose == 1) then
                print *, thickline !========================================
                print *,'Warning: Particle in infinite loop '
                print *, thinline !-----------------------------------------
                print '(A,I7.7,A,I6.6,A,I7.7)', ' ntrac : ', ntrac,     & 
                                          ' niter : ', niter,     &
                                          '    nrj : ', nrj(ntrac,4)
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
             trj(ntrac,1)=x1
             trj(ntrac,2)=y1
             trj(ntrac,3)=z1
             trj(ntrac,4)=tt
             trj(ntrac,5)=subvol
             nrj(ntrac,1)=ib
             nrj(ntrac,2)=jb
             nrj(ntrac,3)=kb
             nrj(ntrac,4)=niter
             nrj(ntrac,5)=idint(ts)
             nrj(ntrac,6) = 1  ! 0=continue trajectory, 1=end trajectory
             nrj(ntrac,7)=1
             nloop=nloop+1             
             errCode = -48
          end if
       case ('bottomError')
          ! if trajectory under bottom of ocean, 
          ! then put in middle of deepest layer 
          ! (this should however be impossible)
           if( z1.le.dble(KM-kmt(ib,jb)) ) then
              print *,'under bottom !!!!!!!',z1,dble(KM-kmt(ib,jb))
              print *,'kmt=',kmt(ia,ja),kmt(ib,jb)
              print *,'ntrac=',ntrac,niter 
              call print_ds
              call print_pos
              call cross_stat(1,ia,ja,ka,x0,dse,dsw) ! zonal
              call cross_stat(2,ia,ja,ka,y0,dsn,dss) ! meridional
              call cross_stat(3,ia,ja,ka,z0,dsu,dsd) ! vertical
              print *,'time step sol:',dse,dsw,dsn,dss,dsu,dsd
              nerror=nerror+1
              nrj(ntrac,6)=1
 !             stop 3957
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
              nrj(ntrac,6)=1
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
#if defined zgrid3Dt || defined zgrid3D
         '    dz(k) : ', dz(kb), '   dzt :  ', dzt(ib,jb,kb)
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
!!$    
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

