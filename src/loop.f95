SUBROUTINE loop
!!------------------------------------------------------------------------------!!
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
!!------------------------------------------------------------------------------          
  USE mod_param
  USE mod_name
  USE mod_time
  USE mod_loopvars
  USE mod_grid
  USE mod_buoyancy
  USE mod_seed
  USE mod_domain
  USE mod_vel
  USE mod_traj
  USE mod_pos
  USE mod_turb
  USE mod_coord
#ifdef tracer
  USE mod_tracer
#endif /*tracer*/
#ifdef streamxy
  USE mod_streamxy
#endif /*streamxy*/
#ifdef streamv
  USE mod_streamv
#endif /*streamv*/
#ifdef streamr
  USE mod_streamr
#endif /*streamr*/
#ifdef stream_thermohaline
  USE mod_stream_thermohaline
#endif /*stream_thermohaline*/
#ifdef tracer
  USE mod_tracer
#endif /*tracer*/
#ifdef sediment
  USE mod_sed
#endif /*sediment*/
#if defined larval_fish
  USE mod_fish
#endif /*fish*/
  
  IMPLICIT none
    
  INTEGER mra,mta,msa
  
#if defined sediment
  INTEGER nsed,nsusp
  logical res
#endif /*sediment*/
 
  INTEGER                                    :: ia, ja, ka, iam
  INTEGER                                    :: ib, jb, kb, ibm
  INTEGER                                    :: i,  j,  k, l, m
  INTEGER                                    :: niter
  INTEGER                                    :: nrh0=0
  
#if defined larval_fish
  ! Specific for fish code
  REAL*8                                     :: rhof    ! density of fish
  REAL*8                                     :: est, group, egg_sg, age
  INTEGER                                    :: istage
#endif /*larval_fish*/

  ! === Variables to interpolate fields ===

  ! Counters
  INTEGER                                    :: nout=0, nloop=0, nerror=0
  INTEGER                                    :: nnorth=0, ndrake=0, ngyre=0
  INTEGER                                    :: nexit(NEND)
  
  REAL*8                                     :: temp, salt, dens
  REAL*8                                     :: temp2, salt2, dens2
  REAL*8                                     :: x0, y0, z0, x1, y1, z1
  REAL*8                                     :: rlon,rlat
  REAL*8                                     :: dt, t0
  REAL*8                                     :: dtreg
  
  
  ! === Error Evaluation ===
  INTEGER                                    :: errCode
  INTEGER                                    :: landError=0 ,boundError=0
  REAL*8                                     :: zz

#if defined sediment
  ! Specific for sediment code
  INTEGER                                    :: nsed,nsusp
  LOGICAL                                    :: res
#endif /*sediment*/


!!------------------------------------------------------------------------------


  iday0=iday
  imon0=imon
  iyear0=iyear
  ! === print some run stats ===
  print *,'------------------------------------------------------'  
  print *,'Files written in directory                  :  ' ,trim(outDataDir)
  print *,'with file names starting with               :  ' ,trim(outDataFile)
  print *,'Time periods (steps) between two GCM fields : ' ,iter
  
  print 999,intstart,intspin,intrun,intend,nff,isec,idir,nqua,num,voltr,&
       tmin0,tmax0,smin0,smax0,rmin0,rmax0
  
999 format(' intstart :',i7,'   intspin :',i7, &
         /,'   intrun :',i7,'   intend  :',i7, &
         /,'      nff :',i2,' isec :',i2,'  idir :',i2,' nqua=',i2,' num=',i7,&
         /,'    voltr : ',f9.0,&
         /,'    tmin0 : ',f7.2,'  tmax0 : ',f7.2, &
         /,'    smin0 : ',f7.2,'  smax0 : ',f7.2,&
         /,'    rmin0 : ',f7.2,'  rmax0 : ',f7.2)

  ! === initialise to zero ===
  nrh0=0
  nexit=0
  ntractot=0
#ifdef sediment
  nsed=0
  nsusp=0
#endif /*sediment*/
  
  nrj=0
  trj=0.d0

  dstep=1.d0/dble(iter)
  dtmin=dstep*tseas
  
  
  !==========================================================
  !=== Read in the end positions from an previous run     === 
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
  if(rlat == 478.) then
     nrj(ntrac,8)=0    ! South
     i=i+1
  elseif(rlon == 100) then
     nrj(ntrac,8)=0    ! East
     j=j+1
  elseif(temp > tmaxe .and. salt < smine .and. tt-t0>200.) then
     nrj(ntrac,8)=1    ! back to the warm pool
     k=k+1
  else
     nrj(ntrac,8)=0 
     l=l+1   
  endif
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#elif defined ifs
  if(rlat == float(jenn(1))) then
     nrj(ntrac,8)=1    ! Southern boundary
  elseif(rlat == float(jens(2))) then
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
  
  print *,'------------------------------------------------------'
  call fancyTimer('initialize dataset','start')
  ff=dble(nff)
!  tstep=dble(intstep) 
  ints=intstart
  call readfields   ! initial dataset
  ntrac=0
  call fancyTimer('initialize dataset','stop')

  !==========================================================
  !==========================================================
  !=== Start main time loop                               ===
  !==========================================================
  !==========================================================
  intsTimeLoop: do ints=intstart+nff,intstart+intrun,nff
!  intsTimeLoop: do ints=intstart+intstep,intstart+intrun,intstep
     call fancyTimer('reading next datafield','start')
     tt = ints*tseas
     call readfields
     call fancyTimer('reading next datafield','stop')
     
     !=======================================================
     !=== write stream functions and "particle tracer"    ===
     !=======================================================
     if(mod(ints,120) == 0) then 
      call writepsi 
      call writetracer
     endif

    intspinCond: if(nff*ints <= nff*(intstart+intspin)) then
        call fancyTimer('seeding','start')
        call seed (tt,ts)
        call fancyTimer('seeding','stop')
        t0 = tt
        dt = 0.d0
     end if intspinCond

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
        
        ! === Write initial data to in.asc file ===
        ! If t0 = tt (first step) 
        
        if(trj(ntrac,4) == trj(ntrac,7)) then
#ifdef tempsalt
        call interp(nrj(ntrac,1),nrj(ntrac,2),nrj(ntrac,3),&
        trj(ntrac,1),trj(ntrac,2),trj(ntrac,3),temp,salt,dens,1)
#endif
        call writedata(10)
        endif

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
#ifdef larval_fish
           ! Find settling velocity for active gridbox ===
!           rhof = fish(ntrac,i_density)
           call fishvel(temp,dens)
#endif /*larval_fish*/
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
           rg=dmod(ts,1.d0) ! time interpolation constant between 0 and 1
           rr=1.d0-rg
           if(rg.lt.0.d0 .or.rg.gt.1.d0) then
              print *,'rg=',rg
              exit intsTimeLoop
           endif
           
#ifndef regional
           ! === Cyclic world ocean/atmosphere === 
           IF (ib == 1 .AND. x1 >= DBLE (IMT)) THEN
              x1 = x1 - DBLE(IMT)
           END IF
#endif
           
           x0=x1
           y0=y1
           z0=z1
           ia=ib
           iam=ia-1
           if(iam == 0)iam=IMT
           ja=jb
           ka=kb

           call calc_dxyz
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
           ! === calculate the turbulent velocities ===
#ifdef turb
           call turbuflux(ia,ja,ka,rr,dt)
#endif /*turb*/
           ! === calculate the vertical velocity ===
           call vertvel(rr,ia,iam,ja,ka)
#ifdef timeanalyt
           ss0=dble(idint(ts))*tseas/dxyz
           call cross_time(1,ia,ja,ka,x0,dse,dsw,ts,tt,dsmin,dxyz,rr) ! zonal
           call cross_time(2,ia,ja,ka,y0,dsn,dss,ts,tt,dsmin,dxyz,rr) ! merid
           call cross_time(3,ia,ja,ka,z0,dsu,dsd,ts,tt,dsmin,dxyz,rr) ! vert
#else
           call cross(1,ia,ja,ka,x0,dse,dsw,rr) ! zonal
           call cross(2,ia,ja,ka,y0,dsn,dss,rr) ! meridional
           call cross(3,ia,ja,ka,z0,dsu,dsd,rr) ! vertical
#endif /*timeanalyt*/
           ds=dmin1(dse,dsw,dsn,dss,dsu,dsd,dsmin)
     
           !if(ds == UNDEF .or.ds == 0.d0)then 
           call errorCheck('dsCrossError', errCode)
           if (errCode.ne.0) cycle ntracLoop

           call calc_time
           
           ! === calculate the new positions ===
           ! === of the trajectory           ===    
           call pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1)
           
#if defined orc || orca025 || orca12
           ! === north fold cyclic for the ORCA grids ===
           ! === only tested with ORCA025 so far !!!  ===
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
            stop 4967
           endif
#endif
           ! === make sure that trajectory ===
           ! === is inside ib,jb,kb box    ===
#ifdef regional
           if (x1.lt.0.d0 .or. x1.gt.(dble(IMT)-0.5) .or. &
     &         y1.lt.0.d0 .or. y1.gt.(dble(JMT)-0.5)) then
                 nexit(NEND)=nexit(NEND)+1
                 exit niterLoop
           endif
           if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! index for correct cell?
#else
           if(x1.lt.0.d0) x1=x1+dble(IMT)           ! east-west cyclic
           if(x1.gt.dble(IMT)) x1=x1-dble(IMT)      ! east-west cyclic
           if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! index for correct cell?
           if(ib.gt.IMT) ib=ib-IMT                  ! east-west cyclic
#endif
           if(y1.ne.dble(idint(y1))) jb=idint(y1)+1 ! index for correct cell?
           
           call errorCheck('boundError', errCode)
           if (errCode.ne.0) cycle ntracLoop

           call errorCheck('landError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           
           call errorCheck('bottomError', errCode)
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
              if(float(ienw(k)) <= x1 .and. x1 <= float(iene(k)) .and. &
                 float(jens(k)) <= y1 .and. y1 <= float(jenn(k))  ) then
                 nexit(k)=nexit(k)+1
                 exit niterLoop                                
              endif
           enddo LBTLOOP
           
#if defined tempsalt
               call interp (ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
!              call statvd(temp, salt, dens, km, depth, lat)
!               if (temp < tmine .or. temp > tmaxe .or. &
!               &   salt < smine .or. salt > smaxe .or. &
!               &   dens < rmine .or. dens > rmaxe      ) then
                if (temp > tmaxe .and. salt < smine .and.  &
               &   (tt-t0)/tday > 200.      ) then
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
        nrj(ntrac,6)=1
     end do ntracLoop
     
#ifdef sediment
     print 599,ints,ntime,ntractot,nout,nloop,nerror,ntractot-nout, & 
          nsed,nsusp,nexit
599  format('ints=',i7,' time=',i10,' ntractot=',i8,' nout=',i8, & 
          ' nloop=',i4,' nerror=',i4,' in ocean/atm=',i8,' nsed=',i8, & 
          ' nsusp=',i8,' nexit=',9i8)
#elif defined ifs || rco || tes || orc || baltix || orca025 || orca12 || AusCOM
     print 799 ,ntime,ints ,ntractot ,nout ,nerror,ntractot-nout
799  format('ntime=',i10,' ints=',i7,' ntractot=',i8,' nout=',i8, & 
          ' nerror=',i4,' in ocean/atm=',i8)
#else
     call fancyTimer('advection','stop') 
     print 799 ,ints ,ntractot-nout ,nout ,nerror,ntractot 
799  format('ints=',i7,' active=',i10,' out=',i10,' err=',i10,' tot=',i10)
#endif
  
   IF (ntractot /= 0 .AND. ntractot - nout - nerror == 0  .AND.                &
   &   seedTime /= 2) THEN
      EXIT intsTimeLoop
   END IF
  
  end do intsTimeLoop
  
  close(56)
  print *,ntractot ,' trajectories calculated'
  print *,nout     ,' trajectories exited the space and time domain'
  print *,nexit    ,' trajectories exited through the boundaries'
#ifdef sediment
  print *,nsed     ,' trajectories sedimented'
  print *,nsusp    ,' trajectories resuspended'
  call writedata(19)
  
#endif
#ifdef tempsalt     
  print *,nrh0,' trajectories outside density range'
#endif
  print *,nloop,' infinite loops'
  print *,nerror,' error loops'
  print *,ntractot-nout-nrh0-nerror,' trajectories in domain'
  
  call writepsi
  
  print *,'The very end of TRACMASS run ',outDataFile,' at'
  call system('tput bel')
  call system('date')
  
return
     
     



























   CONTAINS
     
     subroutine errorCheck(teststr,errCode)
       CHARACTER (len=*)                   :: teststr    
       INTEGER                             :: verbose = 0
       INTEGER                             :: strict  = 0
       INTEGER                             :: errCode

       errCode=0
       select case (trim(teststr))
       case ('ntracGTntracmax')
          if(ntrac.gt.ntracmax) then
             print *,'====================================='
             print *,'ERROR: to many trajectories,'
             print *,'-------------------------------------'             
             print *,'increase ntracmax since'
             print *,'ntrac >',ntrac
             print *,'when ints=',ints,' and ' 
             print *,'intspin=',intspin
             print *,',(intspin-ints)/ints*ntrac='
             print *,(intspin-ints)/ints*ntrac
             print *,'-------------------------------------'
             print *,'The run is terminated'
             print *,'====================================='
             errCode = -38
             stop
          endif
          
       case ('dxyzError')
          if(dxyz == 0.d0) then
             if (verbose == 1) then                 
                print *,'====================================='
                print *,'ERROR: dxyz is zero'
                print *,'-------------------------------------'
                print *,'ntrac=',ntrac,' ints=', ints
                print *,'ib=',ib,'jb=',jb,'kb=',kb
                print *,'kmt=',kmt(ib,jb)
                print *,'dz=',dz(kb)
                print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb)
!                print *,'dztb=',dztb(ib,jb,1),dztb(ib,jb,2)
!                print *,'rg*hs=',rg,hs(ib,jb,NST)
!                print *,'rr*hs=',rr,hs(ib,jb,1)
                print *,'-------------------------------------'
                print *,'The trajectory is killed'
                print *,'====================================='
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
                print *,'====================================='
                print *,'Warning: Trajectory leaving model area'
                print *,'-------------------------------------'
                print *,'iaib',ia,ib,ja,jb,ka,kb
                print *,'xyz',x0,x1,y0,y1,z0,z1
                print *,'ds',dse,dsw,dsn,dss,dsu,dsd
                print *,'dsmin=',ds,dsmin,dtmin,dxyz
                print *,'tt=',tt,ts
                print *,'ntrac=',ntrac
                print *,'-------------------------------------'
                print *,'The trajectory is killed'
                print *,'====================================='
             end if
             call writedata(19)
             nerror=nerror+1
             boundError = boundError +1
             errCode = -50
             if (strict==1) stop
             call writedata(40)
             nrj(ntrac,6)=1
          endif

       case ('landError')
          if(kmt(ib,jb) == 0) then
             if (verbose == 1) then
                print *,'====================================='
                print *,'Warning: Trajectory on land'
                print *,'-------------------------------------'
                print *,'land',ia,ib,ja,jb,ka,kb,kmt(ia,ja)
                print *,'xyz',x0,x1,y0,y1,z0,z1
                print *,'ds',ds,dse,dsw,dsn,dss,dsu,dsd
                print *,'dsmin=',ds,dsmin,dtmin
                print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb),dxdy(ia,ja)
                print *,'hs=',hs(ia,ja,1),hs(ia,ja,2),hs(ib,jb,1),hs(ib,jb,2)
                print *,'tt=',tt,ts,tt/tday,t0/tday
                print *,'ntrac=',ntrac
                print *,'niter=',niter
#ifdef turb
                print *,'upr=',upr
#endif
                print *,'-------------------------------------'
                print *,'The trajectory is killed'
                print *,'====================================='
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
             print *,'========================================'
             print *,'ERROR: Particle overshoot in i direction'
             print *,'----------------------------------------'
             print *,ib-1,x1,ib,ntrac,ib,jb,kb
             x1=dble(ib-1)+0.5d0
             ib=idint(x1)+1
             print *,'error i',ib-1,x1,ib,ntrac,ib,jb,kb
             print *,y1,z1
             print *,'-------------------------------------'
             print *,'The run is terminated'
             print *,'====================================='             
             errCode = -42
             stop
          elseif( dble(jb-1).gt.y1 .or. dble(jb).lt.y1 )  then
             print *,'========================================'
             print *,'ERROR: Particle overshoot in j direction'
             print *,'----------------------------------------'
             print *,'error j',jb-1,y1,jb,ntrac,x1,z1
             print *,'error j',jb-1,y1,jb,ntrac,ib,jb,kb
             print *,'-------------------------------------'
             print *,'The run is terminated'
             print *,'====================================='    
             errCode = -44
             stop
          elseif((dble(kb-1).gt.z1.and.kb.ne.KM).or. & 
               dble(kb).lt.z1 ) then
             print *,'========================================'
             print *,'ERROR: Particle overshoot in k direction'
             print *,'----------------------------------------'
             print *,'error k',kb-1,z1,kb,ntrac,x1,y1
             print *,'error k',kb-1,z1,kb,ntrac,ib,jb,kb
             print *,'-------------------------------------'
             print *,'The run is terminated'
             print *,'====================================='
             errCode = -46
             stop
          end if
       case ('infLoopError')
          if(niter-nrj(ntrac,4).gt.30000) then ! break infinite loops
             nloop=nloop+1             
!             nerror=nerror+1
             if (verbose == 1) then
                print *,'====================================='
                print *,'Warning: Particle in infinite loop '
                print *,'ntrac:',ntrac
                print *,'niter:',niter,'nrj:',nrj(ntrac,4)
                print *,'dxdy:',dxdy(ib,jb),'dxyz:',dxyz
                print *,'kmt:',kmt(ia-1,ja-1),'dz(k):',dz(ka-1)
                print *,'ia=',ia,' ib=',ib,' ja=',ja,' jb=',jb, & 
                     ' ka=',ka,' kb=',kb
                print *,'x1=',x1,' x0=',x0,' y1=',y1,' y0=',y0, & 
                     ' z1=',z1,' z0=',z0
                print *,'u(ia )=',(rbg*uflux(ia ,ja,ka,NST) + &
                     rb*uflux(ia ,ja,ka,1))*ff
                print *,'u(iam)=',(rbg*uflux(iam,ja,ka,NST) + & 
                     rb*uflux(iam,ja,ka,1))*ff
                print *,'v(ja  )=',(rbg*vflux(ia,ja  ,ka,NST) + & 
                     rb*vflux(ia,ja  ,ka,1))*ff
                print *,'v(ja-1)=',(rbg*vflux(ia,ja-1,ka,NST) + & 
                     rb*vflux(ia,ja-1,ka,1))*ff
                print *,'-------------------------------------'
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
             nrj(ntrac,6)=0  ! 0=continue trajectory, 1=end trajectory
             nrj(ntrac,7)=1
             errCode = -48
          end if
       case ('bottomError')
          ! if trajectory under bottom of ocean, 
          ! then put in middle of deepest layer 
          ! (this should however be impossible)
           if( z1.lt.dble(KM-kmt(ib,jb)) ) then
              print *,'under bottom !!!!!!!',z1,dble(KM-kmt(ib,jb)), &
                   kmt(ia,ja),kmt(ib,jb),ntrac
               print *,'ds',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin,dxyz
               print *,'ia=',ia,ib,ja,jb,ka,kb
               print *,'x0=',x0,x1,y0,y1,z0,z1
               call cross(1,ia,ja,ka,x0,dse,dsw,rr) ! zonal
               call cross(2,ia,ja,ka,y0,dsn,dss,rr) ! meridional
               call cross(3,ia,ja,ka,z0,dsu,dsd,rr) ! vertical
               print *,'time step sol:',dse,dsw,dsn,dss,dsu,dsd
               print *,'ntrac=',ntrac
              nerror=nerror+1
              nrj(ntrac,6)=1
               stop 3957
               z1=dble(KM-kmt(ib,jb))+0.5d0
              errCode = -49
           end if
        case ('airborneError')
           ! if trajectory above sea level,
           ! then put back in the middle of shallowest layer (evaporation)
           if( z1.gt.dble(KM) ) then
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
                 print *,'==================================================='
                 print *,'Warning: not find any path for unknown reason '
                 print *, " "
                 write (*,'(A, E9.3, A, E9.3)'), ' uflux= ', &
                      uflux(ia,ja,ka,1),'  vflux= ', vflux(ia,ja,ka,1)

                 write (*,FMT='(A, 5E9.2)'),' ds=',ds,dse,dsw,dsn,dss
                 write (*,FMT='(4E9.2)'), dsu,dsd,dsmin,dxyz
                 print *,'---------------------------------------------------'
                 print *,"   ntrac = ",ntrac
                 write (*,'(A7, I10, A7, I10, A7, I10)'), & 
                      ' ia= ', ia, ' ja= ', ja, ' ka= ', ka
                 write (*,'(A7, I10, A7, I10, A7, I10)'), & 
                      ' ib= ', ib, ' jb= ', jb, ' kb= ', kb
                 write (*,'(A7, F10.3, A7, F10.3, A7, F10.3)'), & 
                      ' x0= ', x0, ' y0= ', y0, ' z0= ', z0
                 write (*,'(A7, F10.3, A7, F10.3, A7, F10.3)'), & 
                      ' x0= ', x0, ' y0= ', y0, ' z0= ', z0
                 write (*,'(A7, I10, A7, I10, A7, I10)'), & 
                      ' k_inv= ', KM+1-kmt(ia,ja), ' kmt= ', kmt(ia,ja), &
                      'lnd= ', mask(ia,ja)
                 print *,'---------------------------------------------------'
                print *,'The trajectory is killed'
                print *,'==================================================='
              end if
              nerror=nerror+1
              nrj(ntrac,6)=1
              errCode = -56
           end if
        end select
      end subroutine errorCheck

  subroutine writedata(sel)
    REAL*8                               :: vort
    INTEGER                              :: sel ,xf ,yf ,zf ,n
    INTEGER, SAVE                        :: recPosIn=0  ,recPosOut=0
    INTEGER, SAVE                        :: recPosRun=0 ,recPosErr=0
    INTEGER, SAVE                        :: recPosKll=0
    REAL*8                               :: x14 ,y14 ,z14
    REAL*8                               :: twrite

#if defined for || sim 
566 format(i8,i7,f7.2,f7.2,f7.1,f10.2,f10.2 &
         ,f10.1,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined rco || defined baltix
566 format(i8,i7,f7.2,f7.2,f7.1,2f12.4 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined tes 
566 format(i8,i7,f8.3,f8.3,f7.3,2f10.2 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined ifs 
566 format(i8,i7,f7.2,f7.2,f7.2,f10.2,f10.2 &
         ,f15.0,f8.2,f8.2,f8.2,f6.0,8e8.1 )
#elif defined orc
!566 format(i8,i7,2f8.2,f6.2,2f10.2 &
!         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#else
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
!566 format(i7,i7,f7.2,f7.2,f7.1,f10.4,f10.4 &
!         ,f13.4,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#endif
    
    xf   = floor(x1)
    yf   = floor(y1)
    zf   = floor(z1)
    
    if ((sel .ne. 19) .and. (sel.ne.40)) then
! this requires too much memory
!       vort = (vvel(xf+1,yf,zf)-vvel(xf-1,yf,zf))/4000 - &
!            (uvel(xf,yf+1,zf)-uvel(xf,yf-1,zf))/4000   
    end if

#if defined textwrite 
    select case (sel)
    case (10)
       write(58,566) ntrac,niter,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens
    case (11)
       if(  (kriva == 1 .AND. nrj(ntrac,4) == niter-1   ) .or. &
            (kriva == 2 .AND. scrivi                    ) .or. &
            (kriva == 3                                 ) .or. &
            (kriva == 4 .AND. niter == 1                ) .or. &
            (kriva == 5 .AND.                                  &
          &  MOD((DBLE(tt)-DBLE(t0))*DBLE(NGCM)/DBLE(ITER), 3600.) == 0.d0 ) .or. &
            (kriva == 6 .AND. .not.scrivi               )        ) then
#if defined tempsalt
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
#if defined biol
          write(56,566) ntrac,ints,x1,y1,z1,tt/3600.,t0/3600.
#else
          write(56,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens
#endif        
       endif
    case (13)
       ! === write sed pos ===
       write(57,566) ntrac,niter,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens 
    case (14)
       write(56,566) ntrac,ints,x1,y1,z1, &
            tt/60.,t0/3600.,subvol,temp,salt,dens
    case (15)
       write(57,566) ntrac,ints,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens
    case (16)
       if(kriva.ne.0 ) then
#if defined tempsalt
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
          write(56,566) ntrac,ints,x1,y1,z1, &
               tt/tday,t0/tday,subvol,temp,salt,dens
       end if
    case (17)
       write(57,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (19)
       ! === write last sedimentation positions ===
       open(34,file=trim(outDataDir)//trim(outDataFile)//'_sed.asc') 
       do n=1,ntracmax
        if(nrj(n,1).ne.0) then
         write(34,566) n,nrj(n,4),trj(n,1),trj(n,2),trj(n,3),trj(n,4)/tday,trj(n,7)/tday
      endif
       enddo
       close(34)
    case (40)
       write(59,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  

    end select
#endif    
#if defined binwrite 

    x14=real(x1,kind=4)
    y14=real(y1,kind=4)
    z14=real(z1,kind=4)
    if (twritetype==1) then
       twrite = tt
    else if (twritetype==2) then
       call updateclock
       twrite = currJDtot
    else
       twrite = real(ints,kind=8)
    end if
    select case (sel)       
    case (10)
       recPosIn = recPosIn+1
       write(unit=78 ,rec=recPosIn) ntrac,ints,x14,y14,z14
       return
    case (11)
       if(  (kriva == 1 .and. nrj(ntrac,4)  ==  niter-1 ) .or. &
            (kriva == 2 .and. scrivi                    ) .or. &
            (kriva == 3                                 ) .or. &
            (kriva == 4 .and. niter == 1                ) .or. &
            (kriva == 5 .and. abs(dmod(tt-t0,9.d0)) < 1e-5 ) .or. &
            (kriva == 6 .and. .not.scrivi               )  ) then
#if defined tempsalt
          call interp(ib,jb,kb,x1,y1,z1,temp, salt,  dens,1)
          call interp(ib,jb,kb,x1,y1,z1,temp2,salt2, dens2,2)
          !z14=real(salt*rb+salt2*(1-rb),kind=4)
#endif
          recPosRun = recPosRun+1
          write(unit=76 ,rec=recPosRun) ntrac,twrite,x14,y14,z14
       end if
    case (13)
       recPosKll = recPosKll+1
       write(unit=77 ,rec=recPosKll) ntrac,twrite,x14,y14,z14   
    case (15)
       recPosRun = recPosRun+1
       write(unit=76 ,rec=recPosRun) ntrac,twrite,x14,y14,z14   
    case (17)
       recPosOut = recPosOut+1
       write(unit=77 ,rec=recPosOut) ntrac,twrite,x14,y14,z14   
    case (19)
       recPosOut = recPosOut+1
       write(unit=75 ,rec=recPosOut) ntrac,twrite,x14,y14,z14
    case (40)
       recPosErr=recPosErr+1    
       write(unit=79 ,rec=recPosErr) ntrac,twrite,x14,y14,z14   
    end select
#endif    

  end subroutine writedata

  subroutine calc_dxyz
    ! T-box volume in m3
#ifdef zgrid3Dt 
    dxyz=rg*dzt(ib,jb,kb,NST)+rr*dzt(ib,jb,kb,1)
#elif  zgrid3D
    dxyz=dzt(ib,jb,kb)
#ifdef freesurface
    if(kb == KM) dxyz=dxyz+rg*hs(ib,jb,NST)+rr*hs(ib,jb,1)
#endif /*freesurface*/
#else
    dxyz=dz(kb)
#ifdef varbottombox
    if(kb == KM+1-kmt(ib,jb) ) dxyz=dztb(ib,jb,1)
#endif /*varbottombox*/
#ifdef freesurface
    if(kb == KM) dxyz=dxyz+rg*hs(ib,jb,NST)+rr*hs(ib,jb,1)
#endif /*freesurface*/
#endif /*zgrid3Dt*/
    dxyz=dxyz*dxdy(ib,jb)
    if (dxyz<0) then
       print *,'=========================================================='
       print *,'ERROR: Negative box volume                                '
       print *,'----------------------------------------------------------'
       print *,'dzt  = ', dxyz/dxdy(ib,jb),dz(kb),hs(ib,jb,:)
       print *,'dxdy = ', dxdy(ib,jb)
       print *,'ib  = ', ib, ' jb  = ', jb, ' kb  = ', kb 
       print *,'----------------------------------------------------------'
       print *,'The run is terminated'
       print *,'=========================================================='
       errCode = -60
       !stop
    end if
  end subroutine calc_dxyz

  subroutine calc_time
#ifdef regulardt
           if(ds == dsmin) then ! transform ds to dt in seconds
!            dt=dt  ! this makes dt more accurate
           else
            dt=ds*dxyz 
           endif
#else
           if(ds == dsmin) then ! transform ds to dt in seconds
              dt=dtmin  ! this makes dt more accurate
           else
              dt=ds*dxyz 
           endif
#endif /*regulardt*/
           if(dt.lt.0.d0) then
              print *,'dt=',dt
              stop 49673
           endif
           ! === if time step makes the integration ===
           ! === exceed the time when fiedls change ===
           if(tss+dt/tseas*dble(iter).ge.dble(iter)) then
              dt=dble(idint(ts)+1)*tseas-tt
              tt=dble(idint(ts)+1)*tseas
              ts=dble(idint(ts)+1)
              tss=dble(iter)
              ds=dt/dxyz
              dsc=ds
           else
              tt=tt+dt
#if defined regulardt
              if(dt == dtmin) then
                 ts=ts+dstep
                 tss=tss+1.d0
              elseif(dt == dtreg) then  
                 ts=nint((ts+dtreg/tseas)*dble(iter))/dble(iter)
!                 ts=ts+dtreg/tseas
                 tss=dble(nint(tss+dt/dtmin))
              else
                 ts=ts+dt/tseas
                 tss=tss+dt/dtmin
              endif
#else
              if(dt == dtmin) then
                 ts=ts+dstep
                 tss=tss+1.d0
              else
                 ts =ts +dt/tseas
                 tss=tss+dt/tseas*dble(iter)
!                 tss=tss+dt/dtmin
              endif
#endif
           end if
           ! === time interpolation constant ===
           rbg=dmod(ts,1.d0) 
           rb =1.d0-rbg
         end subroutine calc_time

  subroutine fancyTimer(timerText ,testStr)
    IMPLICIT NONE

    CHARACTER (len=*)                          :: timerText ,testStr
    REAL*8 ,SAVE                               :: fullstamp1 ,fullstamp2
    REAL*8 ,SAVE ,DIMENSION(2)                 :: timestamp1 ,timestamp2
    REAL*8                                     :: timeDiff
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

