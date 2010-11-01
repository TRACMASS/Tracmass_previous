subroutine loop
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
#ifdef tracer
  USE mod_tracer
#endif /*tracer*/
#ifdef sediment
  USE mod_sed
#endif /*sediment*/
  
  IMPLICIT none

#ifdef streamxy
  REAL sxyy(IMT,JMT),sxyx(IMT,JMT)
#endif /*streamxy*/
#ifdef streamv
  REAL sxz(JMT,KM),syz(JMT,KM)
#endif /*streamv*/
#ifdef streamr
  REAL sxr(IMT,MR,LOV),syr(JMT,MR,LOV)
#endif /*streamr*/
    
  INTEGER mra,mta,msa
  REAL temp,salt,dens
  
#if defined sediment
  INTEGER nsed,nsusp
  logical res
#endif /*sediment*/
  
  INTEGER                                    :: ia, ja, ka, iam
  INTEGER                                    :: ib, jb, kb, ibm
  INTEGER                                    :: i,  j,  k, l, m
  INTEGER                                    :: niter
  INTEGER                                    :: ntrac,nrh0=0

  ! Counters
  INTEGER                                    :: nout=0, nloop=0, nerror=0
  INTEGER                                    :: nnorth=0, ndrake=0, ngyre=0
  INTEGER                                    :: nexit(NEND)
  
  REAL*8                                     :: x0, y0, z0, x1, y1, z1
  REAL*8                                     :: rlon,rlat
  REAL*8                                     :: dt, t0
  REAL*8                                     :: dtreg
  
  
  ! === Error Evaluation ===
  INTEGER                                    :: errCode
  INTEGER                                    :: landError=0 ,boundError=0
  REAL                                       :: zz

!_________________________________________________________________  
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
!  nev=0
  nrh0=0
  nout=0
  nloop=0
  ndrake=0
  nexit=0
  nnorth=0
  ntractot=0
  ngyre=0
  nerror=0
#ifdef sediment
  nsed=0
  nsusp=0
#endif /*sediment*/
  
  nrj=0
  trj=0.d0
  
#ifdef streamxy
  sxyy=0.
  sxyx=0.
#endif /*streamxy*/
#ifdef streamv
  stxz=0.
  styz=0.
#endif /*streamv*/
#ifdef streamr
  stxr=0.
  styr=0.
#endif /*streamr*/


  dstep=1.d0/dble(iter)
  dtmin=dstep*tseas
!  dtmin=tseas/dble(iter)
  
  !==========================================================
  !=== Read in the end positions from an previous run     === 
  !==========================================================
  
#ifdef rerun
  print *,'rerun with initial points from ', & 
       trim(outDataDir)//trim(outDataFile)//'_rerun.asc'
  open(67,file=trim(outDataDir)//trim(outDataFile)//'_rerun.asc')
40 continue
  read(67,566,end=41,err=41) ntrac,niter,rlon,rlat,zz
!  print 566,ntrac,niter,rlon,rlat,zz

#ifdef orc
!  do k=1,LBT
!     if(ienw(k).le.rlon .and. rlon.le.iene(k)) then
!        nrj(ntrac,8)=k                               
!     endif
!  enddo
  nrj(ntrac,8)=1                               
  if(nrj(ntrac,8).eq.0) stop 7395                               
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#elif defined  occ66

  if(rlon.eq.293.) then
     nrj(ntrac,8)=3    ! Drake Passage
  elseif(rlat.eq.-47.00) then
     nrj(ntrac,8)=2    ! Mid Gyre  (47S)
  elseif(rlat.eq.-29.25) then
     nrj(ntrac,8)=1    ! Northern boundary (30S)
  else
     nrj(ntrac,8)=0
     print 566,ntrac,n,rlon,rlat,zz
     stop 4957
  endif

#elif defined ifs
  if(rlat.eq.float(jenn(1))) then
     nrj(ntrac,8)=1    ! Southern boundary
  elseif(rlat.eq.float(jens(2))) then
     nrj(ntrac,8)=2    ! Northern boundary 
  else
     nrj(ntrac,8)=0
     print 566,ntrac,niter,rlon,rlat,zz
     stop 4957
  endif /*orc*/
566 format(i8,i7,2f8.2,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#endif
!  print 566,ntrac,niter,rlon,rlat,zz
!  print *,nrj(ntrac,8)  
  goto 40
41 continue
  print 566,ntrac,niter,rlon,rlat,zz
  do ntrac=1,ntracmax
     if(nrj(ntrac,8).eq.0) nrj(ntrac,6)=1 
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
  tstep=dble(intstep) 
  ints=intstart
  call readfields   ! initial dataset
  ntrac=0
  call fancyTimer('initialize dataset','stop')

  !==========================================================
  !==========================================================
  !=== Start main time loop                               ===
  !==========================================================
  !==========================================================
  
  intsTimeLoop: do ints=intstart+intstep,intstart+intrun,intstep
     call fancyTimer('reading next datafield','start')
     call readfields
     call fancyTimer('reading next datafield','stop')
     
     if(mod(ints,120).eq.0 .and. ints.ne.0) call writepsi ! write psi
     if(mod(ints,120).eq.0) call writetracer

    intspinCond: if(nff*ints <= nff*(intstart+intspin)) then
        call fancyTimer('seeding','start')
        call seed (tt,ts)
        call fancyTimer('seeding','stop')
        t0 = tt
        dt = 0.d0
        arct = 0.d0
     end if intspinCond

     if(ntractot-nout-nerror.eq.0) exit intsTimeLoop
          
     !=======================================================
     !=== Loop over all trajectories and calculate        ===
     !=== a new position for this time step.              ===
     !=======================================================
     call fancyTimer('advection','start')
     ntracLoop: do ntrac=1,ntractot  
        ! === Test if the trajectory is dead   ===
        if(nrj(ntrac,6).eq.1) cycle ntracLoop
        
        ! === Read in the position, etc at the === 
        ! === beginning of new time step       ===
        x1=trj(ntrac,1)
        y1=trj(ntrac,2)
        z1=trj(ntrac,3)
        tt=trj(ntrac,4)
        subvol=trj(ntrac,5)
        arct=trj(ntrac,6)
        t0=trj(ntrac,7)
        ib=nrj(ntrac,1)
        jb=nrj(ntrac,2)
        kb=nrj(ntrac,3)
        niter=nrj(ntrac,4)
        ts=dble(nrj(ntrac,5))
        tss=0.d0
        
        ! === Write initial data to in.asc file ===
        ! If t0 = tt (first step) 
        if(trj(ntrac,4).eq.trj(ntrac,7)) then
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
        if( nrj(ntrac,6).eq.2 ) then
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
           if( niter.ne.1 .and. tss.eq.dble(iter) &
                .and. nrj(ntrac,7).ne.1 ) then
              trj(ntrac,1)=x1
              trj(ntrac,2)=y1
              trj(ntrac,3)=z1
              trj(ntrac,4)=tt
              trj(ntrac,5)=subvol
              trj(ntrac,6)=arct
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
           ! === Cyclic world ocean/atmosphere === 
           if(ib.eq.1.and.x1.eq.dble(IMT)) x1=0.d0
           x0=x1
           y0=y1
           z0=z1
           ia=ib
           iam=ia-1
           if(iam.eq.0)iam=IMT
           ja=jb
           ka=kb

           call calc_dxyz
           call errorCheck('dxyzError'     ,errCode)
           call errorCheck('coordBoxError' ,errCode)
           call errorCheck('infLoopError'  ,errCode)
           if (errCode.ne.0) cycle ntracLoop

           ! === calculate the turbulent velocities ===
#ifdef turb
           call turbuflux(ia,ja,ka,rr)
#endif /*turb*/
           ! === calculate the vertical velocity ===
           call vertvel(rr,ia,iam,ja,ka)
           
           ! === write trajectory ===                       
#ifdef tracer
           if(ts.eq.dble(idint(ts))) then 
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
          
           !if(ds.eq.UNDEF .or.ds.eq.0.d0)then 
           call errorCheck('dsCrossError', errCode)
           if (errCode.ne.0) cycle ntracLoop

           call calc_time

           ! === calculate the new positions ===
           ! === of the trajectory           ===    
           call pos(ia,iam,ja,ka,ib,jb,kb,x0,y0,z0,x1,y1,z1)

#ifdef orc 
! north fold cyclic ORCA grid (only tested with ORCA025)
           if(y1.eq.dble(JMT)) then 
            x1=dble(IMT+4)-x1    
            if(x1.ne.dble(idint(x1))) then
             ib=idint(x1)+1
            else
             ib=idint(x1)
            endif
            jb=JMT-1
            y1=dble(jb-1)
            ! make sure the f-fold is also made for the previous point:
            x0=x1 ; y0=y1 ; ia=ib ; ja=jb  
           endif
#endif           

           ! === make sure that trajectory ===
           ! === is inside ib,jb,kb box    ===
           if(x1.lt.0.d0) x1=x1+dble(IMT)           ! east-west cyclic
           if(x1.gt.dble(IMT)) x1=x1-dble(IMT)      ! east-west cyclic
           if(x1.ne.dble(idint(x1))) ib=idint(x1)+1 ! index for correct cell?
           if(ib.gt.IMT) ib=ib-IMT                  ! east-west cyclic
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
           ! === Calculate arclength of the ===
           ! === trajectory path in the box ===
           call arclength(ia,ja,ka,dt,rr,arc)
           arct=arct+arc*arcscale 
           ! === end trajectory if outside chosen domain ===
#if defined occ
           ! === stop and select stream function ===
           northbound: if( y1.eq.dble(jmt-2) .and. n.ne.1 ) then
              nnorth=nnorth+1
              nexit(1)=nexit(1)+1
           ! or to Drake:
           elseif(x1.eq.dble(idr) .and. ja.lt.94 .and. n.ne.1 ) then 
              ndrake=ndrake+1
              call writedata(15)
              nexit(2)=nexit(2)+1
              cycle niterLoop
              ! or continue trajectory
           else
              cycle niterLoop                                   
           endif northbound
#else     
           LBTloop: do k=1,LBT
!              if(ienw(k).le.ib .and. ib.le.iene(k) .and. &
!                 jens(k).le.jb .and. jb.le.jenn(k)  ) then
              if(float(ienw(k)) <= x1 .and. x1 <= float(iene(k)) .and. &
                 float(jens(k)) <= y1 .and. y1 <= float(jenn(k))  ) then
                 nexit(k)=nexit(k)+1

                 exit niterLoop                                
              endif
           enddo LBTLOOP
           ! === stop trajectory if the choosen time or ===
           ! === water mass properties are exceeded     ===
           if(tt-t0.gt.timax) then
              nexit(NEND)=nexit(NEND)+1
              exit niterLoop
           endif
           call writedata(18)
        end do niterLoop
#endif
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
#elif defined ifs || rco || tes || orc
     print 799 ,ntime,ints ,ntractot ,nout ,nerror,ntractot-nout
799  format('ntime=',i10,' ints=',i7,' ntractot=',i8,' nout=',i8, & 
          ' nerror=',i4,' in ocean/atm=',i8)
#else
     call fancyTimer('advection','stop') 
     print 799 ,ints ,ntractot-nout ,nout ,nerror,ntractot !,nev
799  format('ints=',i7,' active=',i10,' out=',i10,' err=',i10,' tot=',i10)
#endif

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
  call system('date')
  
return
     
     





















































   CONTAINS
     
     subroutine errorCheck(teststr,errCode)
       CHARACTER (len=*)                   :: teststr    
       INTEGER                             :: verboseMess = 0
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
          if(dxyz.eq.0.d0) then
             if (verboseMess == 1) then                 
                print *,'====================================='
                print *,'ERROR: dxyz is zero'
                print *,'-------------------------------------'
                print *,'ntrac=',ntrac,' ints=', ints
                print *,'ib=',ib,'jb=',jb,'kb=',kb
                print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb)
                !print *,'dztb=',dztb(ib,jb,kb,1),dztb(ib,jb,kb,2)
                print *,'rg*hs=',rg*hs(ib,jb,NST)
                print *,'rr*hs=',rr*hs(ib,jb,1)
                print *,'-------------------------------------'
                !             print *,'The run is terminated'
                print *,'The trajectory is killed'
                print *,'====================================='
             end if
             errCode = -39
             call writedata(40)
             nrj(ntrac,6)=1
          endif          

       case ('boundError')
          if(ia>imt .or. ib>imt .or. ja>jmt .or. jb>jmt &
               .or. ia<1 .or. ib<1 .or. ja<1 .or. jb<1) then
             if (verboseMess == 1) then
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
             boundError = boundError +1
             errCode = -50
             stop 40962
             call writedata(40)
             nrj(ntrac,6)=1
          endif

       case ('landError')
          if(kmt(ib,jb).eq.0) then
             if (verboseMess == 1) then
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
             if (verboseMess == 1) then
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
             trj(ntrac,6)=arct
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
           if( z1.le.dble(KM-kmt(ib,jb)) ) then
              print *,'under bottom !!!!!!!',z1,dble(KM-kmt(ib,jb)), &
                   kmt(ia,ja),kmt(ib,jb),ntrac
               print *,'ds',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin,dxyz
               print *,'ia=',ia,ib,ja,jb,ka,kb
               print *,'x0=',x0,x1,y0,y1,z0,z1
               call cross(1,ia,ja,ka,x0,dse,dsw,rr) ! zonal
               call cross(2,ia,ja,ka,y0,dsn,dss,rr) ! meridional
               call cross(3,ia,ja,ka,z0,dsu,dsd,rr) ! vertical
               print *,'time step sol:',dse,dsw,dsn,dss,dsu,dsd
              nerror=nerror+1
              nrj(ntrac,6)=1
!               stop 3957
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
              if(kb.eq.KM+1) kb=KM  ! (should perhaps be removed)
              errCode = -52
           endif
           case ('cornerError')
              ! problems if trajectory is in the exact location of a corner
           if(x1.eq.dble(idint(x1)) .and. y1.eq.dble(idint(y1))) then
              !print *,'corner problem',ntrac,x1,x0,y1,y0,ib,jb
              !print *,'ds=',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
              !stop 34957
              ! corner problems may be solved the following way 
              ! but should really not happen at all
              if(ds.eq.dse .or. ds.eq.dsw) then
                 if(y1.ne.y0) then
                    y1=y0 ; jb=ja
                 else
                    y1=dble(jb)-0.5d0
                 endif
              elseif(ds.eq.dsn .or. ds.eq.dss) then
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
           if(ds.eq.UNDEF .or.ds.eq.0.d0)then 
              if (verboseMess == 0) then
                 print *,'ds cross error',ds,dse,dsw,dsn,dss,dsu
                 print *,dsd,dsmin,dxyz
                 print *,ia,ja,ka,x0,y0,z0,ntrac,niter
                 print *,'k=',ka,kb,KM+1-kmt(ia,ja),kmt(ia,ja)
              end if
              nerror=nerror+1
              nrj(ntrac,6)=1
              errCode = -56
           end if
        end select
      end subroutine errorCheck

  subroutine writedata(sel)
    REAL                                 :: vort
    INTEGER                              :: sel ,xf ,yf ,zf ,n
    INTEGER, SAVE                        :: recPosIn=0  ,recPosOut=0
    INTEGER, SAVE                        :: recPosRun=0 ,recPosErr=0
    INTEGER, SAVE                        :: recPosKll=0
    REAL                                 :: x14 ,y14 ,z14

#if defined for || sim 
566 format(i8,i7,f7.2,f7.2,f7.1,f10.2,f10.2 &
         ,f10.1,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined rco 
566 format(i8,i7,f7.2,f7.2,f7.1,2f10.2 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined tes 
566 format(i8,i7,f8.3,f8.3,f7.3,2f10.2 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined ifs 
566 format(i8,i7,f7.2,f7.2,f7.2,f10.2,f10.2 &
         ,f15.0,f6.1,f6.2,f8.2,f6.0,8e8.1 )
#elif defined orc
!566 format(i8,i7,2f8.2,f6.2,2f10.2 &
!         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#else
566 format(i7,i7,f7.2,f7.2,f7.1,f10.4,f10.4 &
         ,f13.8,f6.2,f6.2,f6.2,f6.0,8e8.1 )
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
       if( (kriva.eq.1 .and. nrj(ntrac,4) .eq. niter-1 ) .or. &
            (scrivi .and. kriva.eq.2)                .or. &
            (kriva.eq.3)                             .or. &
            (kriva.eq.4 .and. niter.eq.1)            .or. &
            (kriva.eq.5 .and. &
            (tt-t0.eq.7.d0*tday.or.tt-t0.eq.14.d0*tday.or.tt-t0.eq.21.d0*tday)) .or. &
            (.not.scrivi .and. kriva.eq.6)           .or. &
            (dmod(tt-t0,1.d0).eq.0.d0 .and. kriva.eq.7)                ) then
#if defined tempsalt
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
#if defined biol
          write(56,566) ntrac,ints,x1,y1,z1,tt/3600.,t0/3600.
#else
          write(56,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens,arct
#endif        
       endif
    case (13)
       ! === write sed pos ===
       write(57,566) ntrac,niter,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens 
    case (14)
       write(56,566) ntrac,ints,x1,y1,z1, &
            tt/60.,t0/3600.,subvol,temp,salt,dens,arct
    case (15)
       write(57,566) ntrac,ints,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens
    case (16)
       if(kriva.ne.0 ) then
#if defined tempsalt
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
          write(56,566) ntrac,ints,x1,y1,z1, &
               tt/tday,t0/tday,subvol,temp,salt,dens,arct
       end if
    case (17)
       write(57,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (18)
       if( kriva.ne.0 .and. ts.eq.dble(idint(ts)) .and. &
            ints.eq.intstart+intrun) then 
#if defined tempsalt
          call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
#endif
          !write(56,566) ntrac,ints,x1,y1,z1, &
          !tt/tday,t0/tday,subvol,temp,salt,dens,arct
!          write(56,566) ntrac,ints,x1,y1,z1, & 
!               uvel(xf,yf,zf),vvel(xf,yf,zf),vort,temp,salt,dens,arct
          ! write(56,566) ntrac,niter,x1,y1,z1,tt/3600.,t0/3600.
          !,subvol,temp,salt,dens,arct
       endif
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
    select case (sel)       
    case (10)
       recPosIn = recPosIn+1
       write(unit=78 ,rec=recPosIn) ntrac,ints,x14,y14,z14
       return
    case (11)
       if( (kriva.eq.1 .and. nrj(ntrac,4) .eq. niter-1 ) .or. &
            (scrivi .and. kriva.eq.2)                .or. &
            (kriva.eq.3)                             .or. &
            (kriva.eq.4 .and. niter.eq.1)            .or. &
            (kriva.eq.5 .and. &
            (tt-t0.eq.7.*tday.or.tt-t0.eq.14.*tday & 
            .or.tt-t0.eq.21.*tday)) ) then
          call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
!         call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)          
          recPosRun = recPosRun+1
          write(unit=76 ,rec=recPosRun) ntrac,ints,x14,y14,z14
       end if
    case (13)
       recPosKll = recPosKll+1
       write(unit=77 ,rec=recPosKll) ntrac,ints,x14,y14,z14   
    case (15)
       recPosRun = recPosRun+1
       write(unit=76 ,rec=recPosRun) ntrac,ints,x14,y14,z14   
    case (17)
       recPosOut = recPosOut+1
       write(unit=77 ,rec=recPosOut) ntrac,ints,x14,y14,z14   
    case (18)
    !   if( kriva.ne.0 .and. ts.eq.dble(idint(ts)) .and. &
    !        ints.eq.intstart+intrun) then 
    !      call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
    !      recPosRun = recPosRun+1
    !      write(unit=76 ,rec=recPosRun) ntrac,ints,x14,y14,z14   
    !   endif
       !   !case (19)
    case (19)
       recPosOut = recPosOut+1
       write(unit=75 ,rec=recPosOut) ntrac,ints,x14,y14,z14
    case (40)
       recPosErr=recPosErr+1    
       write(unit=79 ,rec=recPosErr) ntrac,ints,x14,y14,z14   
    end select
#endif    

  end subroutine writedata

  subroutine calc_dxyz
    ! T-box volume in m3
#ifdef zgrid3Dt 
    dxyz=rg*dzt(ib,jb,kb,NST)+rr*dzt(ib,jb,kb,1)
#elif  zgrid3D
    dxyz=dzt(ib,jb,kb)
#else
    dxyz=dz(kb)
#endif /*zgrid3Dt*/
#ifdef varbottombox
    if(kb.eq.KM+1-kmt(ib,jb) ) dxyz=dztb(ib,jb,1)
#endif /*varbottombox*/
#ifdef freesurface
    if(kb.eq.KM) dxyz=dxyz+rg*hs(ib,jb,NST)+rr*hs(ib,jb,1)
#endif /*freesurface*/
    dxyz=dxyz*dxdy(ib,jb)
  end subroutine calc_dxyz

  subroutine calc_time
#ifdef regulardt
           if(ds.eq.dsmin) then ! transform ds to dt in seconds
!            dt=dt  ! this makes dt more accurate
           else
            dt=ds*dxyz 
           endif
#else
           if(ds.eq.dsmin) then ! transform ds to dt in seconds
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
              if(dt.eq.dtmin) then
                 ts=ts+dstep
                 tss=tss+1.d0
              elseif(dt.eq.dtreg) then  
                 ts=nint((ts+dtreg/tseas)*dble(iter))/dble(iter)
!                 ts=ts+dtreg/tseas
                 tss=dble(nint(tss+dt/dtmin))
              else
                 ts=ts+dt/tseas
                 tss=tss+dt/dtmin
              endif
#else
              if(dt.eq.dtmin) then
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

