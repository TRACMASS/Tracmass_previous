subroutine loop
  USE mod_param
  USE mod_name
  USE mod_time
  USE mod_grid
  USE mod_buoyancy
  USE mod_seed
  USE mod_domain !Only to be used for every project...
  USE mod_vel
  USE mod_turb
#ifdef streamxy
  USE mod_streamxy
#endif
#ifdef streamv
  USE mod_streamv
#endif
#ifdef streamr
  USE mod_streamr
#endif
#ifdef tracer
  USE mod_tracer
#endif
#ifdef sediment
  USE mod_sed
#endif
  
  IMPLICIT none
  
#ifdef stat
#ifdef streamxy
  REAL sxyy(IMT,JMT),sxyx(IMT,JMT)
#endif
#ifdef streamv
  REAL sxz(JMT,KM),syz(JMT,KM)
#endif
#ifdef streamr
  REAL sxr(IMT,MR,LOV),syr(JMT,MR,LOV)
#endif
#endif
  
  INTEGER mra,mta,msa
  REAL temp,salt,dens
  
#if defined sediment
  INTEGER nsed,nsusp
  logical res
#endif
  
  INTEGER nrj(ntracmax,NNRJ)
  REAL*8 dsc,ts,trj(ntracmax,NTRJ)
  INTEGER ist,jst,kst
  INTEGER ib,jb,kb,k,ijt,kkt,ijj,kkk,niter,ia,ja,iam,ibm,ka,i,j,m,l,lbas
  INTEGER ntrac,nev,nrh0,nout,nloop,nerror
  INTEGER nnorth,ndrake,ngyre,ntractot,nexit(NEND)
  
  REAL*8 rlon,rlat,x1,y1,z1,x0,y0,z0,tt,dt,dxyz,t0
  REAL*8 ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
  REAL*8 subvol,vol,arc,arct,rr,rb,rg,rbg,uu
  
  INTEGER                                    :: landError=0
  REAL                                       :: fullstamp1 ,fullstamp2
  REAL, DIMENSION(2)                         :: timestamp1 ,timestamp2
  REAL zz
  
  logical scrivi
  
  ! === Error Evaluation ===
  INTEGER                             :: errCode


  
  iday0=iday
  imon0=imon
  iyear0=iyear
  ! === print some run stats ===
  print *,'------------------------------------------------------'  
  print *,'Traj write dir    :  ' ,trim(outDataDir)
  print *,'Time interp steps : ' ,iter
  
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
  nev=0
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
#endif
  
  nrj=0
  trj=0.d0
  
#ifdef streamxy
  sxyy=0.
  sxyx=0.
#endif
#ifdef streamv
  stxz=0.
  styz=0.
#endif
#ifdef streamr
  stxr=0.
  styr=0.
#endif 
  
  dstep=1.d0/dble(iter)
  dtmin=dstep*tseas
  
  !==========================================================
  !=== Read in the end positions from an previous run     === 
  !==========================================================
  
#ifdef rerun
  open(67,file=trim(outDataDir)//outDataFile//'_rerun.asc')
40 continue
  read(67,566,end=41,err=41) ntrac,n,rlon,rlat,zz
  !if(n.ne.1)       print 566,ntrac,n,rlon,rlat,zz
  
#ifdef orc
  do k=1,LBT
     if(ienw(k).le.rlon .and. rlon.le.iene(k)) then
        nrj(ntrac,8)=k                               
     endif
  enddo
  if(nrj(ntrac,8).eq.0) stop 7395                               
#endif
  
  
  
#ifdef occ66
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
#endif
  
  goto 40
41 continue
  
  do ntrac=1,ntracmax
     if(nrj(ntrac,8).eq.0) nrj(ntrac,6)=1 
  enddo
  
#else
  lbas=1 ! set to 1 if no rerun
#endif
  
  !==========================================================
  !=== read ocean/atmosphere GCM data files               ===
  !==========================================================
  print *,'------------------------------------------------------'
  WRITE (6, FMT="(A)", ADVANCE="NO") ' === Reading initial dataset'
  call etime(timestamp1,fullstamp1)
  ff=dble(nff)
  tstep=dble(intstep) 
  ints=intstart
  call readfields   ! initial dataset
  ntrac=0
  call etime(timestamp2,fullstamp2)
  write (6 , FMT="(A,F5.2,A)") ', done in ' ,(fullstamp2-fullstamp1) ,' sec'

  !==========================================================
  !==========================================================
  !=== Start main time loop                               ===
  !==========================================================
  !==========================================================
  
  intsTimeLoop: do ints=intstart+intstep,intstart+intrun,intstep
     
     call readfields
     if(mod(ints,120).eq.0 .and. ints.ne.0) call writepsi ! write psi
#ifdef tracer 
     if(mod(ints,120).eq.0) call writetracer
#endif
     
     !=======================================================
     ! ===    Seed particles if still in intspin.         ===
     !=======================================================
     intspinCond: if(nff*ints <= nff*(intstart+intspin)) then
        ! === Seed particles ===
        ntrac=ntractot
        ijkstloop: do ijk=1,ijkMax
           ist  = ijkst(ijk,1)
           jst  = ijkst(ijk,2)
           kst  = 23!ijkst(ijk,3)
           idir = ijkst(ijk,4)
           isec = ijkst(ijk,5)
           vol  = 0
           
           ib=ist
           ibm=ib-1
           if(ibm.eq.0) ibm=IMT
           jb=jst
           kb=kst
  
           ! === follow trajectory only if velocity   ===
           ! === in right direction + sets trajectory === 
           ! === transport vol.                       ===
           idirCond: if (idir .ne. 0) then
              select case (isec)
              case (1)
                 vol=uflux(ist,jst,kst,1)
              case (2)
                 vol=vflux(ist,jst,kst,1)
              case default
                 call vertvel(1.d0,ib,ibm,jb,kst)
#ifdef full_wflux
                 vol=wflux(ist,jst,kst,1)
#else
                 vol=wflux(kst)
#endif
              end select
              if (idir*ff*vol.le.0.) cycle ijkstloop
              vol = abs(vol)
           else
              select case (isec)
              case (1)
                 ! WARNING
                 ! Should generate trajectories in 
                 ! both directions but doesn't work 
                 stop 5978 
                 vol=abs(uflux(ist,jst,kst,1))
              case (2)
                 stop 5979
                 vol=abs(vflux(ist,jst,kst,1))
              case (3)
                 call vertvel(1.d0,ib,ibm,jb,kst)
#ifdef full_wflux
                 vol=abs(wflux(ist,jst,kst,1))
#else
                 vol=abs(wflux(kst))
#endif
#ifdef twodim
                 vol = 1
#endif
              case(4)
                 if(KM+1-kmt(ist,jst).gt.kst) cycle ijkstloop 
                 if(uflux(ib,jb,kb,1)+uflux(ibm,jb,kb,1) + & 
                      vflux(ib,jb,kb,1)+vflux(ib,jb-1,kb,1).eq.0.) cycle ijkstloop
              end select
              if(vol.eq.0) cycle ijkstloop
           end if idirCond
           
           ! === trajectory volume in m3 ===
           if(nqua.ge.3 .or. isec.eq.4) then
#if defined ifs || atm
              vol=dztb(ib,jb,kb,1)
#else
#if defined sigma 
              vol=dztb(ib,jb,kb)
#else
              vol=dz(kb)
#endif
#if defined occ66 || orc || for || sim || multcol || jplSCB || eccoSOSE
              if(kb.eq.KM+1-kmt(ib,jb) ) vol=dztb(ib,jb,1)
#endif
              if(kb.eq.KM) vol=vol+hs(ib,jb,1)
#endif
              vol=vol*dxdy(ib,jb)
           endif
           
           ! === number of trajectories for box (ist,jst,kst) ===
           select case (nqua)
           case (1)
              num = partQuant
           case (2)
              num = vol/partQuant
           case (3)
              num = vol/partQuant
           case (4)
              num = ijkst(ijk,6)
           end select
           if(num.eq.0) num=1 ! always at least one trajectory
        
           ijt    = nint(sqrt(float(num)))
           kkt    = nint(float(num)/float(ijt))
           subvol = vol/dble(ijt*kkt)
           
           if(subvol.eq.0.) stop 3956  !?????????????????
           if(subvol.eq.0.) subvol=1.
           
99         format(' ib=',i4,' jb=',i3,' kb=',i2,' vol=',f10.0, &
                ' num=',i6,' ijt=',i4,' kkt=',i7,' subvol=',f12.0) 
           
           ! === loop over the subboxes of box (ist,jst,kst) ===
           ijjLoop: do ijj=1,ijt
              kkkLoop: do kkk=1,kkt
                 
                 ib=ist
                 jb=jst
                 kb=kst
                 
                 ! === Meridional section ===
                 select case(isec)
                 case (1)
                    y1=dble(jb-1) + (dble(ijj)-0.5d0)/dble(ijt) 
                    x1=dble(ist) 
                    if(idir.eq. 1) ib=ist+1
                    if(idir.eq.-1) ib=ist 
                    z1=dble(kb-1) + (dble(kkk)-0.5d0)/dble(kkt)
                    ! === Zonal section      ===
                 case (2)
                    x1=dble(ibm) + (dble(ijj)-0.5d0)/dble(ijt)
                    y1=dble(jst) 
                    if(idir.eq. 1) jb=jst+1
                    if(idir.eq.-1) jb=jst 
                    z1=dble(kb-1) + (dble(kkk)-0.5d0)/dble(kkt)
                    ! === Vertical section   ===
                 case (3)
                    x1=dble(ibm ) + (dble(ijj)-0.5d0)/dble(ijt)
                    y1=dble(jb-1) + (dble(kkk)-0.5d0)/dble(kkt) 
                    z1=dble(kb)
                    ! === Spread even inside T-box ===
                 case (4)
                    x1=dble(ibm ) + (dble(ijj)-0.5d0)/dble(ijt)
                    y1=dble(jb-1) + (dble(kkk)-0.5d0)/dble(kkt) 
                    ! z1=dble(kb-1) + (dble(kkk)-0.5d0)/dble(kkt)
                    z1=dble(kb-1) + 0.5d0
                 end select
                 
                 ibm=ib-1
                 ! === cyclic ocean/atmosphere === 
                 if(ibm.eq.0) ibm=IMT
                 if(ib.eq.1.and.x1.gt.dble(IMT)) x1=x1-dble(IMT)
                 
                 ! === check properties of water- ===
                 ! === mass at initial time       === 
#ifndef ifs || atm
#ifdef tempsalt 
                 !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
                 call interp2(ib,jb,kb,ib,jb,kb,temp,salt,dens,1)
                 if(temp.lt.tmin0 .or. temp.gt.tmax0 .or. &
                      salt.lt.smin0 .or. salt.gt.smax0 .or. &
                      dens.lt.rmin0 .or. dens.gt.rmax0) then
                    !print *,'outside mass range',temp,salt,dens
                    cycle kkkLoop 
                 endif
#endif
#endif
                 
                 ntrac=ntrac+1  ! the trajectory number
                 
                 ! selects only one singe trajectory
#ifdef select
                 if(ntrac.ne.57562) then 
                    nrj(ntrac,6)=1
                     nout=nout+1
                    cycle kkkLoop
                 endif
#endif
                 ! === initialise the trajectory iteration number
                 niter=0 
                 call errorCheck('ntracGTntracmax',errCode)
                 
                 ts=ff*dble(ints-intstep)/tstep !time, fractions of ints
                 tt=ts*tseas !time(sec) rel to start
                 
#ifdef streamv 
                 sxz=0.
                 syz=0.
#endif
#ifdef streamr 
                 syr=0.
                 sxr=0.
#endif
#ifdef streamxy
                 sxyx=0.
                 sxyy=0.
#endif
                 ! === initialise time and ===
                 ! === store trajectory positions ===
                 t0=tt
                 dt=0.d0
                 arct=0.
                 
                 trj(ntrac,1)=x1
                 trj(ntrac,2)=y1
                 trj(ntrac,3)=z1
                 trj(ntrac,4)=tt
                 trj(ntrac,5)=subvol
                 trj(ntrac,6)=arct
                 trj(ntrac,7)=t0
                 nrj(ntrac,1)=ib
                 nrj(ntrac,2)=jb
                 nrj(ntrac,3)=kb
                 nrj(ntrac,4)=niter
                 nrj(ntrac,5)=idint(ts)
                 nrj(ntrac,7)=1
                 call writedata(10)
              end do kkkLoop
           end do ijjLoop
        end do ijkstloop
        
        ! ===  End of seeding part ===
     end if intspinCond
     ! ntracin
     ntractot=ntrac
     ntrac=0

     if(ntractot-nout-nerror.eq.0) exit intsTimeLoop
          
     !=======================================================
     !=== Loop over all trajectories and calculate        ===
     !=== a new position for this time step.              ===
     !=======================================================
     ntracLoop: do ntrac=1,ntractot-1        
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
#ifdef rerun
        lbas=nrj(ntrac,8)
        if(lbas.lt.1 .or.lbas.gt.LBT) then
           print *,'lbas=',lbas,'ntrac=',ntrac
           print *,'trj(ntrac,:)=',trj(ntrac,:)
           print *,'nrj(ntrac,:)=',nrj(ntrac,:)
           goto 1500
        endif
#endif 
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
              z1=z1+0.5
              ! z1=z1+0.1  !resusp l�gre i boxen
              trj(ntrac,3)=z1
              ! === change flag to put trajectory back in circulation ===
              nrj(ntrac,6)=0
              nsed=nsed-1
              nsusp=nsusp+1
           else
              cycle ntracLoop 
           endif
        endif
        
#endif        
        ! ===  start loop for each trajectory ===
        scrivi=.false.
        niterLoop: do        
           niter=niter+1 ! iterative step of trajectory
#ifdef sediment
           ! Find settling velocity for active gridbox ===
           call sedvel(temp,dens) 
#endif
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
              goto 1500
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
           ! T-box volume in m3
#if defined ifs || atm
           dxyz=rg*dztb(ib,jb,kb,NST)+rr*dztb(ib,jb,kb,1)
#else
#ifdef sigma
           dxyz=dztb(ib,jb,kb)
#else
           dxyz=dz(kb)
#endif
#if defined occ66 || orc || for || sim || multcol || jplSCB
           if(kb.eq.KM+1-kmt(ib,jb) ) dxyz=dztb(ib,jb,1)
#endif
           if(kb.eq.KM) dxyz=dxyz+rg*hs(ib,jb,NST)+rr*hs(ib,jb,1)
#endif
           dxyz=dxyz*dxdy(ib,jb)
           call errorCheck('dxyzError'     ,errCode)
           call errorCheck('coordBoxError' ,errCode)
           call errorCheck('infLoopError'  ,errCode)
           if (errCode.ne.0) cycle ntracLoop
           
           ! === calculate the turbulent velocities ===
#ifdef turb
           call turbuflux(ia,ja,ka,rr)
#endif
           ! === calculate the vertical velocity ===
           call vertvel(rr,ia,iam,ja,ka)
           
           ! === write trajectory ===                       
#ifdef tracer
           if(ts.eq.dble(idint(ts))) then 
              tra(ia,ja,ka)=tra(ia,ja,ka)+real(subvol)
           end if
#endif
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
           call cross(1,ia,ja,ka,x0,dse,dsw,rr) ! zonal
           call cross(2,ia,ja,ka,y0,dsn,dss,rr) ! meridional
           call cross(3,ia,ja,ka,z0,dsu,dsd,rr) ! vertical
           
           dsmin=dtmin/dxyz
           ds=min(dse,dsw,dsn,dss,dsu,dsd,dsmin)
           !print *,'min',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin
           if(ds.eq.1.d20 .or.ds.eq.0.d0)then 
              ! === Can not find any path for unknown reasons ===
              print *,'ds cross error',dse,dsw,dsn,dss,dsu,dsd,dsmin,dxyz
              print *,ia,ja,ka,x0,y0,z0,rr,ntrac,niter
              print *,'k=',ka,kb,KM+1-kmt(ia,ja),kmt(ia,ja)
              ! goto 1500
              nerror=nerror+1
              nrj(ntrac,6)=1
              cycle ntracLoop
           endif
           
           dt=ds*dxyz ! transform ds to dt in seconds
           if(ds.eq.dsmin) dt=dtmin  ! this makes dt more accurate
           if(dt.lt.0.d0) then
              print *,'dt=',dt
              goto 1500
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
              if(dt.eq.dtmin) then
                 ts=ts+dstep
                 tss=tss+1.d0
              else
                 ts=ts+dt/tseas
                 tss=tss+dt/tseas*dble(iter)
              endif
           endif
           
           ! === time interpolation constant ===
           rbg=dmod(ts,1.d0) 
           rb =1.d0-rbg
           
           ! === calculate the new positions ===
           ! === of the trajectory           ===    
           if(ds.eq.dse) then ! eastward exit 
              uu=(rbg*uflux(ia,ja,ka,NST)+rb*uflux(ia ,ja,ka,1))*ff
#ifdef turb    
              ! uu=uu+upr(1,2)
#endif
              if(uu.gt.0.d0) then
                 ib=ia+1
                 if(ib.gt.IMT) ib=ib-IMT 
              endif
              x1=dble(ia)
              call pos(2,ia,ja,ka,y0,y1,ds,rr) 
              call pos(3,ia,ja,ka,z0,z1,ds,rr)
              scrivi=.true.
#ifdef streamxy
              ! === zonal component stream function ===
              stxyx(ia,ja,lbas)=stxyx(ia,ja,lbas)+real(subvol*ff)
#endif
#ifdef streamv
              ! === overturning "z" stream function ===
              stxz(ia,ka,lbas)=stxz(ia,ka,lbas)+real(subvol*ff)
#endif
#ifdef streamr
              call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
              mra=nint((dens-rmin)/dr)+1  
              if(mra.lt.1 ) mra=1
              if(mra.gt.MR) mra=MR
              ! === overturn "rho" stream function  ===
              stxr(ia,mra,lbas,1)=stxr(ia,mra,lbas,1)+real(subvol*ff)
#ifdef streamts
              mta=(temp-tmin)/dtemp+1  
              if(mta.lt.1 ) mta=1
              if(mta.gt.MR) mta=MR
              msa=(salt-smin)/dsalt+1  
              if(msa.lt.1 ) msa=1
              if(msa.gt.MR) msa=MR
              ! === overturn "temp" stream function ===
              stxr(ia,mta,lbas,2)=stxr(ia,mta,lbas,2)+real(subvol*ff)
              ! === overturn "salt" stream function ===
              stxr(ia,msa,lbas,3)=stxr(ia,msa,lbas,3)+real(subvol*ff)
#endif
#endif
           elseif(ds.eq.dsw) then ! westward exit
              
              uu=(rbg*uflux(iam,ja,ka,NST)+rb*uflux(iam,ja,ka,1))*ff
#ifdef turb    
              ! uu=uu+upr(2,2)
#endif
              if(uu.lt.0.d0) then
                 ib=iam
              endif
              x1=dble(iam)
              call pos(2,ia,ja,ka,y0,y1,ds,rr) ! meridional position
              call pos(3,ia,ja,ka,z0,z1,ds,rr) ! vertical position
              scrivi=.true.
#ifdef streamxy
              stxyx(iam,ja,lbas)=stxyx(iam,ja,lbas)-real(subvol*ff) 
#endif
#ifdef streamv
              stxz(iam,ka,lbas)=stxz(iam,ka,lbas)-real(subvol*ff)
#endif
#ifdef streamr
              call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
              mra=nint((dens-rmin)/dr)+1  
              if(mra.lt.1) mra=1
              if(mra.gt.MR) mra=MR
              stxr(iam,mra,lbas,1)=stxr(iam,mra,lbas,1)-real(subvol*ff)
#ifdef streamts
              mta=nint((temp-tmin)/dtemp)+1  
              if(mta.lt.1 ) mta=1
              if(mta.gt.MR) mta=MR
              msa=nint((salt-smin)/dsalt)+1  
              if(msa.lt.1 ) msa=1
              if(msa.gt.MR) msa=MR
              stxr(iam,mta,lbas,2)=stxr(iam,mta,lbas,2)-real(subvol*ff)
              stxr(iam,msa,lbas,3)=stxr(iam,msa,lbas,3)-real(subvol*ff)
#endif
#endif
              
           elseif(ds.eq.dsn) then ! northward exit
              
              uu=(rbg*vflux(ia,ja,ka,NST)+rb*vflux(ia,ja,ka,1))*ff
#ifdef turb    
              ! uu=uu+upr(3,2)
#endif
              if(uu.gt.0.d0) then
                 jb=ja+1
              endif
              y1=dble(ja)
              call pos(1,ia,ja,ka,x0,x1,ds,rr)
              call pos(3,ia,ja,ka,z0,z1,ds,rr)
              scrivi=.true.
#ifdef streamxy
              stxyy(ia,ja,lbas)=stxyy(ia,ja,lbas)+real(subvol*ff)
#endif
#ifdef streamv 
              styz(ja,ka,lbas)=styz(ja,ka,lbas)+real(subvol*ff)
              ! print *,'n',styz(ja,ka,lbas),subvol*ff,ja,ka,lbas
#endif
#ifdef streamr
              call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
              mra=nint((dens-rmin)/dr)+1  
              if(mra.lt.1) mra=1
              if(mra.gt.MR) mra=MR
              styr(ja,mra,lbas,1)=styr(ja,mra,lbas,1)+real(subvol*ff)
#ifdef streamts
              mta=nint((temp-tmin)/dtemp)+1  
              if(mta.lt.1 ) mta=1
              if(mta.gt.MR) mta=MR
              msa=nint((salt-smin)/dsalt)+1  
              if(msa.lt.1 ) msa=1
              if(msa.gt.MR) msa=MR
              styr(ja,mta,lbas,2)=styr(ja,mta,lbas,2)+real(subvol*ff)
              styr(ja,msa,lbas,3)=styr(ja,msa,lbas,3)+real(subvol*ff)
#endif
#endif
              
           elseif(ds.eq.dss) then ! southward exit
              
              uu=(rbg*vflux(ia,ja-1,ka,NST)+rb*vflux(ia,ja-1,ka,1))*ff
#ifdef turb    
              ! uu=uu+upr(4,2)
#endif
              if(uu.lt.0.d0) then
                 jb=ja-1
                 if(jb.eq.0) stop 34578
              endif
              y1=dble(ja-1)
              call pos(1,ia,ja,ka,x0,x1,ds,rr)
              call pos(3,ia,ja,ka,z0,z1,ds,rr)
              scrivi=.true.
#ifdef streamxy
              stxyy(ia,ja-1,lbas)=stxyy(ia,ja-1,lbas)-real(subvol*ff)
#endif
#ifdef streamv 
              styz(ja-1,ka,lbas)=styz(ja-1,ka,lbas)-real(subvol*ff)
#endif
#ifdef streamr 
              call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
              mra=nint((dens-rmin)/dr)+1  
              if(mra.lt.1) mra=1
              if(mra.gt.MR) mra=MR
              styr(ja-1,mra,lbas,1)=styr(ja-1,mra,lbas,1)-real(subvol*ff)
#ifdef streamts
              mta=nint((temp-tmin)/dtemp)+1  
              if(mta.lt.1 ) mta=1
              if(mta.gt.MR) mta=MR
              msa=nint((salt-smin)/dsalt)+1  
              if(msa.lt.1 ) msa=1
              if(msa.gt.MR) msa=MR
              styr(ja-1,mta,lbas,2)=styr(ja-1,mta,lbas,2)-real(subvol*ff)
              styr(ja-1,msa,lbas,3)=styr(ja-1,msa,lbas,3)-real(subvol*ff)
#endif
#endif
              
           elseif(ds.eq.dsu) then ! upward exit
              
              scrivi=.false.
              call vertvel(rb,ia,iam,ja,ka)
#ifdef full_wflux
              uu=wflux(ia,ja,ka,1)
#else
              uu=wflux(ka)
#endif
#ifdef turb    
              ! uu=uu+upr(5,2)
#endif
              if(uu.gt.0.d0) then
                 kb=ka+1
              endif
              z1=dble(ka)
              if(kb.eq.KM+1) then
                 nev=nev+1
                 kb=KM
                 z1=dble(KM)-0.5
              endif
              call pos(1,ia,ja,ka,x0,x1,ds,rr)
              call pos(2,ia,ja,ka,y0,y1,ds,rr)
              
           elseif(ds.eq.dsd) then ! downward exit
              
              call vertvel(rb,ia,iam,ja,ka)

#ifdef full_wflux
              if(wflux(ia,ja,ka-1,1).lt.0.d0) kb=ka-1
#else
              if(wflux(ka-1).lt.0.d0) kb=ka-1
#endif              
              z1=dble(ka-1)
              call pos(1,ia,ja,ka,x0,x1,ds,rr)  
              call pos(2,ia,ja,ka,y0,y1,ds,rr) 
#ifdef sediment
              if(kb.eq.KM-kmt(ia,ja)) then
                 nsed=nsed+1
                 nrj(ntrac,6)=2
                 trj(ntrac,1)=x1
                 trj(ntrac,2)=y1
                 trj(ntrac,3)=z1
                 trj(ntrac,4)=tt
                 trj(ntrac,5)=subvol
                 trj(ntrac,6)=arct
                 nrj(ntrac,1)=ib
                 nrj(ntrac,2)=jb
                 nrj(ntrac,3)=ka
                 nrj(ntrac,4)=n
                 nrj(ntrac,5)=idint(ts)
                 nrj(ntrac,7)=1
                 !call writedata(13)
                 !cycle ntracLoop
              endif
#endif
              
           elseif( ds.eq.dsc .or. ds.eq.dsmin ) then  ! inter time steping 
              scrivi=.false.
              call pos(1,ia,ja,ka,x0,x1,ds,rr) ! zonal crossing 
              call pos(2,ia,ja,ka,y0,y1,ds,rr) ! merid. crossing 
              call pos(3,ia,ja,ka,z0,z1,ds,rr) ! vert. crossing 
           endif
           ! === make sure that trajectory ===
           ! === is inside ib,jb,kb box    ===
           
#ifdef orc 
           if(y1.eq.dble(JMT)) then ! north fold cyclic
              x1=722.d0-x1
              y1=dble(JMT-2)
              jb=JMT-2
           endif
#endif
           
           if(x1.lt.0.d0) x1=x1+dble(IMT)      ! east-west cyclic
           if(x1.gt.dble(IMT)) x1=x1-dble(IMT) ! east-west cyclic
           if(x1.ne.dble(idint(x1))) ib=idint(x1)+1
           if(ib.gt.IMT) ib=ib-IMT             ! east-west cyclic
           if(y1.ne.dble(idint(y1))) jb=idint(y1)+1
           
!           print *,z1,km,kmt(ib,jb)
!           print*,'---'
           if( z1.le.dble(KM-kmt(ib,jb)) ) then
!              print *,z1,km,kmt(ib,jb)
              z1=dble(KM-kmt(ib,jb))+0.5d0
!              print *,'==='
           end if
           if( z1.ge.dble(KM) ) then
              z1=dble(KM)-0.5
              kb=KM
           endif
           if(z1.ne.dble(idint(z1))) then
              kb=idint(z1)+1
              if(kb.eq.KM+1) kb=KM  ! ska nog bort
           endif
           
           call errorCheck('landError', errCode)
           if (errCode.ne.0) cycle ntracLoop
           
           ! === diffusion, which adds a random position ===
           ! === position to the new trajectory          ===
#if defined diffusion     
           call diffusion(x1,y1,z1,ib,jb,kb,dt,snew,st0,st1)
           
#endif
           
           ! === Calculate arclength of the ===
           ! === trajectory path in the box ===
           call arclength(ia,ja,ka,dt,rr,arc)
           arct=arct+arc*arcscale  ! orig arc in meters -> 100 km
           ! === end trajectory if outside chosen domain ===
#if defined occam25 || occ66
           ! === stop and select stream function ===
           if( y1.eq.dble(jmt-2) .and. n.ne.1 ) then ! To Northern Boundary
              nnorth=nnorth+1
              nexit(1)=nexit(1)+1
           elseif(x1.eq.dble(idr) .and. ja.lt.94 .and. n.ne.1 ) then ! or to Drake
              ndrake=ndrake+1
              call writedata(15)
              nexit(2)=nexit(2)+1
              cycle niterLoop
              ! or continue trajectory
           else
              cycle niterLoop                                   
           endif
#else     
           do k=1,LBT
              if(ienw(k).le.ib .and. ib.le.iene(k) .and. &
                   jens(k).le.jb .and. jb.le.jenn(k)  ) then
                 nexit(k)=nexit(k)+1
                 exit niterLoop                                
              endif
           enddo
           ! === stop trajectory if the choosen time or ===
           ! === water mass properties are exceeded     ===
           if(tt-t0.ge.timax) then
              nexit(NEND)=nexit(NEND)+1
              exit niterLoop
           endif
           call writedata(18)
        end do niterLoop
#endif
        ! add streamfuction contribution at the end of trajectory for stat
#ifdef streamxy
        stxyx(:,:,lbas)=stxyx(:,:,lbas)+sxyx(:,:)
        stxyy(:,:,lbas)=stxyy(:,:,lbas)+sxyy(:,:)
#endif
#ifdef streamv
        stxz(:,:,lbas)=stxz(:,:,lbas)+sxz(:,:)
        styz(:,:,lbas)=styz(:,:,lbas)+syz(:,:)
#endif
#ifdef streamr
        stxr(:,:,lbas,:)=stxr(:,:,lbas,:)+sxr(:,:,:)
        styr(:,:,lbas,:)=styr(:,:,lbas,:)+syr(:,:,:)
#endif  
        
        nout=nout+1
        
        call writedata(17)
        nrj(ntrac,6)=1
     end do ntracLoop
     
#ifdef sediment
     print 599,ints,ntime,ntractot,nout,nloop,nerror,ntractot-nout,nsed,nsusp,nexit
599  format('ints=',i7,' time=',i10,' ntractot=',i8,' nout=',i8,' nloop=',i4, &
          ' nerror=',i4,' in ocean/atm=',i8,' nsed=',i8, ' nsusp=',i8,' nexit=',9i8)
#else
     print 599,ints,ntime,ntractot,nout,nloop,nerror,ntractot-nout-nerror,nexit
599  format('ints=',i7,' time=',i10,' ntractot=',i8,' nout=',i8,' nloop=',i4, &
         ' nerror=',i4,' in ocean/atm=',i8,' nexit=',9i8)
#endif
     
     
  end do intsTimeLoop
  
1500 close(56)
  
  print *,ntractot ,' trajectories calculated'
  print *,nev      ,' trajectories evaporated'
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
  
  print *,'The very end of tracmass run ',outDataFile,' at'
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
          if(dxyz.eq.0.) then
             print *,'====================================='
             print *,'ERROR: dxyz is zero'
             print *,'-------------------------------------'
             print *,'ntrac=',ntrac,' ints=', ints
             print *,'ib=',ib,'jb=',jb,'kb=',kb
             print *,'dxyz=',dxyz,' dxdy=',dxdy(ib,jb)
             print *,'dztb=',dztb(ib,jb,kb)
             print *,'rg*hs=',rg*hs(ib,jb,NST)
             print *,'rr*hs=',rr*hs(ib,jb,1)
             print *,'-------------------------------------'
             print *,'The run is terminated'
             print *,'====================================='
             errCode = -39
             stop
          endif          
       case ('landError')
          if(kmt(ib,jb).eq.0) then
             if (verboseMess == 1) then
                print *,'====================================='
                print *,'Warning: Trajectory on land'
                print *,'-------------------------------------'
                print *,'land',ia,ib,ja,jb,ka,kb,kmt(ia,ja)
                print *,'xyz',x0,x1,y0,y1,z0,z1
                print *,'ds',dse,dsw,dsn,dss,dsu,dsd
                print *,'dsmin=',ds,dsmin,dtmin,dxyz
                print *,'tt=',tt,ts
                print *,'ntrac=',ntrac
#ifdef turb
                print *,'upr=',upr
#endif
                print *,'-------------------------------------'
                print *,'The trajectory is killed'
                print *,'====================================='
             end if
             call writedata(14)
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
             x1=dble(ib-1)+0.5
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
             print *,'====================================='
             print *,'Warning: Particle in infinite loop '
             print *,'niter:',niter,'nrj:',nrj(ntrac,4)
             print *,'dxdy:',dxdy(ib,jb),'dxyz:',dxyz
             print *,'kmt:',kmt(ia-1,ja-1),'dz(k):',dz(ka-1)
             print *,'x1:',x1,'u:',uflux(ia-1,ja-1,ka-1,2)
             print *,'y1:',y1,'v:',vflux(ia-1,ja-1,ka-1,2)
             print *,'kb:',kb-1,'z1:',z1
             print *,'-------------------------------------'
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
             ! nrj(ntrac,6)=1  ! 1=end trajectory
             nrj(ntrac,7)=1
             errCode = -48
          end if
       end select
  end subroutine errorCheck

  subroutine writedata(sel)
    REAL                                 :: vort
    INTEGER                              :: sel ,xf ,yf ,zf ,n
    INTEGER, SAVE                        :: recPosIn=0  ,recPosOut=0
    INTEGER, SAVE                        :: recPosRun=0 ,recPosErr=0
    REAL                                 :: x14 ,y14 ,z14
#if defined for || sim 
566 format(i8,i7,f7.2,f7.2,f7.1,f10.2,f10.2 &
         ,f10.1,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#elif defined rco 
566 format(i8,i7,f7.2,f7.2,f7.1,f10.0,f10.0 &
         ,f10.0,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#else
566 format(i7,i7,f7.2,f7.2,f7.1,f10.4,f10.4 &
         ,f13.8,f6.2,f6.2,f6.2,f6.0,8e8.1 )
#endif
    
    xf   = floor(x1)
    yf   = floor(y1)
    zf   = floor(z1)
    
    vort = (vvel(xf+1,yf,zf)-vvel(xf-1,yf,zf))/4000 - &
         (uvel(xf,yf+1,zf)-uvel(xf,yf-1,zf))/4000   
    
#if defined textwrite 
    select case (sel)
    case (10)
       write(58,566) ntrac,niter,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens
    case (11)
       if( (kriva.eq.1 .and. ts.eq.dble(idint(ts)) ) .or. &
            (scrivi .and. kriva.eq.2)                .or. &
            (kriva.eq.3)                             .or. &
            (kriva.eq.4 .and. niter.eq.1)            .or. &
            (kriva.eq.5 .and. &
            (tt-t0.eq.7.*tday.or.tt-t0.eq.14.*tday & 
            .or.tt-t0.eq.21.*tday)) ) then
          call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
#if defined biol
          write(56,566) ntrac,ints,x1,y1,z1,tt/3600.,t0/3600.
#else
          !write(56,566) ntrac,ints,x1,y1,z1, & 
          ! tt/tday,t0/tday,subvol,temp,salt,dens,arct
          write(56,566) ntrac,ints,x1,y1,z1, &
               uvel(xf,yf,zf),vvel(xf,yf,zf) &
               ,vort,temp,salt,dens,arct
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
          call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
          write(56,566) ntrac,ints,x1,y1,z1, &
               tt/tday,t0/tday,subvol,temp,salt,dens,arct
       end if
    case (17)
       write(57,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (18)
       if( kriva.ne.0 .and. ts.eq.dble(idint(ts)) .and. &
            ints.eq.intstart+intrun) then 
          call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
          !write(56,566) ntrac,ints,x1,y1,z1, &
          !tt/tday,t0/tday,subvol,temp,salt,dens,arct
          write(56,566) ntrac,ints,x1,y1,z1, & 
               uvel(xf,yf,zf),vvel(xf,yf,zf),vort,temp,salt,dens,arct
          ! write(56,566) ntrac,niter,x1,y1,z1,tt/3600.,t0/3600.
          !,subvol,temp,salt,dens,arct
       endif
    case (19)
       ! === write last sedimentation positions ===
       open(34,file=trim(outDataDir)//name//'_sed.asc') 
       do n=1,ntracmax
          if(nrj(n,1).ne.0) then
             write(34,566) n,nrj(n,4),trj(n,1),trj(n,2), & 
                  trj(n,3),trj(n,4)/tday,trj(n,7)/tday
          endif
       enddo
       close(34)
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
       if( (kriva.eq.1 .and. ts .eq. ints-1) .or. &
            (scrivi .and. kriva.eq.2)                .or. &
            (kriva.eq.3)                             .or. &
            (kriva.eq.4 .and. niter.eq.1)            .or. &
            (kriva.eq.5 .and. &
            (tt-t0.eq.7.*tday.or.tt-t0.eq.14.*tday & 
            .or.tt-t0.eq.21.*tday)) ) then
          call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
          recPosRun = recPosRun+1
          write(unit=76 ,rec=recPosRun) ntrac,ints,x14,y14,z14       
       end if
    case (13)
       recPosOut = recPosOut+1
       write(unit=77 ,rec=recPosOut) ntrac,ints,x14,y14,z14   
    case (14)
       recPosRun = recPosRun+1
       write(unit=76 ,rec=recPosRun) ntrac,ints,x14,y14,z14   
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
    case (40)
       recPosErr=recPosErr+1    
       write(unit=79 ,rec=recPosErr) ntrac,ints,x14,y14,z14   
    end select
#endif    


!!$#if defined binwrite 
!!$    select case (sel)
!!$    case (10)
!!$       write(78) ntrac,ints,x1,y1,z1
!!$    case (11)
!!$       if( (kriva.eq.1 .and. ts.eq.dble(idint(ts)) ) .or. &
!!$            (scrivi .and. kriva.eq.2)                .or. &
!!$            (kriva.eq.3)                             .or. &
!!$            (kriva.eq.4 .and. niter.eq.1)            .or. &
!!$            (kriva.eq.5 .and. &
!!$            (tt-t0.eq.7.*tday.or.tt-t0.eq.14.*tday & 
!!$            .or.tt-t0.eq.21.*tday)) ) then
!!$          call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
!!$          write(76) ntrac,ints,x1,y1,z1
!!$       end if
!!$    case (13)
!!$       write(77) ntrac,ints,x1,y1,z1
!!$    case (14)
!!$       write(76) ntrac,ints,x1,y1,z1§
!!$    case (15)
!!$       write(76) ntrac,ints,x1,y1,z1
!!$    case (16)
!!$       if(kriva.ne.0 ) then
!!$          write(76) ntrac,ints,x1,y1,z1
!!$       end if
!!$    case (17)
!!$       write(77) ntrac,ints,x1,y1,z1
!!$    case (18)
!!$       if( kriva.ne.0 .and. ts.eq.dble(idint(ts)) .and. &
!!$            ints.eq.intstart+intrun) then 
!!$          call interp2(ib,jb,kb,ia,ja,ka,temp,salt,dens,1)
!!$          write(76) ntrac,ints,x1,y1,z1
!!$       endif
!!$       !case (19)
!!$    end select
!!$#endif    


  end subroutine writedata
  
end subroutine loop
