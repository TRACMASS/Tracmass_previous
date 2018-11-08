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
  trajectories(ntrac)%lbas = 1
  if(trajectories(ntrac)%lbas == 0) stop 7395                               
566 format(i8,i7,2f9.3,f6.2,2f10.2 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
#elif defined  occ66

  if(rlon.eq.293.) then
     trajectories(ntrac)%lbas = 3    ! Drake Passage
  elseif(rlat.eq.-47.00) then
     trajectories(ntrac)%lbas = 2    ! Mid Gyre  (47S)
  elseif(rlat.eq.-29.25) then
     trajectories(ntrac)%lbas = 1    ! Northern boundary (30S)
  else
     trajectories(ntrac)%lbas = 0
     print 566,ntrac,n,rlon,rlat,zz
     stop 4957
  endif

#elif defined ifs
  if(rlat.eq.float(jenn(1))) then
     trajectories(ntrac)%lbas = 1    ! Southern boundary
  elseif(rlat.eq.float(jens(2))) then
     trajectories(ntrac)%lbas = 2    ! Northern boundary 
  else
     trajectories(ntrac)%lbas=0
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
     !! if no exit boundary, then kill trajectory
     if(trajectories(ntrac)%lbas == 0) trajectories(ntrac)%active=.false. 
  enddo
  
#else
  lbas=1 ! set to 1 if no rerun
#endif /*rerun*/
