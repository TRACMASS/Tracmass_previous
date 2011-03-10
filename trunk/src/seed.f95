module mod_seed
  
  USE mod_param
  USE mod_time
  USE mod_grid
  USE mod_buoyancy
  USE mod_vel
  USE mod_traj
  
  implicit none
 
  !from mod_seed in modules
  INTEGER                                    :: nff,isec,idir,nqua,num
  INTEGER                                    :: ijk  ,ijkMax
  INTEGER                                    :: seedType ,varSeedFile 
  INTEGER                                    :: ist1 ,ist2   ,jst1 ,jst2
  INTEGER                                    :: kst1, kst2
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: ijkst
  REAL, ALLOCATABLE, DIMENSION(:,:)          :: xyzst
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:)     :: seedMask
  CHARACTER(LEN=200)                         :: seedDir
  CHARACTER(LEN=200)                         :: seedFile
  INTEGER                                    :: ist,jst,kst
  INTEGER                                    :: ijt,kkt,ijj,kkk
  INTEGER                                    :: ntractot

  
CONTAINS


  subroutine seed (tt,ts)
    INTEGER                                  :: errCode
    INTEGER                                  :: ib, jb, kb, ibm
    INTEGER                                  :: i, j, k, l, m
    INTEGER                                  :: ntrac
    REAL*8                                   :: tt,ts
    REAL*8                                   :: x1, y1, z1
    REAL                                     :: temp,salt,dens
    REAL*8                                   :: vol, subvol

    ntrac=ntractot

    ijkstloop: do ijk=1,ijkMax
       if(seedType == 3) then
           ist = int(xyzst(ijk,1))+1
           jst = int(xyzst(ijk,2))+1
           kst = int(xyzst(ijk,3))+1
           idir = int(xyzst(ijk,4))
           isec = int(xyzst(ijk,5))
       else
           ist  = ijkst(ijk,1)
           jst  = ijkst(ijk,2)
           kst  = ijkst(ijk,3)
           idir = ijkst(ijk,4)
           isec = ijkst(ijk,5)
       endif
       vol  = 0
       ib=ist
       ibm=ib-1
       if(ibm.eq.0) ibm=IMT
       jb=jst
       kb=kst

       ! === follow trajectory only if velocity in right direction  ===
       ! === & sets trajectory  transport vol.                      === 
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
             vol=wflux(kst,1)
#endif /*full_wflux*/
          end select
          if (idir*ff*vol.le.0.d0) cycle ijkstloop

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
           !  stop 5979
           vol=abs(vflux(ist,jst,kst,1))
        case (3)
           call vertvel(1.d0,ib,ibm,jb,kst)
#ifdef full_wflux
           vol=abs(wflux(ist,jst,kst,1))
#else
           vol=abs(wflux(kst,1))
#endif /*full_wflux*/
#ifdef twodim
           !                 vol = 1
#endif /*twodim*/
           vol=1
        case(4)
           if(KM+1-kmt(ist,jst).gt.kst) cycle ijkstloop 
           vol=abs(uflux(ib,jb,kb,1))+abs(uflux(ibm,jb  ,kb,1)) + & 
                abs(vflux(ib,jb,kb,1))+abs(vflux(ib ,jb-1,kb,1))
           !              print *,'KM..',ib,jb,kb,vol,uflux(ib,jb,kb,1)
           if(vol.eq.0.d0) cycle ijkstloop
        case(5)
           if(KM+1-kmt(ist,jst).gt.kst) cycle ijkstloop 
           vol=abs(uflux(ib,jb,kb,1))+abs(uflux(ibm,jb  ,kb,1)) + & 
                abs(vflux(ib,jb,kb,1))+abs(vflux(ib ,jb-1,kb,1))
           !                 print *,'KM..',ib,jb,kb,vol,uflux(ib,jb,kb,1)
           if(vol.eq.0.d0) cycle ijkstloop
        end select
        if(vol.eq.0) cycle ijkstloop
     end if idirCond
     
     ! === trajectory volume in m3 ===
     if(nqua.ge.3 .or. isec.ge.4) then
#ifdef zgrid3Dt
        vol=dzt(ib,jb,kb,1)
 
#elif  zgrid3D
        vol=dzt(ib,jb,kb)
#elif  zgrid1D
        vol=dz(kb)
#endif /*zgrid*/
#ifdef varbottombox
        if(kb.eq.KM+1-kmt(ib,jb) ) vol=dztb(ib,jb,1)
#endif /*varbottombox*/
#ifdef freesurface
        if(kb.eq.KM) vol=vol+hs(ib,jb,1)
        vol=vol*dxdy(ib,jb)
#endif /*freesurface*/
     end if
          
     ! === number of trajectories for box (ist,jst,kst) ===
     select case (nqua)
     case (1)
        num = partQuant
     case (2)
        num = vol/partQuant
     case (3)
        num = vol/partQuant
     case (5)
        if(seedType == 3) then
            num = 1
        else
            num = ijkst(ijk,6)
        endif
     end select
     if(num.eq.0 .and. nqua.ne.5) num=1 ! always at least one trajectory
     
     
     ijt    = nint(sqrt(float(num)))
     kkt    = nint(float(num)/float(ijt))
     subvol = vol/dble(ijt*kkt)
     
     if(subvol.eq.0.d0) stop 3956 
99   format(' ib=',i4,' jb=',i3,' kb=',i2,' vol=',f10.0, &
          ' num=',i6,' ijt=',i4,' kkt=',i7,' subvol=',f12.0) 
     
     ! === loop over the subboxes of box (ist,jst,kst) ===
     ijjLoop: do ijj=1,ijt
        kkkLoop: do kkk=1,kkt          
           ib=ist
           jb=jst
           kb=kst
           
           select case(isec)
           case (1)
              ! === Meridional section ===
              y1=dble(jb-1) + (dble(ijj)-0.5d0)/dble(ijt) 
              x1=dble(ist) 
              if(idir.eq. 1) ib=ist+1
              if(idir.eq.-1) ib=ist 
              z1=dble(kb-1) + (dble(kkk)-0.5d0)/dble(kkt)
           case (2)
              ! === Zonal section      ===
              x1=dble(ibm) + (dble(ijj)-0.5d0)/dble(ijt)
              y1=dble(jst) 
              if(idir.eq. 1) jb=jst+1
              if(idir.eq.-1) jb=jst 
              z1=dble(kb-1) + (dble(kkk)-0.5d0)/dble(kkt)
           case (3)
              ! === Vertical section   ===
              x1=dble(ibm ) + (dble(ijj)-0.5d0)/dble(ijt)
              y1=dble(jb-1) + (dble(kkk)-0.5d0)/dble(kkt) 
              z1=dble(kb)
           case (4)
              ! === Spread even inside T-box ===
              x1=dble(ibm ) + 0.25d0*(dble(ijj)-0.5d0)/dble(ijt)
              y1=dble(jb-1) + 0.25d0*(dble(kkk)-0.5d0)/dble(kkt) 
              ! z1=dble(kb-1) + (dble(kkk)-0.5d0)/dble(kkt) ! spread verticaly
              z1=dble(kb-1) + 0.5d0                         ! or on a fixed depth
           case (5)
              if(isec.ne.nqua .or. nqua.ne.5 .or.ijt.ne.1 .or.kkt.ne.1) stop 8461
              ! === Start particles from exact positions set by a file read in readfield 
              ! === (works only for orc at the moment)
              ! trj(ijk,1)=dble(ibm ) + 0.25d0*(dble(ijj)-0.5d0)/dble(ijt)
              ! trj(ijk,2)=dble(jb-1) + 0.25d0*(dble(kkk)-0.5d0)/dble(kkt) 
              ! trj(ijk,3)=dble(kb-1) + 0.5d0   
              x1=trj(ijk,1)
              y1=trj(ijk,2) 
              z1=trj(ijk,3)
           end select
           
           if(seedType == 3) then
              x1 = xyzst(ijk,1)
              y1 = xyzst(ijk,2)
              z1 = xyzst(ijk,3)
           endif
           
           ibm=ib-1
           ! === cyclic ocean/atmosphere === 
           if(ibm.eq.0) ibm=IMT
           if(ib.eq.1.and.x1.gt.dble(IMT)) x1=x1-dble(IMT)
           
           ! === check properties of water- ===
           ! === mass at initial time       === 
#ifndef ifs 
#ifdef tempsalt 
           call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
           !call interp2(ib,jb,kb,ib,jb,kb,temp,salt,dens,1)
           if(temp.lt.tmin0 .or. temp.gt.tmax0 .or. &
                salt.lt.smin0 .or. salt.gt.smax0 .or. &
                dens.lt.rmin0 .or. dens.gt.rmax0) then
              !print *,'outside mass range',temp,salt,dens
              cycle kkkLoop 
           endif
#endif /*tempsalt*/
#endif
           
           ntrac=ntrac+1  ! the trajectory number
           !                 print *,'ntrac=',ntrac,ijk
           
#ifdef select
           ! selects only one singe trajectory
           if(ntrac.ne.57562) then 
              nrj(ntrac,6)=1
              cycle kkkLoop
           endif
#endif /*select*/
           ts=ff*dble(ints-intstep)/tstep !time, fractions of ints
           tt=ts*tseas !time(sec) rel to start
           
           ! === initialise time and ===
           ! === store trajectory positions ===
       
           trj(ntrac,1)=x1
           trj(ntrac,2)=y1
           trj(ntrac,3)=z1
           trj(ntrac,4)=tt
           trj(ntrac,5)=subvol
           trj(ntrac,6)=0
           trj(ntrac,7)=tt
           nrj(ntrac,1)=ib
           nrj(ntrac,2)=jb
           nrj(ntrac,3)=kb
           nrj(ntrac,4) = 0
           nrj(ntrac,5)=idint(ts)
           nrj(ntrac,7)=1
        end do kkkLoop
     end do ijjLoop
  end do ijkstloop

  ntractot = ntrac

end subroutine seed

end module mod_seed
