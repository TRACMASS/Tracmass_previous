
MODULE mod_grid
  REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt
  
  REAL*8                                    :: rmin ,tmin ,smin
  REAL*8                                    :: dr ,dtemp ,dsalt
  REAL*8                                    :: arcscale
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  CHARACTER(LEN=200)                        :: SubGridFile 
ENDMODULE mod_grid

MODULE mod_vel
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: uflux ,vflux
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)    :: wflux
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: hs
  REAL*8                                     :: ff
ENDMODULE mod_vel

MODULE mod_grid
  REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt
  
  REAL*8                                    :: rmin ,tmin ,smin
  REAL*8                                    :: dr ,dtemp ,dsalt
  REAL*8                                    :: arcscale
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  CHARACTER(LEN=200)                        :: SubGridFile 
ENDMODULE mod_grid






MODULE tracmass
  
  USE mod_param
  USE mod_vel
  USE mod_grid
  IMPLICIT NONE
  
  uflux =
  vflux =
  wflux =
  
  

  x1     = trj(ntrac,1)         ! Current grid position  x
  y1     = trj(ntrac,2)         ! Current grid position  y
  z1     = trj(ntrac,3)         ! Current grid position  z
  tt     = trj(ntrac,4)
  subvol = trj(ntrac,5)
  arct   = trj(ntrac,6)
  t0     = trj(ntrac,7)
  ib     = nrj(ntrac,1)         ! Current grid index  x
  jb     = nrj(ntrac,2)         ! Current grid index  y
  kb     = nrj(ntrac,3)         ! Current grid index  z
  niter  = nrj(ntrac,4)
  ts     = dble(nrj(ntrac,5))
  tss    = 0.d0


 x0=x1
 y0=y1
 z0=z1
 ia=ib
 iam=ia-1
 if(iam.eq.0)iam=IMT
 ja=jb
 ka=kb


  nrj(ntrac,7)=0
  rg=dmod(ts,1.d0) ! time interpolation constant between 0 and 1
  rr=1.d0-rg
  if(rg.lt.0.d0 .or.rg.gt.1.d0) then
     print *,'rg=',rg
     goto 1500
  endif
 


  dxyz=dzt(ib,jb,kb)
  dxyz=dxyz*dxdy(ib,jb)


  dsmin=dtmin/dxyz
  ss0=dble(idint(ts))*tseas/dxyz
  call cross_time(1,ia,ja,ka,x0,dse,dsw,ts,tt,dsmin,dxyz,rr) ! zonal
  call cross_time(2,ia,ja,ka,y0,dsn,dss,ts,tt,dsmin,dxyz,rr) ! meridional
  call cross_time(3,ia,ja,ka,z0,dsu,dsd,ts,tt,dsmin,dxyz,rr) ! vertical



  ds=dmin1(dse,dsw,dsn,dss,dsu,dsd,dsmin)
  if(ds.eq.UNDEF .or.ds.eq.0.d0)then 
     ! === Can not find any path for unknown reasons ===
     print *,'ds cross error',ds,dse,dsw,dsn,dss,dsu,dsd,dsmin,dxyz
     print *,ia,ja,ka,x0,y0,z0,ntrac,niter
     print *,'k=',ka,kb,KM+1-kmt(ia,ja),kmt(ia,ja)
     ! goto 1500
     nerror=nerror+1
     nrj(ntrac,6)=1
     cycle ntracLoop
  endif
  

  if(ds.eq.dsmin) then ! transform ds to dt in seconds
     dt=dtmin  ! this makes dt more accurate
  else
     dt=ds*dxyz 
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
        ts =ts +dt/tseas
        tss=tss+dt/tseas*dble(iter)
     endif
  endif
  ! === time interpolation constant ===
  rbg=dmod(ts,1.d0) 
  rb =1.d0-rbg
  
  ! === calculate the new positions ===
  ! === of the trajectory           ===    
  if(ds.eq.dse) then ! eastward grid-cell exit 
     uu=(rbg*uflux(ia,ja,ka,NST)+rb*uflux(ia ,ja,ka,1))*ff
     if(uu.gt.0.d0) then
        ib=ia+1
        if(ib.gt.IMT) ib=ib-IMT 
     endif
     x1=dble(ia)
     call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds)
     call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds)
  elseif(ds.eq.dsw) then ! westward grid-cell exit
     uu=(rbg*uflux(iam,ja,ka,NST)+rb*uflux(iam,ja,ka,1))*ff
     if(uu.lt.0.d0) then
        ib=iam
     endif
     x1=dble(iam)
     call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds)
     call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds)
  elseif(ds.eq.dsn) then ! northward grid-cell exit              
     uu=(rbg*vflux(ia,ja,ka,NST)+rb*vflux(ia,ja,ka,1))*ff
     if(uu.gt.0.d0) then
        jb=ja+1
     endif
     y1=dble(ja)
     call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds)
     call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds)
  elseif(ds.eq.dss) then ! southward grid-cell exit
     uu=(rbg*vflux(ia,ja-1,ka,NST)+rb*vflux(ia,ja-1,ka,1))*ff
     if(uu.lt.0.d0) then
        jb=ja-1
     endif
     y1=dble(ja-1)
     call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds)
     call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds)
  elseif(ds.eq.dsu) then ! upward grid-cell exit
     uu=wflux(ia,ja,ka,1)
     if(uu.gt.0.d0) then
        kb=ka+1
     endif
     z1=dble(ka)
     call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds)
     call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds)
  elseif(ds.eq.dsd) then ! downward grid-cell exit
     call vertvel(rb,ia,iam,ja,ka)
     if(wflux(ia,ja,ka-1,1).lt.0.d0) kb=ka-1
     z1=dble(ka-1)
     call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds)
     call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds)
  elseif( ds.eq.dsc .or. ds.eq.dsmin) then  ! shortest time is the time-steping 
     ! If there is no spatial solution, which should correspond to a convergence zone
     if(dse.eq.UNDEF .and. dsw.eq.UNDEF .and. dsn.eq.UNDEF .and. & 
          dss.eq.UNDEF .and. dsu.eq.UNDEF .and. dsd.eq.UNDEF) then
        x1=x0 ; y1=y0 ; z1=z0 ! let the particle remain in a 
                              !static position from previuos iteration
        ib=ia ; jb=ja ; kb=ka  
        ! If there is at least one spatial solution but the shortest cross time is the time step
     else
        call pos_time(1,ia,ja,ka,x0,x1,ts,tt,dsmin,dxyz,ss0,ds)
        call pos_time(2,ia,ja,ka,y0,y1,ts,tt,dsmin,dxyz,ss0,ds)
        call pos_time(3,ia,ja,ka,z0,z1,ts,tt,dsmin,dxyz,ss0,ds)
     endif
  endif

end MODULE tracmass
