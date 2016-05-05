

MODULE mod_precdef		! Precision definitions
  INTEGER, PARAMETER		            :: QP = selected_real_kind (32)
!  INTEGER, PARAMETER		            :: DP = QP
  INTEGER, PARAMETER		            :: DP = SELECTED_REAL_KIND(15, 307)
ENDMODULE mod_precdef


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_param
  INTEGER                                   :: jmax, ntracmax
  INTEGER, PARAMETER                        :: MR=501 ! or 1001
  INTEGER                                   :: ncoor,kriva,iter,ngcm
  REAL*8, PARAMETER                         :: UNDEF=1.d20 
  REAL*8, PARAMETER                         :: EPS=1.d-7 ! the small epsilon

  REAL*8, PARAMETER                         :: grav = 9.81
  REAL*8, PARAMETER                         :: PI = 3.14159265358979323846d0
  REAL*8, PARAMETER                         :: radius = 6371229.d0 
  REAL*8, PARAMETER                         :: radian = pi/180.d0  
  REAL*8, PARAMETER                         :: deg=radius*radian   
  REAL*8, PARAMETER                         :: tday=24.d0 * 3600.d0
  INTEGER                                   :: idmax(12,1000:3000)
ENDMODULE mod_param
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_loopvars
  REAL*8                                     :: ds, dsmin
  REAL*8                                     :: dse, dsw, dsn, dss
  REAL*8                                     :: dsu, dsd, dsc
  LOGICAL                                    :: scrivi
  INTEGER                                    :: niter
  REAL*8                                     :: ss0
  INTEGER                                    :: lbas
  REAL*8                                     :: subvol
ENDMODULE mod_loopvars
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

MODULE mod_traj

  ! Variables connected to particle positions.
  INTEGER, PARAMETER                         :: NNRJ=8, NTRJ=7
  INTEGER                                    :: nend
  INTEGER                                    :: ntrac, ntractot=0
  ! === Particle arrays ===
  REAL*8, ALLOCATABLE,  DIMENSION(:,:)       :: trj
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: nrj 
  ! === Particle counters ===
  INTEGER                                    :: nout=0, nloop=0, nerror=0, nrh0=0
  INTEGER, ALLOCATABLE,DIMENSION(:)          :: nexit
  ! === Particle positions ===
  INTEGER                                    :: ia, ja, ka, iam
  INTEGER                                    :: ib, jb, kb, ibm
  REAL*8                                     :: x0, y0, z0
  REAL*8                                     :: x1, y1, z1
ENDMODULE mod_traj
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_tempsalt
  REAL*8                                     :: rmin, tmin, smin
  REAL*8                                     :: rmax, tmax, smax
  REAL*8                                     :: dr ,dtemp ,dsalt
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: tem,sal,rho
  REAL*4                                     :: tmin0 ,tmax0
  REAL*4                                     :: smin0 ,smax0
  REAL*4                                     :: rmin0 ,rmax0
  REAL*4                                     :: tmine ,tmaxe
  REAL*4                                     :: smine ,smaxe
  REAL*4                                     :: rmine ,rmaxe
end MODULE mod_tempsalt

  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_grid
  USE mod_param, only: pi, undef, iter
  IMPLICIT NONE

  INTEGER                                   :: imt, jmt, km
  INTEGER                                   :: nst=2
  INTEGER                                   :: nsm=1,  nsp=2
  INTEGER                                   :: wnsm=1, wnsp=2
  REAL*8                                    :: dx,dy
  REAL*8                                    :: dxdeg,dydeg,stlon1,stlat1
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)     :: hs
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)     :: botbox
  REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu, ang
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  REAL*8                                    :: dxyz
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: mask
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: csu,cst,dyt,phi

  ! === Vertical grids ===
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: zlev
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)   :: z_r, z_w
#if defined zgrid3Dt 
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)     :: dzt
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dztb,dzu,dzv
  REAL, ALLOCATABLE, DIMENSION(:,:)         :: abyst,abysu,abysv
#elif zgrid3D
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt, dzu, dzv
  REAL, ALLOCATABLE, DIMENSION(:,:)         :: dzt0surf,dzu0surf,dzv0surf
#endif /*zgrid3Dt*/
#ifdef varbottombox 
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dztb
#endif /*varbottombox*/
#ifdef ifs
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: aa, bb
#endif
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt, kmu, kmv, depth
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  INTEGER                                   :: subGridKmin=1 ,subGridKmax=0
  CHARACTER(LEN=200)                        :: SubGridFile 
  INTEGER                                   :: degrade_space=0

#ifdef ifs
  REAL*8, PARAMETER                         :: R_d = 287.05d0
  REAL*8, PARAMETER                         :: L_v = 2.5d0 * 1e+6   
  REAL*8, PARAMETER                         :: c_d = 1004.d0
#endif

CONTAINS
  function l2d(lon1,lon2,lat1,lat2)
    real                                   :: lon1,lon2,lat1,lat2,l2d
    real                                   :: rlon1,rlon2,rlat1,rlat2
    real                                   :: dlon,dlat,a,c
    dlon = (lon2 - lon1)/180*pi
    rlat1 = lat1 /180.*pi
    rlat2 = lat2 /180.*pi
    dlat = rlat2 - rlat1
    a = (sin(dlat/2))**2 + cos(rlat1) * cos(rlat2) * (sin(dlon/2))**2
    c = 2 * asin(min(1.0,sqrt(a)))
    l2d = 6367 * c * 1000
  end function l2d

  subroutine calc_dxyz(intrpr, intrpg)

    use mod_traj, only: ib,jb,kb
    IMPLICIT NONE
    REAL*8                                    :: intrpr, intrpg

    ! T-box volume in m3
#ifdef zgrid3Dt 
    dxyz = intrpg * dzt(ib,jb,kb,nsp) + intrpr * dzt(ib,jb,kb,nsm)
#elif  zgrid3D
    dxyz = dzt(ib, jb, kb)
#ifdef freesurface
    if(kb == KM) dxyz = dxyz + intrpg * hs(ib,jb,nsp) + intrpr * hs(ib,jb,nsm)
#endif /*freesurface*/
#else
    dxyz =dz(kb)
#ifdef varbottombox
    if(kb == KM+1-kmt(ib,jb) ) dxyz=dztb(ib,jb,1)
#endif /*varbottombox*/
#ifdef freesurface
    if(kb == KM) dxyz=dxyz+intrpg*hs(ib,jb,nsp)+intrpr*hs(ib,jb,nsm)
#endif /*freesurface*/
#endif /*zgrid3Dt*/
    dxyz=dxyz*dxdy(ib,jb)
    if (dxyz<0) then
       print *,'=========================================================='
       print *,'ERROR: Negative box volume                                '
       print *,'----------------------------------------------------------'
       print *,'dxdy = ', dxdy(ib,jb)
       print *,'ib  = ', ib, ' jb  = ', jb, ' kb  = ', kb 
       print *,'----------------------------------------------------------'
       print *,'The run is terminated'
       print *,'=========================================================='
       stop
    end if
  end subroutine calc_dxyz

ENDMODULE mod_grid
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_time
  ! Variables and routines for timekeeping

  !Timestep increasing with one for each new velocity field
  INTEGER                                   :: ints      ,intstart ,intend
  INTEGER                                   :: intrun    ,intspin
  INTEGER                                   :: intmin    ,intmax
  !type for datetimes
  type DATETIME
     REAL*8                                 :: JD=0 
     REAL                                   :: frac=0, yd=0
     INTEGER                                :: Year, Mon, Day
     INTEGER                                :: Hour, Min, Sec
  end type DATETIME
  type(DATETIME), SAVE                      :: basetime, starttime
  type(DATETIME), SAVE                      :: currtime, looptime
  ! === Base for JD (When JD is 1)
  REAL*8                                    :: baseJD=0
  INTEGER                                   :: baseYear  ,baseMon  ,baseDay
  INTEGER                                   :: baseHour  ,baseMin  ,baseSec
  REAL*8                                    :: jdoffset=0
  LOGICAL                                   :: noleap=.false.
  ! === Timerange for velocity fields
  REAL*8                                    :: minvelJD=0,   maxvelJD=0
  INTEGER                                   :: minvelints, maxvelints
  ! === JD when the run starts
  REAL*8                                    :: startJD=-999, ttpart, startFrac
  INTEGER                                   :: startYear, startMon, startDay
  INTEGER                                   :: startHour, startMin, startSec
  REAL*8                                    :: endJD=-999, endFrac
  INTEGER                                   :: endYear=0, endMon,   endDay
  INTEGER                                   :: endHour,   endMin,   endSec
  ! === Current JD
  REAL*8                                    :: currJDtot ,currJDyr,currfrac
  INTEGER                                   :: currYear  ,currMon  ,currDay
  INTEGER                                   :: currHour, currMin, currSec
  INTEGER                                   :: leapoffset=0
  ! === Looping time
  INTEGER                                   :: loopints, loopintstart
  REAL*8                                    :: loopJD, loopJDyr, loopFrac
  INTEGER                                   :: loopYear  ,loopMon  ,loopDay
  INTEGER                                   :: loopHour, loopMin, loopSec 
  ! Old stuff
  INTEGER                                   :: iyear ,imon ,iday ,ihour
  INTEGER                                   :: yearmin ,yearmax

  INTEGER*8                                 :: ntime
  ! Used to figure out when to change file.
  INTEGER                                   :: fieldsPerFile
  ! === Time-interpolation variables in loop ===
  REAL*8                                    :: dt, t0
  REAL*8                                    :: dtreg
  REAL*8                                    :: tseas, tyear, dtmin,voltr
  REAL*8                                    :: tstep, dstep, tss, partQuant
  REAL*8                                    :: ts, tt
  REAL*8                                    :: intrpr, intrpg
  REAL*8                                    :: intrpb, intrpbg
CONTAINS

  subroutine updateClock  
    USE mod_param, only: ngcm
    IMPLICIT NONE
    ttpart = anint((anint(tt,8)/tseas-floor(anint(tt,8)/tseas))*tseas)/tseas 
    !currJDtot = (ints+ttpart)*(dble(ngcm)/24.) !! LD: included startMin and startSec; irrelevant if both are 0; there may be a better way to do this...
    currJDtot = (ints+ttpart)*(dble(ngcm)/24.) + startMin/(60*24.) + startSec/(60*60*24.)
    call  gdate (baseJD+currJDtot-1+jdoffset + leapoffset,  &
                 currYear , currMon ,currDay)
    currJDyr = baseJD + currJDtot - jdate(currYear ,1 ,1) + jdoffset
    if ((mod(currYear, 4) == 0)  .and. (currJDyr>56) .and.     &
         (currJDyr<(56 - leapoffset + ngcm/24.)) .and. noleap) then
       leapoffset = leapoffset + 1
       call  gdate (baseJD+currJDtot-1+jdoffset + leapoffset,  &
            currYear , currMon ,currDay)
    end if
    currJDyr = baseJD + currJDtot - jdate(currYear ,1 ,1) + jdoffset
    currFrac = (currJDtot-dble(int(currJDtot,8)))*24
    currHour = int(currFrac,8)
    currFrac = (currFrac - dble(currHour)) * 60
    CurrMin  = int(currFrac,8)
    currSec  = int((currFrac - dble(currMin)) * 60,8)

    if (ints > (maxvelints-1)) then
       if (minvelints == 0) then
          loopints = ints - intmax * int(real(ints-intstart)/intmax)
       else
          loopints = ints - intmax * int(real(ints-minvelints)/intmax)
          intmax = maxvelints - minvelints
       end if
    else
       loopints = ints
    end if
    !loopJD = (loopints + ttpart)*(dble(ngcm)/24) + 1 !! LD: removed +1 so consistent with interpolation bounds; added startMin and startSec
    loopJD = (loopints + ttpart)*(dble(ngcm)/24) + startHour/(24.) + startMin/(60*24.) + startSec/(60*60*24.)
    call  gdate (baseJD+loopJD-1+jdoffset ,loopYear, loopMon, loopDay)
    loopJDyr = baseJD+loopJD - jdate(loopYear ,1 ,1)
    loopFrac = (loopJD - dble(int(loopJD,8))) * 24
    loopHour = int(loopFrac,8)
    loopFrac = (loopFrac - dble(loopHour)) * 60
    LoopMin  = int(loopFrac,8)
    loopSec  = int((loopFrac - dble(loopMin)) * 60,8)
  end subroutine updateClock
  
  subroutine gdate (rjd, year,month,day)
    !Computes the gregorian calendar date given a julian date (jd).
    !Source: http://aa.usno.navy.mil/faq/docs/JD_Formula.php            
    REAL*8                                   :: rjd
    INTEGER                                  :: jd
    INTEGER                                  :: year ,month ,day
    INTEGER                                  :: i ,j ,k ,l ,n
    
    jd = int(rjd)
    l= jd+68569
    n= 4*l/146097
    l= l-(146097*n+3)/4
    i= 4000*(l+1)/1461001
    l= l-1461*i/4+31
    j= 80*l/2447
    k= l-2447*j/80
    l= j/11
    j= j+2-12*l
    i= 100*(n-49)+i+l
    
    year= i
    month= j
    day= k
    return
  end subroutine gdate

  INTEGER function jd2ints(jd)
    USE mod_param, only: ngcm
    REAL*8                                    :: jd
    jd2ints = int(floor((jd)/(real(ngcm)/24.))) 
    return
  end function jd2ints
  
  INTEGER function jdate (year, month, day)
    !Computes the julian date (JD) given a gregorian calendar date.
    !Source: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    INTEGER                                  :: year, month ,day
    INTEGER                                  :: i, j, k
    i     = year
    j     = month
    k     = day
    jdate = K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12) &
         /12-3*((I+4900+(J-14)/12)/100)/4
    RETURN
  end function jdate

  subroutine calc_time
    USE mod_loopvars, only: ds, dsc, dsmin
    use mod_grid, only: dxyz
    USE mod_param, only: iter
    IMPLICIT NONE

#ifdef regulardt
    if(ds == dsmin) then ! transform ds to dt in seconds
       !            dt=dt  ! this makes dt more accurate
    else
       dt = ds * dxyz 
    endif
#else
    if(ds == dsmin) then ! transform ds to dt in seconds
       dt=dtmin  ! this makes dt more accurate
    else
       dt = ds * dxyz 
    endif
#endif /*regulardt*/
    if(dt.lt.0.d0) then
       Print *,"Error! dt is less than zero."
       print *,'dt=',dt,"ds=",ds,"dxyz=",dxyz,"dsmin=",dsmin
       stop 4968
    endif
    ! === if time step makes the integration ===
    ! === exceed the time when fields change ===
    if(tss+dt/tseas*dble(iter).ge.dble(iter)) then
       dt=dble(int(ts,8)+1)*tseas-tt
       tt=dble(int(ts,8)+1)*tseas
       ts=dble(int(ts,8)+1)
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
          ts=nint((ts+dtreg/tseas)*dble(iter), 8)/dble(iter)
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
#endif /*regulardt*/
    end if
    ! === time interpolation constant ===
    intrpbg=dmod(ts,1.d0) 
    intrpb =1.d0-intrpbg
  end subroutine calc_time

ENDMODULE mod_time
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_buoyancy
ENDMODULE mod_buoyancy
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_domain
  INTEGER, DIMENSION(10)                :: ienw ,iene, jens ,jenn
  REAL*4                                :: timax
ENDMODULE mod_domain
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===        
MODULE mod_dens

ENDMODULE mod_dens
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===        


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_vel
  USE mod_grid, only: nsm, nsp, dzt
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: uflux, vflux
#if defined explicit_w || full_wflux
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)    :: wflux
#else
  REAL*8, ALLOCATABLE, DIMENSION(:,:)        :: wflux
#endif
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)      :: uvel ,vvel ,wvel 
  REAL*8                                     :: ff
  INTEGER                                    :: degrade_time=0
    integer, save                            :: degrade_counter = 0

CONTAINS
 
  subroutine datasetswap

    USE  mod_grid, only      : nsm,nsp,hs
    USE  mod_tempsalt, only  : tem,sal,rho

    IMPLICIT NONE

    hs(:,:,nsm)      = hs(:,:,nsp)
    uflux(:,:,:,nsm) = uflux(:,:,:,nsp)
    vflux(:,:,:,nsm) = vflux(:,:,:,nsp)
!    wflux(:,nsm) = wflux(:,nsp) 
#if defined zgrid3Dt 
    dzt(:,:,:,nsm)   = dzt(:,:,:,nsp)
#endif
#if defined explicit_w || full_wflux
    wflux(:,:,:,nsm) = wflux(:,:,:,nsp)
#endif
#ifdef tempsalt
    tem(:,:,:,nsm)   = tem(:,:,:,nsp)
    sal(:,:,:,nsm)   = sal(:,:,:,nsp)
    rho(:,:,:,nsm)   = rho(:,:,:,nsp)
#endif
  end subroutine datasetswap

#if defined full_wflux
  subroutine calc_implicit_vertvel !This sub needs to be updated to include dzt stuff!!!!

    USE mod_grid, only: imt, jmt, km, nsm, nsp
    IMPLICIT none
    INTEGER                                    :: k
    
    wflux(2:imt,2:jmt,1,nsp)    =  uflux(1:imt-1, 2:jmt,   1,   nsp)  -  &
                                   uflux(2:imt,   2:jmt,   1,   nsp)  +  &
                                   vflux(2:imt,   1:jmt-2, 1,   nsp)  -  & 
                                   vflux(2:imt,   2:jmt,   1,   nsp)
    wflux(1, 2:jmt,  1, nsp)    =  uflux(1,       2:jmt,   1,   nsp)  +  &
                                   vflux(1,       1:jmt-1, 1,   nsp)  -  &
                                   vflux(1,       2:jmt,   1,   nsp)
    wflux(2:imt, 1,  1, nsp)    =  uflux(1:imt-1, 1,       1,   nsp)  -  &
                                   uflux(2:imt,   1,       1,   nsp)  -  &
                                   vflux(2:imt,   1,       1,   nsp) 
    kloop: do k=2,km
       wflux(2:imt,2:jmt,k,nsp) =  wflux(2:imt,   2:jmt,   k-1, nsp)  +  &
                                   uflux(1:imt-1, 2:jmt,   k,   nsp)  -  &
                                   uflux(2:imt,   2:jmt,   k,   nsp)  +  &
                                   vflux(2:imt,   1:jmt-1, k,   nsp)  -  & 
                                   vflux(2:imt,   2:jmt,   k,   nsp) 
       wflux(1,2:jmt,k,nsp)     =  wflux(1,       2:jmt,   k-1, nsp)  -  &
                                   uflux(1,       2:jmt,   k,   nsp)  +  &
                                   vflux(1,       1:jmt-1, k,   nsp)  -  &
                                   vflux(1,       2:jmt,   k,   nsp)
       wflux(2:imt,1,  k, nsp)  =  wflux(2:imt,   1,       k-1, nsp)  +  &
                                   uflux(1:imt-1, 1,       k,   nsp)  -  &
                                   uflux(2:imt,   1,       k,   nsp)  -  &
                                   vflux(2:imt,   1,       k,   nsp) 
    enddo kloop
  end subroutine calc_implicit_vertvel
#endif /*full_wflux*/


ENDMODULE mod_vel
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_name
  CHARACTER(LEN=200)                         :: inDataDir
  CHARACTER(LEN=200)                         :: topoDataDir
  CHARACTER(LEN=200)                         :: projDesc
  CHARACTER(LEN=200)                         :: GCMname   ,GCMsource
  CHARACTER(LEN=200)                         :: gridName  ,gridSource
  CHARACTER(LEN=200)                         :: gridDesc
  CHARACTER (LEN=23)                         :: Project, Case
  CHARACTER(LEN=200)                         :: caseName  ,caseDesc
ENDMODULE mod_name
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_streamfunctions
#ifdef streamxy
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: stxyy, stxyx
#endif
#ifdef streamv
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: stxz, styz
#endif
#ifdef streamr
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)      :: stxr,styr
#endif
#ifdef stream_thermohaline
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)      :: psi_ts
#endif
  INTEGER                        :: INTPSI ! to be read by the xxx.in files in future
#ifdef streamts
  INTEGER, PARAMETER                        :: LOV=3
#else
  INTEGER, PARAMETER                        :: LOV=1
#endif
ENDMODULE mod_streamfunctions
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_tracer
#ifdef tracer
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: tra
#endif
ENDMODULE mod_tracer
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
#if defined diffusion || turb 
!#ifdef diffusion
MODULE mod_diffusion
  REAL                                       :: ah, av
ENDMODULE mod_diffusion
#endif
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_sed
#ifdef sediment
  !  REAL :: wsed,rhos,D,critvel,T,cwamp,kincrit
  REAL                                       :: wsed, partdiam
  REAL                                       :: rhos, cwamp, twave
  REAL                                       :: critvel, kincrit

  INTEGER                                    :: nsed=0, nsusp=0
  LOGICAL                                    :: res
  INTEGER                                    :: numseedsubtimes=0
  REAL, ALLOCATABLE, DIMENSION(:)            :: seedsubtimes


#endif
ENDMODULE mod_sed
MODULE mod_orbital
#ifdef sediment
  REAL, ALLOCATABLE, DIMENSION(:)            :: orb
#endif
ENDMODULE mod_orbital
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
