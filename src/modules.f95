

MODULE mod_precdef		! Precision definitions
  INTEGER, PARAMETER		            :: DP = SELECTED_REAL_KIND(15, 307)
ENDMODULE mod_precdef


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===



MODULE mod_param
  INTEGER                                   :: JMAX, LBT, NTRACMAX
  INTEGER, PARAMETER                        :: MR=501 ! or 1001
#ifdef streamts
  INTEGER, PARAMETER                        :: LOV=3
#else
  INTEGER, PARAMETER                        :: LOV=1
#endif
  INTEGER                                   :: ncoor,kriva,iter,ngcm,ntrac
  INTEGER                                   :: twritetype = 0
  REAL*8                                    :: tseas,tyear,dtmin,voltr
  REAL*8                                    :: tstep,dstep,tss,partQuant
  REAL*8, PARAMETER                         :: UNDEF=1.d20 

  REAL*8, PARAMETER                         :: grav = 9.81
  REAL*8, PARAMETER                         :: PI = 3.14159265358979323846d0
  REAL*8, PARAMETER                         :: radius = 6371229.d0 
  REAL*8, PARAMETER                         :: radian = pi/180.d0  
  REAL*8, PARAMETER                         :: deg=radius*radian   
  REAL*8, PARAMETER                         :: tday=24.d0 * 3600.d0
  INTEGER                                   :: idmax(12,1000:3000)
ENDMODULE mod_param


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

MODULE mod_coord

ENDMODULE mod_coord

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_diff
	INTEGER                             :: dummy	
ENDMODULE mod_diff

MODULE mod_loopvars
  REAL*8                                     :: rr, rb, rg, rbg
  REAL*8                                     :: ds, dsmin
  REAL*8                                     :: dse, dsw, dsn, dss
  REAL*8                                     :: dsu, dsd, dsc
  LOGICAL                                    :: scrivi
  REAL*8                                     :: ts,tt
  REAL*8                                     :: dxyz
  REAL*8                                     :: ss0
  INTEGER                                    :: lbas
  REAL*8                                     :: subvol
ENDMODULE mod_loopvars

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_time
  INTEGER                                   :: ints      ,intstart ,intend
  INTEGER                                   :: intrun    ,intspin  ,intstep
  INTEGER                                   :: intmin    ,intmax
  ! === Calculate julian dates
  INTEGER                                   :: baseYear  ,baseMon  ,baseDay
  INTEGER                                   :: baseHour  ,baseMin  ,baseSec
  INTEGER                                   :: startYear ,startMon ,startDay
  INTEGER                                   :: startHour ,startMin ,startSec
  REAL*8                                    :: currJDtot ,currJDyr,currfrac
  INTEGER                                   :: currYear  ,currMon  ,currDay
  INTEGER                                   :: currHour, currMin, currSec 

  INTEGER                                   :: iyear ,imon ,iday ,ihour
  INTEGER                                   :: iyear0 ,imon0 ,iday0 
  INTEGER                                   :: yearmin ,yearmax

  INTEGER*8                                 :: ntime
 
  REAL*8                                    :: startJD=0 ,baseJD=0,ttpart
  INTEGER                                   :: fieldsPerFile
  ! === Time-interpolation variables in loop ===
  REAL*8                                     :: dt, t0
  REAL*8                                     :: dtreg

CONTAINS

  subroutine updateClock  
    USE mod_param
    USE mod_loopvars
    ttpart = anint((anint(tt)/tseas-floor(anint(tt)/tseas))*tseas)/tseas 
    currJDtot = (ints+ttpart-1)*(real(ngcm)/24) 
    call  gdate (baseJD+currJDtot-1 ,currYear , currMon ,currDay)
    currJDyr = baseJD+currJDtot - jdate(currYear ,1 ,1)
    currFrac = (currJDtot-int(currJDtot))*24
    currHour = int(currFrac)
    currFrac = (currFrac - currHour) * 60
    CurrMin  = int(currFrac)
    currSec  = int((currFrac - currMin) * 60)
  end subroutine updateClock
  
  subroutine gdate (rjd, year,month,day)
    !                                                                      
    !---computes the gregorian calendar date (year,month,day)              
    !   given the julian date (jd).                                        
    !   Source: http://aa.usno.navy.mil/faq/docs/JD_Formula.php            
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

  INTEGER function jdate (YEAR,MONTH,DAY)
    !---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
    !   DATE (YEAR,MONTH,DAY).
    !   Source: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    INTEGER                                  :: YEAR,MONTH,DAY,I,J,K
    i     = year
    j     = month
    k     = day
    jdate = K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12) &
         /12-3*((I+4900+(J-14)/12)/100)/4
    RETURN
  end function jdate
ENDMODULE mod_time


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_grid
  INTEGER                                   :: IMT, JMT, KM
  INTEGER, PARAMETER                        :: NST=2
  INTEGER                                   :: nsm=1     ,nsp=2
  REAL*8                                    :: dx,dy
  REAL*8                                    :: dxdeg,dydeg,stlon1,stlat1
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)     :: botbox
  REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu, ang
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: mask
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: csu,cst,dyt,phi

  ! === Vertical grids ===
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: zw
#ifdef zgrid3Dt 
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)     :: dzt
#elif zgrid3D
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt,dzu,dzv
  REAL, ALLOCATABLE, DIMENSION(:,:)         :: dzt0surf,dzu0surf,dzv0surf
#endif /*zgrid3Dt*/
#ifdef varbottombox 
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dztb
#endif /*varbottombox*/
#ifdef ifs
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: aa,bb
#endif
  REAL*8                                    :: rmin ,tmin ,smin,&
  &                                            rmax ,smax ,tmax
  REAL*8                                    :: dr ,dtemp ,dsalt
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt, kmu, kmv, depth
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  CHARACTER(LEN=200)                        :: SubGridFile 
  INTEGER                                   :: degrade_time=0, degrade_space=0

#ifdef ifs
  REAL*8, PARAMETER                         :: R_d = 287.05d0
  REAL*8, PARAMETER                         :: L_v = 2.5d0 * 1e+6   
  REAL*8, PARAMETER                         :: c_d = 1004.d0
#endif

CONTAINS
  function l2d(lon1,lon2,lat1,lat2)
    real                                   :: lon1,lon2,lat1,lat2,l2d
    real                                   :: dlon,dlat,a,c
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * asin(min(1.0,sqrt(a)))
    l2d = 6367 * c * 1000
  end function l2d

ENDMODULE mod_grid
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_buoyancy
  REAL*4                                    :: tmin0 ,tmax0
  REAL*4                                    :: smin0 ,smax0
  REAL*4                                    :: rmin0 ,rmax0
  REAL*4                                    :: tmine ,tmaxe
  REAL*4                                    :: smine ,smaxe
  REAL*4                                    :: rmine ,rmaxe
ENDMODULE mod_buoyancy
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_domain
  INTEGER, DIMENSION(10)                    :: ienw ,iene
  INTEGER, DIMENSION(10)                    :: jens ,jenn
  REAL*4                                    :: timax
ENDMODULE mod_domain
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_vel
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: uflux ,vflux
#if defined full_wflux
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)    :: wflux
#else
  REAL*8, ALLOCATABLE, DIMENSION(:,:)        :: wflux
#endif
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)      :: uvel ,vvel ,wvel 
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: hs
  REAL*8                                     :: ff

CONTAINS
  
!  subroutine datasetswap
!    hs(:,:,1)      = hs(:,:,2)
!    uflux(:,:,:,1) = uflux(:,:,:,2)
!    vflux(:,:,:,1) = vflux(:,:,:,2)
!#ifdef explicit_w || full_wflux
!    wflux(:,:,:,1) = wflux(:,:,:,2)
!#endif
!#ifdef tempsalt
!    tem(:,:,:,1)   = tem(:,:,:,2)
!    sal(:,:,:,1)   = sal(:,:,:,2)
!    rho(:,:,:,1)   = rho(:,:,:,2)
!#endif
!  end subroutine datasetswap

  subroutine calc_implicit_vertvel
    USE mod_grid
    IMPLICIT none
    ! = Loop variables
    INTEGER                                    :: k
    
    wflux(2:imt,2:jmt,1,2)    =  uflux(1:imt-1, 2:jmt,   1,   2)  -   &
                                 uflux(2:imt,   2:jmt,   1,   2)   +  &
                                 vflux(2:imt,   1:jmt-2, 1,   2)  -   & 
                                 vflux(2:imt,   2:jmt,   1,   2)
    wflux(1, 2:jmt, 1, 2)     =  uflux(1,       2:jmt,   1,   2)   +  &
                                 vflux(1,       1:jmt-1, 1,   2)  -   &
                                 vflux(1,       2:jmt,   1,   2)
    wflux(2:imt, 1, 1, 2)   =    uflux(1:imt-1, 1,       1,   2) -    &
                                 uflux(2:imt,   1,       1,   2)  -  &
                                 vflux(2:imt,   1,       1,   2) 
    kloop: do k=2,km
       wflux(2:imt,2:jmt,k,2) =  wflux(2:imt,   2:jmt,   k-1, 2)   +  &
                                 uflux(1:imt-1, 2:jmt,   k,   2)  -   &
                                 uflux(2:imt,   2:jmt,   k,   2)   +  &
                                 vflux(2:imt,   1:jmt-1, k,   2)  -   & 
                                 vflux(2:imt,   2:jmt,   k,   2) 
       wflux(1,2:jmt,k,2)     =  wflux(1,       2:jmt,   k-1, 2)  -   &
                                 uflux(1,       2:jmt,   k,   2)   +  &
                                 vflux(1,       1:jmt-1, k,   2)  -   &
                                 vflux(1,       2:jmt,   k,   2)
       wflux(2:imt,1,  k,2)   =  wflux(2:imt,   1,       k-1, 2)   +  &
                                 uflux(1:imt-1, 1,       k,   2) -    &
                                 uflux(2:imt,   1,       k,   2)  -  &
                                 vflux(2:imt,   1,       k,   2) 
    enddo kloop
  end subroutine calc_implicit_vertvel

ENDMODULE mod_vel
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_traj
  INTEGER, PARAMETER                         :: NNRJ=8,NTRJ=7
  INTEGER                                    :: NEND

  ! === Particle arrays ===
  REAL*8, ALLOCATABLE, DIMENSION(:,:)        :: trj
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: nrj 

  ! === Particle counters ===
  INTEGER                                    :: nout=0, nloop=0, nerror=0
  INTEGER                                    :: nnorth=0, ndrake=0, ngyre=0
  INTEGER, ALLOCATABLE,DIMENSION(:)          :: nexit

  ! === Particle positions ===
  INTEGER                                    :: ia, ja, ka, iam
  INTEGER                                    :: ib, jb, kb, ibm
  REAL*8                                     :: x0, y0, z0
  REAL*8                                     :: x1, y1, z1
ENDMODULE mod_traj


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


MODULE mod_dens
#ifdef tempsalt
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: tem,sal,rho
#endif
ENDMODULE mod_dens
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_turb
#ifdef turb
  REAL upr(12,2)
#endif
ENDMODULE mod_turb
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_name
  CHARACTER(LEN=200)                         :: outDataFile
  INTEGER                                    :: intminInOutFile
  CHARACTER(LEN=200)                         :: inDataDir ,outDataDir, topoDataDir
  CHARACTER(LEN=200)                         :: projDesc
  CHARACTER(LEN=200)                         :: GCMname   ,GCMsource
  CHARACTER(LEN=200)                         :: gridName  ,gridSource
  CHARACTER(LEN=200)                         :: gridDesc
  CHARACTER(LEN=200)                         :: caseName  ,caseDesc
ENDMODULE mod_name

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
ENDMODULE mod_streamfunctions

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_tracer
#ifdef tracer
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: tra
#endif
ENDMODULE mod_tracer

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
#if defined diffusion || turb 
!#ifdef diffusion
MODULE mod_diffusion
  REAL                                       :: ah, av
ENDMODULE mod_diffusion
#endif

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_sed
#ifdef sediment
  !  REAL :: wsed,rhos,D,critvel,T,cwamp,kincrit
  REAL                                       :: wsed, partdiam
  REAL                                       :: rhos, cwamp, twave
  REAL                                       :: critvel, kincrit

  INTEGER                                    :: nsed=0, nsusp=0
  LOGICAL                                    :: res
#endif
ENDMODULE mod_sed
MODULE mod_orbital
#ifdef sediment
  REAL, ALLOCATABLE, DIMENSION(:)            :: orb
#endif
ENDMODULE mod_orbital

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
