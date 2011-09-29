
MODULE mod_param
  INTEGER                                   :: IMT, JMT, KM
  INTEGER                                   :: JMAX, LBT, NTRACMAX
  INTEGER, PARAMETER                        :: MR=1001
  INTEGER                                   :: NEND
  INTEGER, PARAMETER                        :: NST=2,NNRJ=8,NTRJ=7
#ifdef streamts
  INTEGER, PARAMETER                        :: LOV=3
#else
  INTEGER, PARAMETER                        :: LOV=1
#endif
  INTEGER                                   :: ncoor,kriva,iter,ngcm
  INTEGER                                   :: twritetype = 0
  REAL*8                                    :: tseas,tday,tyear,dtmin,voltr
  REAL*8                                    :: tstep,dstep,tss,partQuant
  REAL*8, PARAMETER                         :: UNDEF=1.d20 
  REAL*8, PARAMETER                         :: PI = 3.14159265358979323846d0
ENDMODULE mod_param

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_precdef		! Precision definitions
  INTEGER, PARAMETER		         :: DP = SELECTED_REAL_KIND(15, 307)
ENDMODULE mod_precdef

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_coord
  REAL*8                                    :: dx,dy
  REAL*8                                    :: deg,stlon1,stlat1
  REAL*8, PARAMETER                         :: grav=9.81
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: zw
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: csu,cst,dyt,phi
  INTEGER idmax(12,1000:3000)
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
  INTEGER                                   :: baseYear  ,baseMon  ,baseDay
  INTEGER                                   :: baseHour  ,baseMin  ,baseSec
  INTEGER                                   :: startYear ,startMon ,startDay
  INTEGER                                   :: startHour ,startMin ,startSec
  INTEGER                                   :: iyear     ,imon     ,iday ,ihour
  INTEGER                                   :: iyear0    ,imon0    ,iday0 
  INTEGER                                   :: yearmin   ,yearmax
  INTEGER*8                                 :: ntime
  INTEGER                                   :: currJDtot ,currJDyr
  INTEGER                                   :: currYear  ,currMon  ,currDay
  INTEGER                                   :: currHour, currMin, currSec
  INTEGER                                   :: startJD   ,baseJD
  INTEGER                                   :: fieldsPerFile
CONTAINS
  subroutine updateClock  
    USE mod_param
    currJDtot = int((ints-1)*(real(ngcm)/24)+1)
    call  gdate (baseJD+currJDtot-1 ,currYear , currMon ,currDay)
    currJDyr = baseJD+currJDtot - jdate(currYear ,1 ,1)
    currHour = int(((ints-intmin)/(24./ngcm)-floor(real(ints-intmin)/ &
         (24./ngcm)))*24)

  end subroutine updateClock
  
  subroutine gdate (jd, year,month,day)
    !                                                                      
    !---computes the gregorian calendar date (year,month,day)              
    !   given the julian date (jd).                                        
    !   Source: http://aa.usno.navy.mil/faq/docs/JD_Formula.php            
    INTEGER                                  :: jd ,year ,month ,day
    INTEGER                                  :: i ,j ,k ,l ,n
    
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
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)     :: botbox
  REAL*4, ALLOCATABLE, DIMENSION(:,:)       :: dxv, dyu, ang
  REAL*8, ALLOCATABLE, DIMENSION(:)         :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: dxdy
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: mask
#ifdef zgrid3Dt 
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)     :: dzt
#elif zgrid3D
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt,dzu,dzv
  REAL, ALLOCATABLE, DIMENSION(:,:)         :: dzt0surf,dzu0surf,dzv0surf
#endif /*zgrid3Dt*/
#ifdef varbottombox 
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dztb
#endif /*varbottombox*/


  REAL*8                                    :: rmin ,tmin ,smin
  REAL*8                                    :: dr ,dtemp ,dsalt
  REAL*8                                    :: arc,arct,arcscale
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt, kmu,kmv, depth
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  CHARACTER(LEN=200)                        :: SubGridFile 
ENDMODULE mod_grid
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
MODULE mod_domain
  INTEGER, DIMENSION(10)                    :: ienw ,iene
  INTEGER, DIMENSION(10)                    :: jens ,jenn
  REAL*4                                    :: timax
ENDMODULE mod_domain
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
ENDMODULE mod_vel
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_traj
  REAL*8, ALLOCATABLE, DIMENSION(:,:)        :: trj
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: nrj 
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
MODULE mod_streamxy
#ifdef streamxy
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: stxyy, stxyx
  REAL                                       :: sxyy(IMT,JMT),sxyx(IMT,JMT)
#endif
ENDMODULE mod_streamxy

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_streamv
#ifdef streamv
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: stxz, styz
  REAL                                       :: sxz(JMT,KM),syz(JMT,KM)
#endif
ENDMODULE mod_streamv

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_streamr
#ifdef streamr
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)      :: stxr,styr
  REAL                                       :: sxr(IMT,MR,LOV)
  REAL                                       :: syr(JMT,MR,LOV)
#endif
ENDMODULE mod_streamr

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
#ifdef sediment
MODULE mod_sed
  !  REAL :: wsed,rhos,D,critvel,T,cwamp,kincrit
  REAL                                       :: wsed, partdiam
  REAL                                       :: rhos, cwamp, twave
  REAL                                       :: critvel, kincrit

  INTEGER                                    :: nsed=0, nsusp=0
  LOGICAL                                    :: res
ENDMODULE mod_sed

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_orbital
  REAL, ALLOCATABLE, DIMENSION(:)            :: orb
ENDMODULE mod_orbital
#endif

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===



































 ! ===   OCCAM  === 
!!$#if defined occ66
!!$#ifdef mod1
!!$  INTEGER, PARAMETER ::  IMT1=1440,JMT1=577
!!$  INTEGER, PARAMETER ::  IMT=IMT1,JMT=JMT1
!!$#endif
!!$#ifdef mod2 
!!$  INTEGER, PARAMETER ::  IMT2=438,JMT2=434
!!$#endif
!!$  INTEGER, PARAMETER ::  KM=66,JMAX=250,LBT=3
!!$  INTEGER, PARAMETER ::  IDR=1171,JEQ=313
!!$  INTEGER, PARAMETER ::  NTRACMAX=3*1000*1000
!!$#endif
!!$  
!!$  
!!$  ! ===   SKB  === 
!!$#ifdef simp
!!$  INTEGER, PARAMETER :: IMT=174, JMT=121
!!$#endif
!!$#ifdef fors
!!$  INTEGER, PARAMETER ::  IMT=241,JMT=241
!!$#endif
!!$#if defined fors || simp 
!!$  INTEGER, PARAMETER ::  KM=39,JMAX=JMT,LBT=4
!!$  INTEGER, PARAMETER ::  NTRACMAX=1*100*1000
!!$#endif
!!$  
!!$  ! ===   TEST  === 
!!$#if defined tes
!!$  !INTEGER, PARAMETER :: IMT=225,JMT=IMT,KM=10,JMAX=JMT,LBT=4
!!$  INTEGER, PARAMETER :: IMT=105,JMT=IMT,KM=10,JMAX=JMT,LBT=4
!!$  INTEGER, PARAMETER :: NTRACMAX=10000
!!$#endif
!!$  ! ===   TUNIS-ROMS  === 
!!$#if defined tun
!!$  INTEGER, PARAMETER ::  (IMT=249,JMT=258,KM=20,JMAX=JMT,LBT=4)
!!$  INTEGER, PARAMETER ::  NTRACMAX=100000
!!$#endif
!!$  
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
