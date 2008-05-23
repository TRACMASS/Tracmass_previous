
MODULE mod_param
  INTEGER  :: IMT, JMT, KM, JMAX, LBT, NTRACMAX
  !__________________________  OCCAM ____________________________________
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
!!$  !__________________________  SKB ______________________________________
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
!!$  !__________________________  TEST ______________________________________
!!$#if defined tes
!!$  !INTEGER, PARAMETER :: IMT=225,JMT=IMT,KM=10,JMAX=JMT,LBT=4
!!$  INTEGER, PARAMETER :: IMT=105,JMT=IMT,KM=10,JMAX=JMT,LBT=4
!!$  INTEGER, PARAMETER :: NTRACMAX=10000
!!$#endif
!!$  !__________________________  TUNIS-ROMS ________________________________
!!$#if defined tun
!!$  INTEGER, PARAMETER ::  (IMT=249,JMT=258,KM=20,JMAX=JMT,LBT=4)
!!$  INTEGER, PARAMETER ::  NTRACMAX=100000
!!$#endif
!!$  !_______________________________________________________________________
!!$  
  INTEGER, PARAMETER :: MR=1001
  INTEGER            :: NEND
#if defined time
  INTEGER, PARAMETER :: NST=2,NNRJ=8,NTRJ=7
#endif
#if defined stat
  INTEGER, PARAMETER ::  (NST=1)
#endif
  
#if defined streamts
  INTEGER, PARAMETER :: LOV=3
#else
  INTEGER, PARAMETER :: LOV=1
#endif
  
#if defined sigma 
  INTEGER, PARAMETER :: KD=KM   ! does this work?
#elif defined atm
  INTEGER, PARAMETER :: KD=11  ! should be KD=KM   
#else
  INTEGER, PARAMETER :: KD=1
#endif
  
  INTEGER            :: nff,isec,idir,nqua,num,ist1,ist2,jst1,jst2,kst1,kst2
  INTEGER            :: ncoor,kriva,iter,ngcm
  REAL*8             :: tseas,tday,tyear,dtmin,voltr,tstep,dstep,tss,partQuant
  
ENDMODULE mod_param
!_______________________________________________________________________________
MODULE mod_coord
  REAL*8 dx,dy,deg,stlon1,stlat1
  REAL*8, ALLOCATABLE, DIMENSION(:) :: zw
  REAL*8, ALLOCATABLE, DIMENSION(:) :: csu,cst,dyt,phi
  INTEGER idmax(12,1000:3000)
ENDMODULE mod_coord
!______________________________________________________________________________
MODULE mod_time
  INTEGER ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
  INTEGER iyear,imon,iday,ihour,iyear0,imon0,iday0,yearmin,yearmax
  INTEGER*8 ntime
ENDMODULE mod_time
!______________________________________________________________________________
MODULE mod_grid
  REAL*8, ALLOCATABLE, DIMENSION(:)      :: dz
  REAL*8, ALLOCATABLE, DIMENSION(:,:)    :: dxdy
#if defined ifs || atm
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)  :: dztb 
#else
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: dztb 
#endif
  REAL*8 rmin,dr,tmin,dtemp,smin,dsalt,arcscale
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: kmt
ENDMODULE mod_grid
!______________________________________________________________________________
MODULE mod_buoyancy
  REAL*4 tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe
ENDMODULE mod_buoyancy
!______________________________________________________________________________
MODULE mod_domain
  INTEGER ienw(10),iene(10),jens(10),jenn(10)
  REAL*4 timax
!  INTEGER, ALLOCATABLE, DIMENSION(:)   :: ienw,iene,jens,jenn
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask
ENDMODULE mod_domain
!______________________________________________________________________________
MODULE mod_vel
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: u,v
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:) ::  hs
  REAL*8, ALLOCATABLE, DIMENSION(:) :: w
  REAL*4, DIMENSION(6) :: rand
  REAL*8 ff
ENDMODULE mod_vel
!______________________________________________________________________________
MODULE mod_dens
#ifdef tempsalt
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: tem,sal,rho
#endif
ENDMODULE mod_dens
!______________________________________________________________________________
MODULE mod_turb
#ifdef turb
  REAL upr(6,2)
#endif
ENDMODULE mod_turb
!______________________________________________________________________________
MODULE mod_name
  CHARACTER(LEN=8) :: name,namep
  CHARACTER(LEN=200) :: directory
ENDMODULE mod_name
!______________________________________________________________________________
MODULE mod_streamxy
#ifdef streamxy
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: stxyy, stxyx
#endif
ENDMODULE mod_streamxy
!______________________________________________________________________________
MODULE mod_streamv
#ifdef streamv
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: stxz, styz
#endif
ENDMODULE mod_streamv
!______________________________________________________________________________
MODULE mod_streamr
#ifdef streamr
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: stxr,styr
#endif
ENDMODULE mod_streamr
!______________________________________________________________________________
MODULE mod_tracer
#ifdef tracer
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tra
#endif
ENDMODULE mod_tracer
!______________________________________________________________________________
#ifdef sediment
MODULE mod_sed
  REAL :: wsed,rhos,D,critvel,T,cwamp,kincrit
ENDMODULE mod_sed

MODULE mod_orbital
  REAL, ALLOCATABLE, DIMENSION(:) :: orb
ENDMODULE mod_orbital
#endif
!______________________________________________________________________________
