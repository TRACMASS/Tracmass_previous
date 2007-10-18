!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x

MODULE mod_param
!__________________________  OCCAM ____________________________________
#if defined occ66
#ifdef mod1
INTEGER, PARAMETER ::  IMT1=1440,JMT1=577
INTEGER, PARAMETER ::  IMT=IMT1,JMT=JMT1
#endif
#ifdef mod2 
INTEGER, PARAMETER ::  IMT2=438,JMT2=434
#endif
INTEGER, PARAMETER ::  KM=66,JMAX=250,LBT=3
INTEGER, PARAMETER ::  IDR=1171,JEQ=313
INTEGER, PARAMETER ::  NTRACMAX=3*1000*1000
#endif

!__________________________  ORCA _____________________________________
#if defined orca
INTEGER, PARAMETER ::  IMT=720,JMT=511,KM=46,JMAX=JMT,LBT=5
INTEGER, PARAMETER ::  NTRACMAX=7*1000*1000
#endif

!__________________________  RCO ______________________________________
#if defined rco
INTEGER, PARAMETER :: IMT=320,JMT=362,KM=41,JMAX=JMT,LBT=1
INTEGER, PARAMETER :: NTRACMAX=20*1000*1000
#endif

!__________________________  SKB ______________________________________
#ifdef simp
INTEGER, PARAMETER :: IMT=174, JMT=121
#endif
#ifdef fors
INTEGER, PARAMETER ::  IMT=241,JMT=241
#endif
#if defined fors || simp 
INTEGER, PARAMETER ::  KM=39,JMAX=JMT,LBT=4
INTEGER, PARAMETER ::  NTRACMAX=1*100*1000
#endif

!__________________________  TEST ______________________________________
#if defined tes
!INTEGER, PARAMETER :: IMT=225,JMT=IMT,KM=10,JMAX=JMT,LBT=4
INTEGER, PARAMETER :: IMT=105,JMT=IMT,KM=10,JMAX=JMT,LBT=4
INTEGER, PARAMETER :: NTRACMAX=10000
#endif
!__________________________  TUNIS-ROMS ________________________________
#if defined tun
INTEGER, PARAMETER ::  (IMT=249,JMT=258,KM=20,JMAX=JMT,LBT=4)
INTEGER, PARAMETER ::  NTRACMAX=100000
#endif
!_______________________________________________________________________

INTEGER, PARAMETER :: MR=1001,NEND=LBT+1

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
INTEGER, PARAMETER :: KD=KM
#else
INTEGER, PARAMETER :: KD=1
#endif

INTEGER nff,isec,idir,nqua,num,ist1,ist2,jst1,jst2,kst1,kst2,ncoor,kriva,iter
REAL*8 tseas,tday,tyear,dtmin,voltr,tstep,dstep,tss



ENDMODULE mod_param

!________________________________________________________________________________________
MODULE mod_coord
USE mod_param
REAL*8 dx,dy,deg,stlon1,stlat1
REAL*8 zw(0:KM)
REAL*8, DIMENSION(JMT) :: csu,cst
INTEGER idmax(12,1000:3000)
ENDMODULE mod_coord
!________________________________________________________________________________________
MODULE mod_time
USE mod_param
INTEGER ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
INTEGER iyear,imon,iday,ihour,iyear0,imon0,iday0,yearmin,yearmax
INTEGER*8 ntime
ENDMODULE mod_time
!________________________________________________________________________________________
MODULE mod_grid
USE mod_param
REAL*8 dxdy(IMT,JMT),dztb(IMT,JMT,KD),dz(KM),rmin,dr,tmin,dtemp,smin,dsalt
INTEGER kmt(IMT,JMT)
ENDMODULE mod_grid
!________________________________________________________________________________________
MODULE mod_buoyancy
REAL*4 tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe,timax
ENDMODULE mod_buoyancy
!________________________________________________________________________________________
MODULE mod_domain
USE mod_param
INTEGER ienw(NEND),iene(NEND),jens(NEND),jenn(NEND),mask(imt,jmt)
ENDMODULE mod_domain
!________________________________________________________________________________________
MODULE mod_vel
USE mod_param
REAL*4, DIMENSION(IMT,0:JMAX,KM,NST) :: u,v
REAL*4 hs(IMT,JMAX,NST),rand(6)
REAL*8 w(0:km)
REAL*8 ff
ENDMODULE mod_vel
!________________________________________________________________________________________
#ifdef tempsalt
MODULE mod_dens
USE mod_param
REAL*4, DIMENSION(IMT,JMAX,KM,NST) :: tem,sal,rho
ENDMODULE mod_dens
#endif
!________________________________________________________________________________________
MODULE mod_turb
REAL upr(6)
ENDMODULE mod_turb
!________________________________________________________________________________________
MODULE mod_name
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=23) :: directory
ENDMODULE mod_name 
!________________________________________________________________________________________
#ifdef streamxy
MODULE mod_streamxy
USE mod_param
REAL stxyy(IMT,JMT,LBT),stxyx(IMT,JMT,LBT)
ENDMODULE mod_streamxy 
#endif
!________________________________________________________________________________________
#ifdef streamv
MODULE mod_streamv
USE mod_param
REAL stxz(IMT,KM,LBT),styz(JMT,KM,LBT)
ENDMODULE mod_streamv 
#endif
!________________________________________________________________________________________
#ifdef streamr
MODULE mod_streamr
USE mod_param
REAL :: stxr(IMT,MR,LBT,LOV)
REAL :: styr(JMT,MR,LBT,LOV)
ENDMODULE mod_streamr 
#endif
!________________________________________________________________________________________
#ifdef tracer
MODULE mod_tracer
USE mod_param
REAL tra(IMT,JMT,KM)
ENDMODULE mod_tracer
#endif
!________________________________________________________________________________________
#ifdef sediment
MODULE mod_sed
REAL wsed,rhos,D,critvel,T,c1,kincrit
ENDMODULE mod_sed

MODULE mod_orbital
REAL orb(KM)
ENDMODULE mod_orbital
#endif
!________________________________________________________________________________________
