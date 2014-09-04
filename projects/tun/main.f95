!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
PROGRAM main

IMPLICIT none
#include "../param.h"

common/buoyancy/tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe
real            tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe

common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
integer ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
 
common/namn/name,namep,directory
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=27) :: directory

common /domain/ienw(NEND),iene(NEND),jens(NEND),jenn(NEND),mask(imt,jmt)
INTEGER ienw,iene,jens,jenn,mask

logical skriva 
integer nff,isec,idir,nqua,num
integer ist1,ist2,jst1,jst2,kst1,kst2,ncoor
real voltr
real*8 tseas,tday,tyear,dtmin

print *,'ORM starts'
call system('date')

tyear=365.25d0 * 24.d0 * 3600.d0
tday=24.d0 * 3600.d0

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????

tseas=15.d0*3600.d0 ! detta är nog inte exakt...fråga Jihène!

!______________________ 5 day datasets
intmin =0      ! 
intmax =10*365 ! maximum length of OGCM fields
intspin=365    ! intspin=trajectory release period
intrun =80*365 ! trajectory run 
intstep=5      ! pos if fw & neg if bw
if(intstep.gt.0) then
! forward
intstart=intmin          
!       intstart=intmin+6*365        
intend  =intmax
elseif(intstep.lt.0) then
! backward 
intstart=intmax
intend  =intmin
intspin =-intspin
intrun  =-intrun   
endif

! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
ncoor = 1
! write trajectory (.true. or .false.)
skriva = .false.
! open the files for the storage of data 
! name of current trajectory run
name='op.fw.b5'   
! name of directory where the files will be written
directory='/Volumes/sjo5/data/tun/orm/'
directory='/Users/doos/data/abcdefghi/'
print *,'writes trajectories in ',directory
#if defined rerun
namep='tu.fw.a0' 
#endif

!-------- initial directions all in MODEL COORDINATES -----------------
! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
 nff=sign(1,intstep)
! isec=1=meridional; isec=2=zonal; isec=3=vertical cross section
isec=2
! if idir= 1 follow positive direction (eastward/northward)  
! if idir=-1 follow negative direction (westward/southward)
! if idir=0  both directions
idir=1

! number of trajectories can be set by
! 1=nqua: "num" which sets a constant number of trajec. in all boxes
! 2=nqua: set by chosen trajectory transport (voltr) 
! 3=nqua: set by chosen trajectory volume (voltr) 
nqua=1
if(nqua.eq.1) then     ! number of trajectories (per time resolution)
 num=1
elseif(nqua.eq.2) then ! volume transport of each trajectory in m3/s
 voltr=500.
elseif(nqua.eq.3) then ! volume of each trajectory in 10^6 m3
 voltr=0.
endif

! ______ defining starting section ____________
kst1=1
kst2=KM

ist1=imt/2
ist2=imt/2

jst1=3
jst2=jmt-2

! iterative timesteping in seconds
#ifdef time
dtmin=0.5*tseas
#endif

! water mass properties with minimum & maximum on temp, salt and density
!
! for starting a trajectory (active only with option tempsalt)
tmin0 = -50.
tmax0 = 400.
smin0 =-500.
smax0 = 400.
rmin0 =-100.
rmax0 = 500.

! for ending a trajectory
tmine =  -50.
tmaxe =  400.
smine = -150.
smaxe =  500.
rmine = -100.
rmaxe =  500.

! defining end sections
! east
ienw(1)=1
iene(1)=1
jens(1)=1
jenn(1)=jmt
! west
ienw(2)=imt
iene(2)=imt
jens(2)=1
jenn(2)=jmt
! south
ienw(3)=1
iene(3)=imt
jens(3)=1
jenn(3)=1
! north
ienw(4)=1
iene(4)=imt
jens(4)=jmt
jenn(4)=jmt


if(skriva) open(56,file=directory//'traj.'//name) ! trajectory path
open(57,file=directory//'traj.ut.'//name)         ! exit position
open(58,file=directory//'traj.in.'//name)         ! entrence position

!??????????????????????????????????? END ???????????????????????????????

! trajectory loops 

call loop(isec,idir,nqua,num,ist1,ist2,jst1,jst2,kst1,kst2,ncoor,nff,&
          tseas,tday,tyear,dtmin,voltr,skriva  )

call writepsi

call system('date')
print *,'the very end of run ',name

stop
end program main

!______________ END OF MAIN PROGRAM _______________________________

# include "../loop.f95"
# include "../vertvel.f95"
# include "../coord.f95"
# include "../cross.f95"
!#include "../interp.f95"
# include "../interp2.f95"
# include "../pos.f95"
!#include "../statocrca.f95"
!#include "../arclength.f95"
# include "../writepsi.f95"
# include "readfield.f95"


 
