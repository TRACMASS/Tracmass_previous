!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
# include "../modules.f95"
PROGRAM main
USE mod_param
USE mod_name
USE mod_time
USE mod_domain
USE mod_buoyancy

IMPLICIT none

INTEGER ngcm,i,j

call coordinat

ngcm=5*24   ! hours between GCM datasets
tseas=dble(ngcm)*3600.d0

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????

intmin =1      ! 
intmax =10*365/5 ! maximum length of OGCM fields
!intmax =3*365/5  ! maximum length of OGCM fields
!intmax =360/5     ! maximum length of OGCM fields
intspin=intmax    ! intspin=trajectory release period
!intspin=10*365/5    ! intspin=trajectory release period
!intspin=2    ! intspin=trajectory release period

intrun =1000*365/5 ! trajectory run 
!intrun =1*365+10 ! trajectory run 
intstep=1      ! pos if fw & neg if bw
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

! starting date
ihour=0 ! alaways zero for orca

iday=2-10
imon=7
iyear=1990

!iday=1-10
!imon=10
!iyear=1992

! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
ncoor = 2
! kriva=0=no writing
!       1=write at time intervals of gcm datasets
!       2=write each grid-crossing and timne change
!       3=write at each iteration (max writing)
!       4=write only start and end positions
kriva=0
! name of current trajectory run
name='op.fw.n0'   
! name of directory where the files will be written
#if defined desk
 directory='/Volumes/hav3/data/orc/'
#elif defined lap
 directory='/Users/doos/data/orc12/'
#else
 stop 49567
#endif
#if defined rerun
namep='op.fw.n0' 
#endif

!-------- initial directions all in MODEL COORDINATES -----------------
! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
 nff=sign(1,intstep)
! isec=1=meridional(y-z); 2=zonal(x-z); 3=horiz(x-y) sections; 4=fill middle of T-box
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
! num=100
! num=10000
 num=1
elseif(nqua.eq.2) then ! volume transport of each trajectory in m3/s
 voltr=1000.
elseif(nqua.eq.3) then ! volume of each trajectory in 10^6 m3
 voltr=0.
endif

timax=1000.*365.*24.*3600.  ! maximum length of a trajectory

! ______ defining starting section ____________
kst1=1
kst2=KM
!kst1=KM

!ist1=150   
!ist2=400

!ist1=120 ! eastern pacific
!ist2=450

!ist1=52
!ist2=415

!ist1=250
!ist2=ist1

ist1=1
ist2=IMT

!jst1=254  ! jeq=250
!jst2=254

!jst1=256  ! jeq=250
jst1=252  ! jeq=250
!jst1=505  ! jeq=250
jst2=jst1

!ist1=300
!ist2=ist1

! iterative timesteping in seconds
#ifdef time
iter=5 ! iteration between two gcm data sets
dstep=1.d0/dble(iter)
dtmin=dstep*tseas
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

! Atlantic
ienw(1)=451
iene(1)=641
jens(1)=0
jenn(1)=jst1  

! western Indian
ienw(2)=641
iene(2)=imt+1
jens(2)=0
jenn(2)=jst1  

! eastern Indian
ienw(3)=0
iene(3)=50
jens(3)=0
jenn(3)=jst1  

! western boundary Pacific
ienw(4)=50
iene(4)=120
jens(4)=0
jenn(4)=jst1  

! interior Pacific
ienw(5)=120
iene(5)=450
jens(5)=0
jenn(5)=jst1  

! Bering strait
ienw(1)=1
iene(1)=IMT
jens(1)=419
jenn(1)=JMT 

if(kriva.ne.0) open(56,file=directory//'orm/traj.'//name) ! trajectory path
open(57,file=directory//'orm/traj.ut.'//name)         ! exit position
open(58,file=directory//'orm/traj.in.'//name)         ! entrence position

! start area
mask=-1
do j=1,JMT
 do i=1,IMT
  if(451.le.i .and. i.le.641) mask(i,j)=0 ! exlude the Atlantic
 enddo
enddo

!??????????????????????????????????? END ???????????????????????????????

! trajectory loops 

call loop

stop
end program main

!______________ END OF MAIN PROGRAM _______________________________

# include "../loop.f95"
# include "../vertvel.f95"
# include "../coord.f95"
# include "../cross.f95"
# include "../interp2.f95"
# include "../pos.f95"
# include "../arclength.f95"
# include "../writepsi.f95"
# include "../writetracer.f95"
# include "readfield.f95"
