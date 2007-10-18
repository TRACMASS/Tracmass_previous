!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
# include "../modules.f95"
PROGRAM main
USE mod_param
USE mod_name
USE mod_time 
USE mod_domain
USE mod_buoyancy

IMPLICIT none

INTEGER ngcm,i,j,l

call coordinat

ngcm=1   ! hours between GCM datasets
tseas=dble(ngcm)*3600.d0 ! time step between data sets

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????

intmin=1

!intmax=24*365   ! max OGCM data simulation
intmax=24*31   ! max OGCM data simulation

!intspin=365*24+1  ! release of  trajectories
!intspin=30*24+1  ! release of  trajectories
intspin=1  ! release of trajectories

!intrun =intspin+15*24 ! ttotal run of trajectories
!intrun =intspin+30*24 
intrun =intmax

intstep=1   ! pos if fw & neg if bw

if(intstep.gt.0) then
! forward
 intstart=intmin          
 intend  =intmax
elseif(intstep.lt.0) then
! backward 
!       intstart=intmax
 intstart=intspin
 intend  =intmin
 intspin =-intspin
 intrun  =-intrun   
endif

! starting date
ihour=0
!ihour=14
iday=1
!iday=15
imon=1
imon=7
if(intstep.lt.0) then
ihour=24
ihour=6
iday=31
iday=4
imon=12
imon=6
endif
#if defined fors
iyear=1988
#elif defined simp
iyear=1981
#endif 

! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
ncoor = 0
! kriva=0=no writing
!       1=write at time intervals of gcm datasets
!       2=write each grid-crossing and timne change
!       3=write at each iteration (all the time)
!       4=write only start and end positions
kriva=1
! open the files for the storage of data 
! name of current trajectory run

name='si.fw.n4'    ! nff= 1,isec=2,idir= 1

! name of directory where the files will be written

#if defined fors && desk
 directory='/Volumes/hav1/data/for/'
#elif defined fors && lap
 directory='/Users/doos/data/for12/'
#elif defined simp && desk
 directory='/Volumes/hav1/data/sim/'
#elif defined simp && lap
 directory='/Users/doos/data/sim12/'
#else
 stop 49567
#endif

#if defined rerun
!      namep='fv.fw.aa' 
#endif 

!-------- initial directions all in MODEL COORDINATES -----------------
! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
nff=sign(1,intstep)
! isec=1=meridional; isec=2=zonal; isec=3=vertical section; isec=4=fill middle of T-box
isec=4
! if idir= 1 follow positive direction (eastward/northward)  
! if idir=-1 follow negative direction (westward/southward)
! if idir=0  both directions
idir=0
! number of trajectories can be set by
! 1=nqua: "num" which sets a !onstant number of trajec. in all boxes
! 2=nqua: set by chosen trajectory transport (voltr) 
! 3=nqua: set by chosen trajectory volume (voltr) 
! 4=Same as 1 but for volume in m3
! 5=number of trajectories set by mask
! 6=only one trajectory per release box defined by mask
nqua=6
if(nqua.eq.1 .or. nqua.eq.4) then     ! number of trajectories (per time resolution)
 num=1
elseif(nqua.eq.2) then ! volume transport of each trajectory in m3/s
 voltr=1.
elseif(nqua.eq.3) then ! volume of each trajectory in 10^6 m3
 voltr=10000.
endif

timax=2.*24.*3600.  ! maximum length of a trajectory

! ______ defining starting section ____________

#if defined fors
ist1=1
ist2=imt
jst1=1
jst2=jmt
kst1=1
kst2=km

#endif
#if defined simp
ist1=1
ist2=imt
jst1=1
jst2=jmt
kst1=1
kst2=km
#endif

! defining end sections
! east
ienw(1)=1
iene(1)=1
jens(1)=1
jenn(1)=jmt
! west
ienw(2)=imt-1
ienw(2)=imt
iene(2)=imt
jens(2)=1
jenn(2)=jmt
! south
ienw(3)=1
iene(3)=imt
jens(3)=1
jenn(3)=2
jenn(3)=1
! north
ienw(4)=1
iene(4)=imt
jens(4)=jmt-1
jenn(4)=jmt

#if defined fors
open(21,file=directory//'topo/mask0',form='unformatted')
read(21) mask
close(21)
#endif
#if defined simp
!open(21,file=directory//'topo/utslappsyd100',form='unformatted')
!open(21,file=directory//'topo/utslappsyd',form='unformatted')
!open(21,file=directory//'topo/utslappgrid',form='unformatted')
open(21,file=directory//'topo/mask',form='unformatted')
read(21) mask
close(21)
!open(21,file=directory//'topo/mask',form='unformatted')
!read(21)izon
!close(21)

do i=1,imt
 do j=1,jmt
!  if(izon(i,j).eq.-99.) mask(i,j)=-99.
!  if(izon(i,j).eq.4) mask(i,j)=-99.
!  if(i.ge.imt-1 .or. j.le.2 .or.j.ge.jmt-1) mask(i,j)=-99.
!  if(j.ge.55) mask(i,j)=0 !eleminate northern point
!  if(i.ne.47 .or. j.ne.45) mask(i,j)=0 !endast en trajektoria
 enddo
enddo

!do j=jmt,1,-1
!print 88,(mask(i,j),i=1,150)
!enddo
!stop 4967
!88 format(300i1)

!print *,'   i    j   vikt'
!l=0

#endif

l=0
do j=1,jmt
do i=1,imt
!if(mask(i,j).ne.0) print 88,i,j,mask(i,j)
l=l+mask(i,j)
enddo
enddo
print *,'total release points=',l
88 format(3i5)


!       print *,(mask(90,j),j=1,jmt)
!       stop 2345

! iterative timesteping in seconds
#ifdef time
iter=5 ! iteration between two gcm data sets
dstep=1.d0/dble(iter)
dtmin=dstep*tseas
#endif

! water mass properties with minimum & maximum on temp, salt and density
! for starting a trajectory (active only with option tempsalt)
tmin0 = -50.
tmax0 = 20.
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

if(kriva.ne.0) open(56,file=directory//'orm/traj.'//name) ! trajectory path
open(57,file=directory//'orm/traj.ut.'//name)         ! exit position
open(58,file=directory//'orm/traj.in.'//name)         ! entrence position

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
!# include "../interp.f95"
# include "../interp2.f95"
# include "../pos.f95"
# include "../arclength.f95"
# include "../writepsi.f95"
# include "../turb.f95"
# include "readfield.f95"
!# include "../stat.f95"

