!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
PROGRAM main

IMPLICIT none
#include "../param.h"

common/buoyancy/tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe,&
timax
real            tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe,&
timax

common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax, &
           ntime,iyear,imon,iday,ihour
integer    ints,intstart,intend,intrun,intspin,intstep,intmin,intmax, &
           ntime,iyear,imon,iday,ihour

common/namn/name,namep,directory
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=23) :: directory

common /domain/ienw(NEND),iene(NEND),jens(NEND),jenn(NEND),mask(imt,jmt)
integer ienw,iene,jens,jenn,mask

#if defined sediment 
common/sed/wsed,rhos,D,critvel,T,c1,kincrit
REAL wsed,rhos,D,critvel,T,c1,kincrit
common/orbital/orb(km)
REAL orb
#endif

logical skriva
integer nff,isec,idir,nqua,num,ist1,ist2,jst1,jst2,kst1,kst2,ncoor,i,j,l,izon(imt,jmt)
real voltr
real*8 tseas,tday,tyear,dtmin

print *,'ORM starts'
call system('date')

! handling quantities concerning coordinates 
!call coordinat

tyear=365.25 * 24. * 3600.
tday=24. * 3600.

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????
! intspin=number of initial datasets

#if defined time
! time step between data sets
tseas=3600.

intmin=1
intmax=24*365   ! max OGCM data simulation

intspin=365*24+1+10  ! release of  trajectories
intrun=3*365*24+10

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
#else
 intend  =intstart        ! end 
#endif

! starting date
ihour=0
iday=1
imon=1
if(intstep.lt.0) then
ihour=24
iday=31
imon=12
endif
#if defined fors
iyear=1988
#elif defined simp
iyear=1981
#endif 

! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
ncoor = 0
! write trajectory (true or false)
skriva = .false.
!skriva = .true.
! open the files for the storage of data 
! name of current trajectory run
name='si.fw.l2'    

! name of directory where the files will be written
#if defined fors
 directory='/Volumes/sjo3/data/for/'
!directory='/Users/doos/data/skb12/'
#endif 
#if defined simp
 directory='/Volumes/sjo4/data/sim/'
!directory='/Users/doos/data/skb12/'
#endif 
print *,'directory=',directory

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
nqua=5
if(nqua.eq.1 .or. nqua.eq.4) then     ! number of trajectories (per time resolution)
 num=1
elseif(nqua.eq.2) then ! volume transport of each trajectory in m3/s
 voltr=1.
elseif(nqua.eq.3) then ! volume of each trajectory in 10^6 m3
 voltr=10000.
endif

timax=1000.*365.*24.*3600.  ! maximum length of a trajectory
#if defined sediment
!---------------sedimentation parameters------------------------------      
! Input section for sedimentation modelling.
! parameters for sedimentation:

! particle diameter in mm: clay 0.0005-0.002, silt 0.002-0.06,
! fine sand 0.06-0.2 (medium sand 0.2-0.6, coarse sand 0.6-2, gravel>2)
D = 0.001   !0.01=silt & 0.001=lera
 
! density of quartz particle: 2600-2650 g/cm^3, mean value 2620
rhos = 2620

! dynamic viscosity of water, value for 10 degrees Celsius
! (if sedvel called from main.f95 = constant viscosity used) 
!      visc = 0.00131
   
! constant for approximating wave amplitude, a = c1*U(surface)
c1=20

! approximative  peak period. Average 4s for Baltic proper
T=8

! critical bottom velocity for resuspension
critvel=0.1

! critical kinetic energy for reducing the sedimentation velocity
kincrit=1500.

!______________________________________________________________________
#endif
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
! west
ienw(1)=1
iene(1)=1
jens(1)=1
jenn(1)=jmt
! east
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

#if defined fors
open(21,file=directory//'topo/mask0',form='unformatted')
read(21) mask
close(21)
#endif
#if defined simp
!open(21,file=directory//'topo/utslappsyd100',form='unformatted')
!open(21,file=directory//'topo/utslappsyd',form='unformatted')
open(21,file=directory//'topo/utslappgrid',form='unformatted')
!open(21,file=directory//'topo/mask',form='unformatted')
read(21) mask
close(21)
!open(21,file=directory//'topo/mask',form='unformatted')
!read(21)izon
!close(21)

!do i=1,imt
! do j=1,jmt
!  if(izon(i,j).eq.-99.) mask(i,j)=-99.
!  if(izon(i,j).eq.4) mask(i,j)=-99.
!  if(i.ge.imt-1 .or. j.le.2 .or.j.ge.jmt-1) mask(i,j)=-99.
! enddo
!enddo

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
if(mask(i,j).ne.0) print 88,i,j,mask(i,j)
l=l+mask(i,j)
enddo
enddo
print *,'total release points=',l
88 format(3i5)


!       print *,(mask(90,j),j=1,jmt)
!       stop 2345

! iterative timesteping in seconds
#ifdef time
dtmin=1.*tseas
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

if(skriva) open(56,file=directory//'orm/traj.'//name) ! trajectory path
open(57,file=directory//'orm/traj.ut.'//name)         ! exit position
open(58,file=directory//'orm/traj.in.'//name)         ! entrence position

!??????????????????????????????????? END ???????????????????????????????

! trajectory loops 

call loop(isec,idir,nqua,num,ist1,ist2,jst1,jst2,kst1,kst2,ncoor,nff,&
          tseas,tday,tyear,dtmin,voltr,skriva  )

call writepsi

print *,'the very end of run ',name
call system('date')

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
#if defined sediment
# include "../sedimentation.f95"
#endif
