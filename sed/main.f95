!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
PROGRAM main
 
IMPLICIT none
#include "../param.h"

common/buoyancy/tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe
real            tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe

common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
integer    ints,intstart,intend,intrun,intspin,intstep,intmin,intmax

common/namn/name,namep,directory
CHARACTER(LEN=8) :: name,namep
CHARACTER(LEN=27) :: directory

common /domain/ienw(NEND),iene(NEND),jens(NEND),jenn(NEND),mask(imt,jmt)
integer ienw,iene,jens,jenn,mask

#if defined sediment
common/sed/wsed,rhos,D,critvel,T,cwamp
REAL wsed,rhos,D,critvel,T,cwamp
common/orbital/orb(km)
REAL orb
#endif

logical skriva
integer nff,isec,idir,nqua,num,ist1,ist2,jst1,jst2,kst1,kst2,ncoor
real voltr
real*8 tseas,tday,tyear,dtmin

print *,'ORM starts'
call system('date')

! handling quantities concerning coordinates
call coordinat

tyear=365.25 * 24. * 3600.
tday=24. * 3600.

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????
! intspin=number of initial datasets

#if defined time
! time step between data sets

tseas=2.*tday
intmin =0  
intmax =365*12    ! maximum length of RCO fields
intspin=365*12    ! intspin=trajectory release period
intrun =100*365   ! trajectory run 
intstep=2         ! pos if fw & neg if bw

if(intstep.gt.0) then ! forward
 intstart=intmin          
 intend  =intmax
elseif(intstep.lt.0) then ! backward
 intstart=intmax
 intend  =intmin
 intspin =-intspin
 intrun  =-intrun   
endif

#else
intend  =intstart        ! end 
#endif

print *,'intspin=',intspin

! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
      ncoor = 1
! write trajectory (true or false)skriva = .false.
skriva = .false.
! open the files for the storage of data 
! name of current trajectory run

name='ne.fw.a0'    ! nff= 1,isec=2,idir= 1

! name of directory where the files will be written

directory='/Volumes/sjo5/data/sed/orm/'
print *,'writes trajectories in ',directory

#if defined rerun
!      namep='fv.fw.aa' 
#endif

!-------- initial directions all in MODEL COORDINATES -----------------

! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
nff=sign(1,intstep)
!nff=1
! isec=1=meridional; isec=2=zonal; isec=3=vertical cross section
isec=1
! if idir= 1 follow positive direction (eastward/northward)  
! if idir=-1 follow negative direction (westward/southward)
! if idir=0  both directions
idir=-1
! number of trajectories can be set by
! 1=nqua: "num" which sets a constant number of trajec. in all boxes
! 2=nqua: set by chosen trajectory transport (voltr) 
! 3=nqua: set by chosen trajectory volume (voltr) 
nqua=1
if(nqua.eq.1) then
! number of trajectories (per time resolution)
num=100
elseif(nqua.eq.2) then
! volume transport of each trajectory in m3/s
!voltr=10. for neva
voltr=100. 
elseif(nqua.eq.3) then
! volume of each trajectory in 10^6 m3
voltr=0.
endif

#if defined sediment
!---------------sedimentation parameters------------------------------      
! Input section for sedimentation modelling.
! parameters for sedimentation:

! particle diameter in mm: clay 0.0005-0.002, silt 0.002-0.06,
! fine sand 0.06-0.2 (medium sand 0.2-0.6, coarse sand 0.6-2, gravel>2)
D = 0.001
 
! density of quartz particle: 2600-2650 g/cm^3, mean value 2620
rhos = 2620

! dynamic viscosity of water, value for 10 degrees Celsius
! (if sedvel called from main.f95 = constant viscosity used) 
!      visc = 0.00131
   
! constant for approximating wave amplitude, a = cwamp*U(surface)
cwamp=20

! approximative  peak period. Average 4s for Baltic proper
T=8

! critical bottom velocity for resuspension
critvel=0.1
!______________________________________________________________________
#endif

! defining starting section 

! trace Gdansk water
!       ist1=140
!       ist2=170
!       jst1=15
!       jst2=27


! Landsort
!      ist1=127
!      ist2=134
!      jst1=139
!      jst2=143
!      kst1=35
!      kst2=km
! Riga
!       ist1=199
!       ist2=imt 
!       jst1=90
!       jst2=139

!  Lilla B?lt
!       ist1=1
!       ist2=20
!       jst1=47
!       jst2=jst1

!  Stora B?lt
!       ist1=26
!       ist2=33
!       jst1=31
!       jst2=jst1

!  ?resund 
!       ist1=53
!       ist2=63
!       jst1=56
!       jst2=jst1

! Wistula
!       ist1=150
!       ist2=ist1
!       jst1=16
!       jst2=jst1

! Oder
!       ist1=80
!       ist2=ist1
!       jst1=4
!       jst2=jst1

! Nemunas ej r?tt
!       ist1=56
!       ist2=ist1
!       jst1=181
!       jst2=jst1

! Neva
       ist1=317
       ist2=ist1
       jst1=183
       jst2=jst1
! Finska viken
!       ist1=200
!       ist2=ist1
!       jst1=153
!       jst2=213
!       jst1=160
!       jst2=jst1
! Entire Baltic
!       ist1=1
!       ist2=imt
!       jst1=1
!       jst2=jmt

      kst1=1
!      kst1=km
      kst2=km

! defining end sections

! Lilla Bält
ienw(1)=1
iene(1)=25
jens(1)=47+1
jenn(1)=jmt
! Stora Bält
ienw(2)=26
iene(2)=36
jens(2)=33+1
jenn(2)=jmt
! Öresund 
ienw(3)=37
iene(3)=63
jens(3)=56+1
jenn(3)=jmt
! utanför Finska viken
!        ienw(1)=0
!        iene(1)=200
!        jens(1)=0
!        jenn(1)=jmt

!        ienw(1)=0
!        iene(1)=ist1
!        jens(1)=0
!        jenn(1)=jmt
! iterative timesteping in seconds
#ifdef time
      dtmin=0.2*tseas
#endif

#ifdef tempsalt
! water mass properties with minimum & maximum on temp, salt and density
!
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
#endif

! trajectory path
if(skriva) open(56,file=directory//'traj.'//name)
! exit position
open(57,file=directory//'traj.ut.'//name)
! entrence position
open(58,file=directory//'traj.in.'//name)


print 999,name,intstart,intspin,intend,intrun,nff,isec,idir,nqua,num,voltr,&
tmin0,tmax0,smin0,smax0,rmin0,rmax0,ist1,ist2,jst1,jst2,kst1,kst2
999 format(' name=',a8,' intstart=',i4,' intspin=',i4,' intend=',i7,' intrun=',i7,/,&
      ' nff=',i2,' isec=',i2,' idir=',i4,' nqua=',i2,' num=',i4,&
      ' voltr=',f9.0,/,&
      ' tmin0=',f7.2,' tmax0=',f7.2,' smin0=',f7.2,' smax0=',f7.2,&
      ' rmin0=',f7.2,' rmax0=',f7.2,/,&
      ' ist1=',i4,' ist2=',i4,' jst1=',i4,' jst2=',i4,' kst1=',i2,&
      ' kst2=',i2)

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
# include "readfield.f95"
# include "stat.f95"

#if defined sediment
# include "../sedimentation.f95"
#endif

!23456789012345678901234567890123456789012345678901234567890123456789012
