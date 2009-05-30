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
tseas=real(ngcm)*3600.d0 ! time step between data sets

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????
! intspin=number of initial datasets

#ifdef time
#ifdef fiveday
!______________________ 5 day datasets
      intmin=1465   ! 5 jan 1985 
      intmax=8395   ! 27 dec 2003  
!      intspin=1*365-5  
      intspin=90
      intrun  =200*365
      intstep=5   ! pos if fw & neg if bw
      if(intstep.gt.0) then
! forward
       intstart=intmin          
!       intstart=intmin+6*365        
       intend  =intmax
      elseif(intstep.lt.0) then
! backward 
       intstart=intmax
!       intstart=intmax-6*365  
       intend  =intmin
       intspin =-intspin
       intrun  =-intrun   
      endif
#endif
!#ifdef year
!______________________ yearly datasets
!      tseas=tyear
! forward
!      intstart=1985            
!      intend  =intstart+200    
!      intstep =1 
!      intspin=2001-1985                  
! backward
!      intstart=2001            
!      intend  =intstart-200    
!      intstep =-1 
!      intspin=1985-2001                  
!#endif
#endif
#if defined stat
      intstart=9102            ! start 
      intend  =intstart        ! end 
#endif

ihour=0-6
iday=15-3           
imon=5
iyear=1994


! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
ncoor = 1
! kriva=0=no writing
!       1=write at time intervals of gcm datasets
!       2=write each grid-crossing and timne change
!       3=write at each iteration (all the time)
!       4=write only start and end positions
kriva=0
! open the files for the storage of data 
! name of current trajectory run

!      name='no.fw.f5'    ! nff= 1,isec=2,idir=-1
!       name='no.fw.z5'    ! nff= 1,isec=2,idir=-1
!      name='no.bw.f5'    ! nff=-1,isec=2,idir=-1
!      name='mn.fw.s2'    ! nff= 1,isec=2,idir= 1
!      name='mn.bw.s2'    ! nff=-1,isec=2,idir= 1
!      name='ms.fw.s2'    ! nff= 1,isec=2,idir=-1 
!      name='ms.bw.s2'    ! nff=-1,isec=2,idir=-1
!       name='dr.bw.f4'    ! nff=-1,isec=1,idir=-1
!      name='dr.fw.f4'    ! nff= 1,isec=1,idir= 1
       name='dr.fw.a0'    ! nff= 1,isec=1,idir= 1
! name of directory where the files will be written

#if defined desk
 directory='/Volumes/hav2/data/occ/'
#elif defined lap
 directory='/Users/doos/data/occ12/'
#else
 stop 49567
#endif

#if defined rerun
!namep='dr.bw.e4' 
#endif

!-------- initial directions all in MODEL COORDINATES -----------------
! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
#if defined time
      nff=sign(1,intstep)
#else
      nff=1
#endif
! isec=1=meridional; isec=2=zonal; isec=3=vertical cross section
      isec=1
! if idir= 1 follow positive direction (eastward/northward)  
! if idir=-1 follow negative direction (westward/southward)
! if idir=0  both directions
      idir=1 
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
!       voltr=100.
        voltr=1.e4
        voltr=1.e5
        voltr=1.e6
      elseif(nqua.eq.3) then
! volume of each trajectory in 10^6 m3
       voltr=0.
      endif

timax=10000.*365.*24.*3600.  ! maximum length of a trajectory


! defining starting section 
       kst1=1
       kst2=km
! northern boundary
!       ist1=1
!       ist2=imt
!       jst1=jmax-2
!       jst2=jmax-2
! middle of Southern Ocean between ACC and Subtropical Gyre 
!        ist1=1
!        ist2=imt
!        jst1=jgyr
!        jst2=jgyr
! drake 
        ist1=idr
        ist2=idr
        jst1=45
        jst2=90

! iterative timesteping in seconds
#ifdef time
iter=5 ! iteration between two gcm data sets
dstep=1.d0/dble(iter)
dtmin=dstep*tseas
#endif

#ifdef tempsalt
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
#endif

mask=-1  ! att fixa om ....


if(kriva.ne.0) open(56,file=directory//'orm/traj.'//name) ! trajectory path
open(57,file=directory//'orm/traj.ut.'//name)         ! exit position
open(58,file=directory//'orm/traj.in.'//name)         ! entrence position


!??????????????????????????????????? END ???????????????????????????????
iday0=iday
imon0=imon
iyear0=iyear

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
# include "readfield.f95"
# include "stat.f95"

