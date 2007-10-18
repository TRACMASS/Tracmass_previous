!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
# include "../modules.f95"
PROGRAM main
USE mod_param
USE mod_name
USE mod_time 
USE mod_domain
USE mod_buoyancy

IMPLICIT none

INTEGER ngcm,i,j,n

call coordinat

ngcm=6   ! hours between GCM datasets
tseas=1.d0 * real(ngcm)*3600.d0 ! time step between data sets

!??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????
! intspin=number of initial datasets 

intmin=1

! maximum length of RCO fields
 intmax  = 1*364*24/ngcm     
! intmax  =(24*365+6)*24/ngcm    ! 1981-2004 including leap days

! intspin=trajectory release period (days*24/ngcm)
#if defined biol
intspin =6  *24/ngcm
#else
!intspin =1*365  *24/ngcm     
intspin =intmax
#endif


! trajectory run   
intrun  =50*365 *24/ngcm 
     
intstep=1   ! pos if fw & neg if bw

if(intstep.gt.0) then ! forward 
 intstart=intmin          
 intend  =intmax
elseif(intstep.lt.0) then ! backward
 intstart=intmax
 intend  =intmin
 intspin =-intspin
 intrun  =-intrun    
endif

#if defined desk
yearmin=1981
yearmax=2004
#elif defined lap
yearmin=1994
yearmax=yearmin 
#else
 stop 49567
#endif

! starting date
#if defined biol
ihour=0-6
iday=15-3           ! musslor
imon=5
iyear=2005

#else
!ihour=0-6
!iday=28           ! estonia
!imon=9
!iyear=1994

ihour=0-6  ! Finska viken overturningen
iday=1
imon=1
iyear=yearmin
#endif


! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
ncoor = 1
! kriva=0=no writing
!       1=write at time intervals of gcm datasets
!       2=write each grid-crossing and timne change
!       3=write at each iteration (all the time)
!       4=write only start and end positions
!       5=write at chosen intervals
kriva=5
! open the files for the storage of data 
! name of current trajectory run

!name='ku.00.a0'        
!name='fi.94.a0'     
!name='fi.94.n0'     
name='ms.05.t1'     

! name of directory where the files will be written

#if defined desk
 directory='/Volumes/hav4/data/rco/'
#elif defined lap
 directory='/Users/doos/data/rco12/'
#else
 stop 49567
#endif

#if defined rerun
! namep='fv.fw.aa' 
#endif

!-------- initial directions all in MODEL COORDINATES -----------------

! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
nff=sign(1,intstep)
! isec=1=meridional(y-z); 2=zonal(x-z); 3=horiz(x-y) sections; 4=fill middle of T-box
isec=4 ! musslor
!isec=1 ! Finska viken
! if idir= 1 follow positive direction (eastward/northward)  
! if idir=-1 follow negative direction (westward/southward)
! if idir=0  both directions
idir=0  ! musslor
!idir=1  ! Finska viken
!idir=-1  ! Neva
! number of trajectories can be set by
! 1=nqua: "num" which sets a constant number of trajec. in all boxes
! 2=nqua: set by chosen trajectory transport (voltr) 
! 3=nqua: set by chosen trajectory volume (voltr) 
nqua=1 ! musslor
!nqua=2
if(nqua.eq.1) then ! number of trajectories (per time resolution)
! num=NTRACMAX
 num=100
! num=1
elseif(nqua.eq.2) then ! volume transport of each trajectory in m3/s
 voltr=1000. 
elseif(nqua.eq.3) then ! volume of each trajectory in 10^6 m3
 voltr=0.
endif

!timax=10.*365.*24.*3600.  ! maximum length of a trajectory
!timax=7. *24.*3600.  ! maximum length of a trajectory
timax=21.*24.*3600.  ! maximum length of a trajectory
!timax=14.*24.*3600.  ! maximum length of a trajectory
!timax=100.*365.*24.*3600.  ! maximum length of a trajectory  ! Finska viken

! defining starting section 


      kst1=1
      kst2=km

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
       ist1=IMT-3
       ist2=ist1
       jst1=183
       jst2=jst1
! Finska viken
!       ist1=204
!       ist2=ist1
!       jst1=145
!       jst2=220
! Entire Baltic
#if defined biol
ist1=1
ist2=imt
jst1=1
jst2=jmt
kst1=km-1
kst2=km
#endif



! defining end sections


! inget stopp alls för musslorna
#if defined biol
       do n=1,NEND
        ienw(n)=IMT+1
        iene(n)=0
        jens(n)=JMT+1
        jenn(n)=0
        enddo
#else
! Lilla Bält
!       ienw(1)=1
!       iene(1)=25
!       jens(1)=47+1
!       jenn(1)=jmt
! Stora Bält
!       ienw(2)=26
!       iene(2)=36
!       jens(2)=33+1
!       jenn(2)=jmt
! Öresund 
!       ienw(1)=37
!       iene(3)=63
!       jens(3)=56+1
!       jenn(3)=jmt
! utanför Finska viken     
        ienw(1)=0
        iene(1)=204
        jens(1)=0
        jenn(1)=jmt
#endif
! iterative timesteping in seconds
#ifdef time
iter=1 ! iteration between two gcm data sets
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

!mask=-1.  ! define start section with ist1,ist2,jst1,jst2
open(21,file=directory//'kmt/kmt',form='unformatted')
!open(21,file=directory//'kmt/mask',form='unformatted')
!open(21,file=directory//'kmt/maskust',form='unformatted')
!open(21,file=directory//'kmt/maskferrysyd',form='unformatted')
!open(21,file=directory//'kmt/maskferrynord',form='unformatted')
read(21)mask
close(21)

do i=1,IMT
 do j=1,JMT
  if(mask(i,j).ne.0 .and. mask(i,j).le.4 .and. j.lt.215) mask(i,j)=-1  ! entire shallow Baltic south of 61N
 enddo
enddo

!mask=-1  ! Finska viken

if(kriva.ne.0) open(56,file=directory//'orm/traj.'//name) ! trajectory path
open(57,file=directory//'orm/traj.ut.'//name)         ! exit position
open(58,file=directory//'orm/traj.in.'//name)         ! entrence position


!??????????????????????????????????? END ???????????????????????????????

! trajectory loops 
call loop

stop
ENDPROGRAM main

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
# include "../turb.f95"
# include "readfield.f95"
# include "stat.f95"
