
PROGRAM main
  
  IMPLICIT none
#include "../param.h"
  
  common/buoyancy/tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe
  real            tmin0,tmax0,smin0,smax0,rmin0,rmax0,tmine,tmaxe,smine,smaxe,rmine,rmaxe
  
  common/tid/ints,intstart,intend,intrun,intspin,intstep,intmin,intmax
  integer ints,intstart,intend,intrun,intspin,intstep,intmin,intmax

  common/namn/name,namep,directory
  CHARACTER(LEN=30) :: name,namep
  CHARACTER(LEN=12) :: directory
  CHARACTER(LEN=40) :: initparticfile


  common /domain/ienw(NEND),iene(NEND),jens(NEND),jenn(NEND),mask(imt,jmt)
  INTEGER ienw,iene,jens,jenn,mask

  logical skriva,fileexists
  integer :: filestat,fileL,ijk
  integer nff,isec,idir,nqua,num,i,j,k
  integer ist1,ist2,jst1,jst2,kst1,kst2,ncoor
  real voltr
  real*8 tseas,tday,tyear,dtmin

  integer, dimension(:,:), ALLOCATABLE :: ijkst
  character (len=*), parameter :: ijkform = "(6I6)"
  character (len=30) :: inparg1, inparg2,inparg3
  integer :: factor
  CHARACTER(LEN=50) :: tablename


  print *,'tracmass starts'
  
  tyear=365.25d0 * 24.d0 * 3600.d0
  tday=24.d0 * 3600.d0
  
  !??????????????????? INPUT/OPTIONS BLOCK ???????????????????????????????
  
  tseas=3600*3
  
  ! -- Check if there is a time argument and if so, use it.
  if (IARGC() .eq. 0)  then
     !    intmin =(365+95+45)*8 +10    ! No 12
     intmin=intmin  +609*8+1
     ! --- Different start times for tracing the bouy...
     !     intmin =3517+590*8  ! No 12
     !     intmin =3517+589*8  ! No 12
     !     intmin =3517+589.5*8  ! No 12
     !     intmin =3517+590.5*8  ! No 12

     !  intrun =40*365*tday ! trajectory run  
     !  intrun =intspin !365*8 ! trajectory run
     intrun =8*90 ! trajectory run                 
  else 
     call getarg(1,inparg1)
     factor=1
     intmin=0
     do i=29,1,-1
        if (ichar(inparg1(i:i)) .ne. 32) then
           intmin=intmin+(ichar(inparg1(i:i))-48)*factor
           factor=factor*10
        end if
     end do
     intmin=(intmin-2)*8-1
     
     call getarg(2,inparg2)
     factor=1
     intrun=0
     do i=29,1,-1
        if (ichar(inparg2(i:i)) .ne. 32) then
           intrun=intrun+(ichar(inparg2(i:i))-48)*factor
           factor=factor*10
        end if
     end do
     intrun=180*8 !10*31      !(intrun-2)*8-1-intmin


     call getarg(3,inparg3)
     name=inparg3
  end if
  
  !Intmax =365*tday ! maximum length of OGCM fields
  intmax =intmin+35*8 ! maximum length of OGCM fields
  intspin=90*8      ! intspin=trajectory release period
  !  intspin=30*8      ! intspin=trajectory release period
  !  intrun =40*365*tday ! trajectory run 
  !  intrun =intspin !365*8 ! trajectory run 
!  intrun =5*8 ! trajectory run 
print *,intmin,intrun,intmax

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

  ! output: long/lat coordinates (ncoor=1) or model coordinates (ncoor=0) 
  ncoor = 0
  ! write trajectory (true or false)
  !      skriva = .false.
  skriva = .true.
  ! open the files for the storage of data 
  ! name of current trajectory run
  !name='go.KB.01'   
  ! name of directory where the files will be written
  directory='/Users/bror/'
  !directory='/Users/doos/data/abcdefghij/'
  print *,'directory=',directory
#if defined rerun
  !      namep='dr.bw.e4' 
#endif
  
  !-------- initial directions all in MODEL COORDINATES -----------------
  ! follow trajectories forward (ff=1), backward (ff=-1), both (ff=0) 
  nff=sign(1,intstep)
  
  ! number of trajectories can be set by
  ! 1=nqua: "num" which sets a constant number of trajec. in all boxes
  ! 2=nqua: set by chosen trajectory transport (voltr) 
  ! 3=nqua: set by chosen trajectory volume (voltr) 
  nqua=2
  if(nqua.eq.1) then     ! number of trajectories (per time resolution)
     num=36
  elseif(nqua.eq.2) then ! volume transport of each trajectory in m3/s
     voltr=100.   !50.
  elseif(nqua.eq.3) then ! volume of each trajectory in 10^6 m3
     voltr=5.
  endif

  !=====================================
  ! ===== define starting section =====
  !=====================================

  ! === Read file with startpoint deinitions if exists === 
  ! initparticfile=  'concmodis_rev.asc'
  !initparticfile=  'concmodis.asc'
  !name='concB.01'
  ! initparticfile=  'onetraj.asc'
  ! initparticfile=  'modis2005_081.txt'
  ! initparticfile=  'modis2005_086.txt'
  !initparticfile= 'randtraj.asc'
  !name='ranwB.01'
  !initparticfile='fullfield10.asc'
  !initparticfile='plumrivseed.asc'
  !name='fullB.01'
 
if (ichar(inparg1(4:4)) .eq. 32) then
   inparg1(2:6)=inparg1(1:5)
   inparg1(1:1)='0'
end if
if (ichar(inparg1(4:4)) .eq. 32) then
   inparg1(2:6)=inparg1(1:5)
   inparg1(1:1)='0'
end if

!  initparticfile='interSeed/s' // trim(name(5:7)) // '.asc' 
!  initparticfile='seedfiles/dist30seed' 
  initparticfile='seedfiles/seed' // trim(inparg1) // '.ascX'

  print *,initparticfile

  INQUIRE(FILE = initparticfile, exist=fileexists)

  if (fileexists) then
     fileL=0
     open(unit=34,file=initparticfile, ACCESS = 'SEQUENTIAL', &
          FORM = 'FORMATTED', ACTION = 'READ')
     findRecl: do
        read (unit=34, fmt=ijkform,iostat=filestat)
        if (filestat < 0) exit findRecl
        fileL=fileL+1
     end do findRecl
     allocate (ijkst(fileL,6))
     
     rewind(34)
     read_ijkst: do ijk=1,fileL
        read (unit=34, fmt=ijkform) ijkst(ijk,1), ijkst(ijk,2), &
             ijkst(ijk,3),  ijkst(ijk,4),  ijkst(ijk,5), ijkst(ijk,6)
     end do read_ijkst
 
     print *,'=== Startfile with seeding info have been used ==='
     ! === Othervise, use old method. === 
  else
     !    if idir= 1 follow positive direction (eastward/northward)  
     !    if idir=-1 follow negative direction (westward/southward)
     !    if idir=0  both directions
     idir=0
     !    isec=1=meridional; isec=2=zonal; isec=3=vertical cross section
     isec=1

!!$     ist1=10 !21            !  21
!!$     ist2=10 !21            !  21
!!$     
!!$     jst1=49 !52            !  52      !5
!!$     jst2=49 !52            !  52      !78
     
!     ist1=10 ! Merimack
!     ist2=10 !
!     jst1=49 !
!     jst2=49 !
!     name='rivmixMeriB'   
!     print *,name

     !ist1=25 ! Saco
     !ist2=25 !  
     !jst1=63 !
     !jst2=63 !
     
     ist1=41 ! Kennebec
     ist2=41 !  
     jst1=65 !
     jst2=65 !
     !name='kennB.01'   

     
     !ist1=64 ! Penobscot
     !ist2=64 !  
     !jst1=76 !
     !jst2=76 !

  

     !ist1=128 ! St John 
     !ist2=128 !  
     !jst1=71 !
     !jst2=71 !
     
     !ist1=17 ! Boundary
     !ist2=112 !  
     !jst1=25 !
     !jst2=25 !   
     !kst1=22
     !kst2=22
     
     kst1=1
     kst2=22

     allocate (ijkst((kst2-kst1+1)*(jst2-jst1+1)*(ist2-ist1+1),6))
     
     ijk=1
     do i=ist1,ist2
        do j=jst1,jst2
           do k=kst1,kst2
              ijkst(ijk,1)=i
              ijkst(ijk,2)=j
              ijkst(ijk,3)=k
              ijkst(ijk,4)=idir
              ijkst(ijk,5)=isec
              ijkst(ijk,6)=1
              ijk=ijk+1
           end do
        end do
     end do
     fileL=(kst2-kst1+1)*(jst2-jst1+1)*(ist2-ist1+1)
  End if

   
  ! iterative timesteping in seconds
#ifdef time
  dtmin=0.02*tseas
#endif

  ! water mass properties with minimum & maximum on temp, salt and density
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

  ienw(1)=163
  iene(1)=164
  jens(1)=1
  jenn(1)=jmt

 ienw(2)=1
 iene(2)=imt
 jens(2)=1
 jenn(2)=2


#ifndef mysqlwrite
 if(skriva) open(56,file=directory//'traj.'//name) ! trajectory path
 if(skriva) open(66,file=directory//'traj.'//name//'kalle') ! trajectory path
 if(skriva) open(76,action='readwrite',file=directory//'traj.rerun.'//name//'') ! trajectory path

 open(57,file=directory//'traj.ut.'//name)         ! exit position
 open(58,file=directory//'traj.in.'//name)         ! entrence position
 open(59,file=directory//'traj.er.'//name)         ! exit position
#endif 
  !??????????????????????????????????? END ???????????????????????????????

  ! trajectory loops 

!#if defined mysqlwrite
! tablename=trim(name(1:5)) // 'run'
! call createmysqltable(tablename)
!#endif

  call loop(isec,idir,nqua,num,ijkst,fileL,ncoor,nff,&
      tseas,tday,tyear,dtmin,voltr,skriva  )

  !call writepsi

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
# include "../readfieldmaine.f95"
!#include "../statocrca.f95"
!#include "../arclength.f95"
# include "../writepsi.f95"



