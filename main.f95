
PROGRAM main
USE mod_param
USE mod_seed
USE mod_name
USE mod_time 
USE mod_domain
USE mod_buoyancy

IMPLICIT none

INTEGER i,j,n


call init_params
call coordinat
call writesetup
call init_seed


do i=1,ijkMax
   print *,ijkst(i,:)
end do
stop 666

tseas=1.d0 * real(ngcm)*3600.d0 ! time step between data sets

modrundirCond: if(intstep.gt.0) then ! forward 
   intstart=intmin          
   intend  =intmax
elseif(intstep.lt.0) then ! backward
   intstart=intmax
   intend  =intmin
   intspin =-intspin
   intrun  =-intrun    
end if modrundirCond

if(nqua.eq.1) then ! number of trajectories (per time resolution)
   ! num=NTRACMAX
   num=partQuant
elseif(nqua.eq.2) then 
   voltr=partQuant 
elseif(nqua.eq.3) then 
   voltr=partQuant
endif

if(kriva.ne.0) open(56,file='/Users/bror/orm/traj.'//name) ! trajectory path
open(57,file='/Users/bror/orm/traj.ut.'//name)             ! exit position
open(58,file='/Users/bror/orm/traj.in.'//name)             ! entrence position


!??????????????????????????????????? END ???????????????????????????????

! trajectory loops 

call loop











































CONTAINS

  subroutine writesetup
    character (len=15)                           :: currDate ,currTime
    
    call date_and_time(currDate, currTime)
    
    print *,'======================================================'
    print *,'=== TRACMASS lagrangian off-line particle tracking ==='
    print *,'------------------------------------------------------'
    print *,'Start date: '//currDate(1:4)//'-'//currDate(5:6)//'-'//currDate(7:8)
    print *,'Start time: '//currTime(1:2)// ':'//currTime(3:4)// ':'//currTime(5:6)
    print *,'Model code: '//trim(GCMname)
    print *,'Data surce: '//trim(gridName)
    print *,'Run name  : '//trim(caseName)
    print *,'          : '//trim(caseDesc)
    print *,'------------------------------------------------------'
#if defined tempsalt
    print *,'with temperature and salinity fields'
#endif
#if defined turb
    print *,'with sub-grid turbulence parameterisation'
#endif
#if defined rerun
    print *,'Rerun in order to store the Lagrangian stream functions in the different basins'
#endif
#if defined twodim                                             
    print *,'Two-dimensional trajectory, which do not change depth'
#endif
#if defined streamxy
    print *,'Lagrangian horizontal stream function stored'
#endif
    
#if defined streamv
    print *,'Lagrangian vertical depth stream function stored'
#endif
    
#if defined streamr
#if defined streamts
    print *,'Lagrangian density, temperature and salinity stream function stored'
#else
    print *,'Lagrangian density stream function stored'
#endif
#endif
    
#if defined streamxy
    print *,'Lagrangian horizontal stream function stored'
#endif
    
#if defined tracer
    print *,'Lagrangian trajectory particle tracer stored'
#endif
  end subroutine writesetup
end PROGRAM main

!______________ END OF MAIN PROGRAM _______________________________

!!$#if defined orc
!!$print *,'ORCA GCM fields'
!!$#elif defined rco
!!$print *,'RCO GCM fields'
!!$#elif defined tes
!!$print *,'ACADEMIC TEST fields'
!!$#elif defined sim
!!$print *,'Simpevarp GCM fields'
!!$#elif defined fors
!!$print *,'Forsmark GCM fields'
!!$#elif defined occam66
!!$print *,'OCCAM 1/4 deg GCM fields'
!!$#elif defined occam083
!!$print *,'OCCAM 1/12 deg GCM fields'
!!$#elif defined sigma
!!$print *,'OGCM with sigma coordinate fields'
!!$#elif defined atm
!!$print *,'IFS (AGCM) with atmospheric sigma coordinate fields'
!!$#endif
