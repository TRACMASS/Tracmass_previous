
PROGRAM main
USE mod_param
USE mod_seed
USE mod_name
USE mod_time 
USE mod_domain
USE mod_buoyancy
#ifdef diffusion
USE mod_diffusion
#endif

IMPLICIT none

INTEGER                                    :: i,j,n
CHARACTER(LEN=200)                         :: fullWritePref
CHARACTER(LEN=8)                           :: WriteStamp

call init_params
call coordinat
call writesetup

!do i=1,ijkMax
!   print *,ijkst(i,:)
!end do
!stop 666

modrundirCond: if(intstep.gt.0) then ! forward 
   intstart=intmin          
   intend  =intmax
elseif(intstep.lt.0) then ! backward
   intstart=intmax
   intend  =intmin
   intspin =-intspin
   intrun  =-intrun    
end if modrundirCond

call init_seed

if(nqua.eq.1) then ! number of trajectories (per time resolution)
   ! num=NTRACMAX
   num=partQuant
elseif(nqua.eq.2) then 
   voltr=partQuant 
elseif(nqua.eq.3) then 
   voltr=partQuant
endif

writeStamp='00000000'
write (writeStamp,'(i8.8)') intstart



if (intminInOutFile.eq.1) then
   fullWritePref =  trim(outDataDir)//trim(outDataFile)//writeStamp
else
   fullWritePref =  trim(outDataDir)//trim(outDataFile)
end if

#if defined textwrite
open(56,file=trim(fullWritePref)//'_run.asc')       ! trajectory path
open(57,file=trim(fullWritePref)//'_out.asc')       ! exit position
open(58,file=trim(fullWritePref)//'_in.asc')        ! entrance position
#endif

#if defined binwrite
open(unit=76 ,file=trim(fullWritePref)//'_run.bin' &  ! Trajectory path
     ,access='direct' ,form='unformatted' ,recl=20 ,status='replace')   !
open(unit=75 ,file=trim(fullWritePref)//'_out.bin' &  ! Exit position
     ,access='direct' ,form='unformatted' ,recl=20 ,status='replace')   !
open(unit=77 ,file=trim(fullWritePref)//'_kll.bin' &  ! Killed position
     ,access='direct' ,form='unformatted' ,recl=20 ,status='replace')   !
open(unit=78 ,file=trim(fullWritePref)//'_in.bin' &  ! Entrance position
     ,access='direct' ,form='unformatted' ,recl=20 ,status='replace')   !
open(unit=79 ,file=trim(fullWritePref)//'_err.bin' &  ! Error position 
     ,access='direct' ,form='unformatted' ,recl=20 ,status='replace')   ! 
#endif



! === Start main loop ===

call loop

CONTAINS

  subroutine writesetup
    character (len=15)                           :: currDate ,currTime
    
    call date_and_time(currDate, currTime)
    
    print *,'======================================================'
    print *,'=== TRACMASS lagrangian off-line particle tracking ==='
    print *,'------------------------------------------------------'
    print *,'Start date  : '//currDate(1:4)//'-'//currDate(5:6)//'-'//currDate(7:8)
    print *,'Start time  : '//currTime(1:2)// ':'//currTime(3:4)// ':'//currTime(5:6)
    print *,'Model code  : '//trim(GCMname)
    print *,'Data surce  : '//trim(gridName)
    print *,'Run name    : '//trim(caseName)
    print *,'Description : '//trim(caseDesc)
    print *,'------------------------------------------------------'
#ifdef timeanalyt 
    print *,'Analytical time scheme used to solve the differential Eqs.'
#elif defined timestep
    print *,'Time steps with analytical stationary scheme used to solve the differential Eqs.'
#endif
#if defined tempsalt
#if defined ifs
    print *,'Temperature and humidity fields included'
#else
    print *,'Temperature and salinity fields included'
#endif
#endif
#if defined turb
    print *,'with sub-grid turbulence parameterisation'
#elif defined diffusion
    print *,'with diffusion parameterisation, Ah=',ah,'m2/s and Av=',av,'m2/s'
#endif
#if defined rerun
    print *,'Rerun in order to store the Lagrangian stream functions in the different basins'
#endif
#if defined twodim                                             
    print *,'Two-dimensional trajectories, which do not change depth'
#endif
#if defined full_wflux
    print *,' 3D vertival volume flux field.'
#endif
#if defined explicit_w
    print *,'Given vertical velocity.'
#endif
#if defined sediment                                             
    print *,'Sedimentation including resuspension activated'
#endif


#if defined streamxy
    print *,'Lagrangian horizontal stream function stored'
#endif
    
#if defined streamv
    print *,'Lagrangian vertical depth stream function stored'
#endif
    
#if defined streamr
#if defined streamts
#if defined ifs
    print *,'Lagrangian density, temperature and humidity stream function stored'
#else
    print *,'Lagrangian density, temperature and salinity stream function stored'
#endif
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
