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

  modrundirCond: if(nff == 1) then ! forward 
     intstart =  intmin          
     intend   =  intmax
  elseif(nff == -1) then ! backward
     intstart =  intmin+intrun
     intend   =  intmin
     intspin  = -intspin
     intrun   = -intrun    
  end if modrundirCond
  call setupgrid
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
  open(59,file=trim(fullWritePref)//'_err.asc')       ! Error position
#endif
  
#if defined binwrite
  ! Trajectory path
  open(unit=76 ,file=trim(fullWritePref)//'_run.bin' &  
       ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
  ! Exit position 
  open(unit=75 ,file=trim(fullWritePref)//'_out.bin' &  
       ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
  ! Killed position
  open(unit=77 ,file=trim(fullWritePref)//'_kll.bin' &
       ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
  ! Entrance position
  open(unit=78 ,file=trim(fullWritePref)//'_in.bin' &  
       ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
  ! Error position
  open(unit=79 ,file=trim(fullWritePref)//'_err.bin' &  
       ,access='direct' ,form='unformatted' ,recl=24 ,status='replace')
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
#elif defined timestat
    print *,'Steady state velocity fields with analytical stationary scheme used to solve the differential Eqs.'
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
#endif
#if defined diffusion
    print *,'with diffusion parameterisation, Ah=',ah,'m2/s and Av=',av,'m2/s'
#if defined anisodiffusion
    print *,'with anisotropic elliptic diffusion along the isopleths'
#endif
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

#if defined stream_thermohaline
    print *,'Lagrangian thermohaline stream function stored'
#endif
    
#if defined tracer
    print *,'Lagrangian trajectory particle tracer stored'
#endif

  end subroutine writesetup
end PROGRAM main


