PROGRAM main
 
  USE mod_param
  USE mod_seed
  USE mod_name
  USE mod_time 
  USE mod_write
  USE mod_domain
  USE mod_buoyancy
#ifdef diffusion
  USE mod_diffusion
#endif
  
  IMPLICIT none
  INTEGER                                    :: i,j,n

  print *,'==================================================================='
  print *,'========= TRACMASS lagrangian off-line particle tracking =========='
  print *,'==================================================================='
  call init_params
  call coordinat
  call writesetup

  modrundirCond: if(intstep.gt.0) then ! forward 
     intstart =  intmin          
     intend   =  intmax
  elseif(intstep.lt.0) then ! backward
     intstart =  intmin+intrun
     intend   =  intmin
     intspin  = -intspin
     intrun   = -intrun    
  end if modrundirCond
 
 call setupgrid
  if (minval(dxv) < 0) then
     print *, " "
     print *, " === Error! === "
     print *, "The array dxv contains negative values."
     print *, "Please check your setupgrid.f95 file."
     stop
  end if
  if (minval(dyu) < 0) then
     print *, " "
     print *, " === Error! === "
     print *, "The array dyu contains negative values."
     print *, "Please check your setupgrid.f95 file."
     stop
  end if

  call init_seed

  if(nqua.eq.1) then ! number of trajectories (per time resolution)
     ! num=NTRACMAX
     num=partQuant
  elseif(nqua.eq.2) then 
     voltr=partQuant 
  elseif(nqua.eq.3) then 
     voltr=partQuant
  endif
  
  call open_outfiles
  call loop
  call close_outfiles

  return

CONTAINS
  
  subroutine writesetup
    character (len=15)                           :: currDate ,currTime
    
    call date_and_time(currDate, currTime)
    
  
    print *,'Start date  : '//currDate(1:4)//'-'//currDate(5:6)//'-'//currDate(7:8)
    print *,'Start time  : '//currTime(1:2)// ':'//currTime(3:4)// ':'//currTime(5:6)
    print *,'Model code  : '//trim(GCMname)
    print *,'Data surce  : '//trim(gridName)
    print *,'Run name    : '//trim(caseName)
    print *,'Description : '//trim(caseDesc)
    print *,'-------------------------------------------------------------------'
    print *,"= Selected compile options:"
#ifdef timeanalyt 
    print *,' - Analytical time scheme used to solve the differential Eqs.'
#elif defined timestep
    print *,' - Time steps with analytical stationary scheme used to solve the differential Eqs.'
#endif
#if defined tempsalt
#if defined ifs
    print *,' - Temperature and humidity fields included'
#else
    print *,' - Temperature and salinity fields included'
#endif
#endif
#if defined turb
    print *,' - Sub-grid turbulence parameterisation'
#endif
#if defined diffusion
    print *,' - Diffusion parameterisation, Ah=',ah,'m2/s and Av=',av,'m2/s'
#if defined anisodiffusion
    print *,' - Anisotropic elliptic diffusion along the isopleths'
#endif
#endif
#if defined rerun
    print *,' - Rerun in order to store the Lagrangian stream functions in the different basins'
#endif
#if defined twodim                                             
    print *,' - Two-dimensional trajectories, which do not change depth'
#endif
#if defined full_wflux
    print *,' - 3D vertival volume flux field.'
#endif
#if defined explicit_w
    print *,' - Explicit vertical velocities from the GCM.'
#endif
#if defined sediment                                             
    print *,' - Sedimentation including resuspension activated'
#endif


#if defined streamxy
    print *,' - Lagrangian horizontal stream function stored'
#endif
    
#if defined streamv
    print *,' - Lagrangian vertical depth stream function stored'
#endif
    
#if defined streamr
#if defined streamts
#if defined ifs
    print *,' - Lagrangian density, temperature and humidity stream function stored'
#else
    print *,' - Lagrangian density, temperature and salinity stream function stored'
#endif
#else
    print *,' - Lagrangian density stream function stored'
#endif
#endif

#if defined stream_thermohaline
    print *,' - Lagrangian thermohaline stream function stored'
#endif
    
#if defined tracer
    print *,' - Lagrangian trajectory particle tracer stored'
#endif
    print *,'-------------------------------------------------------------------'
  end subroutine writesetup
end PROGRAM main


