
module mod_print

  USE mod_param
  USE mod_name
  USE mod_time
  USE mod_loopvars
  USE mod_grid
  USE mod_buoyancy
  USE mod_seed
  USE mod_domain
  USE mod_vel
  USE mod_traj
  USE mod_write
  USE mod_pos

  ! === Selectable moules ===
  USE mod_active_particles
  USE mod_streamfunctions
  USE mod_tracer
  USE mod_sed

  implicit none

  CHARACTER(79)                       :: thinline, thickline
  
 
CONTAINS

  subroutine lines()
    thickline = "==============================================" // &
                "=============================================="
    thinline  = "----------------------------------------------" // &
                "----------------------------------------------"
  end subroutine lines


 subroutine print_header_main
    call lines
    print *, thickline!================================================= 
    print *,'             TRACMASS lagrangian off-line particle tracking '
    print *, thickline!================================================= 
  end subroutine print_header_main



  subroutine writesetup_main
    character (len=15)                           :: currDate ,currTime

    call lines
    call date_and_time(currDate, currTime)
    print *,'Start date  : '//currDate(1:4)//'-'//currDate(5:6)//'-'//currDate(7:8)
    print *,'Start time  : '//currTime(1:2)// ':'//currTime(3:4)// ':'//currTime(5:6)
    print *,'Model code  : '//trim(GCMname)
    print *,'Data surce  : '//trim(gridName)
    print *,'Run name    : '//trim(caseName)
    print *,'Description : '//trim(caseDesc)
    print *, thinline !--------------------------------------------------- 
    print *,'Directory for output files : ' ,trim(outDataDir)
    print *,'Prefix for output files    : ' ,trim(outDataFile)
    print *, thinline !--------------------------------------------------- 
    print *,"Selected compile options:"

#if defined zgrid3Dt
    print *,' - zgrid3Dt has been renamed to zgrid3D. Please update your'
    print *,'   project Makefile.prj accordingly.'
    stop
#endif
    
#ifdef timeanalyt 
    print *,' - Analytical time scheme used to solve diff. Eqs.'
#elif defined timestep
    print *,' - Stationary scheme used to solve  diff. Eqs.'
#endif
#if defined tempsalt
#if defined atmospheric
    print *,' - Temperature and humidity fields included'
#else
    print *,' - Temperature and salinity fields included'
#endif
#endif
#if defined turb
    print *, "FIX ACTPART OUTPUT"
#endif
#if defined diffusion
    print *,' - Diffusion param: Ah=',ah,'m2/s and Av=',av,'m2/s'
#if defined anisodiffusion
    print *,' - Anisotropic elliptic diffusion along the isopleths'
#endif
#endif
#if defined rerun
    print *,' - Rerun for Lagrangian stream functions'
#endif
#if defined twodim                                             
    print *,' - Two-dimensional trajectories, no change in depth'
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
#if defined atmospheric
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
    print *, thinline !--------------------------------------------------- 
  end subroutine writesetup_main

  
  subroutine print_start_loop()
    call lines
    print *, thickline!================================================= 
    write(6,FMT='(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')          &
         ' Start date in model-time     : ' , startYear, '-',  & 
         startMon, '-', startDay,' ' ,startHour, ':', startMin
  write(6,FMT='(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')          &
         ' End date in model-time       : ' , endYear, '-',  & 
         endMon, '-', endDay,' ' ,endHour, ':', endMin
    write(6,FMT='(A,I5)') ' Length of run in timesteps   : ' ,intrun
    write(6,FMT='(A,I5)') ' Number of seeding timesteps  : ' ,intspin
    write(6,FMT='(A,I5)') ' Steps between two GCM fields : ' ,iter


    !write(6,FMT='(A,I7 A,I7 A,I7)') ' intstart : ',intstart, & 
    !                                ' intspin  : ',intspin,  &
    !                                ' intrun   : ',intrun

!    print 999,intstart,intspin,intrun,intend,nff,num,voltr,&
!         tmin0,tmax0,smin0,smax0,rmin0,rmax0
  
!999 format(' intstart :',i7, '   intspin :',i7, &
!         /,'   intrun :',i7, '   intend  :',i7, &
!         /,'   nff :',i2,   
!         /,'    voltr : ',f9.0,&
!         /,'    tmin0 : ',f7.2,'  tmax0 : ',f7.2, &
!         /,'    smin0 : ',f7.2,'  smax0 : ',f7.2,&
!         /,'    rmin0 : ',f7.2,'  rmax0 : ',f7.2)
    print *, thinline !--------------------------------------------------- 
    print *,'t-step        run        out        err '  // & 
            '       tot      dt      model date'
    print *, thinline !--------------------------------------------------- 



  end subroutine print_start_loop

  subroutine print_cycle_loop()
  ! === Timing ===
    INTEGER, dimension(3)                      :: itimearray 
    INTEGER                                    :: sysrate, sysmax
    INTEGER, save                              :: currclock, lastclock=0
    REAL, dimension(2)                         :: wallarray 
    REAL                                       :: walltime, walltot
    INTEGER                                    :: wallmin, wallsec
 
#ifdef sediment
    print 599,ints,ntime,ntractot,nout,nloop,nerror,ntractot-nout, & 
         nsed,nsusp,nexit
599 format('ints=',i7,' time=',i10,' ntractot=',i8,' nout=',i8, & 
         ' nloop=',i4,' nerror=',i4,' in ocean/atm=',i8,' nsed=',i8, & 
         ' nsusp=',i8,' nexit=',9i8)
#elif defined atmospheric || ifs || rco || tes || orc || baltix || orca025  || orca025L75 || AusCOM
    print 799 ,ntime,ints ,ntractot ,nout ,nerror,ntractot-nout
799 format('ntime=',i10,' ints=',i7,' ntractot=',i8,' nout=',i8, & 
         ' nerror=',i4,' in ocean/atm=',i8)
#else
    call SYSTEM_CLOCK(currclock, sysrate, sysmax)
    if (lastclock == 0) lastclock = currclock
    walltime = (currclock - lastclock)/1000
    lastclock = currclock
    !call dtime(wallarray, walltime)
    wallmin = int(walltime/60)
    wallsec = walltime - wallmin*60

    call updateClock
    print 799 ,ints-intstart ,ntractot-nout ,nout ,nerror+nloop,ntractot, &
         wallmin, wallsec, loopYear, loopMon, loopDay, loopHour, loopMin 
799 format(i7, '|', i10,  '|', i10,  '|', i10,  '|', i10, ' | ',  &
         i3.3, ':', i2.2, ' | ', i4.4, '-', i2.2, '-', i2.2, ' ', &
         i2.2, ':', i2.2)
!799 format(i7,' run=',i10,' out=',i10,' err=',i10,' tot=',i10, ' dt=',i2.2,':',i2.2)
#endif
  end subroutine print_cycle_loop


  subroutine print_end_loop()
    character (len=15)                           :: currDate ,currTime
    call lines
    print *, thickline!================================================= 
    print *,ntractot ,  ' particles calculated'
    print *,nout     ,  ' particles exited the space and time domain'
    print *,sum(nexit), ' particles exited through the boundaries'
#ifdef sediment
    print *,nsed     ,' trajectories sedimented'
    print *,nsusp    ,' trajectories resuspended'
    call writedata(19) !end
    
#endif
#ifdef tempsalt     
    print *,nrh0,' particles outside density range'
#endif
    print *,nloop,' infinite loops'
    print *,nerror,' particles flagged with errors'
    print *,ntractot-nout-nrh0-nerror,' particles in domain'
    print *, thinline !--------------------------------------------------- 

#ifdef stremfunction     
    call writepsi
#endif    

    call date_and_time(currDate, currTime)
    print *,'End date  : '//currDate(1:4)//'-'//currDate(5:6)//'-'//currDate(7:8)
    print *,'End time  : '//currTime(1:2)// ':'//currTime(3:4)// ':'//currTime(5:6)


    print *,'The very end of TRACMASS run ',caseName

  end subroutine print_end_loop

end module mod_print
