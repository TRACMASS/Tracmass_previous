SUBROUTINE init_params
! Loads parameters from the namelists (projname).in and (casename).in.
! Allocates matrices and set time variables.

   USE mod_param
   USE mod_seed
   USE mod_grid
   USE mod_name
   USE mod_time 
   USE mod_domain
   USE mod_vel
   USE mod_traj
 !  USE mod_dens
   USE mod_buoyancy
   USE mod_streamfunctions
   USE mod_tracer
   USE mod_getfile
   USE mod_write
   
#if defined diffusion || turb 
   USE mod_diffusion
#endif
#ifdef sediment
   USE mod_orbital
   USE mod_sed
#endif
   IMPLICIT NONE

!!----------------------------------------------------------------------------
   
   INTEGER                                    ::  argint1 ,argint2
   INTEGER                                    ::  dummy ,factor ,i ,dtstep
   INTEGER                                    ::  gridVerNum ,runVerNum
   CHARACTER (LEN=30)                         ::  inparg, argname
   CHARACTER (LEN=200)                        ::  projdir, ormdir
   real*8                                     :: jd

! Setup namelists
   namelist /INIT_NAMELIST_VERSION/ gridVerNum
   namelist /INIT_GRID_DESCRIPTION/ GCMname, GCMsource, gridName, gridSource,&
                                    gridDesc, inDataDir
   namelist /INIT_CASE_DESCRIPTION/ caseName, caseDesc
   namelist /INIT_GRID_SIZE/        imt, jmt, km, nst, subGrid, subGridImin, &
                                    subGridImax, subGridJmin, subGridJmax,   &
                                    SubGridFile, subGridID
   namelist /INIT_BASE_TIME/        baseSec, baseMin, baseHour, baseDay,     &
                                    baseMon, baseYear
   namelist /INIT_GRID_TIME/        fieldsPerFile, ngcm, iter, intmax,       &
                                    minvelJD, maxvelJD
   namelist /INIT_START_DATE/       startSec, startMin, startHour,           & 
                                    startDay, startMon, startYear,           &
                                    startJD, intmin
   namelist /INIT_RUN_TIME/         intspin, intrun
   namelist /INIT_WRITE_TRAJS/      twritetype, kriva, outDataDir,           &
                                    outDataFile, intminInOutFile
   namelist /INIT_SEEDING/          nff, isec, idir, nqua, partQuant,        &
                                    ntracmax, loneparticle, SeedType, ist1,  &
                                    ist2, jst1, jst2, kst1, kst2, tst1, tst2,&
                                    seedDir, seedFile, varSeedFile, seedTime,&
                                    seedAll, seedPos, seedparts, seedpart_id,&
                                    seedsubints
   namelist /INIT_KILLZONES/        nend, ienw, iene, jens, jenn, timax
   namelist /INIT_TEMP_SALT/        tmin0, tmax0, smin0, smax0, rmin0, rmax0,&
                                    tmine, tmaxe, smine, smaxe, rmine, rmaxe
#if defined diffusion || turb 
   namelist /INIT_DIFFUSION/        ah, av
#endif
#ifdef sediment
   namelist /INIT_SEDIMENT/         partdiam, rhos, cwamp, twave, critvel
#endif

!!--------------------------------------------------------------------------  
   
   Project  = PROJECT_NAME
   Case     = CASE_NAME

   IF ((IARGC() > 0) )  THEN
      CALL getarg(1,Case)
   END IF

   CALL getenv('TRMPROJDIR',projdir)
   if (len(trim(projdir)) == 0) then
      CALL getenv('TRMDIR',ormdir)
      if (len(trim(ormdir)) .ne. 0) then
         projdir = trim(ormdir)//'/'//'projects/'//trim(Project)//'/'
      else
         projdir = 'projects/'//trim(Project)
      end if
   end if

   OPEN (8,file=trim(projdir)//'/'//trim(Project)//'.in',    &
        & status='OLD', delim='APOSTROPHE')
   ! -- Check if the namefiles has correct version number. 
   READ (8,nml=INIT_NAMELIST_VERSION)
   IF (gridVerNum < 6) THEN
      PRINT *,'                     ERROR                     '
      PRINT *,'Your namefile out of date. The latest version is described at:'
      PRINT *,'http://docs.tracmass.org/namelist.html'
      PRINT *,'Change gridVerNum to 6 when done.'
      STOP
   END IF
   READ (8,nml=INIT_CASE_DESCRIPTION)
   READ (8,nml=INIT_GRID_SIZE)
   READ (8,nml=INIT_BASE_TIME)
   READ (8,nml=INIT_GRID_TIME)
   READ (8,nml=INIT_START_DATE)
   READ (8,nml=INIT_RUN_TIME)
   READ (8,nml=INIT_WRITE_TRAJS)
   READ (8,nml=INIT_SEEDING)
   READ (8,nml=INIT_KILLZONES)
   READ (8,nml=INIT_TEMP_SALT)
#if defined diffusion || turb 
   READ (8,nml=INIT_DIFFUSION)
#endif
#ifdef sediment
   READ (8,nml=INIT_SEDIMENT)   
#endif
   CLOSE (8)

   print *,'Run file    : ',trim(projdir)//'/'//trim(Case)//'.in'
   OPEN (8,file=trim(projdir)//'/'//trim(Case)//'.in',     &
        & status='OLD', delim='APOSTROPHE')
   READ (8,nml=INIT_NAMELIST_VERSION)
   IF (gridVerNum < 6) THEN
      PRINT *,'                     ERROR                     '
      PRINT *,'Your namefile out of date. The latest version is described at:'
      PRINT *,'http://docs.tracmass.org/namelist.html'
      PRINT *,'Change gridVerNum to 6 when done.'
      STOP
   END IF
   READ (8,nml=INIT_CASE_DESCRIPTION)
   READ (8,nml=INIT_GRID_SIZE)
   READ (8,nml=INIT_BASE_TIME)
   READ (8,nml=INIT_GRID_TIME)
   READ (8,nml=INIT_START_DATE)
   READ (8,nml=INIT_RUN_TIME)
   READ (8,nml=INIT_WRITE_TRAJS)
   READ (8,nml=INIT_SEEDING)
   READ (8,nml=INIT_KILLZONES)
   READ (8,nml=INIT_TEMP_SALT)
#if defined diffusion || turb 
   READ (8,nml=INIT_DIFFUSION)
#endif
#ifdef sediment
   READ (8,nml=INIT_SEDIMENT)   
#endif
   CLOSE (8)

   SELECT CASE (subGrid)
   CASE (0)          
      PRINT *,'Sub-grid    : Use the Full grid.'     
      subGridImin =   1 
      subGridJmin =   1
      subGridKmin =   1
      subGridImax = imt
      subGridJmax = jmt 
      subGridKmax = km 
   CASE (1)
      PRINT *,'Sub-grid    : ', subGridImin ,subGridImax, &
           &   subGridJmin ,subGridJmax
      imt = subGridImax-subGridImin+1
      jmt = subGridJmax-subGridJmin+1
#if !defined(explicit_w) && !defined(twodim)
      if ((subGridKmax-subGridKmin+1) < km) then
         print *, 'ERROR!'
         print *, 'subGridKmin and subGridKmax requires -Dtwodim  or -Dexplicit_w'
         print *, 'to be selected in the project Makefile.'
         stop
      end if
#endif
      km  = subGridKmax-subGridKmin+1
   CASE default
      PRINT *,'==================== ERROR ===================='
      PRINT *,'This subGrid selection is not implemented yet.'
      PRINT *,'subGrid = ' ,subGrid
      STOP
   END SELECT
   start1d  = [subGridKmin]
   count1d  = [subGridKmax]
   start2d  = [1, 1,           subGridImin, subGridJmin]
   count2d  = [1, 1,           imt,         jmt        ]
   start3d  = [1, subGridImin, subGridJmin, subGridKmin]
   count3d  = [1, imt,         jmt,         km         ]
   
   IF ((IARGC() > 1) )  THEN
      ARG_INT1 = 0.1
      CALL getarg(2,inparg)
      if ( ARG_INT1 == 0) then
         read( inparg, '(i15)' ) ARG_INT1
         write( inargstr1, '(A,i9.9 )' ) '_a',ARG_INT1
      else
         read( inparg, '(f15.10)' ) ARG_INT1
         write( inargstr1, '(A,i9.9 )' ) '_a',int(ARG_INT1)
      end if
   END IF
   
   IF ((IARGC() > 2) ) THEN
      ARG_INT2 = 0.1
      CALL getarg(3,inparg)
      if ( ARG_INT2 == 0) then
         read( inparg, '(i15)' ) ARG_INT2
         write( inargstr2, '(A,i9.9)' ) '_b',ARG_INT2
      else
         read( inparg, '(f15.10)' ) ARG_INT2
         write( inargstr2, '(A,i9.9)' ) '_b',int(ARG_INT2)
      end if
   END IF
   
   timax    =  24.*3600.*timax ! convert time lengths from days to seconds
   dstep    =  1.d0/dble(iter)
   dtmin    =  dstep * tseas
   baseJD   =  jdate(baseYear  ,baseMon  ,baseDay)  + &  
           ( dble((baseHour)*3600 + baseMin*60 + baseSec) / 86400 )
   if (startJD < 1) then
      startJD  =  jdate(startYear ,startMon ,startDay) + 1 + &  
           ( dble((startHour)*3600 + startMin*60 + startSec) / 86400 ) -baseJD
   end if

   if (nff == 1) then
      intmin = int(floor((startJD-1)/(real(ngcm)/24.) + 1))
   else
      intmin = int(ceiling((startJD-1)/(real(ngcm)/24.) + 1))
   end if

   if (endJD < 1) then
      endJD  =  jdate(endYear ,endMon ,endDay) + 1 + &  
           ( dble((endHour)*3600 + endMin*60 + endSec) / 86400 ) -baseJD
   end if
   if (endJD < startJD) then
      endJD =  baseJD + startJD + intrun*ngcm/24. -2
   end if
   call  gdate (endJD ,endYear , endMon ,endDay)
   endFrac = (endJD-int(endJD))*24
   endHour = int(endFrac)
   endFrac = (endFrac - endHour) * 60
   EndMin  = int(endFrac)
   endSec  = int((endFrac - currMin) * 60)

   if (nff == 1) then
      intmax = int(floor((endJD-1)/(real(ngcm)/24.) + 1))
   else
      intmax = int(ceiling((endJD-1)/(real(ngcm)/24.) + 1))
   end if
   
   if (maxvelJD > 0) then
      minvelints = (minvelJD)/(real(ngcm)/24.)+1
      maxvelints = (maxvelJD)/(real(ngcm)/24.)+1
      intmax = maxvelints - intmin
   end if

   tseas= dble(ngcm)*3600.d0

   ! --- ist -1 to imt ---
   IF ( ist1 == -1) THEN 
      ist1 = IMT
   END IF
   IF ( ist2 == -1) THEN 
      ist2 = IMT
   END IF
   ! --- jst -1 to jmt ---
   IF ( jst1 == -1) THEN
      jst1=jmt
   END IF
   IF ( jst2 == -1) THEN
      jst2=jmt
   END IF 
   ! --- kst -1 to km ---
   IF ( kst1 == -1) THEN
      kst1 = KM
   END IF
   IF ( kst2 == -1) THEN 
      kst2 = KM
   END IF

   if (len(trim(inDataDir)) == 0) then
      CALL getenv('TRMINDATADIR', projdir)
      if (len(trim(projdir)) .ne. 0) then
         inDataDir = trim(projdir) // trim(Project) // '/'
      end if
   end if

   if (len(trim(outDataDir)) == 0) then
      CALL getenv('TRMOUTDATADIR', projdir)
      if (len(trim(projdir)) .ne. 0) then
         outDataDir = trim(projdir) // trim(Project) // '/'
      end if
   end if

   if (outDataFile == '')  outdataFile = Case

!!---------------------------------------------------------------------------
!!------------------------ A L L O C A T I O N S ----------------------------
!!---------------------------------------------------------------------------

      ! --- Allocate information about the coordinates and grid ---

      ALLOCATE ( csu (0:jmt), cst(jmt)  ) 
      ALLOCATE ( phi(0:jmt),   zw(0:km) ) 
      ALLOCATE ( dyt(jmt), dxv(imt+2,jmt), dyu(imt+2,jmt) ) 
      ALLOCATE ( mask(imt,jmt) )
      dyt = 0
      dxv = 0
      dyu = 0

#ifdef zgrid3Dt
      ALLOCATE ( dzt(imt,jmt,km,nst) )   
#elif  zgrid3D
      ALLOCATE ( dzt(imt,jmt,km) )   
#endif /*zgrid3Dt*/
#ifdef varbottombox
      ALLOCATE ( dztb(imt,jmt,nst) )  !should probably be changed to  dztb(imt,jmt)
#endif /*varbottombox*/
      ALLOCATE ( dxdy(imt,jmt) )   
      ALLOCATE ( kmt(imt,jmt), dz(km) )
    
      ! --- Allocate velocity fields, temperature, salinity, density, --- 
      ! --- sea-surface height, and trajectory data                   ---
      ALLOCATE ( uflux(imt,jmt,km,nst), vflux(imt,0:jmt,km,nst) )
      ALLOCATE ( hs(imt+1,jmt+1,nst) )
#ifdef explicit_w || full_wflux
      ALLOCATE ( wflux(imt+2 ,jmt+2 ,0:km,NST) )
#else
      ALLOCATE ( wflux(0:km,NST) )
#endif
      hs    = 0.
      uflux = 0.
      vflux = 0.
      wflux = 0.d0
      ALLOCATE ( uvel(imt+2,jmt,km) ,vvel(imt+2,jmt,km) ,wvel(imt+2,jmt,km) )
      
      ! === Init mod_traj ===
      ALLOCATE ( trj(NTRJ,ntracmax), nrj(NNRJ,ntracmax) )
      ALLOCATE ( nexit(NEND) ) 
      nrj = 0
      trj = 0.d0
      nexit = 0
      ntractot = 0
      numseedsubints = max(count(seedsubints /= -1), 1)
      if (nqua == 5 .AND. seedsubints(1)==-1) then
         print *,  "Error! "
         print *,  "At least one element in numseedsubints must be " // &
                   "given when nqua=5 is used." 
         stop
      elseif  (nqua /= 5) then
         seedsubints(1) = 0
      end if

#ifdef tempsalt
      ALLOCATE ( tem(imt,jmt,km,nst) ) 
      ALLOCATE ( sal(imt,jmt,km,nst) )
      ALLOCATE ( rho(imt,jmt,km,nst) )
      tem = 0.
      sal = 0.
      rho = 0.
#endif

      ! --- Allocate Lagrangian stream functions ---
#ifdef streamxy
      ALLOCATE ( stxyy(imt,jmt,nend), stxyx(imt,jmt,nend) )
      stxyy=0.
      stxyx=0.
#endif
#ifdef streamv
      ALLOCATE ( stxz(imt,km,nend), styz(jmt,km,nend) )
      stxz=0.
      styz=0.
#endif
#ifdef streamr
      ALLOCATE ( stxr(imt,mr,nend,lov), styr(jmt,mr,nend,lov) )
      stxr=0.
      styr=0
#endif
#ifdef stream_thermohaline
      ALLOCATE ( psi_ts(MR,MR,2,nend) )
      psi_ts=0.
#endif

      ! --- Allocate tracer data ---
#ifdef tracer
      ALLOCATE ( tra(imt,jmt,km) )
      tra=0.
#endif

      ! --- Allocate sedimentation data ---
#ifdef sediment
      ALLOCATE (orb(km) )
      nsed = 0
      nsusp = 0
#endif

END SUBROUTINE init_params

!!----------------------------------------------------------------------------
!!----------------------------------------------------------------------------
