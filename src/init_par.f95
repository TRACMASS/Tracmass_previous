SUBROUTINE init_params
!!----------------------------------------------------------------------------
!!
!!
!!       SUBROUTINE: init_params
!!
!!          Loads parameters from the namelists run.in and grid.in
!!          Allocates matrices for trajectory data, GCM fields,
!!          and Lagrangian stream functions.
!!
!!
!!
!!
!!----------------------------------------------------------------------------
   USE mod_param
   USE mod_seed
   USE mod_coord
   USE mod_grid
   USE mod_name
   USE mod_time 
   USE mod_domain
   USE mod_vel
   USE mod_traj
   USE mod_dens
   USE mod_buoyancy
   USE mod_streamxy
   USE mod_streamv
   USE mod_streamr
   USE mod_stream_thermohaline
   USE mod_tracer
   USE mod_getfile
   
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
   CHARACTER (LEN=23)                         ::  Project, Case
   CHARACTER (LEN=200)                        ::  projdir, ormdir

   real*8                                     :: jd

!!----------------------------------------------------------------------------
!!-------------------- R E A D   F R O M   N A M E L I S T S -----------------
!!----------------------------------------------------------------------------
   ! -------------------------------
   ! --- Parameters from grid.in ---
   ! -------------------------------
   namelist /INITGRIDVER/    gridVerNum
   namelist /INITGRIDDESC/   GCMname, GCMsource, gridName, gridSource, gridDesc
   namelist /INITGRIDGRID/   IMT, JMT, KM, LBT, NEND
   namelist /INITGRIDNTRAC/  NTRACMAX
   namelist /INITGRIDDATE/   yearmin, yearmax, baseSec  ,baseMin  ,baseHour,   &
                         &   baseDay  ,baseMon  ,baseYear
   namelist /INITGRIDTIME/   ngcm, iter, intmax ,fieldsPerFile
   
   ! ------------------------------
   ! --- Parameters from run.in ---
   ! ------------------------------
   namelist /INITRUNVER/     runVerNum
   namelist /INITRUNGRID/    subGrid ,subGridImin ,subGridImax ,subGridJmin,   &
                         &   subGridJmax ,SubGridFile, subGridID
   namelist /INITRUNTIME/    intmin, intspin, intrun, intstep 
   namelist /INITRUNDATE/    startSec ,startMin ,startHour,                    &
                         &   startDay ,startMon ,startYear,                    &
                         &   ihour, iday, imon, iyear
   namelist /INITRUNWRITE/   ncoor, twritetype, kriva,                                     &
                         &   inDataDir ,outDataDir, topoDataDir,               &
                             outDataFile ,intminInOutFile
   namelist /INITRUNSEED/    nff, isec, idir, nqua, partQuant,                 &
                         &   seedType, seedPos, seedTime, seedAll,             &
                         &   ist1, ist2, jst1, jst2, kst1, kst2, tst1, tst2,   &
                             varSeedFile, seedDir, seedFile, timeFile
   namelist /INITRUNDESC/    caseName, caseDesc  
#ifdef tempsalt
   namelist /INITRUNTEMPSALT/ tmin0, tmax0, smin0, smax0, rmin0, rmax0, &
                         &    tmine, tmaxe, smine, smaxe, rmine, rmaxe
#endif
#if defined diffusion || turb 
   namelist /INITRUNDIFFUSION/ ah, av
#endif
#ifdef sediment
   namelist /INITRUNSEDIMENT/  partdiam, rhos, cwamp, twave, critvel
#endif
   namelist /INITRUNEND/ ienw, iene, jens, jenn, timax

!!--------------------------------------------------------------------------  
   
   Project  = PROJECT_NAME
   Case     = CASE_NAME
   
   IF ( (IARGC() == 1 ) .OR. (IARGC() == 4 ) )  then
      CALL getarg(IARGC(),inparg)
      Case = inparg
   END IF

   CALL getenv('ORMPROJDIR',projdir)
   if (len(trim(projdir)) == 0) then
      CALL getenv('ORMDIR',ormdir)
      if (len(trim(ormdir)) .ne. 0) then
         projdir = trim(ormdir)//'/'//'projects/'//trim(Project)//'/'
      else
         projdir = 'projects/'//trim(Project)
      end if
   end if

   OPEN (8,file=trim(projdir)//'/'//trim(Project)//'_grid.in',    &
       & status='OLD', delim='APOSTROPHE')
   
      ! -- Check if the namefiles has correct version number. 
      READ (8,nml=INITGRIDVER)
         IF (gridVerNum < 1) THEN
            PRINT *,'==================== ERROR ===================='
            PRINT *,'Your grid namefile seems to be out of date.'
            PRINT *,'Check the version_grid.txt file for changes.'
            PRINT *,'You have to edit the version number in you grid'
            PRINT *,'manually when done.'
            STOP
         END IF
      READ (8,nml=INITGRIDDESC)
      READ (8,nml=INITGRIDGRID)
      READ (8,nml=INITGRIDNTRAC)
      READ (8,nml=INITGRIDTIME)
      READ (8,nml=INITGRIDDATE)
   
   CLOSE (8)

   print *,' runfile =  ',trim(projdir)//'/'//trim(Case)//'_run.in'

   OPEN (8,file=trim(projdir)//'/'//trim(Case)//'_run.in',     &
        & status='OLD', delim='APOSTROPHE')
   READ (8,nml=INITRUNDESC)
   READ (8,nml=INITRUNGRID)
   SELECT CASE (subGrid)
   CASE (0)          
      PRINT *,'Use the Full grid.'     
      subGridImin =   1
      subGridJmin =   1
      subGridImax = imt
      subGridJmax = jmt 
   CASE (1)
      PRINT *,'Use a subgrid: ', subGridImin ,subGridImax, &
           &   subGridJmin ,subGridJmax
      imt = subGridImax-subGridImin+1
      jmt = subGridJmax-subGridJmin+1
   CASE default
      PRINT *,'==================== ERROR ===================='
      PRINT *,'This subGrid selection is not implemented yet.'
      PRINT *,'subGrid = ' ,subGrid
      STOP
   END SELECT
         start1d  = [  1]
         count1d  = [ km]
         start2d  = [  1 ,  1 ,subGridImin ,subGridJmin]
         count2d  = [  1 ,  1 ,subGridImax ,subGridJmax]
         start3d  = [  1, subGridImin, subGridJmin,  1]
         count3d  = [  1, subGridImax, subGridJmax, km]

         READ (8,nml=INITRUNTIME)
         READ (8,nml=INITRUNDATE)
         READ (8,nml=INITRUNWRITE)  
         READ (8,nml=INITRUNSEED)
#ifdef tempsalt
         READ (8,nml=INITRUNTEMPSALT)
#endif
#if defined diffusion || turb 
         READ (8,nml=INITRUNDIFFUSION)
#endif
#ifdef sediment
         READ (8,nml=INITRUNSEDIMENT)
#endif
         READ(8,nml=INITRUNEND)
         
      CLOSE (8)

      timax    =  24.*3600.*timax ! convert time lengths from days to seconds
      dstep    =  1.d0/dble(iter)
      dtmin    =  dstep * tseas
      baseJD   =  jdate(baseYear  ,baseMon  ,baseDay)
      startJD  =  jdate(startYear ,startMon ,startDay) + 1 + &  
           ( dble((startHour)*3600 + startMin*60 + startSec) / 86400 ) -baseJD
      IF ((IARGC() > 1) )  THEN
         ARG_INT1 = 0.1
         CALL getarg(2,inparg)
         if ( ARG_INT1 == 0) then
            read( inparg, '(i15)' ) ARG_INT1
         else
            read( inparg, '(f15.10)' ) ARG_INT1
         end if
      END IF
    
      IF ((IARGC() > 2) ) THEN
          ARG_INT2 = 0.1
         CALL getarg(3,inparg)
         if ( ARG_INT2 == 0) then
            read( inparg, '(i15)' ) ARG_INT2
         else
            read( inparg, '(f15.10)' ) ARG_INT2
         end if
      END IF

      startYearCond: IF (startYear /= 0) THEN
         IF (ngcm >= 24) THEN 
            intmin      = (startJD)/(ngcm/24)+1
         ELSE ! this needs to be verified
            intmin      = (24*startJD)/ngcm+3-ngcm
         END IF
      END IF startYearCond

      ! tseas - the time step between data sets in [s]
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

!!---------------------------------------------------------------------------
!!------------------------ A L L O C A T I O N S ----------------------------
!!---------------------------------------------------------------------------

      ! --- Allocate information about the coordinates and grid ---

      ALLOCATE ( csu (0:jmt), cst(jmt)  ) 
      ALLOCATE ( phi(0:jmt),   zw(0:km) ) 
      ALLOCATE ( dyt(jmt), dxv(imt+2,jmt), dyu(imt+2,jmt) ) 
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
      hs    = 0.
      uflux = 0.
      vflux = 0.
#ifdef full_wflux
      ALLOCATE ( wflux(imt+2 ,jmt+2 ,0:km,NST) )
#else
      ALLOCATE ( wflux(0:km,NST) )
#endif
      ALLOCATE ( uvel(imt+2,jmt,km) ,vvel(imt+2,jmt,km) ,wvel(imt+2,jmt,km) )
      ALLOCATE ( trj(ntracmax,NTRJ), nrj(ntracmax,NNRJ) )
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
      ALLOCATE ( stxyy(imt,jmt,lbt), stxyx(imt,jmt,lbt) )
      stxyy=0.
      stxyx=0.
#endif
#ifdef streamv
      ALLOCATE ( stxz(imt,km,lbt), styz(jmt,km,lbt) )
      stxz=0.
      styz=0.
#endif
#ifdef streamr
      ALLOCATE ( stxr(imt,mr,lbt,lov), styr(jmt,mr,lbt,lov) )
      stxr=0.
      styr=0
#endif
#ifdef stream_thermohaline
      ALLOCATE ( psi_ts(MR,MR,2,LBT) )
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
#endif

END SUBROUTINE init_params

!!----------------------------------------------------------------------------
!!----------------------------------------------------------------------------
