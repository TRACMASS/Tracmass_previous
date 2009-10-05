subroutine init_params
  USE mod_param
  USE mod_seed
  USE mod_coord
  USE mod_grid
  USE mod_name
  USE mod_time 
  USE mod_domain
  USE mod_vel
  USE mod_dens
  USE mod_buoyancy
  USE mod_streamxy
  USE mod_streamv
  USE mod_streamr
  USE mod_tracer
#ifdef diffusion
  USE mod_diffusion
#endif
#ifdef sediment
  USE mod_orbital
  USE mod_sed
#endif
  implicit none

  INTEGER                                    :: argint1 ,argint2
  INTEGER                                    :: dummy ,factor ,i ,dtstep
  INTEGER                                    :: gridVerNum ,runVerNum
  CHARACTER (LEN=30)                         :: inparg, argname
  CHARACTER (LEN=23)                         ::  Project, Case

  namelist /INITGRIDVER/   gridVerNum
  namelist /INITGRIDGRID/  IMT, JMT, KM, LBT, NEND
  namelist /INITGRIDNTRAC/ NTRACMAX
  namelist /INITGRIDDATE/  yearmin, yearmax, baseSec  ,baseMin  ,baseHour &
                          ,baseDay  ,baseMon  ,baseYear
  namelist /INITGRIDTIME/  ngcm, iter, intmax ,fieldsPerFile
  namelist /INITGRIDARC/   arcscale

  namelist /INITRUNVER/    runVerNum
  namelist /INITRUNGRID/   subGrid ,subGridImin ,subGridImax ,subGridJmin &
                          ,subGridJmax ,SubGridFile, subGridID
  namelist /INITRUNTIME/   intmin, intspin, intrun, intstep 
  namelist /INITRUNDATE/   startSec ,startMin ,startHour &
                          ,startDay ,startMon ,startYear &
                          ,ihour, iday, imon, iyear
  namelist /INITRUNWRITE/  ncoor, kriva ,inDataDir ,outDataDir &
                          ,outDataFile ,intminInOutFile
  namelist /INITRUNSEED/   nff ,isec ,idir ,nqua ,partQuant ,seedType &
                          ,ist1 ,ist2 ,jst1 ,jst2 ,kst1 ,kst2 &
                          ,varSeedFile ,seedDir ,seedFile
  namelist /INITGRIDDESC/  GCMname,GCMsource,gridName,gridSource &
                          ,gridDesc  
  namelist /INITRUNDESC/  caseName ,caseDesc  
#ifdef tempsalt
  namelist /INITRUNTEMPSALT/ tmin0, tmax0, smin0, smax0, rmin0, rmax0, &
                             tmine, tmaxe, smine, smaxe, rmine, rmaxe
#endif
#ifdef diffusion
  namelist /INITRUNDIFFUSION/ ah, av
#endif
#ifdef sediment
  namelist /INITRUNSEDIMENT/ partdiam, rhos, cwamp, twave, critvel
#endif


!  allocate ( ienw (LBT),iene (LBT) )
!  allocate ( jens (LBT),jenn (LBT) )

  namelist /INITRUNEND/ ienw, iene, jens, jenn, timax


  Project  = PROJECT_NAME
  Case     = CASE_NAME
  
  if ( (IARGC() .eq. 1 ) .or. (IARGC() .eq. 4 ) )  then
     call getarg(IARGC(),inparg)
     Case=inparg
  end if
  
  ! -- Check if there is a time argument and if so, use it.
!  print *,trim(Project)//'/'//trim(Case)//'_grid.in'
  open(8,file='projects/'//trim(Project)//'/'//trim(Project)//'_grid.in',  &
       status='OLD', delim='APOSTROPHE')
  ! -- Check if the namefiles has correct version number. 
  read(8,nml=INITGRIDVER)
  if (gridVerNum < 1) then
     print *,'==================== ERROR ===================='
     print *,'Your grid namefile seems to be out of date.'
     print *,'Check the version_grig.txt file for changes.'
     print *,'You have to edit the version number in you grid'
     print *,'manually when done.'
     stop
  end if
  read(8,nml=INITGRIDDESC)
  read(8,nml=INITGRIDGRID)
  read(8,nml=INITGRIDNTRAC)
  read(8,nml=INITGRIDTIME)
  read(8,nml=INITGRIDDATE)
  !=== === ===
  read(8,nml=INITGRIDARC)
  open(8,file='projects/'//trim(Project)//'/'//trim(Case)//'_run.in',  &
       status='OLD', delim='APOSTROPHE')
  read(8,nml=INITRUNDESC)
  read(8,nml=INITRUNGRID)
  select case (subGrid)
  case (0)
     print *,'Use the Full grid.'     
     subGridImin =   1
     subGridJmin =   1
     subGridImax = imt
     subGridJmax = jmt
  case (1)
     print *,'Use a subgrid: ', subGridImin ,subGridImax &
          ,subGridJmin ,subGridJmax
     imt=subGridImax-subGridImin+1
     jmt=subGridJmax-subGridJmin+1
  case default
     print *,'==================== ERROR ===================='
     print *,'This subGrid selection is not implemented yet.'
     print *,'subGrid = ' ,subGrid
     stop
  end select
  !=== === ===
  read(8,nml=INITRUNTIME)
  read(8,nml=INITRUNDATE)
  read(8,nml=INITRUNWRITE)  
  read(8,nml=INITRUNSEED)
#ifdef tempsalt
  read(8,nml=INITRUNTEMPSALT)
#endif
#ifdef diffusion
  read(8,nml=INITRUNDIFFUSION)
#endif
#ifdef sediment
  read(8,nml=INITRUNSEDIMENT)
#endif

  read(8,nml=INITRUNEND)

  timax=24.*3600.*timax ! convert time lengths from days to seconds
  dstep=1.d0/dble(iter)
  dtmin=dtstep*tseas

  if ((IARGC() > 1) .and. (IARGC() < 5) )  then
     call getarg(2,inparg)
     factor=1
     argint1=0
     do i=29,1,-1
        if (ichar(inparg(i:i)) .ne. 32) then
           argint1=argint1+(ichar(inparg(i:i))-48)*factor
           factor=factor*10
        end if
     end do
     ARG_INT1=argint1

!     call getarg(2,inparg)
!     factor=1
!     argint2=0
!     do i=29,1,-1
!        if (ichar(inparg(i:i)) .ne. 32) then
!           argint2=argint2+(ichar(inparg(i:i))-48)*factor
!           factor=factor*10
!        end if
!     end do
!     ARG_INT2=argint2
  
  end if

  baseJD      = jdate(baseYear  ,baseMon  ,baseDay)
  startJD     = jdate(startYear ,startMon ,startDay)
  startYearCond: if (startYear .ne. 0) then
   if(ngcm.ge.24) then 
    intmin      = (startJD-baseJD)/(ngcm/24)+1
   else ! this is a quick fix to avoid division by zero when ngcm < 24
    intmin      = (startJD-baseJD)/(ngcm)+1
   endif
  end if startYearCond
  startHourCond: if ( (startHour .ne. 0)  & 
                 .or. (startMin  .ne. 0)  &
                 .or. (startSec  .ne. 0) ) then
     print *,'------------------------------------------------------'
     print *,'ERROR!'
     print *,'------------------------------------------------------'
     print *,'Fractions of day for start values not implemented yet.'
     print *,'Email Bror (brorfred@gmail.com) to have it fixed.'
     stop 100
  end if startHourCond

  tseas= dble(ngcm)*3600.d0 ! time step between data sets


  ! --ist -1 to imt
  if ( ist1 == -1) then 
     ist1=imt
  end if
  if ( ist2 == -1) then 
     ist2=imt
  end if
  ! --jst -1 to jmt
  if( jst1 == -1) then 
     jst1=jmt
  end if
  if ( jst2 == -1) then 
     jst2=jmt
  end if
  ! --ist -1 to imt  
  if ( kst1 == -1) then 
     kst1=km
  end if
  if ( kst2 == -1) then 
     kst2=km
  end if

  ! mod_coord
  allocate ( csu (0:jmt), cst(jmt)  ) 
  allocate ( phi(0:jmt),   zw(0:km) ) 
  allocate ( dyt(jmt), dxv(imt+2,jmt), dyu(imt+2,jmt), dzt(imt+2,jmt,km) ) 
  ! mod_grid
#if defined ifs || atm
  allocate ( dztb(imt,jmt,km,nst) )   
#else
  allocate ( dztb(imt,jmt,km) )   
#endif
  allocate ( dxdy(imt,jmt) )   
  allocate (kmt(imt,jmt), dz(km) )
  ! mod_domain
!  allocate ( ienw (LBT),iene (LBT) )
!  allocate ( jens (LBT),jenn (LBT) )
  allocate ( uflux(imt,jmt,km,nst), vflux(imt,0:jmt,km,nst) )
  allocate ( hs(imt+1,jmt+1,nst) )
#ifdef full_wflux
  allocate ( wflux(imt+2 ,jmt+2 ,0:km ,2) )
#elif defined timeanalyt || timestep
  allocate ( wflux(0:km,2) )
#else
  allocate ( wflux(0:km) )
#endif
  allocate ( uvel(imt+2,jmt,km) ,vvel(imt+2,jmt,km) ,wvel(imt+2,jmt,km))

  ! mod_dens
#ifdef tempsalt
  allocate ( tem(imt,jmt,km,nst) ,sal(imt,jmt,km,nst), rho(imt,jmt,km,nst) )
#endif
  ! mod_streamxy
#ifdef streamxy
  allocate ( stxyy(imt,jmt,lbt), stxyx(imt,jmt,lbt) )
#endif
  ! mod_streamv
#ifdef streamv
  allocate ( stxz(imt,km,lbt), styz(jmt,km,lbt) )
#endif
  ! mod_streamr
#ifdef streamr
  allocate ( stxr(imt,mr,lbt,lov), styr(jmt,mr,lbt,lov) )
#endif
! mod_tracer
#ifdef tracer
  allocate ( tra(imt,jmt,km) )
#endif
  ! mod_sed
#ifdef sediment
  allocate (orb(km) )
#endif

end subroutine init_params

