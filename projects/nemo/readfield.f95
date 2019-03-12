SUBROUTINE readfields
   
   ! ==========================================================================
   ! 
   ! Purpose
   ! -------
   !
   ! Read NEMO model output to advect trajectories. 
   ! Will be called by loop_pos each time step. 
   !
   ! Method
   ! ------
   !
   ! Read velocities and optionally some tracers from netCDF files and 
   ! update velocity fields for TRACMASS. 
   ! Optionally read temperature, salinity and other tracers, e.g. oxygen, nitrate
   ! 
   ! Updates the variables: 
   !   uflux and vflux (with sub-grid vels if sgsUV)
   !   dzt (if vvl)
   !   temp, salt, dens (if readTS)
   !
   ! There are a bunch of options controlling what model output you have. 
   ! The basic idea is that you set a physPrefixForm e.g. "ORCA1_YYYYMMDD" 
   ! and then this routine will replace YYYY by year, MM with month etc
   ! and then read the data. 
   ! If you have more than one step per file (oneStepPerFile = .false.)
   ! then you can set a "timestamp" e.g. for some runs you have monthly data
   ! with one file per year, so the files are named "ORCA1_19990101_19991231"
   ! In this case, the "timestamp" is "YYYY0101_YYYY1231" and this routine
   ! will replace YYYY by 1999. 
   ! For other runs, e.g. INALT60, there are several step per file, but not 
   ! always the same number. There is a special case for this. 
   !
   ! Some nomenclature:
   ! physDataDir and physPrefixForm refers to physical variables, e.g. u,v,T,S,MLD
   ! bioDataDir and bioPrefixForm refers to biogeochemical variables, N, DO, DIC etc. 
   !
   ! 
   ! Future
   ! ------
   !
   ! We should have a namelist variable that is a 2D array of character strings e.g. 
   ! physTracers(1,1:3) = 'votemper', '3D', 'theta' 
   ! physTracers(2,1:3) = 'vosaline', '3D', 'salinity'
   ! physTracers(3,1:3) = 'somxl010', '2D', 'mld'
   ! bioTracers(1,1:3)  = 'DIC', '3D', 'dic'
   ! and from this, TRACMASS knows what is temperature and salinity and that it can
   ! calculate density, and also knows that MLD is only 2D, and that DIC is a bio tracer
   ! that can be used to make total carbon etc. 
   !
   ! History
   ! -------
   ! 
   ! J. Kjellsson - 02.2019 Unified projects for various NEMO configs
   !
   !
   ! ==========================================================================
   
   USE mod_precdef
   USE mod_calendar
   USE mod_param
   USE mod_vel
   
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_traj
   USE mod_getfile
   use mod_seed
   use mod_tempsalt
   USE mod_deformation
   USE mod_laplacian
   
   USE mod_dens
   USE mod_stat
   
   IMPLICIT none
   
   ! -----------------------------------------------------------------------------
   
   INTEGER                                       :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
   INTEGER                                       :: kbot,ktop, idiag, jdiag
   INTEGER                                       :: ichar
   INTEGER, SAVE                                 :: ntempus=0,ntempusb=0,nread,itime, fieldStep 
   
   ! Variables to set the filenames
   CHARACTER (len=200)                           :: fieldFile, medfieldFile, umFile, vmFile, physPrefix
   CHARACTER (len=200)                           :: tFile, uFile, vFile, wFile, bgcFile
   CHARACTER (len=100)                           :: tmpstr
   CHARACTER (len=200)                           :: dataprefix, timestamp
   INTEGER, ALLOCATABLE, DIMENSION(:,:),SAVE     :: fileMon, fileDay
   CHARACTER (len=200), ALLOCATABLE, DIMENSION(:),SAVE :: file_timestamp
   
   ! Variables to calculate thickness of layers with vvl
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,zstou,zstov,abyst,abysu,abysv
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: xxx
   REAL*4                                        :: dd,hu,hv,uint,vint,zint,hh,h0
   
   ! Time mean u,v if read in
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:),SAVE    :: u_m, v_m
   
   ! Variables to calculate potential density
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: rhozvec, depthzvec, latvec
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: tmpzvec, salzvec
   
   ! Control what to read (should be namelist variables)
   LOGICAL                                       :: around !, useTrmClock = .false.
   
   INTEGER, PARAMETER                            :: NTID=73
   INTEGER, PARAMETER                            :: IJKMAX2=7392 ! for distmax=0.25 and 32 days

   INTEGER,  SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
   REAL*4, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
   
   
   ! -----------------------------------------------------------------------------
   
   if (ints == intstart) then
      
      if (log_level >= 1) then
      print*,' Run ID:                       ',trim(RunID)
      print*,' Prefix for phys files:        ',trim(physPrefixForm)
      print*,' Dir for phys files:           ',trim(physDataDir)
      print*,' Read mean fields:             ',readMean
      print*,' Variable volume layer (vvl):  ',vvl
      print*,' Read SSH:                     ',readSSH
      print*,' Read BGC tracers:             ',readBio
      print*,' Read SGS velocities:          ',sgsUV
      print*,' One step per file:            ',oneStepPerFile
      print*,' Name for gridT:               ',trim(tGridName)
      print*,' File suffix:                  ',trim(fileSuffix)
      
      print*,' Start year,mon,day, hour      ',startYear,startMon,startDay,startHour
      end if
      
      if (.not. useTrmClock) then
         call init_calendar
         !currHour = startHour
         !currDay = startDay
         !currMon = startMon
         !currYear = startYear
      end if
      
      if (readMean) then
         !! Read mean  
         if(.not. allocated (u_m)) then
            allocate ( u_m(imt,jmt,km), v_m(imt,jmt,km) )
            print*,' Read mean fields '
            
            umFile = '/group_workspaces/jasmin2/aopp/joakim/ORCA0083-N001/mean/ORCA0083-N01_1978-2010m00U.nc'
            vmFile = '/group_workspaces/jasmin2/aopp/joakim/ORCA0083-N001/mean/ORCA0083-N01_1978-2010m00V.nc'
            u_m(1:imt,1:jmt,1:km) = get3DfieldNC(umFile,'vozocrtx')
            v_m(1:imt,1:jmt,1:km) = get3DfieldNC(vmFile,'vomecrty')
         end if
      end if
      
      if (.not. oneStepPerFile) then
         allocate( file_timestamp(fieldsperfile), fileDay(fieldsperfile,2), fileMon(fieldsperfile,2) )
         fileDay(:,1) = 1
         fileDay(:,2) = 31
         fileMon(:,1) = 1
         fileMon(:,2) = 12
      end if
      
      if (RunID == '2_INALT60.L120-KRS0020_4h') then
         !
         ! INALT60 data with 4h frequency
         ! Here we set start/stop dates for each file in each year
         !
         fieldsperfile = 16
         deallocate( file_timestamp, fileDay, fileMon )
         allocate( file_timestamp(16), fileDay(16,2), fileMon(16,2) )
         fileDay(1:16,1) = (/  1,  6, 31, 25, 22, 16, 11,  5, 30, 25, 19, 13,  8,  2, 27, 22 /)
         fileDay(1:16,2) = (/  5, 30, 24, 21, 15, 10,  4, 29, 24, 18, 12,  7,  1, 26, 21, 31 /)
         fileMon(1:16,1) = (/  1,  1,  1,  2,  3,  4,  5,  6,  6,  7,  8,  9, 10, 11, 11, 12 /)
         fileMon(1:16,2) = (/  1,  1,  2,  3,  4,  5,  6,  6,  7,  8,  9, 10, 11, 11, 12, 12 /)
         
      end if
      
      if (.not. oneStepPerFile) then   
         ! Now write year to timestamp string
         tmpstr = "YYYYMMDD_YYYYMMDD"
         do ii=1,fieldsperfile
            write(tmpstr(5:8)  ,'(i2.2,i2.2)') fileMon(ii,1),fileDay(ii,1)
            write(tmpstr(14:17),'(i2.2,i2.2)') fileMon(ii,2),fileDay(ii,2)
            file_timestamp(ii) = trim(tmpstr)
         end do      
         
         !  
         ! Find the files for the current step  
         ! 
         fieldStep = 1
         itime = 1
         do while (currMon*100 + currDay - (fileMon(itime,2)*100 + fileDay(itime,2)) > 0)
            itime = itime + 1
         end do  
      end if
      
      if (.not. oneStepPerFile) then
         if (log_level > 0) then
            print*,' First file timestamp:         ',file_timestamp(1)
            print*,' Start with file number:       ',itime
         end if
      end if
      
   end if
   
#ifdef initxyt
   ! 
   ! Allocate variables necessary for drifter simulation
   !
   alloCondGrid: if ( .not. allocated (ntimask) ) then
      allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
   endif alloCondGrid
#endif
   
   !
   ! Allocate temporary array for u,v etc.
   !
   if (.not. allocated(xxx)) then
      allocate ( xxx(imt,jmt,km))
   end if
   
   !
   ! Allocate variables for vvl calculations 
   !
   if(vvl .and. .not. allocated (zstot)) then
      allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
   end if   
   
   !
   ! Allocate variables for density calculations
   !
   if (readTS .and. .not. allocated(tmpzvec)) then
      allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
   end if
   
   !
   ! Swap time levels and update clock
   ! 
   call datasetswap ! Swap between current and previous step
   vort(:,:,:,1) = vort(:,:,:,2)
   hdiv(:,:,:,1) = hdiv(:,:,:,2)
   lapu(:,:,:,1) = lapu(:,:,:,2)
   lapv(:,:,:,1) = lapv(:,:,:,2)
   
   if (useTrmClock) then
      call updateClock 
   else
      if (ints /= intstart) then 
         call update_calendar
      end if
      !if (ints /= intstart) then
      !daysInMonth = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
      !if (ngcm_unit == "minute") then
      !   currStep = ngcm 
      !else if (ngcm_unit == "hour") then
      !   currStep = ngcm * 60
      !else if (ngcm_unit == "day") then
      !   currStep = ngcm * 24 * 60
      !else if (ngcm_unit == "month") then
      !   currStep = daysInMonth(currMon) * 24 * 60 
      !end if
      !! Update the clock manually
      !currHour = currHour + currStep/60
      !do while (currHour >= 24) 
      !   currDay = currDay + 1
      !   currHour = currHour - 24
      !   if (currDay > daysInMonth(currMon)) then
      !      currDay = currDay - daysInMonth(currMon)
      !      currMon = currMon + 1
      !   end if
      !   if (currMon > 12) then
      !      currMon = currMon - 12
      !      currYear = currYear + 1
      !   end if
      !end do
      !end if
   end if
   
   ! === Initialising fields ===
   initFieldcond: if(ints == intstart) then
#ifdef initxyt
      ! Time for individual start positions
      if(IJKMAX2 == 7392) open(84,file=trim(inDataDir)//'topo/masktime_32_025', &
                               form='unformatted')
      read(84) trajinit
      close(84)
      j=0
      do k=1,NTID
         do i=1,IJKMAX2
            if(trajinit(k,i,3) /= 0.) then
               j=j+1
#if orca025l75h6
               trajinit(k,i,3)=float(kst2)-0.5
#endif
            endif
         enddo
      enddo
      ijkst=0
      if(j /= IJKMAX2) then
         stop 4396
      endif
#endif
   
   endif initFieldcond
   
   !
   ! Determine which time level in which file to read
   ! i.e. set itime and fieldStep
   ! If only one step per file, this is always 1. 
   !
   
   ncTpos = fieldStep
   
   if (oneStepPerFile) then
      if (log_level > 0) then
         print*,' Only one step per file '
      end if
      itime = 1
      fieldStep = 1
      
      timestamp = ''
      
   else if (RunID == '2_INALT60.L120-KRS0020_4h') then
      ! 
      ! Special case for INALT60
      !
      fieldStep = fieldStep + 1   
      if (log_level >= 2) then
         print*,' Reading file time stamp, fieldStep ',itime,fieldStep
      end if
      if ( (fieldStep > 150 .and. ( fileMon(itime,1)*100+fileDay(itime,1) > 101 .and. &
                                  & fileMon(itime,1)*100+fileDay(itime,1) < 1222 ) ) .or. & 
         & (fieldStep > 30  .and. ( fileMon(itime,1)*100+fileDay(itime,1) == 101)  ) .or. &
         & (fieldStep > 60  .and. ( fileMon(itime,1)*100+fileDay(itime,1) == 1222) ) ) then
         fieldStep = 1
         itime = itime + 1
         if (itime > 16) then
            itime = 1
         end if
      end if 
      
      timestamp = trim(file_timestamp(itime))
   
   else
      fieldStep = fieldStep + 1
      if (fieldStep > fieldsperfile) then
         fieldStep = 1
      end if
      
   end if   
   
   !
   ! Set the file names that we need to read
   ! We use a physPrefixForm and fill in RunID, currYear etc. 
   !
   
   physPrefix = physPrefixForm
   ichar = INDEX(physPrefix,'RUNID')
   do while (ichar /= 0)
      physPrefix = trim(physPrefix(:ichar-1))//trim(RunID)//trim(physPrefix(ichar+5:))
      ichar = INDEX(physPrefix,'RUNID')
   end do
   
   ichar = INDEX(physPrefix,'YYYY')
   do while (ichar /= 0)
      write(physPrefix(ichar:ichar+3),'(i4)') currYear
      ichar = INDEX(physPrefix,'YYYY')
   end do
   
   ichar = INDEX(timestamp,'YYYY')
   do while (ichar /= 0)
      write(timestamp(ichar:ichar+3),'(i4)') currYear
      ichar = INDEX(timestamp,'YYYY')
   end do
   
   ichar = INDEX(physPrefix,'MM')
   do while (ichar /= 0)
      write(physPrefix(ichar:ichar+1),'(i2.2)') currMon
      ichar = INDEX(physPrefix,'MM')
   end do
      
   ichar = INDEX(physPrefix,'DD')
   do while (ichar /= 0)
      write(physPrefix(ichar:ichar+1),'(i2.2)') currDay   
      ichar = INDEX(physPrefix,'DD')
   end do 
   
   ichar = INDEX(physPrefix,'TSTSTSTS')
   do while (ichar /= 0)
      physPrefix = trim(physPrefix(:ichar-1))//trim(timestamp)//trim(physPrefix(ichar+8:))
      ichar = INDEX(physPrefix,'TSTSTSTS')
   end do
   
   ichar = INDEX(physPrefix,'GRIDX')
   tFile = trim(physDataDir)//trim(physPrefix(:ichar-1))//trim(tGridName)//trim(physPrefix(ichar+5:))//trim(fileSuffix)
   uFile = trim(physDataDir)//trim(physPrefix(:ichar-1))//trim(uGridName)//trim(physPrefix(ichar+5:))//trim(fileSuffix)
   vFile = trim(physDataDir)//trim(physPrefix(:ichar-1))//trim(vGridName)//trim(physPrefix(ichar+5:))//trim(fileSuffix)
   
   if (ints == intstart) then
      if (log_level > 0) then
         print*,' First T file:                 ',trim(tFile)
         print*,' First U file:                 ',trim(uFile)
         print*,' First V file:                 ',trim(vFile)
      end if
   else 
      if (log_level > 1) then
         print*,' Current T file:               ',trim(tFile)
         print*,' Current U file:               ',trim(uFile)
         print*,' Current V file:               ',trim(vFile)
      end if
   end if
    
   ! Read SSH
   if (readSSH .or. vvl) then
      hs(:,     :, nsp) = get2DfieldNC(trim(tFile), ssh_name)
      hs(imt+1, :, nsp) = hs(1,:,nsp)
   end if
   
   ! Depth at U, V, T points as 2D arrays
   allocate ( abyst(imt, jmt) , abysu(imt, jmt) , abysv(imt, jmt) )
   
   abyst = sum(dzt0(:,:,:), dim=3)
   abysu = sum(dzu(:,:,:,1), dim=3)
   abysv = sum(dzv(:,:,:,1), dim=3)
   
   if (vvl) then
      ! Calculate SSH/depth
      where (abyst /= 0)
         zstot = hs(:imt,:jmt,nsp)/abyst + 1
      elsewhere
         zstot = 0.d0
      end where
      
      where (abysu /= 0)
         zstou = 0.5*(hs(:imt,:jmt,nsp)+hs(2:imt+1,:jmt,nsp))/abysu + 1
      elsewhere
         zstou = 0.d0
      end where
   
      where (abysv /= 0)
         zstov = 0.5*(hs(:imt,:jmt,nsp)+hs(:imt,2:jmt+1,nsp))/abysv + 1
      elsewhere
         zstov = 0.d0
      end where
   end if
   
   ! Read temperature 
   if (readTS) then
      xxx(:,:,:) = get3DfieldNC(trim(tFile), temp_name)
      tem(:,:,:,nsp) = xxx(:,:,km:1:-1)
      
      ! Read salinity
      xxx(:,:,:) = get3DfieldNC(trim(tFile), salt_name)
      sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
      ! Calculate potential density
      depthzvec = 0.
      do j=1,jmt
         latvec=-80+1./12.*float(j+subGridJmin-1)
         do i=1,IMT
            tmpzvec = tem(i,j,:,nsp)
            salzvec = sal(i,j,:,nsp)
            call statvd(tmpzvec, salzvec, rhozvec ,km ,depthzvec ,latvec)
            rho(i,j,:,nsp)=rhozvec - 1000.
         end do
      end do
   end if     
   
   !
   ! Read u, v
   !
   xxx(:,:,:) = get3DfieldNC(trim(uFile), trim(ueul_name)) 
   uvel(1:imt,1:jmt,1:km) = xxx(:,:,1:km)
   xxx(:,:,:) = get3DfieldNC(trim(vFile), trim(veul_name))
   vvel(1:imt,1:jmt,1:km) = xxx(:,:,1:km)
   
   if (sgsUV) then
      !
      ! Read subgrid scale velocities (e.g. EIV from NEMO)
      !
      xxx(:,:,:) = get3DfieldNC(trim(uFile), trim(usgs_name))
      uvel(1:imt,1:jmt,1:km) = uvel(1:imt,1:jmt,1:km) + xxx(:,:,1:km)
      xxx(:,:,:) = get3DfieldNC(trim(vFile), trim(vsgs_name))   
      vvel(1:imt,1:jmt,1:km) = vvel(1:imt,1:jmt,1:km) + xxx(:,:,1:km)
   end if
   
   if (readBio) then
      ! Put oxygen in salinity field
      xxx(:,:,:) = get3DfieldNC(trim(medfieldFile)//'D.nc', 'TPP3')
      sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
      
      xxx(:,:,:) = get3DfieldNC(trim(medfieldFile)//'P.nc', 'DIN')
      rho(:,:,:,nsp) = xxx(:,:,km:1:-1)
   end if 
   
   
   !
   ! Calculate zonal and meridional volume flux
   !
   ! Weight by (1 + ssh / depth)
   ! This is only an approximation of what NEMO really does
   ! but is accurate within 1% 
   !
   
   if (vvl) then
      do k = 1, km
      do j = 1, jmt
      do i = 1, imt
         dzt(i,j,k,nsp) = dzt0(i,j,k) * zstot(i,j)
      end do
      end do
      end do
   end if
   
#ifdef nomean   
   !! flip u,v upside down
   !! use uflux, vflux as temporary arrays
   where (uvel == 0)
      u_m = 0
   end where
   where (vvel == 0)
      v_m = 0
   end where
   
   uflux(1:imt,1:jmt,1:km,nsp) = uvel(1:imt,1:jmt,1:km) - u_m(1:imt,1:jmt,1:km)
   vflux(1:imt,1:jmt,1:km,nsp) = vvel(1:imt,1:jmt,1:km) - v_m(1:imt,1:jmt,1:km)

#else
   uflux(1:imt,1:jmt,1:km,nsp) = uvel(1:imt,1:jmt,1:km) 
   vflux(1:imt,1:jmt,1:km,nsp) = vvel(1:imt,1:jmt,1:km) 
#endif
      
   !! calculate volume fluxes
   uflux(:,:,:,nsp) = 0.
   vflux(:,:,:,nsp) = 0.
   !idiag = 844
   !jdiag = 1410
   !print*,'uvel ',uvel(idiag,jdiag,:)
   !print*,'dyu ',dyu(idiag,jdiag)
   !print*,'dzu ',dzu(idiag,jdiag,:,1)
   !print*,'zstou ',zstou(idiag,jdiag)
   !print*,'kmu ',kmu(idiag,jdiag)
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt      
      uflux(i,j,km+1-k,nsp) = uvel(i,j,k) * dyu(i,j) * dzu(i,j,km+1-k,1) * zstou(i,j)  
      vflux(i,j,km+1-k,nsp) = vvel(i,j,k) * dxv(i,j) * dzv(i,j,km+1-k,1) * zstov(i,j)
   enddo
   enddo
   enddo
   
   !! calculate laplacian of u,v
   call laplacian
      
   ! Check that volume fluxes are zero below sea floor
   do i=1,IMT
   do j=1,JMT
   do k=1,KM
   if(k > kmv(i,j) .and. vflux(i,j,km+1-k,nsp) /= 0.) then
      print *,'vflux=',vflux(i,j,km+1-k,nsp),vvel(i,j,k),i,j,k,kmv(i,j),nsp
      stop 4966
   endif
   if(k > kmu(i,j) .and. uflux(i,j,km+1-k,nsp) /= 0.) then
      print *,'uflux=',uflux(i,j,km+1-k,nsp),uvel(i,j,k),i,j,k,kmu(i,j),nsp
      stop 4967
   endif
   enddo
   enddo
   enddo
   
   
#ifdef drifter
   ! average velocity/transport to simulate drifter trajectories
   kbot=65 ; ktop=66 ! number of surface layers to integrate over 
   uint=0. ; vint=0. ; zint=0.
   do k=kbot,ktop
      uint = uint + uflux(:,:,k,nsp) ! integrated transport
      vint = vint + vflux(:,:,k,nsp)
      zint = zint + dz(k)          ! total depth of drougued drifter
   end do
   ! weighted transport for each layer
   do k=kbot,KM
      uflux(:,:,k,nsp) = uint*dz(k)/zint 
      vflux(:,:,k,nsp) = vint*dz(k)/zint
   enddo
#endif

#ifdef initxyt
   ! Set the initial trajectory positions
   !ijkst(:,5)=ntimask(ntempus,:)
#ifdef orca025l75h6
   if( mod(ints,24/ngcm*5) == 1 .or. ints <= 2) ntempus=ntempus+1
   if(ntempus /= ntempusb .and. ntempus <= NTID) then
      ntempusb=ntempus
      !print *,'ints=',ints,' ntempus=',ntempus,' ntempusb=',ntempusb
#else
   if(ints.le.NTID) then
#endif
   
      do ntrac=1,ijkmax
      if(trajinit(ntempus,ntrac,3) /= 0.) then
         ijkst(ntrac,4)=0
         ijkst(ntrac,5)=5
         ijkst(ntrac,6)=ijkst(ntrac,6)+1
         do l=1,3
            ijkst(ntrac,l)=trajinit(ntempus,ntrac,l)+1
            trj(ntrac,l)=trajinit(ntempus,ntrac,l)
            !if(l.eq.1) print *,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
            if(trj(ntrac,l).gt.float(ijkst(ntrac,l)) .or. trj(ntrac,l).lt.float(ijkst(ntrac,l)-1)) then
               print *,l,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
               stop 3946
            endif
         enddo
      else
         ijkst(ntrac,5)=0
         ijkst(ntrac,6)=0
      endif
      enddo
   endif

#ifdef orca025l75h6
#endif
   if( mod(ints,24/ngcm*5).ne.1 .and. ints.gt.2) then
      ijkst=0 
   endif
#endif

   return
   
end subroutine readfields



