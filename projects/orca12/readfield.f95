SUBROUTINE readfields
  
   USE mod_precdef
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
   
#ifdef tempsalt
   USE mod_dens
   USE mod_stat
#endif
   IMPLICIT none
   ! ==========================================================================
   ! === Read velocity, temperature and salinity for ORCA0083 configuration ===
   ! ==========================================================================
   ! Subroutine to read the ocean state from ORCA0083 config
   ! Run each time step
   ! --------------------------------------------------------------------------
   ! The following arrays will be populated:
   !
   ! uflux    - Zonal volume flux (U point)
   ! vflux    - Meridional volume flux (V point)
   !
   ! If run with tempsalt option, the following are also set
   ! tem      - Temperature (T point) 
   ! sal      - Salinity (T point)
   ! rho      - Potential density (T point)
   !
   ! --------------------------------------------------------------------------
   
   ! = Loop variables
   INTEGER                                       :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
   INTEGER                                       :: kbot,ktop
   INTEGER, SAVE                                 :: ntempus=0,ntempusb=0,nread
   ! = Variables used for getfield procedures
   CHARACTER (len=200)                           :: fieldFile, medfieldFile
   ! = Variables for filename generation
   CHARACTER (len=200)                           :: dataprefix
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,zstou,zstov,abyst,abysu,abysv
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: xxx
   REAL*4                                      :: dd,hu,hv,uint,vint,zint,hh,h0
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:),SAVE    :: u_m, v_m
  
#ifdef initxyt
   INTEGER, PARAMETER                            :: NTID=73
   INTEGER, PARAMETER                            :: IJKMAX2=7392 ! for distmax=0.25 and 32 days

   INTEGER,  SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
   REAL*4, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
#endif
   
#ifdef tempsalt
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: rhozvec, depthzvec, latvec
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: tmpzvec, salzvec
#endif

   LOGICAL                                       :: around
 
!---------------------------------------------------------------

#ifdef nomean
if (ints == intstart) then
   !! Read mean  
   if(.not. allocated (u_m)) then
   allocate ( u_m(imt,jmt,km), v_m(imt,jmt,km) )
   print*,' Read mean fields '
   u_m(1:imt,1:jmt,1:km) = get3DfieldNC('/group_workspaces/jasmin2/aopp/joakim/ORCA0083-N001/mean/ORCA0083-N01_1978-2010m00U.nc',&
                            & 'vozocrtx')
   v_m(1:imt,1:jmt,1:km) = get3DfieldNC('/group_workspaces/jasmin2/aopp/joakim/ORCA0083-N001/mean/ORCA0083-N01_1978-2010m00V.nc',&
                            & 'vomecrty')
   end if
end if
#endif

#ifdef initxyt
   ! 
   ! Allocate variables necessary for drifter simulation
   !
   alloCondGrid: if ( .not. allocated (ntimask) ) then
      allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
   endif alloCondGrid
#endif
   
   !
   ! Allocate variables 
   !
   alloCondUVW: if(.not. allocated (zstot)) then
      allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
      allocate ( xxx(imt,jmt,km))
#ifdef tempsalt
      allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
#endif
   endif alloCondUVW
 
   call datasetswap ! Swap between current and previous step
   vort(:,:,:,1) = vort(:,:,:,2)
   hdiv(:,:,:,1) = hdiv(:,:,:,2)
   lapu(:,:,:,1) = lapu(:,:,:,2)
   lapv(:,:,:,1) = lapv(:,:,:,2)
   call updateClock 
 
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
               !   print *,j,trajinit(k,i,:)
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
   ! Find the files for the current step
   ! 
   dataprefix='xxxx/ORCA0083-N06_xxxxxxxx'
   write(dataprefix(1:4),'(i4)')   currYear
   write(dataprefix(19:26),'(i4,i2.2,i2.2)') currYear,currMon,currDay
   fieldFile = trim(inDataDir)//'means/'//trim(dataprefix)//'d05'
   medfieldFile = trim(inDataDir)//'medusa/'//trim(dataprefix)//'d05'
   !fieldFile = trim(inDataDir)//'means/2000/ORCA0083-N01_20000105d05'
   
   ! Read SSH
   hs(:,     :, nsp) = get2DfieldNC(trim(fieldFile)//'T.nc', 'ssh')
   hs(imt+1, :, nsp) = hs(1,:,nsp)
   
   ! Depth at U, V, T points as 2D arrays
   allocate ( abyst(imt, jmt) , abysu(imt, jmt) , abysv(imt, jmt) )
   
   abyst = sum(dzt0(:,:,:), dim=3)
   abysu = sum(dzu(:,:,:,1), dim=3)
   abysv = sum(dzv(:,:,:,1), dim=3)
   
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
 
   ! Read temperature 
#if defined tempsalt 
   xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'potemp')
   tem(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
   ! Read salinity
   xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'salin')
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
#endif     
   
   ! Read u, v
   uvel = get3DfieldNC(trim(fieldFile)//'U.nc', 'uo')
   vvel = get3DfieldNC(trim(fieldFile)//'V.nc', 'vo')
   
   ! Put oxygen in salinity field
   xxx(:,:,:) = get3DfieldNC(trim(medfieldFile)//'D.nc', 'TPP3')
   sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
   xxx(:,:,:) = get3DfieldNC(trim(medfieldFile)//'P.nc', 'DIN')
   rho(:,:,:,nsp) = xxx(:,:,km:1:-1)
   !
   ! Calculate zonal and meridional volume flux
   !
   ! Weight by (1 + ssh / depth)
   ! This is only an approximation of what NEMO really does
   ! but is accurate within 1% 
   !
   
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      dzt(i,j,k,nsp) = dzt0(i,j,k) * zstot(i,j)
   end do
   end do
   end do

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
      
   uvel(:,:,:) = 0.
   vvel(:,:,:) = 0.
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      uvel(i,j,km+1-k) = uflux(i,j,k,nsp) 
      vvel(i,j,km+1-k) = vflux(i,j,k,nsp) 
   enddo
   enddo
   enddo
   
   !! calculate volume fluxes
   uflux(:,:,:,nsp) = 0.
   vflux(:,:,:,nsp) = 0.
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      uflux(i,j,km+1-k,nsp) = uvel(i,j,km+1-k) * dyu(i,j) * dzu(i,j,km+1-k,1) * zstou(i,j)
      vflux(i,j,km+1-k,nsp) = vvel(i,j,km+1-k) * dxv(i,j) * dzv(i,j,km+1-k,1) * zstov(i,j)
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



