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
   CHARACTER (len=200)                           :: fieldFile
   ! = Variables for filename generation
   CHARACTER (len=200)                           :: dataprefix
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,abyst
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_doub
   REAL*4  , ALLOCATABLE, DIMENSION(:,:)         :: temp2d_simp
   REAL*4  , ALLOCATABLE, DIMENSION(:,:,:)       :: temp3d_simp
   REAL*4                                        :: dd,hu,hv,uint,vint,zint,hh,h0
  
   
#ifdef tempsalt
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: rhozvec, depthzvec, latvec
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: tmpzvec, salzvec
#endif

   LOGICAL                                       :: around
   
   !
   ! Allocate variables 
   !
   alloCondUVW: if(.not. allocated (zstot)) then
      allocate ( temp2d_simp(imt,jmt) , zstot(imt,jmt) )
      allocate ( temp3d_simp(imt,jmt,km) )
#ifdef tempsalt
      allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
#endif
   endif alloCondUVW
 
   call datasetswap ! Swap between current and previous step


! Initialise on the first time step
IF (ints == intstart) THEN
   hs    = 0.
   uflux = 0.
   vflux = 0.
#ifdef tempsalt
   tem   = 0.
   sal   = 0.
   rho   = 0.
#endif

   iyear = startYear
   imon  = startMon
   iday  = startDay
   ihour = startHour

END IF

!! Update the time counter

ihour = ihour + int(ff) * ngcm  

IF (ihour >= 24) THEN   !Forward scheme

   ihour = 0
   iday  = iday + 1
   IF (iday > idmax(imon,iyear)) THEN
      iday = 1
      imon = imon + 1
      ncTpos=0
      IF (imon == 13) THEN
         imon  = 1
         iyear = iyear + 1
      END IF
   END IF

ELSE IF (ihour < 0) THEN  !Backward scheme
   
   ihour = 18
   iday  = iday - 1
   IF (iday == 0) THEN
      imon = imon - 1
      IF (imon == 0) THEN
         imon  = 12
         iyear = iyear - 1
      END IF
      iday=idmax(imon,iyear)
   END IF

END IF



ntime = 1000000 * iyear + 10000 * imon + 100 * iday + ihour

ncTpos = ints-113998

! print *,'ints',ints,ncTpos,ntime
 
! stop 4957
 
  start1D  = [ 1]
  count1D  = [KM]
  start2D  = [1 ,1 , ncTpos , 1 ]
  count2D  = [          imt,        jmt ,  1 , 1 ]
  map2D    = [          1 ,          2 ,  3 , 4 ]  
  start3D  = [1 ,1 ,  1 , ncTpos ]
  count3D  = [         imt,        jmt , KM , 1 ]
  map3D    = [          1 ,          2 ,  3 , 4 ] 

!---------------------------------------------------------------


   !
   ! Find the files for the current step
   ! 
!   dataprefix='xxxx/ORCA0083-N01_xxxxxxxx'
!   write(dataprefix(1:4),'(i4)')   currYear
!   write(dataprefix(19:26),'(i4,i2.2,i2.2)') currYear,currMon,currDay
   fieldFile = '/Users/doos/data/nemonordic/fields/NORDIC-NS2_1h_20130101_20131231_GOF_grid_'
   !fieldFile = trim(inDataDir)//'means/2000/ORCA0083-N01_20000105d05'
   
   ! Read SSH
!   hs(1:imt,:, nsp) = get2DfieldNC(trim(fieldFile)//'T.nc', 'SSH')
!   temp2d_simp = get2DfieldNC(trim(fieldFile)//'T.nc', 'SSH')
   
  ierr=NF90_OPEN(trim(fieldFile)//'T.nc',NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 6751
  ierr=NF90_INQ_VARID(ncid,'SSH',varid)
  if(ierr.ne.0) stop 6763
  ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
  if(ierr.ne.0) stop 6799

   do j = 1, jmt
   do i = 1, imt
        if(temp2d_simp(i,j)<10000.) hs(i,j,nsp)=temp2d_simp(i,j)
   enddo
   enddo

!print *,hs(1:imt,1:jmt,nsp)


      ! Depth at U, V, T points as 2D arrays
   allocate ( abyst(imt, jmt)  )
   
   abyst = sum(dzt0(:,:,:), dim=3)
   
!   print *,abyst
   
!   stop 3957


   ! Calculate SSH/depth
   where (abyst /= 0)
      zstot = hs(:imt,:jmt,nsp)/abyst + 1
   elsewhere
      zstot = 0.d0
   end where



   ! Read temperature 
#if defined tempsalt 
   xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'votemper')
   tem(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
   ! Read salinity
   xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'vosaline')
   sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
   ! Calculate potential density
!   depthzvec = 0.
!   do j=1,jmt
!      latvec=-80+1./12.*float(j+subGridJmin-1)
!      do i=1,IMT
!         tmpzvec = tem(i,j,:,nsp)
!         salzvec = sal(i,j,:,nsp)
!         call statvd(tmpzvec, salzvec, rhozvec ,km ,depthzvec ,latvec)
!         rho(i,j,:,nsp)=rhozvec - 1000.
!      end do
!   end do
#endif     
   
   ! Read u, v
   
  ierr=NF90_OPEN(trim(fieldFile)//'U.nc',NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 7751
  ierr=NF90_INQ_VARID(ncid,'u_vtrans',varid)
  if(ierr.ne.0) stop 7763
  ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
  if(ierr.ne.0) stop 7799
  
  uvel(1:imt,:,:)= temp3d_simp(:,:,:)
  
!  print *,nctpos,uvel(1,1,1)
!  print *,'aaa'

!  stop 4964
!do j=10,1,-1
!write (*,"(200i1)") (int(1.e20* uvel(i,j,1)),i=1,imt)
!write (*,"(200i1)") (100*kmu(i,j),i=1,imt)
!enddo


!  stop 4965
!  kmu=99.
!   do j = 1, jmt
!   do i = 1, imt
!        if(uvel(i,j,1)==0.) kmu(i,j)=0.
!   enddo
!   enddo

!   do j=10,1,-1
!!    write (*,"(200i1)") (int(100* temp2d_simp(i,j)),i=1,imt)
!    write (*,"(200i1)") (kmu(i,j),i=1,imt)
!   enddo
  
!  stop 395
  
  ierr=NF90_OPEN(trim(fieldFile)//'V.nc',NF90_NOWRITE,ncid)
  if(ierr.ne.0) stop 8751
  ierr=NF90_INQ_VARID(ncid,'v_vtrans',varid)
  if(ierr.ne.0) stop 8763
  ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
  if(ierr.ne.0) stop 8799

  vvel(1:imt,:,:)= temp3d_simp(:,:,:)

 
   !
   ! Calculate zonal and meridional volume flux
   !
   ! Weight by (1 + ssh / depth)
   ! This is only an approximation of what NEMO really does
   ! but is accurate within 1% 
   !

   
   uflux(:,:,:,nsp) = 0.
   vflux(:,:,:,nsp) = 0.
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      dzt(i,j,k,nsp) = dzt0(i,j,k) * zstot(i,j)
      if(k <= kmu(i,j) ) uflux(i,j,km+1-k,nsp) = uvel(i,j,k)
      if(k <= kmv(i,j) ) vflux(i,j,km+1-k,nsp) = vvel(i,j,k)
   enddo
   enddo
   enddo
   
!     print *,'dzt0 and kmt'
!     print *,dzt0(1,1,:)

do j=10,1,-1
!write (*,"(200f5.0)") (dzt0(i,j,km),i=1,10)
!write (*,"(200i1)") (int(1.e20* dzt0(i,j,km)),i=1,imt)
!write (*,"(200i1)") (int(1.e20* dzt(i,j,km, nsp)),i=1,imt)
!write (*,"(200i1)") (100*kmt(i,j),i=1,imt)
enddo

   
   
   
   ! Check that volume fluxes are zero below sea floor
   do j=1,JMT-1
   do i=1,IMT-1
   do k=KM,1,-1
   if(k > kmv(i,j) .and. vflux(i,j,km+1-k,nsp) /= 0.) then
      print *,'vflux=',vflux(i,j,km+1-k,nsp),vvel(i,j,k),i,j,k,kmv(i,j),nsp
      stop 4966
   endif
   if(k > kmu(i,j) .and. uflux(i,j,km+1-k,nsp) /= 0.) then
      print *,'uflux=',uflux(i,j,km+1-k,nsp),uvel(i,j,k),i,j,k,kmu(i,j),nsp
      stop 4967
   endif
   if(k <= kmt(i,j) .and. dzt(i,j,km+1-k,nsp) == 0.) then
      print *,'dzt =', dzt(i,j,km+1-k,nsp),i,j,k,kmt(i,j),nsp
      stop 4968
   endif
   enddo
   enddo
   enddo
   
   
   return
   
end subroutine readfields



