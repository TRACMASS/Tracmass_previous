SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
!  USE mod_coord   !FC
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  USE mod_seed, only: nff! LD ADDED, for nff    
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2
  INTEGER                                    :: ndates
  INTEGER                                    :: yr1 ,mn1 ,dy1,hr
  INTEGER                                    :: yr2 ,mn2 ,dy2
  
  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  ! = Variables for converting from S to Z
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_r,Cs_r
  REAL*8,       ALLOCATABLE, DIMENSION(:)    :: sc_w,Cs_w
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL*8,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  ! ===   ===   ===

  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt), dzt0(imt,jmt) )
     allocate ( sc_r(km), Cs_r(km) )
     allocate ( sc_w(km), Cs_w(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ
  initFieldcond: if(ints.eq.intstart) then
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif
  end if initFieldcond

  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_r = 0
  Cs_r = 0
  sc_w = 0
  Cs_w = 0


  call datasetswap
  call updateClock

  ! === update the time counting ===
  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  dstamp      = 'coral_avg_XXXXX.nc'

  !write (dstamp(11:15),'(I5.5)') & 
       !int(currJDtot) - 731365
  write (dstamp(11:15),'(I5.5)') & 
       int(currJDtot) - 714777
  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
  print *, "curr JD:", currJDtot, "read file:", dataprefix

  uvel      = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel      = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh       = get2DfieldNC(trim(dataprefix) ,'zeta')
#ifdef explicit_w
  wvel      = get3DfieldNC(trim(dataprefix) ,'omega')
#endif
  where (uvel > 1000)
     uvel = 0
  end where
  where (vvel > 1000)
     vvel = 0
  end where
  where (ssh > 1000)
     ssh = 0
  end where
  !hs(:,:,2) = ssh
  hs(:imt,:jmt,2) = ssh(:imt,:jmt)

#ifdef explicit_w
  wflux(:,:,:,2) = 0.
  do j=1,jmt
    do i=1,imt
      wflux(i,j,0:km-1,2) = wvel(i,j,1:km)*dxdy(i,j)
    end do
  end do
  !print *, 'wflux=', wflux(410,255,50,2), wflux(410,255,49,2), wflux(410,255,48,2)
  !print *, 'wflux=', wflux(410,255,3,2), wflux(410,255,1,2), wflux(410,255,0,2)
#endif

  Cs_w = get1DfieldNC (trim(dataprefix), 'Cs_w')
  sc_w = get1DfieldNC (trim(dataprefix), 's_w')
  Cs_r = get1DfieldNC (trim(dataprefix), 'Cs_r')
  sc_r = get1DfieldNC (trim(dataprefix), 's_rho')
  hc   = getScalarNC (trim(dataprefix), 'hc')

  !do k=1,km
  !   dzt0 = (hc*sc_r(k) + depth*Cs_r(k)) / (hc + depth)
  !   dzt(:,:,k) = ssh + (ssh + depth) / dzt0
  !end do
  !dzt0 = dzt(:,:,km)
  !dzt(:,:,1:km-1)=dzt(:,:,2:km)-dzt(:,:,1:km-1)
  !dzt(:,:,km) = dzt0-ssh

  z_w(:,:,0,2) = depth
  do k=1,km
     dzt0 = (hc*sc_r(k) + depth*Cs_r(k)) / (hc + depth)
     z_r(:,:,k,2) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:)
     dzt0 = (hc*sc_w(k) + depth*Cs_w(k)) / (hc + depth)
#ifdef zgrid3Dt
     z_w(:,:,k,2) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:)
     dzt(:,:,k,2) = z_w(:,:,k,2)
#else
     dzt(:,:,k) = ssh(:imt,:) + (ssh(:imt,:) + depth(:imt,:)) * dzt0(:imt,:)
#endif
  end do
#ifdef zgrid3Dt
  dzt0 = dzt(:,:,km,2)
  dzt(:,:,1:km-1,2)=dzt(:,:,2:km,2)-dzt(:,:,1:km-1,2)
  dzt(:,:,km,2) = ssh(:imt,:) - dzt0
  dzt(:,:,:,1)=dzt(:,:,:,2)
  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5
#else
  dzt0 = dzt(:,:,km)
  dzt(:,:,1:km-1)=dzt(:,:,2:km)-dzt(:,:,1:km-1)
  dzt(:,:,km) = ssh(:imt,:) - dzt0
  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:)*0.5 + dzt(2:imt,:,:)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:)*0.5 + dzt(:,2:jmt,:)*0.5
#endif

  do k=1,km
     !uflux(:,:,k,2)   = uvel(:,:,k) * dzu(:,:,k) * dyu
     !vflux(:,:,k,2)   = vvel(:,:,k) * dzv(:,:,k) * dxv
     uflux(:,:,k,2)   = uvel(:imt,:,k) * dzu(:,:,k) * dyu(:imt,:)
     vflux(:,1:jmt,k,2)   = vvel(:imt,:,k) * dzv(:,:,k) * dxv(:imt,:)
  end do

  if (nff .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if

#ifdef tempsalt
  tem(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'temp')
  sal(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'salt')
  !rho(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'rho')
#ifdef larval_fish
  srflux(:,:,2)     = get2DfieldNC(trim(dataprefix) ,   'swrad')
! Note: this works as long as surface AKt is zero.
  ak2(:,:,:)        = get3DfieldNC(trim(dataprefix) ,   'AKt')
  akt(:,:,0:km-1,2) = ak2(:,:,:)
#endif
#endif

  return

end subroutine readfields
