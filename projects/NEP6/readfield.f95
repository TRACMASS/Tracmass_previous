SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  
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
  REAL*8,       ALLOCATABLE, DIMENSION(:)      :: sc_r,Cs_r
  REAL*8,       ALLOCATABLE, DIMENSION(:)      :: sc_w,Cs_w
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
  dstamp      = 'nep6_avg_XXXXX.nc'

  print *,currJDtot

  write (dstamp(10:14),'(I5.5)') & 
       int(currJDtot) - 731576
  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
  print *,dataprefix

  uvel      = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel      = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh       = get2dfieldNC(trim(dataprefix) ,'zeta')
  where (uvel > 1000)
     uvel = 0
  end where
  where (vvel > 1000)
     vvel = 0
  end where
  where (ssh > 1000)
     ssh = 0
  end where
  hs(:imt,:jmt,2) = ssh


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

  do k=1,km
     dzt0 = (hc*sc_r(k) + depth*Cs_r(k)) / (hc + depth)
     z_r(:,:,k) = ssh + (ssh + depth) * dzt0
     dzt0 = (hc*sc_w(k) + depth*Cs_w(k)) / (hc + depth)
     dzt(:,:,k) = ssh + (ssh + depth) * dzt0
  end do
  dzt0 = dzt(:,:,km)
  dzt(:,:,1:km-1)=dzt(:,:,2:km)-dzt(:,:,1:km-1)
  dzt(:,:,km) = ssh - dzt0
  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:)*0.5 + dzt(2:imt,:,:)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:)*0.5 + dzt(:,2:jmt,:)*0.5

  do k=1,km
     uflux(:,:,k,2)   = uvel(:imt,:,k) * dzu(:,:,k) * dyu(:imt,:)
     vflux(:,1:jmt,k,2)   = vvel(:imt,:,k) * dzv(:,:,k) * dxv(:imt,:)
  end do

  if (intstep .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if

#ifdef tempsalt
  tem(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'temp')
  sal(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'salt')
  rho(:,:,:,2)      = get3DfieldNC(trim(dataprefix) ,   'rho')
#endif

  return

 !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

 
contains

  
  subroutine datasetswap
    hs(:,:,1)      = hs(:,:,2)
    uflux(:,:,:,1) = uflux(:,:,:,2)
    vflux(:,:,:,1) = vflux(:,:,:,2)
#ifdef explicit_w
    wflux(:,:,:,1) = wflux(:,:,:,2)
#endif

#ifdef tempsalt
    tem(:,:,:,1)   = tem(:,:,:,2)
    sal(:,:,:,1)   = sal(:,:,:,2)
    rho(:,:,:,1)   = rho(:,:,:,2)
#endif
  end subroutine datasetswap


end subroutine readfields
