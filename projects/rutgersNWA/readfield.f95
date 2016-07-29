SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel

  USE mod_seed,     only: nff
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
  REAL,       ALLOCATABLE, DIMENSION(:)      :: sc_r,Cs_r
  INTEGER                                    :: hc

  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh,dzt0
  ! ===   ===   ===

  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt), dzt0(imt,jmt) )
     allocate ( sc_r(km), Cs_r(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km,2), dzv(imt,jmt,km,2) )
  end if alloCondDZ
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_r = 0
  Cs_r = 0

  call datasetswap
  call updateClock

  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  dstamp      = 'nwa_avg_XXXXX.nc'

  write (dstamp(9:13),'(I5.5)') & 
       int(currJDtot) - 714782
  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
  uvel        = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel        = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh         = get2dfieldNC(trim(dataprefix) ,'zeta')
  hs(:,:,2)   = ssh

  where (uvel > 1000)
     uvel = 0
  end where
  where (vvel > 1000)
     vvel = 0
  end where

  Cs_r = get1DfieldNC (trim(dataprefix), 'Cs_w')
  sc_r = get1DfieldNC (trim(dataprefix), 's_w')
  hc   = getScalarNC (trim(dataprefix), 'hc')

  do k=1,km
     dzt0 = (sc_r(k)-Cs_r(k))*hc + Cs_r(k) * depth
     dzt(:,:,k,2)= dzt0 + ssh*(1.0 + dzt0/depth)
  end do

  dzt0 = dzt(:,:,km,2)
  dzt(:,:,1:km-1,2)=dzt(:,:,2:km,2)-dzt(:,:,1:km-1,2)
  dzt(:,:,km,2) = ssh-dzt0

  dzu(1:imt-1,:,:,2) = dzt(1:imt-1,:,:,2)*0.5 + dzt(2:imt,:,:,2)*0.5
  dzv(:,1:jmt-1,:,2) = dzt(:,1:jmt-1,:,2)*0.5 + dzt(:,2:jmt,:,2)*0.5

  do k=1,km
     uflux(:,:,k,2)   = uvel(:,:,k) * dzu(:,:,k,2) * dyu
     vflux(:,:,k,2)   = vvel(:,:,k) * dzv(:,:,k,2) * dxv
  end do

  !print *, dzv(:,:,km,:)
  !stop
  
  if (nff .eq. 2) then
     uflux = -uflux
     vflux = -vflux
  end if

  return
end subroutine readfields
