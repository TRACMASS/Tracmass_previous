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

  ! = Variables for filename generation 
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  ! = Loop variables
  INTEGER                                    :: k ,kk
  ! = Variables for converting from S to Z
  REAL,       ALLOCATABLE, DIMENSION(:)      :: stoz
  ! = Input fields from GCM
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: ssh
  
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
     allocate ( stoz(km) )
  end if alloCondUVW
  alloCondDZ: if(.not. allocated (dzu)) then
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  call datasetswap
  call updateClock

  dataprefix  =  trim(indatadir)//"ecom.cdf"
  nctpos = ints-17533082+1

  uvel = get3DfieldNC(trim(dataprefix) ,   'u') * 9.155553e-05
  vvel = get3DfieldNC(trim(dataprefix) ,   'v') * 9.155553e-05
  !wvel = get3DfieldNC(trim(dataprefix) ,   'w')
  ssh  = get2dfieldNC(trim(dataprefix) ,'elev') * 0.000122074
  stoz = get1DfieldNC (trim(dataprefix), 'sigma')

uvel(1:imt-1, :, :) = uvel(2:imt ,:, :)
vvel(:, 1:jmt-1, :) = vvel(:, 2:jmt, :)

  do k=1,km
     dzt(:,:,k) = (ssh + depth) * stoz(km-k+1)
  end do
  dzt(:,:,2:km)=dzt(:,:,2:km)-dzt(:,:,1:km-1)
  dzt(:,:,1) = 0
  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:)*0.5 + dzt(2:imt,:,:)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:)*0.5 + dzt(:,2:jmt,:)*0.5

  do k=1,km
     uflux(:,:,k,2) = uvel(:,:,km-k+1) * dzu(:,:,k) * dyu
     vflux(:,:,k,2) = vvel(:,:,km-k+1) * dzv(:,:,k) * dxv
     !wflux(:,:,k,2) = wvel(:,:,k) * dyu * dxv 
  end do
  hs(:,:,2) = ssh

  return

end subroutine readfields
