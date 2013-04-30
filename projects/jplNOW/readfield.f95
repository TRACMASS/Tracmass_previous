SUBROUTINE readfields

  USE netcdf
  USE mod_time, only: ints, intstep, updateclock, & 
       currYear,currMon,currDay,currHour
  USE mod_name, only: inDataDir
  USE mod_grid, only: imt, jmt, km, kmt, mask, depth
  USE mod_grid, only: dyu, dxv, dzu, dzv, dzt, dxdy
  USE mod_vel, only: uflux, vflux, uvel, vvel, hs, datasetswap
  USE mod_getfile, only: getScalarNC,get1DfieldNC, get2dfieldNC, get3dfieldNC
  
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
  INTEGER                                    :: ncid, ierr

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
     allocate ( dzu(imt,jmt,km), dzv(imt,jmt,km) )
  end if alloCondDZ
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  sc_r = 0
  Cs_r = 0

  call datasetswap
  call updateClock

  ! === update the time counting ===
  intpart1    = mod(ints,24)
  intpart2    = floor((ints)/24.)
  dstamp      = '2011083121_da.nc'

  write (dstamp(1:10),'(i4i2.2i2.2i2.2)') & 
       currYear,currMon,currDay,currHour+3
  dataprefix  = trim(inDataDir) // dstamp
  tpos        = intpart1+1
 
  uvel      = get3DfieldNC(trim(dataprefix) ,   'u')
  vvel      = get3DfieldNC(trim(dataprefix) ,   'v')
  ssh       = get2dfieldNC(trim(dataprefix) ,'zeta')
  hs(:,:,2) = ssh

  ierr = NF90_OPEN(trim(dataprefix) ,NF90_NOWRITE ,ncid)
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'sc_w', sc_r)
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'Cs_w', Cs_r)
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'hc', hc)
  do k=1,km
     dzt0 = (sc_r(k)-Cs_r(k))*hc + Cs_r(k) * depth
     dzt(:,:,k)= dzt0 + ssh*(1.0 + dzt0/depth)
  end do

  dzt0 = dzt(:,:,km)
  dzt(:,:,1:km-1)=dzt(:,:,2:km)-dzt(:,:,1:km-1)
  dzt(:,:,km) = ssh-dzt0

  dzu(1:imt-1,:,:) = dzt(1:imt-1,:,:)*0.5 + dzt(2:imt,:,:)*0.5
  dzv(:,1:jmt-1,:) = dzt(:,1:jmt-1,:)*0.5 + dzt(:,2:jmt,:)*0.5

  do k=1,km
     kk=km+1-k
     uflux(:,:,k,2)   = uvel(:,:,k) * dzu(:,:,k) * dyu
     vflux(:,:,k,2)   = vvel(:,:,k) * dzv(:,:,k) * dxv
  end do

  if (intstep .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if
  return
 
end subroutine readfields
