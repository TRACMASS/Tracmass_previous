SUBROUTINE readfields
  
  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile


#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none

  CHARACTER(len=200), SAVE                   :: ncFile, maskfile
  integer                                    :: i,im,j,jm,k,kk,ints2
  INTEGER, SAVE                              :: nread,ndates
  
  INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: ssh
    
  integer                                    :: year,month,day  
  logical                                    :: around
  character(len=100)                         :: fileName
 
  allocate ( ssh(imt,jmt) )

  call datasetswap

  ! === NetCDF file and fields ===
  fileName = '22450101.ocean_daily.nc'
  ncFile   = trim(inDataDir)//fileName
  ncTpos = ints - 17459
  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556

  ! === Velocities ===
  uvel =  get3DfieldNC(trim(ncFile), 'u')
  vvel =  get3DfieldNC(trim(ncFile), 'v')
  if (intstep .le. 0) then
     uvel = -vvel
     vvel = -vvel
  end if
  where (vvel== -10) 
     uvel = 0
     vvel = 0
  end where 

  !print *,kmt(200:210,32)
  !print *,uvel(200:210,32,1)
  !print *,vvel(200:210,32,1)

  !sal(:,:,:,2) = get3dfieldNC(trim(ncFile) ,'salt') 
  !sal(:,:,:,2) = sal(:,:,:,2) * salt_scale + salt_offset
  !Use  t=1  i=2  j=3  k=4
  map2d    = [2, 3, 1, 4]
  hs(:,:,2) = get2DfieldNC(trim(ncFile) ,'eta_t')

  dzu(:,:,km) = dzt0surf + hs(:,:,2)
  dzv(:,:,km) = dzt0surf + hs(:,:,2)

  do k=1,km
    kk=km-k+1
    uflux(2:imt,2:jmt,k,2)  =  &
         (uvel(1:imt-1,2:jmt,kk) + uvel(1:imt-1,1:jmt-1,kk))/2 * & 
         dyu(1:imt-1,2:jmt) * dzt(1:imt-1,2:jmt,k)
    vflux(2:imt,2:jmt,k,2)  = &
         (vvel(2:imt,1:jmt-1,kk) + vvel(1:imt-1,1:jmt-1,kk))/2 * & 
         dxv(2:imt,1:jmt-1) * dzt(2:imt,1:jmt-1,k)
    uflux(1,2:jmt,k,2) =  &
         (uvel(1,2:jmt,kk) + uvel(imt,1:jmt-1,kk))/2 * & 
         dyu(1,2:jmt) * dzt(imt,2:jmt,k)
    vflux(1,2:jmt,k,2)  = &
         (vvel(1,1:jmt-1,kk) + vvel(1,1:jmt-1,kk))/2 * & 
         dxv(1,1:jmt-1) * dzt(1,1:jmt-1,k)
  !  rho(:,1:80,k,2)=NCfieldr(1,:,:,k)
  enddo

end subroutine readfields
