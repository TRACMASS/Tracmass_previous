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

  CHARACTER(len=200), SAVE                   :: ncFile, maskfile
  integer                                    :: i,im,j,jm,k,kk,ints2
  INTEGER, SAVE                              :: nread,ndates
  integer                                    :: intpart1,intpart2
  
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)     :: ssh
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: dzu,dzv
    
  integer                                    :: year,month,day  
  logical                                    :: around
  character(len=100)                         :: fileName
 
  if ( .NOT. ALLOCATED(dzu) ) then
     allocate ( dzu(imt,jmt,km),dzv(imt,jmt,km) )
  end if
  allocate ( ssh(imt,jmt) )

  call datasetswap

  intpart1 = mod(ints+1,8)
  if (intpart1 .eq. 0) then
     intpart1=8
  endif
  intpart2 = floor((ints)/8.)+1
  ndates   = intpart2
  fileName = '22450101.ocean_daily.nc'
  ncFile   = trim(inDataDir)//fileName

  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556
  nread = mod(ints/5,18) + 1

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

  !sal(:,:,:,2) = get3dfieldNC(trim(ncFile) ,'salt') 
  !sal(:,:,:,2) = sal(:,:,:,2) * salt_scale + salt_offset
  ssh = 0 !get2DfieldNC(trim(ncFile) ,'elev')

  do  k=1,km
     dzt(:,:,km-k+1) = depth(:,:) 
  end do  

  do i=1,imt-1
     dzv(i,:,:)=0.5*(dzt(i,:,:)+dzt(i+1,:,:))
  enddo
  dzv(imt,:,:)=dzv(imt-1,:,:)

  do j=1,jmt-1
     dzu(:,j,:)=0.5*(dzt(:,j,:)+dzt(:,j+1,:))
  enddo
  dzu(:,jmt,:)=dzu(:,jmt-1,:)
  do k=1,km
    kk=km-k+1
    uflux(:,1:80,k,2)  =uvel(:,:,kk) * dyu(:,:) * dzu(:,:,k)
    vflux(:,1:80,k,2)  =vvel(:,:,kk) * dxv(:,:) * dzv(:,:,k)
  !  rho(:,1:80,k,2)=NCfieldr(1,:,:,k)
  enddo

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
