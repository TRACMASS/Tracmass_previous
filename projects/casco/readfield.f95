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

  CHARACTER(len=63), SAVE                    :: ncFile, maskfile
  integer                                    :: i,im,j,jm,k,kk,ints2
  INTEGER, SAVE                              :: nread,ndates
  integer                                    :: intpart1,intpart2
  
  real,         ALLOCATABLE, DIMENSION(:,:)  :: ssh
  real,       ALLOCATABLE, DIMENSION(:,:,:)  :: utemp

  REAL                                       :: x_scale     = 9.1555528e-05
  REAL                                       :: x_offset    = 0
  REAL                                       :: y_scale     = 9.1555528e-05
  REAL                                       :: y_offset    = 0
  REAL                                       :: temp_scale  = 0.0005340739
  REAL                                       :: temp_offset = 12.5
  REAL                                       :: salt_scale  = 0.0006103702
  REAL                                       :: salt_offset = 20.0
  REAL                                       :: ssh_scale   = 0.00024414808
  REAL                                       :: ssh_offset  = 0
  
  integer                                    :: year,month,day  
  logical                                    :: around
  character*18                               :: fileName
  real, dimension(22) :: dS = (/ 0.008, 0.008, 0.017, 0.033, 0.067, &
       0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, &
       0.067, 0.067, 0.067, 0.033, 0.017, 0.008, 0.008, 0.000 /)
   
  allocate ( ssh(imt,jmt), utemp(imt,jmt,km) ) 
  if ( .not. allocated (dzv) ) then
     allocate ( dzv(imt,jmt,km),dzu(imt,jmt,km)  )
  end if

  call datasetswap
  call updateClock



  intpart1 = mod(ints+1,8)
  if (intpart1 .eq. 0) then
     intpart1=8
  endif
  intpart2 = floor((ints)/8.)+1
  ndates   = intpart2
  fileName = 'casco.20000000.cdf'
  write(fileName(7:14),'(i4i2.2i2.2)') currYear,currMon,currDay
  ncFile   = trim(inDataDir)//fileName
  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556
  nread = mod(ints/5,18) + 1
  ncTpos = mod(ints+1,8)+1
  print *,"File and tpos: ", filename, ncTpos

  ! === Velocities ===
  !Use  t=1  i=2  j=3  k=4
  map2d    = [2, 3, 4, 1]   
  map3d    = [2, 3, 4, 1]   
  uvel =  get3DfieldNC(trim(ncFile), 'u') * x_scale + x_offset
  vvel =  get3DfieldNC(trim(ncFile), 'v') * y_scale + y_offset
  ssh  = get2DfieldNC(trim(ncFile) ,'elev') * ssh_scale + ssh_offset

  do k=1,km
     utemp(:,:,k) = ( uvel(:,:,k) * cos(ang) + vvel(:,:,k)*sin(ang) ) * mask
     vvel(:,:,k) = (-uvel(:,:,k) * sin(ang) + vvel(:,:,k)*cos(ang) ) * mask
  end do
  uvel = utemp

  do  k=1,km
     dzt(:,:,km-k+1) = dS(k) * (depth(:,:) + ssh(:,:))
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
    uflux(:,:,k,2)  =uvel(:,:,kk) * dyu(:,:) * dzu(:,:,k)
    vflux(:,:,k,2)  =vvel(:,:,kk) * dxv(:,:) * dzv(:,:,k)
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
