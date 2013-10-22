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
  
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)     :: ssh
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)     :: u0mask,u2mask
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:)     :: v0mask,v2mask
  REAL,    SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: dzu,dzv
  
  REAL                                    :: x_scale     = 9.1556e-5
  REAL                                    :: x_offset    =-3.3845e-19
  REAL                                    :: y_scale     = 9.1556e-5
  REAL                                    :: y_offset    =-5.4561e-18
  REAL                                    :: temp_scale  = 0.0005340739
  REAL                                    :: temp_offset = 12.5
  REAL                                    :: salt_scale  = 0.0006103702
  REAL                                    :: salt_offset = 20.0
  REAL                                    :: ssh_scale   = 0.000244
  REAL                                    :: ssh_offset  = 0
  
  integer                                    :: year,month,day  
  logical                                    :: around
  character*16                               :: fileName
  real, dimension(22) :: dS = (/ 0.008, 0.008, 0.017, 0.033, 0.067, &
       0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, 0.067, &
       0.067, 0.067, 0.067, 0.033, 0.017, 0.008, 0.008, 0.000 /)
  
  if ( .NOT. ALLOCATED(dzu) ) then
     allocate ( dzu(imt,jmt,km),dzv(imt,jmt,km) )
     allocate ( u0mask(imt,jmt),u2mask(imt,jmt) )
     allocate ( v0mask(imt,jmt),v2mask(imt,jmt) )

     maskfile= trim(inDataDir) // 'uvmasks.cdf'
     u0mask = - ( get2DfieldNC(trim(maskfile) ,'u0mask') - 1) * 9999
     u2mask =   ( get2DfieldNC(trim(maskfile) ,'u2mask') - 1) * 9999
     v0mask = - ( get2DfieldNC(trim(maskfile) ,'v0mask') - 1) * 9999
     v2mask =   ( get2DfieldNC(trim(maskfile) ,'v2mask') - 1) * 9999
  end if
  allocate ( ssh(imt,jmt) )

  call datasetswap

  intpart1 = mod(ints+1,8)
  if (intpart1 .eq. 0) then
     intpart1=8
  endif
  intpart2 = floor((ints)/8.)+1
  ndates   = intpart2
  call gdate (2453005+1+ndates, year,month,day)
  fileName = 'unh.20000000.cdf'
  write(fileName(5:12),'(i4i2.2i2.2)') year,month,day
  ncFile   = trim(inDataDir)//fileName

  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556
  nread = mod(ints/5,18) + 1

  ! === Velocities ===
  !Use  t=1  i=2  j=3  k=4
  map3d    = [2, 3, 4, 1]   
  uvel =  get3DfieldNC(trim(ncFile), 'u') * x_scale + x_offset
  vvel =  get3DfieldNC(trim(ncFile), 'v') * y_scale + y_offset

  do k=1,km
     uvel(:,:,k) =  uvel(:,:,k) * cos(ang) + vvel(:,:,k)*sin(ang)
     vvel(:,:,k) = -uvel(:,:,k) * sin(ang) + vvel(:,:,k)*cos(ang)
  end do
  if (intstep .le. 0) then
     uvel = -vvel
     vvel = -vvel
  end if

!!$  mask: do k=1,km
!!$     uvel(:,:,k) = uvel(:,:,k) * u0mask
!!$     uvel(:,:,k) = uvel(:,:,k) * u2mask
!!$     vvel(:,:,k) = vvel(:,:,k) * v0mask
!!$     vvel(:,:,k) = vvel(:,:,k) * v2mask
!!$  end do mask

  uvmask: do k=1,km
     uvel(:,:,k)=min(uvel(:,:,k),real(u0mask))
     uvel(:,:,k)=max(uvel(:,:,k),real(u2mask))
     vvel(:,:,k)=min(vvel(:,:,k),real(v0mask))
     vvel(:,:,k)=max(vvel(:,:,k),real(v2mask))
  end do uvmask

  !sal(:,:,:,2) = get3dfieldNC(trim(ncFile) ,'salt') 
  !sal(:,:,:,2) = sal(:,:,:,2) * salt_scale + salt_offset
  ssh = get2DfieldNC(trim(ncFile) ,'elev') * ssh_scale + ssh_offset

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
    uflux(:,1:80,k,2)  =uvel(:,:,kk) * dyu(:,:) * dzu(:,:,k)
    vflux(:,1:80,k,2)  =vvel(:,:,kk) * dxv(:,:) * dzv(:,:,k)
  !  rho(:,1:80,k,2)=NCfieldr(1,:,:,k)
  enddo

end subroutine readfields
