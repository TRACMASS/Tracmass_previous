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
  
  if ( .NOT. ALLOCATED(dzu) ) then
     allocate ( dzu(imt,jmt,km),dzv(imt,jmt,km) )
  end if


  call datasetswap

  if (ncTpos >= fieldsPerFile) ncTpos=0
  intpart2 = floor((ints)/float(iter))
  ndates   = intpart2
  call gdate (2449944 + ndates, year,month,day)
  fileName = 'oscar_vel0000.nc'
  write(fileName(10:13),'(i4)') year
  ncFile   = trim(inDataDir)//fileName
  inquire(file=ncFile,exist=around)
  if(.not.around) stop 4556
  nread = mod(ints/5,18) + 1
  ncTpos = ncTpos + 1

  ! === Velocities ===
  !Use  t=1  i=2  j=3  k=4
  map2d    = [1, 2, 3, 4]   

  uvel(:,:,1) =  cshift(get2DfieldNC(trim(ncFile), 'u'),481,1)
  vvel(:,:,1) =  cshift(get2DfieldNC(trim(ncFile), 'v'),481,1)

  where (uvel .ne. uvel)
     uvel = 0
  end where
  where(vvel .ne. vvel)
     vvel = 0
  end where

  uflux(:,:,1,2)  =uvel(:,:,1) * dyu(:,:) * dz(1)
  vflux(:,:,1,2)  =vvel(:,:,1) * dxv(:,:) * dz(1)

  uflux(1:imt-1,:,1,2) = (uflux(1:imt-1,:,1,2) + uflux(2:imt,:,1,2)) / 2
  uflux(imt    ,:,1,2) = (uflux(imt    ,:,1,2) + uflux(1    ,:,1,2)) / 2
  vflux(:,1:jmt-1,1,2) = (vflux(:,1:jmt-1,1,2) + uflux(:,2:jmt,1,2)) / 2

end subroutine readfields
