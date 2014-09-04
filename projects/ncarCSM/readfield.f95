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

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk

  ! = Variables for filename generation
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2 ,subYr
  INTEGER                                    :: filePos ,fileJD ,subYrD
  INTEGER                                    :: yr1 ,mn1 ,dy1
  INTEGER                                    :: yr2 ,mn2 ,dy2
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: fieldFile
  CHARACTER (len=50)                         :: varName
  
  call updateClock
  call datasetswap

  ! === update the time counting ===
  filePos = mod((ints-1),fieldsPerFile)+1
  fileJD  = mod(floor(real((ints-1))/real(fieldsPerFile))+2,4)+1
  subYr   = mod(floor(real(ints-1)/real(fieldsPerFile)),4)+1
  dstamp  = 'BEC.gx3.22.pi.cv2.Ar.daily.XXXX-XX.nc'
  write (dstamp(28:31),'(i4.4)') currYear
  write (dstamp(33:34),'(i2.2)') currMon
  ncTpos    = 1 !intpart1
 
  fieldFile = trim(inDataDir)//trim(dstamp)
  fieldFile = trim(inDataDir)//'ccsm3t31savedd.pop.h.0002-03-01.nc'

  uvel = get3DfieldNC(trim(fieldFile), 'UVEL') / 100
  vvel = get3DfieldNC(trim(fieldFile), 'VVEL') / 100
  wvel = get3DfieldNC(trim(fieldFile), 'WVEL') / 100

  where (uvel>1000)
     uvel = 0
  end where
  where (vvel>1000)
     vvel = 0
  end where
  where (wvel>1000)
     wvel = 0
  end where

  do k=1,km
     kk = km + 1 - k
     uflux(:,:,k,2) = uvel(:,:,kk) * dyu(:,:) * dzt(:,:,k)
     vflux(:,:,k,2) = vvel(:,:,kk) * dxv(:,:) * dzt(:,:,k)
     wflux(:,:,k,2) = wvel(:,:,kk) * dxv(:,:) * dyu(:,:)
  end do
  return

end subroutine readfields
