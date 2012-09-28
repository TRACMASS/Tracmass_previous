SUBROUTINE readfields
  !USE mod_param
  USE mod_time
  !USE mod_grid
  USE mod_name
  USE mod_getfile
  USE mod_vel
  
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT NONE

  CHARACTER(LEN=80)                         :: fstamp
  INTEGER                                   :: k, kk
 
  call datasetswap
  call updateClock
 
  fstamp = '.1440x720x50.XXXXXXXX.nc'
  write (fstamp(14:17),'(I4.4)') int(loopYear)
  write (fstamp(18:19),'(I2.2)') int(loopMon)
  write (fstamp(20:21),'(I2.2)') int(loopDay)

  uvel = get3DfieldNC(trim(trim(inDataDir)//'/UVEL'//fstamp), 'UVEL')
  vvel = get3DfieldNC(trim(trim(inDataDir)//'/VVEL'//fstamp), 'VVEL')

  do k=1,km
     kk=km-k+1
     uflux(:,:,k,2)   = uvel(:,:,kk) * dz(kk) * dyu
     vflux(:,:,k,2)   = vvel(:,:,kk) * dz(kk) * dxv
  end do
  if (intstep .le. 0) then
     uflux = -uflux
     vflux = -vflux
  end if

end subroutine readfields
