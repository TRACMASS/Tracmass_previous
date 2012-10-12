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

  CHARACTER(LEN=80)                          :: fstamp
  INTEGER                                    :: k
 
  call datasetswap
  call updateClock
 
  fstamp = '.1440x720x50.XXXXXXXX.nc'
  write (fstamp(14:17),'(I4.4)') int(loopYear)
  write (fstamp(18:19),'(I2.2)') int(loopMon)
  write (fstamp(20:21),'(I2.2)') int(loopDay)

  uvel = get3DfieldNC(trim(trim(inDataDir)//'/UVEL'//fstamp), 'UVEL')
  vvel = get3DfieldNC(trim(trim(inDataDir)//'/VVEL'//fstamp), 'VVEL')

  do k=1,km
     uflux(:,:,km-k+1,2) = uvel(:,:,k) * dz(k) * dyu
     vflux(:,:,km-k+1,2) = vvel(:,:,k) * dz(k) * dxv
  end do
 
end subroutine readfields
