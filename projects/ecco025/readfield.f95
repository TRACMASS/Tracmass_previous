SUBROUTINE readfields
  USE mod_time
  USE mod_name
  USE mod_getfile
  USE mod_vel
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT NONE

  CHARACTER(LEN=80)                          :: fstamp
  INTEGER                                    :: i,j,k
  
  call datasetswap
  call updateClock

  fstamp = '.1440x720x50.XXXXXXXX.nc'
  write (fstamp(14:17),'(I4.4)') int(loopYear)
  write (fstamp(18:19),'(I2.2)') int(loopMon)
  write (fstamp(20:21),'(I2.2)') int(loopDay)
  
  uvel = get3DfieldNC(trim(trim(inDataDir)//'/UVEL'//fstamp), 'UVEL')
  vvel = get3DfieldNC(trim(trim(inDataDir)//'/VVEL'//fstamp), 'VVEL')


  forall(i=1:imt, j=1:jmt) dzt(i,j,:,nsp) = dz
  
  do k=1,km
     uflux(1:imt,1:jmt,km-k+1,nsp) = (uvel(1:imt, 1:jmt, k)/2  + &
                                      uvel(1:imt, 1:jmt, k)/2) * &
                                      dz(km-k+1) * dyu(:,:)
     vflux(1:imt,1:jmt,km-k+1,nsp) = (vvel(1:imt, 1:jmt, k)/2  + &
                                      vvel(1:imt, 1:jmt, k)/2) * & 
                                      dz(km-k+1) * dxv(:,:)
  end do
end subroutine readfields
