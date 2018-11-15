SUBROUTINE readfields
  
  !USE mod_param
  !USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile
  USE mod_seed, only : nff
  
#ifdef tempsalt
  USE mod_dens
#endif
  IMPLICIT none

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos

  ! = Variables for filename generation
  CHARACTER (len=200)                        :: filename, dstamp
    REAL*4, ALLOCATABLE, DIMENSION(:,:)      :: ssh
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: rhof

  REAL, ALLOCATABLE, DIMENSION(:,:)          :: gridLat ,gridLon  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v,gridDX,e2u,gridDY
  !REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: dzu,dzv,dzt
  
  logical around
  
  alloCondGrid: if ( .not. allocated (e1v) ) then
     allocate ( gridLat(IMT,JMT) ,gridLon(IMT,JMT) )
     allocate ( e1v(IMT+2,JMT)    ,gridDX(IMT+2,JMT) )
     allocate ( e2u(IMT+2,JMT)    ,gridDY(IMT+2,JMT) )
     !allocate ( dzu(IMT+2,JMT,KM) ,dzv(IMT+2,JMT,KM),dzt(IMT+2,JMT,KM) )
  end if alloCondGrid
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
     allocate ( rhof(IMT+2,JMT,KM) )
  end if alloCondUVW

  call datasetswap
  call updateClock
  ! === update the time counting ===
  if (currJDyr==366) currJDyr=365
  dstamp   = 'archv.0000_000_00_'
  dstamp   = 'hycom_glb_909_YYYYMMDD00_t000_uv3z.nc'        

  write (dstamp(15:18),'(i4.4)') int(curryear)
  write (dstamp(19:20),'(i2.2)') int(currmon)
  write (dstamp(21:22),'(i2.2)') int(currday)
  
  if (currjdtot > 735096) write (dstamp(11:13),'(A3)') "910"

  filename = trim(inDataDir)//trim(dstamp)
  uvel = get3DfieldNC(trim(filename) ,'water_u') * 0.001 * 2
  where (uvel<-10) uvel=0
  vvel = get3DfieldNC(trim(filename) ,'water_v') * 0.001 * 2 
  where (vvel<-10) vvel=0
  if (nff .eq. -1) then
     uvel = -uvel
     vvel = -vvel
  end if

  hs(:,:,2) = 0.01*ssh
  hs(imt+1,:,2) =hs(1,:,2)
  hs(:,jmt+1,2) =hs(:,1,2)
  
  do k=1,km
     kk=km+1-k
     uflux(:,:,k,2) = uvel(:,:,kk) * dyu(:,:) * dz(k) !dzu(:,:,kk)
     vflux(:,:,k,2) = vvel(:,:,kk) * dxv(:,:) * dz(k) !dzv(:,:,kk)
   !  rho  (:,:,k,2) = rhof(:,:,kk)
  end do

  !uflux(:,:,km,2) = uvel(:,:,1) * e2u * (dzu(:,:,1,1) & 
  !     + 0.5*(hs(:,:,2) + hs(2:imt+1,:,2)))
  !vflux(:,:,km,2) = vvel(:,:,1) * e1v * (dzv(:,:,1,1) & 
  !     + 0.5*(hs(:,:,2) + hs(:,2:jmt+1,2)))
  !rho(:,:,km,2)   = rhof(:,:,1)



  return

end subroutine readfields
