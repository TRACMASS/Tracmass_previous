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

  ! = Loop variables
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos

  integer dimidx,dimidy,dimidz,dimidt       !output ID index of dimension
  integer varid,varidx,varidy,varidz,varidt !output ID index of variable
  integer lenx,leny,lenz,lent,lenz2         !output Length of dimension
  integer p, x1, y1, z1, t1 !?

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
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: rhof
    
  alloCondUVW: if(.not. allocated (rhof)) then
     allocate ( rhof(IMT+2,JMT,KM) )
  end if alloCondUVW

  ! === swap between datasets ===
  hs(:,:,1)=hs(:,:,2)
  uflux(:,:,:,1)=uflux(:,:,:,2)
  vflux(:,:,:,1)=vflux(:,:,:,2)
#ifdef tempsalt 
  tem(:,:,:,1)=tem(:,:,:,2)
  sal(:,:,:,1)=sal(:,:,:,2)
  rho(:,:,:,1)=rho(:,:,:,2)
#endif
  
  call updateClock
  
  ! === update the time counting ===
  filePos = mod((ints-1),fieldsPerFile)+1
  fileJD  = mod(floor(real((ints-1))/real(fieldsPerFile))+2,4)+1
  subYr   = mod(floor(real(ints-1)/real(fieldsPerFile)),4)+1
  dstamp  = 'BEC.gx3.22.pi.cv2.Ar.daily.XXXX-XX.nc'
  write (dstamp(28:31),'(i4.4)') currYear
  write (dstamp(33:34),'(i2.2)') currMon
  tpos    = intpart1
 
  fieldFile = trim(inDataDir)//trim(dstamp)

  uvel = get3DfieldNC(trim(fieldFile), 'UVEL') / 100
  vvel = get3DfieldNC(trim(fieldFile), 'VVEL') / 100

  do k=1,km
     kk=km+1-k
     uflux(:,:,k,2) = uvel(:,:,kk) * dyu(:,:) * dzt(:,:,k)
     vflux(:,:,k,2) = vvel(:,:,kk) * dxv(:,:) * dzt(:,:,k)
  end do
  
  return


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


contains

  real*4 function spherdist(lon1,lat1,lon2,lat2)
    implicit none
    real, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
    !
    ! --- ------------------------------------------------
    ! --- !omputes the distan!e between geo. pos.
    ! --- lon1,lat1 and lon2,lat2. 
    ! --- input is in degrees.
    !
    ! --- output is real*4 for better global consistancy,
    ! --- by truncating double precision roundoff errors.
    ! --- real*4 is not in f90, but is widely supported.
    !
    ! --- Based on m_spherdist.F90 from Geir Evanson.
    ! --- ------------------------------------------------
    !
    double precision, parameter :: invradian=0.017453292d0
    double precision, parameter ::    rearth=6371001.0d0  ! Radius of earth
    !
    double precision  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
    double precision  x1,y1,z1,x2,y2,z2                 ! Cartesian position
    double precision  dr                                ! Arc length
    !
    rlon1=lon1*invradian             !lon1 in rad
    rlat1=(90.d0-lat1)*invradian     !90-lat1 in rad 
    !
    rlon2=lon2*invradian             !lon2 in rad
    rlat2=(90.d0-lat2)*invradian     !90-lat2 in rad 
    !
    x1= sin(rlat1)*cos(rlon1)        !x,y,z of pos 1.
    y1= sin(rlat1)*sin(rlon1)
    z1= cos(rlat1) 
    !
    x2= sin(rlat2)*cos(rlon2)        !x,y,z of pos 2.
    y2= sin(rlat2)*sin(rlon2)
    z2= cos(rlat2) 
    !
    dr=acos(min(1.d0,x1*x2+y1*y2+z1*z2))  ! Arc length
    !
    spherdist=dr*rearth
    
  end function spherdist


end subroutine readfields
