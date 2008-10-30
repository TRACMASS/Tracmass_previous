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

  integer dimidx,dimidy,dimidz,dimidt !output ID index of dimension
  integer varid,varidx,varidy,varidz,varidt !output ID index of variable
  integer lenx,leny,lenz,lent,lenz2 !output Length of dimension
  integer p, x1, y1, z1, t1 !?

  ! = Variables for filename generation
  CHARACTER                                  :: dates(62)*17
  CHARACTER (len=200)                        :: dataprefix, dstamp
  INTEGER                                    :: intpart1 ,intpart2 ,subYr
  INTEGER                                    :: filePos ,fileJD ,subYrD
  INTEGER                                    :: yr1 ,mn1 ,dy1
  INTEGER                                    :: yr2 ,mn2 ,dy2
  

  ! = Variables used for getfield procedures
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  REAL*4, ALLOCATABLE, DIMENSION(:,:)        :: ssh
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: rhof

  REAL, ALLOCATABLE, DIMENSION(:,:)    :: gridLat ,gridLon  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v,e1t,e2u,e2t
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: dzu,dzv,dzt
  
  logical around
  
  alloCondGrid: if ( .not. allocated (e1v) ) then
     allocate ( gridLat(IMT,JMT) ,gridLon(IMT,JMT) )
     allocate ( e1v(IMT+2,JMT)    ,e1t(IMT+2,JMT) )
     allocate ( e2u(IMT+2,JMT)    ,e2t(IMT+2,JMT) )
     allocate ( dzu(IMT+2,JMT,KM) ,dzv(IMT+2,JMT,KM),dzt(IMT+2,JMT,KM) )
  end if alloCondGrid
  alloCondUVW: if(.not. allocated (ssh)) then
     allocate ( ssh(imt,jmt) )
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
  dstamp  = 'archv.0000_000_00_'
  
  write (dstamp( 7:10),'(i4.4)') currYear
  write (dstamp(12:14),'(i3.3)') currJDyr
  dataprefix  = trim(inDataDir) // '/data/'
  tpos        = intpart1
 
!  print *,ints, filePos ,fileJD ,currYear ,currMon ,currDay ,currJDyr,trim(dstamp)
 
  start1d  = [ 1]
  count1d  = [km]
  start2d  = [filePos , 1 ,subGridJmin ,subGridImin]
  count2d  = [      1 , 1 ,subGridJmax ,subGridImax]
  map3D    = [      4 , 3 ,          2 ,          1]  
  start3d  = [filePos , 1 ,          1 ,          1]
  count3d  = [      1 ,km ,        jmt ,        imt]
  map3D    = [      4 , 3 ,          2 ,          1]  
  ! === initialise ===
  initFieldcond: if(ints.eq.intstart) then
     ! call coordinat
     hs     = 0.
     uflux  = 0.
     vflux  = 0.
#ifdef tempsalt
     tem    = 0.
     sal    = 0.
     rho    = 0.
#endif

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================

     gridFile = trim(inDataDir)//'/uvel/'//trim(dstamp)//'3zu.nc'

     print *,' '
     print *,gridFile

     zw(0:km-1) = get1DfieldNC(trim(gridFile) ,'Depth')

     do k=1,km
        kk=km+1-k
        dz(kk)=zw(k)-zw(k-1) 
        ! print *,k,zw(k),kk,dz(kk)
     end do

     dz(1) = dz (2)

     gridLat  = get2DfieldNC(trim(gridFile) ,'Latitude')
     gridLon  = get2DfieldNC(trim(gridFile) ,'Longitude')
     
     print *,shape(gridLat)
     print *,gridLon
     !forall (i=1:imt ,j=1:jmt)
     !   spherdist(lon1,lat1,lon2,lat2)
     


     stop 666
     e2t  = get2DfieldNC(trim(gridFile) ,'e2t')
     e1v  = get2DfieldNC(trim(gridFile) ,'e1v')
     dxdy = e1t * e2t
   
     gridFile = trim(inDataDir)//'topo/mesh_zgr.nc'
     dzt  = get3DfieldNC(trim(gridFile) ,'e3t_ps')
     dzu  = get3DfieldNC(trim(gridFile) ,'e3u_ps')
     dzv  = get3DfieldNC(trim(gridFile) ,'e3v_ps')


     

     fieldFile = trim(inDataDir)//dataprefix
     rhof = get3DfieldNC(trim(fieldFile)//'_sigma.nc'  ,'sigma')
     do j=1,jmt
        do i=1,IMT
           do k=1,km
              kk=km+1-k
              if(rhof(i,j,k).ne.0.) kmt(i,j)=k
              if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
           enddo
           if(kmt(i,j).ne.0) then
              dztb(i,j,1)=dzt(i,j,kmt(i,j))
           else
              dztb(i,j,1)=0.
           endif
        enddo
     enddo
     
  endif initFieldcond
  
  fieldFile = trim(inDataDir)//'/data/'
  
  uvel = get3DfieldNC(trim(fieldFile)//'uvel/'//dataprefix//'3zu.nc' ,'u')
  uvel = get3DfieldNC(trim(fieldFile)//'_grid_V.nc' ,'vomecrty')
  rhof = get3DfieldNC(trim(fieldFile)//'_sigma.nc'  ,'sigma')
  ssh  = get2DfieldNC(trim(fieldFile)//'_SSH.nc'    ,'sossheig')

  hs(:,:,2) = 0.01*ssh
  hs(imt+1,:,2) =hs(1,:,2)
  hs(:,jmt+1,2) =hs(:,1,2)

  do k=1,km-1
     kk=km+1-k
     uflux(:,:,k,2) = uvel(:,:,kk) * e2u(:,:) * dzu(:,:,kk)
     vflux(:,:,k,2) = vvel(:,:,kk) * e1v(:,:) * dzv(:,:,kk)
     rho  (:,:,k,2) = rhof(:,:,kk)
  end do

  uflux(:,:,km,2) = uvel(:,:,1) * e2u * (dzu(:,:,1) & 
       + 0.5*(hs(:,:,2) + hs(2:imt+1,:,2)))
  vflux(:,:,km,2) = vvel(:,:,1) * e1v * (dzv(:,:,1) & 
       + 0.5*(hs(:,:,2) + hs(:,2:jmt+1,2)))
  rho(:,:,km,2)   = rhof(:,:,1)
  
  return


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
