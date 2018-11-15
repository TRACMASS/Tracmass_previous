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
  INTEGER                                    :: t ,i ,j ,k ,kk ,tpos
  INTEGER                                    :: fileInts
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
  CHARACTER (len=200)                        :: gridFile ,fieldFile
  CHARACTER (len=50)                         :: varName

  REAL*4, ALLOCATABLE, DIMENSION(:,:)        :: ssh
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)      :: rhof

  REAL, ALLOCATABLE, DIMENSION(:,:)    :: gridLat ,gridLon  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v,gridDX,e2u,gridDY
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: dzu,dzv,dzt
  
  logical around
  
  alloCondGrid: if ( .not. allocated (e1v) ) then
     allocate ( gridLat(IMT,JMT) ,gridLon(IMT,JMT) )
     allocate ( e1v(IMT+2,JMT)    ,gridDX(IMT+2,JMT) )
     allocate ( e2u(IMT+2,JMT)    ,gridDY(IMT+2,JMT) )
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
  if (currJDyr==366) currJDyr=365

  filePos  = mod((ints),fieldsPerFile)+1
  dstamp   = 'archv.2007_000_00_'

!  write (dstamp( 7:10),'(i4.4)') currYear
  write (dstamp(12:14),'(i3.3)') currJDyr
  dataprefix  = trim(inDataDir) // '/data/'
  tpos        = intpart1
  
  start1D  = [ 1]
  count1D  = [km]
  start2D  = [subGridImin ,subGridJmin , 1 ,1]
  count2D  = [        imt ,        jmt , 1 ,1]
  map2D    = [          1 ,          2 , 3 ,4]  
  start3D  = [subGridImin ,subGridJmin , 1 ,1]
  count3D  = [        imt ,        jmt ,km ,1]
  map3D    = [          1 ,          2 , 3 ,4]  
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

     zw(0:km-1)  = get1DfieldNC(trim(gridFile) ,'Depth')
     do k=1,km
        kk       = km+1-k
        dz(kk)   = zw(k)-zw(k-1) 
     end do
     dz(1)       = dz (2)

     gridLat     = get2DfieldNC(trim(gridFile) ,'Latitude')
     gridLon     = get2DfieldNC(trim(gridFile) ,'Longitude')

     do i=1,imt-1
        do j=1,jmt-1
           gridDY(i,j) = spherdist(gridLon(i,j) ,gridLat(i,j) &
                ,gridLon(i,j+1) ,gridLat(i,j+1))
            gridDX(i,j) = spherdist(gridLon(i,j) ,gridLat(i,j) &
                 ,gridLon(i+1,j) ,gridLat(i+1,j))
         end do
      end do




      gridDX(:,jmt)   = gridDX(:,jmt-1)
      gridDY(:,jmt)   = gridDY(:,jmt-1)
      gridDX(imt,:)   = gridDX(imt-1,:)
      gridDY(imt-1:imt,:)   = gridDY(imt-3:imt-2,:)
      gridDY(imt,jmt) = gridDY(imt-1,jmt-1)
      
      e1v        = gridDX
      e2u        = gridDY
      dxdy       = gridDX*gridDY

      uvel     = get3DfieldNC(trim(gridFile) ,'u')

      do j=1,jmt
         do i=1,IMT
            do k=1,km
               kk=km+1-k
               if(uvel(i,j,k)<10000.) kmt(i,j)=k
               !if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
            enddo
            if(kmt(i,j).ne.0) then
               dztb(i,j,1)=dzt(i,j,kmt(i,j))
            else
               dztb(i,j,1)=0.
            endif
         enddo
      enddo
      
  endif initFieldcond
  

  fieldFile = trim(inDataDir)//'/uvel/'//trim(dstamp)//'3zu.nc'
  uvel = get3DfieldNC(trim(fieldFile) ,'u')
  where (uvel>10000) uvel=0
  fieldFile = trim(inDataDir)//'/vvel/'//trim(dstamp)//'3zv.nc'
  vvel = get3DfieldNC(trim(fieldFile) ,'v')
  where (vvel>10000) vvel=0
 
  hs(:,:,2) = 0.01*ssh
  hs(imt+1,:,2) =hs(1,:,2)
  hs(:,jmt+1,2) =hs(:,1,2)

  do k=1,km
     kk=km+1-k
     uflux(:,:,k,2) = uvel(:,:,kk) * e2u(:,:) * dz(k) !dzu(:,:,kk)
     vflux(:,:,k,2) = vvel(:,:,kk) * e1v(:,:) * dz(k) !dzv(:,:,kk)
   !  rho  (:,:,k,2) = rhof(:,:,kk)
  end do

  uflux(:,:,km,2) = uvel(:,:,1) * e2u * (dzu(:,:,1) & 
       + 0.5*(hs(:,:,2) + hs(2:imt+1,:,2)))
  vflux(:,:,km,2) = vvel(:,:,1) * e1v * (dzv(:,:,1) & 
       + 0.5*(hs(:,:,2) + hs(:,2:jmt+1,2)))
  !rho(:,:,km,2)   = rhof(:,:,1)
  
!  print *,uvel(100,100,1),vvel(100,100,1),gridDX(100,100),gridDY(100,100),dz(33)
!  print *,uvel(100,100,1)*gridDX(
!  stop



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
