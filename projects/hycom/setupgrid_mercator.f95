SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  USE mod_getfile

  IMPLICIT none
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of the GCM. Run once
  ! before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Area of horizontal cell-walls.
  !  dz   - Height of k-cells in 1 dim. |\
  !  dzt  - Height of k-cells i 3 dim.  |- Only one is needed
  !  kmt  - Number of k-cells from surface to seafloor.
  !
  ! The following might be needed to calculate
  ! dxdy, uflux, and vflux
  !
  !  dzu - Height of each u-gridcell.
  !  dzv - Height of each v-gridcell.
  !  dxu -
  !  dyu -
  ! -------------------------------------------------------------



  ! === Init local variables for the subroutine ===
  INTEGER                                    :: i ,j ,k ,kk
  REAL,          ALLOCATABLE, DIMENSION(:)   :: zw
  REAL,       ALLOCATABLE, DIMENSION(:,:)    :: gridLat ,gridLon  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)    :: e1v,gridDX,e2u,gridDY
  CHARACTER (len=200)                        :: gridfile

  

  allocate ( zw(km) )

! alloCondGrid: if ( .not. allocated (e1v) ) then
!     allocate ( gridLat(IMT,JMT)    ,gridLon(IMT,JMT) )
!     allocate ( e1v(IMT+2,JMT)      ,gridDX(IMT+2,JMT) )
!     allocate ( e2u(IMT+2,JMT)      ,gridDY(IMT+2,JMT) )
  allocate ( dzu(IMT+2,JMT,KM,2) ,dzv(IMT+2,JMT,KM,2),dzt(IMT+2,JMT,KM,2) )
  allocate ( gridDY(IMT+2,JMT)      ,gridDX(IMT+2,JMT) )
!  end if alloCondGrid
  

! === Template for setting up grids. Move the code from readfile.f95

! ===


!!$  CHARACTER (len=200)                        :: gridFileXY, gridFileZ
!!$  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: kmask
!!$
!!$  alloCondGrid: if ( .not. allocated (kmask) ) then
!!$     allocate ( kmask(IMT+2,JMT,KM) )
!!$  end if alloCondGrid
!!$  
!!$  start1d  = [  1]
!!$  count1d  = [ km]
!!$  !Order is     t    k            i            j
!!$  start2d  = [  1 ,  1 ,subGridImin ,subGridJmin]
!!$  count2d  = [  1 ,  1 ,subGridImax ,subGridJmax]
!!$  map2d    = [  4 ,  3 ,          1 ,          2]  
!!$  start3d  = [  1 ,  1 ,subGridImin ,subGridJmin]
!!$  count3d  = [  1 , km ,subGridImax ,subGridJmax]
!!$  map3d    = [  4 ,  3 ,          2 ,          1]  
!!$  
!!$  gridFileXY = trim(inDataDir)//'grid_cell_xy.nc'
!!$  gridFileZ  = trim(inDataDir)//'grid_cell_z.nc'
!!$  
!!$  dz   = get1DfieldNC(trim(gridFileZ)  ,'dz')  / 100.
!!$  dxv  = get2DfieldNC(trim(gridFileXY) ,'DXU') / 100.
!!$  dyu  = get2DfieldNC(trim(gridFileXY) ,'DYU') / 100.
!!$  dxdy = dxv * dyu
!!$
!!$  dzt = 0
!!$  kmask  = get3DfieldNC(trim(gridFileZ) ,'SALT')
!!$  do j=1,jmt
!!$     do i=1,imt
!!$        do k=1,km
!!$           kk=km+1-k
!!$           if(kmask(i,j,k) .le. 1000.) then
!!$              kmt(i,j)=k
!!$              dzt(i,j,k) = dz(kk)
!!$           end if
!!$        enddo
!!$     enddo
!!$  enddo



  ! ======================================================
  !    ===  Set up the grid ===
  ! ======================================================
  
  gridFile = trim(inDataDir)//trim('archv.2014_356_00_3zu.nc')
  ncTpos = 1

  zw(0:km-1)  = get1DfieldNC(trim(gridFile) ,'Depth')
  do k=1,km
     kk       = km+1-k
     dz(kk)   = zw(k)-zw(k-1) 
  end do
  dz(1)       = dz (2)

  map2d    = [3, 4, 1, 1]
  map3d    = [2, 3, 4, 1]
  gridLat     = get2DfieldNC(trim(gridFile) ,'Latitude')
  gridLon     = get2DfieldNC(trim(gridFile) ,'Longitude')
  
  do i=1,imt-1
     do j=1,jmt-1
        gridDY(i,j) = spherdist(gridLon(i,j),   gridLat(i,j), &
                                gridLon(i,j+1), gridLat(i,j+1))
        gridDX(i,j) = spherdist(gridLon(i,j),   gridLat(i,j), &
                                gridLon(i+1,j), gridLat(i+1,j))
     end do
  end do
  
  gridDX(:,jmt)   = gridDX(:,jmt-1)
  gridDY(:,jmt)   = gridDY(:,jmt-1)
  gridDX(imt,:)   = gridDX(imt-1,:)
  gridDY(imt-1:imt,:)   = gridDY(imt-3:imt-2,:)
  gridDY(imt,jmt) = gridDY(imt-1,jmt-1)
  
  dxv        = gridDX
  dyu        = gridDY
  dxdy       = gridDX * gridDY


  map2d    = [3, 4, 1, 1]
  map3d    = [2, 3, 4, 1]
  ncTpos = 1
  uvel     = get3DfieldNC(trim(gridFile) ,'u')  
  do j=1,jmt
     do i=1,IMT
        do k=1,km
           kk=km+1-k
           if(uvel(i,j,k)<10000.) kmt(i,j)=k
           !if(k.ne.kmt(i,j)) dz(kk)=dzt(i,j,k)
        enddo
        !if(kmt(i,j).ne.0) then
        !   dztb(i,j,1)=dzt(i,j,kmt(i,j))
        !else
        !   dztb(i,j,1)=0.
        !endif
     enddo
  enddo

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

  
end SUBROUTINE setupgrid
