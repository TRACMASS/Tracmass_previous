SUBROUTINE setupgrid
  
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  
  USE mod_getfile
  USE netcdf
  
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

  integer :: idx, idy, idz, idlon, idlat, iddep, iddx, iddy, iddz
  integer, dimension(1) :: dim1d
  integer, dimension(2) :: dim2d
  integer, dimension(3) :: dim3d
  real, allocatable, dimension(:,:) :: lon, lat
  REAL*4:: rlatt

!  real, allocatable, dimension(:) :: depth
  logical :: lwrite_nc, lread_nc

allocate( lon(imt,jmt), lat(imt,jmt) )
!allocate( depth(km) )

kmt=KM ! flat bottom

!dxdeg=dx*deg
!dydeg=dy*deg

dxv = 250. 
dyu = 250.
!dxdy(1:IMT,1:JMT) = dxv(1:IMT,1:JMT) * dyu(1:IMT,1:JMT)
dxdy = 250. * 250.
dz  = 10.
dzt = 10.

mask(:,:) = 1

do j = 1, jmt
   do i = 1, imt
      lon(i,j) = i * dxv(i,j) / 111000.
      lat(i,j) = j * dyu(i,j) / 111000.
   end do 
end do

!depth(1) = dz(1)/2.
!do k = 2, km
!   depth(k) = SUM(dz(1:k-1)) + dz(k)/2.
!end do



! To mimic a spherical grid
 !!
  !! Grid box sizes
  !!
  dx = 1. ! one degree resolution
  dy = dx  ! Horizontal resolution in degrees [deg]
  print *,'dy=',dy
  dxdeg = dx*deg        
  dydeg = dy*deg        ! Horizontal resolution in radians [rad]
  
  dyu=dydeg
  
  print *,'dxdeg=', dxdeg,dydeg
  
  !! Cosine at each v-point (C-grid) 
  DO j=0,JMT
     phi(j) = 60. + dy * FLOAT(j)    ! Latitude for each v-point (C-grid)
     csu(j) = DCOS (phi(j) * radian) 
     if(j/=0) dxv(:,j)=dxdeg*DCOS (phi(j) * radian)
  END DO
  
  !! Cosine at each u-point (C-grid)
  DO j=1,JMT
     rlatt = 0.5 * ( phi(j) + phi(j-1) )   ! Latitude at each u-point (C-grid)
     cst(j) = DCOS ( rlatt * radian )
  END DO
  
  !! Total horizontal area of grid box
     DO j=1,JMT
  DO i=1,IMT
        dxdy(i,j) = dxdeg * cst(j) * dydeg
     END DO
!     print *,j,dxdy(1,j),phi(j),csu(j),dxv(1,j)
  END DO




!stop 3956

end SUBROUTINE setupgrid
