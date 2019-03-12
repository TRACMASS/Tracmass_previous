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
  real, allocatable, dimension(:) :: depth
  logical :: lwrite_nc, lread_nc

lwrite_nc = .true.
lread_nc  = .true.

allocate( lon(imt,jmt), lat(imt,jmt) )
allocate( depth(km) )

kmt=KM ! flat bottom

!dxdeg=dx*deg
!dydeg=dy*deg

! Nicoletta Fabboni velocities, which have analytical solutions
dxv (:,:) = 250. 
dyu (:,:) = dxv(:,:)
dxdy(:,:) = dxv(:,:) * dyu(:,:)
dz  (:)   = 10.
dzt(:,:,:,:) = 10.

mask(:,:) = 1

do j = 1, jmt
   do i = 1, imt
      lon(i,j) = i * dxv(i,j) / 111000.
      lat(i,j) = j * dyu(i,j) / 111000.
   end do 
end do

depth(1) = dz(1)/2.
do k = 2, km
   depth(k) = SUM(dz(1:k-1)) + dz(k)/2.
end do

if (lwrite_nc) then 
   !! Write grid info to netCDF file
   call check( nf90_create('mesh.nc', nf90_hdf5, ncid) )
   call check( nf90_def_dim(ncid, 'x', imt, idx) )
   call check( nf90_def_dim(ncid, 'y', jmt, idy) )
   call check( nf90_def_dim(ncid, 'z', km, idz) )
   dim1d = (/ idz /)
   dim2d = (/ idx, idy /)
   dim3d = (/ idx, idy, idz /)
   call check( nf90_def_var(ncid, 'lon', nf90_real, dim2d, idlon) )
   call check( nf90_def_var(ncid, 'lat', nf90_real, dim2d, idlat) )
   call check( nf90_def_var(ncid, 'depth', nf90_real, dim1d, iddep) )
   call check( nf90_def_var(ncid, 'dx', nf90_real, dim2d, iddx) )
   call check( nf90_def_var(ncid, 'dy', nf90_real, dim2d, iddy) )
   call check( nf90_def_var(ncid, 'dz', nf90_real, dim3d, iddz) )
   call check( nf90_enddef(ncid) )
   
   call check( nf90_put_var(ncid, idlon, lon) )
   call check( nf90_put_var(ncid, idlat, lat) )
   call check( nf90_put_var(ncid, iddep, depth) )
   call check( nf90_put_var(ncid, iddx, dxv(1:jmt,1:imt)) )
   call check( nf90_put_var(ncid, iddy, dyu(1:jmt,1:imt)) )
   call check( nf90_put_var(ncid, iddz, dzt(1:imt,1:jmt,1:km,2)) )
   call check( nf90_close(ncid) )
end if

if (lread_nc) then 
  !! Read grid info from netCDF file
   ncTpos = 1
   map2d = [3, 4, 1, 1]
   map3d = [2, 3, 4, 1]
   lon = get2DfieldNC('mesh.nc', 'lon')
   lat = get2DfieldNC('mesh.nc', 'lat')
   depth = get1DfieldNC('mesh.nc', 'depth')
   dxv(1:imt,1:jmt) = get2DfieldNC('mesh.nc', 'dx')
   dyu(1:imt,1:jmt) = get2DfieldNC('mesh.nc', 'dy')
   dzt(1:imt,1:jmt,1:km,2) = get3DfieldNC('mesh.nc', 'dz')
   dzt(:,:,:,1) = dzt(:,:,:,2)
end if


! ===

contains
   SUBROUTINE check(ierr2)
   
      INTEGER :: ierr2
   
      IF( ierr2 /= 0) THEN
         PRINT*, NF90_STRERROR(ierr2)
         STOP
      END IF
   
   END SUBROUTINE 




end SUBROUTINE setupgrid
