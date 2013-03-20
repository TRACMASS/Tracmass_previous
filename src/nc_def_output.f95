subroutine nc_def_output

USE mod_param
USE mod_coord
USE mod_grid
USE mod_name
USE mod_netcdf
IMPLICIT none

INTEGER i,j,k,kk
      integer timedim, nfltdim
      integer t2dgrd(2) 

!=======================================================================
!  Create a new NetCDF file.
!=======================================================================

      rcode = nf90_create(TRIM(gridfile), nf90_clobber, ncgridid)
      if (rcode.ne.0) print *, trim(nf90_strerror(rcode))

!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------

      rcode = nf90_def_dim(ncgridid,'nfloat',L,nfltdim)
      if (rcode.ne.0) print *, trim(nf90_strerror(rcode))
      rcode = nf90_def_dim(ncgridid,'time',nf90_unlimited,timedim)
      if (rcode.ne.0) print *, trim(nf90_strerror(rcode))

      t2dgrd = (/ nfltdim, timedim /)

return
end subroutine nc_def_output

