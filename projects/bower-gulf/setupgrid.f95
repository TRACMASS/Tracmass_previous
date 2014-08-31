subroutine setupgrid
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
  
   use mod_param
   use mod_vel
   use mod_coord
   use mod_time
   use mod_grid
   use mod_name
   use mod_vel
  
   implicit none
   
   !! Domain length, width and depth
   real(4)                                  ::  xmin, xmax, ymin, ymax, zmin, zmax
   
   
   !!
   !! Read in data from namelist and set up grid
   !!                       
   
   namelist /GRIDSIZE/  xmin, xmax, ymin, ymax, zmin, zmax
   
   open (unit=11,file='projects/bower-gulf/namelist',status='old',delim='apostrophe')
   read (11,nml=GRIDSIZE) 
   close (11)
   
   !! x grid
   dx   = (xmax - xmin) / (imt - 1)
   dxv(1:imt,1:jmt) = dx
   
   !! y grid
   dy   = (ymax - ymin) / (jmt - 1)
   dyu(1:imt,1:jmt) = dy
   
   !! depth of each grid box
   dz(1:km) = (zmax - zmin) / (km - 1)
   
   !! Area of horizontal grid faces
   dxdy(1:imt,1:jmt) = dx*dy
  
   !! Mask
   mask(1:imt,1:jmt) = 1 !! All cells active
   
   !! Number of depth levels
   kmt(1:imt,1:jmt) = km
    
   return

end subroutine setupgrid
