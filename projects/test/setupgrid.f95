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
   
   !!-------------------------------------------------------------
    
   !!
   !! Set name of testcase
   !!
   
   testcase = trim(GCMsource)
   
   ! -------------------------------------------------------------
   
   if (trim(testcase) == 'circular-stationary'   .or. &
   &   trim(testcase) == 'circular-evolving'     .or. &
   &   trim(testcase) == 'time-oscillation1'     .or. &
   &   trim(testcase) == 'time-oscillation2'     .or. &
   &   trim(testcase) == 'time-oscillation3'     .or. &
   &   trim(testcase) == 'nicoletta-fabboni'          ) then
      
      dx = 250. 
      dy = 250.
      dz(1:km) = 1.
      
      dxdy(1:imt,1:jmt) = dx*dy
      
      mask(1:imt,1:jmt) = 1
      kmt(1:imt,1:jmt) = km
      
      
   else if (trim(testcase) == 'bower-gulf') then
      
      !! === Set up a grid with x = 0 to the left and y = 0 in the middle
      
      dx   = (xmax - xmin) / (imt - 1)
      dxv(1:imt,1:jmt) = dx
      
      dy   = (ymax - ymin) / (jmt - 1)
      dyu(1:imt,1:jmt) = dy
   
      dz(1:km) = (zmax - zmin) / (km - 1)
   
      dxdy(1:imt,1:jmt) = dx*dy
      
      !! All cells active
      mask(1:imt,1:jmt) = 1 
      kmt(1:imt,1:jmt) = km
   
   else
      
      print*, trim(testcase)//' is not a valid testcase '
      print*,' Stopping in setupgrid '
      stop
      
   end if
    
   return

end subroutine setupgrid
