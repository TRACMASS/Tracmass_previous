SUBROUTINE setupgrid
  
  USE netcdf
  USE mod_param
  USE mod_vel
  USE mod_coord
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
  INTEGER                                    :: ji,jj

  
  !!
  !! Vertical levels
  !!
  dz(41)=   3.d0
  dz(40)=   3.d0
  dz(39)=   3.d0
  dz(38)=   3.d0
  dz(37)=   3.d0
  dz(36)=   3.d0
  dz(35)=   3.d0
  dz(34)=   3.d0
  dz(33)=   3.d0
  dz(32)=   3.d0
  dz(31)=   3.d0
  dz(30)=   3.d0
  dz(29)=   3.d0
  dz(28)=   3.007080d0
  dz(27)=   3.063581d0
  dz(26)=   3.175872d0
  dz(25)=   3.342542d0
  dz(24)=   3.561495d0
  dz(23)=   3.829976d0
  dz(22)=   4.144610d0
  dz(21)=   4.501440d0
  dz(20)=   4.895979d0
  dz(19)=   5.323265d0
  dz(18)=   5.777925d0
  dz(17)=   6.254241d0
  dz(16)=   6.746222d0
  dz(15)=   7.247683d0
  dz(14)=   7.752317d0
  dz(13)=   8.253778d0
  dz(12)=   8.745760d0
  dz(11)=   9.222075d0
  dz(10)=   9.676735d0
  dz( 9)=   10.10402d0
  dz( 8)=   10.49856d0
  dz( 7)=   10.85539d0
  dz( 6)=   11.17002d0
  dz( 5)=   11.43851d0
  dz( 4)=   11.65746d0
  dz( 3)=   11.82413d0
  dz( 2)=   11.93642d0
  dz( 1)=   11.99292d0
  
  !!
  !! Grid box size
  !!
  dx = 1./15.d0
  dy = 1./30.d0
  DO ji=1,IMT
     DO jj=1,JMT
        dxdy(ji,jj)=dx * cst(jj) * dy * deg**2
     END DO
  END DO
  
  !!
  !! Min/max of density, temperature, and salinity
  !!
  rmin    =  -2.d0
  rmax    =  10.d0
  tmin    =  -2.d0
  tmax    =  25.d0 
  smin    =   0.d0
  smax    =  15.d0
  dr      = (rmax-rmin)/DBLE(MR-1)
  dtemp   = (tmax-tmin)/DBLE(MR-1)
  dsalt   = (smax-smin)/DBLE(MR-1)
  
  
END SUBROUTINE setupgrid
