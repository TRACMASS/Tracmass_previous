SUBROUTINE readfields
!!---------------------------------------------------------------------------------
!!
!!
!!   Set up mass fluxes and tracer fields from idealised flows 
!!
!!   NOTE: These test cases do not need any netCDF modules or any input fields. 
!!
!!
!!   The following examples exist: 
!!
!!      1) circular-stationary - Circular stationary field
!!      2) circular-evolving   - Time-evolving circular field
!!      3) time-oscillation1   - 3 different fields with oscillations  
!!         time-oscillation2
!!         time-oscillation3 
!!      4) Nicoletta-Fabboni   - Oscillatory field with analytical solutions
!!      5) bower-gulf          - Idealised gulf stream model from Bower (JPO,1991) 
!!   
!!
!!   History: 
!!      
!!      Sep 2014, JK: Examples 1-4 by K. Doos and example 5 merged
!!
!!   
!!   Joakim Kjellsson
!!   British Antarctic Survey
!!   August 2014
!!
!!
!!
!!---------------------------------------------------------------------------------

  use mod_param
  use mod_vel
  use mod_coord
  use mod_time
  use mod_grid
  use mod_name
  use mod_vel
  use mod_traj
  use mod_seed

#ifdef tempsalt
  use mod_stat
#endif
  
  implicit none
  
   !!------------------------------------------------------------------------------
   
   !! Loop/dummy variables
   integer                                  ::  ji, jj, jk, jn, jim, jip, jjm, jjp, &
   &                                            jkm, jkp, jnm, jnp
   real                                     ::  x, y, z, t, alpha, beta, yc
   real(4), allocatable, dimension(:,:,:)   ::  dPsidx, dPsidy
   
   !! Grid 
   real(4)                                  ::  xmin, xmax, ymin, ymax, zmin, zmax
   character(len=4), dimension(4)           ::  hour
   character(len=2), dimension(12)          ::  month
   character(len=2), dimension(31)          ::  date
   
   !! For 
   real(8), save                            ::  cox, coy, om, co, uwe, dl
   real(8), save                            ::  ug, u0, gamma, gammag, fcor
   
   !! For gulf stream function
   real(4)                                  ::  Lx, kx, Ay, Hz, lambda, cx, sc, psi0
   real(4), allocatable, dimension(:,:,:,:) ::  psi
   !!                                           km/d -> m/s
   real(4), parameter                       ::  kmd = 1000.0 / (60.0*60.0*24.0)
#ifdef tempsalt
   real(4), parameter                       ::  thigh = 25.0 !![C]
   real(4), parameter                       ::  shigh = 35.0 !![PSU]
   real(4), parameter                       ::  tlow = 10.0  !![C]
   real(4), parameter                       ::  slow = 33.0  !![PSU]
#endif
      
   
   !!------------------------------------------------------------------------------
   
   namelist /GRIDSIZE/  xmin, xmax, ymin, ymax, zmin, zmax
   
   !!------------------------------------------------------------------------------
   
   data month /'01','02','03','04','05','06','07','08','09','10','11','12'/
   data date  /'01','02','03','04','05','06','07','08','09','10',&
               '11','12','13','14','15','16','17','18','19','20',&
               '21','22','23','24','25','26','27','28','29','30','31'/
   data hour  /'0000','0600','1200','1800'/
   
   !!------------------------------------------------------------------------------
   
   !!
   !! Put time step 2 at step 1 so we can put the new stuff at step 2
   !!
   
   call datasetswap
   
   !!------------------------------------------------------------------------------
   !! === Grid is determined from a namelist and is read here and in setupgrid.f95
   !!------------------------------------------------------------------------------
   
   open (unit=11,file='projects/bower-gulf/namelist',status='old',delim='apostrophe')
   read (11,nml=GRIDSIZE) 
   close (11)
   
   !! We assume a horizontal grid that is 800 x 500 km 
   !! and 1000 m depth                                 
   
   !! x grid
   dx   = (xmax - xmin) / (imt - 1)
   
   !! y grid
   dy   = (ymax - ymin) / (jmt - 1)
   
   
   !!------------------------------------------------------------------------------
   
   
   !!------------------------------------------------------------------------------
   
   !! 
   !! Determine the time in seconds since start
   !! (We are using a very simple calendar with 30-day months)
   !!
   
   if (ints == intstart) then !! if first time step
      
      t = 0.
      
      currHour = startHour
      currDay  = startDay
      currMon  = startMon
      currYear = startYear
      
   else 
      
      t = t + nff * ngcm * 60.0 * 60.0 !![s]
      
      currHour = currHour + nff * ngcm
      
      if (currHour >= 24) then
         
         currHour = currHour - 24
         currDay  = currDay + 1
         
      elseif (currHour < 0) then 
         
         currHour = currHour + 24
         currDay  = currDay - 1
         
      endif
      
      if (currDay > 30) then 
         
         currDay = currDay - 30
         currMon = currMon + 1
         
         if (currMon == 13) then
            
            currMon  = 1
            currYear = currYear + 1
         
         end if
            
      else if (currDay <= 0) then
         
         currDay  = currDay + 30   
         currMon  = currMon - 1
            
         if(currMon == 0) then
               
            currMon  = 12
            currYear = currYear - 1
               
         end if
            
      end if
   
   end if   
   
   100 format(' Current date and time: ',i4,'-',i2,'-',i2,' ',i2,':00')
   !print 100,currYear,currMon,currDay,currHour 
   
   
   !!------------------------------------------------------------------------------
   
   !!
   !! Set parameters specific for each test case
   !!
   
   if ( testcase == 'circular-stationary' .or. &
      & testcase == 'circular-evolving'   .or. &
      & testcase == 'time-oscillation3'        ) then
      
      t = dble(ints) / dble(5) + 1.
      
      cox = 0.5d0 + 0.5d0 * dcos(t)
      coy = 0.5d0 + 0.5d0 * dcos(t+pi)
      
   end if
   
   if ( testcase == 'nicoletta-fabboni') then
      
      ug     = 0.04  
      u0     = 0.5 
      fcor   = 2. * 2.*PI / (24.*3600.) * cos(45.*PI/180.)
      
      gamma  = 1./(2.89 *24. * 3600.)
      gammag = 1./(28.9 *24. * 3600.)
      t = dble(ints-intstart) * dble(ngcm * 3600)

   end if
   
   if (testcase == 'circular-evolving') then
      
      !! === Parameters for circular time evolving field
      uwe = -0.4d0
      dl  = dble(ints) * 0.01d0 * pi
   
   else if (testcase == 'bower-gulf') then
      
      !! === Parameters for Amy Bowers gulf stream ===
      
      Lx = 400.0      !! Wave length in [km]
      kx = 2.0*PI/Lx  !! Wave number [m-1]
      
      Ay = 50.0       !! wave amplitude in y direction [km]
      
      Hz = 100.0      !! e-folding depth scale, psi decreases with depth [m]
      
      lambda = 40.0   !! width of jet in y [km]
      
      cx = 10.0       !! phase speed of wave in x [km/day] 
                      !! 10 km/day approx. 0.1 m/s
      
      sc = 200.0      !! speed at jet center [km/day]
                      !! 200 km/day approx. 2.3 m/s
      
      !! Convert to SI units
      Lx     = Lx * 1000.0 !![m]
      Ay     = Ay * 1000.0 !![m]
      lambda = lambda * 1000.0 !![m]
      cx     = cx * kmd  !![m/s]
      sc     = sc * kmd !![m/s]
      
      psi0   = sc * lambda !! Max amplitude of stream function [m2/s]
   
   
   end if
   
   
   !!------------------------------------------------------------------------------
   
   if (testcase == 'bower-gulf') then
         
      !!
      !! Calculate stream function for this time step
      !!
      
      allocate ( psi(imt,jmt,km,1) )
      
      do jk=1,km
         
         do jj=1,jmt
            
            do ji=1,imt
               
               x = float(ji) * dx     + xmin
               y = float(jj) * dy     + ymin
               z = float(jk) * dz(jk) + zmin
               
               yc = Ay * sin( kx * x )
               alpha = atan( Ay * kx * cos(kx * x) )
               beta  = (y - yc) / (lambda/cos(alpha))
               
               !! Equation (2) from Bower (JPO, 1991)
               psi(ji,jj,jk,1) = psi0 * exp(-z/Hz) * (1.0 - tanh(beta)) + cx*y 
            
            end do
            
         end do
         
      end do
      
      !!
      !!  Calculate u and v from 
      !!  u = - dPsi/dy
      !!  v =   dPsi/dx
      !!
      
      allocate ( dPsidx(imt,jmt,km), dPsidy(imt,jmt,km) )
      
      dPsidy(1:imt,2:jmt,1:km) = (psi(1:imt,2:jmt,1:km,1) - psi(1:imt,1:jmt-1,1:km,1)) / dy
      dPsidy(1:imt,1,1:km)     = dPsidy(1:imt,2,1:km)
      
      dPsidx(2:imt,1:jmt,1:km) = (psi(2:imt,1:jmt,1:km,1) - psi(1:imt-1,1:jmt,1:km,1)) / dx
      dPsidx(1,    1:jmt,1:km) = (psi(1,    1:jmt,1:km,1) - psi(imt,    1:jmt,1:km,1)) / dx
   
      
   end if
   
   !!------------------------------------------------------------------------------
   
   !!
   !!  Calculate uflux and vflux
   !!
   
   do jk=1,km
      
      if (testcase == 'circular-stationary') then
         
         !!----------------------------------
         !! === Circular stationary field ===
         !!----------------------------------
         do jj = 1,jmt
            
            do ji = 1,imt
               
               uflux(ji,jj,jk,nsp) = -cox * dy * deg * dz(k) * dble(jj-1-jmt/2) / dble(jmt-1) 
               vflux(ji,jj,jk,nsp) =  cox * dx * deg * dz(k) * dble(ji-1-imt/2) / dble(imt-1)
            
            end do
            
         end do
         
      else if (testcase == 'circular-evolving') then
          
         !!----------------------------------------
         !! === Time evolving oscillating field ===
         !!----------------------------------------
         do jj=1,jmt
            
            do ji=1,imt
         
               uflux(ji,jj,jk,nsp) = dy * deg * dz(k) * cox * &
               &                     ( dcos( PI * dble(ji-1-imt/2)/dble(imt-1)   + dl) *  &
               &                       dsin(-PI * dble(jj-1-jmt/2)/dble(jmt-1) ) + uwe    &
               &                       + dcos(t) )
               
               vflux(ji,jj,jk,nsp) = dx * deg * dz(k) * coy * &
               &                     ( dsin( PI * dble(i-1-imt/2)/dble(imt-1)    + dl) *  & 
               &                       dcos( PI * dble(j-1-jmt/2)/dble(jmt-1) )           &
               &                       + dsin(t) )
               
            end do
         
         end do
      
         
      else if (testcase == 'time-oscillation1') then   
         
         !!---------------------------------------------------------
         !! === Spatially constant field that oscillates in time ===
         !!---------------------------------------------------------
         do jj=1,jmt
            
            do ji=1,imt
               
               uflux(ji,jj,jk,nsp) = dy * deg * dz(k) * (dcos(t)-0.01  + dcos(t/pi))
               vflux(ji,jj,jk,nsp) = dx * deg * dz(k) * (dsin(t)-0.001 + dsin(t/pi))
               
            end do
            
         end do
         
         
      else if (testcase == 'time-oscillation2') then
         
         !!---------------------------------------------------------
         !! === Spatially constant field that oscillates in time ===
         !!---------------------------------------------------------
         do jj=1,jmt
            
            do ji=1,imt
               
               uflux(ji,jj,jk,nsp) = dy * deg * dz(k) * (-0.01  + dcos(t/pi))
               vflux(ji,jj,jk,nsp) = dx * deg * dz(k) * (-0.001 + dsin(t/pi))
               
            end do
            
         end do
      
         
      else if (testcase == 'time-oscillation3') then
         
         !!---------------------------------------------------------
         !! === Spatially constant field that oscillates in time ===
         !!---------------------------------------------------------
         do jj=1,jmt
            
            do ji=1,imt
               
               uflux(ji,jj,jk,nsp) = dy * deg * dz(k) * cox *                &
               &                     ( -0.05 + t/5000.    + dcos(t) +        &
               &                        0.5  * dsin(t*2.) + 2. * dsin(t*10.) )
               
               vflux(ji,jj,jk,nsp) = dx * deg * dz(k) * coy *                &
               &                     ( -0.025 + dsin(t) + 0.5 * dcos(t*2.) + &
               &                        0.5 * dcos(t/2)                      )
               
            end do
            
         end do
         
         
      else if (testcase == 'nicoletta-fabboni') then
         
         !!----------------------------------------------------------------------
         !! === Nicoletta Fabboni velocities, which have analytical solutions ===
         !!----------------------------------------------------------------------
         uflux(1:imt,1:jmt,jk,nsp) = dy * dz(k) * ( ug * dexp(-gammag*t) + &
         &                           (u0-ug) * dexp(-gamma*t) * cos(fcor*t) )
         
         vflux(1:imt,1:jmt,jk,nsp) = dx * dz(k) * (-(u0-ug) * dexp(-gamma*t) * &
         &                           sin(fcor*t) )
      
         
      else if (testcase == 'bower-gulf') then
         
         !!-------------------------------------------------
         !! === Idealised gulf stream model by Amy Bower ===
         !!-------------------------------------------------
         
         !! Interpolate dPsidy to u points and multiply by dy and dz
         uflux(1:imt-1,1:jmt,jk,nsp) = -0.5 * ( dPsidy(2:imt  ,1:jmt,jk) + & 
         &                                      dPsidy(1:imt-1,1:jmt,jk) ) * &
         &                                      dyu   (1:imt-1,1:jmt) * dz(jk)
         
         !! Special case for cyclic point in x
         uflux(imt,    1:jmt,jk,nsp) = -0.5 * ( dPsidy(imt,    1:jmt,jk) + & 
         &                                      dPsidy(1,      1:jmt,jk) ) * &
         &                                      dyu   (1,      1:jmt) * dz(jk)
         
         !! Interpolate dPsidx to v points and multiply by dx and dz
         vflux(1:imt,1:jmt-1,jk,nsp) =  0.5 * ( dPsidx(1:imt,2:jmt,  jk) + & 
         &                                      dPsidx(1:imt,1:jmt-1,jk) ) * &
         &                                      dxv   (1:imt,1:jmt-1) * dz(jk)
         
         !! Set zero flux at y boundaries
         vflux(1:imt,  0,jk,nsp) = 0.0
         vflux(1:imt,jmt,jk,nsp) = 0.0
         
         
      end if 
   
   
   end do
   
   
   
   !!------------------------------------------------------------------------------
   
   
#ifdef tempsalt 

   thigh = 30.
   tlow  = 10.
   shigh = 35.
   slow  = 33.
   
   do jk=1,km
      
      do jj=1,jmt
         
         do ji=1,imt
            
            x = float(ji) * dx     + xmin
            y = float(jj) * dy     + ymin
            z = float(jk) * dz(jk) + zmin
            
            if (testcase == 'circular-stationary'   .or. &
                testcase == 'circular-evolving'     .or. &
                testcase == 'time-oscillation1'     .or. & 
                testcase == 'time-oscillation2'     .or. & 
                testcase == 'time-oscillation3'     .or. & 
                testcase == 'nicoletta-fabboni'          ) then
               
               !! 
               !! Simplified distribution of T,S and rho
               !!
               tem(ji,jj,jk,nsp) = 20. * float(k) / float(km)
               sal(ji,jj,jk,nsp) = 30.
               rho(ji,jj,jk,nsp) = (28. - 20.) * float(km-k) / float(km) + 20. 
            
            else if (testcase == 'bower-gulf') then
               
               !!
               !! Set up a simple distribution of temperature and salinity 
               !! where temperature and salt in low at high y and high and low y
               !! Meridional variations are set by a arctan function
               !! Uniform with depth
               !! 
               yc = Ay * sin( kx * (x-cx*t) )
               
               alpha = atan( Ay * kx * cos(kx * (x-cx*t)) )
               beta  = (y - yc) / (lambda/cos(alpha))
               
               toff = 0.5 * (tmax + tmin)
               tscale = tmax-tmin
               
               tem(ji,jj,jk,nsp) = (-1.0) * atan(beta) / PI * (thigh-tlow) + 0.5*(thigh+tlow)
               sal(ji,jj,jk,nsp) = (-1.0) * atan(beta) / PI * (shigh-slow) + 0.5*(shigh+slow)
            
            end if
            
         end do
      
      end do
   
   end do
   
   if (testcase == 'bower-gulf') then
      
      do jj=1,jmt
         
         do ji=1,imt
            
            call statvd(tem(ji,jj,:,nsp),sal(ji,jj,:,nsp),rho(ji,jj,:,nsp),km)
         
         end do
         
      end do
   
   end if
   
#endif
   
   !!------------------------------------------------------------------------------
   !!------------------------------------------------------------------------------  
   
   !!
   !! Print the fields to binary files for testing
   !!
   
   !open(unit=11,file='psi.bin',form='unformatted')
   !write(11) psi,dPsidy,dPsidx
   !close(11)
   !stop
   
   !!------------------------------------------------------------------------------  
   
   if (testcase == 'bower-gulf') then
      deallocate ( psi, dPsidx, dPsidy )
   end if
   
   return
   
   !!------------------------------------------------------------------------------
   
end subroutine readfields

