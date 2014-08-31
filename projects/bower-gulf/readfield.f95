SUBROUTINE readfields
!!---------------------------------------------------------------------------------
!!
!!
!!   Calculate stream function and velocities from the idealised 
!!   Gulf Stream model by Amy Bower. 
!!   
!!   The function is based on a cosine wave in x and arctan in y that produces
!!   a meandering jet. 
!!   The jet is calulcated in a frame of reference that moves with a given 
!!   phase speed, i.e. the function is time-independent.
!!
!!   Reference: Bower, A.S., 1991, 
!!              A Simple Kinematic Mechanism for Mixing
!!              Fluid Particles Across a Meandering Jet. 
!!              J. Phys. Oceanogr., 21, 173-180
!!
!!   This test case does not need any netCDF modules or any
!!   input fields. 
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
   
   !! Dummy variables
   integer                                  :: ji, jj, jk, jn, jim, jip, jjm, jjp, &
   &                                           jkm, jkp, jnm, jnp
   real                                     :: x, y, z, t, alpha, beta, yc
   real(4), allocatable, dimension(:,:,:)   :: dPsidx, dPsidy
   
   !! Grid 
   
   !! Domain length, width and depth
   real(4)                                  ::  xmin, xmax, ymin, ymax, zmin, zmax
   
   !! For stream function
   real(4)                                  ::  Lx, kx, Ay, Hz, lambda, cx, sc, psi0
   
   !!                                           km/d -> m/s
   real(4), parameter                       ::  kmd = 1000.0 / (60.0*60.0*24.0)
   
                                              
   !! Local variables for stream function
   real(4), allocatable, dimension(:,:,:,:) ::  psi
   
#ifdef tempsalt
   !!                                           high values of T and S 
   real(4), parameter                       ::  thigh = 25.0 !![C]
   real(4), parameter                       ::  shigh = 35.0 !![PSU]
   
   !!                                           low values of T and S
   real(4), parameter                       ::  tlow = 10.0  !![C]
   real(4), parameter                       ::  slow = 33.0  !![PSU]
#endif
   
   !!------------------------------------------------------------------------------
   
   namelist /GRIDSIZE/  xmin, xmax, ymin, ymax, zmin, zmax
   namelist /WAVE/      Lx, Ay, Hz, lambda, cx, sc
   
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
   
   open (unit=11,file='projects/bower-gulf/namelist',status='old',delim='apostrophe')
   read (11,nml=WAVE)
   close(11)
   
   Lx     = Lx * 1000.0 !![m]
   Ay     = Ay * 1000.0 !![m]
   lambda = lambda * 1000.0 !![m]
   kx     = 2.0*PI/Lx !![m-1]
   cx     = cx * kmd  !![m/s]
   sc     = sc * kmd !![m/s]
   psi0   = sc * lambda !![m2/s]
   
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
   !! Calculate stream function for this time step
   !!
   
   allocate ( psi(imt,jmt,km,1) )
   
   do jk=1,km
      
      do jj=1,jmt
         
         do ji=1,imt
            
            x = float(ji) * dx     + xmin
            y = float(jj) * dy     + ymin
            z = float(jk) * dz(jk) + zmin
            
            
            !   !! y position of jet center
            !   yc = Ay * sin( kx * (x-cx*t) )
            !   
            !   alpha = atan( Ay * kx * cos(kx * (x-cx*t)) )
            !   beta  = (y - yc) / (lambda/cos(alpha))
            !
            !   !! Equation (1) from Bower (JPO,1991)  
            !   psi(ji,jj,jk,1) = psi0 * exp(-z/Hz) * (1.0 - tanh(beta))
            
            
               yc = Ay * sin( kx * x )
               alpha = atan( Ay * kx * cos(kx * x) )
               beta  = (y - yc) / (lambda/cos(alpha))
               
               !! Equation (2) from Bower (JPO, 1991)
               psi(ji,jj,jk,1) = psi0 * exp(-z/Hz) * (1.0 - tanh(beta)) + cx*y 
            
         end do
         
      end do
      
   end do

   
   !!------------------------------------------------------------------------------
   
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
   
   
   !!
   !!  Calculate uflux and vflux
   !!
   
   do jk=1,km
      
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
      
   end do
   
   !! Set zero flux at y boundaries
   vflux(1:imt,  0,1:km,nsp) = 0.0
   vflux(1:imt,jmt,1:km,nsp) = 0.0
   
   !!------------------------------------------------------------------------------
   
   !!
   !! Set up a simple distribution of temperature and salinity 
   !! where temperature and salt in low at high y and high and low y
   !! Meridional variations are set by a arctan function
   !! Uniform with depth
   !! 
   
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
            
            !! y position of jet center
            yc = Ay * sin( kx * (x-cx*t) )
            
            alpha = atan( Ay * kx * cos(kx * (x-cx*t)) )
            beta  = (y - yc) / (lambda/cos(alpha))
            
            toff = 0.5 * (tmax + tmin)
            tscale = tmax-tmin
            
            tem(ji,jj,jk,nsp) = (-1.0) * atan(beta) / PI * (thigh-tlow) + 0.5*(thigh+tlow)
            sal(ji,jj,jk,nsp) = (-1.0) * atan(beta) / PI * (shigh-slow) + 0.5*(shigh+slow)
            
         end do
      
      end do
   
   end do
   
   
   do jj=1,jmt
      
      do ji=1,imt
         
         call statvd(tem(ji,jj,:,nsp),sal(ji,jj,:,nsp),rho(ji,jj,:,nsp),km)
      
      end do
      
   end do
   
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
   
   deallocate ( psi, dPsidx, dPsidy )
   
   return
   
   !!------------------------------------------------------------------------------
   
end subroutine readfields

