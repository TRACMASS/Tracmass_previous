! This subroutine is not updated anymore and will be replaced by the TRACMASS manual
!                                  
!                                TRACMASS
!                       Lagrangian trajectory code 
!                             in fortran 95
!                                 by
!      Kristofer Doos, Bror Jonsson, Hanna Corell, Andrew Coward, 
!             Pedro de Vries, Donatella Faggioli, etc.
!
!  Department of Meteorology, Stockholm University (doos@misu.su.se)
!
! Doos,K., 1995: Inter-ocean exchange of water masses. JGR Vol. 100, No. C7, 13499-13514.
!
! Vries,P.de and K.Doos, 2001: Calculating Lagrangian trajectories using time-dependent
!                              velocity fields. JAOT Vol. 18, No. 6, 1092-1101.
!
! Units in SI
! uflux zonal transport in m3/s (zonal velocity * dy * dz)
! vvlux meridional transport in m3/s (meridional velocity * dx * dz)
! h sea surface elevation in m
! grid coordinates in model index coordinates
!
! the OGCM velocities are interpolated on a C-grid box in and multiplied by the box wall
! in order to have transports in m3/s in the zonal, meridional and vertical directions
!
!   ja  -----------------vflux(ia,ja)-----------------
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!    (x0,y0)                                         |
!       |          *                                 |
!       |                  *                         |
!       |                         *                  |
!  uflux(ia-1,ja)                x         *         uflux(ia,ja)
!       |                  h(ia,ja)         *        |
!       |                                       *    |
!       |                                          * |
!       |                                         (x1,y1)
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!  ja-1 ----------------vflux(ia,ja-1)----------------
!     ia-1                                           ia
!
! Inside a grid box the zonal transport can be obtinaed by interpolating linearly
! uflux(x)=uflux(ia-1)+(x-x(ia-1))(uflux(i)-uflux(ia-1))  where x is in model index coordinates
! Since u=dx/ds => dx/ds + beta*x + delta = 0  where 
! beta=uflux(ia-1)-uflux(ia) 
! delta=-uflux(ia-1)-beta*x(ia-1)
! ds=dt/(dx*dy*dz) 
!
! The differntial equation can be solved analytically with initial and boundary conditions
! x(s0)=x0 =>
! x(s)=(x0+delta/beta)*exp(-beta(s-s0)) - delta/beta
! The time s1 when the trajectory reaches a zonal wall is
! s1=s0-1/beta*( log(x1+delta/beta) - log(x0+delta/beta) )
!
! ds=dt/(dx*dy*dz) Units in (s/m**3)
! tt time in seconds of the trajectory relative to the code start
! dt time step of trajectory iteration in seconds
! tseas (sec) velocity fields time step = time between velocity fields
! dtmin iterative time steping in seconds
! dsmin=dtmin/(dx*dy*dz) Units in (s/m**3)
!
! NNTRJ and NTRJ number of informations to be stored for each 
! trajectory between the dataset time steps (ints)
! nrj(ntrac,1)=ib  zonal index
! nrj(ntrac,2)=jb  meridional index
! nrj(ntrac,3)=kb  vertical index
! nrj(ntrac,4)=n   trajectory iterative step
! nrj(ntrac,5)=ts  model dataset time step of trajectory
! nrj(ntrac,6)=0 if trajectory  active & =1 if trajectory terminated & =2 if sedimented
! nrj(ntrac,7)=1 if trajectory active in data set time step else =0
! nrj(ntrac,8)=1,2,...,NEND index of end section for rerun in order to calculate
!                       Lagrangian stream functions
! trj(ntrac,1)=x1 zonal      model index coordinate of trajectory
! trj(ntrac,2)=y1 meridional model index coordinate of trajectory
! trj(ntrac,3)=z1 vertical   model index coordinate of trajectory
! trj(ntrac,4)=tt time of trajectory in seconds
! trj(ntrac,5)=subvol volume transport of trajctory in m3/s
! trj(ntrac,7)=t0 initial time in seconds of trajectory
!
!__________________________ Array dimensions defined in modules.f95_______________________
! IMT= zonal model array dimension
! JMT= meridional    "
! KM=  vertical      "
! MR= number of density/pressure, temperature & salinity/specific humidity levels
! NTRACMAX= maximum number of trajectories
! NST=number of time levels to be stored in centrel memory 1=stationary, 2=time dependent
! LBT=stream function separation dimension. If LBT>1 one must use Drerun if time dependent
! NEND=number of end sections that are defined in main
! LOV=Stream function dimension with LOV=3 for Dstreamts else LOV=1
!_________________________________________________________________________________________
! Run the programme with gotraj with the following precompilation options:
!
!__________ Ocean and atmospheric GCM possibilities
!
! -Dorca        NEMO with one of the following ORCA grids
! -Dorca2       NEMO with 2    degree ORCA2 grid 
! -Dorca1       NEMO with 1    degree ORCA1 grid
! -Dorca05      NEMO with 0.5  degree ORCA05 grid
! -Dorca025     NEMO with 0.25 degree ORCA025 grid
! -Dorca12      NEMO with 1/12 degree ORCA12 grid
! -Docc         OCCAM 66 levles 1/4 deree resolution
! -Drco         RCO model 2 nm
! -Dfors        Forsmark model
! -Dsimp        Simpevarp model
! -Dtes         Test basin with analytical time depenent velocity fields
! -Dtun         ROMS GCM for the Med 
! -Datm         The AGCM (IFS) from ECMWF, which is part of EC-Earth
!
!__________ Time possibilities
! -Dtime        Time changing velocity fields
! -Dstat        Stationary velocity fields (does not work at the moment)
! -Dtimeanalyt  Analytical time dependent option is based on Vries and Doos (JAOT,2001)  
!               It has to be reimplemented from old TRACMASS code 
!_________________________________________________
! -Dsediment    Sediment code developed for RCO
! -Dtwodim      Two-dimensional trajectories, which do not change depth. 
!               The vertical velocity is simply set to zero.
! -Drerun       Stores the Lagrangian stream functions as a function of the end
!               positions that has been calculated in an identical previous trajectory run
! -Dturb        Adds parameterised sub-grid turbulent velocities u' and v'
!
!_______________ Output possibilities
! -Dstreamxy    Calculates the barotropic stream function.  
! -Dstreamv     Calculates the vertical stream function as a function of depth
! -Dstreamr     ------------------------ " ----------------------------- density 
!               the latter two have to be used with
! -Dtempsalt    Calculates the temperature, salinity and density for the ocean
!                          and temperature, humidity and pressure for the atmosphere
! -Ddensity     Calculates only the density along the trajectory.  
!
! -Dtracer      Stores trajectory particle positions as a simulated tracer
!
!
!OOOOOOOOOOOOOOOOOOOOOOOOOOO  OLD precompilation options OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

! -Dmean 
! -Dfiveday
! -Dsigma 
