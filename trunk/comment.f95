!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890x
!                                  
!                                TRACMASS
!                       Lagrangian trajectory code 
!                             in fortran 95
!                                 by
!      Kristofer Döös, Bror Jönsson, Hanna Kling, Andrew Coward, 
!             Pedro de Vries, Donatella Faggioli, etc.
!
!  Department of Meteorology, Stockholm University (doos@misu.su.se)
!
! Döös,K., 1995: Inter-ocean exchange of water masses. JGR Vol. 100, No. C7, 13499-13514.
!
! Vries,P.de and K.Döös, 2001: Calculating Lagrangian trajectories using time-dependent
!                              velocity fields. JAOT Vol. 18, No. 6, 1092-1101.
!
! Units in SI
! u zonal transport in m3/s (zonal velocity * dy * dz)
! v meridional transport in m3/s (meridional velocity * dx * dz)
! h sea surface elevation in m
! grid coordinates in model index coordinates
!
! the OGCM velocities are interpolated on a C-grid box in and multiplied by the box wall
! in order to have transports in m3/s in the zonal, meridional and vertical directions
!
!   ja  -------------------v(ia,ja)-------------------
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!    (x0,y0)                                         |
!       |          *                                 |
!       |                  *                         |
!       |                         *                  |
!   u(ia-1,ja)                x         *         u(ia,ja)
!       |                  h(ia,ja)         *        |
!       |                                       *    |
!       |                                          * |
!       |                                         (x1,y1)
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!       |                                            |
!  ja-1 ------------------v(ia,ja-1)------------------
!     ia-1                                          ia
!
! Inside a grid box the zonal transport can be obtinaed by interpolating linearly
! u(x)=u(ia-1)+(x-x(ia-1))(u(i)-u(ia-1))  where x is in model index coordinates
! Since u=dx/ds => dx/ds + beta*x + delta = 0  where 
! beta=u(ia-1)-u(ia) 
! delta=-u(ia-1)-beta*x(ia-1)
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
! trj(ntrac,6)=arct length of the trajectory computed by arclength.f
! trj(ntrac,7)=t0 initial time in seconds of trajectory
!
!__________________________ Array dimensions defined in param.h __________________________
! IMT= zonal model array dimension
! JMT= meridional    "
! KM=  vertical      "
! JMAX= maximum meridional orm array dimension =JMT if entire model domain is to be used
! MR= number od density, temperature & salinity levels
! NTRACMAX= maximum number of trajectories
! NST=number of time levels to be stored in centrel memory 1=stationary, 2=time dependent
! LBT=stream function separation dimension. If LBT>1 one must use Drerun if time dependent
! NEND=number of end sections that are defined in main
! Stream function dimension with LOV=1 for Ddensity and LOV=3 for Dstreamts



!_________________________________________________________________________________________
! Run the programme with gotraj with the following precompilation options:
!
!__________ OGCM possibilities
! -Dorca        ORCA5/OPA model
! -Docc66       Occam 66 levles 1/4 deree resolution
! -Drco         RCO model 2 nm
! -Dfors        Forsmark model
! -Dsimp        Simpevarp model
! -Dtes         Test basin with analytical time depenent velocity fields
! -Dtun         ROMS GCM for the Med 
! 
! -Dtime        Time changing velocity fields
! -Dstat        Stationary velocity fields
!
! -Dsediment    Sediment code developed for RCO
!
! -Dstreamxy    Calculates the barotropic stream function.  
! -Dstreamv     Calculates the vertical stream function as a function of depth
! -Dstreamr     ------------------------ " ----------------------------- density 
!               the latter two have to be used with
! -Dtempsalt    Calculates the temperature, salinity and the density 
! -Ddensity     Calculates only the density along the trajectory.  
!
! -Dmean 
!
! -Dtracer      Stores a simulated tracer
!
! -Drerun       Stores the Lagrangian stream functions as a function of the end
!               positions that has been calculated in an identical previous trajectory run
!
!
!OOOOOOOOOOOOOOOOOOOOOOOOOOO  OLD precompilation options OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
!
! -Dtwodim: Uses a two-dimensional velocity field. 
!           surface layer velocity. Trajectory stops when converged.
!
! -Deqatl: release trajectories only from the Equator in the Atlantic
!          Ocean. To be used only with -Dmod1
!
! -Dinitone: Trajectories are only generated at the very first time step
!
! 
! -Dinterp1: used for working with NADW. The depth of the 
! isopycnal 27.625 is it stored in two files 
! 
! -Dnome: It stops all the trajectories enetring the Mediterranean Sea. 
! The variable  imed counts the number of trajectories entering.
! 
! -Dupwelling: It is used to calculate the upwelling of NADW. The result
! is stored in binary format in unit 58, file data.upw
! 
! -Dpoints: It is used for release trajectories from a fixed section.
! 
! -Dseason: To be used in conjuction with  -Dtime. It will check 
!           only the change in season.
! 

!
! -Deqatl: release trajectories only from the Equator in the Atlantic
!          Ocean. To be used only with -Dmod1
!
! -Dline: it is used for considering only trajectories below a certain value in 
! latitude (only for model 2)
! 
! -Dinterp1: used for working with NADW. The depth of the 
! isopycnal 27.625 is it stored in two files 
! (/mare/faggioli/depthtime1.dat, 
! /mare/faggioli/depthtime2.dat) and it is assigned to unit 85. 
! The variables read are called  depth1(imt1,jmt1) and  depth2(imt2,jmt2).
! 
! -Dnome: It stops all the trajectories enetring the Mediterranean Sea. 
! The variable  imed counts the number of trajectories entering.
! 
! -Dpoints: It is used for release trajectories from a fixed section.
! 

