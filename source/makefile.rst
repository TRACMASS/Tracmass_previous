

Project Makefile options
========================


Variables set from command line
-------------------------------

 .. make:var:: INPUT_INT1    

    Define first Variable set from command line. Use 'dummy' if not
    used.

 .. make:var:: INPUT_INT2    
    
    Define second Variable set from command line. Use 'dummy' if not
    used.


Time interpolation
------------------

 .. make:var::  -Dtimestep   

    Time steps with analytical stationary scheme differential Eqs.

 .. make:var::  -Dtimeanalyt 

    Analytical time scheme used to solve the

 .. make:var::  -Dregulardt      

    Regular time steps to be used with -Dtimestep


.. _makeoutput:

Output formats
--------------

 .. make:var:: -Dtextwrite

    Write results to textfile

 .. make:var:: -Dcsvwrite

    Write results to comma-separated file

 .. make:var:: -Dbinwrite

    Write results to binaryfile

 .. make:var:: -Dmysqlwrite

    Write results to mysql


Stream functions
----------------

 .. make:var:: -Dstreamxy      
    
    Calculates the barotropic stream function.

 .. make:var:: -Dstreamr     

    vertical stream function, z=density

 .. make:var::  -Dstreamts

    Lagrangian thermohaline/hydrothermal stream function

 .. make:var:: -Dstreamv      

    vertical stream function, z=depth


 .. make:var:: -Drerun

    Stores the Lagrangian stream functions as a function of the end
    positions that has been calculated in an identical previous run.   

 .. make:var:: -Dinitxyt

    Start trajectories at given positions and times


Vertical advection
------------------

 .. make:var:: -Dtwodim

    Turn off vertical velocities.

 .. make:var:: -Dfull_wflux 

    Use a full 3D wflux field.

 .. make:var:: -Dexplicit_w 

    Use a given vertical velocity.


Grid shape
----------

 .. make:var::  -Dvarbottombox

    Variable bottom box to match actual depth

 .. make:var::  -Dfreesurface

    Variable bottom box to match actual depth

 .. make:var::  -Dzgrid1D

    Cell depths defined as individual column vector for z-coordinates

 .. make:var::  -Dzgrid3D

    Cell depths defined as 3D grid (for sigma)

 .. make:var::  -Dzgrid3Dt

    Cell depths 3D and time interp. (for atm)


Diffusion
---------

 .. make:var::  -Dturb

    Adds subgrid turbulent velocities 

 .. make:var::  -Ddiffusion

    Adds a diffusion on trajectory

 .. make:var::  -Danisodiffusion

    Adds an anisotropic diffusion on trajectory


Special functions
-----------------
 .. make:var::  -Dtempsalt

    Include temp and salt

 .. make:var::  -Dtracer

    Stores a simulated tracer

 .. make:var::  -Dsediment

    Sediment code developed for RCO

 .. make:var::  -Dorca025

    orca1 or orca025 or orca025l75h6

 .. make:var:: -Ddrifter 

    surface drifter depth average of uflux/vflux
