

Lagrangian trajectory code for ocean and atmospheric general circulation models

TRACMASS makes it possible to estimate water paths, Lagrangian stream functions (barotropic and overturning), exchange times, particle sedimentation, etc.

TRACMASS has been used in studies of the global ocean circulation, of sea circulation in The Baltic, The Mediterranean and in coastal regions.

The code is written in FORTRAN 95 with modules and runs on UNIX platforms such as MAC OS X and Linux.

TRACMASS has been set up to run with velocities integrated with models such as NEMO, ORCA, ROMS, POM, MICOM, IFS-ECMWF and EC-Earth.


Installation on Mac Os X
========================

Download and install macports (www.macports.org)
Install the following ports:

 sudo port install openmpi +gcc45 +valgrind
 sudo port install netcdf +openmpi +netcdf4 
 sudo port install netcdf-fortran +gcc45 +openmpi
 sudo port install git-core 
 sudo port install git-extras


Compile
=======
Copy Makefile_tmpl to Makefile, change PROJECT and CASE in Makefile to your projectname.

run:

 make
 ./runtrm
