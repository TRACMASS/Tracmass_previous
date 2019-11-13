TRACMASS is a Lagrangian trajectory code for ocean and atmospheric general circulation models. The code makes it possible to estimate water paths, Lagrangian stream functions (barotropic and overturning), exchange times, particle sedimentation, etc. TRACMASS has been used in studies of the global ocean circulation, of sea circulation in The Baltic, The Mediterranean and in coastal regions.

The code is written in FORTRAN 95 with modules and runs on UNIX platforms such as MAC OS X and Linux.

TRACMASS has been set up to run with velocities integrated with models such as NEMO, ORCA, ROMS, POM, MICOM, IFS-ECMWF and EC-Earth.


Quickstart
==========

Prerequisites to run on Mac OS X
--------------------------------

Download and install macports (www.macports.org)
Install the following ports:

```sh
sudo port install openmpi +gcc45 +valgrind
sudo port install netcdf +openmpi +netcdf4 
sudo port install netcdf-fortran +gcc45 +openmpi
sudo port install git-core 
sudo port install git-extras
```

Download test data
------------------

You can find some input data for testing the code on 

```bash
https://www.dropbox.com/sh/b6leim7tb99i3hm/AAA7M0mxuYTwJLKjqmLnoJuca?dl=0
```

Before doing any analysis we recommend to download some of the test data and make sure TRACMASS is working properly. 

Compile the code
----------------

Copy Makefile_tmpl to Makefile. 
Choose a change PROJECT and CASE in Makefile to your projectname.
If you are using netCDF version >4.4 you can set NETCDFLIBS to "automatic-44", in which case netCDF paths will be taken from nc-config and nf-config. 

run:

```bash
make
./runtrm
```


Change log
==========


6.0.0 First public release