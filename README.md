TRACMASS is a Lagrangian trajectory code for ocean and atmospheric general circulation models. The code makes it possible to estimate water paths, Lagrangian stream functions (barotropic and overturning), exchange times, particle sedimentation, etc. TRACMASS has been used in studies of the global ocean circulation, of sea circulation in The Baltic, The Mediterranean and in coastal regions.

The code is written in FORTRAN 95 with modules and runs on UNIX platforms such as MAC OS X and Linux.

TRACMASS has been set up to run with velocities integrated with models such as NEMO, ORCA, ROMS, POM, MICOM, IFS-ECMWF and EC-Earth.

WARNING
=======

This TRACMASS version is an alpha version. 
It is still in development. 

Quickstart
==========

Download the code
-----------------

```bash
git clone https://github.com/TRACMASS/tracmass.git
```

Compile
-------

Enter the tracmass directory and copy the template Makefile

```bash
cd tracmass
cp Makefile_tmpl Makefile
```

Modify the Makefile to fit your system. 
You will need to set ARCH, which is the name of your system, i.e. macports, tetralith etc. 
You will also need to configure how TRACMASS should find the netCDF libraries, if at all. 
For most systems, we recommend the option automatic-44.
Then you can run the make command. 

```bash
make
```

Running a first test case
-------------------------

We recommend testing that TRACMASS was properly compiled by letting PROJECT and CASE be "theoretical" in the Makefile (which is the default). 
In this case, TRACMASS will use a simple oscillating velocity field to trace trajectories. 
You can run this case by setting

```bash
./runtracmass
```

Download test data
------------------

You can find some input data for testing the code on 

```bash
https://www.dropbox.com/sh/b6leim7tb99i3hm/AAA7M0mxuYTwJLKjqmLnoJuca?dl=0
```

This test data includes data from NEMO, IFS and ROMS. 
Before doing any analysis we recommend to download some of the test data and make sure TRACMASS is working properly. 

In order to set up TRACMASS to run trajectories on e.g. NEMO data, you will need to change PROJECT and CASE to NEMO, and then re-compile the code. 

```bash 
make clean
make
./runtracmass
```

Run your analysis
-----------------

If you wish to run TRACMASS using fields from NEMO, ROMS, IFS or CMEMS satellite altimetry, you may use these projects as is (see projects directory). 
If you wish to run another case or a very specific case of the above models, you may create your own project in the projects directory. 
To run with e.g. your own NEMO data, you will need to modify the projects/nemo/nemo.in namelist to suit your needs. 


Change log
==========


6.0.0 First public release
7.0.0_alpha New features for parametrising sub-grid scale motions, tracing water in the atmosphere, simplified set-up for users + code clean–up + bug-fixes.  


Prerequisites to run on Mac OS X
================================

Download and install macports (www.macports.org)
Install the following ports:

```sh
sudo port install openmpi +gcc45 +valgrind
sudo port install netcdf +openmpi +netcdf4 
sudo port install netcdf-fortran +gcc45 +openmpi
sudo port install git-core 
sudo port install git-extras
```
