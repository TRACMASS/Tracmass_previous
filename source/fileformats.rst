
Saving particle postions
========================

TRACMASS can save calculated particle postion at different frequencies. This timing and the location of files is controlled by the case namelist (:ref:`INITRUNWRITE`). The resulting output files can also have different formats, which is selectes in each project's Makefile.prj file (:ref:`makeoutput`). They all save the same type of data:

**ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens**
where

:ntrac:      Trajectory number
:niter:      TRACMASS code iteration (only important for TRACMASS modellers)
:x1:         Zoonal position of the trajectory particle
:y1:         Meridional position of the trajectory particle
:z1:         Vertical position of the trajectory particle
:tt:         Time of the trajectory particle (in days)
:t0:         Initial time of the trajectory particle
:subvol:     The "volume" transport in m3/s of the trajectory

The following columns are added if tempsalt is selected:

:temp:       Temperature of the trajectory particle
:salt:       Salinity/specific humidity of the trajectory particle
:dens:       Density of the trajectory particle



Ascii (fortran) format
----------------------
This is a native fortran text format where each column is right aligned. Best used when the data is postprocessed by other fortran softwares (e.g. Kristofers stuff). Each file is very large and take a long time to read. example::


CSV (Comma-separated values) format
-----------------------------------
CSV files are great if you want to look at small datasets quickly. Works well with python, R,  and matlab. Pretty safe when transferring data between computers with different processor architecture and for longterm storage. CSV files are also very large and slow to be processed. Avoid this format if you have large datasets. Example of a file::

  1,27.00000,9.12500,49.12500,4.50000
  1,27.00000,9.13220,49.12500,4.50000
  1,27.00000,9.13940,49.12500,4.50000
  1,27.00000,9.14660,49.12499,4.50000


Binary format
------------- 
The most efficient file format to use is binary. These files are not directly readble and need to be imported using either numpy or matlab. They are also endian depended, which means that data genreated on one computer might not be directly readable on a computer of a different architecture. 

We don't provide the capacity to save to netcdf since this format are less optimal for his usage. The best self contained format to us would be HDF which we  plan to support in the future. 


Filenames
---------

Each run generates the following files

:casename_ini.ext: Where and when particles are seeded
:casename_kll.ext: Where and when particles are killed
:casename_out.ext: Where and when particles leave the domain
:casename_run.ext: All particle trajectories.

Where *ext* is either asc, csv, or bin depending on format.

It is the run-file that you normally want for postprocessing, the others are mainly for diagnostics. The filename might also include a time value and information about chunks if very large runs are split (see XXX).
