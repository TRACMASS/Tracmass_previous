.. TRACMASS documentation master file, created by
   sphinx-quickstart on Tue Sep  2 11:54:43 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the TRACMASS documentation
=====================================


TRACMASS is a Lagrangian trajectory code for ocean and atmospheric
general circulation models. The code makes it possible to estimate
water paths, Lagrangian stream functions (barotropic and overturning),
exchange times, particle sedimentation, etc. TRACMASS has been used in
studies of the global ocean circulation, of sea circulation in The
Baltic, The Mediterranean and in coastal regions. The code is written
in FORTRAN 95 and runs on UNIX platforms such as MAC OS X and
Linux. TRACMASS has been set up to run with velocity fields from
general circulation models such as NEMO, ORCA, ROMS, POM, MICOM,
IFS-ECMWF, MOM, HYCOM, MITGCM, ECCO, and EC-Earth.


Quickstart
----------

You need gfortran (gcc) or iFort to compile the code and run the testcases. We
suggest macports (http://www.macports.org) if you use Mac
OS X. Install netcdf as well if you are going use TRACMASS with real data.

Download the source code from
https://github.com/TRACMASS/tracmass/archive/master.zip 

or clone the project::

  git clone https://github.com/TRACMASS/tracmass.git

Make a copy of Makefile_tmpl named Makefile, edit the new file to
select your compiler, compile, and run:

.. code-block:: bash

  cp Makefile_tmpl Makefile
  nano Makefile
  # (Uncomment your compiler, save and close with ctrl-x y)
  make
  ./runtrm

There should now be a couple of new files in the folder. You can
use python to plot the advected trajectories:

.. code-block:: python
   
   import numpy as np
   import pylab as pl

   data = np.loadtxt('analytical_run.csv', delimiter=",")
   pl.plot(data[:,2], data[:,3])

Or matlab:

.. code-block:: matlab

   Someone write this!


Contents:

.. toctree::
   :maxdepth: 2
            
   time
   fileformats
   envvar
   makefile
   namelist
   sourcecode

.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. #* :ref:`search`

