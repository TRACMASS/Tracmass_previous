#
#
# Basic script to read and plot TRACMASS trajectories
#
# Depends on:
#   * numpy       - Array handling in python
#   * matplotlib  - For plotting
#   * pandas      - Data processing
#   * netcdf4     - To read a topography/bathymetry and plot as background
#
# Author : Joakim Kjellsson, 2019 (at 35000 feet over the Indian Ocean)
#
#

import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as mpl
from netCDF4 import Dataset

import matplotlib
print(matplotlib.matplotlib_fname())

def read_txt(filename,mode='old'):
   if mode == 'old':
      intrac = 0
      inn    = 1
      ix     = 2
      iy     = 3
      iz     = 4
      itt    = 5
      itt0   = 6
      ivol   = 7
      item   = 8
      isal   = 9
      idens  = 10
      df = pd.read_table(filename, delim_whitespace=True)
      print(df)
   
   #ntracmax = df[intrac].max() + 1    
   #nsteps = 0
   #for ntrac in range(1,ntracmax+1):
   #   steps = pd.where( df[intrac] == ntrac )[0]
   #   nsteps = max(nsteps,steps.shape[0])
   
   #df.rename({0: "ntrac", 1: "nn", 2: "x", 3: "y", 4:"z", 5:"tt",6:"t0",7:"vol",8:"tem",9:"sal",10:"dens"},axis=1,errors="raise")
   #print(df['ntrac'])
   #print(df)
   
   return df
   

def read_bin(filename,undef=-999.,mode='old'):   
   if mode == 'new':
      dtype = np.dtype([('ntrac','i4'), ('ints','f8'), 
                        ('x','f4'), ('y','f4'), ('z','f4'),
                        ('tem','f4'), ('sal','f4')])
   else:
       dtype = np.dtype([('ntrac','i4'), ('ints','f8'),('x','f4'), ('y','f4'), ('z','f4')])
   
   data = np.fromfile(open(filename), dtype)
   
   ntracmax = data['ntrac'].max() + 1    
   nsteps = 0
   for ntrac in range(1,ntracmax+1):
      steps = np.where( data['ntrac'] == ntrac )[0]
      nsteps = max(nsteps,steps.shape[0])
   
   # array to store trajectories in
   traj = np.ma.ones((ntracmax,nsteps,10)) * undef
   
   for ntrac in range(1,ntracmax+1):
      steps = np.where( data['ntrac'] == ntrac )[0]
      steps = steps[:]
      traj[ntrac-1,0:steps.shape[0],0] = data['x'][steps]
      traj[ntrac-1,0:steps.shape[0],1] = data['y'][steps]
      traj[ntrac-1,0:steps.shape[0],2] = data['z'][steps]
      traj[ntrac-1,0:steps.shape[0],3] = np.arange(0,steps.shape[0])
      if mode == 'new':
         traj[ntrac-1,0:steps.shape[0],4] = data['tem'][steps]
         traj[ntrac-1,0:steps.shape[0],5] = data['sal'][steps]

   # mask all invalid data points
   traj = np.ma.masked_where(traj == -999.,traj)
   
   return traj 
   
def read_lagrpsi(filename,mode='stxy'):
   if mode == 'stxy':
      dtype = np.dtype('(292,362)f4')
   elif mode == 'styz':
      #dtype = np.dtype(['(5,75,292)f4','(5,75,362)f4'])
      #dtype = np.dtype([('styz','(5,75,292)f4'), ('stxz','(5,75,362)f4')])
      dtype = np.float32
      
   data = np.fromfile(open(filename), dtype)
   print(data)
   ii = 5*75*292
   styz = data[0:ii].reshape(5,75,292)
   print(styz[:,5,194])
   
   #print(data.shape)
   #stxyy = data[:,:,:]
   #stxyx = data[5:10,:,:]
   return styz

#
# Directories etc.
#

tracmass_dir = '/Users/joakim/Downloads/tracmass_out/'
trm_run      = 'orca1_test'

runfile = '%s/%s_run.bin' % (tracmass_dir,trm_run) 
inifile  = '%s/%s_ini.bin'  % (tracmass_dir,trm_run) 
outfile = '%s/%s_out.bin' % (tracmass_dir,trm_run) 

runfile_asc = '%s/%s_run.asc' % (tracmass_dir,trm_run) 
inifile_asc  = '%s/%s_ini.asc'  % (tracmass_dir,trm_run) 
outfile_asc = '%s/%s_out.asc' % (tracmass_dir,trm_run) 
   
psi_xz_yz = '%s/%s_psi_yz_xz.bin' % (tracmass_dir,trm_run)
   
plot_traj_map = True
plot_psi      = False

topofile = '/Users/joakim/data/tracmass_test_data/orca1/topo/bathy_level.nc'
nc = Dataset(topofile,'r')
topo = nc.variables['Bathy_level'][0,:,:]
nc.close()


if plot_traj_map:
   #
   # Read data into Pandas DataFrame
   #
   #traj = read_bin(runfile,mode='old')
   #traj = read_txt(runfile_asc,mode='old')
   traj1 = read_txt('/Users/joakim/Downloads/tracmass_out/orca1_test_run.asc',mode='old')
   traj2 = read_txt('/Users/joakim/Downloads/tracmass_out/orca1_test_newtracers_run.asc',mode='old')
   print(traj1['sal']-traj2['vosaline'])
   sys.exit()
   # Set min max for plot
   imin = traj['x'].min()
   imax = traj['x'].max()
   jmin = traj['y'].min()
   jmax = traj['y'].max()
   # Set some minimum size of plot
   imin -= 2
   jmin -= 2
   imax = max(imax, imin+10)
   jmax = max(jmax, jmin+10)
      
   #
   # Plot trajectories on map
   #  
   fig1 = plt.figure()
   ax1  = fig1.add_subplot(111)
      
   cf = ax1.contourf(topo)
   ax1.set_xlim([imin,imax])
   ax1.set_ylim([jmin,jmax])
      
   vmin1 = 0 
   vmax1 = 120
   vmin2 = 0
   vmax2 = 75
   cmap1 = plt.cm.winter
   cmap2 = plt.cm.winter
   norm1 = mpl.colors.Normalize(vmin=vmin1, vmax=vmax1)
   norm2 = mpl.colors.Normalize(vmin=vmin2, vmax=vmax2)
   
   # Group data by trajectory ID
   # and scatter plot in x,y 
   traj.groupby('ntrac').plot.scatter(x='x',y='y',ax=ax1)
   
   figname1 = 'traj_map_dots.png' 
   fig1.savefig(figname1,format='png')
   
   fig2 = plt.figure()
   ax2  = fig2.add_subplot(111)
   traj.groupby('ntrac').plot.scatter(x='y',y='tem',ax=ax2)
   print(traj['tem'][0:100])
   figname2 = 'traj_tem_dots.png' 
   fig2.savefig(figname2,format='png')
   

if plot_psi:
   
   psi = read_lagrpsi(psi_xz_yz,mode='styz')
   
   zplot = psi[0,:,:] * 1e-6
   for k in range(1,zplot.shape[0]):
      zplot[k,:] += zplot[k-1,:]
   
   fig1,ax1 = plt.subplots(1,1)
   cf = ax1.contourf(zplot)
   plt.colorbar(cf,ax=ax1)
      
   
plt.show()
   
   
