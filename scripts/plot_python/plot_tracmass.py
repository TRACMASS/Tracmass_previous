#
#
# Basic script to read and plot TRACMASS trajectories
#
# Depends on:
#   * numpy
#   * matplotlib 
#   * cython
#
# Author : Joakim Kjellsson, 2019 (at 35000 feet over the Indian Ocean)
#
#

import os,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib as mpl
import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)
import tracmass_functions

def read_bin(filename,undef=-999.,mode='old',cython=True):   
   if mode == 'new':
      dtype = np.dtype([('ntrac','i4'), ('ints','f8'), 
                        ('x','f4'), ('y','f4'), ('z','f4'),
                        ('tem','f4'), ('sal','f4')])
   else:
       dtype = np.dtype([('ntrac','i4'), ('ints','f8'),('x','f4'), ('y','f4'), ('z','f4')])
   
   data = np.fromfile(open(filename), dtype)
   
   ntracmax = data['ntrac'].max() + 1    
   
   if not cython:
      nsteps = 0
      for ntrac in xrange(1,ntracmax+1):
         steps = np.where( data['ntrac'] == ntrac )[0]
         nsteps = max(nsteps,steps.shape[0])
      
      # array to store trajectories in
      traj = np.ma.ones((ntracmax,nsteps,10)) * undef
      
      for ntrac in xrange(1,ntracmax+1):
         steps = np.where( data['ntrac'] == ntrac )[0]
         steps = steps[:]
         traj[ntrac-1,0:steps.shape[0],0] = data['x'][steps]
         traj[ntrac-1,0:steps.shape[0],1] = data['y'][steps]
         traj[ntrac-1,0:steps.shape[0],2] = data['z'][steps]
         traj[ntrac-1,0:steps.shape[0],3] = np.arange(0,steps.shape[0])
         if np.mod(ntrac,100) == 0: 
            print('Processing : ',ntrac)
         if mode == 'new':
            traj[ntrac-1,0:steps.shape[0],4] = data['tem'][steps]
            traj[ntrac-1,0:steps.shape[0],5] = data['sal'][steps]
   
   elif cython:
      rawdata = np.ones((data['ntrac'].shape[0],7)) * undef
      rawdata[:,0] = data['ntrac'][:]
      print('x min max ',data['x'][:].min(),data['x'][:].max())
      rawdata[:,1] = data['ints'][:]
      rawdata[:,2] = data['x'][:]
      rawdata[:,3] = data['y'][:]
      rawdata[:,4] = data['z'][:]
      if mode == 'new':
         rawdata[:,5] = data['tem'][:]
         rawdata[:,6] = data['sal'][:]
         print(rawdata[:,5].min(),rawdata[:,5].max())
      
      rawdata = np.array(rawdata,dtype='float32')      
      nsteps = tracmass_functions.find_numsteps(rawdata,ntracmax=ntracmax)
      print(nsteps)
      nstepsmax = nsteps.max()
      traj   = tracmass_functions.sort_particles(rawdata,ntracmax=ntracmax,nsteps=nstepsmax)
   
   # mask all invalid data points
   traj = np.ma.masked_where(traj == -999.,traj)
   
   return traj 
   
def read_lagrpsi(filename,mode='stxy'):
   if mode == 'stxy':
      dtype = np.dtype('(292,362)f4')
   
   data = np.fromfile(open(filename), dtype)
   
   print(data.shape)
   stxyy = data[:,:,:]
   #stxyx = data[5:10,:,:]
   return stxyy#,stxyx

#
# Directories etc.
#

tracmass_dir = '/gws/nopw/j04/aopp/joakim/tracmass_out/nemo/19920105-0000/' #'/Users/jkjellsson/Downloads/nemo/20000101-0000/'
trm_run      = 'orca12_overflow-backward'

runfile = '%s/%s_run.bin' % (tracmass_dir,trm_run) 
infile  = '%s/%s_in.bin'  % (tracmass_dir,trm_run) 
outfile = '%s/%s_out.bin' % (tracmass_dir,trm_run) 
   
plot_traj_map = True

if 1:   
   if plot_traj_map:
      traj = read_bin(runfile,mode='old')
      
      ##
      ## Plot trajectories on map
      ##  
      fig1 = plt.figure()
      ax1  = fig1.add_subplot(111)
      
      fig2 = plt.figure()
      ax2  = fig2.add_subplot(111)
      
      vmin1 = 0 
      vmax1 = 120
      vmin2 = 0
      vmax2 = 75
      cmap1 = plt.cm.winter
      cmap2 = plt.cm.winter
      norm1 = mpl.colors.Normalize(vmin=vmin1, vmax=vmax1)
      norm2 = mpl.colors.Normalize(vmin=vmin2, vmax=vmax2)
      label1 = 'Time step' 
      label2 = 'Vertical level' 
      xin = np.array([])
      xut = np.array([])
      yin = np.array([])
      yut = np.array([])
      
      for ntrac in range(0,traj.shape[0]):
         x = traj[ntrac,:,0].compressed() 
         y = traj[ntrac,:,1].compressed()
         z = traj[ntrac,:,2].compressed()
         z1 = traj[ntrac,:,3].compressed()
         z2 = np.array(z,dtype='float32')         
         if x.shape[0]>0:
            print('Plot ntrac ',ntrac)
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=cmap2,norm=plt.Normalize(vmin2, vmax2),linewidth=0.2)
            ax2.add_collection(lc)
            lc.set_array(z2)
            xin = np.append(xin,x[0])
            xut = np.append(xut,x[-1])
            yin = np.append(yin,y[0])
            yut = np.append(yut,y[-1])
            
         if x.shape[0]>0:         
            sc = ax1.scatter(x,y,s=5,c=z1,cmap=cmap1,vmin=vmin1,vmax=vmax1,edgecolors=None,zorder=1)
            #ax2.plot(x,y,color=cmap(norm(z[0])),lw=0.5)
            #ax1.scatter([x[0]],[y[0]],c='w',marker='o',edgecolors='k')
            #ax1.scatter([x[-1]],[y[-1]],c='w',marker='d',edgecolors='k')
      
      for ax in [ax1,ax2]:
         ax.scatter(xin,yin,c='w',marker='o',edgecolors='k',zorder=2)                                              
         ax.scatter(xut,yut,c='w',marker='d',edgecolors='k',zorder=2)      
         ax.set_title('Trajectories from run: %s ' % (trm_run,))      
         ax.set_xlabel('Model zonal index')
         ax.set_ylabel('Model meridional index')
      
      cb1 = plt.colorbar(sc, ax=ax1, orientation='vertical', pad=0.) 
      cb1.set_label(label1)      
      cb2 = plt.colorbar(lc, ax=ax2, orientation='vertical', pad=0.)
      cb2.set_label(label2)
      for fig in [fig1,fig2]:
         fig.tight_layout()
      #figname1 = 'traj_map_dots_%s_%s.png' % (trm_runs[jf],trm_time[jf])
      #figname2 = 'traj_map_lines_%s_%s.png' % (trm_runs[jf],trm_time[jf])
      #fig1.savefig(figname1,format='png')
      #fig2.savefig(figname2,format='png')
   
plt.show()
   
   
