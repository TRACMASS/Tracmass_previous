import sys, os
import datetime
import glob
from datetime import datetime as dtm
from itertools import izip
import csv

import numpy as np
import pylab as pl
import scipy.io
from scipy.stats import nanmean, linregress
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree
import matplotlib.cm as cm

from hitta import GBRY
import projmaps, anim
import trm
import batch
import figpref
import mycolor

miv = np.ma.masked_invalid

class Discs(trm.Trm):
    """ Class to calculate residence times in circular discs.

    This class generates connectivity matrices from defined regions
    using the output from TRACMASS. The reagions are by default discs
    with a prescribed radius. The discs are packed in a semi-optimal
    fashion within the part of the grid defined by self.mask.

    Example to calculate the residence timesnd dt:

    cn = connect.Matrix('rutgersNWA','rutgersNWA')
    cn.load(jdstart=730120)
    cn.calc_conmat()
    """

    def __init__(self,projname,casename="", datadir="", datafile="",
                 ormdir="", griddir="",radius=2):
        super(Discs,self).__init__(projname, casename,
                                    datadir, datafile, ormdir)
        self.radius = radius
        self.add_default_regmask()
                
    def generate_regions(self, mask=[],regcoord='ij'):
        """Create discs defining the regions used for connectivities"""
        if len(mask)==0: mask = self.mask
        ncnt = 1
        nvec = []
        xvec = []
        yvec = []
        r = self.radius
        self.regmat = self.llon * 0
        if regcoord == "ij":
            regi = np.arange(r, self.imt) #
            regj = np.arange(r, self.jmt) #
        
        for n,(i,j) in enumerate([(i,j) for i in regi for j in regj]):
            if (mask[j-r:j+r, i-r:i+r]).all():
                xvec.append(i)
                yvec.append(j)
                nvec.append(ncnt)
                ncnt += 1
        self.disci = np.array(xvec)
        self.discj = np.array(yvec)
    
        duse = self.disci > -1
        for p in range(0,len(self.discj)):
            if duse[p]:
                mask = np.sqrt((self.disci[p+1:] - self.disci[p])**2 +
                               (self.discj[p+1:] - self.discj[p])**2) < r*2
                duse[p+1:][mask] = False
        self.disci = self.disci[duse]
        self.discj = self.discj[duse]
        self.discn = np.arange(1,len(self.discj)+1)
        self.discKD = cKDTree(list(np.vstack((np.ravel(self.disci),
                                              np.ravel(self.discj))).T))
        self.nreg = self.discn.max()+1

    @trm.Traj.trajsloaded
    def regvec_from_regions(self,mask=False):
        """Generate a vector with region IDs"""
        if not mask:
            self.add_default_regmask()
            mask = self.mask
        self.generate_regions(mask)
        dist,ij = self.discKD.query(list(np.vstack((self.x,self.y)).T), 1)
        self.reg = self.discn[ij]
        self.reg[dist>self.radius] = 0

    @trm.Traj.trajsloaded
    def calc_numpart(self,jd=None):
        """Calculate # of particles remaining in all regions at time jd."""
        jdstart = self.jd.min()
        if not hasattr(self,'reg'):
             self.regvec_from_regions()
        tmask1 = self.jd == jdstart
        tmask2 = self.jd == jd        
        ntracmax = max(self.ntrac[tmask1].max(), self.ntrac[tmask2].max())
        convec = np.zeros((2, ntracmax+1))
        convec[0,self.ntrac[tmask1]] = self.reg[tmask1]
        convec[1,self.ntrac[tmask2]] = self.reg[tmask2]
        convec = convec.astype(np.int)
        convec = convec[:,convec[0,:] == convec[1,:]]
        self.numpart = np.bincount(convec[0,:],minlength=self.nreg)

    @trm.Traj.trajsloaded
    def calc_decaymatrix(self):
        """Calculate the change of particles #'s in all regions over time"""
        self.jdvec = np.unique(self.jd)
        self.regvec_from_regions()
        self.decaymat = np.zeros((self.nreg, len(self.jdvec)))
        for n,jd in enumerate(self.jdvec):
            print n,jd
            self.calc_numpart(jd=jd)
            self.decaymat[:,n] = self.numpart
            
    def all_residencetimes(self):
        """Calculate residence times for all runs and regions"""
        self.listfiles()
        self.generate_regions()
        self.restime_mat = np.zeros((len(self.runfiles), self.nreg, 240))
        for n,f in enumerate(self.runfiles):
            self.load(filename=os.path.basename(f))
            self.calc_decaymatrix()
            self.restime_mat[n,1:,:] = self.decaymat[1:,:240]
        np.savez('restime_mat.npz', restime_mat=self.restime_mat)


    @trm.Traj.trajsloaded
    def calc_halftimes(self):
        if not hasattr(self, 'restime_mat'):
            print "loading restime_mat from file"
            self.restime_mat = np.load('restime_mat.npz')['restime_mat']
        self.jdvec = np.unique(self.jd)
        mat = self.restime_mat - self.restime_mat[:,:,-1][:,:,np.newaxis]
        mat = mat / mat[:,:,0][:,:,np.newaxis]
        meanmat = nanmean(mat,axis=0)
        xi = self.jdvec-self.jdvec[0]
        nreg = meanmat.shape[0]
        self.regtaus = np.zeros((nreg,))
        for i in np.arange(nreg):
            yi = meanmat[i,:]
            try:
                ipos = np.nonzero(yi>0.1)[0].max()
            except:
                continue
            k,m,_,_,_ = linregress(xi[:ipos], np.log(yi[:ipos]))
            tau = np.interp(-1, (xi*k+m)[::-1], xi[::-1])
            self.regtaus[i] = tau
        
    @trm.Traj.trajsloaded
    def tauregmap(self):
        if not hasattr(self,'regtaus'): self.calc_halftimes()
        if not hasattr(self,'reg'): self.regvec_from_regions()
        if not hasattr(self,'lon'): self.ijll()

        figpref.current()
        pl.close(1)
        pl.figure(1)
        mask = (self.jd==self.jd.min()) & (self.reg>0)
        x,y = self.gcm.mp(self.lon[mask], self.lat[mask])
        regs = self.reg[mask]
        taus = self.reg[mask].astype(np.float)
        for n,t in enumerate(self.regtaus):
            taus[regs==n] = t
        self.gcm.mp.scatter(x, y, 5, taus*24)
        pl.clim(0,30)
        cb = pl.colorbar(aspect=30,pad=0.02,ticks=[0,6,12,18,24,30])
        cb.ax.set_ylabel(r'$\tau$ (Hours)')
        self.gcm.mp.nice()
        pl.title('Residence times defined as e-folding half-life ')
        pl.savefig('residencetimes.png')

    def add_default_regmask(self):
        """Generate a mask to remove unwanted parts of the domain"""
        self.mask = ~self.gcm.landmask
        #self.mask = (self.gcm.depth<200) & (self.gcm.depth>10)
        #self.mask[:,:250] = False
        #self.mask[:160,:] = False

    def export(self,filename,type='csv'):
        np.savetxt(filename,co.conmat,fmt="%f",delimiter=',')
        
    

class Rectangles(trm.Trm):
    """ Class to  calculate residence times in rectangles.

    """
    def __init__(self,projname,casename="", datadir="", datafile="",
                 ormdir="", griddir="",radius=2):
        super(Rectangles,self).__init__(projname, casename,
                                    datadir, datafile, ormdir)
        self.radius = radius
        self.add_default_regmask()
                
    def generate_regions(self, di=20, dj=20, mask=[]):
        """Create region matrix defining the regions used for connectiv."""
        if len(mask)==0: mask = self.llat>-9999
        ncnt = 1
        self.regmat = self.llon * 0
        regi = np.arange(0, self.imt, di)
        regj = np.arange(0, self.jmt, dj)
        for n,(i,j) in enumerate([(i,j) for i in regi for j in regj]):
            if (mask[j:j+dj, i:i+di]).any():
                self.regmat[j:j+dj, i:i+di] = ncnt
                ncnt += 1
        self.regmat[~mask] = 0
        self.nreg = ncnt - 1




    

