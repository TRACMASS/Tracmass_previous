import sys
import datetime
import glob
from datetime import datetime as dtm
from itertools import izip

import numpy as np
import pylab as pl
import scipy.io
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree

from hitta import GBRY
import projmaps, anim
from trm import trm
import batch

miv = np.ma.masked_invalid

class Matrix(trm):

    def __init__(self,projname,casename="", datadir="", datafile="",
                 ormdir="", griddir=""):
        trm.__init__(self,projname,casename,datadir,datafile,ormdir)
                
    def create_regmat(self, di=20, dj=20, mask=[]):
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

    def create_regdiscs(self, r=5, mask=[]):
        """Create discs defining the regions used for connectivities"""
        if len(mask)==0: mask = self.llat>-9999
        ncnt = 1
        nvec = []
        ivec = []
        jvec = []
        self.regmat = self.llon * 0
        regi = np.arange(r, self.imt)
        regj = np.arange(r, self.jmt)
        for n,(i,j) in enumerate([(i,j) for i in regi for j in regj]):
            if (mask[j-r:j+r, i-r:i+r]).all():
                ivec.append(i)
                jvec.append(j)
                nvec.append(ncnt)
                ncnt += 1
        self.disci = np.array(ivec)
        self.discj = np.array(jvec)
    
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

    def create_conmat_from_regmat(self,jd=None,dt=20):
        """Create connectivity matrix using a regions matrix."""
        if not jd: jd = self.jd.min()
        tmask1 = self.jd == jd
        tmask2 = self.jd == jd + dt
        ntracmax = max(self.ntrac[tmask1].max(), self.ntrac[tmask2].max())
        convec = np.zeros((2, ntracmax+1))

        self.reg = self.regmat[self.y.astype(int),self.x.astype(int)]
        convec[0,self.ntrac[tmask1]] = self.reg[tmask1]
        convec[1,self.ntrac[tmask2]] = self.reg[tmask2]
        convec = convec.astype(np.int)

        weights = np.ones(convec.shape[1])
        flat_coord = np.ravel_multi_index(convec, (self.nreg, self.nreg))
        sums = np.bincount(flat_coord, weights)
        self.conmat = np.zeros((self.nreg,self.nreg))
        self.conmat.flat[:len(sums)] = sums
        
    def trajsloaded( aFunc ):
        """Trace entry, exit and exceptions."""
        def bFunc( *args, **kw ):
            if not "x" in dir(args[0]):
                raise NameError, "Trajectory data not loaded."
            return aFunc( *args, **kw )
        bFunc.__name__= aFunc.__name__
        bFunc.__doc__= aFunc.__doc__
        return bFunc

        
    @trajsloaded
    def create_conmat_from_discs(self,jd=None,dt=20, r=5):
        """Create connectivity matrix using a regions matrix."""
        if not jd: jd = self.jd.min()
        tmask1 = self.jd == jd
        tmask2 = self.jd == jd + dt
        ntracmax = max(self.ntrac[tmask1].max(), self.ntrac[tmask2].max())
        convec = np.zeros((2, ntracmax+1))

        dist,ij = self.discKD.query(list(np.vstack(
            (self.x,self.y)).T),1)
        ij[dist>r]=0
        self.reg = self.discn[ij]
        

        convec[0,self.ntrac[tmask1]] = self.reg[tmask1]
        convec[1,self.ntrac[tmask2]] = self.reg[tmask2]
        convec = convec.astype(np.int)

        weights = np.ones(convec.shape[1])
        flat_coord = np.ravel_multi_index(convec, (self.nreg, self.nreg))
        sums = np.bincount(flat_coord, weights)
        self.conmat = np.zeros((self.nreg,self.nreg))
        self.conmat.flat[:len(sums)] = sums
