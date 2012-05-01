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
import figpref

miv = np.ma.masked_invalid

class Matrix(trm):

    def __init__(self,projname,casename="", datadir="", datafile="",
                 ormdir="", griddir="",radius=2):
        trm.__init__(self,projname,casename,datadir,datafile,ormdir)
        self.radius = radius
        self.add_default_regmask()
                
    def generate_regmat(self, di=20, dj=20, mask=[]):
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

    def generate_regdiscs(self, mask=[]):
        """Create discs defining the regions used for connectivities"""
        if len(mask)==0: mask = self.mask
        ncnt = 1
        nvec = []
        ivec = []
        jvec = []
        r = self.radius
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

    def conmat_from_regmat(self,jd=None,dt=20):
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
    def reg_from_discs(self,mask):
        """Generate a vector with region IDs"""
        self.generate_regdiscs(mask)
        dist,ij = self.discKD.query(list(np.vstack(
            (self.x,self.y)).T),1)
        ij[dist>self.radius]=0
        self.reg = self.discn[ij]
        
    @trajsloaded
    def calc_conmat(self,jd=None,dt=20):
        """Create connectivity matrix using a regions matrix."""
        if not jd: jd = self.jd.min()
        tmask1 = self.jd == jd
        tmask2 = self.jd == jd + dt
        ntracmax = max(self.ntrac[tmask1].max(),
                       self.ntrac[tmask2].max())

        convec = np.zeros((2, ntracmax+1))
        convec[0,self.ntrac[tmask1]] = self.reg[tmask1]
        convec[1,self.ntrac[tmask2]] = self.reg[tmask2]
        convec = convec.astype(np.int)

        weights = np.ones(convec.shape[1])
        flat_coord = np.ravel_multi_index(convec, (self.nreg, self.nreg))
        sums = np.bincount(flat_coord, weights)
        self.conmat = np.zeros((self.nreg,self.nreg))
        self.conmat.flat[:len(sums)] = sums

    def __getitem__(self,val):
        djd = 5
        if isinstance(val[0], slice):
            jd1 = val[0].start; jd2 = val[0].stop
            if val[0].step: djd = val[0].step
        else:
            jd1 = val[0]; jd2 = jd1 + 1
        if isinstance(val[1], slice):
            dt1 = val[1].start; dt2 = val[1].stop
        else:
            dt1 = val[1]; dt2 = dt1 + 1

        for n1,jd in enumerate(np.arange(jd1, jd2, djd)):
            for n2,dt in enumerate(np.arange(dt1,dt2)):
                cmobj = np.load('conmatfiles/conmat_%04i_%04i.npz' %
                                (jd,dt))
                try:
                    conmat += cmobj['conmat']
                except UnboundLocalError:
                    conmat = cmobj['conmat']
                cmobj.close()
        conmat[0:2,:] = 0
        conmat[:,0:2] = 0
        conmat[conmat==0] = np.nan
        return conmat


    def multiplot(self,jd1=730120.0, djd=60, dt=20):

        if not hasattr(self,'disci'):
            self.generate_regdiscs()
            self.x = self.disci
            self.y = self.discj
            self.ijll()

        figpref.manuscript()
        pl.close(1)
        pl.figure(1,(10,10))

        conmat = self[jd1-730120.0:jd1-730120.0+60, dt:dt+10]
        x,y = self.gcm.mp(self.lon, self.lat)
        self.gcm.mp.merid = []
        self.gcm.mp.paral = []

        pl.subplots_adjust(wspace=0,hspace=0,top=0.95)

        pl.subplot(2,2,1)
        pl.pcolormesh(miv(conmat))
        pl.clim(0,50)
        pl.plot([0,800],[0,800],'g',lw=2)
        pl.gca().set_aspect(1)
        pl.setp(pl.gca(),yticklabels=[])
        pl.setp(pl.gca(),xticklabels=[])

        pl.subplot(2,2,2)
        colorvec = (np.nansum(conmat,axis=1)-np.nansum(conmat,axis=0))[1:]
        self.gcm.mp.scatter(x, y, 10, 'w', edgecolor='k')
        self.gcm.mp.scatter(x, y, 10, colorvec)
        self.gcm.mp.nice()
        pl.clim(0,10000)

        
        pl.subplot(2,2,3)
        colorvec = np.nansum(conmat,axis=1)[1:]
        self.gcm.mp.scatter(x, y, 10, 'w', edgecolor='k')
        self.gcm.mp.scatter(x, y, 10, colorvec)
        self.gcm.mp.nice()
        pl.clim(0,10000)

        pl.subplot(2,2,4)
        colorvec = np.nansum(conmat,axis=0)[1:]
        self.gcm.mp.scatter(x, y, 10, 'w', edgecolor='k')
        self.gcm.mp.scatter(x, y, 10, colorvec)
        self.gcm.mp.nice()
        pl.clim(0,10000)

        pl.suptitle("Trajectories seeded from %s to %s, Duration: %i-%i days" %
                    (pl.num2date(jd1).strftime("%Y-%m-%d"),
                     pl.num2date(jd1+djd).strftime("%Y-%m-%d"), dt,dt+10))

        pl.savefig('multplot_%i_%03i.png' % (jd1,dt))

    def add_default_regmask(self):
        self.mask = (self.gcm.depth<200) & (self.gcm.depth>10)
        self.mask[:,:250] = False
        self.mask[:160,:] = False


def all_conmats(mask,jd1=730120.0,jd2=730360.0, djd=5):
    """Generate  nonmats for all trm runs and dt's"""
    co = Matrix('rutgersNWA')
    for jd in np.arange(jd1,jd2,djd):
        co.load(jd)
        co.reg_from_discs(mask)
        for dt in np.arange(1,120):
            co.calc_conmat(dt=dt)
            np.load('conmatfiles/conmat_%04i_%04i.npz' % (jd-jd1,dt),
                     conmat=co.conmat)
            print jd-jd1,dt
        















def generate_netcdf(co, dtlen=120,jd1=730120.0,jd2=730360.0, djd=5):
    from scipy.io import netcdf
    f = netcdf.netcdf_file('conmat.cdf', 'w')
    f.history = 'Connectivity matrices for NWA'

    f.createDimension('jd', None)
    jdvec = f.createVariable('jd', 'i', ('jd',))
    jdvec.units = 'Julian days from 0001-01-01 as defined by scipy'

    f.createDimension('dt', dtlen-1)
    dtvec = f.createVariable('dtvec', 'i', ('dt',))
    dtvec.units = 'Time since start of trajecories (days)'
    dtvec[:] = np.arange(1,dtlen)

    f.createDimension('reg', co.nreg)
    regid = f.createVariable('regid', 'i', ('reg',))
    regid.units = 'ID for the different regions.'
    regx = f.createVariable('regx', 'f', ('reg',))
    regx.units = 'X-pos of the different regions.'
    regy = f.createVariable('regy', 'f', ('reg',))
    regy.units = 'Y-pos of the different regions.'
    regid[1:] = co.discn
    regx[1:] = co.disci
    regy[1:] = co.discj

    conmat = f.createVariable('conmat', 'i', ('jd','dt','reg','reg'))
    conmat.units = 'Connectivity matrix (number of particles).'


    for n1,jd in enumerate(np.arange(jd1,jd2,djd)):
        for n2,dt in enumerate(np.arange(2,3)):
            cmobj = np.load('conmatfiles/conmat_%04i_%04i.npz' % (jd-jd1,dt))
            conmat[n1,n2,:,:] = cmobj['conmat']
            print jd-jd1,dt
            cmobj.close()
        jdvec[n1] = jd
        return
    
    f.close()

def rsquared(dt, jd=0):
    mask = ~np.isnan(np.ravel(co[jd,dt]))
    return linregress(ravel(imat)[mask], ravel(jmat)[mask])[2]        
    
