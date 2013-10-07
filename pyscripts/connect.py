import sys,os
import datetime
import glob
from datetime import datetime as dtm
from itertools import izip
import csv

import numpy as np
import pylab as pl
import tables as td
import scipy.io
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree
from mp_kdtree import mpKDTree
import matplotlib.cm as cm

from hitta import GBRY
import projmaps, anim
import trm
import batch
import figpref
import mycolor

miv = np.ma.masked_invalid

class Matrix(trm.Trm):
    """ Class to handle connectivity matrices from trajectories

    This class generates connectivity matrices from defined regions
    using the output from TRACMASS. The reagions are by default discs
    with a prescribed radius. The discs are packed in a semi-optimal
    fasion withing the part of the grid defined by self.mask.

    Example to calculate a connectivity matrix for a give jd and dt:

    cn = connect.Matrix('rutgersNWA','rutgersNWA')
    cn.load(jdstart=730120)
    cn.calc_conmat()
    """

    
    def __init__(self,projname,casename="", radius=2, **kwargs):
        super(Matrix,self).__init__(projname, casename, **kwargs)
        self.radius = radius
        self.filetype = "hdf"
        
        self.add_default_regmask()
        if not hasattr(self, 'conmatdir'):
            self.conmatdir = os.path.join(os.getcwd(), 'conmatfiles')
        if not os.path.exists(self.conmatdir): os.makedirs(self.conmatdir)
        
                
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
        #self.discKD = cKDTree(list(np.vstack((np.ravel(self.disci),
        #                                      np.ravel(self.discj))).T))
        self.discKD = mpKDTree(self.zip(self.disci, self.discj))
        self.nreg = self.discn.max()+1

    @trm.Traj.trajsloaded
    def regvec_from_discs(self,mask=False, vecmask=False):
        """Generate a vector with region IDs"""
        if not mask:
            self.add_default_regmask()
            mask = self.mask
        if not hasattr(self, 'nreg'): self.generate_regdiscs(mask)
        if vecmask is False:
            xvec = self.x; yvec = self.y
        else:
            xvec = self.x[vecmask]; yvec = self.y[vecmask]

        #dist,ij = self.discKD.query(list(np.vstack((xvec, yvec)).T), 1)
        #self.reg = self.discn[ij]
        #self.reg[dist>self.radius] = 0

        dist,ij = self.discKD.parallel_query(self.zip(xvec, yvec), 1)
        self.reg = self.discn[ij]
        self.reg[dist>self.radius] = 0

    @property
    def dtmax(self):
        if not 'dtmax' in self.__dict__.keys():
            if not hasattr(self, 'x'):
                raise AttributeError, "No particle trajectories loaded."
            self.__dict__['dtmax'] = int(self.jd.max() - self.jd.min()) + 1
        return self.__dict__['dtmax']
    @dtmax.setter
    def dtmax(self,val):
        self.__dict__['dtmax'] = val


    @trm.Traj.trajsloaded
    def calc_conmat(self, jd=None, dt=20, write=False):
        """Create connectivity matrix using a regions matrix."""
        if not jd: jd = self.jd.min()
        if not hasattr(self,'reg'): self.regvec_from_discs()
        tmask1 = self.jd == jd
        tmask2 = self.jd == jd + dt
        try:
            ntracmax = max(self.ntrac[tmask1].max(), self.ntrac[tmask2].max())
        except ValueError:
            return False

        convec = np.zeros((2, ntracmax+1))
        convec[0,self.ntrac[tmask1]] = self.reg[tmask1]
        convec[1,self.ntrac[tmask2]] = self.reg[tmask2]
        convec = convec.astype(np.int)

        flat_coord = np.ravel_multi_index(convec, (self.nreg, self.nreg))
        sums = np.bincount(flat_coord, minlength=self.nreg*self.nreg)
        self.conmat = np.zeros((self.nreg,self.nreg))
        self.conmat.flat[:len(sums)] = sums

        def write_hdf(jd, dt):
            filename = os.path.join(self.conmatdir,"conmat_%s_%s_%06i.h5" %
                                    (self.projname, self.casename, jd))
            shape = (self.dtmax, self.nreg, self.nreg)
            atom = td.UInt32Atom()
            #filters = td.Filters(complevel=5, complib='zlib')
            with td.openFile(filename, 'a') as h5f:
                if hasattr(h5f.root, 'conmat'):
                    ca = h5f.root.conmat
                else:
                    ca = h5f.createCArray(h5f.root, 'conmat', atom, shape)
                ca[dt,:,:] = self.conmat.astype(np.uint32)

        def write_npz(jd, dt):
            conmatfile = ("conmat_%s_%s_%06i_%04i.npz" %
                          (self.projname, self.casename, jd, dt))
            np.savez(os.path.join(self.conmatdir, conmatfile),
                     conmat=self.conmat.astype(np.uint32), jd=jd, dt=dt)

        if write is True:
            write_hdf(jd, dt)


    def __getitem__(self,val):
        if isinstance(val[0], slice):
            jd1 = val[0].start; jd2 = val[0].stop
            if val[0].step: self.djd = val[0].step
        else:
            jd1 = val[0]; jd2 = jd1 + 1
        if isinstance(val[1], slice):
            dt1 = val[1].start; dt2 = val[1].stop
        else:
            dt1 = val[1]; dt2 = dt1 + 1

        for n1,jd in enumerate(np.arange(jd1+1, jd2+1, self.djd)):
            for n2,dt in enumerate(np.arange(dt1,dt2)):
                prefix = "conmat_%s_%s_%i" % (self.projname,self.casename,jd)
                predir = os.path.join(self.conmatdir,prefix)
                if self.filetype is "npz":
                    cmobj = np.load('%s_%04i.npz' % (predir,dt))['conmat']
                else:
                    with td.openFile( "%s.h5" % predir, 'r') as h5f:
                        cmobj = h5f.root.conmat[dt,:,:]            
                try:
                    conmat += cmobj
                except UnboundLocalError:
                    conmat = cmobj
        #conmat[0,:] = 0
        #conmat[:,0] = 0
        #conmat[conmat==0] = np.nan
        return conmat

    def open_hdfconmatfile(self):
        if not hasattr(self,'reg'): self.regvec_from_discs()
        self.h5filename = os.path.join(self.conmatdir,"conmat_%s_%s.h5" %
                                       (self.projname, self.casename))
        jdvec = int((self.jdmax-self.jdmin+1)/self.djd)+1
        shape = (jdvec, self.dtmax, self.nreg, self.nreg)
        iatom = td.UInt32Atom()
        fatom = td.FloatCol()
        batom = td.BoolAtom()
        filtr = td.Filters(complevel=5, complib='zlib')
        self.h5f = h5f = td.openFile(self.h5filename, 'a')
        if not hasattr(h5f.root, 'conmat'):
            crc = h5f.createCArray
            cnmat = crc(h5f.root, 'conmat', iatom,  shape, filters=filtr)
            jdvec = crc(h5f.root, 'jdvec',  fatom, (shape[0],))
            exist = crc(h5f.root, 'exist',  batom, (shape[0],shape[1]))
            jdvec[:] = np.arange(self.jdmin, self.jdmax+1, self.djd)
            exist[:] = False
        else:
            cnmat = h5f.root.conmat
            jdvec = h5f.root.jdvec
            exist = h5f.root.exist
        return cnmat, jdvec, exist
        """
        hf =  td.openFile('test.h5', 'a')
        for n,jd in enumerate(arange(tr.jdmin,tr.jdmax+1,tr.djd)):
            for dt in arange(120):
                try:
                    c[n,dt,:,:] = tr[jd,dt,:]
                    e[n,dt] = True
                except IOError:
                    print "Nope"
                print jd,dt
        """
    
    @trm.Traj.trajsloaded
    def multiplot(self,jd1=730120.0, djd=60, dt=20):

        if not hasattr(self,'disci'):
            self.generate_regdiscs()
            self.x = self.disci
            self.y = self.discj
        if not hasattr(self,'lon'):
            self.ijll()

        figpref.presentation()
        pl.close(1)
        pl.figure(1,(10,10))

        conmat = self[jd1-730120.0:jd1-730120.0+60, dt:dt+10]
        x,y = self.gcm.mp(self.lon, self.lat)
        self.gcm.mp.merid = []
        self.gcm.mp.paral = []

        pl.subplots_adjust(wspace=0,hspace=0,top=0.95)

        pl.subplot(2,2,1)
        pl.pcolormesh(miv(conmat),cmap=cm.hot)
        pl.clim(0,250)
        pl.plot([0,800],[0,800],'g',lw=2)
        pl.gca().set_aspect(1)
        pl.setp(pl.gca(),yticklabels=[])
        pl.setp(pl.gca(),xticklabels=[])
        pl.colorbar(aspect=40,orientation='horizontal',
                    pad=0,shrink=.8,fraction=0.05,ticks=[0,50,100,150,200])

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

        mycolor.freecbar([0.2,.06,0.6,0.020],[2000,4000,6000,8000])

        pl.suptitle("Trajectories seeded from %s to %s, Duration: %i-%i days" %
                    (pl.num2date(jd1).strftime("%Y-%m-%d"),
                     pl.num2date(jd1+djd).strftime("%Y-%m-%d"), dt,dt+10))

        pl.savefig('multplot_%i_%03i.png' % (jd1,dt),transparent=True)
    def all_multiplots(self):
        for jd in np.arange(0,235,60):
            for dt in [10,20,40,60,90]:
                self.multiplot(730120+jd,dt=dt)


    def add_default_regmask(self):
        self.mask = (self.gcm.depth<200) & (self.gcm.depth>10)
        self.mask[:,:250] = False
        self.mask[:160,:] = False

    def export(self,filename,type='csv'):
        np.savetxt(filename,co.conmat,fmt="%f",delimiter=',')
        
    
def ncfile(co):
    nc = Netcdf()
  
    nc.write_conmat(co.conmat,0,0)
    nc.close()


from scipy.io import netcdf

class Netcdf(object):
    """Class to create and populate necdf files for connectivity mats"""

    def __init__(self):
        nc.create_file('test.cdf')
        nc.create_jdvar()
        nc.create_dtvar(np.arange(1,120))
        nc.create_regions(co.discn,co.disci,co.discj)
        nc.create_conmat()

    def create_file(self,filename):
        self.f = netcdf.netcdf_file(filename, 'w')
        self.f.history = 'Connectivity matrices for NWA'

    def create_jdvar(self):
        self.f.createDimension('jd', None)
        self.jdvec = self.f.createVariable('seed_time', 'i', ('jd',))
        self.jdvec.units = 'Julian days from 0001-01-01 (scipy)'

    def create_dtvar(self, pldvec):
        self.f.createDimension('dt', len(pldvec))
        self.dtvec = self.f.createVariable('dtvec', 'i', ('dt',))
        self.dtvec.units = 'Time since start of trajecories (days)'
        self.dtvec[:] = pldvec

    def create_regions(self,regid, regi, regj):
        self.f.createDimension('reg', len(regid)+1)
        ncregid = self.f.createVariable('regid', 'i', ('reg',))
        ncregid.units = 'ID for the different regions.'
        ncregx = self.f.createVariable('regx', 'f', ('reg',))
        ncregx.units = 'X-pos of the region centers'
        ncregy = self.f.createVariable('regy', 'f', ('reg',))
        ncregy.units = 'Y-pos of the region centers.'
        ncregid[1:] = regid
        ncregx[1:] =  regi
        ncregy[1:] =  regj

    def create_conmat(self):
        self.conmat = self.f.createVariable('conmat', 'i',
                                       ('jd','dt','reg','reg'))
        self.conmat.units = 'Connectivity matrix (number of particles).'


    def write_conmat(self,conmat,jdpos=None,dtpos=None):
        self.conmat[jdpos,dtpos,:,:] = conmat

    
    def close(self):
        self.f.close()

def rsquared(dt, jd=0):
    mask = ~np.isnan(np.ravel(co[jd,dt]))
    return linregress(ravel(imat)[mask], ravel(jmat)[mask])[2]        
    

def npz_to_csv(datadir):
    """Convert conmat npz files to csv"""
    
    files = glob.glob(datadir + '/conmat*.npz')
    for fname in files:
        print fname
        conmat = np.load(fname)['conmat']
        fH = open (fname[:-3] + "csv", 'w')
        csw = csv.writer(fH)
        for row in conmat:
            csw.writerow(row)
        fH.close()
    

