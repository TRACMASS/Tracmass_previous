import sys, os
import ConfigParser

import numpy as np
import pylab as pl
import matplotlib as mpl

import anim
from hitta import GrGr
import figpref

import lldist
from postgresql import DB

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

class Traj(object):
    """Main class for trajectory post-processing"""
    def __init__(self, projname, casename=None, region=None):
        """ Setup variables and make sure everything needed is there """
        self.projname = projname
        if casename is None:
            self.casename = projname
        else:
            self.casename = casename
        self.region = region
        self.setup_njord()

    def setup_njord(self):
        cfg = ConfigParser.ConfigParser()
        cfg.read(PROJECT_ROOT + "/projects.cfg")
        if not self.projname in cfg.sections():
            print "Project not included in config file"
            print "Geo-reference functionality will be limited"
            return
        nmod = cfg.get(self.projname, 'njord_module')
        ncls = cfg.get(self.projname, 'njord_class')
        if self.region is None:
            self.region = cfg.get(self.projname, 'map_region')

        #try:
        self.gcm = (__import__("njord." + nmod,
                               fromlist=["njord"]).__dict__[ncls])()
        self.gcm.add_landmask()
        self.landmask = self.gcm.landmask
        self.llon = self.gcm.llon
        self.llat = self.gcm.llat
        #except:
        #print "Couldn't load the NJORD module." 
        #print "Geo-reference functionality will be limited"""
            
    def trajsloaded(aFunc):
        """Decorator function to check if trajs are loaded."""
        def bFunc( *args, **kw ):
            if not "x" in dir(args[0]):
                raise NameError, "Trajectory data not loaded."
            return aFunc( *args, **kw )
        bFunc.__name__ = aFunc.__name__
        bFunc.__doc__ = aFunc.__doc__
        return bFunc

    def add_mp(self, map_region=None):
        if map_region is not None:
            self.region = map_region
            if hasattr(self,'mp'): del self.mp
        try:
            import projmap
        except:
            raise ImportError("Module PROJMAP not available")            
        if not hasattr(self,'mp'):
            self.mp = projmap.Projmap(self.region)
            self.mpxll,self.mpyll = self.mp(self.llon,self.llat)

    @trajsloaded
    def ijll(self,ps=None):
        from scipy.ndimage.interpolation import map_coordinates
        self.lon = map_coordinates(self.llon, [self.y,self.x])
        self.lat = map_coordinates(self.llat, [self.y,self.x])
        self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
        self.lon[self.lon> 180] = self.lon[self.lon> 180] - 360
        self.lon[self.lon==0] = np.nan
        self.lat[self.lat==0] = np.nan

    @trajsloaded
    def ij2utm(self,ps=None):
        from scipy.ndimage.interpolation import map_coordinates
        self.gcm.add_utmxy()
        self.utmx = map_coordinates(self.gcm.utmx, [self.y,self.x])
        self.utmy = map_coordinates(self.gcm.utmy, [self.y,self.x])

    @trajsloaded
    def map(self,weights=[]):
        """Create a field with number of particles in each gridcell"""
        xy = np.vstack( (self.x.astype(np.int), self.y.astype(np.int)) )
        if len(weights) == 0: weights = np.ones(xy.shape[1])
        flat_coord = np.ravel_multi_index(xy, (self.imt, self.jmt))
        sums = np.bincount(flat_coord, weights)
        fld = np.zeros((self.imt, self.jmt))
        fld.flat[:len(sums)] = sums
        return fld.T

    @trajsloaded
    def calc_dists(self):
        if not hasattr(self, 'lon'): self.ijll()
        self.dists = self.lon * np.nan
        self.speed = self.lon * np.nan
        for jd1,jd2 in zip(self.jdvec[:-1], self.jdvec[1:]):
            mask1 = self.jd == jd1
            mask2 = self.jd == jd2
            print jd1,jd2
            ntracmax = max(self.ntrac[mask1].max(), self.ntrac[mask2].max())
            lonvecs = np.zeros((2, ntracmax+1)) * np.nan
            latvecs = np.zeros((2, ntracmax+1)) * np.nan
            lonvecs[0,self.ntrac[mask1]] = self.lon[mask1]
            lonvecs[1,self.ntrac[mask2]] = self.lon[mask2]
            latvecs[0,self.ntrac[mask1]] = self.lat[mask1]
            latvecs[1,self.ntrac[mask2]] = self.lat[mask2]
            self.dists[mask2] = lldist.ll2dist(lonvecs,latvecs)[self.ntrac[mask2]]
            self.speed[mask2] = self.dists[mask2]/(jd2-jd1)/24/60/60



    @trajsloaded
    def field(self,fieldname):
        t = self.ints - self.ints.min()
        ifloor = np.floor(self.i).astype(int)
        jfloor = np.floor(self.y).astype(int)
        iceil  = np.ceil(self.x).astype(int)
        jceil  = np.ceil(self.y).astype(int)
        b1 = fld[t,ifloor,jfloor]
        b2 = fld[t,iceil,jfloor] - b1
        b3 = fld[t,ifloor,jceil] - b1
        b4 = b1 - fld[t,iceil,jfloor] - fld[t,ifloor,jceil] + fld[t,iceil,jceil]
        x = self.x - ifloor
        y = self.y - jfloor
        self.__dict__[fieldname] = b1 + b2*x + b3*y + b4*x*y

    @trajsloaded
    def interp(self,fieldname):
        pass

    @trajsloaded
    def permanent_mask(self, mask):
        """Delete all masked data"""
        for v in vars(self):
            try:
                if len(self.__dict__[v]) == len(mask):
                    self.__dict__[v] = self.__dict__[v][mask]
            except:
                pass
        self.jdvec = np.unique(self.jd)

    @trajsloaded
    def insert(self, database="traj"):
        """Insert current trajectories into database"""
        DB.insert(db)
  
    def select(self,jd=None, runid=0, ints=0, ntrac=0, database=None):
        """ Retrive trajectories from database """
        if not jd: jd = ints
        if not hasattr(self, 'db'):
            self.db = DB(self.projname, self.casename,database=database)
        res = self.db.select(jd, runid, ints, ntrac)
        if len(res) > 0:
            for n,a in enumerate(['runid','ints','ntrac','x','y','z']):
                self.__dict__[a] = np.array(res[n])

    @trajsloaded
    def scatter(self,mask=None, ntrac=None, jd=None, k1=None, k2=None,
                c="g", clf=True, coord="latlon", land="nice", map_region=None):
        if (not hasattr(self,'mp'))  & (coord=="latlon"):
            self.add_mp(map_region)
        if (not hasattr(self,'lon')) & (coord=="latlon"):
            self.ijll()
        if mask is None:
            mask = self.ntrac==self.ntrac
            if jd is not None:
                mask = mask & (self.jd==jd)
            if ntrac is not None:
                mask = mask & (self.ntrac==ntrac)

        figpref.current()
        if clf: pl.clf()
    
        if coord is "latlon":
            x,y = self.mp(self.lon[mask],self.lat[mask])
            scH = self.mp.scatter(x, y, 5, c)
        else:
            scH = pl.scatter(self.x[mask], self.y[mask], 5, c)
    
        if land is "nice":
            land = self.gcm.mp.nice()
        elif coord is "latlon":
            self.mp.pcolormesh(self.mpxll,
                               self.mpyll,
                               np.ma.masked_equal(self.landmask, False),
                               cmap=GrGr())
        else:
            pl.pcolormesh(np.arange(self.gcm.i1, self.gcm.i2+1),
                          np.arange(self.gcm.j1, self.gcm.j2+1),
                          np.ma.masked_equal(self.landmask, False),
                          cmap=GrGr())
        if (ntrac is not None) & (coord is "latlon"):
            self.mp.plot(x,y,'-w',lw=0.5)

        """ 
        xl,yl = self.mp(
            [self.llon[0,0], self.llon[0,-1], self.llon[-1,-1],
             self.llon[-1,0],self.llon[0,0]],
            [self.llat[0,0], self.llat[0,-1], self.llat[-1,-1],
             self.llat[-1,0], self.llat[0,0]]
             )
        #self.mp.plot(xl,yl,'0.5')
        """
        if jd: pl.title(pl.num2date(jd).strftime("%Y-%m-%d %H:%M"))
     
        return scH
        
    @trajsloaded
    def movie(self,di=10, coord='latlon',land="nice"):
        mv = anim.Movie()
        for jd in self.jdvec:
            print self.jdvec[-1] - jd
            if jd/di == float(jd)/di:
                self.scatter(jd=jd, coord=coord, land=land)
                mv.image()
        mv.video(self.projname+self.casename+"_mov.mp4")

    @trajsloaded
    def double_movie(tr1,tr2,di=10):
        mv = anim.Movie()
        ints = np.intersect1d(tr1.ints,tr2.ints).sort()
        for i in ints:
            if i/di == float(i)/di:
                tr1.scatter(ints=i,c="b",clf=True)
                tr2.scatter(ints=i,c="r",clf=False)
                mv.image()
        mv.video(tr1.projname + tr1.casename + "_" +
                 tr2.projname + tr2.casename + "_mov.mp4")

    @trajsloaded
    def export(self,filename='',filetype="mat"):
        import scipy.io as sio

        if type(self.intstart) == int:
            intstart = self.intstart
        else:
            instart = self.intstart.min()
            
        if not filename:
            filename = ( "%s_%s_%i_%i_%i_%i_%i.%s" %
                         (self.projname,self.casename,
                          intstart,
                          self.ints.min(),self.ints.max(),
                          self.ntrac.min(),self.ntrac.max(),filetype)
                         )
        sio.savemat(filename, {'ints':self.ints,
                               'jd':self.jd,
                               'ntrac':self.ntrac,
                               'x':self.x,

                               'y':self.y})

            
        

 

    trajsloaded = staticmethod(trajsloaded)
