import datetime
from datetime import datetime as dtm
import glob
import os
from itertools import izip
import cStringIO
import subprocess as spr

import numpy as np
import pylab as pl
import matplotlib as mpl
#import MySQLdb,oursql
import psycopg2

import pycdf
#from pyhdf.SD import SD,SDC

import anim
import namelist as nlt
import lldist
import projmaps
from hitta import GrGr

class trm:
    """ main class for TRACMASS data manipulation"""
    def __init__(self,projname,casename="",datadir="/Users/bror/ormOut/", 
                 datafile="", ormdir="/Users/bror/git/orm"):
        self.projname = projname
        if len(casename) == 0:
            self.casename = projname
        else:
            self.casename = casename
        self.datadir = datadir
        self.ormdir = ormdir
        self.isobase = datetime.datetime(2004,1,1)
        
        self.conn = psycopg2.connect (host="localhost", database="partsat")
        self.c = self.conn.cursor()
        self.tablename = "%s%s" % (projname ,casename)

        self.nlgrid = nlt.parse('/%s/projects/%s/%s_grid.in' %
                            (self.ormdir,self.projname,self.projname))
        self.nlrun = nlt.parse('/%s/projects/%s/%s_run.in' %
                            (self.ormdir,self.projname,self.casename))
        if datafile:
            self.datafile = datafile
        else:
            self.datafile=self.nlrun.outDataFile

        if projname == 'oscar':
            import oscar
            self.gcm = oscar.Oscar()
            self.llat = self.gcm.llat
            self.llon = self.gcm.llon
            self.imt = 1080
            self.jmt = 480
        elif projname=="topaz":
            griddir  = '/projData/TOPAZ/1yr_1d/'
            gridname = '/22450101.ocean_daily.nc' 
            g = pycdf.CDF(griddir + gridname)
            lon = g.var('xu_ocean')[:]
            lat = g.var('yu_ocean')[:]
            #self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
            self.llon,self.llat = np.meshgrid(lon, lat)
            self.imt = 90
            self.jmt = 360
        elif projname=="casco":
            import casco
            self.gcm = casco.GCM()
            self.llon = self.gcm.llon
            self.llat = self.gcm.llat
            self.imt = 285
            self.jmt = 274
            self.region = "casco"
            self.landmask = self.gcm.get_landmask()
            self.base_iso = pl.date2num(dtm(2004,1,1))
        elif projname=="gompom":
            n = pycdf.CDF(griddir + 'grid.cdf')
            self.llon = n.var('x')[:]    
            self.llat = n.var('y')[:]
            self.base_iso = pl.date2num(dtm(2004,1,1))
        elif projname=="jplSCB":
            import jpl
            self.gcm = jpl.SCB()
            self.gcm.add_landmask()
            self.landmask = self.gcm.landmask
            self.llon = self.gcm.llon
            self.llat = self.gcm.llat
            self.imt = 211
            self.jmt = 111
            self.region = "scb"
            self.base_iso = pl.date2num(dtm(2001,1,1))-3./24
        elif projname=="jplNOW":
            import jpl
            self.gcm = jpl.NOW()
            self.gcm.add_landmask()
            self.landmask = self.gcm.landmask
            self.llon = self.gcm.llon
            self.llat = self.gcm.llat
            self.imt = self.nlgrid.IMT
            self.jmt = self.nlgrid.JMT
            self.region = "scb"
            self.base_iso = pl.date2num(dtm(
                self.nlgrid.baseYear,
                self.nlgrid.baseMon,
                self.nlgrid.baseDay))-1
            

    def ijll(self,ps=None):
        from scipy.ndimage.interpolation import map_coordinates
        self.lon = map_coordinates(self.llon, [self.y,self.x])
        self.lat = map_coordinates(self.llat, [self.y,self.x])
        self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
        self.lon[self.lon> 180] = self.lon[self.lon> 180] - 360
        self.lon[self.lon==0] = np.nan
        self.lat[self.lat==0] = np.nan

    def dist(self):
        if not hasattr(self, 'lon'):
            self.ijll()
        self.dist = lldist.lldist(self.lon,self.lat)
        msk = np.zeros([len(self.x)])
        msk[1:] = self.ntrac[1:]-self.ntrac[:-1]
        self.dist[msk != 0] = 0

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

    def db_addfield(self,tp,fieldname):
        pass

    def read_bin(self, filename):
        """ Read binary output from TRACMASS """
        with open(filename) as fd:
            runvec = np.fromfile(fd,np.dtype([
                ('ntrac','>i4'), ('ints','>f8'), 
                ('x','>f4'), ('y','>f4'), ('z','>f4')
                ]))
        return runvec

    def create_table(self):
        """Create a postgres  table """
        itp = " INT  "
        ftp = " REAL "
        CT1 = ( "CREATE TABLE %s (" % self.tablename )
        CT2 = "   runid " + itp + "DEFAULT -999,ints " + ftp + ",ntrac " + itp
        CT3 = "   ,x " + ftp + " ,y " + ftp + ",z " + ftp
        CT4 = "   )"
        CT  = CT1 + CT2 + CT3 + CT4
        try:
            self.c.execute(CT)
        except:
            pass
        finally:
            self.conn.commit()

        """CREATE FUNCTION trm_bl_filter(int, float8, int, real, real, real) RETURNS record AS $$ SELECT -999, $2, $1, $3, $4, $5 $$ LANGUAGE SQL;
        """



        
    def generate_runid(self):
        """Check if run exists in runs table. If not insert run info.
        Return runid"""
        sql = "SELECT id FROM runs WHERE jd1=%s AND jd2=%s AND tablename=%s"
        self.c.execute(sql, (self.jd.min(),self.jd.max(),self.tablename) )
        id = self.c.fetchall()
        if len(id) == 1:
            return id[0][0]
        elif len(id) == 0:
            sql = ("INSERT INTO runs (jd1,jd2,tablename) " +
                   " values (%s,%s,%s) RETURNING id" )
            self.c.execute(sql, (self.jd.min(),self.jd.max(),self.tablename) )
            self.conn.commit()
            return self.c.fetchone()[0]
        else:
            raise ValueError,"More than one runid in database"
        
    def disable_indexes(self,table):
        self.create_table()
        """Disable indexes to speed up large inserts."""
        sql = ("ALTER TABLE %s DISABLE KEYS;" % table)
        self.c.execute(sql)

    def enable_indexes(self,table):
        """Enable indexes again after a large insert."""
        sql = ("ALTER TABLE %s ENABLE KEYS;" % table)
        self.c.execute(sql)

    def remove_earlier_data_from_table(self, runid):
        self.c.execute("SELECT DISTINCT(runid) FROM %s;" %self.tablename)
        if self.c.rowcount > -2:
            DL = ( "DELETE FROM %s WHERE runid=%s;" % 
                   (self.tablename,runid) )
            print "Any old posts with runid=%s deleted." % (runid)
        else :
            DL = "TRUNCATE  TABLE %s;" % self.tablename
            print "The table %s was truncated." % self.tablename
            #self.create_indexes()
        self.c.execute(DL)
        self.conn.commit()

    def load(self, jdstart=0, intstart=0,
             ftype="run", stype='bin', filename='',rt=False):
        """Load a tracmass output file. Add data to class instance."""
        if jdstart != 0:
            ints = (jdstart+self.base_iso) * 24./self.nlgrid.ngcm  + 1
            filename = ("%s%08i_%s.%s" % (self.datafile,ints,ftype,stype))
            self.jd = jdstart
        elif intstart != 0:
            filename = ("%s%08i_%s.%s" % (self.datafile,intstart,ftype,stype))
        elif filename == '':
            print (self.datadir + self.casename +
                                 "*" + ftype)
        if filename[-3:] == "bin":
            runtraj = self.read_bin(self.datadir + filename)
        elif filename[-3:] == "asc":
            print "Not implemented yet."
            raise
        else:
            print "Unknown file format, data file should be bin or asc"
            raise
        if rt: self.runtraj=runtraj
        tvec = ['ntrac', 'ints', 'x', 'y', 'z']
        for tv in tvec:
            self.__dict__[tv] = runtraj[:][tv]
        #self.ints = self.ints.astype(np.int64)
        self.x = self.x - 1
        self.y = self.y - 1
        self.x[self.x<0] = self.x[self.x<0] + self.imt
        #assert self.x.min() >= 0
        #assert self.y.min() >= 0
        assert self.x.min() <= self.imt
        assert self.y.min() <= self.jmt
        self.intstart = intstart
        if self.nlrun.twritetype == 1:
            self.jd = (self.ints.astype(float)/60/60/24 + self.base_iso)
        elif self.nlrun.twritetype == 2:
            self.jd = (self.ints + self.base_iso)
        else:
            self.jd = (self.ints * self.nlgrid.ngcm/24. +self.base_iso) 
        
    def db_insert(self, ftype="run", stype='bin', filename='', db="trm"):
        """Load a tracmass output file and add to mysql table"""
        self.create_table()
        id = self.generate_runid()
        print id
        self.remove_earlier_data_from_table(id)
        sql = ("INSERT INTO " + self.tablename +
               " (runid, ntrac, ints, x, y, z) " + 
               " values (%s,%s,%s,%s,%s,%s)")
        vals = izip((self.x*0+id).astype('float'), self.ntrac.astype('float'),
                    self.ints.astype('float'),   self.x.astype('float'),
                    self.y.astype('float'),      self.z.astype('float'))
        self.c.executemany(sql,vals)

    def db_bulkinsert(self,datafile=None):
        """Insert trm bin-files data using pg_bulkload"""
        pg_bulkload = "/opt/local/lib/postgresql90/bin/pg_bulkload"
        ctl_file = "load_trm.ctl"
        db = "-dpartsat"
        outtable = "-O" + self.tablename

        def run_command(datafile):
            t1 = dtm.now()
            print datafile
            infile = "-i" + datafile
            spr.call([pg_bulkload,ctl_file,db,infile,outtable])
            t2 = dtm.now()
            print "Elapsed time: " + str(t2-t1)

        if datafile:
            run_command(datafile)
        else:
            flist = glob.glob(self.datadir + "/" +
                              self.projname + self.casename + "*_run.bin")
            for f in flist[:10]: run_command(f)

    def db_copy(self):
        if len(self.x) == 0: return False
        self.create_table()
        id = self.generate_runid()
        self.remove_earlier_data_from_table(id)
        vf = cStringIO.StringIO()
        mat = np.vstack( ((self.jd*0+1).astype(np.int),self.jd,self.ntrac,
                           self.x, self.y, self.z) ).T
        np.savetxt(vf,mat,fmt=('%i %f %i %f %f %f') )
        vf.seek(0)
        self.c.copy_from(vf,self.tablename,sep=' ')
        self.conn.commit()

    def ints2iso(self,ints):
        base_iso = mpl.dates.date2num(self.isobase)
        return base_iso + float(ints)/6-1

    def trajs(self,intstart=0, ints=0, ntrac=0):
        """ Retrive trajectories from database """
        whstr = ""
        if intstart != 0:
            whstr += " intstart = %i AND" % intstart
        if ints != 0:
            whstr += " ints = %i AND" % ints
        if ntrac != 0:
            whstr += " ntrac = %i " % ntrac
        whstr = whstr.rstrip("AND")

        self.c.execute('SELECT * FROM %s WHERE %s' %
                       (self.tablename,whstr) )
        res = zip(*self.c.fetchall())
        if len(res) > 0:
            for n,a in enumerate(['intstart','ints','ntrac','x','y','z']):
                self.__dict__[a] = np.array(res[n])
        #self.ijll()

    def scatter(self,ntrac=None,ints=None,k1=None,k2=None,c="g",clf=True):
        if not hasattr(self,'mp'):
            self.mp = projmaps.Projmap(self.region)
            self.xll,self.yll = self.mp(self.llon,self.llat)
        mask = self.ntrac==self.ntrac
        if ints:
            mask = mask & (self.ints==ints)
        if ntrac:
            mask = mask & (self.ntrac==ntrac)

        if clf: pl.clf()
        self.ijll()
        x,y = self.mp(self.lon[mask],self.lat[mask])

        self.mp.pcolormesh(self.xll,self.yll,
                           np.ma.masked_equal(self.landmask,1),cmap=GrGr())
        xl,yl = self.mp(
            [self.llon[0,0], self.llon[0,-1], self.llon[-1,-1],
             self.llon[-1,0],self.llon[0,0]],
            [self.llat[0,0], self.llat[0,-1], self.llat[-1,-1],
             self.llat[-1,0], self.llat[0,0]]
             )
        self.mp.plot(xl,yl,'0.5')

        if ntrac: self.mp.plot(x,y,'-w',lw=0.5)
        self.mp.scatter(x,y,5,c)
        if ints:
            jd = self.jd[self.ints==ints][0]
            pl.title(pl.num2date(jd).strftime("%Y-%m-%d %H:%M"))
        print len(x)

    def movie(self,di=10):
        mv = anim.Movie()
        ints = np.sort(np.unique(self.ints))
        for i in ints:
            print i-ints[0]
            if i/di == float(i)/di:
                self.scatter(ints=i)
                mv.image()
        mv.video(self.projname+self.casename+"_mov.mp4")

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
    def ls(self):
        flist = glob.glob(self.datadir + "/" +
                          self.projname + self.casename + "*")
        for f in flist: print f



def movie2(tr1,tr2,di=10):
    mv = anim.Movie()
    ints = np.intersect1d(tr1.ints,tr2.ints).sort()
    for i in ints:
        if i/di == float(i)/di:
            tr1.scatter(ints=i,c="b",clf=True)
            tr2.scatter(ints=i,c="r",clf=False)
            mv.image()
    mv.video(tr1.projname + tr1.casename + "_" +
             tr2.projname + tr2.casename + "_mov.mp4")



















































            
def grid2ll(latMat ,lonMat ,xVec ,yVec):
    
    xMax = a=latMat.shape[1]-2
    yMax = a=latMat.shape[0]-2

    yVec[yVec>xMax] = xMax
    xVec[xVec>yMax] = yMax

    xFlr  = xVec.astype('int')
    xFrac = xVec-xFlr

    yFlr  = yVec.astype('int')
    yFrac = yVec-yFlr

    latVec1 = latMat[xFlr,yFlr]*(1-xFrac)   + latMat[xFlr+1,yFlr]*xFrac
    latVec2 = latMat[xFlr,yFlr+1]*(1-xFrac) + latMat[xFlr+1,yFlr+1]*xFrac
    latVec  = latVec1*(1-yFrac)             + latVec2*yFrac
    del latVec1 ,latVec2

    lonVec1 = lonMat[xFlr,yFlr]*(1-xFrac)   + lonMat[xFlr+1,yFlr]*xFrac
    lonVec2 = lonMat[xFlr,yFlr+1]*(1-xFrac) + lonMat[xFlr+1,yFlr+1]*xFrac
    lonVec  = lonVec1*(1-yFrac)             + lonVec2*yFrac

    return latVec ,lonVec

def grid2z(k2zMat ,kVec):
    k2zVec = squeeze(k2zMat)
    kVec   = k2zVec.size-kVec-1;
    gC     = ceil(kVec).astype('int')
    gF     = floor(kVec).astype('int')
    gM     = kVec-gF;
    
    zC     = k2zVec[gC];
    zF     = k2zVec[gF];

    return zC*gM+zF*(1-gM);
    





"""
CREATE TABLE runs (id INT NOT NULL AUTO_INCREMENT, PRIMARY KEY(id),
jd1 FLOAT, jd2 FLOAT, tablename VARCHAR(50));
"""
