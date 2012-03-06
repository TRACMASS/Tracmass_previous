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

import batch
import anim
import namelist as nlt
import lldist
import projmaps
from hitta import GrGr

class trm:
    """ main class for TRACMASS data manipulation"""
    def __init__(self,projname, casename="", datadir="", datafile="", ormdir=""):
        self.projname = projname
        if len(casename) == 0:
            self.casename = projname
        else:
            self.casename = casename
        self.datadir = datadir
        if ormdir:
            self.ormdir = ormdir
        else:
            self.ormdir = os.getenv('ORMDIR')

        self.isobase = datetime.datetime(2004,1,1)
        
        self.conn = psycopg2.connect (host="localhost", database="partsat")
        self.c = self.conn.cursor()
        self.tablename = ("%s%s" % (projname ,casename)).lower()

        self.nlgrid = nlt.parse('/%s/projects/%s/%s_grid.in' %
                            (self.ormdir,self.projname,self.projname))
        self.nlrun = nlt.parse('/%s/projects/%s/%s_run.in' %
                            (self.ormdir,self.projname,self.casename))
        if datadir:
            self.datadir = datadir
        else:
            self.datadir = self.nlrun.outDataDir
        if datafile:
            self.datafile = datafile
        else:
            self.datafile=self.nlrun.outDataFile

        self.base_iso = pl.date2num(dtm(
             self.nlgrid.baseYear,
             self.nlgrid.baseMon,
             self.nlgrid.baseDay))-1
        self.imt = self.nlgrid.IMT
        self.jmt = self.nlgrid.JMT

        if projname == 'oscar':
            import oscar
            self.gcm = oscar.Oscar()
        elif projname=="topaz":
            griddir  = '/projData/TOPAZ/1yr_1d/'
            gridname = '/22450101.ocean_daily.nc' 
            g = pycdf.CDF(griddir + gridname)
            lon = g.var('xu_ocean')[:]
            lat = g.var('yu_ocean')[:]
            #self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
            self.llon,self.llat = np.meshgrid(lon, lat)
        elif projname=="casco":
            import casco
            self.gcm = casco.GCM()
            self.region = "casco"
            self.base_iso = pl.date2num(dtm(2004,1,1))
        elif projname=="gompom":
            n = pycdf.CDF(griddir + 'grid.cdf')
            self.llon = n.var('x')[:]    
            self.llat = n.var('y')[:]
            self.base_iso = pl.date2num(dtm(2004,1,1))
        elif projname=="jplSCB":
            import jpl
            self.gcm = jpl.SCB()
            self.region = "scb"
            self.base_iso = pl.date2num(dtm(2001,1,1))-3./24
        elif projname=="jplNow":
            import jpl
            self.gcm = jpl.NOW()
            self.region = "scb"
        elif projname=="rutgersNWA":
            import rutgers
            self.gcm = rutgers.NWA()
            self.region = "nwa"
            
        if hasattr(self,'gcm'):
            self.gcm.add_landmask()
            self.landmask = self.gcm.landmask
            self.llon = self.gcm.llon
            self.llat = self.gcm.llat
            

    def ijll(self,ps=None):
        from scipy.ndimage.interpolation import map_coordinates
        self.lon = map_coordinates(self.llon, [self.y,self.x])
        self.lat = map_coordinates(self.llat, [self.y,self.x])
        self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
        self.lon[self.lon> 180] = self.lon[self.lon> 180] - 360
        self.lon[self.lon==0] = np.nan
        self.lat[self.lat==0] = np.nan

    def ij2utm(self,ps=None):
        from scipy.ndimage.interpolation import map_coordinates
        self.gcm.add_utmxy()
        self.utmx = map_coordinates(self.gcm.utmx, [self.y,self.x])
        self.utmy = map_coordinates(self.gcm.utmy, [self.y,self.x])



    def map(self,weights=[]):
        #xy = np.array(zip(self.x.astype(np.int), self.y.astype(np.int)))
        xy = np.vstack( (self.x.astype(np.int), self.y.astype(np.int)) )
        if len(weights) == 0: weights = np.ones(xy.shape[1])
        flat_coord = np.ravel_multi_index(xy, (self.imt, self.jmt))
        sums = np.bincount(flat_coord, weights)
        fld = np.zeros((self.imt, self.jmt))
        fld.flat[:len(sums)] = sums
        return fld.T

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
                ('ntrac','i4'), ('ints','f8'), 
                ('x','f4'), ('y','f4'), ('z','f4')
                ]))
        return runvec

    def db_drop_all_tables(self, tablename, really=False):
        """Drop master and inherited tables in a partition"""
        if really:
            sql = "SELECT tablename FROM pg_tables WHERE tablename LIKE '%s%%';"
            self.c.execute(sql % tablename)
            tablelist = self.c.fetchall()
            for table in tablelist:
                print table
                self.c.execute("DROP TABLE %s;" % table)
                self.conn.commit()

    def table_exists(self, table_name):
        """ Check if an index exists in the DB """
        sql = "SELECT tablename FROM pg_tables WHERE tablename LIKE '%s';"
        self.c.execute(sql % table_name.lower())
        if self.c.rowcount > 0:
            return True
        else:
            return False

    def db_create_table(self,tablename=None):
        """Create a postgres table for tracmass data"""
        if not tablename: tablename = self.tablename
        if self.table_exists(tablename): return
        CT1 = "CREATE TABLE %s " % tablename
        CT2 = "( runid INT, ints float8, ntrac INT, x REAL ,y REAL, z REAL )"
        self.c.execute(CT1 + CT2)
        self.conn.commit()

    def db_create_partition(self,tablename,partition):
        if self.table_exists(partition): return
        self.db_create_table(tablename)
        sql= "CREATE TABLE %s ( ) INHERITS (%s);"
        self.c.execute(sql % (partition, tablename))
        self.conn.commit()

    def db_create_bulkload_table(self):
        """Create a temporary postgres table for bulkload of trm data """
        if self.table_exists("temp_bulkload"): return
        CT1 = "CREATE TABLE temp_bulkload "
        CT2 = "( ntrac INT, ints float8, x REAL, y REAL, z REAL)"
        self.c.execute(CT1 + CT2)
        self.conn.commit()
        
    def generate_runid(self, jd1=None, jd2=None, temp=False,
                       tablename
                       =None,filename=''):
        """Find or generate row for current run in runs."""
        if not tablename: tablename = self.tablename

        def insert_runid(jd1,jd2):
            sql = ("INSERT INTO runs (jd1,jd2,tablename,filename) " +
                   " values (%s,%s,%s,%s) RETURNING id" )
            self.c.execute(sql, (jd1, jd2, tablename, filename) )
            self.conn.commit()
            return self.c.fetchone()[0]

        if temp:
            return insert_runid(-999, -998)
        else:
            if not jd1: jd1,jd2 = (self.jd.min(),self.jd.max())
            sql = ("SELECT id FROM runs " +
                   " WHERE jd1=%s AND jd2=%s AND tablename='%s'" %
                (jd1, jd2, tablename))
            self.c.execute(sql )
            id = self.c.fetchall()
            if len(id) == 1:
                return id[0][0]
            elif len(id) == 0:
                return insert_runid(jd1,jd2)
            else:
                raise ValueError,"More than one runid in database"

    def index_exists(self, index_name):
        """ Check if an index exists in the DB """
        sql = "SELECT * FROM pg_indexes WHERE indexname LIKE '%s';"
        self.c.execute(sql % index_name.lower())
        if self.c.rowcount > 0:
            return True
        else:
            return False

    def db_add_primary_keys(self, tablename):
        if not self.index_exists("%s_pkey" % tablename):
            sql = "ALTER TABLE %s ADD PRIMARY KEY (runid,ints,ntrac);"
            self.c.execute(sql % ("%s" % tablename) )
            self.conn.commit()

    def db_add_index(self, tablename, indexname, rows):
        index = "%s_%s_idx" % (indexname, tablename)
        if not self.index_exists(index):
            sql = "CREATE INDEX %s ON %s USING btree (%s)"
            self.c.execute(sql % (index, tablename, rows ))
            self.conn.commit()
            
    def db_create_indexes(self):
        """ Create all missing indexes """
        sql = "SELECT distinct(tablename) FROM runs;"
        self.c.execute(sql)
        rowlist = self.c.fetchall()
        for row in rowlist:
            t1 = dtm.now()
            print row[0]
            db_add_primary_keys(row[0])
            print "Time passed: ",dtm.now()-t1
            db_add_index(self, row[0], 'ints', 'ints')
            print "Time passed: ",dtm.now()-t1
            db_add_index(self, row[0], 'runtrac', 'runid,ntrac')
            print "Time passed: ",dtm.now()-t1
            batch.purge()

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
            filename = self.currfile()
         
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

    def tablejds(self,jd):
        return "%s_%i_%i" % (self.tablename, int(jd/10)*10, int(jd/10)*10+10)

    def db_bulkinsert(self,datafile=None):
        """Insert trm bin-files data using pg_bulkload"""
        pg_bulkload = "/opt/local/lib/postgresql90/bin/pg_bulkload"
        ctl_file = "load_trm.ctl"
        db = "-dpartsat"
        outtable = "-O" + "temp_bulkload" # self.tablename

        def run_command(datafile):
            t1 = dtm.now()
            sql = "truncate table temp_bulkload;"
            self.c.execute(sql)
            self.conn.commit()      

            infile = "-i%s/%s" % (self.datadir, datafile)
            spr.call([pg_bulkload,ctl_file,db,infile,outtable])
            print "Elapsed time: " + str(dtm.now()-t1)            

            sql = "SELECT min(ints),max(ints) FROM temp_bulkload;"
            self.c.execute(sql)
            jd1,jd2 = self.c.fetchall()[0]
            tablename = self.tablejds(jd1)
            runid = self.generate_runid(jd1=jd1, jd2=jd2,filename = datafile,
                                        tablename=tablename)
            print "Elapsed time: " + str(dtm.now()-t1)            

         
            self.db_create_partition(self.tablename, tablename)
            sql1 = "INSERT INTO %s (runid,ints,ntrac,x,y,z) " % tablename
            sql2 = "   SELECT %i as runid,ints,ntrac,x,y,z " % runid
            sql3 = "      FROM temp_bulkload;"
            self.c.execute(sql1 + sql2 + sql3)
            self.conn.commit()
            print "Elapsed time: " + str(dtm.now()-t1)

            batch.purge()
            
        if datafile:
            run_command(datafile)
        else:
            flist = glob.glob( "%s/%s*_run.bin" % (self.datadir, self.datafile))
            for f in flist: run_command(os.path.basename(f))

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

    def select(self,jd=None, runid=0, ints=0, ntrac=0):
        """ Retrive trajectories from database """
        if not jd: jd = ints
        whstr = ""
        if runid != 0: whstr += " runid = %i AND" % intstart
        if ints != 0:  whstr += " ints = %i AND" % ints
        if ntrac != 0: whstr += " ntrac = %i " % ntrac
        whstr = whstr.rstrip("AND")
        self.c.execute('SELECT * FROM %s WHERE %s' % (self.tablename,whstr) )
        res = zip(*self.c.fetchall())
        if len(res) > 0:
            for n,a in enumerate(['runid','ints','ntrac','x','y','z']):
                self.__dict__[a] = np.array(res[n])



    def create_seedfile(self,filename, k, mask):
        """Create a seed file based on a 2D mask for TRACMASS """
        ii,jj = np.meshgrid(np.arange(mask.shape[1]),
                            np.arange(mask.shape[0]))
        f = open(filename,'w')
        for i,j,m in zip(np.ravel(ii), np.ravel(jj), np.ravel(mask)):
            if m: f.writelines("% 6i% 6i% 6i% 6i% 6i% 6i\n" %
                               (i+1, j+1, k+1,3,0,50))
        f.close()



    def scatter(self,ntrac=None,ints=None,k1=None,k2=None,c="g",clf=True):
        self.add_mp()
        mask = self.ntrac==self.ntrac
        if ints:
            mask = mask & (self.ints==ints)
        if ntrac:
            mask = mask & (self.ntrac==ntrac)

        if clf: pl.clf()
        self.ijll()
        x,y = self.mp(self.lon[mask],self.lat[mask])

        self.mp.pcolormesh(self.mpxll,self.mpyll,
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

    def currfile(self, ftype='run', stype='bin'):
        flist = glob.glob("%s/%s*%s.%s" %
                          (self.datadir,self.nlrun.outDataFile,
                           ftype,stype))
        datearr = np.array([ os.path.getmtime(f) for f in flist])
        listpos = np.nonzero(datearr == datearr.max())[0][0]
        return os.path.basename(flist[listpos])
            
    
    def ls(self):
        flist = glob.glob("%s/%s*"% (self.datadir,self.nlrun.outDataFile))
        for f in flist: print f

    def add_mp(self):
        if not hasattr(self,'mp'):
            self.mp = projmaps.Projmap(self.region)
            self.mpxll,self.mpyll = self.mp(self.llon,self.llat)


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

"""
  def db_create_loadfilter(self,runid=-999):

        sql = "DROP FUNCTION  trm_bl_filter(int, float8, real, real, real)";
        try:
            self.c.execute(sql )
        except:
            pass
        finally:
            self.conn.commit()

        sql = \"""CREATE FUNCTION trm_bl_filter(int, float8, real, real, real)
                   RETURNS record AS $$
                   SELECT -999, $2, $1, $3, $4, $5
                   $$ LANGUAGE SQL;\"""
        self.c.execute(sql, (runid,) )
        self.conn.commit()
"""
