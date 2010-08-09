import datetime
import glob

import numpy as np
import pylab as pl
import matplotlib as mpl
import MySQLdb

import pycdf
from pyhdf.SD import SD,SDC

class trm:
    """ main class for TRACMASS data manipulation"""
    def __init__(self,projname,casename="",datadir="/Users/bror/ormOut/", 
                 ormdir="/Users/bror/svn/orm"):
        self.projname = projname
        if len(casename) == 0:
            self.casename = projname
        else:
            self.casename = casename
        self.datadir = datadir
        self.ormdir = ormdir
        self.isobase = datetime.datetime(2004,1,1)
        
        conn = MySQLdb.connect (host = "localhost",
                                user = "root", passwd = "",
                                db = "partsat")
        self.c = conn.cursor ()
        self.tablename = "partsat.%s%s" % (projname ,casename)

        griddir  = '/projData/TOPAZ/1yr_1d/'
        gridname = '/22450101.ocean_daily.nc' 
        g = pycdf.CDF(griddir + gridname)
        lon = g.var('xu_ocean')[:]
        lat = g.var('yu_ocean')[:]
        #self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
        self.llon,self.llat = np.meshgrid(lon, lat)
  

    def ijll(self):
        def interp(M):
            ifloor = np.floor(self.y).astype(int)
            jfloor = np.floor(self.x).astype(int)
            iceil  = np.ceil(self.y).astype(int)
            jceil  = np.ceil(self.x).astype(int)
            iceil[iceil==80] = 79
            jceil[jceil==164] = 163
            i1j1 = M[ifloor,jfloor]
            i2j1 = M[iceil, jfloor]
            i1j2 = M[ifloor, jceil]
            
            idf = (i2j1 - i1j1) * (self.y-np.floor(self.y))
            jdf = (i1j2 - i1j1) * (self.x-np.floor(self.x))
            return i1j1 + (idf+jdf)

        self.lon = interp(self.llon)
        self.lat = interp(self.llat)
        self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
        self.lon[self.lon> 180] = self.lon[self.lon> 180] - 360

    def dist(self):
        import lldist

        if not hasattr(self, 'lon'):
            self.ijll()
        self.dist = lldist.lldist(self.lon,self.lat)
        msk = np.zeros([len(self.x)])
        msk[1:] = self.ntrac[1:]-self.ntrac[:-1]
        self.dist[msk != 0] = 0

    def read_bin(self, filename):
        """ Load binary output from TRACMASS """
        fd = open(filename)
        runvec = np.fromfile(fd,np.dtype([
                    ('ntrac','>i4'), ('ints','>i4'), 
                    ('x','>f4'), ('y','>f4'), ('z','>f4')
                    ]))
        return runvec

    def create_table(self):
        """Create a mysql table  table """
        itp = " MEDIUMINT NOT NULL "
        ftp = " FLOAT NOT NULL "
        CT1 = ( "CREATE TABLE IF NOT EXISTS %s (" % self.tablename )
        CT2 = "   intstart " + itp + ",ints " + itp + ",ntrac " + itp
        CT3 = "   ,x " + ftp + " ,y " + ftp + ",z " + ftp
        CT4 = "   )"
        CT  = CT1 + CT2 + CT3 + CT4
        self.c.execute(CT)

    def disable_indexes(self,table):
        self.create_table()
        """Disable indexes to speed up large inserts."""
        sql = ("ALTER TABLE %s DISABLE KEYS;" % table)
        self.c.execute(sql)

    def enable_indexes(self,table):
        """Enable indexes again after a large insert."""
        sql = ("ALTER TABLE %s ENABLE KEYS;" % table)
        self.c.execute(sql)

    def remove_earlier_data_from_table(self, intstart):
        self.c.execute("SELECT DISTINCT(intstart) FROM %s;" %self.tablename)
        if self.c.rowcount > 0:
            DL = ( "DELETE FROM %s WHERE intstart=%s;" % 
                   (self.tablename,intstart) )
            print "Any old posts with intstart=%s deleted." % (intstart)
        else :
            DL = "TRUNCATE  TABLE %s;" % self.tablename
            print "The table %s was truncated." % self.tablename
            self.create_indexes()
        self.c.execute(DL)

    def load(self, intstart=0, ftype="run", stype='bin', filename=''):
        """Load a tracmass output file. Add data to class instance."""
        def unzip(seq): return zip(*seq)
        
        if intstart != 0:
            filename = ("%s%08i_%s.%s" %
                        (self.casename,intstart,ftype,stype) )
        elif filename == '':
            print (self.datadir + self.casename +
                                 "*" + ftype)
            print filename
            #intstart = int(filename[-16:-8])
        if filename[-3:] == "bin":
            runtraj = self.read_bin(self.datadir + filename)
        elif filename[-3:] == "asc":
            print "Not implemented yet."
            raise
        else:
            print "Unknown file format, data file should be bin or asc"
            raise
        tvec = ['ntrac', 'ints', 'x', 'y', 'z']
        self.ntrac, self.ints, self.x, self.y, self.z  = unzip(runtraj)
        for tv in tvec:
            self.__dict__[tv] = np.array(self.__dict__[tv])
        ind = np.lexsort((self.ints,self.ntrac))
        for tv in tvec:
            self.__dict__[tv] = self.__dict__[tv][ind]
 
    def load_to_mysql(self, intstart=0, ftype="run", stype='bin',
                      filename='', db="trm"):        
        """Load a tracmass output file and add to mysql table"""
        self.create_table()
        self.remove_earlier_data_from_table(intstart)
        if intstart != 0:
            filename = "%s%08i_%s.%s" % (self.casename,intstart,ftype,stype)
        else:
            intstart = int(filename[-16:-8])
        if filename[-3:] == "bin":
            runtraj = self.read_bin(self.datadir + filename)
        elif filename[-3:] == "asc":
            print "Not implemented yet."
            raise
        else:
            print "Unknown file format, data file should be .bin or .asc"
            raise
        
        for r in runtraj:
            self.c.execute(
                "INSERT INTO " + self.tablename +
                " (intstart, ntrac, ints, x, y, z) " + 
                " VALUES (-999, %s, %s, %s, %s, %s)", tuple(r) )
        self.c.execute("UPDATE %s SET intstart=%i WHERE intstart=-999" % 
                       (self.tablename, intstart) )

    def ints2iso(self,ints):
        base_iso = mpl.dates.date2num(self.isobase)
        return base_iso + float(ints)/6-1

 
    def trajs(self,intstart=0, ints=0, ntrac=0):
        """ Retrive trajectories from database """
        class traj: pass

        whstr = ""
        if intstart != 0:
            whstr += " intstart = %i AND" % intstart
        if ints != 0:
            whstr += " ints = %i AND" % ints
        if ntrac != 0:
            whstr += " ntrac = %i " % ntrac
        whstr = whstr.rstrip("AND")

        self.c.execute('SELECT * FROM %s WHERE %s' %(self.tablename,whstr) )
        res = zip(*self.c.fetchall())
        if len(res) > 0:
            for n,a in enumerate(['intstart','ints','ntrac','x','y','z']):
                traj.__dict__[a] = np.array(res[n])
        traj = self.ijll(traj)
        return traj

def import_batchrun(batchfile='batch_ints_start.asc'):
    tr = trm('gompom')
    file = tr.ormdir + "/projects/gomoos/" + batchfile
    for t in open(file):
        tr.load_to_mysql(int(t))
    tr._enable_indexes()

#def field_recorder(tr,'fld',datadir=

























































            
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
    





