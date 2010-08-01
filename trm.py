import datetime

import numpy as np
import pylab as pl
import matplotlib as mpl
import MySQLdb

import pycdf
from pyhdf.SD import SD,SDC

class trm:
    """ main class for TRACMASS data manipulation"""
    def __init__(self,projname,casename="",datadir="/Users/bror/ormOut/", 
                 ormdir="/Users/bror/svn/orm",
                 griddir='/projData/GOMPOM/'):
        self.projname = projname
        if len(casename) == 0:
            self.casename = projname
        else:
            self.casename = casename
        self.datadir = datadir
        self.ormdir = ormdir
        conn = MySQLdb.connect (host = "localhost",
                                user = "root", passwd = "",
                                db = "partsat")
        self.c = conn.cursor ()
        self.tablename = "partsat.%s%s" % (projname ,casename)
        n = pycdf.CDF(griddir + 'grid.cdf')
        self.llon = n.var('x')[:]    
        self.llat = n.var('y')[:]


    def ijll(self,tS):
        def interp(M):
            ifloor = np.floor(tS.y).astype(int)
            jfloor = np.floor(tS.x).astype(int)
            iceil  = np.ceil(tS.y).astype(int)
            jceil  = np.ceil(tS.x).astype(int)
            iceil[iceil==80] = 79
            jceil[jceil==164] = 163
            i1j1 = M[ifloor,jfloor]
            i2j1 = M[iceil, jfloor]
            i1j2 = M[ifloor, jceil]
            
            idf = (i2j1 - i1j1) * (tS.y-np.floor(tS.y))
            jdf = (i1j2 - i1j1) * (tS.x-np.floor(tS.x))
            return i1j1 + (idf+jdf)
        tS.lon = interp(self.llon)
        tS.lat = interp(self.llat)
        return tS

    def read_bin(self, filename):
        """ Load binary output from TRACMASS """
        fd = open(filename)
        runvec = np.fromfile(fd,np.dtype([
                    ('ints','>i4'), ('ntrac','>i4'), 
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
        base_iso = mpl.dates.date2num(datetime.datetime(2004,1,1))
        return base_iso + float(ints)/6-1

    def sat_field(self,field,ints):
        ds = mpl.dates.num2date(self.ints2iso(ints))
        box8dir = '/projData/JCOOT_sat_Box8/'
        filepre = "A%i%03i" % (ds.year, ds.timetuple().tm_yday)
        print filepre
        try:
            satfile = glob.glob(box8dir + filepre + "*")[0]
        except IndexError:
            print "Satellite file missing"
            raise
        print satfile
        h = SD(satfile, SDC.READ)
        return h.select(field)[:]

    def sat_to_db(self,field,intstart,ints,batch=False):
        """Load field data to a table. """
        def create_table():
            """Create a mysql table  table """
            itp = " MEDIUMINT NOT NULL "
            ftp = " FLOAT NOT NULL "
            CT  = ( """CREATE TABLE IF NOT EXISTS %s 
                       (intstart %s, ints %s, ntrac %s ,val %s 
                        ,PRIMARY KEY (intstart,ints,ntrac) )"""
                    % (self.tablename + field, itp, itp, itp, ftp) )
            self.c.execute(CT)
            self.disable_indexes(self.tablename + field)

        create_table()
        traj = self.trajs(intstart=intstart,ints=ints)
        if not hasattr(traj,"x"):
            return False
        sati,satj = self.gcmij_to_satij(traj)
        val = self.sat_field(field,ints)[satj,sati]
        for p in zip(traj.ints[val>=0], traj.ntrac[val>=0], val[val>=0]):
            self.c.execute(
                "INSERT INTO " + self.tablename + field +
                " (intstart, ints, ntrac, val) VALUES (-999, %s, %s, %s)", p )
            self.c.execute("UPDATE %s SET intstart=%i WHERE intstart=-999" % 
                           (self.tablename + field, intstart) )
        if not batch:
            self.enable_indexes(self.tablename + field)
            
    def gcmij_to_satij(self,traj):
        n = pycdf.CDF('/Users/bror/svn/modtraj/box8_gompom.cdf')
        igompom = n.var('igompom')[:]
        jgompom = n.var('jgompom')[:]
        ibox8 = n.var('ibox8')[:]
        jbox8 = n.var('jbox8')[:]
        
        k = KDTree(zip(igompom,jgompom))
        ann = k.query(zip(traj.x,traj.y),1,1,1,10)
        return ibox8[ann[1]], jbox8[ann[1]]

    def field_jds(self,field):
        table = self.tablename + field
        class ints: pass
        self.c.execute("SELECT intstart, min(ints), max(ints)" +
                       "FROM %s GROUP BY intstart" % table)
        res = zip(*self.c.fetchall())
        if len(res) > 0:
            for n,a in enumerate(['t0','t1','t2']):
                ints.__dict__[a] = np.array(res[n])
        return ints

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

    def sat_trajs(self,intstart,field,pos='start'):
        """Retrive start- and end-values together with x-y pos for trajs""" 
        if pos == 'start':
            t1str = " AND t.ints=t1.ints "; t2str = ""
        else:
            t2str = " AND t.ints=t2.ints "; t1str = ""
        class traj: pass
        table = self.tablename + field
        self.c.execute("SELECT DISTINCT(ints) FROM %s WHERE intstart=%i" %
                       (table,intstart) )
        ints_t1t2 = self.c.fetchall()
        if len(ints_t1t2) < 2:
            return False 
        sql = """
        SELECT t.x as x , t.y as y, t.ntrac as n,
               t1.val as t1, t2.val as t2
           FROM gompom t
           INNER JOIN %s t1 ON
              t.intstart=t1.intstart %s AND t.ntrac=t1.ntrac
           INNER JOIN %s t2 ON
              t.intstart=t2.intstart %s AND t.ntrac=t2.ntrac
           WHERE t.intstart=%i AND t1.ints=%i AND t2.ints=%i
           """ % (table, t1str, table, t2str, intstart,
                  ints_t1t2[0][0], ints_t1t2[1][0])
        self.c.execute(sql)
        res = zip(*self.c.fetchall())
        if len(res) > 0:
            for n,a in enumerate(['x','y','ntrac','t1','t2']):
                traj.__dict__[a] = np.array(res[n])
            traj = self.ijll(traj)
        return traj

    def traj_ncp(self,t):
        class traj: pass
        t = int(t)
        tr_chl  = self.sat_trajs(t, 'chlor_a')
        tr_chl2 = self.sat_trajs(t, 'chlor_a', pos='end')
        if not tr_chl: return False
        tr_k49 = self.sat_trajs(t, 'K_490')
        msk = np.intersect1d(tr_chl.ntrac,tr_k49.ntrac)
        tr_k49.msk = [np.flatnonzero(tr_k49.ntrac==m).item() for m in msk]
        tr_chl.msk = [np.flatnonzero(tr_chl.ntrac==m).item() for m in msk]
        tr_k49.eu1 = float(4.6) / tr_k49.t1
        tr_k49.eu2 = float(4.6) / tr_k49.t2

        jds =self.field_jds('chlor_a')

        traj.pc2 = tr_chl.t2[tr_chl.msk] * tr_k49.eu2[tr_k49.msk]*60.
        traj.pc1 = tr_chl.t1[tr_chl.msk] * tr_k49.eu1[tr_k49.msk]*60.
        dt = (jds.t2[jds.t0==t]/6. - jds.t1[jds.t0==t]/6.).item()

        print dt

        traj.ncp = (traj.pc2-traj.pc1) / dt 
        traj.ntrac = tr_chl.ntrac[tr_chl.msk]
        traj.x1 = tr_chl.x[tr_chl.msk]
        traj.y1 = tr_chl.y[tr_chl.msk]
        traj.x2 = tr_chl2.x[tr_chl.msk]
        traj.y2 = tr_chl2.y[tr_chl.msk]
        traj.lon1 = tr_chl.lon[tr_chl.msk]
        traj.lat1 = tr_chl.lat[tr_chl.msk]
        traj.lon2 = tr_chl2.lon[tr_chl.msk]
        traj.lat2 = tr_chl2.lat[tr_chl.msk]
        return traj

def import_batchrun(batchfile='batch_ints_start.asc'):
    tr = trm('gompom')
    file = tr.ormdir + "/projects/gomoos/" + batchfile
    for t in open(file):
        tr.load_to_mysql(int(t))
    tr._enable_indexes()
    #t.load_to_mysql(563)
    

def batch_sat_to_db(field='chlor_a',batchprefix='batch_ints'):
    tr = trm('gompom')
    file1 = tr.ormdir + "/projects/gomoos/" + batchprefix + "_start.asc"
    file2 = tr.ormdir + "/projects/gomoos/" + batchprefix + "_end.asc"
    for t1,t2 in zip(open(file1),open(file2)):
        print int(t1),int(t2), int(t2)-int(t1) 
        try:
            tr.sat_to_db(field,int(t1),int(t1)+1,batch=True)
        except IOError:
            print "*** File missing! ***"
        try:
            tr.sat_to_db(field,int(t1),int(t2),batch=True)
        except IOError:
            print "*** File missing! ***"
    tr.enable_indexes()

def calc_ncp_time(batchfile='batch_ints_start.asc'): 
    class tS: pass
    tr = trm('gompom')
    tS.ncp = []
    tS.iso = []
    tS.cnt = []
    for t in open(tr.ormdir + "/projects/gomoos/" + batchfile):
        print t
        ncp = tr.traj_ncp(int(t))
        if not ncp: continue
        tS.ncp.append(np.median(ncp))
        tS.iso.append(tr.ints2iso(int(t)))
        tS.cnt.append(len(ncp))
    return tS
            
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
    


def loadtrajs (trajFile):

    fd   = open(trajFile, mode='rb')
    traj = fromfile(file=fd, dtype='>i4', count=-1)

    ntrac = traj[0::5]
    ints  = traj[1::5]
    x1    = traj[2::5]
    y1    = traj[3::5]
    z1    = traj[4::5]
    x1.dtype='>f4'
    y1.dtype='>f4'
    z1.dtype='>f4'

    return ntrac ,ints ,x1 ,y1 ,z1

#def diagnostics (gridName,runName):

    


