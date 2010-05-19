import datetime

import numpy as np
import pylab as pl
import matplotlib as mpl
import MySQLdb

class trm:
    """ main class for TRACMASS data manipulation"""
    def __init__(self,projname,casename="",datadir="/Users/bror/ormOut/", 
                 ormdir="Users/bror/svn/orm"):
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

    def create_indexes(self):
        """Create necessary indexes for the traj table."""
        indexlist=[("allints" ,"intstart,ints,ntrac"),("ints" ,"ints"),
                   ("ntrac" ,"ntrac")]
        for i in indexlist:
            sql = ("ALTER TABLE %s ADD INDEX %s (%s);" % 
                   (self.tablename ,i[0],i[1]) )
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

    def ints_to_iso(self,ints):
        base_iso = mpl.dates.date2num(datetime.datetime(2004,1,1))
        

    def sat_to_mysql(self,intstart,field):
        """Load field data to a table. """
        def create_table():
            """Create a mysql table  table """
            itp = " MEDIUMINT NOT NULL "
            ftp = " FLOAT NOT NULL "
            CT  = ( """CREATE TABLE IF NOT EXISTS %s 
                       (ints %s,ntrac %s ,x %s 
                        ,INDEX allints (ints,ntrac) )"""
                    % (self.tablename + field, itp, itp, ftp) )
            self.c.execute(CT)

        create_table()
        

#mpl.dates.num2date(jd30+jd0-1)[0]


def test_insert():
    t = trm('gompom')
    #t.load_to_mysql(33)
    #t.load_to_mysql(563)
    t.sat_to_mysql(33,'chlor_a')

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

    


