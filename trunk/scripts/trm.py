
import numpy as np
import pylab as pl
import MySQLdb

class trm:
    """ main class for TRACMASS data manipulation"""
    def __init__(self,casename,datadir="/Users/bror/ormOut/", 
                 ormdir="Users/bror/svn/orm"):
        self.casename = casename
        self.datadir = datadir
        self.ormdir = ormdir

    def read_bin(self, filename):
        """ Load binary output from TRACMASS """
        fd = open(filename)
        runvec = np.fromfile(fd,np.dtype([
                    ('ints','>i4'), ('ntrac','>i4'), 
                    ('x','>f4'), ('y','>f4'), ('z','>f4')
                    ]))
        return runvec

    def load_to_mysql(self, ints0=0, ftype="run", stype='bin', filename='', 
                      db="trm", table="traj"):        
        """Load a tracmass output file and add to mysql table"""
        if ints0 != 0:
            filename = "%s%08i_%s.%s" % (self.casename,ints0,ftype,stype)
        else:
            ints0 = int(filename[-16:-8])
        if filename[-3:] == "bin":
            runtraj = self.read_bin(self.datadir + filename)
        elif filename[-3:] == "asc":
            print "Not implemented yet."
            raise
        else:
            print "Unknown file format, data file should be .bin or .asc"
            raise
        
        conn = MySQLdb.connect (host = "localhost",
                                user = "root",
                                passwd = "",
                                db = "trm")        
        c = conn.cursor ()
        for r in runtraj:
            c.execute(
                """INSERT INTO traj (ints0, ntrac, ints, x, y, z)
                         VALUES (-999, %s, %s, %s, %s, %s)""",
                tuple(r) )
        c.execute("UPDATE traj SET ints0=%i WHERE ints0=-999" % ints0)
        

        """ Sql commands to create database:
        CREATE DATABASE trm;
        CREATE TABLE traj (ints0 INT, ntrac INT, ints INT, 
                           x FLOAT, y FLOAT, z FLOAT);
        """

    def sat_to_trajs(self,ints0,field, db="trm", table="traj"):
        pass
        

#mpl.dates.num2date(jd30+jd0-1)[0]


def test_insert():
    t = trm('gompom')
    #t.load_to_mysql(filename="/Users/bror/ormOut/gompom00000033_run.bin")
    t.load_to_mysql(33)
    t.load_to_mysql(563)


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

    


