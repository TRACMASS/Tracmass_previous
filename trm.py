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

from traj import Traj

import batch
import namelist as nlt
import lldist

class Trm(Traj):
    """Class for TRACMASS specific functionality"""
    def __init__(self,projname, casename=None, region=None,
                 datadir="", datafile="", trmdir=""):
        super(Trm, self).__init__(projname, casename, region)
        if trmdir:
            self.trmdir = trmdir
        else:
            self.trmdir = os.getenv('TRMDIR')
            if self.trmdir is None:
                raise EnvironmentError, """ Trmdir is not set.
                Add TRMDIR=/path/to/tracmass to your local environment
                or specify trmdir when calling Trm."""

        def parse(od,pn,cn,sfx):
            gridfile = '/%s/projects/%s/%s_%s.in' % (od, pn, cn, sfx)
            if not os.path.isfile(gridfile):
                raise IOError("Can't find the file %s" % gridfile)
            return nlt.parse(gridfile)
        self.nlgrid = parse(self.trmdir,self.projname,self.projname,"grid")
        self.nlrun =  parse(self.trmdir,self.projname,self.casename,"run")

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

        """
        elif projname=="topaz":
            griddir  = '/projData/TOPAZ/1yr_1d/'
            gridname = '/22450101.ocean_daily.nc' 
            g = pycdf.CDF(griddir + gridname)
            lon = g.var('xu_ocean')[:]
            lat = g.var('yu_ocean')[:]
            #self.lon[self.lon<-180] = self.lon[self.lon<-180] + 360
            self.llon,self.llat = np.meshgrid(lon, lat)
        elif projname=="gompom":
            n = pycdf.CDF(griddir + 'grid.cdf')
            self.llon = n.var('x')[:]    
            self.llat = n.var('y')[:]
            self.base_iso = pl.date2num(dtm(2004,1,1))
        """

    def read_bin(self, filename):
        """Read binary output from TRACMASS"""
        with open(filename) as fd:
            runvec = np.fromfile(fd,np.dtype([
                ('ntrac','i4'), ('ints','f8'), 
                ('x','f4'), ('y','f4'), ('z','f4')
                ]))
        return runvec

    def load(self, jdstart=0, intstart=0,
             ftype="run", stype='bin', filename='',rt=False):
        """Load a tracmass output file. Add data to class instance."""
        if jdstart != 0:
            ints = (jdstart+self.base_iso) * 24./self.nlgrid.ngcm  + 1
            filename = ("%s%08i_%s.%s" % (self.datafile,ints,ftype,stype))
            self.jd = jdstart
        elif intstart != 0:
            filename = ("%s%08i_%s.%s" % (self.datafile,intstart,
                                          ftype,stype))
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
            runid = self.generate_runid(jd1=jd1, jd2=jd2,
                                        filename=datafile,
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
            flist = glob.glob( "%s/%s*_run.bin" % (self.datadir,
                                                   self.datafile))
            for f in flist: run_command(os.path.basename(f))

    def create_seedfile(self,filename, k, mask):
        """Create a seed file based on a 2D mask for TRACMASS """
        ii,jj = np.meshgrid(np.arange(mask.shape[1]),
                            np.arange(mask.shape[0]))
        f = open(filename,'w')
        for i,j,m in zip(np.ravel(ii), np.ravel(jj), np.ravel(mask)):
            if m: f.writelines("% 6i% 6i% 6i% 6i% 6i% 6i\n" %
                               (i+1, j+1, k+1,3,0,50))
        f.close()

    def currfile(self, ftype='run', stype='bin'):
        flist = glob.glob("%s/%s*%s.%s" %
                          (self.datadir,self.nlrun.outDataFile,
                           ftype,stype))
        datearr = np.array([ os.path.getmtime(f) for f in flist])
        listpos = np.nonzero(datearr == datearr.max())[0][0]
        return os.path.basename(flist[listpos])

    @property
    def ls(self):
        flist = glob.glob("%s/%s*"% (self.datadir,self.nlrun.outDataFile))
        for f in flist: print f

