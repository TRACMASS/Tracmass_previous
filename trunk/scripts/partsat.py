import sys
import datetime
import glob
from datetime import datetime as dtm
from itertools import izip

import numpy as np
import pylab as pl
import matplotlib as mpl
import scipy.io
from scipy.spatial import KDTree

import MySQLdb

import pycdf
from pyhdf.SD import SD,SDC

from trm import trm

sys.path.append('/Users/bror/svn/chlcarb/')
import chl2c

miv = np.ma.masked_invalid

class traj(trm):

    def __init__(self,projname,casename="",
                 datadir="/Users/bror/ormOut/", 
                 datafile="", ormdir="/Users/bror/git/orm",
                 griddir='/projData/GOMPOM/'):
        trm.__init__(self,projname,casename,datadir,datafile,ormdir)
        self.flddict = {'par':('L3',),'chl':('box8',)}

        if projname == 'oscar':
            import pysea.NASA
            self.sat = pysea.NASA.nasa(res='4km',ijarea=(700,1700,2000,4000))
            def calc_jd(ints,intstart):
                return self.base_iso + float(ints)/6-1
        elif projname=="casco":
            self.sat = casco.Sat(res='500m')
            def calc_jd(ints,intstart):
                return (self.base_iso +(ints-(intstart)*10800)/150+intstart/8)
        elif projname=="gompom":
            n = pycdf.CDF('/Users/bror/svn/modtraj/box8_gompom.cdf')
            self.gomi = n.var('igompom')[:]
            self.gomj = n.var('jgompom')[:]
            self.sati = n.var('ibox8')[:]
            self.satj = n.var('jbox8')[:]
        elif projname=="jplSCB":
            import mati
            self.sat = mati.Cal()
        elif projname=="jplNOW":
            import mati
            self.sat = mati.Cal()

    def sat_trajs(self,intstart,field,pos='start'):
        """Retrive start- and end-values together
        with x-y pos for trajs""" 
        if pos == 'start':
            t1str = " AND t.ints=t1.ints "; t2str = ""
        else:
            t2str = " AND t.ints=t2.ints "; t1str = ""
        table = self.tablename + '__' + field
        self.c.execute(
            "SELECT DISTINCT(ints) FROM %s WHERE intstart=%i" %
            (table,intstart) )
        ints_t1t2 = self.c.fetchall()
        if len(ints_t1t2) < 2:
            self.empty = True
            return
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
                self.__dict__[a] = np.array(res[n])
            self.empty = False
        else:
            self.empty = True
        self.ijll()

    def sat_conc(self,intstart,field,pos='start'):
        """Retrive fields of start- and end-values""" 
        jdS = self.field_jds(field)
        if intstart in jdS.t0:
            if pos == 'start':
                ints = jdS.t1[jdS.t0==intstart].item()
            else:
                ints = jdS.t2[jdS.t0==intstart].item()
        else:            
            return self.llat * 0
        table = self.tablename + field
        sql = """
               SELECT round(t.x) x, round(t.y) y, avg(n.val) val
                  FROM gompom t
                  INNER JOIN %s n ON
                     t.intstart=n.intstart AND t.ntrac=n.ntrac
                  WHERE t.intstart=%i AND n.ints=%s AND t.ints=%s
                  GROUP BY round(t.x),round(t.y)
           """ % (table, intstart, intstart+1, ints)
        n = self.c.execute(sql)
        fld = self.llat * 0
        if n == 0: return fld
        x,y,val = zip(*self.c.fetchall())
        fld[np.array(y).astype(int),np.array(x).astype(int)-1] = val
        return fld

    def trajs(self,intstart=0, ints=0, ntrac=0,fld=''):
        """ Retrive trajectories from database """
        whstr = ""
        if intstart != 0:
            whstr += " t.intstart = %i AND" % intstart
        if ints != 0:
            whstr += " t.ints = %i AND" % ints
        if ntrac != 0:
            whstr += " t.ntrac = %i " % ntrac
        whstr = whstr.rstrip("AND")
        if fld:
            valstr = " ,n.val val "
            table2 = self.tablename + fld
            valwhere = " AND n.ints=%s " % ints
            valjoin = (" INNER JOIN %s n ON " % table2 +
                       " t.intstart=n.intstart AND t.ntrac=n.ntrac ")
        else:
            valstr = " ,t.intstart val "
            table2 = valwhere = valjoin = ''
        sql = ("SELECT t.ints ints, t.ntrac ntrac, t.x x, t.y y %s FROM %s t %s WHERE %s %s "
               % (valstr, self.tablename, valjoin, whstr, valwhere) )
        n = self.c.execute(sql)
        res = zip(*self.c.fetchall())
        if len(res) > 0:
               for n,a in enumerate(['ints','ntrac','x','y','val']):
                   self.__dict__[a] = np.array(res[n])
        #self.ijll()
  
    def load_sat(self,field,ints=0,intstart=0,jd=0):
        """ Load satellite field for a given jd or ints"""
        self.sat.load(field,jd=jd)
        return self.sat.__dict__[field]
        """ GOMPOM
            ds = mpl.dates.num2date(self.ints2iso(ints))
            box8dir = '/projData/JCOOT_sat_Box8/'
            filepre = "A%i%03i" % (ds.year, ds.timetuple().tm_yday)
            try:
                satfile = glob.glob(box8dir + filepre + "*")[0]
            except IndexError:
                print "Satellite file missing"
                raise
            h = SD(satfile, SDC.READ)
            return h.select(field)[:]
        """

    def gcmij_to_satij(self,mask=[]):
        if len(mask)==0: mask = tr.x==tr.x
        self.ijll()
        self.si, self.sj = self.sat.ll2ij(self.lon[mask],self.lat[mask])

    def remove_satnans(self,jd=False):
        if not jd: jd = self.jd.min()
        self.gcmij_to_satij(self.jd==jd)
        chl = self.load_sat('chl',jd=jd)
        chl = chl[self.si,self.sj]
        ntrac = self.ntrac[self.jd==jd]
        ntrac = ntrac[~np.isnan(chl)]
        mask = np.in1d(self.ntrac,ntrac)
        self.x = self.x[mask]
        self.y = self.y[mask]
        self.z = self.z[mask]
        self.jd = self.jd[mask]
        self.ntrac = self.ntrac[mask]
        self.ints = self.ints[mask]


    def sat_to_db(self,field,intstart,ints,batch=False,traj=False):
        """Load field data to a table. """
        table = self.tablename + '__' + field
        def create_table():
            """Create a mysql table  table """
            itp = " INT NOT NULL "
            ftp = " FLOAT NOT NULL "
            CT  = ( """CREATE TABLE IF NOT EXISTS %s 
                       (intstart %s, ints %s, ntrac %s ,val %s 
                        ,PRIMARY KEY (intstart,ints,ntrac) )"""
                    % (table, itp, itp, itp, ftp) )
            self.c.execute(CT)
            self.disable_indexes(table)
        create_table()
        if not hasattr(self,'gomi'): self.create_ijvecs(field)
        if traj:
            val = traj.val
        else:
            self.trajs(intstart=intstart,ints=ints)
            if not hasattr(self,"x"):
                return False
            self.gcmij_to_satij()
            self.val = self.load_sat(field,ints)[self.sj,self.si]
            print len(self.val[self.val>-1e9])
            print ("min:", min(self.ints[self.val>-1e9]*0+intstart),
                   min(self.ints[self.val>-1e9]))
        plist = zip(
                    self.ints[self.val>-1e9]*0+intstart,
                    self.ints[self.val>-1e9],
                    self.ntrac[self.val>-1e9],
                    self.val[self.val>-1e9])
        sql = ("INSERT INTO " + table + " (intstart,ints,ntrac,val) " +
               " VALUES (%s,%s,%s,%s)")
        self.c.executemany(sql,plist)
        if not batch:
            self.enable_indexes(table)
            
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

    def median_ncp(self):
        class svec: pass
        self.c.execute(
            """CREATE TEMPORARY TABLE IF NOT EXISTS temp_median 
                    (id INT AUTO_INCREMENT PRIMARY KEY) 
                SELECT ints, val FROM gompomncp 
                ORDER BY ints, val;
             """)
        self.c.execute(
            """CREATE TEMPORARY TABLE temp_median_ids 
               SELECT ROUND(AVG(id)) AS id FROM temp_median 
               GROUP BY ints;
            """)
        self.c.execute(
            """SELECT ints, val FROM temp_median_ids 
               LEFT JOIN temp_median USING (id) ORDER BY ints;
            """)
        res = zip(*self.c.fetchall())
        svec.t1 = res[0]
        svec.medianNCP = res[1]
        self.c.execute(
            """SELECT ints, count(val),max(val),min(val),avg(val) 
                   FROM gompomncp GROUP BY ints
                   ORDER BY ints;
            """)
        res = zip(*self.c.fetchall())
        svec.sum   = res[1]
        svec.max   = res[2]
        svec.min   = res[3]
        svec.mean  = res[4]
        
        return svec


#####################################################################
def ncp(t, db=False,kpar=False):
        class tr: pass
        t = int(t)
        tr_chl = traj('gompom')
        tr_chl.sat_trajs(t, 'chlor_a')
        tr_chl2 = traj('gompom') 
        tr_chl2.sat_trajs(t, 'chlor_a', pos='end')
        if tr_chl.empty: return False
        tr_k49 = traj('gompom')
        tr_k49.sat_trajs(t, 'K_490')
        if kpar:
            """
http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=2997
            """
            tr_k49.t1 = 0.0864 + 0.884*tr_k49.t1 - 0.00137/tr_k49.t1
            tr_k49.t2 = 0.0864 + 0.884*tr_k49.t2 - 0.00137/tr_k49.t2
        msk = np.intersect1d(tr_chl.ntrac,tr_k49.ntrac)
        tr_k49.msk = [np.flatnonzero(tr_k49.ntrac==m).item()
                      for m in msk]
        tr_chl.msk = [np.flatnonzero(tr_chl.ntrac==m).item()
                      for m in msk]
        tr_k49.eu1 = float(4.6) / tr_k49.t1
        tr_k49.eu2 = float(4.6) / tr_k49.t2
        if len(tr_chl.msk) == 0: return False
        jds =tr_chl.field_jds('chlor_a')

        tr.pc2 = tr_chl.t2[tr_chl.msk] * tr_k49.eu2[tr_k49.msk]*60.
        tr.pc1 = tr_chl.t1[tr_chl.msk] * tr_k49.eu1[tr_k49.msk]*60.
        dt = (jds.t2[jds.t0==t]/6. - jds.t1[jds.t0==t]/6.).item()

        tr.ncp = (tr.pc2-tr.pc1) / dt
        tr.ntrac = tr_chl.ntrac[tr_chl.msk]
        tr.x1 = tr_chl.x[tr_chl.msk]
        tr.y1 = tr_chl.y[tr_chl.msk]
        tr.x2 = tr_chl2.x[tr_chl.msk]
        tr.y2 = tr_chl2.y[tr_chl.msk]
        tr.lon1 = tr_chl.lon[tr_chl.msk]
        tr.lat1 = tr_chl.lat[tr_chl.msk]
        tr.lon2 = tr_chl2.lon[tr_chl.msk]
        tr.lat2 = tr_chl2.lat[tr_chl.msk]

        if db:
            tr.ints=tr.x1*0+t
            tr.val = tr.ncp
            tr_chl.sat_to_db('ncp',intstart=t,ints=t,batch=False,traj=tr)
        return tr

def batch_sat_to_db(field='chlor_a',batchprefix='batch_ints'):
    tr = traj('gompom')
    file1 = tr.ormdir + "/projects/gomoos/" + batchprefix + "_start.asc"
    file2 = tr.ormdir + "/projects/gomoos/" + batchprefix + "_end.asc"
    for t1,t2 in zip(open(file1),open(file2)):
        print int(t1), int(t1)+1, int(t1) ,int(t2), int(t2)-int(t1)
        if int(t1) > 4000:
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
    import partsat
    class tS: pass
    tr = partsat.traj('gompom')
    tS.ncp = []
    tS.iso = []
    tS.cnt = []
    for t in open(tr.ormdir + "/projects/gomoos/" + batchfile):
        ncp = tr.traj_ncp(int(t))
        if not ncp: continue
        tS.ncp.append(np.median(ncp.ncp))
        tS.iso.append(tr.ints2iso(int(t)))
        tS.cnt.append(len(ncp.ncp))
    return tS

def test(ps2):
    ps = traj(projname='casco')
    ps.gomi,ps.gomj,ps.sati,ps.satj,ps.sat = (
        ps2.gomi,ps2.gomj,ps2.sati,ps2.satj,ps2.sat)
    print ps.sat_field('chlor_a',68612400,6535)

def add_all(ps):
    intstart = [6353,6361,7665,7865,7873,7993,8001]
    ints = [68612400, 68698800, 68785200, 84250800,
            84942000, 85028400, 86324400, 86410800] 

    for s,i in zip(intstart,ints)[4:]:
        print s,i
        ps.sat_to_db('chlor_a',intstart=s,ints=i,batch=True)
        print s,i+24*3600
        ps.sat_to_db('chlor_a',intstart=s,ints=(i+24*3600),batch=True)



def interp(ps,intstart=6353):
    import figpref
    import projmaps
    figpref.current()
    
    ps.c.execute("SELECT distinct(ints) FROM casco__chlor_a " + 
                 "   WHERE intstart=%i" % intstart) 
    ints1,ints2 = ps.c.fetchall()

    sql = """SELECT t1.ntrac,t1.val,t2.val,p.x,p.y,p.ints
               FROM casco__chlor_a t1
               INNER JOIN casco__chlor_a t2 ON t1.ntrac=t2.ntrac
               INNER JOIN casco p ON t1.ntrac=p.ntrac
               WHERE t1.intstart=%i
                  AND t2.intstart=%i
                  AND p.intstart=%i
                  AND t1.ints=%i
                  AND t2.ints=%i;
                  """ % (intstart,intstart,intstart,ints1[0],ints2[0])
    ps.c.execute(sql)
    res = zip(*ps.c.fetchall())
    class trj: pass
    if len(res) > 0:
        for n,a in enumerate(['ntrac','t1val','t2val','x','y','ints']):
            trj.__dict__[a] = np.array(res[n])

    mask = (trj.ints.astype(np.float)/100-trj.ints/100)>0.5
    trj.ints[mask]=trj.ints[mask]+1
    tvec = np.unique(trj.ints)
    tvec = tvec[tvec<=ints2[0]]
    itvec = tvec-tvec.min()
    itvec = itvec.astype(np.float)/itvec.max()

    """
    for n,t in enumerate(tvec):
        pl.clf()
        pl.scatter(trj.x[trj.ints==t],
                   trj.y[trj.ints==t],5,
                   np.log(trj.t1val[trj.ints==t])*(1-itvec[n]) +
                   np.log(trj.t2val[trj.ints==t])*(itvec[n]) )
        print itvec[n]
        pl.savefig("interp_%i.png" % t, dpi=100)
    """
    mp = projmaps.Projmap('casco')
    xl,yl = mp(ps.llon,ps.llat)
    for n,t in enumerate(tvec):
        fld = ps.cs.llat * 0
        cnt = ps.cs.llat * 0

        xvc = (trj.x[trj.ints==t].astype(np.int))
        yvc = (trj.y[trj.ints==t].astype(np.int))
        val = ( np.log(trj.t1val[trj.ints==t])*(1-itvec[n]) +
                np.log(trj.t2val[trj.ints==t])*(itvec[n]) )
        for x,y,v in zip(xvc,yvc,val):
            fld[y,x] += v
            cnt[y,x] += 1

        pl.clf()
        mp.pcolormesh(xl,yl,miv(fld/cnt))
        mp.nice()
        jd = (ps.base_iso +
              (float(t)-(intstart)/8*3600*24)/3600/24 + intstart/8) + 0.583
        pl.title(pl.num2date(jd).strftime("log(Chl) %Y-%m-%d %H:%M"))
        pl.clim(-2,2.5)
        pl.savefig("interp_%i_%03i.png" % (intstart,n), dpi=100)


def test():
    tr = traj('jplNOW','ftp','/Volumes/keronHD3/ormOut/')
    tr.load(733833.0)
    tr.remove_satnans()
    tr.db_insert()

def batch_insert():
    tr = traj('jplNOW','ftp','/Volumes/keronHD3/ormOut/')
    for jd in np.arange(733773.0, 734138.0):
        print pl.num2date(jd), jd
        tr.load(jd)
        tr.remove_satnans()
        tr.db_insert()

def profile():
    import cProfile

    cProfile.run('test()')
