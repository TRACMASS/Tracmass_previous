import datetime
import glob

import numpy as np
import pylab as pl
import matplotlib as mpl
import MySQLdb

import pycdf
from pyhdf.SD import SD,SDC

from trm import trm

class traj(trm):

    def __init__(self,projname="gmopom",
                 griddir='/projData/GOMPOM/'):
        trm.__init__(self,projname)
        n = pycdf.CDF(griddir + 'grid.cdf')
        self.llon = n.var('x')[:]    
        self.llat = n.var('y')[:]

    def sat_trajs(self,intstart,field,pos='start'):
        """Retrive start- and end-values together
        with x-y pos for trajs""" 
        if pos == 'start':
            t1str = " AND t.ints=t1.ints "; t2str = ""
        else:
            t2str = " AND t.ints=t2.ints "; t1str = ""
        table = self.tablename + field
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

    
    def sat_field(self,field,ints):
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

    def sat_to_db(self,field,intstart,ints,batch=False,traj=False):
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
        if traj:
            val = traj.val
        else:
            traj = self.trajs(intstart=intstart,ints=ints)
            if not hasattr(traj,"x"):
                return False
            sati,satj = self.gcmij_to_satij(traj)
            val = self.sat_field(field,ints)[satj,sati]
        print intstart
        for p in zip(traj.ints[val>-1e9],traj.ntrac[val>-1e9],val[val>-1e9]):
            self.c.execute(
                "INSERT INTO " + self.tablename + field +
                " (intstart,ints,ntrac,val) VALUES (-999,%s,%s,%s)", p )
        self.c.execute("UPDATE %s SET intstart=%i WHERE intstart=-999;" 
                       % (self.tablename + field, intstart) )
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
def ncp(t, db=False):
        class tr: pass
        t = int(t)
        tr_chl = traj('gompom')
        tr_chl.sat_trajs(t, 'chlor_a')
        tr_chl2 = traj('gompom') 
        tr_chl2.sat_trajs(t, 'chlor_a', pos='end')
        if tr_chl.empty: return False
        tr_k49 = traj('gompom') 
        tr_k49.sat_trajs(t, 'K_490')
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
