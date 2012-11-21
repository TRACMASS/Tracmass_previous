
import psycopg2


class DB(object):
    """Helper class to Traj for database handling"""
    def __init__(self,projname,casename,
                 host="localhost",database="trajs",user="postgres"):
        """ Initialize connection to database """
        self.conn = psycopg2.connect(database=database,host=host,
                                     user=user)
        self.c = self.conn.cursor()
        self.tablename = ("%s%s" % (projname ,casename)).lower()

    def drop_all_tables(self, tablename, really=False):
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
        """ Check if an index already exists in the DB """
        sql = "SELECT tablename FROM pg_tables WHERE tablename LIKE '%s';"
        self.c.execute(sql % table_name.lower())
        if self.c.rowcount > 0:
            return True
        else:
            return False

    def create_table(self,tablename=None):
        """Create a postgres table for tracmass data"""
        if not tablename: tablename = self.tablename
        if self.table_exists(tablename): return
        CT1 = "CREATE TABLE %s " % tablename
        CT2 = "( runid INT, ints float8, ntrac INT, x REAL ,y REAL, z REAL )"
        self.c.execute(CT1 + CT2)
        self.conn.commit()

    def create_partition(self,tablename,partition):
        if self.table_exists(partition): return
        self.create_table(tablename)
        sql= "CREATE TABLE %s ( ) INHERITS (%s);"
        self.c.execute(sql % (partition, tablename))
        self.conn.commit()

    def create_bulkload_table(self):
        """Create a temporary postgres table for bulkload of data """
        if self.table_exists("temp_bulkload"): return
        CT1 = "CREATE TABLE temp_bulkload "
        CT2 = "( ntrac INT, ints float8, x REAL, y REAL, z REAL)"
        self.c.execute(CT1 + CT2)
        self.conn.commit()

    def generate_runid(self, jd1=None, jd2=None, temp=False,
                       tablename=None,filename=''):
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

    def add_primary_keys(self, tablename):
        if not self.index_exists("%s_pkey" % tablename):
            sql = "ALTER TABLE %s ADD PRIMARY KEY (runid,ints,ntrac);"
            self.c.execute(sql % ("%s" % tablename) )
            self.conn.commit()

    def add_index(self, tablename, indexname, rows):
        index = "%s_%s_idx" % (indexname, tablename)
        if not self.index_exists(index):
            sql = "CREATE INDEX %s ON %s USING btree (%s)"
            self.c.execute(sql % (index, tablename, rows ))
            self.conn.commit()
            
    def create_indexes(self):
        """ Create all missing indexes """
        sql = "SELECT distinct(tablename) FROM runs;"
        self.c.execute(sql)
        rowlist = self.c.fetchall()
        for row in rowlist:
            t1 = dtm.now()
            print row[0]
            add_primary_keys(row[0])
            print "Time passed: ",dtm.now()-t1
            add_index(self, row[0], 'ints', 'ints')
            print "Time passed: ",dtm.now()-t1
            add_index(self, row[0], 'runtrac', 'runid,ntrac')
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


    def insert(self, db="trm"):
        """Insert trajectories from Traj into database"""
        self.create_table()
        id = self.generate_runid()
        print id
        self.remove_earlier_data_from_table(id)
        sql = ("INSERT INTO " + self.tablename +
               " (runid, ntrac, ints, x, y, z) " + 
               " values (%s,%s,%s,%s,%s,%s)")
        vals = izip((self.x*0+id).astype('float'),
                    self.ntrac.astype('float'),
                    self.ints.astype('float'),   self.x.astype('float'),
                    self.y.astype('float'),      self.z.astype('float'))
        self.c.executemany(sql,vals)

    def tablejds(self,jd):
        return "%s_%i_%i" % (self.tablename, int(jd/10)*10, int(jd/10)*10+10)

    def bulkinsert(self,datafile=None):
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

         
            self.create_partition(self.tablename, tablename)
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

    def copy(self):
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

    def select(self,jd=None, runid=0, ints=0, ntrac=0):
        """ Retrive trajectories from database """
        if not jd: jd = ints
        whstr = ""
        if runid != 0: whstr += " runid = %i AND" % runid
        if ints != 0:  whstr += " ints = %i AND" % ints
        if ntrac != 0: whstr += " ntrac = %i " % ntrac
        whstr = whstr.rstrip("AND")
        self.c.execute('SELECT * FROM %s WHERE %s' %
                       (self.tablename,whstr) )
        return zip(*self.c.fetchall())
