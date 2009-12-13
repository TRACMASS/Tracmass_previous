#!/sw/bin/python2.5
import numpy 
import _mysql

blockSize = 100000

projName     = raw_input("Which project name? ")
caseName     = raw_input("Which case name?    ")
intStartStr  = raw_input("What is intstart?   ")

intstart     = int(intStartStr)

db  = _mysql.connect("localhost","root","","partsat")
CT1 = "CREATE TABLE IF NOT EXISTS partsat.%s%s (" % (projName ,caseName)
CT2 = "    intstart MEDIUMINT UNSIGNED NOT NULL "
CT3 = "    ,ints    MEDIUMINT UNSIGNED NOT NULL "
CT4 = "    ,ntrac   MEDIUMINT UNSIGNED NOT NULL "
CT5 = "    ,x1 FLOAT NOT NULL ,y1 FLOaT NOT NULL ,z1 FLOAT NOT NULL "
CT6 = "   )"
CT  = CT1 + CT2 + CT3 + CT4 + CT5 + CT6
db.query(CT)

chNumSql = "SELECT DISTINCT(intstart) FROM %s%s;" % (projName ,caseName )
db.query(chNumSql)
chNumRes = db.store_result()
numInts  = len(chNumRes.fetch_row(maxrows=0 ,))
if (numInts>1):
    DL = "DELETE FROM partsat.%s%s WHERE intstart=%s;" % \
        (projName ,caseName ,intstart)
    db.query(DL)
    print "Any old posts with intstart=%s deleted." % (intstart)
else :
    DL = "TRUNCATE  TABLE partsat.%s%s;" % (projName ,caseName)
    db.query(DL)
    print "The table %s%s was truncated." % (projName ,caseName)


trajFile  = "ormOut/%s%s_run.bin" % (projName ,caseName)
print trajFile
fd    = open(trajFile, mode='rb')


while (1 == 1) :
    traj  = numpy.fromfile(file=fd, dtype='>i4', count=100000000)
    ntrac = traj[0::5]
    ints  = traj[1::5]
    x1    = traj[2::5]
    y1    = traj[3::5]
    z1    = traj[4::5]
    x1.dtype='>f4' 
    y1.dtype='>f4' 
    z1.dtype='>f4' 

    rows=list()
    for block in range(0,len(x1),blockSize) :
        print "%i %i" % (intstart ,block)
        tt=0
    
        rowStr=''
        blockEnd = min(block+blockSize-1,len(x1))
        for bl in range(block,blockEnd):
            if (x1[bl]+y1[bl]+z1[bl]==x1[bl]+y1[bl]+z1[bl]):
                rowStr += "(%s,%s,%s,%s,%s,%s)," % (
                    intstart ,ints[bl] ,ntrac[bl] ,x1[bl] ,y1[bl] ,z1[bl])
        rowStr = rowStr.rstrip(',')
        if (len(rowStr)>0):
            insertSql = " %s%s%s %s %s " % (
                "INSERT INTO partsat." ,projName ,caseName
                ,"(intstart ,ints ,ntrac ,x1 ,y1 ,z1) VALUES "
                ,rowStr)
            db.query(insertSql)       
                    

print "Create index for all ints"
indexSql = "ALTER TABLE %s%s ADD INDEX %s (%s);" % (
    projName  ,caseName ,"allints" ,"intstart,ints,ntrac")
db.query(indexSql)
print "Create index for ints"
indexSql = "ALTER TABLE %s%s ADD INDEX %s (%s);" % (
    projName  ,caseName ,"ints" ,"ints")
db.query(indexSql)
print "Create index for ntrac"
indexSql = "ALTER TABLE %s%s ADD INDEX %s (%s);" % (
    projName  ,caseName ,"ntrac" ,"ntrac")
db.query(indexSql)
