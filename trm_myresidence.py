#!/sw/bin/python2.5
import numpy 
import MySQLdb

projName = raw_input("Which project name? ")
caseName = raw_input("Which case name?    ")

db    = MySQLdb.connect("localhost","root","","partsat")
myCur = db.cursor()

xmin = 1848
xmax = 2200
ymin = 1098
ymax = 1199

intsList  = []
k1List    = []
decayList = []

sql =       "SELECT ints,z1,count(ntrac) "
sql = sql + "  FROM %s%s " % (projName,caseName)
sql = sql + "  WHERE x1>%i AND x1<%i " % (xmin ,xmax)
sql = sql + "    AND y1>%i AND y1<%i " % (ymin ,ymax)
sql = sql + "  GROUP BY ints,floor(z1) ORDER BY ints;"

myCur.execute(sql)
resTup = myCur.fetchall()

for i in range(len(resTup)): 
    intsList.append(resTup[i][0])
    k1List.append(resTup[i][1])
    decayList.append(resTup[i][2])

ints   = numpy.array(intsList)
k1     = numpy.array(k1List)
decay  = numpy.array(decayList)
