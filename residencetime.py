import numpy ,time
from pycdf import *
from pylab import *
from mpl_toolkits.basemap import Basemap


gr    = CDF('/projData/or5/topo/mask.nc')
tmask = gr.var('tmask')
lat=gr.var('nav_lat')
lon=gr.var('nav_lon')

minLat =  -20
maxLat =  -10
minLon = -125
maxLon = -105

a=numpy.greater(lat[:,:] ,minLat) * numpy.greater(maxLat ,lat[:,:])
b=numpy.greater(lat[:,:] ,minLon) * numpy.greater(maxLon ,lon[:,:])

seedMask=a*b

iGrid ,jGrid = numpy.meshgrid(range(1,511), range(1,722))

minI = iGrid[:,:,seedMask].min()
maxI = iGrid[:,:,seedMask].max()

minJ = jGrid[:,:,seedMask].min()
maxJ = jGrid[:,:,seedMask].max()

trajFile = "ormOut/KAB042j_run.bin"

fd   = open(trajFile, mode='rb')
traj = numpy.fromfile(file=fd, dtype='>i4', count=-1)

ntrac = traj[0::5]
ints  = traj[1::5]
x1    = traj[2::5]
y1    = traj[3::5]
z1    = traj[4::5]
x1.dtype='>f4'
y1.dtype='>f4'
z1.dtype='>f4'

partMask  = ints
partInBox = []

for t in range(2,partMask.max(),1):
   xVec=x1[:,partMask.T==t] 
   yVec=y1[:,partMask.T==t] 
   zVec=z1[:,partMask.T==t] 
   scatter(xVec ,yVec ,s=20 ,c=zVec ,linewidths=0)
   print t
 
   xMask = numpy.greater(xVec ,minI) * numpy.greater(maxI ,xVec)
   yMask = numpy.greater(yVec ,minJ) * numpy.greater(maxJ ,yVec)

   partInBox.append((xMask*yMask).sum())
   print partInBox[-1]


   contour(tmask[0,0,:,:])
   show()
   time.sleep(1)
   clf()
