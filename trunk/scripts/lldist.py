from numpy import *

pi180 =  pi/180
earth_radius = 6378.137
 

def lldist(lon,lat):

    lat = lat * pi180
    dlon = (lon[1:] - lon[:-1])*pi180 
    dlat = (lat[1:] - lat[:-1])

    a = (sin(dlat/2))**2 + cos(lat[:-1]) * cos(lat[1:]) * (sin(dlon/2))**2
    angles = lon * 0
    angles[1:] = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles
