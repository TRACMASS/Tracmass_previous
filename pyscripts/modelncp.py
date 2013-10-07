import sys, os
import datetime
import glob
from datetime import datetime as dtm
from itertools import izip
import csv

import numpy as np
import pylab as pl
import scipy.io
from scipy.stats import nanmean, linregress
from matplotlib.colors import LogNorm
from matplotlib import dates as mdates
from scipy.spatial import cKDTree
import matplotlib.cm as cm

from hitta import GBRY
import projmaps, anim
import trm
import batch
import figpref
import mycolor

miv = np.ma.masked_invalid
mf = mdates.DateFormatter('%m-%d')

class Discs(trm.Trm):
    def __init__(self,projname,casename="", datadir="", datafile="",
                 ormdir="", griddir="",radius=2):
        super(Discs,self).__init__(projname, casename,
                                    datadir, datafile, ormdir)

    @trm.Traj.trajsloaded
    def lagDfieldDt(self,fieldname='chlo', jdstart=None, tpos=None, dt=1):
        """Calculate change in field between two timesteps along trajs."""
        if not hasattr(self, fieldname):
            self.fld2trajs(fieldname)
        steps = 24 / self.nlgrid.ngcm
        self.jdvec = np.unique(self.jd)
        ntracmax = self.ntrac.max()
        convec = np.zeros((2, ntracmax+1))
        self.sumvec = np.zeros((1, ntracmax+1))
        self.cntvec = np.zeros((1, ntracmax+1))

        if jdstart is not None:
            tpos = np.nonzero(self.jdvec == jdstart)[0][0]

        for tp in np.arange(0,steps,dt):
            tmask1 = self.jd == self.jdvec[tpos+tp]
            tmask2 = self.jd == self.jdvec[tpos+tp+dt]
            convec[0,self.ntrac[tmask1]] = self.__dict__[fieldname][tmask1]
            convec[1,self.ntrac[tmask2]] = self.__dict__[fieldname][tmask2]
            self.sumvec = np.nansum(np.vstack((convec[1,:]-convec[0,:],
                                            self.sumvec)),axis=0)
            self.cntvec = np.nansum(np.vstack(((convec[1,:]-convec[0,:])*0+1,
                                            self.cntvec)),axis=0)
        return (float(steps)/dt) * (self.sumvec/self.cntvec)

    @trm.Traj.trajsloaded
    def find_live_trajs(self, jd=None, tpos=None, dt=1, prec=1e-5):
        """Find live (moving) trajcories between two timesteps."""
        if jd is not None: tpos = np.nonzero(self.jdvec == jd)[0][0]
        ntracmax = self.ntrac.max()
        convec = np.zeros((2, ntracmax+1))
        tmask1 = self.jd == self.jdvec[tpos]
        tmask2 = self.jd == self.jdvec[tpos + dt]
        convec[0,self.ntrac[tmask1]] = self.x[tmask1]  +  self.y[tmask1]
        convec[1,self.ntrac[tmask2]] = self.x[tmask2]  +  self.y[tmask2]
        #return convec
        return np.nonzero(abs(convec[0,:]-convec[1,:]) > prec)[0]


    @trm.Traj.trajsloaded
    def modelncp(self,fieldname='chlo', dt=1):
        """Calculate model ncp for all loaded trajectories"""
        self.jdvec = np.unique(self.jd)
        steps = 24 / self.nlgrid.ngcm
        tposvec =  np.arange(0,len(self.jdvec),steps)
        histmat = np.zeros((100, len(tposvec)))
        hstyvec = np.linspace(np.log(0.001), np.log(10), 101)

        for n,tpos in enumerate(tposvec[:-1]):
            print n,tpos
            vec = self.lagDfieldDt(fieldname=fieldname, tpos=tpos, dt=dt)
            vec = vec[~np.isnan(vec)]
            histmat[:,n],_ = np.histogram(np.log(vec[vec>0]), hstyvec)
        return histmat
      
    @trm.Traj.trajsloaded
    def calc_decaymatrix(self):
        """Calculate the change of particles #'s in all regions over time"""
        self.regvec_from_regions()
        self.decaymat = np.zeros((self.nreg, len(self.jdvec)))
        for n,jd in enumerate(self.jdvec):
            print n,jd
            self.calc_numpart(jd=jd)
            self.decaymat[:,n] = self.numpart


def all_histmats():
    tr = Discs('bem','restime')
    tr.listfiles()
    
    for dt in [1,6,12,18,24]:
        tr.__dict__['totmat%02i' % dt] = np.zeros((len(tr.runfiles),100,10))

    tr.jdstarts = []
    for n,file in enumerate(tr.runfiles):
        tr.load(filename=os.path.basename(file))
        tr.fld2trajs('chlo')
        for dt in [1,6,12,18,24]:
            tr.__dict__['totmat%02i' % dt][n,:,:] = tr.modelncp(dt=dt)
        tr.jdstarts.append(tr.jdvec[0])
    np.savez('totmat.npz', mat01=tr.totmat01, mat06=tr.totmat06,
                           mat12=tr.totmat12, mat18=tr.totmat18,
                           mat24=tr.totmat24, jdstarts=tr.jdstarts)
    return tr

def stitch_totmats(totmat, jdstarts):
    tvec = (jdstarts-jdstarts[0]).astype(np.int)
    outmat = np.zeros((100,30))
    for n,t in enumerate(tvec):
        outmat[:,0+t:10+t] = outmat[:,0+t:10+t] + totmat[n,:,:]
    return outmat

def totmatplot(mat,dt):
    """Plot a histmoeller diagram of totmat"""
    xi = np.linspace(np.log(0.001), np.log(10), 101)
    figpref.current()
    pl.clf()
    pl.contourf(np.arange(30),np.exp((xi[1:]+xi[:-1])/2),
                miv(mat/mat.max(axis=0)[np.newaxis,:]), np.arange(0,1.1,0.125),
                cmap=cm.OrRd)
    pl.clim(0,1)
    pl.setp(pl.gca(),yscale='log')
    cb = pl.colorbar(aspect=40, pad=0.001, ticks=[0,0.25,0.5,0.75,1])
    cb.ax.set_ylabel('Relative number of particles')
    pl.xlim(0,26)
    pl.xlabel("Time (days)")
    pl.ylabel(r"Positive change in Chlorophyll (mg l$^{-1}$ d$^{-1}$)")
    pl.title("Biological production from particle tracking. Dt=%i hour(s)" % dt)



def hr_vs_dy_figs(ntrac=200):

    figpref.current()
    tr = Discs('bem','bem')
    tr.load()
    tr.fld2trajs('chlo')
    tr.ijll()
    mask = tr.ntrac == ntrac
    chl = tr.chlo[mask][:200]
    jd = tr.jd[mask][:200]
    dchl = chl[1:]-chl[:-1]

    #plot_date(jd,chl,'-')
    #plot(jd[1:], cumsum(dchl)+chl[0])
    plot_date(jd[1:], dchl,'-')
    plot(jd[1::24], dchl[::24])

    pl.figure(1)
    pl.clf()
    x,y = tr.gcm.mp(tr.lon[mask],tr.lat[mask])
    tr.gcm.mp.scatter(x,y)
    tr.gcm.mp.nice()
    pl.savefig('chlhrdy_map_%i.png' % ntrac)

    pl.figure(2)
    pl.clf()
    pl.subplot(2,1,1)
    pl.plot_date(jd,chl,'-')
    pl.plot(jd[1::24], (cumsum(dchl)+chl[0])[::24],'o-')
    pl.ylabel('Chl (mg C / ml)')

    pl.gca().xaxis.set_major_formatter(mf)

    pl.subplot(2,1,2)
    pl.plot_date(jd[1:], dchl,'-')
    pl.plot(jd[1::24], dchl[::24],'o-')
    pl.gca().xaxis.set_major_formatter(mf)
    pl.ylabel('Change in Chl (mg C / ml)')
    pl.xlabel('Time (days)')
    pl.savefig('chlhrdy_dchl_%i.png' % ntrac)
