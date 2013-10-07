
import numpy as np
import pylab as pl
import matplotlib.cm as cm
from scipy.stats import nanmean

import casco
import trm
import figpref
import anim
import stload

miv = np.ma.masked_invalid

def hrvec(tr):
    """Return a vector with hourly time postions in a jdvec"""
    jdvec = np.unique(tr.jd)
    dtvec = jdvec -tr.jdmin
    hours = np.round((dtvec-np.round(dtvec))*24*1000)/1000
    return (jdvec[hours==hours.round()])

####
def tracer_rms(field='salt'):
    """ Make a nrmsd vs dt plot"""

    cs = casco.GCM()
    cs.fatten_landmask()
    cs.fatten_landmask()

    def load_tr(filename):
        tr = trm.trm('casco',"full")
        trl = stload.load(filename)
        tr.ints = trl.ints
        tr.ntrac = trl.ntrac
        tr.x = trl.x
        tr.y = trl.y
        tr.z = trl.z
        tr.jd = trl.jd
        tr.jdmin = tr.jd.min()
        return tr

    def get_field(jd):
        field = tr.traj2grid('z',tr.jd==jd)
        field[cs.landmask==0]=np.nan
        return  np.ravel(field)

    pl.clf()

    lHlist = []
    for file in ['tr_salt_0301.npz','tr_temp_0301.npz',
                 'tr_active_0301.npz']:
        tr = load_tr(file)
        tr.z[cs.landmask[tr.y.astype(int), tr.x.astype(int)]==0] = np.nan
        dtvec = hrvec(tr)-tr.jdmin
        
        msk0 = tr.jd==tr.jdmin
        fld0 = np.ravel(cs.get_field(field,tr.jdmin))
        lagvec = []

        for dt in dtvec:
            diffmat = tr.lagr_diff('z',msk0,tr.jd==(tr.jdmin+dt))
            rmsd,nrmsd = tr.rmsd(diffmat[:,3],diffmat[:,2],False)
            lagvec.append(nrmsd)

        lHlist.append(pl.plot(dtvec*24,np.array(lagvec)*100))

    pl.xlabel('Time (hours)')
    pl.ylabel('NRMSD (%)')
    pl.xlim(0,48)
    pl.ylim(0,14)
    pl.legend(lHlist,('Salinity','Temperature','Reactive tracer'),
              'upper left' )
    pl.title('Temporal decorrelation of tracers')
    #pl.savefig('figs/nrmsd_tracer.png')
    pl.xlim(0,12)
    pl.ylim(0,10)
    #pl.savefig('figs/nrmsd_tracer_detail.png')
    pl.savefig('figs/nrmsd_tracers.png')


###
def rms(field='salt'):
    """ Make a nrmsd vs dt plot"""
    tr = trm.trm('casco',"full")
    cs = casco.GCM()
    cs.fatten_landmask()
    cs.fatten_landmask()
    tr.load(intstart=5858969)
    #tr.load(intstart=5859217)
    tr.z[cs.landmask[tr.y.astype(int), tr.x.astype(int)]==0] = np.nan
    dtvec = hrvec(tr)-tr.jdmin

    def get_field(jd):
        field = tr.traj2grid('z',tr.jd==jd)
        field[cs.landmask==0]=np.nan
        return  np.ravel(field)
        
    msk0 = tr.jd==tr.jdmin
    fld0 = np.ravel(cs.get_field(field,tr.jdmin))
    fld0 = get_field(tr.jdmin)
    lagvec = []
    eulvec = []

    def eulmsk(diffmat):
        eulmsk = cs.landmask*np.nan
        eulmsk[diffmat[:,1].astype(np.int),
               diffmat[:,0].astype(np.int)] = diffmat[:,3]-diffmat[:,2]
        return eulmsk*0+1

 
    
    for dt in dtvec:
        diffmat = tr.lagr_diff('z',msk0,tr.jd==(tr.jdmin+dt))
        rmsd,nrmsd = tr.rmsd(diffmat[:,3],diffmat[:,2],False)
        lagvec.append(nrmsd)
        fld1 = np.ravel(cs.get_field(field,tr.jdmin+dt) * eulmsk(diffmat))
        fld1 = get_field(tr.jdmin+dt)
        rmsd,nrmsd = tr.rmsd(fld0,fld1,False)
        eulvec.append(nrmsd)

    pl.clf()
    pl.plot(dtvec*24,np.array(lagvec)*100)
    pl.plot((dtvec*24)[0:-2:2],(np.array(eulvec)*100)[0:-2:2])
    pl.xlabel('Time (hours)')
    pl.ylabel('NRMSD (%)')
    pl.xlim(0,48)
    pl.ylim(0,14)
    pl.legend(('Lagrangian frame','Eularian frame'),'upper left')
    pl.title('Temporal decorrelation of %s' % field)
    pl.savefig('figs/nrmsd_%s.png' % field)
    pl.xlim(0,12)
    pl.ylim(0,10)
    pl.savefig('figs/nrmsd_%s_detail.png' % field)


def eulrms():
    pl.clf()
    cs = casco.GCM()
    tr = trm.trm('casco',"full")
    tr.load(intstart=5858969)
    dtvec = hrvec(tr)-tr.jdmin

    def get_salt(jd): return np.ravel(tr.traj2grid('z',tr.jd==jd))

    fld0 = get_salt(tr.jdmin)
    rmsvec = []
    for dt in dtvec:
        fld1 = get_salt(tr.jdmin+dt)
        rmsd,nrmsd = tr.rmsd(fld0,fld1)
        rmsvec.append(rmsd)
    pl.clf()
    pl.plot(dtvec*24,np.array(rmsvec)*100,'*-')
    pl.xlabel('Time (hours)')
    pl.ylabel('NRMSD (%)')
    pl.xlim(0,24)
    
def eul(tr):
    pl.clf()
    cs = casco.GCM()
    jdmin = tr.jd.min()
    dtvec = np.arange(0,1.125,0.125)

    def get_salt(jd):
        cs.load('temp',jd=jd)
        return cs.temp
    fld0 = get_salt(jdmin)

    def gen_dvec(dt):
        fld1 = get_salt(jdmin+dt)
        dvec = np.ravel((fld1[0,:,:]-fld0)[0,:,:])
        return dvec[~np.isnan(dvec)]

    dlist = []
    for dt in dtvec:
        dlist.append(np.abs(gen_dvec(dt)))
    pl.boxplot(dlist,sym='',positions=dtvec*24)

def eultraj(tr):
    pl.clf()
    dtvec = hrvec(tr)
    def get_salt(jd): return tr.traj2grid('z',tr.jd==jd)

    fld0 = get_salt(jdmin)
    def gen_dvec(dt):
        fld1 = get_salt(jdmin+dt)
        dvec = np.ravel(fld1-fld0)
        return dvec[~np.isnan(dvec)]

    dlist = []
    for dt in dtvec:
        dlist.append(np.abs(gen_dvec(dt)))
    pl.boxplot(dlist,sym='',positions=np.round(dtvec*24).astype(np.int))
    pl.xlabel('Time (Hours)')
    pl.ylabel('Difference in Salinity')
    pl.ylim(0,0.6)
    pl.xlim(0,24)
    pl.savefig('figs/boxplot_salt.png')



def patch(tr):
    
    figpref.current()
    pl.figure(1)
    sH0,_ = tr.scatter(jd=732371.125+0,s=3,c='g',clf=False)
    sH1,_ = tr.scatter(jd=732371.125+1,s=3,c='r',clf=False)
    sH2,_ = tr.scatter(jd=732371.125+2,s=3,c='b',clf=False)

    pl.legend((sH0,sH1,sH2),('0 Hours','24 Hours','48 Hours') )
    pl.savefig('figs/patch_days.png')


def trajs(tr):

    pl.clf()
    jdmin = tr.jd.min()

    hmask = (tr.jd-np.floor(tr.jd))*24 == np.floor((tr.jd-np.floor(tr.jd))*24)
    dmask = (tr.jd-np.floor(tr.jd)) == np.floor((tr.jd-np.floor(tr.jd)))

    tr.scatter(jd=jdmin, s=5, c='g',clf=False)

    nvec = [1,  20, 40, 220,240,260,400,420,440]
    cvec = ['k','b','r','g','c','y','m','b','r']
    nvec = [1,  20, 40, 220,400]
    cvec = ['k','b','r','g','m']

    for n,c in zip(nvec,cvec):
        tr.scatter(ntrac=n, mask=dmask, s=30, c=c,clf=False)
    pl.title ('Positions every 24 hours')
    pl.savefig('figs/trajs_days.png')


    for n,c in zip(nvec,cvec):
        tr.scatter(ntrac=n, mask=hmask, s=10, c=c,clf=False)
    pl.title ('Positions every hour')
    pl.savefig('figs/trajs_hours.png')

    #for dt in np.arange(0.125,1.125,0.125):
    #    y,x = np.histogram(np.abs(gen_dvec(dt)),250)
    #    pl.plot(x[:-1],y)



    #pl.pcolormesh(cs.llon,cs.llat,miv(cs.salt[0,:,:]),cmap=cm.Paired)


def saltmovie(tr):
    mv = anim.Movie()
    ints = np.sort(np.unique(tr.ints))
    for i in ints[:40]:
        pl.pcolormesh(miv(tr.traj2grid('z',tr.jd==i)),cmap=cm.Paired)
        pl.clim(28,33)
        print i
        mv.image()
    mv.video(tr.projname+tr.casename+"_salt_mov.mp4")




def degrade():
    import subprocess as sp
    import os

    tr = trm.trm('casco',"decorr")
    tr.load(intstart=5859329) # 04-15

    pl.clf()


    dslist = []
    
    for t in np.arange(1,9):
        #proc = sp.Popen(['/Users/bror/git/orm/runtraj', 'degrade', str(t)] )
        #proc.wait()

        os.system('/Users/bror/git/orm/runtraj degrade '+ str(t))
        td = trm.trm('casco',"degrade")
        td.load(intstart=5859329) # 04-15
        dist = ([nanmean(tr.dist(td,tr.jd==tr.jdmin+h/24,
                                    td.jd==td.jdmin+h/24))
                                    for h in np.arange(48.)])
        pl.plot(np.array(dist)*0.2)
        dslist.append(np.array(dist))
    return dslist
        #stload.save(td,'data/degrade_%i' %t)

def degrade_plot(dslist):
    dsmat = np.array(dslist)
    pl.clf()
    lH0 = pl.plot(dsmat[0,:]*0.2)
    lH3 = pl.plot(dsmat[3,:]*0.2)
    pl.legend( ('Velocity fields every 6 hours',
                 'Velocity fields every 12 hours') ,'upper left')
    pl.xlabel('Time')
    pl.xlabel('Time (hours)')
    pl.ylabel('Mean distance (km)')
    pl.title('Particle positions compared to velocity fields every 3 hours')
    pl.savefig('figs/degrade.png')
