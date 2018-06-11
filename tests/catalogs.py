'''

Tests for mnucosmomap.catalogs


'''
from mnucosmomap import util as UT 
from mnucosmomap import catalogs as mNuCat 

# -- plotting --
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['axes.xmargin'] = 1
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['legend.frameon'] = False


def _mNuParticles_subbox(): 
    '''check that the function mnucosmomap.catalogs.mNuParticles_subbox is 
    properly reading in the snapshot particles 
    '''
    mneut = 0.0 
    nreal = 1 
    nzbin = 2 # z=1 
    sim = 'paco' # only paco is supported anyway. 
    nside = 8 
    
    # read/write subbox 
    subboxes = [] 
    for i in range(5):
        sb = mNuCat.mNuParticles_subbox(mneut, nreal, nzbin, i, sim=sim, nside=nside, verbose=False)
        subboxes.append(sb) 
    
    # read in the particles
    box = mNuCat.mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=False)

    # now run some common sense checks 
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.scatter(box['Position'][:,0][::10000], box['Position'][:,1][::10000], c='k', s=0.1) 
    for i, subbox in enumerate(subboxes): 
        sub.scatter(subbox['Position'][:,0], subbox['Position'][:,1], c='C'+str(i), s=0.1) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([0., 1000.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 1000.]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuParticles_subbox.png']), bbox_inches='tight') 
    return None 


def _mNuParticles(): 
    '''check that the function mnucosmomap.catalogs.mNuParticles is 
    properly reading in the snapshot particles 
    '''
    mneut = 0.0 
    nreal = 1 
    nzbin = 2 # z=1 
    sim = 'paco' # only paco is supported anyway. 
    
    # read in the particles
    cat = mNuCat.mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=False)
    for x, i in zip(['x', 'y', 'z'], range(3)): 
        print('%s axis ranges from %f to %f' % (x, cat['Position'][:,i].min(), cat['Position'][:,i].max())) 

    # now run some common sense checks 
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.scatter(cat['Position'][:,0][::10000], cat['Position'][:,1][::10000], c='k', s=0.1) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([0., 1000.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 1000.]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuParticles.png']), bbox_inches='tight') 
    return None 


def _mNuICs_subbox(): 
    ''' examine particles from catalog.mNuICs within
    a 10^3 Mpc subbox.
    '''
    nreal = 1 
    sim = 'paco' # only paco is supported anyway. 
    
    # read in the particles
    cat = mNuCat.mNuICs(nreal, sim=sim, verbose=False)
    insubbox = ((cat['Position'][:,0] < 10.) & (cat['Position'][:,1] < 10.) & (cat['Position'][:,2] < 10.)) 
    print(cat['ID'][insubbox])

    fig = plt.figure()
    sub = fig.add_subplot(131) # x vs y 
    sub.scatter(cat['Position'][:,0][insubbox], cat['Position'][:,1][insubbox], c='k', s=3) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([0., 10.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 10.]) 

    sub = fig.add_subplot(132) # z vs y 
    sub.scatter(cat['Position'][:,2][insubbox], cat['Position'][:,1][insubbox], c='k', s=3) 
    sub.set_xlabel('Z', fontsize=20) 
    sub.set_xlim([0., 10.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 10.]) 
    
    sub = fig.add_subplot(133) # z vs y 
    sub.scatter(cat['Position'][:,2][insubbox], cat['Position'][:,0][insubbox], c='k', s=3) 
    sub.set_xlabel('Z', fontsize=20) 
    sub.set_xlim([0., 10.]) 
    sub.set_ylabel('X', fontsize=20) 
    sub.set_ylim([0., 10.]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuICs_subbox.png']), bbox_inches='tight') 
    return None 


def _mNuICs(): 
    ''' check that catalog.mNuICs reads intial conditions properly  
    '''
    nreal = 1 
    sim = 'paco' # only paco is supported anyway. 
    
    # read in the particles
    cat = mNuCat.mNuICs(nreal, sim=sim, verbose=False)
    for x, i in zip(['x', 'y', 'z'], range(3)): 
        print('%s axis ranges from %f to %f' % (x, cat['Position'][:,i].min(), cat['Position'][:,i].max())) 
    
    # now run some common sense checks 
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.scatter(cat['Position'][:,0][::10000], cat['Position'][:,1][::10000], c='k', s=0.1) 
    sub.scatter(cat['Position'][:100,0], cat['Position'][:100,1], c='C1', s=5) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([0., 1000.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 1000.]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuICs.png']), bbox_inches='tight') 
    return None 


if __name__=="__main__": 
    _mNuICs()
    _mNuICs_subbox()
