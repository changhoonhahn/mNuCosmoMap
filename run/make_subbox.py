import sys 
import numpy as np 
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


def mNuDispField_subbox(mneut, nreal, nzbin, sim='paco', nside=8, nsubbox=None): 
    '''
    '''
    if nsubbox is None: nsubboxes = range(nside**3)
    else: nsubboxes = [nsubbox]
    for i in nsubboxes: # write subbox 
        print('writing subbox %i of %i' % (i, nside**3))
        dfield = mNuCat.mNuDispField_subbox(i, mneut, nreal, nzbin, 
                sim=sim, boundary_correct=True, verbose=False) 
    return None


def mNuParticles_subbox(mneut, nreal, nzbin=2, sim='paco', nside=8, nsubbox=None): 
    ''' Generate suboxes for given mneut, nreal, nzbin, sim, and nside using 
    mnucosmomap.catalogs.mNuParticles_subbox  
    '''
    if nsubbox is None: nsubbox = nside**3

    for i in range(nsubbox): # write subbox 
        print('writing subbox %i of %i' % (i, nside**3))
        sb = mNuCat.mNuParticles_subbox(mneut, nreal, nzbin, i, sim=sim, nside=nside, verbose=False)
    return None 


def _check_subbox(mneut, nreal, i_subbox, nzbin=2, sim='paco', nside=8): 
    ''' Check a subbox
    '''
    # read in the particles
    #box = mNuCat.mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=False)
    # read in subbox 
    sb = mNuCat.mNuParticles_subbox(mneut, nreal, nzbin, i_subbox, sim=sim, nside=nside, 
            verbose=False)

    # now run some common sense checks 
    fig = plt.figure(figsize=(12,4))
    coord = ['X', 'Y', 'Z'] 
    pairs = [(0, 1), (2, 1), (0, 2)]  
    for i in range(len(pairs)): 
        sub = fig.add_subplot(1,3,i+1)
        #sub.scatter(box['Position'][:,pairs[i][0]][::10000], box['Position'][:,pairs[i][1]][::10000], c='k', s=0.1) 
        sub.scatter(sb['Position'][:,pairs[i][0]][::100], sb['Position'][:,pairs[i][1]][::100], c='C1', s=0.1) 
        sub.set_xlabel(coord[pairs[i][0]], fontsize=20) 
        sub.set_xlim([0., 1000.]) 
        sub.set_ylabel(coord[pairs[i][1]], fontsize=20) 
        sub.set_ylim([0., 1000.]) 

    fig.subplots_adjust(wspace=0.15)
    fig.savefig(''.join([UT.fig_dir(), 
        '_check_subbox.', str(mneut), 'eV.', str(nreal), '.subbox', str(i_subbox), '.png']), 
        bbox_inches='tight') 
    plt.close()
    return None 


if __name__=="__main__": 
    mneut = float(sys.argv[1].strip('eV')) 
    nreal = int(sys.argv[2]) 
    nzbin = int(sys.argv[3]) 
    nside = int(sys.argv[4]) 
    nsubbox = int(sys.argv[5]) 
    mNuDispField_subbox(mneut, nreal, nzbin=nzbin, nside=nside, nsubbox=nsubbox)
    
    # randomly check a few subboxes to make sure. 
    i_rand = np.random.choice(range(nsubbox), size=3, replace=False) 
    for i in i_rand: 
        _check_subbox(mneut, nreal, i, nzbin=nzbin, nside=nside)
