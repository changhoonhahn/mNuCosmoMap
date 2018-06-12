'''

Tests for mnucosmomap.catalogs


'''
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

    # read in particles ID within subbox
    subbox = mNuCat.mNuICs_subbox(nreal, 1, sim=sim, nside=8, verbose=False)
    sub_id = subbox['ID'].flatten()
    # read in initial conditions  
    ics = mNuCat.mNuICs(nreal, sim=sim, verbose=False)
    
    isort_ics = np.argsort(ics['ID']) 
    assert np.array_equal(ics['ID'][isort_ics][(sub_id-1).astype('int')], sub_id)
    x_sub = ics['Position'][:,0][isort_ics][(sub_id-1).astype('int')]
    y_sub = ics['Position'][:,1][isort_ics][(sub_id-1).astype('int')]
    z_sub = ics['Position'][:,2][isort_ics][(sub_id-1).astype('int')]

    fig = plt.figure(figsize=(12,4))
    sub = fig.add_subplot(131) # x vs y 
    sub.scatter(x_sub, y_sub, c='k', s=3) 
    sub.set_xlabel('X', fontsize=20) 
    #sub.set_xlim([-0.1, 10.]) 
    sub.set_ylabel('Y', fontsize=20) 
    #sub.set_ylim([-0.1, 10.]) 

    sub = fig.add_subplot(132) # z vs y 
    sub.scatter(z_sub, y_sub, c='k', s=3) 
    sub.set_xlabel('Z', fontsize=20) 
    #sub.set_xlim([-0.1, 10.]) 
    sub.set_ylabel('Y', fontsize=20) 
    #sub.set_ylim([-0.1, 10.]) 
    
    sub = fig.add_subplot(133) # z vs y 
    sub.scatter(z_sub, x_sub, c='k', s=3) 
    sub.set_xlabel('Z', fontsize=20) 
    #sub.set_xlim([-0.1, 10.]) 
    sub.set_ylabel('X', fontsize=20) 
    #sub.set_ylim([-0.1, 10.]) 
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
    sub.scatter(cat['Position'][:1000,0], cat['Position'][:1000,1], c='C1', s=1) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([0., 1000.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 1000.]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuICs.png']), bbox_inches='tight') 
    return None 


def _mNuDispField_subbox(): 
    '''
    '''
    dfield = mNuCat.mNuDispField_subbox(1, 0.0, 1, 2, sim='paco', boundary_correct=True, verbose=True) 
    fig = plt.figure(figsize=(15,5))
    X, Y = np.meshgrid(np.arange(dfield['dispfield'].shape[1]), np.arange(dfield['dispfield'].shape[2]))
    for i in range(3): 
        sub = fig.add_subplot(1,3,i+1)
        print('%f < d_x < %f' % 
                (dfield['dispfield'][0,:,:,10*i].min(), dfield['dispfield'][0,:,:,10*i].max()))
        print('%f < d_y < %f' % 
                (dfield['dispfield'][1,:,:,10*i].min(), dfield['dispfield'][1,:,:,10*i].max()))
        sub.quiver(X, Y, dfield['dispfield'][0,:,:,10*i], dfield['dispfield'][1,:,:,10*i]) 
        sub.set_xlabel('X', fontsize=25) 
        sub.set_xlim([0., dfield['dispfield'].shape[1]])
        if i == 0: sub.set_ylabel('Y', fontsize=25) 
        sub.set_ylim([0., dfield['dispfield'].shape[2]])
        sub.set_title(r'$i_z =$ %i' % (10*i), fontsize=20) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuDispField_subbox.png']), bbox_inches='tight') 
    return None


def _mNuDispField(): 
    ''' Check that the displacement fields are properly calculated.  
    '''
    dfield = mNuCat.mNuDispField(0.0, 1, 2, boundary_correct=False, sim='paco', verbose=False) 
    print('%i objects with cross boundary displacements' % np.sum(np.abs(dfield['dispfield']) > 900.))
    print('they account for %f of all the objects' % 
            (float(np.sum(np.abs(dfield['dispfield']) > 900.))/float(len(dfield['ID']))))

    dfield_corr = mNuCat.mNuDispField(0.0, 1, 2, boundary_correct=True, sim='paco', verbose=False) 

    #crossbound = (np.abs(dfield['dispfield']) > 250.)
    #crossbound_corr = (np.abs(dfield_corr['dispfield']) > 250.)
    fig = plt.figure()
    sub = fig.add_subplot(111)
    #_ = sub.hist(dfield['dispfield'], range=(-1000., 1000.), bins=20, normed=True)
    for i in range(3): 
        _ = sub.hist(dfield['dispfield'][:,i], range=(-1000., 1000.), bins=200, color='k', normed=True)
    for i in range(3): 
        _ = sub.hist(dfield_corr['dispfield'][:,i], range=(-1000., 1000.), bins=200, color='C1', normed=True)
    sub.set_xlim([-1000., 1000.]) 
    sub.set_ylim([0., 0.001]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuDispField.png']), bbox_inches='tight') 
    plt.close() 
    return None 
 

if __name__=="__main__": 
    _mNuDispField_subbox()
    #_mNuDispField()
    #_mNuICs_subbox()
    #_mNuICs()
