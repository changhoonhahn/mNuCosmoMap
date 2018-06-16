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


def _mNuDispField_subbox_mneut(): 
    '''
    '''
    dfields = [] 
    ics = [] 
    mneuts = [0.0, 0.06, 0.10, 0.15, 0.6]
    for mneut in mneuts: 
        dfield = mNuCat.mNuDispField_subbox(1, mneut, 1, 2, sim='paco', boundary_correct=True, overwrite=True, verbose=True) 
        dfields.append(dfield) 
        ics_subbox = mNuCat.mNuICs_subbox(1, 1, sim='paco', verbose=False)
        ics.append(ics_subbox)

    for i, dfield, ic in zip(range(len(mneuts)), dfields, ics): 
        fig = plt.figure(figsize=(15,5))
        for ii in range(3):
            sub = fig.add_subplot(1,3,ii+1)
            print('%f < d_x < %f' % 
                    (dfield['DispField'][0,:,:,10].min(), dfield['DispField'][0,:,:,10].max()))
            print('%f < d_y < %f' % 
                    (dfield['DispField'][1,:,:,10].min(), dfield['DispField'][1,:,:,10].max()))

            sub.quiver(ic['Position'][0,:,:,5*ii], ic['Position'][1,:,:,5*ii],
                    dfield['DispField'][0,:,:,5*ii], dfield['DispField'][1,:,:,5*ii], 
                    scale=1, scale_units='xy')

            sub.set_xlabel('X', fontsize=20) 
            sub.set_xlim([120., 255.]) 
            if ii == 0: sub.set_ylabel('Y', fontsize=20) 
            sub.set_ylim([-0.5, 130.]) 
            if ii == 1: sub.set_title(r'$\sum m_\nu = $ %f eV' % round(mneuts[i],2), fontsize=20) 
        fig.savefig(''.join([UT.fig_dir(), '_mNuDispField_subbox_mneut', str(round(mneuts[i],2)), '.png']), 
                bbox_inches='tight') 
    return None


def _mNuDispField_subbox(): 
    '''
    '''
    subbox = mNuCat.mNuParticles_subbox(1, 0.0, 1, 2, sim='paco', nside=8, verbose=False) 
    ics_subbox = mNuCat.mNuICs_subbox(1, 1, sim='paco', nside=8, verbose=False)
    dfield = mNuCat.mNuDispField_subbox(1, 0.0, 1, 2, sim='paco', boundary_correct=False, verbose=True) 

    fig = plt.figure(figsize=(15,5))
    #X, Y = np.meshgrid(np.arange(dfield['dispfield'].shape[1]), np.arange(dfield['dispfield'].shape[2]))
    for i in range(3): 
        sub = fig.add_subplot(1,3,i+1)
        #sub.quiver(X, Y, dfield['dispfield'][0,:,:,10*i], dfield['dispfield'][1,:,:,10*i]) 
        sub.scatter(subbox['Position'][0,:,:,5*i], subbox['Position'][1,:,:,5*i], 
                c='C0', s=0.4) 
        sub.scatter(
                ics_subbox['Position'][0,:,:,5*i] + dfield['DispField'][0,:,:,5*i], 
                ics_subbox['Position'][1,:,:,5*i] + dfield['DispField'][1,:,:,5*i],  
                c='C1', s=0.2) 
        sub.set_xlabel('X', fontsize=25) 
        sub.set_xlim([120., 255.]) 
        if i == 0: sub.set_ylabel('Y', fontsize=20) 
        sub.set_ylim([-0.5, 130.]) 
        sub.set_title(r'$i_z =$ %i' % (10*i), fontsize=20) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuDispField_subbox.png']), bbox_inches='tight') 
    return None


def _mNuParticles_subbox_mneut(): 
    '''check that the function mnucosmomap.catalogs.mNuParticles_subbox is 
    properly reading in the snapshot particles 
    '''
    subboxes = [] 
    mneuts = [0.0, 0.06, 0.10, 0.15, 0.6]
    for mneut in mneuts: 
        subbox = mNuCat.mNuParticles_subbox(1, mneut, 1, 2, sim='paco', nside=8, 
                overwrite=False, verbose=False) 
        subboxes.append(subbox)
    
    ics_subbox = mNuCat.mNuICs_subbox(1, 1, sim='paco', nside=8, overwrite=False, verbose=False)

    for i in range(3): 
        print('%i : %f - %f' % (i, 
            ics_subbox['Position'][i,:,:,:].flatten().min(), 
            ics_subbox['Position'][i,:,:,:].flatten().max()))

    for i, subbox in enumerate(subboxes): 
        print('subbox %i' % i)
        for ii in range(3): 
            print('%i : %f - %f' % (ii, 
                ics_subbox['Position'][ii,:,:,:].flatten().min(), 
                ics_subbox['Position'][ii,:,:,:].flatten().max()))
            print('%i : %f - %f' % (ii, 
                subbox['Position'][ii,:,:,:].flatten().min(), 
                subbox['Position'][ii,:,:,:].flatten().max()))
        fig = plt.figure(figsize=(15,5))
        for ii in range(3): 
            sub = fig.add_subplot(1,3,ii+1)
            sub.scatter(ics_subbox['Position'][0,:,:,5*ii], ics_subbox['Position'][1,:,:,5*ii], 
                    c='k', s=0.2) 
            sub.scatter(subbox['Position'][0,:,:,5*ii], subbox['Position'][1,:,:,5*ii], 
                    c='C1', s=0.2) 
            sub.set_xlabel('X', fontsize=20) 
            sub.set_xlim([120., 255.]) 
            if ii == 0: sub.set_ylabel('Y', fontsize=20) 
            sub.set_ylim([-0.5, 130.]) 
            if ii == 1: sub.set_title(r'$\sum m_\nu = $ %f eV' % round(mneuts[i],2), fontsize=20) 
        fig.savefig(''.join([UT.fig_dir(), '_mNuParticles_subbox_mneut', str(round(mneuts[i],2)), '.png']), 
                bbox_inches='tight') 
    return None 


def mNuParticles_subbox(mneut=0.0, nreal=1, nzbin=2, sim='paco', nside=8): 
    '''check that the function mnucosmomap.catalogs.mNuParticles_subbox is 
    properly reading in the snapshot particles 
    '''
    # read/write subbox 
    subboxes = mNuCat.mNuParticles_subbox(range(5), mneut, nreal, nzbin, sim=sim, nside=nside, verbose=False)
    
    # read in the particles
    box = mNuCat.mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=False)

    # now run some common sense checks 
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.scatter(box['Position'][:,0][::10000], box['Position'][:,1][::10000], c='k', s=0.1) 
    for i, subbox in enumerate(subboxes): 
        sub.scatter(subbox['Position'][0,:,:,0], subbox['Position'][1,:,:,0], c='C'+str(i), s=0.1) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([0., 1000.]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([0., 1000.]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuParticles_subbox.png']), bbox_inches='tight') 
    return None 


def mNuICs_subbox(nsubbox=1, nreal=1, sim='paco'): 
    ''' examine particles from catalog.mNuICs within
    a 10^3 Mpc subbox.
    '''
    # read in initial condition particles within subbox
    subbox = mNuCat.mNuICs_subbox(nsubbox, nreal, sim=sim, nside=8, verbose=False)
    x_sub = subbox['Position'][0,:,:,:].flatten()
    y_sub = subbox['Position'][1,:,:,:].flatten()
    z_sub = subbox['Position'][2,:,:,:].flatten()

    i_x, i_y, i_z = subbox['meta']['subbox_ijk'] 

    fig = plt.figure(figsize=(12,4))
    sub = fig.add_subplot(131) # x vs y 
    sub.scatter(x_sub, y_sub, c='k', s=1) 
    sub.set_xlabel('X', fontsize=20) 
    sub.set_xlim([(1000/8)*i_x-5, (1000/8)*(i_x+1)+5]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([(1000/8)*i_y-5, (1000/8)*(i_y+1)+5]) 

    sub = fig.add_subplot(132) # z vs y 
    sub.scatter(z_sub, y_sub, c='k', s=1) 
    sub.set_xlabel('Z', fontsize=20) 
    sub.set_xlim([(1000/8)*i_z-5, (1000/8)*(i_z+1)+5]) 
    sub.set_ylabel('Y', fontsize=20) 
    sub.set_ylim([(1000/8)*i_y-5, (1000/8)*(i_y+1)+5]) 
    
    sub = fig.add_subplot(133) # z vs y 
    sub.scatter(z_sub, x_sub, c='k', s=1) 
    sub.set_xlabel('Z', fontsize=20) 
    sub.set_xlim([(1000/8)*i_z-5, (1000/8)*(i_z+1)+5]) 
    sub.set_ylabel('X', fontsize=20) 
    sub.set_ylim([(1000/8)*i_x-5, (1000/8)*(i_x+1)+5]) 
    fig.savefig(''.join([UT.fig_dir(), '_mNuICs_subbox', str(nsubbox), '.png']), 
            bbox_inches='tight') 
    return None 


def mNuParticles(mneut=0.0, nreal=1, nzbin=2, sim='paco'): 
    '''check that the function mnucosmomap.catalogs.mNuParticles is 
    properly reading in the snapshot particles 
    '''
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


def mNuICs(nreal=1, sim='paco'): 
    ''' check that catalog.mNuICs reads intial conditions properly  
    '''
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


if __name__=="__main__": 
    #mNuICs()
    #mNuParticles()
    #mNuICs_subbox(nsubbox=0)
    mNuParticles_subbox()
    #_mNuParticles_subbox_mneut() 
    #_mNuDispField_subbox_mneut()
    #_mNuDispField_subbox()
