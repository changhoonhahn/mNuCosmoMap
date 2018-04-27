'''

Deal with the different catalogs associated with Paco's massive
neutrino simulations 


'''
import os 
import h5py 
import numpy as np 
# -- local -- 
from mnucosmomap import util as UT 
from mnucosmomap import readsnap as ReadSnap


def mNuParticles(mneut, nreal, nzbin, sim='paco'): 
    ''' Read in snapshot generated from Paco's code 

    parameters
    ----------
    mneut : float, 
        total neutrino mass 

    nreal : int,
        realization number 

    nzbin : int, 
        integer specifying the redshift of the snapshot. 
        nzbin = 0 --> z=3
        nzbin = 1 --> z=2
        nzbin = 2 --> z=1
        nzbin = 3 --> z=0.5
        nzbin = 4 --> z=0

    sim : str
        option kwarg to specify which simulation. At this moment
        this doesn't really do much since we only have paco's 
        simulatin 
    '''
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 
    _dir = ''.join([UT.dat_dir(), 'sims/paco/', # directory 
        str(mneut), 'eV/', str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'snap_', str(nzbin).zfill(3)]) # snapshot file 
    if not os.path.isfile(f): raise ValueError("file %s not found" % f)

    # read in Gadget header
    header = ReadSnap.read_gadget_header(f)

    # read in CDM particles (parttype = 1) and create catalogue
    read_keys = ['POS ', 'VEL ', 'ID ', 'MASS '] # currently only reading POS, VEL, ID, and MASS
    save_keys = ['Position', 'Velocity', 'ID', 'Mass'] 
    
    particle_data = {} 
    particle_data['meta'] = header # store meta data
    for k_s, k_r in zip(save_keys, read_keys): 
        particle_data[k_s] = ReadSnap.read_block(f, k_r, parttype=1)
        if k_s == 'Position': 
            particle_data[k_s] /= 1000. # convert ot Mpc/h 
    return particle_data 

