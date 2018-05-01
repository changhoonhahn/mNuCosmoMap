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


def mNuParticles_subbox(mneut, nreal, nzbin, nsubbox, sim='paco', nside=8, overwrite=False, verbose=False): 
    ''' Read in (and write out if it doesn't exist) subbox of snapshots generated from 
    '''
    if verbose: print('reading in %i of %i^3 subboxes' % (nsubbox, nside)) 
    if nsubbox > nside**3: 
        raise ValueError('%i exceeds number of subboxes specifed' % nsubbox) 
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 
    
    if mneut == 0.1: str_mneut = '0.10'
    else: str_mneut = str(mneut) 
    _dir = ''.join([UT.dat_dir(), 'sims/paco/', # directory 
        str_mneut, 'eV/', str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    # snapshot subbox file 
    f = ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.nside', str(nside), '.', str(nsubbox), '.hdf5']) 

    if os.path.isfile(f) and not overwrite: # read in the file 
        subbox = {} 
        
        fsub = h5py.File(f, 'r') # read in hdf5 file 
        subbox['meta'] = {} 
        for k in fsub.attrs.keys(): 
            subbox['meta'][k] = fsub.attrs[k]
        # read in datasets
        for k in fsub.keys(): 
            subbox[k] = fsub[k].value 
        fsub.close() 

    else: # write file 
        fullbox = mNuParticles(mneut, nreal, nzbin, sim='paco', verbose=verbose) # read in full box

        # append extra metadata
        fullbox['meta']['n_side'] = nside
        fullbox['meta']['n_subbox'] = nsubbox

        x, y, z = fullbox['Position'].T

        L_subbox = 1000. / float(nside) # L_subbox

        i_x = ((nsubbox % nside**2) % nside) 
        i_y = ((nsubbox % nside**2) // nside) 
        i_z = (nsubbox // nside**2) 

        xlim = (x >= L_subbox * float(i_x)) & (x < L_subbox * float(i_x + 1))
        ylim = (y >= L_subbox * float(i_y)) & (y < L_subbox * float(i_y + 1))
        zlim = (z >= L_subbox * float(i_z)) & (z < L_subbox * float(i_z + 1))
    
        subbox = {}  
        subbox['meta'] = fullbox['meta'].copy() 

        fsub = h5py.File(f, 'w') # save to hdf5 file 
        for k in fullbox.keys(): 
            if k == 'meta': continue 
            subbox[k] = fullbox[k][xlim & ylim & zlim]
            fsub.create_dataset(k, data=subbox[k]) 

        for k in fullbox['meta'].keys(): # store meta data 
            fsub.attrs.create(k, fullbox['meta'][k]) 
        fsub.close() 

    return subbox 


def mNuParticles(mneut, nreal, nzbin, sim='paco', verbose=False): 
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
    if mneut == 0.1: str_mneut = '0.10'
    else: str_mneut = str(mneut) 
    _dir = ''.join([UT.dat_dir(), 'sims/paco/', # directory 
        str_mneut, 'eV/', str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'snap_', str(nzbin).zfill(3)]) # snapshot 
    #if not os.path.isfile(f): raise ValueError("file %s not found" % f)

    # read in Gadget header
    header = ReadSnap.read_gadget_header(f)

    # read in CDM particles (parttype = 1) and create catalogue
    read_keys = ['POS ', 'VEL ', 'ID  ', 'MASS'] # currently only reading POS, VEL, ID, and MASS
    save_keys = ['Position', 'Velocity', 'ID', 'Mass'] 
    
    particle_data = {} 
    particle_data['meta'] = header # store meta data
    for k_s, k_r in zip(save_keys, read_keys): 
        particle_data[k_s] = ReadSnap.read_block(f, k_r, parttype=1, verbose=verbose)
        if k_s == 'Position': 
            particle_data[k_s] /= 1000. # convert ot Mpc/h 
    return particle_data 
