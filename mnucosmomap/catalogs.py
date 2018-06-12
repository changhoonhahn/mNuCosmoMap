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


def mNuDispField_subbox(nsubbox, mneut, nreal, nzbin, sim='paco', nside=8, 
        overwrite=False, verbose=False): 
    ''' Read in (or construct) displacement field of particles that are within
    the subbox in the initial conditions.
    '''
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 
    if mneut == 0.1: str_mneut = '0.10'
    else: str_mneut = str(mneut) 
    _dir = ''.join([UT.dat_dir(), 'sims/paco/', # directory 
        str_mneut, 'eV/', str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.dispfield.nside', 
        str(nside), '.', str(nsubbox), '.hdf5']) # snapshot 

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
        # read in displacement of all the particles
        dispfield = mNuDispField(mneut, nreal, nzbin, sim=sim, 
                overwrite=overwrite, verbose=verbose)
        # read in subbox indicies
        subb = mNuICs_subbox(nreal, nsubbox, sim=sim, nside=nside, verbose=verbose)
        sub_shape = subb['ID'].shape
        sub_id = subb['ID'].flatten()

        subbox = {}  
        subbox['meta'] = dispfield['meta'].copy() 
        subbox['ID'] = sub_id  
        subbox['dispfield'] = np.array([
            dispfield['dispfield'][:,0][(sub_id - 1).astype('int')].reshape(sub_shape), 
            dispfield['dispfield'][:,1][(sub_id - 1).astype('int')].reshape(sub_shape), 
            dispfield['dispfield'][:,2][(sub_id - 1).astype('int')].reshape(sub_shape)])
        if verbose: print(subbox['dispfield'].shape)

        fsub = h5py.File(f, 'w') # save ID's to hdf5 file 
        for k in dispfield['meta'].keys(): # store meta data 
            fsub.attrs.create(k, dispfield['meta'][k]) 
        fsub.create_dataset('ID', data=subbox['ID']) 
        fsub.create_dataset('dispfield', data=subbox['dispfield']) 
        fsub.close() 
    return subbox


def mNuICs_subbox(nreal, nsubbox, sim='paco', nside=8, overwrite=False, verbose=False): 
    ''' Read in (and write out if it doesn't exist) particle ID of initial condition
    particles within subbox  
    '''
    assert 512 % nside == 0 
    if verbose: print('reading in %i of %i^3 subboxes' % (nsubbox, nside)) 
    if nsubbox > nside**3: 
        raise ValueError('%i exceeds number of subboxes specifed' % nsubbox) 
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 

    _dir = ''.join([UT.dat_dir(), 'sims/paco/0.0eV/', str(nreal), '/ICs/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'ics.nside', str(nside), '.', str(nsubbox), '.hdf5']) # snapshot 

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
        fullbox = mNuICs(nreal, sim='paco', verbose=verbose) # read in full IC 

        # append extra metadata
        fullbox['meta']['n_side'] = nside
        fullbox['meta']['n_subbox'] = nsubbox

        x, y, z = fullbox['Position'].T

        L_subbox = 1000. / float(nside) # L_subbox
        L_res = 1000./512.
        L_halfres = 0.5 * L_res
        N_partside = 512/nside
        N_subbox = (N_partside)**3

        i_x = ((nsubbox % nside**2) % nside) 
        i_y = ((nsubbox % nside**2) // nside) 
        i_z = (nsubbox // nside**2) 

        xlim = ((x > L_subbox * float(i_x) + L_halfres) & 
                (x < L_subbox * float(i_x + 1) + L_halfres))
        ylim = ((y > L_subbox * float(i_y) + L_halfres) & 
                (y < L_subbox * float(i_y + 1) + L_halfres))
        zlim = ((z > L_subbox * float(i_z) + L_halfres) & 
                (z < L_subbox * float(i_z + 1) + L_halfres)) 
        in_subbox = (xlim & ylim & zlim)
        assert np.sum(in_subbox) == N_subbox
        
        ID_sub = fullbox['ID'][in_subbox]
        x_sub = x[in_subbox] - i_x * L_subbox
        y_sub = y[in_subbox] - i_y * L_subbox
        z_sub = z[in_subbox] - i_z * L_subbox
        if verbose: 
            print('%f < x_sub < %f' % (x_sub.min(), x_sub.max()))
            print('%f < y_sub < %f' % (y_sub.min(), y_sub.max()))
            print('%f < z_sub < %f' % (z_sub.min(), z_sub.max()))
        subbox_ID = np.zeros((N_partside, N_partside, N_partside))
        for j_z in range(N_partside): 
            #print('j_z = %i , %f < z < %f' % (j_z, L_res* float(j_z) + L_halfres, L_res * float(j_z + 1) + L_halfres))
            zlim_sub = ((z_sub > L_res* float(j_z) + L_halfres) & 
                    (z_sub < L_res * float(j_z + 1) + L_halfres))
            for j_y in range(N_partside): 
                #print('j_y = %i , %f < y < %f' % (j_y, L_res* float(j_y) + L_halfres, L_res * float(j_y + 1) + L_halfres))
                ylim_sub = ((y_sub > L_res * float(j_y) + L_halfres) & 
                        (y_sub < L_res * float(j_y + 1) + L_halfres))
                for j_x in range(N_partside): 
                    j_x_sorted = np.argsort(x_sub[ylim_sub & zlim_sub])
                    subbox_ID[:,j_y,j_z] = ID_sub[ylim_sub & zlim_sub][j_x_sorted]
    
        subbox = {}  
        subbox['meta'] = fullbox['meta'].copy() 
        subbox['ID'] = subbox_ID

        fsub = h5py.File(f, 'w') # save ID's to hdf5 file 
        for k in fullbox['meta'].keys(): # store meta data 
            fsub.attrs.create(k, fullbox['meta'][k]) 
        fsub.create_dataset('ID', data=subbox['ID']) 
        fsub.close() 

    return subbox 


def mNuDispField(mneut, nreal, nzbin, sim='paco', overwrite=False, verbose=False): 
    ''' Read in displacement field of the snapshot generated from Paco's code 

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
    f = ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.dispfield.hdf5']) # snapshot 
    
    if os.path.isfile(f) and not overwrite: # read in the file 
        dispfield = {} 
        
        fdisp = h5py.File(f, 'r') # read in hdf5 file 
        dispfield['meta'] = {} 
        for k in fdisp.attrs.keys(): 
            dispfield['meta'][k] = fdisp.attrs[k]
        # read in datasets
        for k in fdisp.keys(): 
            dispfield[k] = fdisp[k].value 
        fdisp.close() 

    else: # write file 
        # read in particles of snapshot 
        par = mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=verbose)
        # read in initial conditions 
        ics = mNuICs(nreal, sim=sim, verbose=verbose) 

        isort_ics = np.argsort(ics['ID']) 
        isort_par = np.argsort(par['ID']) 
        assert np.array_equal(ics['ID'][isort_ics], par['ID'][isort_par])

        dispfield = {} 
        dispfield['meta'] = par['meta'].copy() 
        dispfield['ID'] = ics['ID'][isort_ics]
        dispfield['dispfield'] = ics['Position'][isort_ics] - par['Position'][isort_par]
        
        fdisp = h5py.File(f, 'w') # save ID's to hdf5 file 
        for k in dispfield.keys(): 
            if k == 'meta': continue 
            fdisp.create_dataset(k, data=dispfield[k]) 

        for k in par['meta'].keys(): # store meta data 
            fdisp.attrs.create(k, dispfield['meta'][k]) 
        fdisp.close() 

    return dispfield 


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


def mNuICs(nreal, sim='paco', verbose=False): 
    ''' Read in initial condition of  Paco's simulations. According to 
    Paco, each of the realizations (regardless of neutrino mass) should 
    have the same initial conditions. 

    parameters
    ----------
    nreal : int,
        realization number 

    sim : str
        option kwarg to specify which simulation. At this moment
        this doesn't really do much since we only have paco's 
        simulatin 
    '''
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 
    _dir = ''.join([UT.dat_dir(), 'sims/paco/0.0eV/', str(nreal), '/ICs/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'ics']) # snapshot 
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
