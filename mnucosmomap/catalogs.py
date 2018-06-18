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


def mNuDispField_subbox(nsubbox, mneut, nreal, nzbin, sim='paco', nside=8, boundary_correct=True, 
        overwrite=False, verbose=False): 
    ''' Read in (or construct) displacement field of particles that are within
    the subbox in the initial conditions.
    '''
    try: 
        Nsub = len(nsubbox) 
    except TypeError: 
        Nsub = 1 
        nsubbox = [nsubbox]
    for isubbox in nsubbox: 
        if verbose: print('reading in %i of %i^3 subboxes' % (isubbox, nside)) 
        if isubbox > nside**3: raise ValueError('%i exceeds number of subboxes specifed' % isubbox) 
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 

    _dir = ''.join([UT.mNuDir(mneut, sim=sim), str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    F = lambda isub: ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.dispfield.nside', 
        str(nside), '.', str(isub), '.hdf5']) # snapshot 

    subboxes = [] 
    for isubbox in nsubbox: 
        subbox = {} 
        if os.path.isfile(F(isubbox)) and not overwrite: # read in the file 
            fsub = h5py.File(F(isubbox), 'r') # read in hdf5 file 
            subbox['meta'] = {} 
            for k in fsub.attrs.keys(): 
                subbox['meta'][k] = fsub.attrs[k]
            # read in datasets
            for k in fsub.keys(): 
                subbox[k] = fsub[k].value 
            fsub.close() 
        else: # write file 
            print('Constructing %s ......' % F(isubbox)) 
            # read in displacement of all the particles 
            subpar = mNuParticles_subbox(isubbox, mneut, nreal, nzbin, sim=sim, nside=nside, verbose=verbose) 
            # read in subbox indicies
            subics = mNuICs_subbox(isubbox, nreal, sim=sim, nside=nside, verbose=verbose)

            subbox['meta'] = subpar['meta'].copy() 
            subbox['ID'] = subics['ID'].copy()  
            subbox['DispField'] = subpar['Position'] - subics['Position']

            fsub = h5py.File(F(isubbox), 'w') # save ID's to hdf5 file 
            for k in subbox['meta'].keys(): # store meta data 
                fsub.attrs.create(k, subbox['meta'][k]) 
            fsub.create_dataset('ID', data=subbox['ID']) 
            fsub.create_dataset('DispField', data=subbox['DispField']) 
            fsub.close() 
        subboxes.append(subbox)

    if boundary_correct: 
        # correct for the cases where particles cross the periodic boundaries 
        # the correction at the moment is brute forced. 
        for subbox in subboxes: 
            for i in range(3): 
                cross_neg = (subbox['DispField'][i,:,:,:].flatten() < -900.) 
                subbox['DispField'][i,cross_neg.reshape(subbox['DispField'].shape[1:])] += 1000.
                cross_pos = (subbox['DispField'][i,:,:,:].flatten() > 900.) 
                subbox['DispField'][i,cross_pos.reshape(subbox['DispField'].shape[1:])] = \
                        1000. - subbox['DispField'][i,cross_pos.reshape(subbox['DispField'].shape[1:])]
    if Nsub == 1: 
        return subboxes[0]
    else: 
        return subboxes


def mNuParticles_subbox(nsubbox, mneut, nreal, nzbin, sim='paco', nside=8, 
        overwrite=False, verbose=False): 
    ''' Read in (and write out if it doesn't exist) subbox of snapshots generated from 
    Paco's neutrino simulations. The positions and velocity of the CDM particles
    have dimensions 3xNxNxN
    '''
    try: 
        Nsub = len(nsubbox) 
    except TypeError: 
        Nsub = 1 
        nsubbox = [nsubbox]
    for isubbox in nsubbox: 
        if verbose: print('reading in %i of %i^3 subboxes' % (isubbox, nside)) 
        if isubbox > nside**3: raise ValueError('%i exceeds number of subboxes specifed' % isubbox) 
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 
    
    _dir = ''.join([UT.mNuDir(mneut, sim=sim), str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    # snapshot subbox file 
    F = lambda isub: ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.nside', str(nside), '.', str(isub), '.hdf5']) 

    subboxes = [] 
    for isubbox in nsubbox: 
        subbox = {} 
        if os.path.isfile(F(isubbox)) and not overwrite: # read in the file 
            fsub = h5py.File(F(isubbox), 'r') # read in hdf5 file 
            subbox['meta'] = {} 
            for k in fsub.attrs.keys(): 
                subbox['meta'][k] = fsub.attrs[k]
            for k in fsub.keys(): # read in datasets
                subbox[k] = fsub[k].value 
            fsub.close() 
        else: # write file 
            print('Constructing %s ......' % f(isubbox)) 
            if isubbox == nsubbox[0]: 
                # read in full box of Particles
                fullbox = mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=verbose)
                isort_box = np.argsort(fullbox['ID'])

            # read in subbox indicies
            subb = mNuICs_subbox(isubbox, nreal, sim=sim, nside=nside, verbose=verbose)
            N_partside, N_partside, N_partside = subb['ID'].shape
            sub_id = subb['ID'].flatten()

            subbox['meta'] = fullbox['meta'].copy() 
            subbox['meta']['n_side'] = nside # append extra metadata
            subbox['meta']['n_subbox'] = nsubbox
            subbox['meta']['subbox_ijk'] = subb['meta']['subbox_ijk']
            subbox['ID'] = subb['ID'].copy()  
            for k in ['Position', 'Velocity']: 
                subbox[k] = np.zeros((3, N_partside, N_partside, N_partside)) 
                for i in range(3): 
                    subbox[k][i,:,:,:] = \
                            fullbox[k][:,i][isort_box][(sub_id - 1).astype('int')].reshape((N_partside, N_partside, N_partside))
        
            fsub = h5py.File(F(isubbox), 'w') # save to hdf5 file 
            for k in ['ID', 'Position', 'Velocity']:  
                fsub.create_dataset(k, data=subbox[k]) 

            for k in subbox['meta'].keys(): # store meta data 
                fsub.attrs.create(k, subbox['meta'][k]) 
            fsub.close() 
        subboxes.append(subbox) 

    if Nsub == 1: 
        return subboxes[0]
    else: 
        return subboxes


def mNuICs_subbox(nsubbox, nreal, sim='paco', nside=8, overwrite=False, verbose=False): 
    ''' Read in (and write out if it doesn't exist) particle ID of initial condition
    particles within subbox  
    '''
    assert 512 % nside == 0 
    try: 
        Nsub = len(nsubbox) 
    except TypeError: 
        Nsub = 1 
        nsubbox = [nsubbox]
    for isubbox in nsubbox: 
        if isubbox > nside**3: raise ValueError('%i exceeds number of subboxes specifed' % isubbox) 
    if sim not in ['paco']: raise NotImplementedError('%s simulation not supported yet') 

    _dir = ''.join([UT.mNuDir(0.0, sim=sim), str(nreal), '/ICs/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    F = lambda isub: ''.join([_dir, 'ics.nside', str(nside), '.', str(isub), '.hdf5']) # snapshot 
    
    iread = 0 
    subboxes = [] 
    for isubbox in nsubbox: 
        subbox = {} 
        if os.path.isfile(F(isubbox)) and not overwrite: # read in the file 
            if verbose: print('reading in %i of %i^3 subboxes' % (isubbox, nside)) 
            fsub = h5py.File(F(isubbox), 'r') # read in hdf5 file 
            subbox['meta'] = {} 
            for k in fsub.attrs.keys(): 
                subbox['meta'][k] = fsub.attrs[k]
            # read in datasets
            for k in fsub.keys(): 
                subbox[k] = fsub[k].value 
            fsub.close() 
        else: # write file 
            if verbose: print('Constructing %s ......' % F(isubbox)) 
            if iread == 0: # read in full IC 
                fullbox = mNuICs(nreal, sim=sim, verbose=verbose)
                # append extra metadata
                fullbox['meta']['n_side'] = nside
                fullbox['meta']['n_subbox'] = nsubbox

                x, y, z = fullbox['Position'].T
                vx, vy, vz = fullbox['Velocity'].T

                L_subbox = 1000./float(nside) # L_subbox
                L_res = 1000./512.
                L_halfres = 0.5 * L_res
                N_partside = 512/nside
                N_subbox = (N_partside)**3
                iread += 1

            i_x = ((isubbox % nside**2) % nside) 
            i_y = ((isubbox % nside**2) // nside) 
            i_z = (isubbox // nside**2) 
            if verbose: print('%i, %i, %i' % (i_x, i_y, i_z))

            xmin = L_subbox * float(i_x) + L_halfres
            xmax = (L_subbox * float(i_x+1) + L_halfres) % 1000.
            ymin = L_subbox * float(i_y) + L_halfres
            ymax = (L_subbox * float(i_y+1) + L_halfres) % 1000.
            zmin = L_subbox * float(i_z) + L_halfres
            zmax = (L_subbox * float(i_z+1) + L_halfres) % 1000.
            if xmin <= xmax: xlim = ((x >= xmin) & (x < xmax))
            else: xlim = ((x >= xmin) | (x < xmax))
            if ymin <= ymax: ylim = ((y >= ymin) & (y < ymax))
            else: ylim = ((y >= ymin) | (y < ymax))
            if zmin <= zmax: zlim = ((z >= zmin) & (z < zmax))
            else: zlim = ((z >= zmin) | (z < zmax))
            in_subbox = (xlim & ylim & zlim)
            assert np.sum(in_subbox) == N_subbox
            
            ID_sub = fullbox['ID'][in_subbox]
            x_subbox = x[in_subbox]
            y_subbox = y[in_subbox]
            z_subbox = z[in_subbox]
            x_sub = (x_subbox - i_x * L_subbox) % 1000.
            y_sub = (y_subbox - i_y * L_subbox) % 1000.
            z_sub = (z_subbox - i_z * L_subbox) % 1000.

            vx_subbox = vx[in_subbox]
            vy_subbox = vy[in_subbox]
            vz_subbox = vz[in_subbox]
            if verbose: 
                print('%f < x_sub < %f' % (x_sub.min(), x_sub.max()))
                print('%f < y_sub < %f' % (y_sub.min(), y_sub.max()))
                print('%f < z_sub < %f' % (z_sub.min(), z_sub.max()))
            subbox_ID = np.zeros((N_partside, N_partside, N_partside))
            subbox_pos = np.zeros((3, N_partside, N_partside, N_partside))
            subbox_vel = np.zeros((3, N_partside, N_partside, N_partside))
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
                        subbox_pos[0,:,j_y,j_z] = x_subbox[ylim_sub & zlim_sub][j_x_sorted]
                        subbox_pos[1,:,j_y,j_z] = y_subbox[ylim_sub & zlim_sub][j_x_sorted]
                        subbox_pos[2,:,j_y,j_z] = z_subbox[ylim_sub & zlim_sub][j_x_sorted]
                        subbox_vel[0,:,j_y,j_z] = vx_subbox[ylim_sub & zlim_sub][j_x_sorted]
                        subbox_vel[1,:,j_y,j_z] = vy_subbox[ylim_sub & zlim_sub][j_x_sorted]
                        subbox_vel[2,:,j_y,j_z] = vz_subbox[ylim_sub & zlim_sub][j_x_sorted]
            subbox['meta'] = fullbox['meta'].copy() 
            subbox['meta']['subbox_ijk'] = [i_x, i_y, i_z]
            subbox['ID'] = subbox_ID
            subbox['Position'] = subbox_pos
            subbox['Velocity'] = subbox_vel

            fsub = h5py.File(F(isubbox), 'w') # save ID's to hdf5 file 
            for k in subbox['meta'].keys(): # store meta data 
                fsub.attrs.create(k, subbox['meta'][k]) 
            fsub.create_dataset('ID', data=subbox['ID']) 
            fsub.create_dataset('Position', data=subbox['Position']) 
            fsub.create_dataset('Velocity', data=subbox['Velocity']) 
            fsub.close() 
        subboxes.append(subbox) 
    if Nsub == 1: 
        return subboxes[0]
    else: 
        return subboxes


def mNuParticles(mneut, nreal, nzbin, sim='paco', overwrite=False, verbose=False): 
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
    _dir = ''.join([UT.mNuDir(mneut, sim=sim), str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.hdf5']) # snapshot 

    particle_data = {} 
    if os.path.isfile(f) and not overwrite: # read in the file 
        f_par = h5py.File(f, 'r') # read in hdf5 file 
        particle_data['meta'] = {} 
        for k in f_par.attrs.keys(): 
            particle_data['meta'][k] = f_par.attrs[k]
        # read in datasets
        for k in f_par.keys(): 
            particle_data[k] = f_par[k].value 
        f_par.close() 
    else: # write file 
        f_gadget = ''.join([_dir, 'snap_', str(nzbin).zfill(3)]) # snapshot 
        # read in Gadget header
        header = ReadSnap.read_gadget_header(f_gadget)

        # read in CDM particles (parttype = 1) and create catalogue
        read_keys = ['POS ', 'VEL ', 'ID  ', 'MASS'] # currently only reading POS, VEL, ID, and MASS
        save_keys = ['Position', 'Velocity', 'ID', 'Mass'] 
        
        particle_data['meta'] = header # store meta data
        for k_s, k_r in zip(save_keys, read_keys): 
            particle_data[k_s] = ReadSnap.read_block(f_gadget, k_r, parttype=1, verbose=verbose)
            if k_s == 'Position': 
                particle_data[k_s] /= 1000. # convert ot Mpc/h 
        
        f_par = h5py.File(f, 'w') # save ID's to hdf5 file 
        for k in particle_data.keys(): 
            if k == 'meta': continue 
            f_par.create_dataset(k, data=particle_data[k]) 

        for k in particle_data['meta'].keys(): # store meta data 
            f_par.attrs.create(k, particle_data['meta'][k]) 
        f_par.close() 
    return particle_data 


def mNuICs(nreal, sim='paco', overwrite=False, verbose=False): 
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
    _dir = ''.join([UT.mNuDir(0.0, sim=sim), str(nreal), '/ICs/'])
    if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
    f = ''.join([_dir, 'ics.hdf5']) # snapshot 
    
    particle_data = {} 
    if os.path.isfile(f) and not overwrite: # read in the file 
        f_ics = h5py.File(f, 'r') # read in hdf5 file 
        particle_data['meta'] = {} 
        for k in f_ics.attrs.keys(): 
            particle_data['meta'][k] = f_ics.attrs[k]
        # read in datasets
        for k in f_ics.keys(): 
            particle_data[k] = f_ics[k].value 
        f_ics.close() 
    else: # write file 
        f_gadget = ''.join([_dir, 'ics']) # read in Gadget header
        header = ReadSnap.read_gadget_header(f_gadget)

        # read in CDM particles (parttype = 1) and create catalogue
        read_keys = ['POS ', 'VEL ', 'ID  ', 'MASS'] # currently only reading POS, VEL, ID, and MASS
        save_keys = ['Position', 'Velocity', 'ID', 'Mass'] 
        
        particle_data['meta'] = header # store meta data
        for k_s, k_r in zip(save_keys, read_keys): 
            particle_data[k_s] = ReadSnap.read_block(f_gadget, k_r, parttype=1, verbose=verbose)
            if k_s == 'Position': 
                particle_data[k_s] /= 1000. # convert ot Mpc/h 

        f_ics = h5py.File(f, 'w') # save ID's to hdf5 file 
        for k in particle_data.keys(): 
            if k == 'meta': continue 
            f_ics.create_dataset(k, data=particle_data[k]) 

        for k in particle_data['meta'].keys(): # store meta data 
            f_ics.attrs.create(k, particle_data['meta'][k]) 
        f_ics.close() 
    return particle_data 


"""
    def mNuDispField(mneut, nreal, nzbin, boundary_correct=True, sim='paco', overwrite=False, verbose=False): 
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
            print('Constructing %s ......' % f) 
            # read in particles of snapshot 
            par = mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=verbose)
            # read in initial conditions 
            ics = mNuICs(nreal, sim=sim, verbose=verbose) 
            assert len(ics['ID']) == len(par['ID']) 

            isort_ics = np.argsort(ics['ID']) 
            isort_par = np.argsort(par['ID']) 
            #assert np.array_equal(ics['ID'][isort_ics], par['ID'][isort_par])

            dispfield = {} 
            dispfield['meta'] = par['meta'].copy() 
            dispfield['ID'] = ics['ID'][isort_ics]
            # position @ snapshot - position @ initial condition
            dispfield['dispfield'] = par['Position'][isort_par] - ics['Position'][isort_ics]
            
            fdisp = h5py.File(f, 'w') # save ID's to hdf5 file 
            for k in dispfield.keys(): 
                if k == 'meta': continue 
                fdisp.create_dataset(k, data=dispfield[k]) 

            for k in par['meta'].keys(): # store meta data 
                fdisp.attrs.create(k, dispfield['meta'][k]) 
            fdisp.close() 

        if boundary_correct: 
            # correct for the cases where particles cross the periodic boundaries 
            # the correction at the moment is brute forced. 
            for i in range(3): 
                cross_neg = (dispfield['dispfield'][:,i] < -900.) 
                dispfield['dispfield'][cross_neg,i] += 1000.
                
                cross_pos = (dispfield['dispfield'][:,i] > 900.) 
                dispfield['dispfield'][cross_pos,i] = 1000. - dispfield['dispfield'][cross_pos,i]
        return dispfield 
"""
