import os
import sys 
import time
import h5py
import numpy as np 
import multiprocessing as MP 
from mnucosmomap import util as UT 
from mnucosmomap import catalogs as mNuCat 


if __name__=="__main__": 
    tt = sys.argv[1] 
    if tt == 'ics': 
        nreal = int(sys.argv[2])
        nside = int(sys.argv[3]) 
        nsubbox = range(int(sys.argv[4]), int(sys.argv[5])+1)
        n_cpu = int(sys.argv[6]) 
        sim = 'paco' 
        
        #mNuICs_subboxes(range(nsubbox0, nsubbox1+1), nreal, sim='paco', nside=nside, n_cpu=n_cpu)
        t00 = time.time() 
        assert np.array(nsubbox).max() <= nside**3 

        _dir = ''.join([UT.mNuDir(0.0, sim=sim), str(nreal), '/ICs/'])
        if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
        F = lambda isub: ''.join([_dir, 'ics.nside', str(nside), '.', str(isub), '.hdf5'])

        # read in full IC box. 
        print('Reading in full IC box...') 
        t0 = time.time() 
        fullbox = mNuCat.mNuICs(nreal, sim=sim, verbose=True)
        print('full box read in takes: %f mins' % ((time.time() - t0)/60.))
        # append extra metadata
        fullbox['meta']['n_side'] = nside
        fullbox['meta']['n_subbox'] = nside**3

        id_full = fullbox['ID']
        x, y, z = fullbox['Position'].T
        vx, vy, vz = fullbox['Velocity'].T
        
        L_subbox = 1000./float(nside) # L_subbox
        L_res = 1000./512.
        L_halfres = 0.5 * L_res
        N_partside = 512/nside
        N_subbox = (N_partside)**3

        def _make_subbox(isubbox): 
            t0 = time.time()  
            i_x = ((isubbox % nside**2) % nside) 
            i_y = ((isubbox % nside**2) // nside) 
            i_z = (isubbox // nside**2) 
            print('%i, %i, %i' % (i_x, i_y, i_z))

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
            
            ID_sub = id_full[in_subbox] #fullbox['ID'][in_subbox]
            x_subbox = x[in_subbox]
            y_subbox = y[in_subbox]
            z_subbox = z[in_subbox]
            x_sub = (x_subbox - i_x * L_subbox) % 1000.
            y_sub = (y_subbox - i_y * L_subbox) % 1000.
            z_sub = (z_subbox - i_z * L_subbox) % 1000.

            vx_subbox = vx[in_subbox]
            vy_subbox = vy[in_subbox]
            vz_subbox = vz[in_subbox]
            print('%f < x_sub < %f' % (x_sub.min(), x_sub.max()))
            print('%f < y_sub < %f' % (y_sub.min(), y_sub.max()))
            print('%f < z_sub < %f' % (z_sub.min(), z_sub.max()))
            subbox_ID = np.zeros((N_partside, N_partside, N_partside))
            subbox_pos = np.zeros((3, N_partside, N_partside, N_partside))
            subbox_vel = np.zeros((3, N_partside, N_partside, N_partside))

            j_x = ((x_sub - L_halfres) // L_res).astype(int) 
            j_y = ((y_sub - L_halfres) // L_res).astype(int) 
            j_z = ((z_sub - L_halfres) // L_res).astype(int) 
            subbox_ID[j_x,j_y,j_z] = ID_sub
            subbox_pos[0,j_x,j_y,j_z] = x_subbox
            subbox_pos[1,j_x,j_y,j_z] = y_subbox
            subbox_pos[2,j_x,j_y,j_z] = z_subbox
            subbox_vel[0,j_x,j_y,j_z] = vx_subbox
            subbox_vel[1,j_x,j_y,j_z] = vy_subbox
            subbox_vel[2,j_x,j_y,j_z] = vz_subbox
            
            subbox = {} 
            subbox['meta'] = fullbox['meta'].copy()
            subbox['meta']['i_subbox'] = isubbox
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
            print('making subbox %i took: %f mins' % (isubbox, (time.time() - t0)/60.))
            return None 

        print('Constructing subboxes using %i processes...' % n_cpu)
        pewl = MP.Pool(processes=n_cpu)
        pewl.map(_make_subbox, [(isubbox) for isubbox in nsubbox]) 
        pewl.close()
        pewl.terminate()
        pewl.join() 
        print('Everything took %f mins' % ((time.time() - t00)/60.))
    elif tt == 'particles': 
        mneut = float(sys.argv[2].strip('eV')) 
        nreal = int(sys.argv[3]) 
        nzbin = int(sys.argv[4]) 
        nside = int(sys.argv[5]) 
        nsubbox = range(int(sys.argv[6]), int(sys.argv[7])+1)
        sim = 'paco'
        t00 = time.time() 
    
        _dir = ''.join([UT.mNuDir(mneut, sim=sim), str(nreal), '/snapdir_', str(nzbin).zfill(3), '/'])
        if not os.path.isdir(_dir): raise ValueError("directory %s not found" % _dir)
        # snapshot subbox file 
        F = lambda isub: ''.join([_dir, 'snap_', str(nzbin).zfill(3), '.nside', str(nside), '.', str(isub), '.hdf5']) 
                
        # read in full box of Particles
        t0 = time.time() 
        fullbox = mNuCat.mNuParticles(mneut, nreal, nzbin, sim=sim, verbose=True)
        print('full box read in takes: %f mins' % ((time.time() - t0)/60.))
        isort_box = np.argsort(fullbox['ID'])

        def _make_particle_subbox(isubbox): 
            t0 = time.time()  
            # read in subbox indicies
            subb = mNuCat.mNuICs_subbox(isubbox, nreal, nside=nside, sim=sim, verbose=True)
            N_partside, N_partside, N_partside = subb['ID'].shape
            sub_id = subb['ID'].flatten()

            subbox['meta'] = fullbox['meta'].copy() 
            subbox['meta']['n_side'] = nside # append extra metadata
            subbox['meta']['n_subbox'] = nside**3 
            subbox['meta']['i_subbox'] = isubbox 
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
            print('making subbox %i took: %f mins' % (isubbox, (time.time() - t0)/60.))
            return None 

        print('Constructing subboxes using %i processes...' % n_cpu)
        pewl = MP.Pool(processes=n_cpu)
        pewl.map(_make_particle_subbox, [(isubbox) for isubbox in nsubbox]) 
        pewl.close()
        pewl.terminate()
        pewl.join() 
        print('Everything took %f mins' % ((time.time() - t00)/60.))

    # randomly check a few subboxes to make sure. 
    #i_rand = np.random.choice(range(nsubbox), size=3, replace=False) 
    #for i in i_rand: 
    #    _check_subbox(mneut, nreal, i, nzbin=nzbin, nside=nside)
