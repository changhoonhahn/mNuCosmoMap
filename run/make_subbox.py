import os
import sys 
import time
import h5py
import numpy as np 
import multiprocessing as MP 
from mnucosmomap import util as UT 
from mnucosmomap import catalogs as mNuCat 


def mNuParticles_subbox(nsubbox, mneut, nreal, nzbin=2, sim='paco', nside=8): 
    ''' Generate suboxes for given mneut, nreal, nzbin, sim, and nside using 
    mnucosmomap.catalogs.mNuParticles_subbox  
    '''
    sb = mNuCat.mNuParticles_subbox(nsubbox, mneut, nreal, nzbin, 
            sim=sim, nside=nside, overwrite=True, verbose=False)
    return None 


"""
    def mNuICs_subbox(nsubbox, nreal, sim='paco', nside=8): 
        ''' Generate subboxes of ICs 
        '''
        sb = mNuCat.mNuICs_subbox(nsubbox, nreal, sim=sim, nside=nside, 
                overwrite=False, verbose=True)
        return None

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
"""


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
        nsubbox0 = int(sys.argv[6])
        nsubbox1 = int(sys.argv[7])
        mNuParticles_subbox(range(nsubbox0, nsubbox1+1), mneut, nreal, nzbin=nzbin, sim='paco', nside=nside)
    # randomly check a few subboxes to make sure. 
    #i_rand = np.random.choice(range(nsubbox), size=3, replace=False) 
    #for i in i_rand: 
    #    _check_subbox(mneut, nreal, i, nzbin=nzbin, nside=nside)
