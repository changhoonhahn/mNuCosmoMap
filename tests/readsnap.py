'''

Some tests of the readsnap.py module 

'''
import numpy as np 
from mnucosmomap import util as UT 
from mnucosmomap import readsnap as ReadSnap

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


def ParticleTypes_ID(): 
    ''' parttype specifies the type of particle. Try to figure out
    the ordering of the IDs for the different particle types
    '''
    # read in Gadget header
    nzbin = 2
    _dir = ''.join([UT.dat_dir(), 'sims/paco/0.06eV/1/snapdir_', str(nzbin).zfill(3), '/'])
    f_gadget = ''.join([_dir, 'snap_', str(nzbin).zfill(3)]) # snapshot 
    header = ReadSnap.read_gadget_header(f_gadget)
    print('----------------') 
    for k in header.keys(): 
        print(k, header[k])
    print('----------------') 

    # read in CDM particles (parttype = 1) and create catalogue
    read_keys = ['ID  ']#['POS ', 'VEL ', 'ID  ', 'MASS'] # currently only reading POS, VEL, ID, and MASS
    save_keys = ['ID']#['Position', 'Velocity', 'ID', 'Mass'] 
    
    # CDM particles
    particle_data = {} 
    for k_s, k_r in zip(save_keys, read_keys): 
        particle_data[k_s] = ReadSnap.read_block(f_gadget, k_r, parttype=1, verbose=False)
        if k_s == 'Position': 
            particle_data[k_s] /= 1000. # convert ot Mpc/h 
    print('%i CDM particle IDs' % len(particle_data['ID']))
    isort = np.argsort(particle_data['ID'])
    print(particle_data['ID'][isort][:15]) 
    print('----------------') 

    # neutrino particles
    particle_data = {} 
    for k_s, k_r in zip(save_keys, read_keys): 
        particle_data[k_s] = ReadSnap.read_block(f_gadget, k_r, parttype=2, verbose=False)
        if k_s == 'Position': 
            particle_data[k_s] /= 1000. # convert ot Mpc/h 
    print('%i Nu particle IDs' % len(particle_data['ID']))
    isort = np.argsort(particle_data['ID'])
    print(particle_data['ID'][isort][:15]) 
    return None 


if __name__=="__main__": 
    ParticleTypes_ID()
