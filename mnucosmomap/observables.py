'''

code for calculating observables for mNuCosmoMap.catalogs 
objects

author(s) : ChangHoon Hahn 

'''
import numpy as np 
# -- nbodykit -- 
from nbodykit.base.catalog import CatalogSource
from nbodykit.algorithms.fftpower import FFTPower
# -- local -- 
from mnucosmomap import util as UT
from mnucosmomap import catalogs as mNuCat


def Pk(pos, mode='1d', Nmesh=None, BoxSize=None,  **kwargs): 
    ''' Wrapper for nbodykit.algorithms.fftpower.FFTPower in the 
    nbodykit package. Given xyz positions of objects in a **periodic
    box**, this code will calculate the powerspectrum. 
    
    pos : (3, N_particles)
    '''
    n_part = pos.shape[1] # number of particles
    # generate CatalogSource object  
    cat = CatalogSource()
    cat._size = n_part
    cat._csize = cat.comm.allreduce(cat.size)
    cat['Position'] = pos.T
    # measure powerspectrum
    pique = FFTPower(cat, mode, Nmesh=Nmesh, BoxSize=BoxSize, **kwargs)
    return pique.power
