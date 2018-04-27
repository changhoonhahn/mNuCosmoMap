'''

Tests for mnucosmomap.catalogs


'''
from mnucosmomap import catalogs as mNuCat 


def _mNuParticles(): 
    '''check that the function mnucosmomap.catalogs.mNuParticles is 
    properly reading in the snapshot particles 
    '''
    mneut = 0.0 
    nreal = 1 
    nzbin = 2 # z=1 
    sim = 'paco' # only paco is supported anyway. 
    
    # read in the particles
    cat = mNuCat.mNuParticles(mneut, nreal, nzbin, sim=sim)
    
    print(cat['Position'].shape)
    # now run some common sense checks 
    #fig = plt.figure()
    #sub = fig.add_subplot(111)
    return None 


if __name__=="__main__": 
    _mNuParticles()
