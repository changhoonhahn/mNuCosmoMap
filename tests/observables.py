import time
import h5py
import numpy as np 
from mnucosmomap import util as UT 
from mnucosmomap import catalogs as mNuCat 
from mnucosmomap import observables as mNuObvs 
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

def Pk_fullbox(): 
    mneuts = [0.0]#, 0.06, 0.10, 0.15, 0.6]
    fig = plt.figure(figsize=(6*len(mneuts), 6))
    for mneut in mneuts: 
        box = mNuCat.mNuParticles(mneut, 1, 2) 
        pos = box['Position'].T
        pk = mNuObvs.Pk(pos, Nmesh=256, BoxSize=1e3) 
        for k in pk.attrs: 
            print("%s = %s" %(k, str(pk.attrs[k])))

        sub = fig.add_subplot(111)
        sub.plot(pk['k'], pk['power']) 
        sub.set_xlabel('$k$', fontsize=25)
        sub.set_xscale('log') 
        sub.set_xlim([8e-3,0.5]) 
        sub.set_ylabel('$P(k)$', fontsize=25)
        sub.set_ylim([1e1, 4e4]) 
        sub.set_yscale('log') 
        sub.set_title(r'$\sum m_\nu = '+str(round(mneut,2))+'$', fontsize=25)
    fig.savefig(''.join([UT.fig_dir(), '_box_pos_pk.png']), bbox_inches='tight') 
    return None 


def Pk_subbox(): 
    mneuts = [0.0]#, 0.06, 0.10, 0.15, 0.6]
    fig = plt.figure(figsize=(6*len(mneuts), 6))
    for i_nu, mneut in enumerate(mneuts): 

        box = mNuCat.mNuParticles(mneut, 1, 2) 
        pk_full = mNuObvs.Pk(box['Position'].T, Nmesh=256, BoxSize=1e3) 

        subs = mNuCat.mNuParticles_subbox(range(10), mneut, 1, 2, nside=8) 
        pks = [] 
        for sub in subs: 
            pos = np.array([
                sub['Position'][0,:,:,:].flatten(), 
                sub['Position'][1,:,:,:].flatten(), 
                sub['Position'][2,:,:,:].flatten()]) 
            pk = mNuObvs.Pk(pos, Nmesh=64, BoxSize=125) 
            pks.append(pk) 

        sub = fig.add_subplot(1,len(mneuts),1+i_nu)
        sub.plot(pk_full['k'], pk_full['power'], c='k', ls='--', label='Full box') 
        for i, pk in enumerate(pks): 
            sub.plot(pk['k'], pk['power'], label=('subbox %i' % i)) 
        sub.legend(loc='upper right', prop={'size':15}) 
        sub.set_xlabel('$k$', fontsize=25)
        sub.set_xscale('log') 
        #sub.set_xlim([1e-2,5.0]) 
        sub.set_xlim([8e-3,2.0]) 
        sub.set_ylabel('$P(k)$', fontsize=25)
        sub.set_yscale('log') 
        sub.set_title(r'$\sum m_\nu = '+str(round(mneut,2))+'$', fontsize=25)
    fig.savefig(''.join([UT.fig_dir(), '_subbox_pos_pk.png']), bbox_inches='tight') 
    return None 


if __name__=="__main__": 
    #Pk_fullbox() 
    Pk_subbox()
