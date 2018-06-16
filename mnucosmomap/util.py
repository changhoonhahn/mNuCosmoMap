'''

General utility functions 

'''
import os
import numpy as np


def check_env(): 
    if os.environ.get('MNUCOSMOMAP_DIR') is None: 
        raise ValueError("set $MNUCOSMOMAP_DIR environment varaible!") 
    if os.environ.get('MNUCOSMOMAP_CODEDIR') is None: 
        raise ValueError("set $MNUCOSMOMAP_DIR environment varaible!") 
    return None


def mNuDir(mneut, sim='paco'): 
    ''' Directory where massive neutrino simulations are 
    '''
    if mneut == 0.1: 
        str_mneut = '0.10'
    else: 
        str_mneut = str(mneut) 
    return ''.join([dat_dir(), 'sims/', sim, '/', str_mneut, 'eV/'])


def dat_dir(): 
    ''' directory that contains all the data files, defined by environment 
    variable $IQUENCH_DIR
    '''
    return os.environ.get('MNUCOSMOMAP_DIR') 


def code_dir(): 
    if os.environ.get('MNUCOSMOMAP_CODEDIR') is None: 
        raise ValueError("set $MNUCOSMOMAP_CODEDIR environment varaible!") 
    return os.environ.get('MNUCOSMOMAP_CODEDIR') 
    

def fig_dir(): 
    ''' directory to dump all the figure files 
    '''
    if os.path.isdir(code_dir()+'figs/'):
        return code_dir()+'figs/'
    else: 
        raise ValueError("create figs/ folder in $MNUCOSMOMAP_CODEDIR directory for figures")


def doc_dir(): 
    ''' directory for paper related stuff 
    '''
    if os.path.isdir(code_dir()+'doc/'):
        return code_dir()+'doc/'
    else: 
        raise ValueError("create doc/ folder in $MNUCOSMOMAP_CODEDIR directory for documntation")
