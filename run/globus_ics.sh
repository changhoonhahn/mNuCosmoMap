#!/bin/bash/
# transfer ICs from Paco
source_ep='8ecb2e28-da42-11e5-976d-22000b9da45e'
dest_ep='9d6d994a-6d04-11e5-ba46-22000b92c6ec'
   
nersc_dir="/global/cscratch1/sd/chahah/mNuCosmoMap/"

for snap in {21..100}; do 
    dir_nersc=$nersc_dir"sims/paco/0.0eV/"$snap"/ICs/"
    dir_cca="/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.0eV/"$snap"/ICs/"
    globus transfer $source_ep:$dir_cca $dest_ep:$dir_nersc --recursive
done 
