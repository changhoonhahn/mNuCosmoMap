#!/bin/bash/
# transfer snapdir_002 from Paco to Nersc using Globus
source_ep='8ecb2e28-da42-11e5-976d-22000b9da45e'
dest_ep='9d6d994a-6d04-11e5-ba46-22000b92c6ec'
   
nersc_dir="/global/cscratch1/sd/chahah/mNuCosmoMap/"

echo "password:" 
read -s pwd 

#0.0eV 
#for mneut in 0.06eV 0.10eV 0.15eV 0.6eV; do 
#done 

mneut=0.6eV
for snap in {7..100}; do 
    dir_nersc=$nersc_dir"sims/paco/"$mneut"/"$snap"/snapdir_002/"
    dir_cca="/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/"$mneut"/"$snap"/snapdir_002/"

    if sshpass -p $pwd ssh edison '[ ! -d '$dir_nersc' ]'; then
        echo "transfering ... "$dir_nersc
        globus transfer $source_ep:$dir_cca $dest_ep:$dir_nersc --recursive
    fi
done 
