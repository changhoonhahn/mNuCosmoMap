#!/bin/bash
echo "password" 
read -s pwd 

edison_dir=/scratch2/scratchdirs/chahah/mNuCosmoMap/sims/paco
local_dir=/Volumes/chang_eHDD/projects/mNuCosmoMap/sims/paco

for mneut in 0.0eV 0.06eV 0.10eV 0.15eV 0.6eV
do
    for snap in 1
    do
        # check that the neutrino directory exists
        mneut_dir=$edison_dir/$mneut/$snap/snapdir_002/
        if sshpass -p $pwd ssh edison '[ ! -d '$mneut_dir' ]'; then
            sshpass -p $pwd ssh edison "mkdir -p "$mneut_dir
        fi
        for isub in {0..9}; do
            subbox_file=$mneut_dir"snap_002.nside8."$isub".hdf5"
            if sshpass -p $pwd ssh edison '[ ! -f '$subbox_file' ]'; then
                echo "transfering ..."$subbox_file
                sshpass -p $pwd scp $local_dir/$mneut/$snap/snapdir_002/snap_002.nside8."$isub".hdf5 edison:$mneut_dir
            fi
        done 
    done 

done 
