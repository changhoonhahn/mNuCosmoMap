#!bin/bash/
source_ep='8ecb2e28-da42-11e5-976d-22000b9da45e'
dest_ep='1c93ed5e-b8da-11e7-b143-22000a92523b'

for mneut in 0.0eV 0.06eV 0.10eV 0.15eV 0.6eV
do 
    for snap in 1 #{1..5}
    do 
        dir_local=$MNUCOSMOMAP_DIR"sims/paco/"$mneut/$snap/snapdir_002/
        if [ ! -d $dir_local ]; then
            echo $dir_local 
            mkdir -p $dir_local
        fi 
        if [ ! -f $dir_local"snap_002.0" ]; then 
            echo $dir_local"snap_002.0 does not exist!"
            dir_source=/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/$mneut/$snap/snapdir_002/
            globus transfer $source_ep:$dir_source $dest_ep:$dir_local --recursive
        fi 
    done 
done 
