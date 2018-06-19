#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 03:00:00 
#SBATCH -J subbox_400_512
#SBATCH -o subbox_400_512.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nreal=1 
nzbin=2 
nside=8
    
python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 400 511

now=$(date +"%T") 
echo "end time ... $now"
