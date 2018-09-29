#!/bin/bash -l 
#SBATCH -q debug 
#SBATCH -N 1 
#SBATCH -t 00:10:00 
#SBATCH -J ics_nside32_nreal1
#SBATCH -o _ics_nside32_nreal1.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nzbin=2 
nside=32
n_cpu=24

nreal=1 

#for nreal in {1..100}; do 
# generate ICs subboxes ~100 takes 1 mins 
srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 0 511 $n_cpu
#done 

now=$(date +"%T") 
echo "end time ... $now"
