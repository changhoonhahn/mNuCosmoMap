#!/bin/bash -l 
#SBATCH -q debug 
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -J ics_nside16_nreal1_5
#SBATCH -o _ics_nside16_nreal1_5.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nzbin=2 
n_cpu=24
for nreal in {1..5}; do 
    # generate ICs subboxes 512 takes ~2.5 mins 
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal 16 0 511 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
