#!/bin/bash -l 
#SBATCH -q debug 
#SBATCH -N 1 
#SBATCH -t 00:05:00 
#SBATCH -C haswell
#SBATCH -J _make_subbox_mptest
#SBATCH -o _make_subbox_mptest.o 
#SBATCH -L SCRATCH,project 


now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nreal=1 
nzbin=2 
nside=8

n_cpu=2
srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics_mp" $nreal $nside 0 1 $n_cpu

now=$(date +"%T") 
echo "end time ... $now"
