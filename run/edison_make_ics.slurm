#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 01:30:00 
#SBATCH -J ics_11_100
#SBATCH -o _ics_11_100.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nzbin=2 
n_cpu=24
for nreal in {20..100}; do 
    # generate ICs subboxes ~100 takes 1 mins 
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal 16 0 0 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
