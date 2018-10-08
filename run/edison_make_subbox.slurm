#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 01:00:00 
#SBATCH -J part_nside32_0.6eV_nreal1_10
#SBATCH -o _part_nside32_0.6eV_nreal1_10.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nzbin=2 
nside=32
n_cpu=24
mneut=0.6eV

# particles 
for nreal in {1..10}; do 
    echo "neutrino mass = "$mneut
    # each 512 takes ~6mins 
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside 0 511 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
