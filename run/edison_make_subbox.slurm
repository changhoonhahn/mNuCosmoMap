#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -J part_nside16_test
#SBATCH -o _part_nside16_test.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nzbin=2 
nside=16
n_cpu=8
mneut=0.0eV

# particles 
for nreal in {11..11}; do 
    echo "neutrino mass = "$mneut
    # each 512 takes ~6mins 
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside 0 511 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
