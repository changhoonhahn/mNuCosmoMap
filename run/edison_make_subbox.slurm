#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -J _part_subbox_nreal2
#SBATCH -o _part_subbox_nreal2.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nreal=5 
nzbin=2 
nside=8
n_cpu=24

# ICs
#srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 0 511 $n_cpu

# particles 
for menut in 0.0eV 0.06eV 0.10eV 0.15eV 0.6eV; 
do 
    echo "neutrino mass = "$mneut
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 0 511 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
