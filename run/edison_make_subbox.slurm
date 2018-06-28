#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 01:00:00 
#SBATCH -J _ic_subbox_nreal5
#SBATCH -o _ic_subbox_nreal5.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nreal=5 
nzbin=2 
nside=8

# ICs
n_cpu=24
srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 0 511 $n_cpu

# particles 
#mneut=0.6eV
#echo "neutrino mass = "$mneut
#python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside  0 511 


now=$(date +"%T") 
echo "end time ... $now"
