#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 01:30:00 
#SBATCH -J subbox_pars_0.6eV
#SBATCH -o subbox_pars_0.6eV.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nreal=1 
nzbin=2 
nside=8

# ICs
#python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 500 511 

# particles 
mneut=0.6eV
echo "neutrino mass = "$mneut
python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside  0 511 


now=$(date +"%T") 
echo "end time ... $now"
