#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -C haswell
#SBATCH -J _part_subbox_nreal5
#SBATCH -o _part_subbox_nreal5.o 
#SBATCH -L SCRATCH,project 


now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nreal=5 
nzbin=2 
nside=8
n_cpu=32

#srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "ics_mp" $nreal $nside 0 1 $n_cpu

# particles 
for mneut in 0.0eV 0.06eV 0.10eV 0.15eV 0.6eV
do 
    echo "neutrino mass = "$mneut
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside 0 31 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
