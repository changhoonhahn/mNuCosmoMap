#!/bin/bash -l 
#SBATCH -q debug 
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -J part_nside32_nreal1
#SBATCH -o _part_nside32_nreal1.o 
#SBATCH -L SCRATCH,project 

now=$(date +"%T") 
echo "start time ... $now"

module load python/2.7-anaconda
source activate myenv0 

nzbin=2 
nside=32
n_cpu=24

nreal=1

# particles 
#for mneut in 0.0eV; do # 0.06eV 0.10eV 0.15eV 0.6eV; do 
for mneut in 0.06eV 0.10eV 0.15eV 0.6eV; do 
    echo "neutrino mass = "$mneut
    srun -n 1 -c $n_cpu python /global/homes/c/chahah/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside 0 511 $n_cpu
done 

now=$(date +"%T") 
echo "end time ... $now"
