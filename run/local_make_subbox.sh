#!/bin/bash
nreal=1 
nzbin=2
nside=8

# ICs 
#python /Users/chang/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside 0 99 

# particles 
for mneut in 0.0eV; do #0.06eV 0.10eV 0.15eV 0.6eV; do 
    echo "neutrino mass = "$mneut
    python /Users/chang/projects/mNuCosmoMap/run/make_subbox.py "particles" $mneut $nreal $nzbin $nside  0 5 
done 
