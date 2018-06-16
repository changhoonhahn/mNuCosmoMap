#!/bin/bash
nreal=1 
nzbin=2
nside=8
nsubbox=10

    
python /Users/chang/projects/mNuCosmoMap/run/make_subbox.py "ics" $nreal $nside

#for mneut in 0.0eV 0.06eV 0.10eV 0.15eV 0.6eV; do 
#    echo "neutrino mass = "$mneut
#    python /Users/chang/projects/mNuCosmoMap/run/make_subbox.py "ics" $mneut $nreal $nzbin $nside $nsubbox
#done 
