#!/bin/zsh

#NX=8
#NY=1

OUTFILE="energy_XXZ.txt"

echo "Delta\t\tEnergy" > "$OUTFILE"

# looping over Delta values

for Delta in $(seq -1.5 0.05 1.5); do
    #echo -e "Nx=$NX\nNy=$NY\nNsites=Nx*Ny\nDelta=$Delta" > parameter.py
    echo "$Delta"
    sed -i '' "s/^Delta=.*/Delta=$Delta/" parameter.py
    python3 XXZ_energy.py

done 

echo "Completed!"