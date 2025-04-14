#!/bin/zsh

OUTFILE="energy_XXZ_1D_12_JW_parity.txt"

echo "Delta\t\tEnergy" > "$OUTFILE"

# looping over Delta values

for Delta in $(seq -2.0 0.05 2.0); do
    #echo -e "Nx=$NX\nNy=$NY\nNsites=Nx*Ny\nDelta=$Delta" > parameter.py
    echo "$Delta"
    sed -i '' "s/^ *Delta *=.*/Delta=$Delta/" parameter.py  # * ensures if there is space in the parameter.py
    python3 XXZ_energy.py

done 

echo "Completed!"