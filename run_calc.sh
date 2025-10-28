#!/bin/zsh

OUTFILE="energy_XXZ_1D_12_JW_parity_OBC.txt"

# Create/overwrite with header 
echo -e "Delta\tEnergy" > "$OUTFILE"

# Sweep Delta from -1.0 to 1.5 in steps of 0.05
for Delta in $(seq -1.0 0.05 1.5); do
  echo "Running Delta = $Delta"

  # Update Delta in parameter.py (handles optional leading spaces)
  sed -i '' "s/^ *Delta *=.*/Delta = $Delta/" parameter.py

  # Run and extract the energy from the specific line
  ENERGY=$($(which python) Main.py | awk '/^The Energy for/ {print $NF; exit}')

  # Append Delta and Energy to the outfile
  echo -e "$Delta\t$ENERGY" >> "$OUTFILE"
done

echo "Completed! Results saved to $OUTFILE"
