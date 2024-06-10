#!/bin/bash

# Initialize total number of atoms
total_atoms=54
elements=("V" "Cr" "Ba" "W")
increment=3

# Function to create alloy with given composition
create_alloy () {
    local -a composition=("$@")
    local filename="alloy_${composition[*]}.xsf"

    # Check if the composition values are integers
    for value in "${composition[@]}"; do
        if ! [[ "$value" =~ ^[0-9]+$ ]]; then
            echo "Error: Composition values must be integers."
            exit 1
        fi
    done

    # Create an initial supercell with just the first element
    atomsk --create bcc 2.87 "${elements[0]}" -duplicate 3 3 3 supercell.xsf

    # Loop over the elements and their respective counts
    for i in "${!elements[@]}"; do
        if [ "${composition[i]}" -gt 0 ]; then
            # Randomly replace some atoms with the current element
            atomsk supercell.xsf -sub "${elements[0]}" "${elements[i]}" "${composition[i]}" "$filename"
            mv "$filename" supercell.xsf
        fi
    done

    # Check the final number of atoms
    final_atoms=$(atomsk supercell.xsf -count | grep "Total number of atoms" | awk '{print $6}')
    if [ "$final_atoms" -ne "$total_atoms" ]; then
        echo "Error: The final structure has $final_atoms atoms, but it should have $total_atoms atoms."
        exit 1
    fi

    # Output the final structure in the XSF format
    atomsk supercell.xsf xsf "$filename"
    echo "Alloy structure created successfully with $total_atoms atoms."
}

# Loop over each element and vary its count in increments of increment
for (( i=0; i<=$total_atoms; i+=$increment )); do
    for (( j=0; j<=$total_atoms-i; j+=$increment )); do
        for (( k=0; k<=$total_atoms-i-j; k+=$increment )); do
            # Calculate the remaining atoms for the last element
            l=$((total_atoms - i - j - k))
            if [ $l -ge 0 ] && [ $((i+j+k+l)) -eq $total_atoms ]; then
                create_alloy $i $j $k $l
            fi
        done
    done
done
