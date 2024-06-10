import os
import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.io import write

# Function to create alloys using ASE
def create_alloy(composition, elements, total_atoms, filename):
    # Create a base structure (e.g., BCC V)
    base_element = elements[0]
    base_structure = bulk(base_element, 'bcc', a=2.87)

    # Create a supercell
    supercell = make_supercell(base_structure, [[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    supercell *= int(total_atoms / len(supercell))

    # Ensure the supercell has the correct number of atoms
    while len(supercell) > total_atoms:
        supercell.pop(-1)

    # Replace atoms to create the alloy
    indices = list(range(len(supercell)))
    np.random.shuffle(indices)
    
    current_index = 0
    for i, count in enumerate(composition):
        for _ in range(count):
            if current_index < len(indices):
                supercell[indices[current_index]].symbol = elements[i]
                current_index += 1

    # Write the structure to a file
    write(filename, supercell)

# Function to generate alloy compositions with one element varied and others equal
def generate_alloy_compositions(elements, total_atoms, increment):
    compositions = []
    num_elements = len(elements)
    
    # Generate compositions with one element varied and others equal
    for varying_element in range(num_elements):
        for i in range(0, total_atoms + 1, increment):
            remaining_atoms = total_atoms - i
            if remaining_atoms % (num_elements - 1) == 0:
                equal_share = remaining_atoms // (num_elements - 1)
                composition = [equal_share] * num_elements
                composition[varying_element] = i
                compositions.append(composition)
    
    # Generate any remaining combinations
    for i in range(0, total_atoms + 1, increment):
        for j in range(0, total_atoms - i + 1, increment):
            for k in range(0, total_atoms - i - j + 1, increment):
                l = total_atoms - i - j - k
                if l >= 0:
                    composition = [i, j, k, l]
                    if composition not in compositions:
                        compositions.append(composition)
    
    return compositions

# Isolated part to generate alloy files
def generate_alloy_files():
    total_atoms = 54
    elements = ["V", "Cr", "Ba", "W"]
    increment = 3
    
    # Generate alloy compositions
    compositions = generate_alloy_compositions(elements, total_atoms, increment)

    for composition in compositions:
        alloy_name = f"alloy_{'_'.join(map(str, composition))}"
        filename = f"{alloy_name}.xsf"
        create_alloy(composition, elements, total_atoms, filename)

        # Check if the file was created successfully
        if not os.path.exists(filename):
            print(f"Error: File {filename} was not created.")
        else:
            print(f"Successfully created {filename}")

if __name__ == "__main__":
    generate_alloy_files()