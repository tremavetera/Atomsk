import os
import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.io import write, read
from ase.calculators.vasp import Vasp
from py4vasp import Calculation
import make_alloys



# Function to run VASP calculation
def run_vasp_calculation(structure, alloy_name):
    calc = Vasp(
        directory=f'vasp_calculations/{alloy_name}',
        xc='PBE',
        kpts=[4, 4, 4],
        encut=350,
        ibrion=6,
        nsw=50,
        isif=3,
        atoms=structure
    )
# Trigger the calculation by requesting the potential energy
    try:
        energy = structure.get_potential_energy()
        return energy
    except Exception as e:
        print(f"Error during VASP calculation: {e}")
        return None

# Function to analyze VASP output using Py4VASP
def analyze_vasp_output(alloy_name):
    calc = Calculation.from_path(f'vasp_calculations/{alloy_name}')
    energy = calc.get_potential_energy()
    print(f"Energy for {alloy_name}: {energy}")
    return energy

# Isolated part to generate alloy files and analyze them
def generate_alloy_files_and_analyze():
    total_atoms = 54
    elements = ["V", "Cr", "Ti", "W"]
    increment = 3
    
    # Generate alloy compositions
    compositions = make_alloys.generate_alloy_compositions(elements, total_atoms, increment)

    results = []

    for composition in compositions:
        alloy_name = f"alloy_{'_'.join(map(str, composition))}"
        filename = f"{alloy_name}.xsf"
        make_alloys.create_alloy(composition, elements, total_atoms, filename,iscubic=False)

        # Check if the file was created successfully
        if not os.path.exists(filename):
            print(f"Error: File {filename} was not created.")
            continue

        structure = read(filename)
        print("this is the structure type in the main",type(structure))
        
        # Run VASP calculation and analyze output
        run_vasp_calculation(structure, alloy_name)
        energy = analyze_vasp_output(alloy_name)
        results.append((alloy_name, energy))

    # Save results
    with open('alloy_analysis_results.txt', 'w') as f:
        for alloy, energy in results:
            f.write(f"{alloy}: {energy}\n")

    print("Alloy analysis completed.")

if __name__ == "__main__":
    generate_alloy_files_and_analyze()
