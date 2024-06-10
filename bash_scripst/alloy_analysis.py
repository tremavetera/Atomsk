import os
import subprocess
from ase.io import read, write
from ase.calculators.vasp import Vasp
from py4vasp import Calculation

# Function to generate alloy compositions using bash script
def generate_alloy_compositions():
    subprocess.run(["./generate_alloys.sh"], check=True)

# Function to create atomic structures using ASE
def create_structure(alloy_name):
    structure = read(f"{alloy_name}.xsf")
    return structure

# Function to run VASP calculation
def run_vasp_calculation(structure, alloy_name):
    calc = Vasp(
        directory=f'vasp_calculations/{alloy_name}',
        xc='PBE',
        kpts=[4, 4, 4],
        encut=350,
        ibrion=2,
        nsw=50,
        atoms=structur
    )
    structure.set_calculator(calc)
    calc.calculate()

# Function to analyze VASP output using Py4VASP
def analyze_vasp_output(alloy_name):
    calc = Calculation.from_path(f'vasp_calculations/{alloy_name}')
    energy = calc.energy
    print(f"Energy for {alloy_name}: {energy}")
    return energy

# Main workflow
def main():
    # Generate alloy compositions using the bash script
    generate_alloy_compositions()

    # List all generated alloy files
    alloy_files = [f for f in os.listdir() if f.startswith('alloy_') and f.endswith('.xsf')]

    results = []

    for alloy_file in alloy_files:
        alloy_name = alloy_file.split('.')[0]
        structure = create_structure(alloy_name)
        run_vasp_calculation(structure, alloy_name)
        energy = analyze_vasp_output(alloy_name)
        results.append((alloy_name, energy))

    # Save results
    with open('alloy_analysis_results.txt', 'w') as f:
        for alloy, energy in results:
            f.write(f"{alloy}: {energy}\n")

    print("Alloy analysis completed.")

if __name__ == "__main__":
    main()