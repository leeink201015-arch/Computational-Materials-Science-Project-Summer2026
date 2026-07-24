# In order to generalize the molecules that are undergoing metal-catalysed reactions, we make a library

import os
import pymatgen.core as mg
# In order to create a Slab from Bulk Fe, we need to perform the Slab operation, therefore,
from pymatgen.core.surface import SlabGenerator

bcc_lattice = mg.Lattice.cubic(2.87) # Uppercase Lattice, because it's an object?

species_pure = ["Fe", "Fe"]

bcc_coordinates =([
    [0.00, 0.00, 0.00], 
    [0.50, 0.50, 0.50]
])

bulk_iron = mg.Structure(bcc_lattice, species_pure, bcc_coordinates)

# Generating a slab from bulky iron so that we can create a surface for catalysing reactions

slab_gen = SlabGenerator(
    initial_structure = bulk_iron,
    miller_index      = (1, 0, 0),
    min_slab_size     = 10.0,
    min_vacuum_size   = 10.0
)

slabs = slab_gen.get_slabs()
fe_surface = slabs[0]
print(f"Number of atoms in surface slab: {len(fe_surface)}")

# Generate the structure and find the number of atoms

fe_surface.to(fmt="poscar", filename="fe_100_surface.vasp")

print(f"File saved in {os.getcwd()}")

# Trying to find the uppermost iron atom that's going to react with CO first, then we arrange them in order

print(fe_surface) # First getting a better visual of our surface atoms
rearranged_sites = sorted(fe_surface, key = lambda site: site.coords[2], reverse = True)

print(rearranged_sites)

# Reacting gas molecules generalization

reacting_molecules = {
    "CO":   {"bottom": "C", "top": "O", "bond_length": 1.13, "floating_distance": 2.0},
    "CO2":  {"bottom": "C", "top": "O", "bond_length": 1.16, "floating_distance": 2.2}, # Note: CO2 is linear, treated as C-O vector here
    
    "NO":   {"bottom": "N", "top": "O", "bond_length": 1.15, "floating_distance": 1.8},
    "NO2":  {"bottom": "N", "top": "O", "bond_length": 1.20, "floating_distance": 1.9},
    "N2":   {"bottom": "N", "top": "N", "bond_length": 1.10, "floating_distance": 2.1},
    "N2O":  {"bottom": "N", "top": "O", "bond_length": 1.19, "floating_distance": 2.0},
    
    "CH":   {"bottom": "C", "top": "H", "bond_length": 1.09, "floating_distance": 1.5},
    "OH":   {"bottom": "O", "top": "H", "bond_length": 0.96, "floating_distance": 1.4},
    "NH":   {"bottom": "N", "top": "H", "bond_length": 1.04, "floating_distance": 1.5},
    "HF":   {"bottom": "F", "top": "H", "bond_length": 0.92, "floating_distance": 1.3},
    
    "O2":   {"bottom": "O", "top": "O", "bond_length": 1.21, "floating_distance": 1.9},
    "H2":   {"bottom": "H", "top": "H", "bond_length": 0.74, "floating_distance": 1.2},
    "Cl2":  {"bottom": "Cl", "top": "Cl", "bond_length": 1.99, "floating_distance": 2.3}
}

# Creating the atomic_species block

involved_elements_informations = {
    "Fe": {"mass": 55.845, "pseudopotential": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"},
    "C":  {"mass": 12.011, "pseudopotential": "C.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "O":  {"mass": 15.999, "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "N":  {"mass": 14.007, "pseudopotential": "N.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "H":  {"mass": 1.008,  "pseudopotential": "H.pbe-rrkjus_psl.1.0.0.UPF"},
    "F":  {"mass": 18.998, "pseudopotential": "F.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "Cl": {"mass": 35.453, "pseudopotential": "Cl.pbe-n-kjpaw_psl.1.0.0.UPF"}
}

# Automation (Quantum Espresso) converting to a readable format

def write_qespresso_input (structure, molecular_name, output_file):
    element_present = structure.composition.elements
    num_atoms = len(structure)
    num_species = len(element_present)

    with open(output_file, "w") as file:
        # Control Block
        file.write("&CONTROL\n")
        file.write("calculation = 'scf'\n")
        file.write("restart_mode = 'from_scratch'\n")
        file.write(f"prefix = 'fe_{molecular_name}'\n")
        file.write("pseudo_dir = './'\n")
        file.write("outdir = './tmp/'\n")
        file.write("/\n\n")

        # System Block
        file.write("&SYSTEM\n")
        file.write("ibrav = 0\n")
        file.write(f"nat = {num_atoms}\n")
        file.write(f"ntyp = {num_species}\n")
        file.write("ecutwfc = 45.0\n")
        file.write("ecutrho = 360.0\n")
        file.write("occupations = 'smearing'\n")
        file.write("smearing = 'mv'\n")
        file.write("degauss = 0.02\n")
        file.write("/\n\n")

        # Electrons Block
        file.write("&ELECTRONS\n")
        file.write("conv_thr = 1.0d-8\n")
        file.write("mixing_beta = 0.7\n")
        file.write("/\n\n")

        # Atomic Species
        file.write("ATOMIC_SPECIES\n")
        for elem in element_present:
            sym = elem.symbol
            data = involved_elements_informations[sym]
            file.write(f"{sym:2s}  {data['mass']:7.3f}  {data['pseudopotential']}\n")
        file.write("\n")

        # Cell Parameters
        file.write("CELL_PARAMETERS angstrom\n")
        for vec in structure.lattice.matrix:
            file.write(f"{vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
        file.write("\n")

        # Atomic Positions
        file.write("ATOMIC_POSITIONS angstrom\n")
        for site in structure.sites:
            file.write(f"{site.species_string:2s}  {site.coords[0]:12.6f} {site.coords[1]:12.6f} {site.coords[2]:12.6f}\n")
        file.write("\n")

        # K-Points
        file.write("K_POINTS automatic\n")
        file.write("3 3 1 0 0 0\n")

    print(f"Quantum ESPRESSO input file: {output_file}")

def write_slurm_script(molecular_name, output_file):
    with open(output_file, "w") as file:
        file.write("#!/bin/bash\n")
        file.write(f"#SBATCH --job-name=fe_{molecular_name}\n")
        file.write("#SBATCH --output=qe_job.out\n")
        file.write("#SBATCH --error=qe_job.err\n")
        file.write("#SBATCH --nodes=1\n")
        file.write("#SBATCH --ntasks-per-node=16\n")
        file.write("#SBATCH --time=04:00:00\n")
        file.write("#SBATCH --partition=standard\n\n")

        file.write("module load quantum-espresso\n\n")
        file.write("mpirun -np 16 pw.x -in qespresso_input.in > qe_output.out\n")

# Following last project, instead of making it interactive by having the input function, we create a for loop so that everything's being saved in file

for reacting_molecule in reacting_molecules:
    folder_path = os.path.join("calculations", reacting_molecule)
    os.makedirs(folder_path, exist_ok=True)

    molecular_data = reacting_molecules[reacting_molecule]

    bottom_atom_coordinates = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + reacting_molecules[reacting_molecule]["floating_distance"]]
    top_atom_coordinates = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + reacting_molecules[reacting_molecule]["floating_distance"] + reacting_molecules[reacting_molecule]["bond_length"]]

    fe_surface_with_gas_molecule = fe_surface.copy()

    fe_surface_with_gas_molecule.append(species=reacting_molecules[reacting_molecule]["bottom"], coords=bottom_atom_coordinates, coords_are_cartesian=True)
    fe_surface_with_gas_molecule.append(species=reacting_molecules[reacting_molecule]["top"], coords=top_atom_coordinates, coords_are_cartesian=True)

    vasp_path = os.path.join(folder_path, "fe_surface_with_gas.vasp")
    qe_path = os.path.join(folder_path, "qespresso_input.in")
    slurm_path = os.path.join(folder_path, "run_job.sh")

    fe_surface_with_gas_molecule.to(fmt="poscar", filename=vasp_path)

    write_qespresso_input(
        structure=fe_surface_with_gas_molecule,
        molecular_name=reacting_molecule,
        output_file=qe_path
    )

    write_slurm_script(
        molecular_name=reacting_molecule,
        output_file=slurm_path
    )

    print(f"File saved in {os.getcwd()}")