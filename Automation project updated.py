# In order to generalize the molecules that are undergoing metal-catalysed reactions, we make a library

import os
import numpy as np
import pymatgen.core as mg
# In order to create a Slab from Bulk Fe, we need to perform the Slab operation, therefore,
from pymatgen.core.surface import SlabGenerator

VSEPR_VECTORS = {
    "linear_diatomic": np.array([
        [1.0, 0.0, 0.0]
        #[-1.0, 0.0, 0.0]
    ]),

    "linear_triatomic": np.array([
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0]
    ]),

    "trigonal_planar": np.array([
        [ 1.0, 0.0, 0.0],
        [-0.5, np.sqrt(3)/2, 0.0],
        [-0.5, -np.sqrt(3)/2, 0.0]
    ]),
    
    "tetrahedral": np.array([
        [0.0, 0.0, 1.0],
        [2*np.sqrt(2)/3, 0.0, -1/3],
        [-np.sqrt(2)/3,  np.sqrt(2/3), -1/3],
        [-np.sqrt(2)/3, -np.sqrt(2/3), -1/3]
    ]),

    "square_planar": np.array([
        [np.sqrt(2)/2, np.sqrt(2)/2, 0.0],
        [-np.sqrt(2)/2, np.sqrt(2)/2, 0.0],
        [np.sqrt(2)/2, -np.sqrt(2)/2, 0.0],
        [-np.sqrt(2)/2, -np.sqrt(2)/2, 0.0],
    ]),

    "t_shaped": np.array([
        [0.0, 1.0, 0.0],
        [np.cos(np.radians(86.2)), np.sin(np.radians(86.2)), 0.0],
        [-np.cos(np.radians(86.2)), np.sin(np.radians(86.2)), 0.0]
    ]),

    "trigonal_bipyramidal": np.array([
        [0.0, 0.0, 1.05],
        [0.0, 0.0, -1.05],
        [1.0, 0.0, 0.0],
        [np.cos(np.radians(120)), np.sin(np.radians(120)), 0.0],
        [np.cos(np.radians(240)), np.sin(np.radians(240)), 0.0]
    ]),

    "octahedral": np.array([
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0]
    ]),
}

def bent_vectors (bond_angle):
    half_bond_angle = np.radians(bond_angle/2)

    return np.array([
        [np.sin(half_bond_angle), np.cos(half_bond_angle), 0.0],
        [-np.sin(half_bond_angle), np.cos(half_bond_angle), 0.0],
    ])

def trigonal_pyramidal_vectors (lp_centre_ligand_angle):
    supplementary_angle = np.radians(180)-np.radians(lp_centre_ligand_angle)

    return np.array([
        [1.0 * np.sin(supplementary_angle), 0.0, np.cos(supplementary_angle)],
        [-0.5 * np.sin(supplementary_angle), np.sqrt(3)/2 * np.sin(supplementary_angle), np.cos(supplementary_angle)],
        [-0.5 * np.sin(supplementary_angle), -np.sqrt(3)/2 * np.sin(supplementary_angle), np.cos(supplementary_angle)]
    ])

def molecule_generated (centre, ligands, geometry, bond_length, bond_angle=None, lp_centre_ligand_angle=None):

    if geometry == "bent":
        molecular_geometry = bent_vectors(bond_angle)

    elif geometry == "trigonal_pyramidal":
        molecular_geometry = trigonal_pyramidal_vectors(lp_centre_ligand_angle)

    else:
        molecular_geometry = VSEPR_VECTORS[geometry]

    symbols = [centre]
    coordinates = [[0.0, 0.0, 0.0]]

    for element, vector, length in zip(ligands, molecular_geometry, bond_length):
        symbols.append(element)
        coordinates.append((vector * length).tolist())
        
    return symbols, np.array(coordinates)

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
    "CO": {"centre": "C", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.128], "floating_distance": 2.0},
    "NO": {"centre": "N", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.150], "floating_distance": 1.8},

    "H2": {"centre": "H", "ligands": ["H"], "geometry": "linear_diatomic", "bond_length": [0.741], "floating_distance": 1.2},
    "O2": {"centre": "O", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.208], "floating_distance": 1.9},
    "N2": {"centre": "N", "ligands": ["N"], "geometry": "linear_diatomic", "bond_length": [1.098], "floating_distance": 2.1},
    
    "CO2": {"centre": "C", "ligands": ["O", "O"], "geometry": "linear_triatomic", "bond_length": [1.162, 1.162], "floating_distance": 2.2},
    "N2O": {"centre": "N", "ligands": ["N", "O"], "geometry": "linear_triatomic", "bond_length": [1.120, 1.190], "floating_distance": 2.0},

    "NO2": {"centre": "N", "ligands": ["O", "O"], "geometry": "bent", "bond_length": [1.197, 1.197], "bond_angle": 134.1, "floating_distance": 1.9},
    "SO2": {"centre": "S", "ligands": ["O", "O"], "geometry": "bent", "bond_length": [1.430, 1.430], "bond_angle": 119.5, "floating_distance": 2.2},
    "H2O": {"centre": "O", "ligands": ["H", "H"], "geometry": "bent", "bond_length": [0.958, 0.958], "bond_angle": 104.5, "floating_distance": 2.2},
    
    "CH4": {"centre": "C", "ligands": ["H", "H", "H", "H"], "geometry": "tetrahedral", "bond_length": [1.087, 1.087, 1.087, 1.087], "floating_distance": 3.5},
    "NH3": {"centre": "N", "ligands": ["H", "H", "H"], "geometry": "trigonal_pyramidal", "bond_length": [1.012, 1.012, 1.012], "lp_centre_ligand_angle": 112.0, "floating_distance": 2.0}
}

# Creating the atomic_species block

involved_elements_informations = {
    "Fe": {"mass": 55.845, "pseudopotential": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"},
    "C":  {"mass": 12.011, "pseudopotential": "C.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "O":  {"mass": 15.999, "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "N":  {"mass": 14.007, "pseudopotential": "N.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "H":  {"mass": 1.008,  "pseudopotential": "H.pbe-rrkjus_psl.1.0.0.UPF"},
    "S":  {"mass": 32.059, "pseudopotential": "S.pbe-n-kjpaw_psl.1.0.0.UPF"}
}

# Converting to a "readable" format (readable by QEspresso)

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

for reacting_molecule, molecular_data in reacting_molecules.items():
    folder_path = os.path.join("calculations", reacting_molecule)
    os.makedirs(folder_path, exist_ok=True)

    symbols, molecule_coordinates_relative = molecule_generated(
        centre = molecular_data["centre"],
        ligands = molecular_data["ligands"],
        geometry = molecular_data["geometry"],
        bond_length = molecular_data["bond_length"],
        bond_angle = molecular_data.get("bond_angle"),
        lp_centre_ligand_angle = molecular_data.get("lp_centre_ligand_angle")
    )

    adsorption_site = rearranged_sites[0]

    anchor_position = np.array([adsorption_site.coords[0], adsorption_site.coords[1], adsorption_site.coords[2]+ molecular_data["floating_distance"]])

    individual_atomic_positions = (anchor_position + molecule_coordinates_relative)

    fe_surface_with_gas_molecule = fe_surface.copy()

    for element, coords in zip(symbols, individual_atomic_positions):
        fe_surface_with_gas_molecule.append(
            species=element,
            coords=coords,
            coords_are_cartesian=True
        )

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