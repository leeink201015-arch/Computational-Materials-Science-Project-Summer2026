# So far, we've modelled all the possible adsorption geometry and positions for the adsorbate on the Fe(100) surface. Now we need to evaluate the thermodynamic favorability of each adsorption reaction.
# To calculate the adsorption energy, we use the formula: E_ads = E_total - (E_surface + E_adsorbate), where E_total is the total energy of the surface with the adsorbate, E_surface is the energy of the clean surface, and E_adsorbate is the energy of the isolated adsorbate molecule.
# We already know E_total from the relaxation calculations
# So we need to calculate E_surface and E_adsorbate. For E_surface, we can use the clean Fe(100) surface, and for E_adsorbate, we can use the isolated gas molecule in a vacuum.

import os
import numpy as np
import pymatgen.core as mg

from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

from scipy.spatial.transform import Rotation as R

VSEPR_VECTORS = {
    "linear_diatomic": np.array([
        [1.0, 0.0, 0.0]
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

# New lines:
tilt_angles = np.array([0, 30, 60, 90, 120, 150, 180])
azimuth_angles = np.array([0, 60, 120, 180, 240, 300])

a = 2.8326
c = a*8.2958

lattice = mg.Lattice.tetragonal(2*a,c)

relaxed_layers = [
    (0.0000000000, 0.0000000000, 0.3189370858),
    (0.5000000000, 0.5000000000, 0.3778786719),
    (0.0000000000, 0.0000000000, 0.4398482736),
    (0.5000000000, 0.5000000000, 0.5000000000),
    (0.0000000000, 0.0000000000, 0.5601517264),
    (0.5000000000, 0.5000000000, 0.6221213281),
    (0.0000000000, 0.0000000000, 0.6810629142),
]

positions_2x2 = []

for x, y, z in relaxed_layers:
    for i in range(2):
        for j in range(2):
            positions_2x2.append([(x+i)/2, (y+j)/2, z])

fe_surface = mg.Structure(lattice, ["Fe"]*len(positions_2x2), positions_2x2)
print("# of Fe atoms:", len(fe_surface))
print("Lattice parameters:", fe_surface.lattice.abc)
assert len(fe_surface) == 28
# Ima run it but only to test the new slab geometry after modifications...
# raise SystemExit("Geometry check complete")

print(f"Number of atoms in surface slab: {len(fe_surface)}")

# Generate the structure and find the number of atoms

fe_surface.to(fmt="poscar", filename="fe_100_surface.vasp")

print(f"File saved in {os.getcwd()}")

# Just created a slab surface, now let AdsobateSiteFinder deal with the cases with different adsorption positions

adsorption_finder = AdsorbateSiteFinder(fe_surface)
adsorption_sites = adsorption_finder.find_adsorption_sites(positions=["ontop", "bridge", "hollow"], symm_reduce=False)

for site_type, adsorbing_positions in adsorption_sites.items():

    print(f"\n{site_type}")

    for position in adsorbing_positions:
        print(position)

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

involved_elements_informations = {
    "Fe": {"mass": 55.845, "pseudopotential": "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"},
    "C":  {"mass": 12.011, "pseudopotential": "C.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "O":  {"mass": 15.999, "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "N":  {"mass": 14.007, "pseudopotential": "N.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "H":  {"mass": 1.008,  "pseudopotential": "H.pbe-rrkjus_psl.1.0.0.UPF"},
    "S":  {"mass": 32.059, "pseudopotential": "S.pbe-n-kjpaw_psl.1.0.0.UPF"}
}

def write_qespresso_input (structure, molecular_name, output_file, calculation = "relax"):
    element_present = structure.composition.elements
    num_atoms = len(structure)
    num_species = len(element_present)

    with open(output_file, "w") as file:
        # Control Block
        file.write("&CONTROL\n")
        file.write(f"calculation = '{calculation}'\n")
        file.write("restart_mode = 'from_scratch'\n")
        file.write(f"prefix = 'fe_{molecular_name}'\n")
        file.write("pseudo_dir = '/usr/share/espresso/pseudo/'\n") # Made a change here, because pseudopotentials are actually installed in /usr/share/espresso/pseudo/
        file.write("outdir = './tmp/'\n")
        file.write("/\n\n")

        # System Block
        file.write("&SYSTEM\n")
        file.write("ibrav = 0\n")
        file.write(f"nat = {num_atoms}\n")
        file.write(f"ntyp = {num_species}\n")
        file.write("nosym = .true.\n")  # Added this line to disable symmetry operations, which is important for surface calculations
        file.write("ecutwfc = 65.0\n")
        file.write("ecutrho = 782.0\n")
        file.write("nspin = 2\n")  # Change, added spin-polarized calculation for Fe
        file.write("starting_magnetization(1) = 0.5\n")  # Initial magnetic moment for Fe. Added these two lines I didn't have them before lol
        file.write("occupations = 'smearing'\n")
        file.write("smearing = 'mv'\n")
        file.write("degauss = 0.02\n")
        file.write("/\n\n")

        # Electrons Block
        file.write("&ELECTRONS\n")
        file.write("conv_thr = 1.0d-8\n")
        file.write("mixing_beta = 0.3\n")  # Changed from 0.7 to 0.3 so that the updates for electron density are more conservation. Important because the complex Fe+CO system can make convergence challenging
        file.write("/\n\n")

        # Ions Block (for relaxation, so block is skippable if calculation is not "relax")
        if calculation == "relax":
            file.write("&IONS\n")
            file.write("ion_dynamics = 'bfgs'\n")
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

        # Atomic Positions (modified, added mobility to the top layer Fe atoms)
        file.write("ATOMIC_POSITIONS angstrom\n")

        fe_sites = [
            site for site in structure.sites
            if site.species_string == "Fe"
        ]

        if fe_sites:
            top_fe_z_coord = max(
                [s.coords[2] for s in fe_sites]
            )

        tolerance = 0.1  # Define a small tolerance for floating-point comparison

        for site in structure.sites:
            if site.species_string != "Fe":
                flags = "1 1 1"  # Basically the adsorbates

            elif abs(site.coords[2] - top_fe_z_coord) < tolerance:
                flags = "1 1 1"  # Top layer Fe atoms (giving mobility)

            else:
                flags = "0 0 0"  # Lower layers of Fe atoms (fixed, no relaxation)

            file.write(f"{site.species_string:2s}  {site.coords[0]:12.6f} {site.coords[1]:12.6f} {site.coords[2]:12.6f} {flags}\n")
        file.write("\n")

        # K-Points
        file.write("K_POINTS automatic\n")
        file.write("5 5 1 0 0 0\n")  # Changed from 3x3x1 to 5x5x1 because we used 5x5x1 for the clean Fe(100) surface calculations so it's important to keep that consistency. 5x5x1 and 4x4x1 have differen Brillouin-zone sampling errors

    print(f"Quantum ESPRESSO input file: {output_file}")

# Test one CO moleucle on an ontop site
co_data = reacting_molecules["CO"]

symbols, molecule_coordinates = molecule_generated(
    centre=co_data["centre"],
    ligands=co_data["ligands"],
    geometry=co_data["geometry"],
    bond_length=[1.1411]
)

ontop_site = adsorption_sites["ontop"][0].copy()

top_fe_z = max(site.coords[2] for site in fe_surface)
ontop_site[2] = top_fe_z + 2.0

co_coordinates = molecule_coordinates @ R.from_euler(
    "y", -90, degrees=True
).as_matrix().T

fe_co_system = fe_surface.copy()

for element, coords in zip(symbols, co_coordinates):
    fe_co_system.append(element, coords + ontop_site, coords_are_cartesian=True)

print("Fe + CO atoms: ", len(fe_co_system))
print("Carbon: ", fe_co_system[-2].coords)
print("Oxygen: ", fe_co_system[-1].coords)

fe_co_system.to(fmt="cif", filename="fe100_2x2_co_ontop.cif")

write_qespresso_input(
    structure=fe_co_system,
    molecular_name="CO_ontop",
    output_file="fe100_2x2_co_ontop.in",
    calculation="relax"
)

raise SystemExit("Fe + CO geometry and QE input generated.")

print("Top Fe z:", max(site.coords[2] for site in fe_surface))
print("C-Fe vertical distance:", fe_co_system[-2].coords[2] - max(site.coords[2] for site in fe_surface))
print("C-O distance:", fe_co_system.get_distance(28, 29))

ontop_sites = adsorption_finder.find_adsorption_sites(
    positions = ["ontop"]
)

bridge_sites = adsorption_finder.find_adsorption_sites(
    positions = ["bridge"]
)

hollow_sites = adsorption_finder.find_adsorption_sites(
    positions = ["hollow"]
)

# Trying to find the uppermost iron atom that's going to react with CO first, then we arrange them in order

# print(fe_surface) First getting a better visual of our surface atoms
# rearranged_sites = sorted(fe_surface, key = lambda site: site.coords[2], reverse = True)

# print(rearranged_sites)

# Converting to a "readable" format (readable by QEspresso)

# Deleted the first write_qespresso_input function because it was redundant and not used in the code. The second write_qespresso_input function is the one that will be used to generate the input files for QE calculations.

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

# So how to create files for all three kinds of adsorption positions...
# Trying to use a nested for loop idk:

for reacting_molecule, molecular_data in reacting_molecules.items():

    symbols, molecule_coordinates_relative = molecule_generated(
        centre = molecular_data["centre"],
        ligands = molecular_data["ligands"],
        geometry = molecular_data["geometry"],
        bond_length = molecular_data["bond_length"],
        bond_angle = molecular_data.get("bond_angle"),
        lp_centre_ligand_angle = molecular_data.get("lp_centre_ligand_angle")
    )

    for site_type in ["ontop", "bridge", "hollow"]:
        
        for site_number, adsorption_site in enumerate (adsorption_sites[site_type]):

            for tilt_angle in tilt_angles:
                R_tilted = np.array([
                    [ np.cos(np.radians(tilt_angle)), 0, np.sin(np.radians(tilt_angle))],
                    [0, 1, 0],
                    [-np.sin(np.radians(tilt_angle)), 0, np.cos(np.radians(tilt_angle))]
                ])

                for azimuth_angle in azimuth_angles:
                    R_azimuth = np.array([
                        [np.cos(np.radians(azimuth_angle)), -np.sin(np.radians(azimuth_angle)), 0],
                        [np.sin(np.radians(azimuth_angle)), np.cos(np.radians(azimuth_angle)), 0],
                        [0, 0, 1]
                    ])

                    folder_path = os.path.join("calculations", reacting_molecule, f"{site_type}_{site_number}", f"tilt_{tilt_angle}_azimuth_{azimuth_angle}")
                    os.makedirs(folder_path, exist_ok=True)

                    # Matrix multiplication (it's not communtative), Final Point P_f = R_azimuth x (R_tilt x P_i) (I'll convert this to LaTeX later whatever)

                    molecule_coordinates_tilted = molecule_coordinates_relative @ R_tilted.T
                    molecule_coordinates_rotated = molecule_coordinates_tilted @ R_azimuth.T

                    anchor_position = adsorption_site.copy()
                    anchor_position[2] += molecular_data["floating_distance"]
                    individual_atomic_positions = (anchor_position + molecule_coordinates_rotated)

                    fe_surface_with_gas_molecule = fe_surface.copy()

                    for element, coords in zip(symbols, individual_atomic_positions):
                        fe_surface_with_gas_molecule.append(
                            species = element,
                            coords = coords,
                            coords_are_cartesian = True
                        )

                    qe_path = os.path.join(folder_path, "qespresso_input.in")
                    cif_path = os.path.join(folder_path, f"{reacting_molecule}_{site_type}_{site_number}_tilt{tilt_angle}_az{azimuth_angle}.cif")
                    slurm_path = os.path.join(folder_path, "run_job.sh")

                    fe_surface_with_gas_molecule.to(fmt = "cif", filename = cif_path)

                    write_qespresso_input(
                        structure = fe_surface_with_gas_molecule,
                        molecular_name = reacting_molecule,
                        output_file = qe_path
                    )

                    write_slurm_script(
                        molecular_name = reacting_molecule,
                        output_file = slurm_path
                    )

                    print(f"File saved in {os.getcwd()}")

# Deleted that entire section at the bottom becasue they were written to generate a clean Fe(100) realaxation input
# I'm not doing that rn, I don't need it