# TRYING TO RUN CONVERGENCE TESTS BUT SO FAR I'M ONLY PARTIALLY FINISHED, i'VE ALREADY CONVERGED THE CUTOFF ENERGIES, BUT NOT THE K-POINTS YET, WHICH I WILL DO NEXT IN THE FILE "K POINTS CONVERGENCE TEST.PY"

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

# Just created a slab surface, now let AdsobateSiteFinder deal with the cases with different adsorption positions

adsorption_finder = AdsorbateSiteFinder(fe_surface)
adsorption_sites = adsorption_finder.find_adsorption_sites()

# Manually defining the hollow site, because AdsorbateSiteFinder doesn't find it for some reason (maybe because it's a 4-fold hollow site)
# It's hollow (x,y) = top Fe (x,y) + (1/2, 1/2), so yea:

top_fe_site = max(
    (site for site in fe_surface if site.species_string == "Fe"),
    key=lambda site: site.coords[2]
)

top_fe_fractional = top_fe_site.frac_coords.copy()

hollow_fractional = top_fe_fractional.copy()
hollow_fractional[0] = (hollow_fractional[0] + 0.5) % 1.0
hollow_fractional[1] = (hollow_fractional[1] + 0.5) % 1.0

hollow_position = fe_surface.lattice.get_cartesian_coords(
    hollow_fractional
)

adsorption_sites["hollow"] = [hollow_position]

for site_type, adsorbing_positions in adsorption_sites.items():

    print(f"\n{site_type}")

    for position in adsorbing_positions:
        print(position)

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

# Reacting gas molecules generalization
# Updates from Aug 30, 2026: Added the anchor_index key to the dictionary, doing this allow us to specify which atom in the molecule is going to be anchored to the surface, shifting the molecule around, instead of fixing the anchor to [0, 0, 0] (the centre atom of the molecule).

reacting_molecules = {
    "CO": {"centre": "C", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.128], "floating_distance": 2.0, "anchor_index": 0},
    "NO": {"centre": "N", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.150], "floating_distance": 1.8, "anchor_index": 0},

    "H2": {"centre": "H", "ligands": ["H"], "geometry": "linear_diatomic", "bond_length": [0.741], "floating_distance": 1.2, "anchor_index": 0},
    "O2": {"centre": "O", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.208], "floating_distance": 1.9, "anchor_index": 0},
    "N2": {"centre": "N", "ligands": ["N"], "geometry": "linear_diatomic", "bond_length": [1.098], "floating_distance": 2.1, "anchor_index": 0},
    
    "CO2": {"centre": "C", "ligands": ["O", "O"], "geometry": "linear_triatomic", "bond_length": [1.162, 1.162], "floating_distance": 2.2, "anchor_index": 0},
    "N2O": {"centre": "N", "ligands": ["N", "O"], "geometry": "linear_triatomic", "bond_length": [1.120, 1.190], "floating_distance": 2.0, "anchor_index": 0},

    "NO2": {"centre": "N", "ligands": ["O", "O"], "geometry": "bent", "bond_length": [1.197, 1.197], "bond_angle": 134.1, "floating_distance": 1.9, "anchor_index": 0},
    "SO2": {"centre": "S", "ligands": ["O", "O"], "geometry": "bent", "bond_length": [1.430, 1.430], "bond_angle": 119.5, "floating_distance": 2.2, "anchor_index": 0},
    "H2O": {"centre": "O", "ligands": ["H", "H"], "geometry": "bent", "bond_length": [0.958, 0.958], "bond_angle": 104.5, "floating_distance": 2.2, "anchor_index": 0},
    
    "CH4": {"centre": "C", "ligands": ["H", "H", "H", "H"], "geometry": "tetrahedral", "bond_length": [1.087, 1.087, 1.087, 1.087], "floating_distance": 3.5, "anchor_index": 0},
    "NH3": {"centre": "N", "ligands": ["H", "H", "H"], "geometry": "trigonal_pyramidal", "bond_length": [1.012, 1.012, 1.012], "lp_centre_ligand_angle": 112.0, "floating_distance": 2.0, "anchor_index": 0}
}

# Creating the atomic_species block

involved_elements_informations = {
    "Fe": {"mass": 55.845, "pseudopotential": "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"},
    "C":  {"mass": 12.011, "pseudopotential": "C.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "O":  {"mass": 15.999, "pseudopotential": "O.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "N":  {"mass": 14.007, "pseudopotential": "N.pbe-n-kjpaw_psl.1.0.0.UPF"},
    "H":  {"mass": 1.008,  "pseudopotential": "H.pbe-rrkjus_psl.1.0.0.UPF"},
    "S":  {"mass": 32.059, "pseudopotential": "S.pbe-n-kjpaw_psl.1.0.0.UPF"}
}

# I got these results from QEspresso after running a relaxation calculation yesterday (Aug 29, 2026), now I'm using these results to create a new structure with the relaxed coordinates, so that it can be used for convergence tests

final_relaxed_coords = np.array([
    [0.0000000000,  0.0000000000,  2.1525000000],
    [1.4350000000,  1.4350000000,  0.7175000000],
    [0.0000000000,  0.0000000000,  5.0225000000],
    [1.4350000000,  1.4350000000,  3.5875000000],
    [-0.0000000223, -0.0000000316, 7.8467151517],
    [1.4349999930,  1.4349999972,  6.4171424594],
    [-0.0000000481, 0.0000000207, 10.6837505334],
    [1.4349999766,  1.4349999800,  9.3029766345],
])

# Letting Pymatgen create a new structure with the relaxed coordinates (used for convergence tests)

relaxed_fe_surface = mg.Structure(
    fe_surface.lattice,  # Using the same lattice as the original Fe surface
    [site.species_string for site in fe_surface],
    final_relaxed_coords,
    coords_are_cartesian=True
)

# Converting to a "readable" format (readable by QEspresso)

# (Update from Aug 30, 2026): I basically cleaned up the QE input file generation function, because I had two write_qespresso_input functions, one for the clean Fe surface and one for the Fe surface with adsorbates. Now I combined them into one function.
# This is so far the biggest duplication bug risk ever encountered in this project

# Wondered if this would cause inaccuracy in QE calculations (but we'll see later on)

# Function header modified (Aug 30, 2026): while doing actual convergence tests, the ecutwfc, ecutrho, and k_grid values will be changed, so I added them as function parameters with default values. This way, we can easily change them when doing convergence tests without having to modify the function itself.
def write_qespresso_input (structure, molecular_name, output_file, calculation = "relax", ecutwfc=45.0, ecutrho=360.0, k_grid=(4, 4, 1)):
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
        file.write(f"ecutwfc = {ecutwfc}\n")  # Added this line to set the plane-wave cutoff energy for wavefunctions (in Ry) (August 30, 2026)
        file.write(f"ecutrho = {ecutrho}\n")  # Added this line to set the plane-wave cutoff energy for charge density and potential (in Ry) (August 30, 2026)
        file.write("occupations = 'smearing'\n")
        file.write("smearing = 'mv'\n")
        file.write("nspin = 2\n")  # Added this line to enable spin-polarized calculations for ferromagnetic Fe
        file.write("starting_magnetization(1) = 0.5\n")  # Initial magnetization for Fe (can be adjusted based on the system)
        file.write("degauss = 0.02\n")
        file.write("/\n\n")

        # Electrons Block
        file.write("&ELECTRONS\n")
        file.write("conv_thr = 1.0d-8\n")
        file.write("mixing_beta = 0.7\n")
        file.write("/\n\n")

        # Ions Block (if statement no longer true)
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
        file.write("CELL_PARAMETERS (angstrom)\n")
        for vec in structure.lattice.matrix:
            file.write(f"{vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")
        file.write("\n")

        # Atomic Positions (modified, added mobility to the top layer Fe atoms)
        # Update: From previous QEspresso calculations, we found that the top layer contraction was way larger than expected (0.39332 Å inwards, 27%)
        # This shows bad accuracy I think, because I only relaxed the top layer of Fe atoms, which was not enough to get a good relaxation result.
        # So now we are going to relax the top 4 layers of Fe atoms, and fix the lower 4 layers of Fe atoms. This should give a more accurate relaxation result.
        
        file.write("ATOMIC_POSITIONS (angstrom)\n")

        fe_sites = [
            site for site in structure.sites
            if site.species_string == "Fe"
        ]

        # Moved out of the for loop cuz it's not necessary to recalculate for every site

        fe_z_coords = sorted(
            set(round(site.coords[2], 6) for site in fe_sites)
        )

        relaxed_fe_layers = fe_z_coords[-4:]  # Top 4 layers of Fe atoms will be relaxed

        tolerance = 0.1  # Define a small tolerance for floating-point comparison

        for site in structure.sites:

            if site.species_string != "Fe":
                flags = "1 1 1"  # Basically the adsorbates

            elif any(
                abs(site.coords[2] - z_layer) < tolerance
                for z_layer in relaxed_fe_layers
            ):
                flags = "1 1 1"  # Relaxing all Fe atoms in the top 4 layers

            else:
                flags = "0 0 0"  # The bottom 4 layers of Fe atoms (fixed, no relaxation)

            if calculation == "relax":
                file.write(f"{site.species_string:2s}  {site.coords[0]:12.6f} {site.coords[1]:12.6f} {site.coords[2]:12.6f} {flags}\n")  # Flags are only necessary if relaxed
            else:
                file.write(f"{site.species_string:2s}  {site.coords[0]:12.6f} {site.coords[1]:12.6f} {site.coords[2]:12.6f}\n")  # Flags aren't necessary if not relaxed, in which all atoms are fixed
        file.write("\n")

        # K-Points
        file.write("K_POINTS automatic\n")
        file.write(f"{k_grid[0]} {k_grid[1]} {k_grid[2]} 0 0 0\n")  # Changed from fixed 4x4x1 to a variable k_grid, so that we can easily change it for convergence tests (Aug 30, 2026)

    print(f"Quantum ESPRESSO input file: {output_file}")

cutoff_tests = [
    (45.0, 360.0), 
    (55.0, 440.0),
    (65.0, 520.0),
]

# Adding for loop to generate input files for convergence tests with different cutoff energies

for ecutwfc, ecutrho in cutoff_tests:
    output_folder = os.path.join("calculations", "convergence_tests", "cutoff",f"ecut_{int(ecutwfc)}_ecutrho_{int(ecutrho)}")
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, "scf.in")

    write_qespresso_input(
        structure = relaxed_fe_surface,
        molecular_name = "clean_fe_surface",
        output_file = output_file,
        calculation = "scf",
        ecutwfc = ecutwfc,
        ecutrho = ecutrho,
        k_grid = (4, 4, 1)
    )

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

# Took this out of the for loop because we only need to calculate the top Fe atom once, not for every adsorption site
top_fe_z = max(
    site.coords[2]
    for site in fe_surface
    if site.species_string == "Fe"
)

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

                    anchor_index = molecular_data["anchor_index"]

                    molecule_coordinates_anchored = (molecule_coordinates_relative - molecule_coordinates_relative[anchor_index])
                    
                    molecule_coordinates_tilted = molecule_coordinates_anchored @ R_tilted.T
                    molecule_coordinates_rotated = molecule_coordinates_tilted @ R_azimuth.T

                    # Previously had                anchor_position = adsorption_site.copy()
                    #                               anchor_position[2] += molecular_data["floating_distance"]
                    # Made a change here because AbsorbateSiteFinder already returns the adsorption site coordinates positioned above the surface by its own default floating distance
                    # Shouldn't be adding another moleule-specific height on top of that
                    # Otherwise floating distance would be higher than expected

                    anchor_position = adsorption_site.copy()
                    anchor_position[2] = top_fe_z + molecular_data["floating_distance"]
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
                    # slurm_path = os.path.join(folder_path, "run_job.sh")
                    # Commented out the slurm_path because don't need to generate a slurm script for every single adsorption configuration, just one is enough (prob why GitHub thinks my repository is mostly Shell instead of Python)

                    fe_surface_with_gas_molecule.to(fmt = "cif", filename = cif_path)

                    write_qespresso_input(
                        structure = fe_surface_with_gas_molecule,
                        molecular_name = reacting_molecule,
                        output_file = qe_path
                    )

                    # write_slurm_script(
                    #     molecular_name = reacting_molecule,
                    #     output_file = slurm_path
                    # )
                    # Commented this out reason same as above

                    print(f"File saved in {os.getcwd()}")

# Now we are trying to generate a relaxation input file for the clean Fe(100) surface
clean_surface_folder = os.path.join("calculations", "references", "clean_fe_surface")
os.makedirs(clean_surface_folder, exist_ok=True)

clean_surface_relaxation_input_path = os.path.join(clean_surface_folder, "relax.in")

print("Generating clean surface relax.in now...")
print("Target path:", os.path.abspath(clean_surface_relaxation_input_path))

write_qespresso_input(
    structure=fe_surface,
    molecular_name="clean_fe_surface",
    output_file=clean_surface_relaxation_input_path,
    calculation="relax"
)

print("Finished write_qespresso_input")
print("File exists:", os.path.exists(clean_surface_relaxation_input_path))

print(f"File saved in {os.getcwd()}")

# Verifying the relaxed Fe surface coordinates

print("\nRelaxed Fe surface coordinates:")
for site in relaxed_fe_surface:
    print(site.species_string, site.coords)
