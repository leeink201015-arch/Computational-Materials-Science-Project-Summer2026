# In order to generalize the molecules that are undergoing metal-catalysed reactions, we make a library

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

fe_surface.to(fmt="cif", filename="fe_100_surface.cif")
import os; os.system("explorer .")
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

reacting_molecule_input = input("Which gaseous molecule should we simulate: ")
if reacting_molecule_input in reacting_molecules:
    print("Molecule found, generating reaction model")

    molecular_data = reacting_molecules[reacting_molecule_input]

    bottom_atom_coordinates = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + reacting_molecules[reacting_molecule_input]["floating_distance"]]
    top_atom_coordinates = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + reacting_molecules[reacting_molecule_input]["floating_distance"] + reacting_molecules[reacting_molecule_input]["bond_length"]]

    fe_surface_with_gas_molecule = fe_surface.copy()

    fe_surface_with_gas_molecule.append(species=reacting_molecules[reacting_molecule_input]["bottom"], coords=bottom_atom_coordinates, coords_are_cartesian=True)
    fe_surface_with_gas_molecule.append(species=reacting_molecules[reacting_molecule_input]["top"], coords=top_atom_coordinates, coords_are_cartesian=True)

    fe_surface_with_gas_molecule.to(fmt="cif", filename="fe_surface_with_gas_molecule.cif")
    import os; os.system("explorer .")
    print(f"File saved in {os.getcwd()}")

else:
    print(f"The molecule {reacting_molecule_input} is not in included.")