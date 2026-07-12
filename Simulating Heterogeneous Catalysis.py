# Take pure BCC Iron (Ferrite) crystal, cut a (100) surface (creating a 2D slab)
# Visualize it to inspect the exposed reactive surface sites
# Then coordinates-wise "adsorb" a CO molecule to a dangling surface iron atom.

import pymatgen.core as mg # Pymatgen.core is where all the lattices and stuff are being stored
# In order to create a Slab from Bulk Fe, we need to perform the Slab operation, therefore,
from pymatgen.core.surface import SlabGenerator
# from pymatgen.visualization.plotter import StructurePlotter | old, unupdated, can't use it anymore
import matplotlib.pyplot as plt

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
# rearranged_sites[0] are the coordinates for the first Fe atom to react with CO
# print(rearranged_sites[0].coords[2]) is the z coordinate for the first Fe atom to react with CO

# Have our 8 Fe surface atoms located, now we add the CO gas to start a rxn

# It should absorb like this:
#       ..
#       O
#      |||
#       C
#       ..
# Fe Fe Fe Fe Fe
# Fe Fe Fe Fe Fe

co_floating_distance = 2.00
co_triple_bond_length = 1.13

c_coordinate = [rearranged_sites[0].coords[0], rearranged_sites[0].coords[1], rearranged_sites[0].coords[2] + co_floating_distance]
o_coordinate = [rearranged_sites[0].coords[0], rearranged_sites[0].coords[1], rearranged_sites[0].coords[2] + co_floating_distance + co_triple_bond_length]

fe_surface_with_co = fe_surface.append(species="O", coords=o_coordinate, coords_are_cartesian=True)
fe_surface_with_co = fe_surface.append(species="C", coords=c_coordinate, coords_are_cartesian=True)

fe_surface_with_co.to(fmt="cif", filename="fe_100_with_CO.cif")
import os; os.system("explorer .")
print(f"File saved in {os.getcwd()}")

# Reaction going on
# First, the CO needs to shift from directly above an Fe atom to the hollow site

co_hollow_floating_distance = 0.80

c_hollow_coordinate = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + co_hollow_floating_distance]
o_hollow_coordinate = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + co_hollow_floating_distance + co_triple_bond_length]

################################################################################

slab_reaction_initial = slab_gen.get_slabs()
fe_surface_reaction_initial = slab_reaction_initial[0]

fe_surface_reaction_initial.append(species="O", coords=o_hollow_coordinate, coords_are_cartesian=True)
fe_surface_reaction_initial.append(species="C", coords=c_hollow_coordinate, coords_are_cartesian=True)

fe_surface_reaction_initial.append(species="O", coords=[2.87, 0.00, rearranged_sites[1].coords[2] + co_hollow_floating_distance], coords_are_cartesian=True)
fe_surface_reaction_initial.append(species="C", coords=[2.87, 0.00, rearranged_sites[1].coords[2] + co_hollow_floating_distance + co_triple_bond_length], coords_are_cartesian=True)


fe_surface_reaction_initial.to(fmt="cif", filename="reaction_process_initial.cif") # Not fe_surface reaction as the variables hve all been wiped out
import os; os.system("explorer .")
print(f"File saved in {os.getcwd()}")

# Simulate bond dissociation and formation - disproportionation of CO

# Reaction is CO(g) --> CO2(g) + C(ad)

slab_reaction = slab_gen.get_slabs()
fe_surface_reaction = slab_reaction[0]

c_hollow_floating_distance = 0.50
carbon_atom_adsorbed_coordinates = [rearranged_sites[1].coords[0], rearranged_sites[1].coords[1], rearranged_sites[1].coords[2] + c_hollow_floating_distance]

fe_surface_reaction.append(species="C", coords=carbon_atom_adsorbed_coordinates, coords_are_cartesian=True)

carbon_dioxide_floating_distance = 4.5

co2_center_z = rearranged_sites[1].coords[2] + carbon_dioxide_floating_distance  # Floating safely above the surface

co2_c_coordinates = [rearranged_sites[0].coords[0] + 2.87, rearranged_sites[0].coords[1], co2_center_z]                             # Carbon center
co2_o1_coordinates = [rearranged_sites[0].coords[0] + 2.87, rearranged_sites[0].coords[1], co2_center_z - co_triple_bond_length]    # Bottom Oxygen
co2_o2_coordinates = [rearranged_sites[0].coords[0] + 2.87, rearranged_sites[0].coords[1], co2_center_z + co_triple_bond_length]    # Top Oxygen

# I set the (x,y) coordinates of the CO2 molecule to be (2.87, 0.0) because that's going to keep the movie running smoothly like a real simulation, setting it to (0,0) is a mismatch

fe_surface_reaction.append(species="C", coords=co2_c_coordinates, coords_are_cartesian=True)
fe_surface_reaction.append(species="O", coords=co2_o1_coordinates, coords_are_cartesian=True)
fe_surface_reaction.append(species="O", coords=co2_o2_coordinates, coords_are_cartesian=True)

fe_surface_reaction.to(fmt="cif", filename="reaction_process.cif")
import os; os.system("explorer .")
print(f"File saved in {os.getcwd()}")

# Making the movie #

def make_movie(start, end, steps = 5):
    import os
    
    folder = "./reaction_movie"
    os.makedirs(folder, exist_ok=True)
    total = steps + 2
    for idx in range(total):
        frac = idx/(total - 1)
        frame = start.copy()
        for atom_idx in range(len(frame)):
            frame [atom_idx].coords = start[atom_idx].coords + frac * (end[atom_idx].coords - start[atom_idx].coords)
        frame.to("cif", filename=f"{folder}/frame_{idx}.cif")
    print("Movie saved in {folder}")

make_movie(fe_surface_reaction_initial, fe_surface_reaction)