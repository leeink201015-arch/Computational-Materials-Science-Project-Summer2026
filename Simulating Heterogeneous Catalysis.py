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

# Adsorption - basically add carbon and oxygen coordinates to our rearranged_sites
rearranged_sites.append(c_coordinate[0])